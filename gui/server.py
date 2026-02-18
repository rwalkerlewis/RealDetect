"""
RealDetect Web GUI Server

FastAPI application that serves a tiled map interface showing
seismic stations and detected events in real time.
"""

import argparse
import asyncio
import json
import math
import os
import random
import re
import sqlite3
import subprocess
import time
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.requests import Request


@asynccontextmanager
async def lifespan(application: FastAPI):
    task = asyncio.create_task(poll_events())
    yield
    task.cancel()


app = FastAPI(title="RealDetect", description="Real-time Seismic Monitoring GUI",
              lifespan=lifespan)

BASE_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = BASE_DIR / "templates"
STATIC_DIR = BASE_DIR / "static"

templates = Jinja2Templates(directory=str(TEMPLATES_DIR))
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

# Runtime configuration (set from CLI args or env)
CONFIG = {
    "stations_file": os.environ.get("REALDETECT_STATIONS", "config/stations.txt"),
    "config_file": os.environ.get("REALDETECT_CONFIG", "config/realdetect.conf"),
    "db_file": os.environ.get("REALDETECT_DB", "data/realdetect_catalog.db"),
    "build_dir": os.environ.get("REALDETECT_BUILD_DIR", "build"),
}


# ---------------------------------------------------------------------------
# Station file parser
# ---------------------------------------------------------------------------

def parse_stations(filepath: str) -> list[dict]:
    """Parse the station inventory file into a list of dicts."""
    stations = []
    try:
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    stations.append({
                        "network": parts[0],
                        "station": parts[1],
                        "latitude": float(parts[2]),
                        "longitude": float(parts[3]),
                        "elevation": float(parts[4]) if len(parts) >= 5 else 0.0,
                    })
    except FileNotFoundError:
        pass
    return stations


def parse_config(filepath: str) -> dict:
    """Parse the realdetect.conf into a nested dict."""
    config: dict = {}
    section = "general"
    try:
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                m = re.match(r"^\[(\w+)]$", line)
                if m:
                    section = m.group(1)
                    continue
                m = re.match(r"^(\w+)\s*=\s*(.+)$", line)
                if m:
                    key = f"{section}.{m.group(1)}"
                    config[key] = m.group(2).strip()
    except FileNotFoundError:
        pass
    return config


# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------

def query_events_from_db(db_path: str, limit: int = 200) -> list[dict]:
    """Query recent events from the CSS3.0 SQLite database."""
    events = []
    if not os.path.exists(db_path):
        return events
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row

        cursor = conn.execute("""
            SELECT o.orid, o.lat, o.lon, o.depth, o.time AS otime,
                   o.nass, o.ndef, o.dtype,
                   e.evid, e.prefor
            FROM origin o
            LEFT JOIN event e ON e.prefor = o.orid
            ORDER BY o.time DESC
            LIMIT ?
        """, (limit,))
        for row in cursor:
            mag_val = None
            mag_type = None
            try:
                mag_row = conn.execute(
                    "SELECT magnitude, magtype FROM netmag WHERE orid = ? LIMIT 1",
                    (row["orid"],)
                ).fetchone()
                if mag_row:
                    mag_val = mag_row["magnitude"]
                    mag_type = mag_row["magtype"]
            except Exception:
                pass

            arrivals = []
            try:
                arr_cursor = conn.execute("""
                    SELECT a.sta, a.time AS atime, a.iphase, a.deltim, a.snr,
                           assoc.delta, assoc.seaz, assoc.timeres
                    FROM arrival a
                    JOIN assoc ON assoc.arid = a.arid AND assoc.orid = ?
                    ORDER BY a.time
                """, (row["orid"],))
                for ar in arr_cursor:
                    arrivals.append({
                        "station": ar["sta"],
                        "time": ar["atime"],
                        "phase": ar["iphase"],
                        "residual": ar["timeres"],
                        "distance": ar["delta"],
                        "azimuth": ar["seaz"],
                        "snr": ar["snr"],
                    })
            except Exception:
                pass

            events.append({
                "evid": row["evid"] or row["orid"],
                "orid": row["orid"],
                "latitude": row["lat"],
                "longitude": row["lon"],
                "depth": row["depth"],
                "time": row["otime"],
                "num_phases": row["nass"] or 0,
                "num_stations": row["ndef"] or 0,
                "magnitude": mag_val,
                "magnitude_type": mag_type,
                "arrivals": arrivals,
            })
        conn.close()
    except Exception as exc:
        print(f"DB query error: {exc}")
    return events


def get_db_stats(db_path: str) -> dict:
    """Get summary statistics from the database."""
    stats = {"event_count": 0, "arrival_count": 0, "origin_count": 0}
    if not os.path.exists(db_path):
        return stats
    try:
        conn = sqlite3.connect(db_path)
        for table, key in [("event", "event_count"), ("arrival", "arrival_count"),
                           ("origin", "origin_count")]:
            try:
                row = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()
                stats[key] = row[0] if row else 0
            except Exception:
                pass
        conn.close()
    except Exception:
        pass
    return stats


# ---------------------------------------------------------------------------
# WebSocket manager for real-time push
# ---------------------------------------------------------------------------

class ConnectionManager:
    def __init__(self):
        self.active: list[WebSocket] = []

    async def connect(self, ws: WebSocket):
        await ws.accept()
        self.active.append(ws)

    def disconnect(self, ws: WebSocket):
        if ws in self.active:
            self.active.remove(ws)

    async def broadcast(self, message: dict):
        dead = []
        for ws in self.active:
            try:
                await ws.send_json(message)
            except Exception:
                dead.append(ws)
        for ws in dead:
            self.disconnect(ws)


manager = ConnectionManager()


# ---------------------------------------------------------------------------
# Background task: poll DB for new events
# ---------------------------------------------------------------------------

async def poll_events():
    """Periodically check the database for new events and push updates."""
    last_count = 0
    while True:
        await asyncio.sleep(2)
        db_path = CONFIG["db_file"]
        stats = get_db_stats(db_path)
        current_count = stats.get("event_count", 0)
        if current_count != last_count:
            last_count = current_count
            events = query_events_from_db(db_path, limit=50)
            await manager.broadcast({
                "type": "events_update",
                "events": events,
                "stats": stats,
            })
        elif manager.active:
            await manager.broadcast({
                "type": "heartbeat",
                "stats": stats,
                "timestamp": time.time(),
            })


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.get("/api/stations")
async def api_stations():
    stations = parse_stations(CONFIG["stations_file"])
    return {"stations": stations, "count": len(stations)}


@app.get("/api/events")
async def api_events(limit: int = 200):
    events = query_events_from_db(CONFIG["db_file"], limit=limit)
    return {"events": events, "count": len(events)}


@app.get("/api/stats")
async def api_stats():
    stats = get_db_stats(CONFIG["db_file"])
    stations = parse_stations(CONFIG["stations_file"])
    stats["station_count"] = len(stations)
    return stats


@app.get("/api/config")
async def api_config():
    return parse_config(CONFIG["config_file"])


@app.websocket("/ws")
async def websocket_endpoint(ws: WebSocket):
    await manager.connect(ws)
    # Send initial data
    stations = parse_stations(CONFIG["stations_file"])
    events = query_events_from_db(CONFIG["db_file"], limit=50)
    stats = get_db_stats(CONFIG["db_file"])
    stats["station_count"] = len(stations)
    await ws.send_json({
        "type": "init",
        "stations": stations,
        "events": events,
        "stats": stats,
    })
    try:
        while True:
            data = await ws.receive_text()
            msg = json.loads(data)
            if msg.get("type") == "ping":
                await ws.send_json({"type": "pong"})
            elif msg.get("type") == "refresh":
                events = query_events_from_db(CONFIG["db_file"], limit=50)
                stats = get_db_stats(CONFIG["db_file"])
                stats["station_count"] = len(stations)
                await ws.send_json({
                    "type": "events_update",
                    "events": events,
                    "stats": stats,
                })
    except WebSocketDisconnect:
        manager.disconnect(ws)


# ---------------------------------------------------------------------------
# Demo data generation — populates the DB with realistic sample events
# ---------------------------------------------------------------------------

DEMO_EVENTS = [
    {
        "name": "DPRK Nuclear Test (2017-09-03)",
        "lat": 41.343, "lon": 129.036, "depth": 0.5,
        "time": 1504405801.6, "mag": 6.3, "magtype": "mb",
        "stations_file": "config/stations_korea.txt",
    },
    {
        "name": "Ridgecrest M7.1 (2019-07-06)",
        "lat": 35.770, "lon": -117.599, "depth": 8.0,
        "time": 1562383619.0, "mag": 7.1, "magtype": "ML",
        "stations_file": "config/stations.txt",
    },
    {
        "name": "Southern California M4.2",
        "lat": 33.95, "lon": -117.75, "depth": 12.3,
        "time": None, "mag": 4.2, "magtype": "ML",
        "stations_file": "config/stations.txt",
    },
    {
        "name": "Anza M3.5",
        "lat": 33.55, "lon": -116.57, "depth": 14.1,
        "time": None, "mag": 3.5, "magtype": "ML",
        "stations_file": "config/stations.txt",
    },
    {
        "name": "San Bernardino M2.8",
        "lat": 34.10, "lon": -117.30, "depth": 6.7,
        "time": None, "mag": 2.8, "magtype": "ML",
        "stations_file": "config/stations.txt",
    },
    {
        "name": "Baja California M5.0",
        "lat": 32.45, "lon": -115.30, "depth": 10.0,
        "time": None, "mag": 5.0, "magtype": "ML",
        "stations_file": "config/stations.txt",
    },
]


def _ensure_css30_schema(db_path: str):
    """Create CSS3.0 tables if they don't exist."""
    os.makedirs(os.path.dirname(db_path) or ".", exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS event (
            evid INTEGER PRIMARY KEY,
            evname TEXT DEFAULT '-',
            prefor INTEGER DEFAULT -1,
            commid INTEGER DEFAULT -1,
            lddate REAL DEFAULT 0
        );
        CREATE TABLE IF NOT EXISTS origin (
            orid INTEGER PRIMARY KEY,
            evid INTEGER DEFAULT -1,
            jdate INTEGER DEFAULT -1,
            grn INTEGER DEFAULT -1,
            srn INTEGER DEFAULT -1,
            etype TEXT DEFAULT '-',
            lat REAL DEFAULT -999.0,
            lon REAL DEFAULT -999.0,
            depth REAL DEFAULT -999.0,
            time REAL DEFAULT -9999999999.0,
            nass INTEGER DEFAULT -1,
            ndef INTEGER DEFAULT -1,
            ndp INTEGER DEFAULT -1,
            dtype TEXT DEFAULT '-',
            algorithm TEXT DEFAULT '-',
            auth TEXT DEFAULT '-',
            commid2 INTEGER DEFAULT -1,
            lddate REAL DEFAULT 0
        );
        CREATE TABLE IF NOT EXISTS origerr (
            orid INTEGER PRIMARY KEY,
            sxx REAL DEFAULT -1, syy REAL DEFAULT -1, szz REAL DEFAULT -1,
            stt REAL DEFAULT -1, sxy REAL DEFAULT -1, sxz REAL DEFAULT -1,
            syz REAL DEFAULT -1, stx REAL DEFAULT -1, sty REAL DEFAULT -1,
            stz REAL DEFAULT -1, sdobs REAL DEFAULT -1, smajax REAL DEFAULT -1,
            sminax REAL DEFAULT -1, strike REAL DEFAULT -1,
            sdepth REAL DEFAULT -1, stime REAL DEFAULT -1,
            conf REAL DEFAULT 0.90, commid INTEGER DEFAULT -1,
            lddate REAL DEFAULT 0
        );
        CREATE TABLE IF NOT EXISTS arrival (
            arid INTEGER PRIMARY KEY,
            sta TEXT DEFAULT '-',
            time REAL DEFAULT -9999999999.0,
            endtime REAL DEFAULT -9999999999.0,
            nsamp INTEGER DEFAULT -1,
            samprate REAL DEFAULT -1.0,
            chan TEXT DEFAULT '-',
            iphase TEXT DEFAULT '-',
            stype TEXT DEFAULT '-',
            deltim REAL DEFAULT -1.0,
            azimuth REAL DEFAULT -1.0,
            delaz REAL DEFAULT -1.0,
            slow REAL DEFAULT -1.0,
            delslo REAL DEFAULT -1.0,
            ema REAL DEFAULT -1.0,
            rect REAL DEFAULT -1.0,
            amp REAL DEFAULT -1.0,
            per REAL DEFAULT -1.0,
            logat REAL DEFAULT -1.0,
            clip TEXT DEFAULT '-',
            fm TEXT DEFAULT '-',
            snr REAL DEFAULT -1.0,
            qual TEXT DEFAULT '-',
            auth TEXT DEFAULT '-',
            commid INTEGER DEFAULT -1,
            lddate REAL DEFAULT 0
        );
        CREATE TABLE IF NOT EXISTS assoc (
            arid INTEGER,
            orid INTEGER,
            sta TEXT DEFAULT '-',
            phase TEXT DEFAULT '-',
            belief REAL DEFAULT -1.0,
            delta REAL DEFAULT -1.0,
            seaz REAL DEFAULT -1.0,
            esaz REAL DEFAULT -1.0,
            timeres REAL DEFAULT -999.0,
            timedef TEXT DEFAULT 'd',
            azres REAL DEFAULT -999.0,
            azdef TEXT DEFAULT '-',
            slores REAL DEFAULT -999.0,
            slodef TEXT DEFAULT '-',
            emares REAL DEFAULT -999.0,
            wgt REAL DEFAULT -1.0,
            vmodel TEXT DEFAULT '-',
            commid INTEGER DEFAULT -1,
            lddate REAL DEFAULT 0,
            PRIMARY KEY (arid, orid)
        );
        CREATE TABLE IF NOT EXISTS netmag (
            magid INTEGER PRIMARY KEY,
            net TEXT DEFAULT '-',
            orid INTEGER DEFAULT -1,
            evid INTEGER DEFAULT -1,
            magtype TEXT DEFAULT '-',
            nsta INTEGER DEFAULT -1,
            magnitude REAL DEFAULT -999.0,
            uncertainty REAL DEFAULT -1.0,
            auth TEXT DEFAULT '-',
            commid INTEGER DEFAULT -1,
            lddate REAL DEFAULT 0
        );
        CREATE TABLE IF NOT EXISTS site (
            sta TEXT,
            ondate REAL DEFAULT -1,
            offdate REAL DEFAULT -1,
            lat REAL DEFAULT -999.0,
            lon REAL DEFAULT -999.0,
            elev REAL DEFAULT -999.0,
            staname TEXT DEFAULT '-',
            statype TEXT DEFAULT '-',
            refsta TEXT DEFAULT '-',
            dnorth REAL DEFAULT 0.0,
            deast REAL DEFAULT 0.0,
            lddate REAL DEFAULT 0,
            PRIMARY KEY (sta, ondate)
        );
    """)
    conn.commit()
    conn.close()


def _haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    rlat1, rlat2 = math.radians(lat1), math.radians(lat2)
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(rlat1) * math.cos(rlat2) * math.sin(dlon / 2) ** 2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def _azimuth(lat1, lon1, lat2, lon2):
    rlat1, rlat2 = math.radians(lat1), math.radians(lat2)
    dlon = math.radians(lon2 - lon1)
    y = math.sin(dlon) * math.cos(rlat2)
    x = math.cos(rlat1) * math.sin(rlat2) - math.sin(rlat1) * math.cos(rlat2) * math.cos(dlon)
    return (math.degrees(math.atan2(y, x)) + 360) % 360


def generate_demo_db(db_path: str, stations_file: str = "config/stations.txt"):
    """Insert demo events into the CSS3.0 database."""
    _ensure_css30_schema(db_path)
    conn = sqlite3.connect(db_path)
    now = time.time()

    evid = int(conn.execute("SELECT COALESCE(MAX(evid),0) FROM event").fetchone()[0]) + 1
    orid = int(conn.execute("SELECT COALESCE(MAX(orid),0) FROM origin").fetchone()[0]) + 1
    arid = int(conn.execute("SELECT COALESCE(MAX(arid),0) FROM arrival").fetchone()[0]) + 1
    magid = int(conn.execute("SELECT COALESCE(MAX(magid),0) FROM netmag").fetchone()[0]) + 1

    for i, demo in enumerate(DEMO_EVENTS):
        ev_time = demo["time"] if demo["time"] else now - random.uniform(600, 86400)
        st_file = demo.get("stations_file", stations_file)
        stations = parse_stations(st_file)
        if not stations:
            stations = parse_stations(stations_file)

        # Sort stations by distance, pick closest ones
        for s in stations:
            s["dist"] = _haversine(demo["lat"], demo["lon"], s["latitude"], s["longitude"])
        stations.sort(key=lambda s: s["dist"])
        used = stations[:min(12, len(stations))]

        conn.execute(
            "INSERT INTO event (evid, evname, prefor, lddate) VALUES (?,?,?,?)",
            (evid, demo["name"], orid, now),
        )
        conn.execute(
            "INSERT INTO origin (orid, evid, lat, lon, depth, time, nass, ndef, "
            "dtype, algorithm, auth, lddate) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (orid, evid, demo["lat"], demo["lon"], demo["depth"], ev_time,
             len(used), len(used), "f", "geiger", "demo", now),
        )
        conn.execute(
            "INSERT INTO origerr (orid, smajax, sminax, sdepth, stime, conf, lddate) "
            "VALUES (?,?,?,?,?,?,?)",
            (orid, 2.5, 1.8, 3.0, 0.5, 0.90, now),
        )
        conn.execute(
            "INSERT INTO netmag (magid, net, orid, evid, magtype, nsta, magnitude, "
            "uncertainty, auth, lddate) VALUES (?,?,?,?,?,?,?,?,?,?)",
            (magid, "XX", orid, evid, demo["magtype"], len(used),
             demo["mag"], 0.2, "demo", now),
        )

        for s in used:
            dist_km = s["dist"]
            az = _azimuth(demo["lat"], demo["lon"], s["latitude"], s["longitude"])
            p_tt = dist_km / 6.0 + random.gauss(0, 0.3)
            arr_time = ev_time + p_tt
            residual = random.gauss(0, 0.15)
            snr = max(3, demo["mag"] * 5 - math.log10(max(dist_km, 1)) * 3 + random.gauss(0, 2))

            conn.execute(
                "INSERT INTO arrival (arid, sta, time, chan, iphase, snr, auth, lddate) "
                "VALUES (?,?,?,?,?,?,?,?)",
                (arid, s["station"], arr_time, "BHZ", "P", snr, "demo", now),
            )
            conn.execute(
                "INSERT INTO assoc (arid, orid, sta, phase, delta, seaz, timeres, "
                "vmodel, lddate) VALUES (?,?,?,?,?,?,?,?,?)",
                (arid, orid, s["station"], "P", dist_km, az, residual, "iasp91", now),
            )
            arid += 1

        evid += 1
        orid += 1
        magid += 1

    conn.commit()
    conn.close()
    return len(DEMO_EVENTS)


@app.post("/api/demo")
async def api_load_demo():
    """Populate the database with demo events for map testing."""
    db_path = CONFIG["db_file"]
    st_file = CONFIG["stations_file"]
    count = generate_demo_db(db_path, st_file)
    events = query_events_from_db(db_path, limit=50)
    stats = get_db_stats(db_path)
    stations = parse_stations(st_file)
    stats["station_count"] = len(stations)
    await manager.broadcast({
        "type": "events_update",
        "events": events,
        "stats": stats,
    })
    return {"status": "ok", "events_added": count}


@app.post("/api/run-simulation")
async def api_run_simulation():
    """Run the C++ simulator binary if available."""
    build_dir = CONFIG["build_dir"]
    sim_path = os.path.join(build_dir, "realdetect_sim")
    if not os.path.isfile(sim_path):
        sim_path = "realdetect_sim"
    try:
        result = subprocess.run(
            [sim_path], capture_output=True, text=True, timeout=30,
            cwd=build_dir if os.path.isdir(build_dir) else ".",
        )
        return {
            "status": "ok",
            "returncode": result.returncode,
            "stdout": result.stdout[-2000:] if result.stdout else "",
            "stderr": result.stderr[-1000:] if result.stderr else "",
        }
    except FileNotFoundError:
        return {"status": "error", "message": "Simulator binary not found. Build the project first."}
    except subprocess.TimeoutExpired:
        return {"status": "error", "message": "Simulation timed out."}


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="RealDetect Web GUI")
    parser.add_argument("--stations", default=CONFIG["stations_file"],
                        help="Station inventory file")
    parser.add_argument("--config", default=CONFIG["config_file"],
                        help="RealDetect configuration file")
    parser.add_argument("--db", default=CONFIG["db_file"],
                        help="CSS3.0 SQLite database file")
    parser.add_argument("--build-dir", default=CONFIG["build_dir"],
                        help="C++ build directory")
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8080)
    args = parser.parse_args()

    CONFIG["stations_file"] = args.stations
    CONFIG["config_file"] = args.config
    CONFIG["db_file"] = args.db
    CONFIG["build_dir"] = args.build_dir

    import uvicorn
    uvicorn.run(app, host=args.host, port=args.port, log_level="info")


if __name__ == "__main__":
    main()
