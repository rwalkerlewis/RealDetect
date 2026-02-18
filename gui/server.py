"""
RealDetect Web GUI Server

FastAPI application that serves a tiled map interface showing
seismic stations and detected events in real time.
"""

import argparse
import asyncio
import json
import os
import re
import sqlite3
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
