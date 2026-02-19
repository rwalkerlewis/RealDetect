#!/usr/bin/env python3
"""
Fetch a small real data fixture for unit tests.

Uses a subset of Ridgecrest M7.1 data: 5 closest stations, 60s window.
Small enough for test fixtures but real enough for proper validation.
"""

import os
import sys
from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

# ── Event parameters (Ridgecrest M7.1) ──────────────────────────
ORIGIN_TIME = UTCDateTime("2019-07-06T03:19:53.04")
EVENT_LAT = 35.770
EVENT_LON = -117.599

# ── Output ───────────────────────────────────────────────────────
OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "data/test"
os.makedirs(OUT_DIR, exist_ok=True)

# ── 5 closest stations for a compact fixture ─────────────────────
STATIONS = [
    ("CI", "CLC"),   # ~27 km
    ("CI", "SLA"),   # ~40 km
    ("CI", "SRT"),   # ~15 km
    ("CI", "MPM"),   # ~32 km
    ("CI", "WRC2"),  # ~36 km
]


def main():
    print("Fetching test fixture data from EarthScope...")

    scedc = Client("SCEDC")
    iris = Client("IRIS")

    t_start = ORIGIN_TIME - 10.0
    t_end = ORIGIN_TIME + 60.0

    all_streams = Stream()
    station_info = []

    for net, sta in STATIONS:
        client = scedc if net == "CI" else iris
        try:
            inv = client.get_stations(
                network=net, station=sta,
                starttime=t_start, endtime=t_end,
                channel="BHZ", level="station"
            )
            info = None
            for n in inv:
                for s in n:
                    info = {
                        "network": n.code,
                        "station": s.code,
                        "latitude": s.latitude,
                        "longitude": s.longitude,
                        "elevation": s.elevation,
                    }
                    break
            if info is None:
                print(f"  {net}.{sta}: no metadata, skipping")
                continue

            st = client.get_waveforms(net, sta, "*", "BHZ", t_start, t_end)
            if len(st) == 0:
                print(f"  {net}.{sta}: no waveform data, skipping")
                continue

            dist_m, _, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON,
                                             info["latitude"], info["longitude"])
            print(f"  {net}.{sta}: {st[0].stats.npts} samples, "
                  f"{dist_m / 1000:.1f} km")
            all_streams += st
            station_info.append(info)

        except Exception as e:
            print(f"  {net}.{sta}: failed ({e})")

    if len(all_streams) == 0:
        print("ERROR: No data fetched!")
        sys.exit(1)

    all_streams.merge(method=1, fill_value=0)

    mseed_path = os.path.join(OUT_DIR, "waveforms.mseed")
    all_streams.write(mseed_path, format="MSEED")
    print(f"Saved {len(all_streams)} traces to {mseed_path}")

    sta_path = os.path.join(OUT_DIR, "stations.txt")
    with open(sta_path, "w") as f:
        f.write("# Test fixture stations (Ridgecrest M7.1 subset)\n")
        f.write("# Network Station Latitude Longitude Elevation(m)\n")
        for info in station_info:
            f.write(f"{info['network']} {info['station']} "
                    f"{info['latitude']:.4f} {info['longitude']:.4f} "
                    f"{info['elevation']:.1f}\n")
    print(f"Saved {len(station_info)} stations to {sta_path}")
    print("Test fixture data ready.")


if __name__ == "__main__":
    main()
