#!/usr/bin/env python3
"""
Fetch real seismic data for the 2017-09-03 DPRK 6th nuclear test
from EarthScope (IRIS) FDSN web services.

USGS/CTBTO Parameters:
  Date:      2017-09-03
  Time:      03:30:01.6 UTC
  Location:  41.343°N, 129.036°E
  Depth:     ~0 km (shallow explosion)
  Magnitude: mb 6.3
"""

import os
import sys
from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

# ── Event parameters ─────────────────────────────────────────────
ORIGIN_TIME = UTCDateTime("2017-09-03T03:30:01.6")
EVENT_LAT = 41.343
EVENT_LON = 129.036
EVENT_DEPTH_KM = 0.5
EVENT_MAG = 6.3

# ── Data windows ─────────────────────────────────────────────────
PRE_ORIGIN = 15.0
POST_ORIGIN = 180.0

# ── Output ───────────────────────────────────────────────────────
OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "data/dprk"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Target stations (verified available from IRIS) ───────────────
TARGET_STATIONS = [
    ("IU", "INCN"),   # Incheon, South Korea
    ("IC", "BJT"),    # Baijiatuan, Beijing, China
    ("IC", "MDJ"),    # Mudanjiang, China
    ("IC", "SSE"),    # Shanghai, China
    ("IU", "MAJO"),   # Matsushiro, Japan
    ("IU", "YSS"),    # Yuzhno-Sakhalinsk, Russia
    ("II", "ERM"),    # Erimo, Japan
    ("KS", "BUS2"),   # Busan, South Korea
    ("KS", "SEO2"),   # Seoul, South Korea
    ("IU", "ULN"),    # Ulaanbaatar, Mongolia
    ("IU", "ADK"),    # Adak, Alaska
    ("IU", "TATO"),   # Taipei, Taiwan
    ("JP", "JTU"),    # Tsushima, Japan
    ("JP", "JMM"),    # Matsumoto, Japan
    ("JP", "JTM"),    # Towada, Japan
]


def fetch_station_info(client, network, station, t):
    """Fetch station metadata from FDSN."""
    try:
        inv = client.get_stations(
            network=network, station=station,
            starttime=t - 60, endtime=t + 60,
            channel="BHZ", level="station"
        )
        for net in inv:
            for sta in net:
                return {
                    "network": net.code,
                    "station": sta.code,
                    "latitude": sta.latitude,
                    "longitude": sta.longitude,
                    "elevation": sta.elevation,
                }
    except Exception:
        pass
    return None


def fetch_waveform(client, network, station, t_start, t_end):
    """Fetch 3-component waveform data (BH? = BHZ+BHN+BHE or BH1+BH2)."""
    for chan in ["BH?", "HH?"]:
        try:
            st = client.get_waveforms(network, station, "*", chan, t_start, t_end)
            if len(st) > 0:
                return st
        except Exception:
            continue
    return None


def main():
    print("=" * 64)
    print("Fetching DPRK Nuclear Test real data from EarthScope")
    print("=" * 64)
    print(f"Origin: {ORIGIN_TIME}")
    print(f"Location: {EVENT_LAT}°N, {EVENT_LON}°E")
    print(f"Depth: {EVENT_DEPTH_KM} km, Magnitude: mb {EVENT_MAG}")
    print()

    iris = Client("IRIS")

    t_start = ORIGIN_TIME - PRE_ORIGIN
    t_end = ORIGIN_TIME + POST_ORIGIN

    all_streams = Stream()
    station_info = []

    print("Fetching regional stations from IRIS...")
    for net_code, sta_code in TARGET_STATIONS:
        info = fetch_station_info(iris, net_code, sta_code, ORIGIN_TIME)
        if info is None:
            print(f"  {net_code}.{sta_code}: station metadata not found, skipping")
            continue

        st = fetch_waveform(iris, net_code, sta_code, t_start, t_end)
        if st is None:
            print(f"  {net_code}.{sta_code}: no waveform data, skipping")
            continue

        dist_m, _, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON,
                                         info["latitude"], info["longitude"])
        dist_km = dist_m / 1000.0

        print(f"  {net_code}.{sta_code}: {st[0].stats.npts} samples, "
              f"{st[0].stats.sampling_rate} Hz, {dist_km:.1f} km")
        all_streams += st
        station_info.append(info)

    if len(all_streams) == 0:
        print("\nERROR: No waveform data fetched!")
        sys.exit(1)

    # Merge overlapping traces
    all_streams.merge(method=1, fill_value=0)

    # Write MiniSEED
    mseed_path = os.path.join(OUT_DIR, "waveforms.mseed")
    all_streams.write(mseed_path, format="MSEED")
    print(f"\nSaved {len(all_streams)} traces to {mseed_path}")

    # Write station file
    sta_path = os.path.join(OUT_DIR, "stations.txt")
    with open(sta_path, "w") as f:
        f.write("# Station inventory fetched from EarthScope FDSN\n")
        f.write("# Network Station Latitude Longitude Elevation(m)\n")
        for info in station_info:
            f.write(f"{info['network']} {info['station']} "
                    f"{info['latitude']:.4f} {info['longitude']:.4f} "
                    f"{info['elevation']:.1f}\n")
    print(f"Saved {len(station_info)} stations to {sta_path}")

    # Summary
    print(f"\n{'=' * 64}")
    print(f"DPRK data fetch complete")
    print(f"  Stations: {len(station_info)}")
    print(f"  Traces:   {len(all_streams)}")
    print(f"  Window:   {t_start} to {t_end}")
    print(f"  Output:   {OUT_DIR}/")
    print(f"{'=' * 64}")


if __name__ == "__main__":
    main()
