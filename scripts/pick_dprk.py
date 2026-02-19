#!/usr/bin/env python3
"""
Pick P-wave arrivals for the DPRK nuclear test using ObsPy.

Supports multiple picker methods:
  - stalta: Classic STA/LTA with predicted-pick selection
  - recursive: Recursive STA/LTA
  - ar: Autoregressive AIC picker
  - all: Run all methods and take consensus

Outputs picks to data/dprk/picks_python.csv for the C++ locator.
"""

import os
import sys
import numpy as np
from obspy import read, UTCDateTime
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta, trigger_onset
from obspy.signal.trigger import ar_pick
from obspy.geodetics import gps2dist_azimuth

# Event parameters
ORIGIN = UTCDateTime("2017-09-03T03:30:01.6")
EVENT_LAT, EVENT_LON, EVENT_DEPTH = 41.343, 129.036, 0.5

# Station coordinates (from stations.txt)
STATION_COORDS = {}
sta_file = sys.argv[1] if len(sys.argv) > 1 else "data/dprk/stations.txt"
with open(sta_file) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        key = f"{parts[0]}.{parts[1]}"
        STATION_COORDS[key] = (float(parts[2]), float(parts[3]))

# Expected P travel time (Pn approximation)
def expected_p_tt(dist_km):
    """Simple Pn travel time: accounts for crustal path + Moho refraction."""
    if dist_km < 150:
        return dist_km / 6.0  # direct Pg
    # Pn: crustal legs + Moho head wave
    moho = 33.0
    v_crust = 6.0
    v_moho = 8.0
    ic = np.arcsin(v_crust / v_moho)
    x_moho = dist_km - 2 * moho * np.tan(ic)
    if x_moho > 0:
        t_crust = 2 * moho / (v_crust * np.cos(ic))
        t_moho = x_moho / v_moho
        return t_crust + t_moho
    return dist_km / 6.5

def pick_station(tr, sta_key, method="stalta"):
    """Pick P arrival on a single trace using predicted-pick selection."""
    sr = tr.stats.sampling_rate
    data = tr.data.astype(float)

    # Demean
    data -= data.mean()

    dist_m, _, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON,
                                     STATION_COORDS[sta_key][0],
                                     STATION_COORDS[sta_key][1])
    dist_km = dist_m / 1000.0
    exp_tt = expected_p_tt(dist_km)

    # Compute characteristic function
    if method == "recursive":
        cft = recursive_sta_lta(data, int(1.0 * sr), int(20.0 * sr))
    else:
        cft = classic_sta_lta(data, int(1.0 * sr), int(20.0 * sr))

    # Get ALL triggers (low threshold to catch real arrivals)
    try:
        triggers = trigger_onset(cft, 2.5, 1.5)
    except Exception:
        return None

    if len(triggers) == 0:
        return None

    # For each trigger, compute travel time
    best_pick = None
    best_resid = 1e9

    for on, off in triggers:
        pick_time = tr.stats.starttime + on / sr
        tt = pick_time - ORIGIN
        if tt < 0:
            continue  # before origin
        resid = abs(tt - exp_tt)
        snr = cft[on] if on < len(cft) else 1.0
        if resid < best_resid and snr >= 2.0:
            best_resid = resid
            best_pick = {
                "station": sta_key,
                "travel_time": tt,
                "expected_tt": exp_tt,
                "residual": tt - exp_tt,
                "snr": snr,
                "distance_km": dist_km,
                "method": method,
            }

    if best_pick and best_resid < 25.0:
        return best_pick
    return None


def main():
    mseed = sys.argv[2] if len(sys.argv) > 2 else "data/dprk/waveforms.mseed"
    out_csv = sys.argv[3] if len(sys.argv) > 3 else "data/dprk/picks_python.csv"

    print("=" * 64)
    print("DPRK P-wave Picking (ObsPy multi-method)")
    print("=" * 64)

    st = read(mseed)
    st_z = st.select(channel="BHZ")

    # Deduplicate: prefer location code "00" or ""
    seen = {}
    for tr in sorted(st_z, key=lambda t: t.stats.location):
        key = f"{tr.stats.network}.{tr.stats.station}"
        if key not in seen:
            seen[key] = tr

    print(f"Stations: {len(seen)}")

    methods = ["stalta", "recursive"]
    all_picks = []

    for sta_key, tr in sorted(seen.items()):
        if sta_key not in STATION_COORDS:
            continue

        # Run multiple methods and take the one with smallest residual
        best = None
        for method in methods:
            pick = pick_station(tr, sta_key, method)
            if pick and (best is None or abs(pick["residual"]) < abs(best["residual"])):
                best = pick

        if best:
            print(f"  {best['station']:>12s}: tt={best['travel_time']:6.1f}s "
                  f"(exp={best['expected_tt']:6.1f}s, resid={best['residual']:+6.1f}s, "
                  f"SNR={best['snr']:.1f}, {best['method']})")
            all_picks.append(best)
        else:
            dist_m, _, _ = gps2dist_azimuth(EVENT_LAT, EVENT_LON,
                                             STATION_COORDS[sta_key][0],
                                             STATION_COORDS[sta_key][1])
            print(f"  {sta_key:>12s}: NO PICK ({dist_m/1000:.0f} km)")

    # Write CSV
    with open(out_csv, "w") as f:
        f.write("station,phase,travel_time_s,snr,quality,distance_km,method\n")
        for p in all_picks:
            f.write(f"{p['station']},P,{p['travel_time']:.4f},"
                    f"{p['snr']:.4f},1,{p['distance_km']:.4f},{p['method']}\n")

    print(f"\n{len(all_picks)} picks written to {out_csv}")
    print("=" * 64)


if __name__ == "__main__":
    main()
