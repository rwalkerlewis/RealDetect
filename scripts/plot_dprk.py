#!/usr/bin/env python3
"""
Generate publication-quality plots from the DPRK nuclear test example output.

Uses PyGMT (GMT backend) for map figures with real coastlines.
Uses matplotlib for non-map data plots.

Produces five figures:
  1. Regional station map (PyGMT — real coastlines, country borders)
  2. Record section (matplotlib)
  3. STA/LTA characteristic function (matplotlib)
  4. Travel-time curve (matplotlib)
  5. Summary dashboard (matplotlib + PyGMT inset)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pygmt

# ── Configuration ───────────────────────────────────────────────────
DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "output/dprk"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "output/plots/dprk"
os.makedirs(OUT_DIR, exist_ok=True)

DARK_BG = "#0f1117"
PANEL_BG = "#1a1d26"
TEXT_CLR = "#e0e4f0"
DIM_CLR = "#8890a8"
ACCENT = "#4e9af1"
ACCENT2 = "#7ab8ff"
RED = "#f05050"
GREEN = "#3ddc84"
ORANGE = "#ffb833"
GRID_CLR = "#2e3345"

plt.rcParams.update({
    "figure.facecolor": DARK_BG,
    "axes.facecolor": PANEL_BG,
    "axes.edgecolor": GRID_CLR,
    "axes.labelcolor": TEXT_CLR,
    "text.color": TEXT_CLR,
    "xtick.color": DIM_CLR,
    "ytick.color": DIM_CLR,
    "grid.color": GRID_CLR,
    "grid.alpha": 0.5,
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
})

TRUTH_LAT = 41.343
TRUTH_LON = 129.036
TRUTH_DEPTH = 0.5

# ── Load data ───────────────────────────────────────────────────────
stations = pd.read_csv(os.path.join(DATA_DIR, "stations.csv"))
picks = pd.read_csv(os.path.join(DATA_DIR, "picks.csv"))
waveforms = pd.read_csv(os.path.join(DATA_DIR, "waveforms.csv"))
summary = pd.read_csv(os.path.join(DATA_DIR, "summary.csv"))
vmodel = pd.read_csv(os.path.join(DATA_DIR, "velocity_model.csv"))

try:
    stalta = pd.read_csv(os.path.join(DATA_DIR, "stalta.csv"))
except Exception:
    stalta = pd.DataFrame()

try:
    arrivals = pd.read_csv(os.path.join(DATA_DIR, "arrivals.csv"))
except Exception:
    arrivals = pd.DataFrame()

comp_lat = float(summary[summary.parameter == "latitude"].computed.values[0])
comp_lon = float(summary[summary.parameter == "longitude"].computed.values[0])
epi_err = float(summary[summary.parameter == "epicentral_error_km"].computed.values[0])
rms = float(summary[summary.parameter == "rms_s"].computed.values[0])

print(f"Data loaded from {DATA_DIR}")
print(f"  {len(stations)} stations, {len(picks)} picks")
print(f"  Epicentral error: {epi_err:.1f} km,  RMS: {rms:.2f} s")

# =====================================================================
# FIGURE 1: Regional Station Map (PyGMT)
# =====================================================================
print("  Generating 01_station_map.png (PyGMT)...")

# Focus on East Asia region (exclude very distant stations like ADK)
nearby = stations[stations.distance_km < 2500]
if len(nearby) < 5:
    nearby = stations
lon_min = min(nearby.longitude.min(), TRUTH_LON) - 5
lon_max = max(nearby.longitude.max(), TRUTH_LON) + 5
lat_min = min(nearby.latitude.min(), TRUTH_LAT) - 3
lat_max = max(nearby.latitude.max(), TRUTH_LAT) + 3

fig = pygmt.Figure()
fig.coast(
    region=[lon_min, lon_max, lat_min, lat_max],
    projection="M18c",
    land="gray90",
    water="lightblue",
    shorelines="0.5p,gray40",
    borders=["1/0.8p,gray30", "2/0.3p,gray50"],
    frame=["af", "+tDPRK Nuclear Test — Regional Station Network (EarthScope Data)"],
    resolution="h",
)

# Ray paths
for _, r in stations.iterrows():
    fig.plot(
        x=[TRUTH_LON, r.longitude],
        y=[TRUTH_LAT, r.latitude],
        pen="0.3p,dodgerblue,dashed",
    )

# Station triangles
fig.plot(
    x=stations.longitude.values,
    y=stations.latitude.values,
    style="t0.35c",
    fill="dodgerblue",
    pen="0.5p,black",
)

# Station labels
for _, r in stations.iterrows():
    fig.text(
        x=r.longitude,
        y=r.latitude,
        text=f"{r.network}.{r.station}",
        font="5p,Helvetica,gray30",
        offset="0.2c/0.15c",
    )

# USGS/CTBTO epicenter (red star)
fig.plot(
    x=[TRUTH_LON],
    y=[TRUTH_LAT],
    style="a0.6c",
    fill="red",
    pen="0.8p,darkred",
)
fig.text(
    x=TRUTH_LON,
    y=TRUTH_LAT,
    text="Punggye-ri (USGS)",
    font="8p,Helvetica-Bold,red",
    offset="0.4c/-0.4c",
)

# Computed epicenter
fig.plot(
    x=[comp_lon],
    y=[comp_lat],
    style="d0.5c",
    fill="green3",
    pen="0.8p,darkgreen",
)

fig.savefig(os.path.join(OUT_DIR, "01_station_map.png"))
print(f"  Saved 01_station_map.png")

# =====================================================================
# FIGURE 2: Record Section
# =====================================================================
fig2, ax = plt.subplots(figsize=(12, 10))
ax.set_title("DPRK Nuclear Test — Record Section (real EarthScope waveforms)", pad=14)

sta_dist = waveforms.groupby("station").distance_km.first().sort_values()
p_picks_df = picks[picks.phase == "P"] if "phase" in picks.columns else pd.DataFrame()

plotted = 0
for i, (sta_name, dist) in enumerate(sta_dist.items()):
    if plotted >= 12:
        break
    subset = waveforms[waveforms.station == sta_name]
    t = subset.time_s.values
    amp = subset.amplitude.values
    mx = np.max(np.abs(amp)) if np.max(np.abs(amp)) > 0 else 1
    norm_amp = amp / mx * 0.4

    ax.fill_between(t, dist + norm_amp, dist, where=norm_amp > 0,
                    facecolor=ACCENT, alpha=0.25, linewidth=0)
    ax.fill_between(t, dist + norm_amp, dist, where=norm_amp < 0,
                    facecolor=RED, alpha=0.15, linewidth=0)
    ax.plot(t, dist + norm_amp, color=ACCENT, linewidth=0.4, alpha=0.8)

    ax.text(t.min() - 1, dist, sta_name.split(".")[-1], fontsize=7,
            color=TEXT_CLR, va="center", ha="right")

    if not p_picks_df.empty:
        p_row = p_picks_df[p_picks_df.station == sta_name]
        if not p_row.empty:
            ax.plot(p_row.travel_time_s.values[0], dist, "|",
                    color=GREEN, markersize=10, markeredgewidth=1.5, zorder=6)
    plotted += 1

ax.set_xlabel("Time relative to origin (s)")
ax.set_ylabel("Epicentral distance (km)")
ax.grid(True, linewidth=0.3)
fig2.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "02_record_section.png"), dpi=180)
print(f"  Saved 02_record_section.png")

# =====================================================================
# FIGURE 3: STA/LTA
# =====================================================================
if not stalta.empty and stalta.station.nunique() > 0:
    n_panels = min(6, stalta.station.nunique())
    fig3, axes3 = plt.subplots(nrows=n_panels, ncols=1,
                                figsize=(12, 8), sharex=True)
    fig3.suptitle("DPRK Nuclear Test — STA/LTA Trigger (real data)", fontweight="bold", y=0.98)

    sta_groups = list(stalta.groupby("station"))
    sta_groups.sort(key=lambda x: x[1].distance_km.iloc[0])

    if n_panels == 1:
        axes3 = [axes3]

    for idx, (sta_name, grp) in enumerate(sta_groups[:n_panels]):
        ax = axes3[idx]
        t = grp["sample"].values
        ratio = grp["ratio"].values
        dist = grp.distance_km.iloc[0]

        ax.fill_between(t, ratio, 0, alpha=0.3, color=ACCENT)
        ax.plot(t, ratio, color=ACCENT, linewidth=0.8)
        ax.axhline(y=3.0, color=RED, linewidth=1, linestyle="--", alpha=0.8)
        ax.set_ylabel("Ratio", fontsize=8)
        ax.text(0.02, 0.85, f"{sta_name}  ({dist:.0f} km)",
                transform=ax.transAxes, fontsize=8, color=TEXT_CLR,
                bbox=dict(facecolor=PANEL_BG, edgecolor=GRID_CLR, alpha=0.9))
        finite_ratio = ratio[np.isfinite(ratio)]
        if len(finite_ratio) > 0:
            ax.set_ylim(0, max(8, np.percentile(finite_ratio, 99) * 1.2))
        ax.grid(True, linewidth=0.3)

    axes3[-1].set_xlabel("Time relative to origin (s)")
    fig3.tight_layout(rect=[0, 0, 1, 0.96])
    fig3.savefig(os.path.join(OUT_DIR, "03_stalta_function.png"), dpi=180)
    print(f"  Saved 03_stalta_function.png")
else:
    print("  Skipping 03_stalta_function.png")

# =====================================================================
# FIGURE 4: Travel-Time Curve
# =====================================================================
fig4, ax4 = plt.subplots(figsize=(10, 7))
ax4.set_title("DPRK Nuclear Test — Travel-Time Curve (real data)", pad=14)

if not p_picks_df.empty:
    ax4.scatter(p_picks_df.distance_km, p_picks_df.travel_time_s, s=40, c=GREEN,
                edgecolors="white", linewidths=0.4, zorder=5, label="P picks (real)")

d = np.linspace(0, max(picks.distance_km) * 1.1, 200) if len(picks) > 0 else np.linspace(0, 2000, 200)
hypo = np.sqrt(d**2 + TRUTH_DEPTH**2)
ax4.plot(d, hypo / 6.0, "--", color=GREEN, alpha=0.6, linewidth=1.5, label="Pg (6.0 km/s)")
ax4.plot(d, hypo / 7.8, ":", color=GREEN, alpha=0.4, linewidth=1, label="Pn (7.8 km/s)")

ax4.set_xlabel("Epicentral Distance (km)")
ax4.set_ylabel("Travel Time (s)")
ax4.legend(fontsize=8, framealpha=0.8, facecolor=PANEL_BG, edgecolor=GRID_CLR)
ax4.grid(True, linewidth=0.3)
fig4.tight_layout()
fig4.savefig(os.path.join(OUT_DIR, "04_travel_time.png"), dpi=180)
print(f"  Saved 04_travel_time.png")

# =====================================================================
# FIGURE 5: Summary Dashboard
# =====================================================================
fig5 = plt.figure(figsize=(14, 9))
gs = GridSpec(2, 3, figure=fig5, hspace=0.35, wspace=0.35)

# 5a: Velocity model
ax5a = fig5.add_subplot(gs[0, 0])
ax5a.set_title("Korean Velocity Model")
ax5a.step(vmodel.vp_km_s.values, vmodel.depth_km.values, color=ACCENT, linewidth=2, where="post", label="Vp")
ax5a.step(vmodel.vs_km_s.values, vmodel.depth_km.values, color=ORANGE, linewidth=2, where="post", label="Vs")
ax5a.set_xlabel("Velocity (km/s)")
ax5a.set_ylabel("Depth (km)")
ax5a.invert_yaxis()
ax5a.legend(fontsize=7, framealpha=0.8, facecolor=PANEL_BG, edgecolor=GRID_CLR)
ax5a.grid(True, linewidth=0.3)

# 5b: Residual histogram
ax5b = fig5.add_subplot(gs[0, 1])
ax5b.set_title("Arrival Residuals")
if not arrivals.empty:
    res = arrivals.residual_s.values
    ax5b.hist(res, bins=15, color=ACCENT, alpha=0.7, edgecolor=ACCENT2, linewidth=0.5)
    ax5b.axvline(x=0, color=RED, linewidth=1.5, linestyle="--")
    ax5b.set_xlabel("Residual (s)")
    ax5b.set_ylabel("Count")
ax5b.grid(True, linewidth=0.3)

# 5c: Azimuthal coverage
ax5c = fig5.add_subplot(gs[0, 2], projection="polar")
ax5c.set_title("Azimuthal Coverage", pad=16)
ax5c.set_facecolor(PANEL_BG)
for _, r in stations.iterrows():
    dlat = r.latitude - TRUTH_LAT
    dlon = r.longitude - TRUTH_LON
    az = np.arctan2(dlon * np.cos(np.radians(TRUTH_LAT)), dlat)
    dist = r.distance_km
    ax5c.scatter(az, dist, s=30, c=ACCENT, edgecolors=ACCENT2, linewidths=0.4, zorder=5)
    ax5c.plot([az, az], [0, dist], color=ACCENT, alpha=0.2, linewidth=0.5)
ax5c.set_rmax(max(stations.distance_km) * 1.1)
ax5c.tick_params(colors=DIM_CLR, labelsize=7)
ax5c.grid(True, linewidth=0.3, color=GRID_CLR, alpha=0.5)

# 5d: SNR vs distance
ax5d = fig5.add_subplot(gs[1, 0])
ax5d.set_title("P-pick SNR vs Distance")
if not p_picks_df.empty:
    ax5d.scatter(p_picks_df.distance_km, p_picks_df.snr, s=30, c=GREEN,
                 edgecolors="white", linewidths=0.3)
    ax5d.set_xlabel("Distance (km)")
    ax5d.set_ylabel("SNR")
ax5d.grid(True, linewidth=0.3)

# 5e: Location map (PyGMT)
zoom_margin = max(2.0, epi_err / 111.0 * 2)
zoom_region = [
    min(TRUTH_LON, comp_lon) - zoom_margin,
    max(TRUTH_LON, comp_lon) + zoom_margin,
    min(TRUTH_LAT, comp_lat) - zoom_margin,
    max(TRUTH_LAT, comp_lat) + zoom_margin,
]
fig_zoom = pygmt.Figure()
fig_zoom.coast(
    region=zoom_region,
    projection="M6c",
    land="gray90",
    water="lightblue",
    shorelines="0.3p,gray50",
    borders=["1/0.5p,gray30"],
    frame=["af", "+tLocation Comparison"],
    resolution="h",
)
fig_zoom.plot(x=[TRUTH_LON], y=[TRUTH_LAT], style="a0.4c", fill="red", pen="0.5p,darkred")
fig_zoom.plot(x=[comp_lon], y=[comp_lat], style="d0.35c", fill="green3", pen="0.5p,darkgreen")
zoom_path = os.path.join(OUT_DIR, "05e_location_zoom.png")
fig_zoom.savefig(zoom_path)

ax5e = fig5.add_subplot(gs[1, 1])
ax5e.set_title("Location Comparison (GMT)")
try:
    img = plt.imread(zoom_path)
    ax5e.imshow(img)
except Exception:
    pass
ax5e.axis("off")

# 5f: Results text
ax5f = fig5.add_subplot(gs[1, 2])
ax5f.axis("off")
ax5f.set_title("Pipeline Summary")

depth_val = float(summary[summary.parameter == 'depth'].computed.values[0])
mag_val = float(summary[summary.parameter == 'magnitude'].computed.values[0])

results_text = (
    f"DPRK 6th Nuclear Test\n"
    f"2017-09-03 03:30:01.6 UTC\n"
    f"Real data from EarthScope\n"
    f"{'─' * 32}\n"
    f"Stations:          {len(stations)}\n"
    f"P picks:           {len(p_picks_df)}\n"
    f"{'─' * 32}\n"
    f"USGS Location:\n"
    f"  {TRUTH_LAT:.3f}°N, {TRUTH_LON:.3f}°E\n"
    f"  Depth: {TRUTH_DEPTH} km\n"
    f"RealDetect Location:\n"
    f"  {comp_lat:.3f}°N, {comp_lon:.3f}°E\n"
    f"  Depth: {depth_val:.1f} km\n"
    f"{'─' * 32}\n"
    f"Epicentral error:  {epi_err:.1f} km\n"
    f"RMS residual:      {rms:.2f} s\n"
    f"Magnitude:         {mag_val:.2f}\n"
)
ax5f.text(0.05, 0.95, results_text, transform=ax5f.transAxes,
          fontsize=9, color=TEXT_CLR, va="top", family="monospace",
          bbox=dict(facecolor=PANEL_BG, edgecolor=GRID_CLR, pad=10))

fig5.savefig(os.path.join(OUT_DIR, "05_summary_dashboard.png"), dpi=180)
print(f"  Saved 05_summary_dashboard.png")

plt.close("all")
print(f"\nAll plots saved to {OUT_DIR}/")
