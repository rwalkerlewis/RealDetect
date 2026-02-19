#!/usr/bin/env python3
"""
Generate publication-quality plots from the Ridgecrest M7.1 example output.

Produces five figures:
  1. Station & event map
  2. Record section (waveforms sorted by distance with P/S picks)
  3. STA/LTA characteristic function
  4. Travel-time curve
  5. Summary dashboard (velocity model, residuals, location comparison)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.gridspec import GridSpec

# ── Configuration ───────────────────────────────────────────────────
DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "build/output/ridgecrest"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "output/plots"
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

TRUTH_LAT = 35.770
TRUTH_LON = -117.599
TRUTH_DEPTH = 8.0

# ── Load data ───────────────────────────────────────────────────────
stations = pd.read_csv(os.path.join(DATA_DIR, "stations.csv"))
picks = pd.read_csv(os.path.join(DATA_DIR, "picks.csv"))
waveforms = pd.read_csv(os.path.join(DATA_DIR, "waveforms.csv"))
summary = pd.read_csv(os.path.join(DATA_DIR, "summary.csv"))
vmodel = pd.read_csv(os.path.join(DATA_DIR, "velocity_model.csv"))
stalta = pd.read_csv(os.path.join(DATA_DIR, "stalta.csv"))

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
# FIGURE 1: Station & Event Map
# =====================================================================
fig1, ax = plt.subplots(figsize=(10, 8))
ax.set_title("Ridgecrest M7.1 — Station Network & Location", pad=14)

# Stations
ax.scatter(stations.longitude, stations.latitude, s=60, marker="^",
           c=ACCENT, edgecolors=ACCENT2, linewidths=0.6, zorder=5,
           label="Stations")
for _, r in stations.iterrows():
    ax.annotate(r.station, (r.longitude, r.latitude),
                fontsize=6, color=DIM_CLR, xytext=(4, 4),
                textcoords="offset points")

# True location
ax.scatter(TRUTH_LON, TRUTH_LAT, s=250, marker="*", c=RED,
           edgecolors="white", linewidths=0.8, zorder=10,
           label=f"USGS epicentre ({TRUTH_LAT:.2f}°, {TRUTH_LON:.2f}°)")

# Computed location
ax.scatter(comp_lon, comp_lat, s=180, marker="D", c=GREEN,
           edgecolors="white", linewidths=0.8, zorder=10,
           label=f"RealDetect ({comp_lat:.3f}°, {comp_lon:.3f}°)\n  error = {epi_err:.1f} km")

# Lines from event to stations
for _, r in stations.iterrows():
    ax.plot([comp_lon, r.longitude], [comp_lat, r.latitude],
            color=ACCENT, alpha=0.15, linewidth=0.5, zorder=1)

ax.set_xlabel("Longitude (°)")
ax.set_ylabel("Latitude (°)")
ax.legend(loc="lower left", fontsize=8, framealpha=0.8,
          facecolor=PANEL_BG, edgecolor=GRID_CLR)
ax.grid(True, linewidth=0.3)
ax.set_aspect("equal")
fig1.tight_layout()
fig1.savefig(os.path.join(OUT_DIR, "01_station_map.png"), dpi=180)
print(f"  Saved 01_station_map.png")

# =====================================================================
# FIGURE 2: Record Section (waveforms + picks)
# =====================================================================
fig2, ax = plt.subplots(figsize=(12, 10))
ax.set_title("Ridgecrest M7.1 — Record Section (waveforms sorted by distance)", pad=14)

stas_in_wf = waveforms.station.unique()
# sort by distance
sta_dist = waveforms.groupby("station").distance_km.first().sort_values()

p_picks = picks[picks.phase == "P"]
s_picks = picks[picks.phase == "S"]

plotted = 0
for i, (sta_name, dist) in enumerate(sta_dist.items()):
    if plotted >= 12:
        break
    subset = waveforms[waveforms.station == sta_name]
    t = subset.time_s.values
    amp = subset.amplitude.values
    # normalize per trace
    mx = np.max(np.abs(amp)) if np.max(np.abs(amp)) > 0 else 1
    norm_amp = amp / mx * 0.4

    ax.fill_between(t, dist + norm_amp, dist, where=norm_amp > 0,
                    facecolor=ACCENT, alpha=0.25, linewidth=0)
    ax.fill_between(t, dist + norm_amp, dist, where=norm_amp < 0,
                    facecolor=RED, alpha=0.15, linewidth=0)
    ax.plot(t, dist + norm_amp, color=ACCENT, linewidth=0.4, alpha=0.8)

    ax.text(-9.5, dist, sta_name.split(".")[-1], fontsize=7,
            color=TEXT_CLR, va="center", ha="right")

    # P pick markers
    p_row = p_picks[p_picks.station == sta_name]
    if not p_row.empty:
        p_t = p_row.travel_time_s.values[0]
        ax.axvline(x=p_t - dist / 6.0 + dist / 6.0,  # pick is already travel-time
                   ymin=0, ymax=1, color="white", alpha=0)
        # plot as offset from expected
        ax.plot(p_t, dist, "|", color=GREEN, markersize=10, markeredgewidth=1.5, zorder=6)

    # S pick markers
    s_row = s_picks[s_picks.station == sta_name]
    if not s_row.empty:
        s_t = s_row.travel_time_s.values[0]
        ax.plot(s_t, dist, "|", color=ORANGE, markersize=10, markeredgewidth=1.5, zorder=6)

    plotted += 1

# Theoretical travel-time curves
dists_th = np.linspace(0, 400, 100)
ax.plot(dists_th / 6.0, dists_th, "--", color=GREEN, linewidth=1.2,
        alpha=0.7, label="P (Vp≈6.0 km/s)")
ax.plot(dists_th / 3.46, dists_th, "--", color=ORANGE, linewidth=1.2,
        alpha=0.7, label="S (Vs≈3.5 km/s)")

ax.set_xlabel("Time relative to P arrival (s)")
ax.set_ylabel("Epicentral distance (km)")
ax.set_xlim(-10, 55)
ax.legend(loc="lower right", fontsize=8, framealpha=0.8,
          facecolor=PANEL_BG, edgecolor=GRID_CLR)
ax.grid(True, linewidth=0.3)
ax.invert_yaxis()
fig2.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "02_record_section.png"), dpi=180)
print(f"  Saved 02_record_section.png")

# =====================================================================
# FIGURE 3: STA/LTA Characteristic Function
# =====================================================================
fig3, axes3 = plt.subplots(nrows=min(6, stalta.station.nunique()), ncols=1,
                            figsize=(12, 8), sharex=True)
fig3.suptitle("Ridgecrest M7.1 — STA/LTA Trigger Function", fontweight="bold", y=0.98)

sta_groups = list(stalta.groupby("station"))
sta_groups.sort(key=lambda x: x[1].distance_km.iloc[0])

for idx, (sta_name, grp) in enumerate(sta_groups[:min(6, len(sta_groups))]):
    ax = axes3[idx] if hasattr(axes3, "__len__") else axes3
    t = grp["sample"].values
    ratio = grp["ratio"].values
    dist = grp.distance_km.iloc[0]

    ax.fill_between(t, ratio, 0, alpha=0.3, color=ACCENT)
    ax.plot(t, ratio, color=ACCENT, linewidth=0.8)
    ax.axhline(y=3.5, color=RED, linewidth=1, linestyle="--", alpha=0.8)
    ax.set_ylabel("Ratio", fontsize=8)
    ax.text(0.02, 0.85, f"{sta_name.split('.')[-1]}  ({dist:.0f} km)",
            transform=ax.transAxes, fontsize=8, color=TEXT_CLR,
            bbox=dict(facecolor=PANEL_BG, edgecolor=GRID_CLR, alpha=0.9))
    ax.set_ylim(0, max(8, np.percentile(ratio[np.isfinite(ratio)], 99) * 1.2))
    ax.grid(True, linewidth=0.3)

if hasattr(axes3, "__len__"):
    axes3[-1].set_xlabel("Time relative to expected P (s)")
else:
    axes3.set_xlabel("Time relative to expected P (s)")

fig3.tight_layout(rect=[0, 0, 1, 0.96])
fig3.savefig(os.path.join(OUT_DIR, "03_stalta_function.png"), dpi=180)
print(f"  Saved 03_stalta_function.png")

# =====================================================================
# FIGURE 4: Travel-Time Curve
# =====================================================================
fig4, ax4 = plt.subplots(figsize=(10, 7))
ax4.set_title("Ridgecrest M7.1 — Travel-Time Curve", pad=14)

p_data = picks[picks.phase == "P"]
s_data = picks[picks.phase == "S"]

ax4.scatter(p_data.distance_km, p_data.travel_time_s, s=40, c=GREEN,
            edgecolors="white", linewidths=0.4, zorder=5, label="P picks")
ax4.scatter(s_data.distance_km, s_data.travel_time_s, s=40, c=ORANGE,
            edgecolors="white", linewidths=0.4, zorder=5, label="S picks")

# Theoretical curves
d = np.linspace(0, max(picks.distance_km) * 1.1, 200)
hypo = np.sqrt(d**2 + TRUTH_DEPTH**2)
ax4.plot(d, hypo / 6.0, "--", color=GREEN, alpha=0.6, linewidth=1.5, label="P (6.0 km/s)")
ax4.plot(d, hypo / 6.7, ":", color=GREEN, alpha=0.4, linewidth=1, label="P (6.7 km/s)")
ax4.plot(d, hypo / 3.46, "--", color=ORANGE, alpha=0.6, linewidth=1.5, label="S (3.46 km/s)")
ax4.plot(d, hypo / 3.87, ":", color=ORANGE, alpha=0.4, linewidth=1, label="S (3.87 km/s)")

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
ax5a.set_title("SoCal Velocity Model")
depths = vmodel.depth_km.values
vp = vmodel.vp_km_s.values
vs = vmodel.vs_km_s.values
ax5a.step(vp, depths, color=ACCENT, linewidth=2, where="post", label="Vp")
ax5a.step(vs, depths, color=ORANGE, linewidth=2, where="post", label="Vs")
ax5a.axhline(y=TRUTH_DEPTH, color=RED, linewidth=1, linestyle="--", alpha=0.6,
             label=f"Event depth ({TRUTH_DEPTH} km)")
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
    ax5b.hist(res, bins=20, color=ACCENT, alpha=0.7, edgecolor=ACCENT2, linewidth=0.5)
    ax5b.axvline(x=0, color=RED, linewidth=1.5, linestyle="--")
    ax5b.set_xlabel("Residual (s)")
    ax5b.set_ylabel("Count")
    ax5b.text(0.95, 0.9, f"RMS = {rms:.2f} s",
              transform=ax5b.transAxes, fontsize=9, ha="right",
              color=TEXT_CLR, bbox=dict(facecolor=PANEL_BG, edgecolor=GRID_CLR))
else:
    ax5b.text(0.5, 0.5, "No arrival data", transform=ax5b.transAxes,
              ha="center", va="center", color=DIM_CLR)
ax5b.grid(True, linewidth=0.3)

# 5c: Azimuthal coverage
ax5c = fig5.add_subplot(gs[0, 2], projection="polar")
ax5c.set_title("Azimuthal Coverage", pad=16)
ax5c.set_facecolor(PANEL_BG)
for _, r in stations.iterrows():
    dlat = r.latitude - comp_lat
    dlon = r.longitude - comp_lon
    az = np.arctan2(dlon * np.cos(np.radians(comp_lat)), dlat)
    dist = r.distance_km
    ax5c.scatter(az, dist, s=30, c=ACCENT, edgecolors=ACCENT2, linewidths=0.4, zorder=5)
    ax5c.plot([az, az], [0, dist], color=ACCENT, alpha=0.2, linewidth=0.5)
ax5c.set_rmax(max(stations.distance_km) * 1.1)
ax5c.set_rticks([100, 300, 500, 700, 1000])
ax5c.tick_params(colors=DIM_CLR, labelsize=7)
ax5c.grid(True, linewidth=0.3, color=GRID_CLR, alpha=0.5)

# 5d: Distance vs SNR
ax5d = fig5.add_subplot(gs[1, 0])
ax5d.set_title("P-pick SNR vs Distance")
if not p_data.empty:
    ax5d.scatter(p_data.distance_km, p_data.snr, s=30, c=GREEN,
                 edgecolors="white", linewidths=0.3)
    ax5d.set_xlabel("Distance (km)")
    ax5d.set_ylabel("Signal-to-Noise Ratio")
ax5d.grid(True, linewidth=0.3)

# 5e: Location comparison (zoomed)
ax5e = fig5.add_subplot(gs[1, 1])
ax5e.set_title("Location Comparison (zoomed)")
ax5e.scatter(TRUTH_LON, TRUTH_LAT, s=200, marker="*", c=RED,
             edgecolors="white", linewidths=0.8, zorder=10,
             label=f"USGS")
ax5e.scatter(comp_lon, comp_lat, s=120, marker="D", c=GREEN,
             edgecolors="white", linewidths=0.8, zorder=10,
             label=f"RealDetect")
# error circle
theta = np.linspace(0, 2 * np.pi, 100)
err_r = epi_err / 111.0
ax5e.plot(comp_lon + err_r * np.cos(theta), comp_lat + err_r * np.sin(theta),
          "--", color=GREEN, alpha=0.5, linewidth=1)
ax5e.set_xlabel("Longitude (°)")
ax5e.set_ylabel("Latitude (°)")
ax5e.legend(fontsize=8, framealpha=0.8, facecolor=PANEL_BG, edgecolor=GRID_CLR)
ax5e.set_aspect("equal")
margin = max(0.05, epi_err / 111.0 * 2)
ax5e.set_xlim(min(TRUTH_LON, comp_lon) - margin, max(TRUTH_LON, comp_lon) + margin)
ax5e.set_ylim(min(TRUTH_LAT, comp_lat) - margin, max(TRUTH_LAT, comp_lat) + margin)
ax5e.grid(True, linewidth=0.3)

# 5f: Results text box
ax5f = fig5.add_subplot(gs[1, 2])
ax5f.axis("off")
ax5f.set_title("Pipeline Summary")

results_text = (
    f"Ridgecrest M7.1 Earthquake\n"
    f"2019-07-06 03:19:53 UTC\n"
    f"{'─' * 32}\n"
    f"Stations:          {len(stations)}\n"
    f"P picks:           {len(p_data)}\n"
    f"S picks:           {len(s_data)}\n"
    f"{'─' * 32}\n"
    f"USGS Location:\n"
    f"  {TRUTH_LAT:.3f}°N, {TRUTH_LON:.3f}°E\n"
    f"  Depth: {TRUTH_DEPTH} km\n"
    f"RealDetect Location:\n"
    f"  {comp_lat:.3f}°N, {comp_lon:.3f}°E\n"
    f"  Depth: {float(summary[summary.parameter=='depth'].computed.values[0]):.1f} km\n"
    f"{'─' * 32}\n"
    f"Epicentral error:  {epi_err:.1f} km\n"
    f"RMS residual:      {rms:.2f} s\n"
    f"Magnitude (ML):    {float(summary[summary.parameter=='magnitude'].computed.values[0]):.2f}\n"
)
ax5f.text(0.05, 0.95, results_text, transform=ax5f.transAxes,
          fontsize=9, color=TEXT_CLR, va="top", family="monospace",
          bbox=dict(facecolor=PANEL_BG, edgecolor=GRID_CLR, pad=10))

fig5.savefig(os.path.join(OUT_DIR, "05_summary_dashboard.png"), dpi=180)
print(f"  Saved 05_summary_dashboard.png")

plt.close("all")
print(f"\nAll plots saved to {OUT_DIR}/")
