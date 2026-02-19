#!/usr/bin/env python3
"""
Generate ensemble comparison plots from pipeline comparison CSVs.

Usage:
  python3 scripts/plot_ensemble.py output/ensemble_ridgecrest ridgecrest
  python3 scripts/plot_ensemble.py output/ensemble_dprk dprk
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

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "output/ensemble_ridgecrest"
EVENT = sys.argv[2] if len(sys.argv) > 2 else "ridgecrest"
OUT_DIR = DATA_DIR
os.makedirs(OUT_DIR, exist_ok=True)

# Dark theme
DARK_BG = "#0f1117"
PANEL_BG = "#1a1d26"
TEXT_CLR = "#e0e4f0"
DIM_CLR = "#8890a8"
ACCENT = "#4e9af1"
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
})

# Event parameters
if EVENT == "dprk":
    TRUTH_LAT, TRUTH_LON = 41.343, 129.036
    TITLE = "DPRK Nuclear Test"
else:
    TRUTH_LAT, TRUTH_LON = 35.770, -117.599
    TITLE = "Ridgecrest M7.1"

# Load comparison data
df = pd.read_csv(os.path.join(DATA_DIR, "comparison.csv"))
truth = df[df.config == "USGS_Truth"].iloc[0] if "USGS_Truth" in df.config.values else None
configs = df[df.config != "USGS_Truth"].copy()
configs = configs[configs.error_km >= 0]  # Filter out failed configs

if len(configs) == 0:
    print("No successful configurations to plot")
    sys.exit(0)

print(f"Plotting {len(configs)} configurations for {TITLE}")

# =====================================================================
# FIGURE 1: Location Comparison Map (PyGMT)
# =====================================================================
print("  Generating 01_location_map.png (PyGMT)...")

margin = max(2.0, configs.error_km.max() / 111.0 * 1.5) if len(configs) > 0 else 2.0
all_lons = list(configs.longitude.values) + [TRUTH_LON]
all_lats = list(configs.latitude.values) + [TRUTH_LAT]

region = [
    min(all_lons) - margin,
    max(all_lons) + margin,
    min(all_lats) - margin / 1.5,
    max(all_lats) + margin / 1.5,
]

fig = pygmt.Figure()
fig.coast(
    region=region,
    projection="M15c",
    land="gray90",
    water="lightblue",
    shorelines="0.5p,gray40",
    borders=["1/0.5p,gray50"],
    frame=["af", f"+t{TITLE} — Ensemble Location Comparison"],
    resolution="h",
)

# Plot each configuration's location
colors = ["dodgerblue", "orange", "purple", "cyan", "magenta", "yellow", "green3"]
for i, (_, row) in enumerate(configs.iterrows()):
    c = colors[i % len(colors)]
    fig.plot(x=[row.longitude], y=[row.latitude], style="c0.25c", fill=c, pen="0.5p,black")
    fig.text(x=row.longitude, y=row.latitude,
             text=f"{row.error_km:.0f}km", font="6p,Helvetica,gray20",
             offset="0.2c/0.15c")

# Truth
fig.plot(x=[TRUTH_LON], y=[TRUTH_LAT], style="a0.5c", fill="red", pen="0.8p,darkred")
fig.text(x=TRUTH_LON, y=TRUTH_LAT, text="USGS Truth",
         font="8p,Helvetica-Bold,red", offset="0.3c/-0.3c")

fig.savefig(os.path.join(OUT_DIR, "01_location_map.png"))
print(f"  Saved 01_location_map.png")

# =====================================================================
# FIGURE 2: Error Bar Chart + RMS Comparison
# =====================================================================
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
fig2.suptitle(f"{TITLE} — Ensemble Comparison", fontweight="bold", color=TEXT_CLR, fontsize=14)

names = [c[:25] for c in configs.config.values]
errors = configs.error_km.values
rms_vals = configs.rms_s.values
picks = configs.num_picks.values

# Error bar chart
bars = ax1.barh(range(len(names)), errors, color=ACCENT, alpha=0.8, edgecolor=ACCENT)
ax1.set_yticks(range(len(names)))
ax1.set_yticklabels(names, fontsize=8)
ax1.set_xlabel("Epicentral Error (km)")
ax1.set_title("Location Error by Configuration")
ax1.grid(True, linewidth=0.3, axis='x')
ax1.invert_yaxis()

# Annotate with pick counts
for i, (err, n) in enumerate(zip(errors, picks)):
    ax1.text(err + max(errors) * 0.02, i, f"{err:.1f} km ({n} picks)",
             va='center', fontsize=8, color=TEXT_CLR)

# RMS comparison
bars2 = ax2.barh(range(len(names)), rms_vals, color=ORANGE, alpha=0.8, edgecolor=ORANGE)
ax2.set_yticks(range(len(names)))
ax2.set_yticklabels(names, fontsize=8)
ax2.set_xlabel("RMS Residual (s)")
ax2.set_title("RMS Residual by Configuration")
ax2.grid(True, linewidth=0.3, axis='x')
ax2.invert_yaxis()

for i, rms in enumerate(rms_vals):
    ax2.text(rms + max(rms_vals) * 0.02, i, f"{rms:.2f} s",
             va='center', fontsize=8, color=TEXT_CLR)

fig2.tight_layout(rect=[0, 0, 1, 0.95])
fig2.savefig(os.path.join(OUT_DIR, "02_error_comparison.png"), dpi=180)
print(f"  Saved 02_error_comparison.png")

# =====================================================================
# FIGURE 3: Summary Table
# =====================================================================
fig3, ax3 = plt.subplots(figsize=(14, 6))
ax3.axis("off")
ax3.set_title(f"{TITLE} — Ensemble Results Table", fontweight="bold", pad=20)

header = ["Configuration", "Lat (°N)", "Lon (°E)", "Depth (km)", "RMS (s)", "Error (km)", "Picks"]
rows = []
for _, row in configs.iterrows():
    rows.append([
        row.config[:30],
        f"{row.latitude:.3f}",
        f"{row.longitude:.3f}",
        f"{row.depth_km:.1f}",
        f"{row.rms_s:.2f}",
        f"{row.error_km:.1f}",
        str(int(row.num_picks)),
    ])
rows.append(["USGS Truth", f"{TRUTH_LAT:.3f}", f"{TRUTH_LON:.3f}", "8.0" if EVENT == "ridgecrest" else "0.5", "—", "0.0", "—"])

table = ax3.table(cellText=rows, colLabels=header, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.5)

# Style header
for j in range(len(header)):
    table[0, j].set_facecolor(ACCENT)
    table[0, j].set_text_props(color="white", fontweight="bold")

# Style data rows
for i in range(len(rows)):
    for j in range(len(header)):
        table[i+1, j].set_facecolor(PANEL_BG)
        table[i+1, j].set_text_props(color=TEXT_CLR)
    # Highlight best result
    if i < len(rows) - 1:
        err = configs.iloc[i].error_km
        if err == configs.error_km.min():
            for j in range(len(header)):
                table[i+1, j].set_facecolor("#1a3d26")
    # Highlight truth row
    if i == len(rows) - 1:
        for j in range(len(header)):
            table[i+1, j].set_facecolor("#3d1a1a")

fig3.tight_layout()
fig3.savefig(os.path.join(OUT_DIR, "03_results_table.png"), dpi=180)
print(f"  Saved 03_results_table.png")

plt.close("all")
print(f"\nAll ensemble plots saved to {OUT_DIR}/")
