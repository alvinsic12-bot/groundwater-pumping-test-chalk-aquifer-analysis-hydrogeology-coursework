# RESET MATPLOTLIB
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import atan2, degrees
from mpl_toolkits.mplot3d import Axes3D


# WELL DATA — OSGB36 + RL + dipper → REAL groundwater head
wells = pd.DataFrame({
    "Site": ["PL10A",  "PL10B",   "PL10E",    "BBA"],
    "E":    [451304.8, 451320.4,  451343.7,   451330.0],
    "N":    [175029.6, 175063.1,  175044.1,   175010.0],
    "RL":   [107.358,  110.102,   107.0315,   105.760],
    "dipper":[21.050,  23.815,    20.825,     19.550]
})

# compute groundwater head
wells["Head"] = wells["RL"] - wells["dipper"]

# extract arrays
E = wells["E"].values
N = wells["N"].values
H = wells["Head"].values

# FIT PLANE: h = a0 + a1*x + a2*y

X = np.column_stack([np.ones_like(E), E, N])
a0, a1, a2 = np.linalg.lstsq(X, H, rcond=None)[0]

# fitted heads and residuals
H_fit = X @ np.array([a0, a1, a2])
wells["H_fit"] = H_fit
wells["Error"] = wells["Head"] - wells["H_fit"]


# HYDRAULIC GRADIENT

grad_mag = np.sqrt(a1*a1 + a2*a2)
flow_bearing = (degrees(atan2(-a1, -a2)) + 360) % 360

cardinal_dirs = ["N","NE","E","SE","S","SW","W","NW"]
cardinal = cardinal_dirs[int((flow_bearing + 22.5)//45) % 8]

# GRID FOR SURFACE + CONTOUR MAP

E_grid, N_grid = np.meshgrid(
    np.linspace(E.min()-20, E.max()+20, 60),
    np.linspace(N.min()-20, N.max()+20, 60)
)
H_grid = a0 + a1*E_grid + a2*N_grid


# FIGURE CANVAS

fig = plt.figure(figsize=(18, 12))


# (a) 3D POTENTIOMETRIC SURFACE (Improved Version)

# coordinate offset for cleaner axes
E0 = E - np.mean(E)
N0 = N - np.mean(N)

# grid offset
E_grid0 = E_grid - np.mean(E)
N_grid0 = N_grid - np.mean(N)

ax1 = fig.add_subplot(2, 2, 1, projection='3d')

ax1.plot_surface(E_grid0, N_grid0, H_grid,
                 color="#66c2a5", alpha=0.75, edgecolor="none")

# observed wells
ax1.scatter(E0, N0, H, color="#fc8d62", s=160, marker="v")

for _, r in wells.iterrows():
    ax1.text(r["E"] - np.mean(E),
             r["N"] - np.mean(N),
             r["Head"] + 0.05,
             f"{r['Site']} ({r['Head']:.2f})",
             fontsize=9)

ax1.view_init(elev=25, azim=215)        # better viewpoint
ax1.set_box_aspect([1, 1, 0.4])         # reduce Z stretch

ax1.set_title("(a) 3D Potentiometric Surface")
ax1.set_xlabel("Easting offset (m)")
ax1.set_ylabel("Northing offset (m)")
ax1.set_zlabel("Head (mAOD)")


# (b) PLAN VIEW

ax2 = fig.add_subplot(2,2,2)

cont = ax2.contour(E_grid, N_grid, H_grid, 12, linewidths=2, colors="#66c2a5")
ax2.clabel(cont, inline=True, fontsize=8)

ax2.scatter(E, N, color="#fc8d62", s=140, marker="v")
for _, r in wells.iterrows():
    ax2.text(r["E"]+1.5, r["N"]+1.5, r["Site"], fontsize=11)

cx, cy = np.mean(E), np.mean(N)
scale = 80

ax2.arrow(cx, cy, -a1*scale*10, -a2*scale*10, width=0.5, color="#66c2a5") 
ax2.arrow(cx, cy,  a1*scale*10,  a2*scale*10, width=0.5, color="#fc8d62") 

ax2.text(
    cx+10, cy+10,
    f"|∇h| = {grad_mag:.4f}\nFlow = {flow_bearing:.1f}° ({cardinal})",
    fontsize=12, bbox=dict(facecolor="white", alpha=0.8)
)

ax2.set_title("(b) Plan View – Contours & Flow Direction")
ax2.set_xlabel("Easting (m)")
ax2.set_ylabel("Northing (m)")
ax2.grid(True, linestyle="--", alpha=0.5)


# (c) EAST–WEST CROSS SECTION  —— FITTED PLANE POINTS ADDED

ax3 = fig.add_subplot(2,2,3)

# sort by Easting
idx_E = np.argsort(E)
E_sorted = E[idx_E]
H_sorted = H[idx_E]
Err_sorted = wells["Error"].iloc[idx_E]

# fitted plane values along EW
EW_fit = a0 + a1*E_sorted + a2*np.mean(N)

# observed points
ax3.scatter(E_sorted, H_sorted, color="#fc8d62", marker="v", s=120,
            label="Observed Head")

# fitted plane points (your request)
ax3.scatter(E_sorted, EW_fit, color="#1b9e77", marker="o", s=100,
            label="Fitted Head (on plane)")

# error band
ax3.fill_between(E_sorted,
                 H_sorted - abs(Err_sorted),
                 H_sorted + abs(Err_sorted),
                 alpha=0.25, color="#fc8d62")

# fitted plane line
ax3.plot(E_sorted, EW_fit, "-", color="#1b9e77", linewidth=2.5,
         label="Fitted Plane")

ax3.set_title("(c) East–West Cross Section")
ax3.set_xlabel("Easting (m)")
ax3.set_ylabel("Head (mAOD)")
ax3.grid(True, linestyle="--", alpha=0.5)
ax3.legend()


# (d) NORTH–SOUTH CROSS SECTION —— FITTED PLANE POINTS ADDED

ax4 = fig.add_subplot(2,2,4)

idx_N = np.argsort(N)
N_sorted = N[idx_N]
H_sorted2 = H[idx_N]
Err_sorted2 = wells["Error"].iloc[idx_N]

# fitted plane values along NS
NS_fit = a0 + a1*np.mean(E) + a2*N_sorted

# observed points
ax4.scatter(N_sorted, H_sorted2, color="#fc8d62", marker="v", s=120,
            label="Observed Head")

# fitted plane head points (your request)
ax4.scatter(N_sorted, NS_fit, color="#1b9e77", marker="o", s=100,
            label="Fitted Head (on plane)")

# error band
ax4.fill_between(N_sorted,
                 H_sorted2 - abs(Err_sorted2),
                 H_sorted2 + abs(Err_sorted2),
                 alpha=0.25, color="#fc8d62")

# fitted plane line
ax4.plot(N_sorted, NS_fit, "-", color="#1b9e77", linewidth=2.5,
         label="Fitted Plane")

ax4.set_title("(d) North–South Cross Section")
ax4.set_xlabel("Northing (m)")
ax4.set_ylabel("Head (mAOD)")
ax4.grid(True, linestyle="--", alpha=0.5)
ax4.legend()


plt.savefig(r"E:\course work\hydrogeology coursework\2x2_hydraulic_gradient.png", dpi=300, bbox_inches="tight")
plt.tight_layout()
plt.show()
