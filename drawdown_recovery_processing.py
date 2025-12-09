import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# 1. Load Excel (MultiIndex header)

file = r"E:\course work\hydrogeology coursework\工作簿1.xlsx"
df = pd.read_excel(file, header=[0,1])     # multi-level header


# 2. Well column keys

# --- time ---
time_col = {
    "PL10A": ("PL10A","Time"),
    "PL10B": ("PL10B","Time"),
    "PL10E": ("PL10E","Time"),
    "BBA_dip": ("BBA ","Time"),
    "BBA_log": ("BBA (logged)","Time"),
    "BBA_obs": ("BBA Obs (logged)","Time")
}

# --- WL or Dip ---
level_col = {
    "PL10A": ("Group 2","Dip (m)"),
    "PL10B": ("Group 1","Dip (m)"),
    "PL10E": ("Groups 1 & 2","Dip (m)"),
    "BBA_dip": ("BBA ","Dip (m)"),
    "BBA_log": ("BBA (logged)","Level (mAOD)"),
    "BBA_obs": ("BBA Obs (logged)","Level (mAOD)")
}

# --- casing elevation (None = logged well that already uses mAOD) ---
casing = {
    "PL10A": 107.358,
    "PL10B": 110.102,
    "PL10E": 107.0315,
    "BBA_dip": 105.757,
    "BBA_log": None,
    "BBA_obs": None
}

# global baseline for non-logged wells
baseline_target = pd.to_datetime("2025-10-22 12:00:00")

draw_color = "#65B2FA"
reco_color = "#FCC24C"



# 3. Processing function

def process_manual_well(well):
    """
    For PL10A / PL10B / PL10E / BBA_dip:
    Use baseline = nearest to 12:00.
    """
    t_col = time_col[well]
    w_col = level_col[well]

    t_raw = df[t_col].dropna()
    wl_raw = df[w_col].loc[t_raw.index]

    # convert time
    t_dt = pd.to_datetime("2025-10-22 " + t_raw.astype(str), errors="coerce")

    # dip → WL
    WL = casing[well] - wl_raw

    d = pd.DataFrame({"time_dt": t_dt, "WL": WL})
    d = d.dropna()

    # find baseline nearest to 12:00
    d["abs_diff"] = (d["time_dt"] - baseline_target).abs()
    t0 = d.loc[d["abs_diff"].idxmin(), "time_dt"]

    # keep only after baseline
    d = d[d["time_dt"] >= t0].copy()
    d["t_min"] = (d["time_dt"] - t0).dt.total_seconds()/60

    # compute drawdown
    static = d["WL"].iloc[0]
    d["s"] = static - d["WL"]
    d["ds"] = d["s"].diff().fillna(0)

    # remove flat / no-response region
    mask = (d["s"].abs() < 0.002) & (d["ds"].abs() < 0.001)
    first_valid = d.index[~mask][0]
    d = d.loc[first_valid:].copy()

    # split by maximum drawdown
    idx_max = d["s"].idxmax()
    return d.loc[:idx_max].copy(), d.loc[idx_max+1:].copy()


def process_logged_well(well):
    """
    For BBA_log & BBA_obs:
    t = cumulative minutes from own start (original sampling interval)
    NO baseline, NO trimming!
    """
    t_col = time_col[well]
    w_col = level_col[well]

    t_raw = df[t_col].dropna()
    wl_raw = pd.to_numeric(df[w_col].loc[t_raw.index], errors="coerce")

    # convert to datetime for computing deltas
    t_dt = pd.to_datetime("2025-10-22 " + t_raw.astype(str), errors="coerce")
    d = pd.DataFrame({"time_dt": t_dt, "WL": wl_raw})
    d = d.dropna().reset_index(drop=True)

    # time since first record
    d["t_min"] = (d["time_dt"] - d["time_dt"].iloc[0]).dt.total_seconds()/60

    # drawdown: relative to first WL
    static = d["WL"].iloc[0]
    d["s"] = static - d["WL"]

    # split at maximum s
    idx_max = d["s"].idxmax()
    return d.loc[:idx_max].copy(), d.loc[idx_max+1:].copy()



# 4. Process all wells

all_data = {}

for w in ["PL10A","PL10B","PL10E","BBA_dip"]:
    all_data[w] = process_manual_well(w)

for w in ["BBA_log","BBA_obs"]:
    all_data[w] = process_logged_well(w)



# 5. 2×3 PANEL PLOT (final perfect output)

fig, axes = plt.subplots(2, 3, figsize=(22, 13))
axes = axes.flatten()
fig.suptitle("Drawdown & Recovery (All 6 Wells)", fontsize=18)

for ax, (well, (dfD, dfR)) in zip(axes, all_data.items()):

    # --- drawdown ---
    ax.plot(dfD["t_min"], dfD["s"], "o", markerfacecolor="none",
            markeredgecolor=draw_color, label="Drawdown")

    # --- recovery ---
    if not dfR.empty:
        ax.plot(dfR["t_min"], dfR["s"], "s", markerfacecolor="none",
                markeredgecolor=reco_color, label="Recovery")

    # tp / sp
    tp = dfD["t_min"].iloc[-1]
    sp = dfD["s"].iloc[-1]

    # extend dashed lines to panel borders (not axes)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

      # horizontal line (sp → left border)
    ax.hlines(sp, xmin, tp, colors="black", linestyles="--")

    # vertical line (tp → bottom)
    ax.vlines(tp, ymin, sp, colors="black", linestyles="--")
    ax.plot(tp, sp, "v", markersize=7, color="black")

    ax.set_title(well, fontsize=16)
    ax.set_xlabel("t (min)", fontsize=16)
    ax.set_ylabel("Drawdown (m)", fontsize=16)
    ax.grid(False)
    ax.legend()
plt.savefig("Drawdown & Recovery.png", dpi=300, bbox_inches="tight")
plt.tight_layout()
plt.show()
