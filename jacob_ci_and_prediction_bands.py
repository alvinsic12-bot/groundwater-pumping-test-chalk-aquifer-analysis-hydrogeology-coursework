import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, t


#SETTINGS

file = r"E:\course work\hydrogeology coursework\工作簿1.xlsx"
pumping_target = pd.to_datetime("2025-10-22 12:00:00")

Q = 58.0 / 1000.0   # pumping rate (m³/s)

# casing elevations (mAOD)
casing = {
    "PL10A": 107.358,
    "PL10B": 110.102,
    "PL10E": 107.0315
}

# distances to pumping well BBA (m)
r_dict = {
    "PL10A": 37.2,
    "PL10B": 53.5,
    "PL10E": 34.6
}

draw_color = "#65B2FA"
reco_color = "#FCC24C"


# Load Excel

df = pd.read_excel(
    file,
    skiprows=2,
    names=[
        "Time_PL10A", "PL10A_dip",
        "Time_PL10B", "PL10B_dip",
        "Time_PL10E", "PL10E_dip",
        "Time_BBA_dip", "BBA_dip",
        "col8",
        "Time_BBA_log", "BBA_log_WL",
        "Time_BBA_obs2", "BBA_obs2_WL"
    ]
)

def convert_time(series):
    return pd.to_datetime("2025-10-22 " + series.astype(str), errors="coerce")



# Function: 95% CI for slope & prediction band

def regression_CI(x, y, slope, intercept):
    n = len(x)
    y_pred = slope * x + intercept
    residuals = y - y_pred
    dof = n - 2

    tval = t.ppf(0.975, dof)
    s_err = np.sqrt(np.sum(residuals**2) / dof)
    x_mean = np.mean(x)
    Sxx = np.sum((x - x_mean)**2)

    # slope CI
    slope_se = s_err / np.sqrt(Sxx)
    slope_lo = slope - tval * slope_se
    slope_hi = slope + tval * slope_se

    # prediction band function
    def pred_band(x_new):
        return (
            slope * x_new + intercept,
            slope * x_new + intercept + tval * s_err * np.sqrt(1/n + ((x_new - x_mean)**2 / Sxx)),
            slope * x_new + intercept - tval * s_err * np.sqrt(1/n + ((x_new - x_mean)**2 / Sxx))
        )

    return slope_lo, slope_hi, pred_band



# Process wells

def process_well(time_col, dip_col, casing_elev):

    dfw = df[[time_col, dip_col]].dropna().copy()
    dfw["time_dt"] = convert_time(dfw[time_col])

    dfw["abs_diff"] = (dfw["time_dt"] - pumping_target).abs()
    idx0 = dfw["abs_diff"].idxmin()
    t0 = dfw.loc[idx0, "time_dt"]

    dfw = dfw[dfw["time_dt"] >= t0].copy()
    dfw["t_min"] = (dfw["time_dt"] - t0).dt.total_seconds() / 60

    dfw["WL"] = casing_elev - dfw[dip_col]
    static_WL = dfw["WL"].iloc[0]
    dfw["s"] = static_WL - dfw["WL"]
    dfw["ds"] = dfw["s"].diff().fillna(0)

    mask_flat = (dfw["s"].abs() < 0.002) & (dfw["ds"].abs() < 0.001)
    idx_start = dfw.index[~mask_flat][0]
    dfw = dfw.loc[idx_start:].copy()

    idx_max = dfw["s"].idxmax()
    dfD = dfw.loc[:idx_max].copy()
    dfR = dfw.loc[idx_max+1:].copy()

    dfD = dfD[dfD["t_min"] > 0]
    dfR = dfR[dfR["t_min"] > dfD["t_min"].max()]

    return dfD, dfR


draw_reco = {
    w: process_well(f"Time_{w}", f"{w}_dip", casing[w])
    for w in ["PL10A", "PL10B", "PL10E"]
}

# Compute Ts, S, CI, and prediction bands

results = {}

for well in ["PL10A", "PL10B", "PL10E"]:

    dfD, dfR = draw_reco[well]
    r = r_dict[well]

    # ---------- Drawdown ----------
    xD = np.log(dfD["t_min"])
    yD = dfD["s"]
    slopeD, interceptD, *_ = linregress(xD, yD)

    slopeD_lo, slopeD_hi, predD = regression_CI(xD.values, yD.values, slopeD, interceptD)

    T_draw = Q / (4 * np.pi * abs(slopeD))
    T_draw_lo = Q / (4 * np.pi * abs(slopeD_hi))
    T_draw_hi = Q / (4 * np.pi * abs(slopeD_lo))

    S_draw = 4 * T_draw / (r**2 * np.exp(interceptD/slopeD + 0.5772 ))

    # ---------- Recovery ----------
    tp = dfD["t_min"].iloc[-1]
    sp = dfD["s"].iloc[-1]

    tR = dfR["t_min"].values
    sR = dfR["s"].values

    mask = tR > tp
    tR = tR[mask]
    sR = sR[mask]

    sr = sR - sp
    xR = np.log(tR / (tR - tp))

    slopeR, interceptR, *_ = linregress(xR, sr)
    slopeR_lo, slopeR_hi, predR = regression_CI(xR, sr, slopeR, interceptR)

    T_reco = Q / (4 * np.pi * abs(slopeR))
    T_reco_lo = Q / (4 * np.pi * abs(slopeR_hi))
    T_reco_hi = Q / (4 * np.pi * abs(slopeR_lo))

    # Recovery S 
    S_reco = (4*T_reco*tp) / (r**2 * np.exp(0.5772 + (4*np.pi*T_reco/Q)*sp))

    results[well] = dict(
        xD=xD, yD=yD, predD=predD,
        xR=xR, yR=sr, predR=predR,
        slopeD=slopeD, interceptD=interceptD,
        slopeR=slopeR, interceptR=interceptR,
        T_draw=T_draw, T_draw_lo=T_draw_lo, T_draw_hi=T_draw_hi,
        T_reco=T_reco, T_reco_lo=T_reco_lo, T_reco_hi=T_reco_hi,
        S_draw=S_draw, S_reco=S_reco
    )



# Plot 1×3 panel WITH CONFIDENCE BANDS

fig, axes = plt.subplots(1, 3, figsize=(22, 7))
fig.suptitle("Jacobs Drawdown & Recovery (with 95% CI bands)", fontsize=17)

for ax, well in zip(axes, ["PL10A", "PL10B", "PL10E"]):

    res = results[well]

    # ---------------- Drawdown ----------------
    ax.plot(res["xD"], res["yD"], "o", markerfacecolor="none",
            markeredgecolor=draw_color, label="Drawdown")

    xd = np.linspace(res["xD"].min(), res["xD"].max(), 200)
    yfit, yhi, ylo = res["predD"](xd)
    ax.plot(xd, yfit, color=draw_color, lw=1.6)
    ax.fill_between(xd, ylo, yhi, color=draw_color, alpha=0.4)

    # ---------------- Recovery ----------------
    ax.plot(res["xR"], res["yR"], "s", markerfacecolor="none",
            markeredgecolor=reco_color, label="Recovery")

    xr = np.linspace(res["xR"].min(), res["xR"].max(), 200)
    yfit_r, yhi_r, ylo_r = res["predR"](xr)
    ax.plot(xr, yfit_r, color=reco_color, lw=1.6)
    ax.fill_between(xr, ylo_r, yhi_r, color=reco_color, alpha=0.4)

    # ---------------- Label ----------------
    ax.text(
        0.98, 0.02,
        f"T(d) = {res['T_draw']:.2e}\n"
        f"  CI = [{res['T_draw_lo']:.2e}, {res['T_draw_hi']:.2e}]\n"
        f"S(d) = {res['S_draw']:.2e}\n"
        f"T(r) = {res['T_reco']:.2e}\n"
        f"  CI = [{res['T_reco_lo']:.2e}, {res['T_reco_hi']:.2e}]\n"
        f"S(r) = {res['S_reco']:.2e}",
        transform=ax.transAxes,
        ha="right", va="bottom",
        fontsize=9,
        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none")
    )

    ax.set_title(well, fontsize=15)
    ax.set_xlabel("ln(t) or ln[t/(t−tp)]",fontsize=15)
    ax.set_ylabel("s or s_r",fontsize=15)
    ax.grid(False)
    ax.legend(loc="upper left")
plt.savefig("Jacobs Drawdown & Recovery.png", dpi=300, bbox_inches="tight")
plt.tight_layout(rect=[0,0,1,0.94])
plt.show()
