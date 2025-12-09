import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# INPUT

file = r"E:\course work\hydrogeology coursework\工作簿1.xlsx"

casing = {"PL10A": 107.358, "PL10B": 110.102, "PL10E": 107.0315}
r_dict = {"PL10A": 37.2, "PL10B": 53.5, "PL10E": 34.6}

Q = 58/1000         # m³/s
pumping_target = pd.to_datetime("2025-10-22 12:00:00")

draw_color = "#65B2FA"
reco_color = "#FCC24C"



# Load data

df = pd.read_excel(
    file,
    skiprows=2,
    names=[
        "Time_PL10A","PL10A_dip",
        "Time_PL10B","PL10B_dip",
        "Time_PL10E","PL10E_dip",
        "Time_BBA_dip","BBA_dip",
        "col8",
        "Time_BBA_log","BBA_log_WL",
        "Time_BBA_obs2","BBA_obs2_WL"
    ]
)

def convert(series):
    return pd.to_datetime("2025-10-22 " + series.astype(str), errors="coerce")



# Process wells (baseline = nearest to 12:00)

def process_well(time_col, dip_col, casing_elev):

    dfw = df[[time_col, dip_col]].dropna().copy()
    dfw["time_dt"] = convert(dfw[time_col])

    # find baseline
    dfw["abs_diff"] = (dfw["time_dt"] - pumping_target).abs()
    idx_base = dfw["abs_diff"].idxmin()
    t0 = dfw.loc[idx_base, "time_dt"]

    dfw = dfw[dfw["time_dt"] >= t0].copy()
    dfw["t_min"] = (dfw["time_dt"] - t0).dt.total_seconds()/60

    # WL → drawdown
    dfw["WL"] = casing_elev - dfw[dip_col]
    static = dfw["WL"].iloc[0]
    dfw["s"] = static - dfw["WL"]

    # remove early no-response region
    dfw["ds"] = dfw["s"].diff().fillna(0)
    mask = (dfw["s"].abs()<0.002)&(dfw["ds"].abs()<0.001)
    idx = dfw.index[~mask][0]
    dfw = dfw.loc[idx:].copy()

    # split at max s
    idx_max = dfw["s"].idxmax()
    dfD = dfw.loc[:idx_max].copy()
    dfR = dfw.loc[idx_max+1:].copy()

    tp = dfD["t_min"].iloc[-1]

    return dfD, dfR, tp



# Jacob functions

def jacob_draw(t, T, S, r):
    return (Q/(4*np.pi*T)) * ( np.log(4*T*t/(r*r*S)) - 0.5772 )

def jacob_recovery(t, T, tp):
    return (Q/(4*np.pi*T)) * np.log(t/(t-tp))



# Fit each well

results = {}
wells = ["PL10A","PL10B","PL10E"]

for w in wells:

    dfD, dfR, tp = process_well(f"Time_{w}", f"{w}_dip", casing[w])
    r = r_dict[w]

    tD = dfD["t_min"].values
    sD = dfD["s"].values

    tR = dfR["t_min"].values
    sR = dfR["s"].values

    # --- Drawdown fit ---
    poptD, _ = curve_fit(lambda t,T,S: jacob_draw(t,T,S,r),
                         tD, sD, p0=[0.002, 0.01], maxfev=20000)
    T_draw, S_draw = poptD

    # --- Recovery fit ---
    poptR, _ = curve_fit(lambda t,T: jacob_recovery(t,T,tp),
                         tR, sR, p0=[T_draw], maxfev=20000)
    T_reco = poptR[0]

    # --- Both fit ---
    t_all = np.concatenate([tD, tR])
    s_all = np.concatenate([sD, sR])

    def both_model(t, T, S):
        out = np.zeros_like(t)
        maskD = t <= tp
        maskR = t > tp
        out[maskD] = jacob_draw(t[maskD], T, S, r)
        out[maskR] = jacob_recovery(t[maskR], T, tp)
        return out

    poptB, _ = curve_fit(both_model, t_all, s_all,
                         p0=[T_draw, S_draw], maxfev=20000)
    T_both, S_both = poptB

    results[w] = dict(T_draw=T_draw, S_draw=S_draw,
                      T_reco=T_reco,
                      T_both=T_both, S_both=S_both,
                      tD=tD, sD=sD, tR=tR, sR=sR, tp=tp, r=r)



# PLOT 1 × 3 PANEL

fig, axes = plt.subplots(1,3, figsize=(22,7))
fig.suptitle("Simultaneous Jacob Fitting (Drawdown / Recovery / Both)", fontsize=20)

for ax, w in zip(axes, wells):

    res = results[w]

    # ---------------- data point----------------
    ax.scatter(res["tD"], res["sD"], edgecolor=draw_color,
               facecolor="none", label="Drawdown data")
    ax.scatter(res["tR"], res["sR"], edgecolor=reco_color,
               facecolor="none", marker="s", label="Recovery data")

    # ---------------- Drawdown fit ----------------
    tfitD = np.linspace(res["tD"].min(), res["tD"].max(), 200)
    ax.plot(tfitD, jacob_draw(tfitD, res["T_draw"], res["S_draw"], res["r"]),
            draw_color, lw=2.5, alpha=1, label="Drawdown fit")

    # ---------------- Recovery fit ----------------
    tfitR = np.linspace(res["tR"].min(), res["tR"].max(), 200)
    ax.plot(tfitR, jacob_recovery(tfitR, res["T_reco"], res["tp"]),
            reco_color, lw=2.5, alpha=1, label="Recovery fit")

    
    tfit_left = np.linspace(res["tD"].min(), res["tp"], 200)
    tfit_right = np.linspace(res["tp"], res["tR"].max(), 200)

    y_left = jacob_draw(tfit_left, res["T_both"], res["S_both"], res["r"])
    y_right = jacob_recovery(tfit_right, res["T_both"], res["tp"])

    ax.plot(tfit_left, y_left, "k--", lw=1.5, alpha=1)
    ax.plot(tfit_right, y_right, "k--", lw=1.5, alpha=1, label="Both fit")


    sp = res["sD"][-1]
    ax.scatter(res["tp"], sp, s=80, marker="v", color="k", zorder=5)

    ax.set_title(w, fontsize=15)
    ax.set_xlabel("Time after pumping started (min)", fontsize=15)
    ax.set_ylabel("Drawdown (m)", fontsize=15)
    ax.grid(False)
    ax.legend()
plt.savefig("Simultaneous Jacob Fitting.png", dpi=300, bbox_inches="tight")
plt.tight_layout(rect=[0,0,1,0.92])
plt.show()



# PRINT RESULTS TABLE

print("\n================ PARAMETERS TABLE ================\n")
for w in wells:
    r=results[w]
    print(f"--- {w} ---")
    print(f"T_draw = {r['T_draw']:.4e}   S_draw = {r['S_draw']:.4e}")
    print(f"T_both = {r['T_both']:.4e}   S_both = {r['S_both']:.4e}")
    print(f"T_reco = {r['T_reco']:.4e}\n")
