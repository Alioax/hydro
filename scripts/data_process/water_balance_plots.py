# -*- coding: utf-8 -*-
"""
Water balance and hydroclimate figures for the report.

- Time series snippet: flow rate, rainfall, PET, temperature (2017-2020), stacked subplots, shared time axis.
- Monthly climatology: 12-month average P, Q, PET.
- Flow duration curve: % time exceeded vs discharge.
- Annual water balance: bar chart P, Q, PET, residual by year.

All outputs to reports/report 1/figs/
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["figure.dpi"] = 800
plt.rcParams["font.family"] = "Times New Roman"
mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
    color=["#FF5F05", "#13294B", "#009FD4", "#FCB316",
           "#006230", "#007E8E", "#5C0E41", "#7D3E13"]
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed_data"
FIG_DIR = PROJECT_ROOT / "reports" / "report 1" / "figs"
STATION_CODE = "K287191001"
START_YEAR = 1990
END_YEAR = 2020
SNIPPET_START = "2017-01-01"
SNIPPET_END = "2020-12-31"


def load_daily():
    path = PROCESSED_DIR / f"daily_{STATION_CODE}_{START_YEAR}_{END_YEAR}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Run timeseries_qc_aggregate.py first: {path}")
    df = pd.read_csv(path)
    df["date"] = pd.to_datetime(df["date"])
    return df


def load_monthly():
    path = PROCESSED_DIR / f"monthly_{STATION_CODE}_{START_YEAR}_{END_YEAR}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Run timeseries_qc_aggregate.py first: {path}")
    df = pd.read_csv(path)
    df["date"] = pd.to_datetime(df["date"])
    return df


def load_annual():
    path = PROCESSED_DIR / f"annual_{STATION_CODE}_{START_YEAR}_{END_YEAR}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Run timeseries_qc_aggregate.py first: {path}")
    return pd.read_csv(path)


def run():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    daily = load_daily()
    monthly = load_monthly()
    annual = load_annual()

    # --- 1) Time series snippet: flow, rainfall, PET, temperature (2017-2020), aligned ---
    snippet = daily[(daily["date"] >= SNIPPET_START) & (daily["date"] <= SNIPPET_END)].copy()
    if len(snippet) == 0:
        snippet = daily.tail(365 * 3)  # fallback: last 3 years
    fig, axes = plt.subplots(4, 1, figsize=(9.5, 7.25), sharex=True)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    # Discharge Q: navy (C1); Precipitation P: blue (C2); PET: orange (C0); Temperature: orange (C0)
    axes[0].fill_between(snippet["date"], snippet["Q_mm"], alpha=0.6, color="C1")
    axes[0].set_ylabel("Discharge (mm/day)")
    axes[0].set_title("Flow rate")
    axes[0].grid(True, alpha=0.3)
    axes[1].bar(snippet["date"], snippet["P_mm"], width=1, color="C2", edgecolor="none")
    axes[1].set_ylabel("Precipitation (mm)")
    axes[1].set_title("Rainfall")
    axes[1].grid(True, alpha=0.3)
    axes[2].plot(snippet["date"], snippet["PET_mm"], color="C0", linewidth=0.8)
    axes[2].set_ylabel("PET (mm/day)")
    axes[2].set_title("PET")
    axes[2].grid(True, alpha=0.3)
    axes[3].plot(snippet["date"], snippet["T_C"], color="C0", linewidth=0.8)
    axes[3].set_ylabel("Temperature (°C)")
    axes[3].set_xlabel("Date")
    axes[3].set_title("Temperature")
    axes[3].grid(True, alpha=0.3)
    # Minimal symmetric x-margin so plot uses full width (shared x across subplots)
    for ax in axes:
        ax.margins(x=0.005)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "timeseries_flow_precip_temp_2017_2020.pdf", bbox_inches="tight", dpi=800)
    plt.close()

    # --- 2) Monthly climatology (12-month average) ---
    monthly["month"] = monthly["date"].dt.month
    clim = monthly.groupby("month").agg(
        P_mm=("P_mm", "mean"),
        Q_mm=("Q_mm", "mean"),
        PET_mm=("PET_mm", "mean"),
    ).reset_index()
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 2.5))
    ax.spines[["top", "right"]].set_visible(False)
    x = np.arange(12)
    w = 0.25
    # P: blue (C2), Q: navy (C1), PET: orange (C0)
    ax.bar(x - w, clim["P_mm"], width=w, label="P (mm)", color="C2")
    ax.bar(x, clim["Q_mm"], width=w, label="Q (mm)", color="C1")
    ax.bar(x + w, clim["PET_mm"], width=w, label="PET (mm)", color="C0", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(months)
    ax.set_ylabel("Monthly mean (mm)")
    fig.legend(loc="upper center", frameon=False, ncol=3, bbox_to_anchor=(0.5, 1.075))
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "monthly_climatology.pdf", bbox_inches="tight", dpi=800)
    plt.close()

    # --- 3) Flow duration curve ---
    q = daily["Q_mm"].dropna()
    q = q[q >= 0]
    q_sorted = np.sort(q)[::-1]  # descending
    n = len(q_sorted)
    pct = 100 * np.arange(1, n + 1) / (n + 1)  # % time flow is exceeded
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.75))
    ax.spines[["top", "right"]].set_visible(False)
    ax.semilogy(pct, q_sorted, color="C1", linewidth=1)  # Q: navy, same as discharge elsewhere
    ax.set_xlabel("Percent of time flow is exceeded")
    ax.set_ylabel("Discharge (mm/day)")
    ax.grid(True, alpha=0.3, which="both")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "flow_duration_curve.pdf", bbox_inches="tight", dpi=800)
    plt.close()

    # --- 4) Annual water balance ---
    ann = annual.copy()
    ann["residual"] = ann["P_mm"] - ann["Q_mm"] - ann["PET_mm"]
    x = np.arange(len(ann))
    w = 0.2
    fig, ax = plt.subplots(1, 1, figsize=(9.5, 4))
    ax.spines[["top", "right"]].set_visible(False)
    # P: blue (C2), Q: navy (C1), PET: orange (C0), Residual: green (C4)
    ax.bar(x - 1.5 * w, ann["P_mm"], width=w, label="P", color="C2")
    ax.bar(x - 0.5 * w, ann["Q_mm"], width=w, label="Q", color="C1")
    ax.bar(x + 0.5 * w, ann["PET_mm"], width=w, label="PET", color="C0", alpha=0.8)
    ax.bar(x + 1.5 * w, ann["residual"], width=w, label="Residual (P−Q−PET)", color="C4", alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(ann["year"].astype(int), rotation=45, ha="right")
    ax.set_ylabel("mm/yr")
    ax.margins(x=0.025)
    fig.legend(loc="upper center", frameon=False, ncol=4, bbox_to_anchor=(0.5, 1.05))
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "annual_water_balance.pdf", bbox_inches="tight", dpi=800)
    plt.close()

    # Summary numbers (for report)
    p_mean = ann["P_mm"].mean()
    q_mean = ann["Q_mm"].mean()
    pet_mean = ann["PET_mm"].mean()
    res_mean = ann["residual"].mean()
    runoff_ratio = q_mean / p_mean if p_mean else np.nan
    res_frac = (p_mean - q_mean - pet_mean) / p_mean if p_mean else np.nan
    summary_path = PROCESSED_DIR / "water_balance_summary.txt"
    with open(summary_path, "w") as f:
        f.write(f"Mean annual P (mm/yr): {p_mean:.1f}\n")
        f.write(f"Mean annual Q (mm/yr): {q_mean:.1f}\n")
        f.write(f"Mean annual PET (mm/yr): {pet_mean:.1f}\n")
        f.write(f"Mean residual (mm/yr): {res_mean:.1f}\n")
        f.write(f"Runoff ratio (Q/P): {runoff_ratio:.3f}\n")
        f.write(f"Residual fraction ((P-Q-PET)/P): {res_frac:.3f}\n")
    print(f"Water balance summary -> {summary_path}")
    print(f"Figures saved to {FIG_DIR}")
    return None


if __name__ == "__main__":
    run()
