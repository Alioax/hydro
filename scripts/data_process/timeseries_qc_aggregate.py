# -*- coding: utf-8 -*-
"""
Time series QC and aggregation for CAMELS-FR daily data.

- Load daily P, Q, PET (Oudin), temperature; skip comment lines.
- Select analysis period (configurable, default 1990-2020 for water balance).
- Missing-data policy: drop days with missing or invalid Q; aggregate to monthly
  only for months with at least 90% valid days.
- Output: monthly and annual CSV (date, P_mm, Q_mm, PET_mm, T_C) for water balance and plots.
"""
import csv
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "camels_required"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed_data"
STATION_CODE = "K287191001"
# Analysis period for water balance and FDC
START_YEAR = 1990
END_YEAR = 2020
# Minimum fraction of valid days in a month to keep that month
MIN_DAYS_FRAC = 0.9


def load_daily(station_code):
    path = DATA_DIR / "CAMELS_FR_time_series" / "daily" / f"CAMELS_FR_tsd_{station_code}.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    lines = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip().startswith("#"):
                continue
            lines.append(line)
    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep=";")
    df["date"] = pd.to_datetime(df["tsd_date"].astype(str), format="%Y%m%d")
    df = df.rename(columns={
        "tsd_prec": "P_mm",
        "tsd_q_mm": "Q_mm",
        "tsd_pet_ou": "PET_mm",
        "tsd_temp": "T_C",
    })
    return df[["date", "P_mm", "Q_mm", "PET_mm", "T_C"]].copy()


def run(start_year=START_YEAR, end_year=END_YEAR, min_days_frac=MIN_DAYS_FRAC):
    df = load_daily(STATION_CODE)
    # Replace comma by dot for numeric columns
    for col in ["P_mm", "Q_mm", "PET_mm", "T_C"]:
        if df[col].dtype == object:
            df[col] = pd.to_numeric(df[col].astype(str).str.replace(",", "."), errors="coerce")
    df = df.sort_values("date").reset_index(drop=True)

    # Restrict to analysis period
    t0 = f"{start_year}-01-01"
    t1 = f"{end_year}-12-31"
    df = df[(df["date"] >= t0) & (df["date"] <= t1)].copy()

    # Missing-data: drop days with missing or invalid Q (negative or NaN)
    valid_q = df["Q_mm"].notna() & (df["Q_mm"] >= 0)
    df = df[valid_q].copy()

    # Monthly: only months with >= min_days_frac of valid days (we already dropped invalid Q days)
    df["year"] = df["date"].dt.year
    df["month"] = df["date"].dt.month
    days_per_month = df.groupby(["year", "month"]).size()
    expected_days = pd.Series(
        index=days_per_month.index,
        data=[pd.Period(f"{y}-{m:02d}").days_in_month for y, m in days_per_month.index],
    )
    keep_month = (days_per_month >= (min_days_frac * expected_days)).to_dict()
    df["keep_month"] = df.apply(lambda r: keep_month.get((r["year"], r["month"]), False), axis=1)
    df_month = df[df["keep_month"]].groupby(["year", "month"]).agg(
        P_mm=("P_mm", "sum"),
        Q_mm=("Q_mm", "sum"),
        PET_mm=("PET_mm", "sum"),
        T_C=("T_C", "mean"),
    ).reset_index()
    df_month["date"] = pd.to_datetime(df_month["year"].astype(str) + "-" + df_month["month"].astype(str) + "-01")
    df_month = df_month[["date", "P_mm", "Q_mm", "PET_mm", "T_C"]]

    # Annual
    df_annual = df_month.groupby(df_month["date"].dt.year).agg(
        P_mm=("P_mm", "sum"),
        Q_mm=("Q_mm", "sum"),
        PET_mm=("PET_mm", "sum"),
        T_C=("T_C", "mean"),
    ).reset_index()
    df_annual = df_annual.rename(columns={"date": "year"})

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    out_month = PROCESSED_DIR / f"monthly_{STATION_CODE}_{start_year}_{end_year}.csv"
    out_annual = PROCESSED_DIR / f"annual_{STATION_CODE}_{start_year}_{end_year}.csv"
    df_month.to_csv(out_month, index=False)
    df_annual.to_csv(out_annual, index=False)

    # Also save daily (filtered) for FDC and time-series snippet
    out_daily = PROCESSED_DIR / f"daily_{STATION_CODE}_{start_year}_{end_year}.csv"
    df.to_csv(out_daily, index=False)

    print(f"Monthly: {len(df_month)} rows -> {out_month}")
    print(f"Annual: {len(df_annual)} rows -> {out_annual}")
    print(f"Daily: {len(df)} rows -> {out_daily}")
    return df_month, df_annual, df


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--start-year", type=int, default=START_YEAR)
    p.add_argument("--end-year", type=int, default=END_YEAR)
    args = p.parse_args()
    run(start_year=args.start_year, end_year=args.end_year)
