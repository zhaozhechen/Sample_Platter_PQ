"""Audit source datasets used by the generalized Sampler Platter TE workflow.

This script is read-only. It streams large CSV files and writes compact metadata
summaries into 04_Results/Data_Audit.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT.parents[1] / "Data"
OUTPUT_ROOT = REPO_ROOT / "04_Results" / "Data_Audit"


def stream_group_counts(path: Path, columns: list[str], group_columns: list[str]) -> pd.DataFrame:
    totals: Counter[tuple] = Counter()
    complete: Counter[tuple] = Counter()
    minima: dict[tuple, pd.Timestamp] = {}
    maxima: dict[tuple, pd.Timestamp] = {}

    for chunk in pd.read_csv(path, usecols=columns, chunksize=250_000, low_memory=False):
        time_column = next(name for name in columns if "date" in name.lower())
        chunk[time_column] = pd.to_datetime(chunk[time_column], errors="coerce")
        value_columns = [name for name in columns if name not in group_columns + [time_column]]
        for key, group in chunk.groupby(group_columns, dropna=False, sort=False):
            if not isinstance(key, tuple):
                key = (key,)
            totals[key] += len(group)
            complete[key] += int(group[value_columns + [time_column]].notna().all(axis=1).sum())
            valid_time = group[time_column].dropna()
            if not valid_time.empty:
                group_min = valid_time.min()
                group_max = valid_time.max()
                minima[key] = min(minima.get(key, group_min), group_min)
                maxima[key] = max(maxima.get(key, group_max), group_max)

    records = []
    for key in sorted(totals, key=lambda item: tuple(str(value) for value in item)):
        record = {name: value for name, value in zip(group_columns, key)}
        record.update(
            rows=totals[key],
            complete_rows=complete[key],
            start=minima.get(key),
            end=maxima.get(key),
        )
        records.append(record)
    return pd.DataFrame(records)


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    hja = stream_group_counts(
        DATA_ROOT / "Hjandrews_Keith" / "hja_q_ppt_30min_processed.csv",
        ["sitecode_q", "sitecode_ppt", "datetime_pst", "q_30m_avg_m3s", "ppt_30m_total_mm"],
        ["sitecode_q", "sitecode_ppt"],
    )
    hja.to_csv(OUTPUT_ROOT / "HJA_site_pair_inventory.csv", index=False)

    casper_q = pd.read_csv(
        DATA_ROOT / "Casper_Salli" / "cas_daily_q_01to23.csv",
        usecols=["Date", "WS", "DailyTotQ_mm", "MeanDailyQ_cms"],
    )
    casper_p = pd.read_csv(
        DATA_ROOT / "Casper_Salli" / "cas_sf_daily_ppt_01to23.csv",
        usecols=["Year", "Month", "Day", "Site", "PPT_mm"],
    )
    casper_q["Date"] = pd.to_datetime(casper_q["Date"], errors="coerce")
    casper_p["Date"] = pd.to_datetime(casper_p[["Year", "Month", "Day"]], errors="coerce")

    q_summary = (
        casper_q.groupby("WS", dropna=False)
        .agg(rows=("Date", "size"), complete_rows=("MeanDailyQ_cms", "count"), start=("Date", "min"), end=("Date", "max"))
        .reset_index()
    )
    p_summary = (
        casper_p.groupby("Site", dropna=False)
        .agg(rows=("Date", "size"), complete_rows=("PPT_mm", "count"), start=("Date", "min"), end=("Date", "max"))
        .reset_index()
    )
    q_summary.to_csv(OUTPUT_ROOT / "Casper_discharge_site_inventory.csv", index=False)
    p_summary.to_csv(OUTPUT_ROOT / "Casper_precip_site_inventory.csv", index=False)

    print(f"Wrote data audit summaries to {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()

