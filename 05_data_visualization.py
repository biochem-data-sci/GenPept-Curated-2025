#!/usr/bin/env python3
"""Reproduce the current manuscript length-stratified EDA from the frozen v1.1 dataset.

Source basis (no model/results code):
- DataVisualization_exact.ipynb: Wilson CI calculation and length-bin count/proportion plots.
- main (1).pdf: current Figure 2 is a three-panel length-bin composition figure.

The script reads the frozen canonical dataset and generates, at run time:
  - Figure2_Class_Composition_by_Length.png/pdf/svg
  - Table1_Length_Bin_Summary.csv

No published result table is bundled; the CSV is generated from the canonical dataset.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

BINS = [(10, 50), (51, 100), (101, 150), (151, 200)]
BIN_LABELS = ["10–50", "51–100", "101–150", "151–200"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce the manuscript length-stratified EDA from the frozen GenPept dataset."
    )
    parser.add_argument(
        "--dataset-csv",
        default="data/GenPept_Curated_2025_primary_split_v1.1.csv",
        help="Frozen canonical GenPept-Curated-2025 v1.1 CSV.",
    )
    parser.add_argument("--outdir", default="outputs/eda")
    return parser.parse_args()


def load_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    required = {"label", "sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if "length" not in df.columns:
        df["length"] = df["sequence"].astype(str).str.len()

    df["length"] = pd.to_numeric(df["length"], errors="raise").astype(int)
    labels = set(df["label"].astype(str).unique())
    if labels != {"AMP", "non-AMP"}:
        raise ValueError(f"Unexpected label set: {sorted(labels)}")
    return df


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (lo, hi), label in zip(BINS, BIN_LABELS):
        sub = df[(df["length"] >= lo) & (df["length"] <= hi)]
        total = int(len(sub))
        if total == 0:
            raise ValueError(f"No records in length bin {label}")
        amp = int((sub["label"] == "AMP").sum())
        non_amp = total - amp
        low, high = proportion_confint(amp, total, alpha=0.05, method="wilson")
        rows.append(
            {
                "length_interval_aa": label,
                "total": total,
                "AMP": amp,
                "non-AMP": non_amp,
                "AMP_proportion_percent": round(100.0 * amp / total, 1),
                "Wilson_95_CI_low_percent": round(100.0 * float(low), 1),
                "Wilson_95_CI_high_percent": round(100.0 * float(high), 1),
            }
        )

    out = pd.DataFrame(rows)
    total = int(len(df))
    amp = int((df["label"] == "AMP").sum())
    low, high = proportion_confint(amp, total, alpha=0.05, method="wilson")
    total_row = pd.DataFrame(
        [
            {
                "length_interval_aa": "Total",
                "total": total,
                "AMP": amp,
                "non-AMP": total - amp,
                "AMP_proportion_percent": round(100.0 * amp / total, 1),
                "Wilson_95_CI_low_percent": round(100.0 * float(low), 1),
                "Wilson_95_CI_high_percent": round(100.0 * float(high), 1),
            }
        ]
    )
    return pd.concat([out, total_row], ignore_index=True)


def plot_figure(summary: pd.DataFrame, outdir: Path) -> None:
    bins = summary.iloc[:4].copy()
    x = np.arange(len(bins))

    amp = bins["AMP"].to_numpy(dtype=float)
    non_amp = bins["non-AMP"].to_numpy(dtype=float)
    total = bins["total"].to_numpy(dtype=float)
    amp_pct = bins["AMP_proportion_percent"].to_numpy(dtype=float)
    non_pct = 100.0 - amp_pct
    ci_low = bins["Wilson_95_CI_low_percent"].to_numpy(dtype=float)
    ci_high = bins["Wilson_95_CI_high_percent"].to_numpy(dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    # Panel A: AMP proportion with Wilson 95% CI.
    yerr = np.vstack([amp_pct - ci_low, ci_high - amp_pct])
    axes[0].errorbar(x, amp_pct, yerr=yerr, fmt="o", capsize=4)
    axes[0].axhline(50.0, linestyle="--", linewidth=1.0)
    axes[0].set_xticks(x, BIN_LABELS)
    axes[0].set_ylim(0, 100)
    axes[0].set_ylabel("AMP proportion (%)")
    axes[0].set_title("A")
    for i, (n, k, p) in enumerate(zip(total.astype(int), amp.astype(int), amp_pct)):
        axes[0].annotate(f"n={n:,}; AMP={k:,}\n{p:.1f}%", (i, p), xytext=(0, 10), textcoords="offset points", ha="center", fontsize=8)

    # Panel B: absolute counts.
    width = 0.36
    axes[1].bar(x - width / 2, amp, width=width, label="AMP")
    axes[1].bar(x + width / 2, non_amp, width=width, label="non-AMP")
    axes[1].set_xticks(x, BIN_LABELS)
    axes[1].set_ylabel("Count")
    axes[1].set_title("B")
    axes[1].legend(frameon=False)

    # Panel C: normalized proportions.
    axes[2].bar(x, amp_pct, label="AMP")
    axes[2].bar(x, non_pct, bottom=amp_pct, label="non-AMP")
    axes[2].axhline(50.0, linestyle="--", linewidth=1.0)
    axes[2].set_xticks(x, BIN_LABELS)
    axes[2].set_ylim(0, 100)
    axes[2].set_ylabel("Proportion (%)")
    axes[2].set_title("C")
    axes[2].legend(frameon=False)

    for ax in axes:
        ax.set_xlabel("Length interval (aa)")

    fig.tight_layout()
    for suffix, kwargs in [
        ("png", {"dpi": 300}),
        ("pdf", {}),
        ("svg", {}),
    ]:
        fig.savefig(outdir / f"Figure2_Class_Composition_by_Length.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    dataset = Path(args.dataset_csv)
    outdir = Path(args.outdir)
    if not dataset.is_file():
        raise FileNotFoundError(dataset)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_dataset(dataset)
    summary = build_summary(df)

    expected = {
        "10–50": (2618, 372, 2246),
        "51–100": (2379, 1602, 777),
        "101–150": (2281, 1759, 522),
        "151–200": (3722, 1767, 1955),
    }
    for _, row in summary.iloc[:4].iterrows():
        observed = (int(row["total"]), int(row["AMP"]), int(row["non-AMP"]))
        if observed != expected[str(row["length_interval_aa"])]:
            raise RuntimeError(
                f"Frozen manuscript count mismatch for {row['length_interval_aa']}: "
                f"expected {expected[str(row['length_interval_aa'])]}, observed {observed}"
            )

    summary.to_csv(outdir / "Table1_Length_Bin_Summary.csv", index=False)
    plot_figure(summary, outdir)
    print(summary.to_string(index=False))
    print(f"Wrote manuscript EDA outputs to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
