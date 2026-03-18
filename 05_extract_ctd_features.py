#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from collections import Counter
from typing import Dict, Optional

import pandas as pd
from Bio import SeqIO

try:
    import numpy as np
    import matplotlib.pyplot as plt
except Exception:
    plt = None

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
PROPERTIES = {
    "hydrophobicity": {"polar": "RKEDQN", "neutral": "GASTPHY", "hydrophobic": "CLVIMFW"},
    "vdw_volume": {"small": "GASTPDC", "medium": "NVEQIL", "large": "MHKFRYW"},
    "polarity": {"polar": "EDQNKR", "neutral": "GASTPHY", "nonpolar": "CLVIMFW"},
    "polarizability": {"low": "GASDT", "medium": "CPNVEQIL", "high": "MHKFRYW"},
    "charge": {"negative": "DE", "neutral": "ACFGHILMNPQSTVWY", "positive": "KR"},
    "secondary_structure": {"helix": "EALMQKRH", "sheet": "VIYCWT", "coil": "GNPSD"},
    "solvent_accessibility": {"buried": "ALFCGIVW", "intermediate": "MPSTHY", "exposed": "DEKNQR"},
}
PERCENTILES = [0.0, 0.25, 0.5, 0.75, 1.0]
PERCENTILE_LABELS = ["p0", "p25", "p50", "p75", "p100"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Extract CTD-like distribution features from peptide FASTA or dataset CSV.")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--fasta", help="Input FASTA file")
    group.add_argument("--dataset-csv", help="Input dataset CSV with sequence and optional label/id columns")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--label", default=None, help="Optional label filter when using --dataset-csv")
    ap.add_argument("--min-len", type=int, default=10)
    ap.add_argument("--max-len", type=int, default=200)
    ap.add_argument("--no-plots", action="store_true")
    return ap.parse_args()


def parse_label(description: str) -> Optional[str]:
    desc = description.strip()
    m = re.search(r"(?:^|\s)(?:label|Label)\s*=\s*([A-Za-z0-9\-_/]+)", desc)
    if m:
        return m.group(1)
    m = re.search(r"(?:^|\s)(?:label|Label)\s*:\s*([A-Za-z0-9\-_/]+)", desc)
    if m:
        return m.group(1)
    tokens = desc.split()
    if tokens and tokens[-1].upper() in {"AMP", "NON-AMP", "NONAMP", "NAMP"}:
        return tokens[-1]
    return None


def normalize_label(value: object) -> Optional[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip().replace("_", "-").replace(" ", "").lower()
    if s == "amp":
        return "AMP"
    if s in {"nonamp", "non-amp", "namp"}:
        return "non-AMP"
    return str(value)


def validate_sequence(seq: str) -> bool:
    s = str(seq).strip().upper().replace("*", "")
    return len(s) > 0 and set(s).issubset(VALID_AA)


def distribution_for_class(sequence: str, amino_class: str) -> Dict[str, float]:
    seq = sequence.upper()
    idxs = [i for i, aa in enumerate(seq) if aa in amino_class]
    total = len(idxs)
    if total == 0:
        return {lbl: float("nan") for lbl in PERCENTILE_LABELS}
    out: Dict[str, float] = {}
    for pct, lbl in zip(PERCENTILES, PERCENTILE_LABELS):
        pos_in_class = max(0, math.ceil(total * pct) - 1)
        pos_in_seq = idxs[pos_in_class]
        out[lbl] = (pos_in_seq + 1) / len(seq) * 100.0
    return out


def ctd_distribution_features(sequence: str) -> Dict[str, float]:
    features: Dict[str, float] = {}
    for prop, classes in PROPERTIES.items():
        for cls, residues in classes.items():
            dist = distribution_for_class(sequence, residues)
            for lbl, val in dist.items():
                features[f"{prop}_{cls}_{lbl}"] = val
    return features


def fasta_to_dataframe(fasta_path: str, min_len: int, max_len: int) -> pd.DataFrame:
    rows = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).strip().upper().replace("*", "")
        if not validate_sequence(seq):
            continue
        if not (min_len <= len(seq) <= max_len):
            continue
        rows.append({"ID": rec.id.split()[0], "Label": parse_label(rec.description), "Sequence": seq, "Length": len(seq)})
    return pd.DataFrame(rows)


def dataset_to_dataframe(csv_path: str, min_len: int, max_len: int, label_filter: Optional[str]) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    cols = {c.lower(): c for c in df.columns}
    if "sequence" not in cols:
        raise ValueError("dataset CSV must contain a 'sequence' column")
    seq_col = cols["sequence"]
    label_col = cols.get("label")
    id_col = cols.get("accession_version") or cols.get("accession") or cols.get("id")
    out = pd.DataFrame()
    out["Sequence"] = df[seq_col].astype(str).str.strip().str.upper().str.replace("*", "", regex=False)
    out = out[out["Sequence"].map(validate_sequence)].copy()
    out["Length"] = out["Sequence"].str.len()
    out = out[(out["Length"] >= min_len) & (out["Length"] <= max_len)].copy()
    if id_col:
        out["ID"] = df.loc[out.index, id_col].astype(str).values
    else:
        out["ID"] = [f"seq_{i+1:06d}" for i in range(len(out))]
    if label_col:
        out["Label"] = df.loc[out.index, label_col].map(normalize_label).values
    else:
        out["Label"] = None
    if label_filter is not None:
        wanted = normalize_label(label_filter)
        out = out[out["Label"] == wanted].copy()
    return out.reset_index(drop=True)


def ensure_plotting() -> None:
    if plt is None:
        raise RuntimeError("Matplotlib and NumPy are required for plotting")


def plot_length_hist(lengths: pd.Series, out_dir: str) -> str:
    ensure_plotting()
    os.makedirs(out_dir, exist_ok=True)
    bins = range(int(lengths.min()), int(lengths.max()) + 2)
    plt.figure(figsize=(6, 4))
    plt.hist(lengths, bins=bins)
    plt.title(f"Peptide length distribution (n={len(lengths)})")
    plt.xlabel("Length (aa)")
    plt.ylabel("Count")
    plt.tight_layout()
    path = os.path.join(out_dir, "length_hist.png")
    plt.savefig(path, dpi=150)
    plt.close()
    return path


def plot_aa_frequency(seqs: pd.Series, out_dir: str) -> str:
    ensure_plotting()
    os.makedirs(out_dir, exist_ok=True)
    aa_counter = Counter()
    for s in seqs:
        aa_counter.update(ch for ch in s if ch in VALID_AA)
    total = sum(aa_counter.values()) or 1
    aas = sorted(VALID_AA)
    freqs = [aa_counter[a] / total for a in aas]
    plt.figure(figsize=(7.5, 4))
    plt.bar(aas, freqs)
    plt.title("Amino acid frequency")
    plt.xlabel("AA")
    plt.ylabel("Fraction")
    plt.tight_layout()
    path = os.path.join(out_dir, "aa_frequency.png")
    plt.savefig(path, dpi=150)
    plt.close()
    return path


def plot_property_p50_hists(df_features: pd.DataFrame, out_dir: str) -> list[str]:
    ensure_plotting()
    os.makedirs(out_dir, exist_ok=True)
    saved = []
    for prop, classes in PROPERTIES.items():
        cls_list = list(classes.keys())
        fig, axes = plt.subplots(1, len(cls_list), figsize=(4 * len(cls_list), 3), squeeze=False)
        for i, cls in enumerate(cls_list):
            col = f"{prop}_{cls}_p50"
            series = pd.to_numeric(df_features.get(col), errors="coerce").dropna() if col in df_features.columns else pd.Series(dtype=float)
            ax = axes[0, i]
            if len(series) == 0:
                ax.text(0.5, 0.5, "No data", ha="center", va="center")
            else:
                ax.hist(series, bins=20, range=(0, 100))
            ax.set_title(f"{cls} (p50)")
            ax.set_xlabel("% position")
            ax.set_ylabel("Count")
        plt.suptitle(f"{prop} p50 distributions", y=1.02)
        plt.tight_layout()
        path = os.path.join(out_dir, f"{prop}_p50_hist.png")
        plt.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(path)
    return saved


def main() -> int:
    args = parse_args()
    if args.dataset_csv:
        df = dataset_to_dataframe(args.dataset_csv, args.min_len, args.max_len, args.label)
    else:
        df = fasta_to_dataframe(args.fasta, args.min_len, args.max_len)
    if df.empty:
        print("No valid sequences after filtering", file=sys.stderr)
        return 2
    feat = df["Sequence"].apply(ctd_distribution_features).apply(pd.Series)
    out_df = pd.concat([df, feat], axis=1)
    out_dir = os.path.dirname(args.out) or "."
    os.makedirs(out_dir, exist_ok=True)
    out_df.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")
    if not args.no_plots:
        plots_dir = os.path.join(out_dir, "plots")
        for path in [plot_length_hist(out_df["Length"], plots_dir), plot_aa_frequency(out_df["Sequence"], plots_dir), *plot_property_p50_hists(out_df, plots_dir)]:
            print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
