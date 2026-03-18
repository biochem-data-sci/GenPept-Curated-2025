#!/usr/bin/env python3
from __future__ import annotations

import argparse
import platform
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BINS = [10, 50, 100, 150, 200]
BIN_LABELS = ["10-50", "51-100", "101-150", "151-200"]
TARGET_SPLIT = {"train": 7700, "val": 990, "test": 2310}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Run CD-HIT and create cluster-intact train/val/test splits.")
    ap.add_argument("--dataset-csv", required=True, help="dataset.csv from step 02")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cdhit-bin", default="cd-hit", help="CD-HIT executable name")
    ap.add_argument("--cdhit-id", type=float, default=0.90)
    ap.add_argument("--cdhit-word", type=int, default=5)
    ap.add_argument("--cdhit-cov-short", type=float, default=0.80, help="Coverage threshold on the shorter sequence")
    ap.add_argument("--seed", type=int, default=42)
    return ap.parse_args()


def normalize_label(x: str) -> str:
    s = str(x).strip().replace("_", "-").replace(" ", "").lower()
    if s == "amp":
        return "AMP"
    if s in {"non-amp", "nonamp", "namp"}:
        return "non-AMP"
    return str(x)


def load_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    req = {"sequence", "label"}
    missing = req - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in dataset: {sorted(missing)}")
    if "accession_version" not in df.columns:
        df["accession_version"] = [f"seq_{i+1:06d}" for i in range(len(df))]
    if "length" not in df.columns:
        df["length"] = df["sequence"].astype(str).str.len()
    if "length_bin" not in df.columns:
        df["length_bin"] = pd.cut(df["length"], bins=BINS, labels=BIN_LABELS, include_lowest=True)
    df["label"] = df["label"].map(normalize_label)
    expected_total = sum(TARGET_SPLIT.values())
    if len(df) != expected_total:
        raise ValueError(f"Expected {expected_total} sequences for the official split workflow, found {len(df)}")
    return df.reset_index(drop=True)


def write_cdhit_input_fasta(df: pd.DataFrame, path: Path) -> None:
    records = []
    for _, row in df.iterrows():
        header = f"{row['accession_version']}|{row['label']}|len:{row['length']}"
        records.append(SeqRecord(Seq(str(row["sequence"])), id=header, description=""))
    SeqIO.write(records, str(path), "fasta")


def resolve_wsl_path(path: Path) -> str:
    out = subprocess.run(["wsl", "wslpath", "-a", str(path.resolve())], capture_output=True, text=True, check=True)
    return out.stdout.strip()


def run_cdhit(fasta_in: Path, out_prefix: Path, executable: str, cdhit_id: float, cdhit_word: int, cdhit_cov_short: float) -> None:
    native = shutil.which(executable)
    cmd_native = [
        executable,
        "-i",
        str(fasta_in),
        "-o",
        str(out_prefix),
        "-c",
        str(cdhit_id),
        "-n",
        str(cdhit_word),
        "-aS",
        str(cdhit_cov_short),
        "-G",
        "0",
        "-g",
        "1",
        "-d",
        "0",
    ]
    if native:
        subprocess.run(cmd_native, check=True)
        return
    if platform.system().lower().startswith("win"):
        cmd_wsl = [
            "wsl",
            executable,
            "-i",
            resolve_wsl_path(fasta_in),
            "-o",
            resolve_wsl_path(out_prefix),
            "-c",
            str(cdhit_id),
            "-n",
            str(cdhit_word),
            "-aS",
            str(cdhit_cov_short),
            "-G",
            "0",
            "-g",
            "1",
            "-d",
            "0",
        ]
        subprocess.run(cmd_wsl, check=True)
        return
    raise RuntimeError("CD-HIT not found in PATH")


def parse_cdhit_clusters(clstr_path: Path) -> list[list[str]]:
    clusters: list[list[str]] = []
    current: list[str] = []
    with clstr_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">Cluster"):
                if current:
                    clusters.append(current)
                current = []
                continue
            m = re.search(r">([^|>\s]+)\|", line)
            if m:
                current.append(m.group(1))
    if current:
        clusters.append(current)
    return clusters


def build_component_table(df: pd.DataFrame, acc2cluster: dict[str, int]) -> pd.DataFrame:
    df2 = df.copy()
    df2["cluster"] = df2["accession_version"].astype(str).map(acc2cluster)
    if df2["cluster"].isna().any():
        missing = df2[df2["cluster"].isna()]["accession_version"].head(10).tolist()
        raise ValueError(f"Some accessions were not mapped to CD-HIT clusters, e.g. {missing}")
    comp = (
        df2.groupby(["cluster", "label", "length_bin"]).size().unstack(["label", "length_bin"]).fillna(0).astype(int)
    )
    full_cols = pd.MultiIndex.from_product([["AMP", "non-AMP"], BIN_LABELS])
    comp = comp.reindex(columns=full_cols, fill_value=0)
    return comp


def infer_targets(df: pd.DataFrame) -> dict[str, dict[tuple[str, str], int]]:
    frac = {k: v / len(df) for k, v in TARGET_SPLIT.items()}
    bin_tot = df["length_bin"].value_counts().to_dict()
    split_bin_target: dict[str, dict[str, int]] = {s: {} for s in TARGET_SPLIT}
    for b in BIN_LABELS:
        raw = {s: bin_tot.get(b, 0) * frac[s] for s in TARGET_SPLIT}
        alloc = {s: int(2 * round(raw[s] / 2)) for s in TARGET_SPLIT}
        diff = int(bin_tot.get(b, 0) - sum(alloc.values()))
        while diff != 0:
            moved = False
            for s in sorted(TARGET_SPLIT, key=lambda x: raw[x] - alloc[x], reverse=(diff > 0)):
                if diff > 0:
                    alloc[s] += 2
                    diff -= 2
                    moved = True
                elif alloc[s] >= 2:
                    alloc[s] -= 2
                    diff += 2
                    moved = True
                if diff == 0:
                    break
            if not moved:
                break
        for s in TARGET_SPLIT:
            split_bin_target[s][b] = max(0, alloc[s])
    needs: dict[str, dict[tuple[str, str], int]] = {s: {} for s in TARGET_SPLIT}
    for s in TARGET_SPLIT:
        for b in BIN_LABELS:
            per_label = split_bin_target[s][b] // 2
            needs[s][("AMP", b)] = per_label
            needs[s][("non-AMP", b)] = per_label
    return needs


def score_put(split: str, row: pd.Series, needs: dict[str, dict[tuple[str, str], int]]) -> tuple[int, int]:
    overflow_penalty = 0
    fit_score = 0
    for lab in ["AMP", "non-AMP"]:
        for b in BIN_LABELS:
            r = int(row.get((lab, b), 0))
            deficit = max(0, needs[split][(lab, b)])
            fit_score += min(r, deficit)
            overflow_penalty += max(0, r - deficit)
    split_current = sum(TARGET_SPLIT[split] for _ in [0]) - (sum(needs[split].values()))
    size_penalty = max(0, split_current + int(row.sum()) - TARGET_SPLIT[split])
    return (overflow_penalty * 100 + size_penalty * 10, -fit_score)


def assign_clusters(comp: pd.DataFrame, needs: dict[str, dict[tuple[str, str], int]], seed: int) -> dict[int, str]:
    assign: dict[int, str] = {}
    cluster_order = list(comp.index)
    rng = random.Random(seed)
    rng.shuffle(cluster_order)
    cluster_order.sort(key=lambda cid: tuple(-int(comp.loc[cid].get((lab, b), 0)) for lab in ["AMP", "non-AMP"] for b in BIN_LABELS))
    for cluster_id in cluster_order:
        row = comp.loc[cluster_id]
        best = min(TARGET_SPLIT, key=lambda s: score_put(s, row, needs))
        assign[int(cluster_id)] = best
        for lab in ["AMP", "non-AMP"]:
            for b in BIN_LABELS:
                r = int(row.get((lab, b), 0))
                if r > 0:
                    needs[best][(lab, b)] = max(0, needs[best][(lab, b)] - r)
    return assign


def repair_with_singletons(df: pd.DataFrame, comp: pd.DataFrame, assign: dict[int, str], targets: dict[str, dict[tuple[str, str], int]]) -> dict[int, str]:
    # Most clusters are singletons; move singleton clusters to fix small off-by-one mismatches.
    def counts(assign_map: dict[int, str]) -> dict[str, dict[tuple[str, str], int]]:
        out = {s: {(lab, b): 0 for lab in ["AMP", "non-AMP"] for b in BIN_LABELS} for s in TARGET_SPLIT}
        for cluster_id, split in assign_map.items():
            row = comp.loc[cluster_id]
            for lab in ["AMP", "non-AMP"]:
                for b in BIN_LABELS:
                    out[split][(lab, b)] += int(row.get((lab, b), 0))
        return out

    singleton_clusters = {cid for cid, row in comp.iterrows() if int(row.sum()) == 1}
    current = counts(assign)
    for lab in ["AMP", "non-AMP"]:
        for b in BIN_LABELS:
            for need_split in TARGET_SPLIT:
                target = targets[need_split][(lab, b)]
                while current[need_split][(lab, b)] < target:
                    donor = next((s for s in TARGET_SPLIT if current[s][(lab, b)] > targets[s][(lab, b)]), None)
                    if donor is None:
                        break
                    moved = False
                    for cluster_id in singleton_clusters:
                        if assign.get(int(cluster_id)) != donor:
                            continue
                        row = comp.loc[int(cluster_id)]
                        if int(row.get((lab, b), 0)) != 1:
                            continue
                        assign[int(cluster_id)] = need_split
                        current[donor][(lab, b)] -= 1
                        current[need_split][(lab, b)] += 1
                        moved = True
                        break
                    if not moved:
                        break
    return assign


def validate_split_output(df: pd.DataFrame) -> None:
    observed = df["split"].value_counts().to_dict()
    for split, expected in TARGET_SPLIT.items():
        actual = int(observed.get(split, 0))
        if actual != expected:
            raise ValueError(f"Split {split} has {actual} sequences; expected {expected}.")

    expected_per_label = {
        "train": {"AMP": 3850, "non-AMP": 3850},
        "val": {"AMP": 495, "non-AMP": 495},
        "test": {"AMP": 1155, "non-AMP": 1155},
    }
    observed_labels = df.groupby(["split", "label"]).size().to_dict()
    for split, label_targets in expected_per_label.items():
        for label, expected in label_targets.items():
            actual = int(observed_labels.get((split, label), 0))
            if actual != expected:
                raise ValueError(
                    f"Split {split} / label {label} has {actual} sequences; expected {expected}."
                )


def export_split_files(df: pd.DataFrame, outdir: Path) -> None:
    release_csv = outdir / "release_cluster_split.csv"
    df.to_csv(release_csv, index=False)
    for split in ["train", "val", "test"]:
        sub = df[df["split"] == split].copy()
        sub.to_csv(outdir / f"{split}.csv", index=False)
        records = [
            SeqRecord(Seq(str(r.sequence)), id=str(r.accession_version), description=f"label={r.label}")
            for r in sub.itertuples(index=False)
        ]
        SeqIO.write(records, str(outdir / f"{split}.fasta"), "fasta")
    stats = df.groupby(["split", "label", "length_bin"]).size().rename("count").reset_index()
    stats.to_csv(outdir / "split_summary_by_bin.csv", index=False)
    cluster_map = df[["accession_version", "cluster", "split"]].copy()
    cluster_map.to_csv(outdir / "cluster_assignment.csv", index=False)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_dataset(Path(args.dataset_csv))
    fasta_in = outdir / "all_11000_for_cdhit.fasta"
    cdhit_prefix = outdir / f"cdhit_c{args.cdhit_id}_n{args.cdhit_word}"
    clstr_path = outdir / f"cdhit_c{args.cdhit_id}_n{args.cdhit_word}.clstr"

    write_cdhit_input_fasta(df, fasta_in)
    run_cdhit(fasta_in, cdhit_prefix, args.cdhit_bin, args.cdhit_id, args.cdhit_word, args.cdhit_cov_short)
    if not clstr_path.exists():
        print(f"ERROR: missing CD-HIT cluster file: {clstr_path}", file=sys.stderr)
        return 2

    clusters = parse_cdhit_clusters(clstr_path)
    acc2cluster = {acc: i for i, members in enumerate(clusters) for acc in members}
    comp = build_component_table(df, acc2cluster)
    targets = infer_targets(df)
    assign = assign_clusters(comp, needs={s: dict(v) for s, v in targets.items()}, seed=args.seed)
    assign = repair_with_singletons(df, comp, assign, targets)

    out = df.copy()
    out["cluster"] = out["accession_version"].map(acc2cluster)
    out["split"] = out["cluster"].map(assign)
    validate_split_output(out)
    export_split_files(out, outdir)

    print(f"Wrote {outdir / 'release_cluster_split.csv'}")
    print(f"Wrote {outdir / 'train.csv'}")
    print(f"Wrote {outdir / 'val.csv'}")
    print(f"Wrote {outdir / 'test.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
