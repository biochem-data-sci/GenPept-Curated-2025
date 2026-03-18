#!/usr/bin/env python3
from __future__ import annotations

import argparse
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REQUIRED_BASE = {"sequence", "split"}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Check cross-split homology leakage with MMseqs2.")
    ap.add_argument("--release-csv", required=True, help="release_cluster_split.csv from step 03")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mmseqs-bin", default="mmseqs")
    ap.add_argument("--min-pident", type=float, default=90.0)
    ap.add_argument("--min-cov", type=float, default=0.80, help="Coverage threshold on qcov or tcov")
    return ap.parse_args()


def resolve_wsl_path(path: Path) -> str:
    out = subprocess.run(["wsl", "wslpath", "-a", str(path.resolve())], capture_output=True, text=True, check=True)
    return out.stdout.strip()


def run_mmseqs(fasta_path: Path, out_m8: Path, tmp_dir: Path, mmseqs_bin: str, min_pident: float, min_cov: float) -> None:
    native = shutil.which(mmseqs_bin)
    if native:
        cmd = [
            mmseqs_bin,
            "easy-search",
            str(fasta_path),
            str(fasta_path),
            str(out_m8),
            str(tmp_dir),
            "--min-seq-id",
            str(min_pident / 100.0),
            "-c",
            str(min_cov),
            "--cov-mode",
            "5",
            "--format-output",
            "query,target,pident,qcov,tcov,qlen,tlen",
        ]
    elif platform.system().lower().startswith("win"):
        cmd = [
            "wsl",
            mmseqs_bin,
            "easy-search",
            resolve_wsl_path(fasta_path),
            resolve_wsl_path(fasta_path),
            resolve_wsl_path(out_m8),
            resolve_wsl_path(tmp_dir),
            "--min-seq-id",
            str(min_pident / 100.0),
            "-c",
            str(min_cov),
            "--cov-mode",
            "5",
            "--format-output",
            "query,target,pident,qcov,tcov,qlen,tlen",
        ]
    else:
        raise RuntimeError("MMseqs2 not found in PATH")
    subprocess.run(cmd, check=True)


def write_fasta(df: pd.DataFrame, path: Path, acc_col: str) -> None:
    records = [
        SeqRecord(Seq(str(r.sequence)), id=f"{getattr(r, acc_col)}|{r.split}", description="")
        for r in df.itertuples(index=False)
    ]
    SeqIO.write(records, str(path), "fasta")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.release_csv)
    missing = REQUIRED_BASE - set(df.columns)
    if missing:
        print(f"ERROR: release CSV missing columns: {sorted(missing)}", file=sys.stderr)
        return 2
    acc_col = "accession_version" if "accession_version" in df.columns else "accession"
    if acc_col not in df.columns:
        print("ERROR: release CSV requires accession_version or accession", file=sys.stderr)
        return 2

    fasta_path = outdir / "mmseqs_input.fasta"
    out_m8 = outdir / "mmseqs_all_vs_all.tsv"
    tmp_dir = outdir / "mmseqs_tmp"
    write_fasta(df, fasta_path, acc_col=acc_col)
    run_mmseqs(fasta_path, out_m8, tmp_dir, args.mmseqs_bin, args.min_pident, args.min_cov)

    hits = pd.read_csv(out_m8, sep="	", header=None, names=["query", "target", "pident", "qcov", "tcov", "qlen", "tlen"])
    hits = hits[hits["query"] != hits["target"]].copy()
    hits["pair"] = hits.apply(lambda r: tuple(sorted((r["query"], r["target"]))), axis=1)
    hits = hits.drop_duplicates(subset=["pair"]).copy()
    hits[["query_accession", "query_split"]] = hits["query"].str.split("|", n=1, expand=True)
    hits[["target_accession", "target_split"]] = hits["target"].str.split("|", n=1, expand=True)

    hits["shorter_cov"] = hits.apply(
        lambda r: float(r["qcov"]) if float(r["qlen"]) <= float(r["tlen"]) else float(r["tcov"]),
        axis=1,
    )

    leaks = hits[
        (hits["pident"] >= args.min_pident)
        & (hits["shorter_cov"] >= args.min_cov)
        & (hits["query_split"] != hits["target_split"])
    ].copy()

    all_hits_csv = outdir / "mmseqs_all_pairs.csv"
    leak_csv = outdir / "mmseqs_cross_split_leaks.csv"
    summary_txt = outdir / "mmseqs_leak_summary.txt"
    hits.to_csv(all_hits_csv, index=False)
    leaks.to_csv(leak_csv, index=False)
    with summary_txt.open("w", encoding="utf-8") as f:
        f.write(f"all_pairs={len(hits)}\n")
        f.write(f"cross_split_leaks={len(leaks)}\n")
        f.write(f"min_pident={args.min_pident}\n")
        f.write(f"min_cov={args.min_cov}\n")
    print(f"Wrote {all_hits_csv}")
    print(f"Wrote {leak_csv}")
    print(f"Wrote {summary_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
