#!/usr/bin/env python3
# PUBLIC-UPDATE SOURCE NOTE:
# Preserved from the user-supplied data-pipeline source. It implements sequence QC,
# annotation-derived labels, precursor handling, and exact-sequence deduplication.
# It does not by itself reproduce the complete IPG representative-selection provenance
# described by the current manuscript. Therefore the frozen Zenodo v1.1 CSV is the
# authoritative input for reproducing the paper; this script is a procedural audit aid.

"""Clean, annotate, deduplicate, and audit downloaded GenPept records.

Important label semantics:
  AMP     = annotation-supported AMP according to the explicit lexicon below.
  non-AMP = candidate non-AMP / no AMP-associated annotation under this rule.

This script does not claim experimental inactivity for non-AMP records.
It preserves pre-balance counts and exclusion reasons for reviewer audit.
"""
from __future__ import annotations
import argparse, gzip, json, re
from collections import Counter
from pathlib import Path
import pandas as pd

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
LOW_QUALITY = re.compile(r"\b(fragment|hypothetical|putative|possible|predicted)\b|LOW\s+QUALITY\s+PROTEIN", re.I)
FALSE_AMP_GENES = re.compile(r"(?:^|[; ,])(ampC|ampD|ampR|ampS)(?:$|[; ,])", re.I)
AMP_TERMS = [
    r"\bantimicrobial[ -]+peptide\b",
    r"\banti[ -]?bacterial[ -]+peptide\b",
    r"\banti[ -]?fungal[ -]+peptide\b",
    r"\banti[ -]?parasitic[ -]+peptide\b",
    r"\bbacteriocin\b",
    r"\blantibiotic\b",
    r"\bmicrocin\b",
    r"\bthiopeptide\b",
    r"\bsactipeptide\b",
]
AMP_RE = re.compile("|".join(AMP_TERMS), re.I)
PRECURSOR_RE = re.compile(r"\b(precursor|prepropeptide|propeptide)\b", re.I)
MIN_LEN, MAX_LEN = 10, 200


def clean_sequence(seq: str) -> str:
    # Remove formatting whitespace only. Do not silently map/remove non-standard residues.
    return re.sub(r"\s+", "", str(seq).upper())


def iter_raw(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                yield json.loads(line)


def classify(row):
    text = " ".join([str(row.get("description", "")), str(row.get("cds_product", ""))])
    genes = str(row.get("cds_gene", ""))
    if FALSE_AMP_GENES.search(genes):
        return "non-AMP", "false_amp_gene_guard"
    m = AMP_RE.search(text)
    if m:
        return "AMP", m.group(0).lower()
    return "non-AMP", "no_amp_associated_annotation"


def parse_master_sequences(path: Path):
    seqs = []
    cur = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur: seqs.append("".join(cur)); cur=[]
            elif line:
                cur.append(line.upper())
    if cur: seqs.append("".join(cur))
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="outputs/00_retrieval/raw_genpept_records.jsonl.gz")
    ap.add_argument("--out-dir", default="outputs/00_curated")
    ap.add_argument("--exclude-precursors", action="store_true")
    ap.add_argument("--master-fasta", default="input/MASTER_DATASET_11000.fasta",
                    help="Optional release manifest by exact sequence. If present, audit whether all final 11k sequences are recoverable.")
    args = ap.parse_args()

    raw = Path(args.raw); out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    if not raw.exists():
        raise SystemExit(f"Missing raw retrieval file: {raw}")

    counts = Counter(); kept = []
    seen_acc = set()
    for r in iter_raw(raw):
        counts["raw_records"] += 1
        acc = str(r.get("accession_version", ""))
        if acc and acc in seen_acc:
            counts["excluded_duplicate_accession_version"] += 1; continue
        if acc: seen_acc.add(acc)
        seq = clean_sequence(r.get("sequence", ""))
        if not (MIN_LEN <= len(seq) <= MAX_LEN):
            counts["excluded_length"] += 1; continue
        if not seq or any(a not in AA20 for a in seq):
            counts["excluded_noncanonical_aa"] += 1; continue
        text = " ".join([str(r.get("description", "")), str(r.get("cds_product", ""))])
        if LOW_QUALITY.search(text):
            counts["excluded_annotation_quality"] += 1; continue
        precursor = bool(PRECURSOR_RE.search(text))
        if precursor and args.exclude_precursors:
            counts["excluded_precursor"] += 1; continue
        label, evidence = classify(r)
        x = dict(r)
        x.update(sequence=seq, length=len(seq), label=label, label_rule=evidence,
                 precursor_flag=precursor,
                 label_semantics=("annotation-supported AMP" if label == "AMP" else "candidate non-AMP; absence of AMP-associated annotation"))
        kept.append(x)
        counts["after_qc_before_exact_sequence_dedup"] += 1

    df = pd.DataFrame(kept)
    if df.empty: raise SystemExit("No records remain after QC.")

    # Exact-sequence deduplication, deterministic representative: smallest accession.version.
    df = df.sort_values(["sequence", "accession_version"], kind="mergesort")
    dup_sizes = df.groupby("sequence").size().rename("exact_sequence_group_size")
    df = df.drop_duplicates("sequence", keep="first").copy()
    df = df.join(dup_sizes, on="sequence")
    counts["after_exact_sequence_dedup"] = len(df)
    counts["amp_prebalance"] = int((df.label == "AMP").sum())
    counts["nonamp_prebalance"] = int((df.label == "non-AMP").sum())
    counts["precursor_flagged_after_qc"] = int(df.precursor_flag.sum())

    df.to_csv(out / "candidate_pool_prebalance.csv.gz", index=False, compression="gzip")
    pd.DataFrame([{"stage": k, "count": v} for k, v in counts.items()]).to_csv(out / "stage_counts.csv", index=False)

    master = Path(args.master_fasta)
    if master.exists():
        target = parse_master_sequences(master)
        target_set = set(target)
        recovered = df[df.sequence.isin(target_set)].copy()
        missing = sorted(target_set - set(recovered.sequence))
        recovered.to_csv(out / "master_11000_provenance_recovered.csv.gz", index=False, compression="gzip")
        (out / "master_recovery_summary.json").write_text(json.dumps({
            "master_sequences": len(target),
            "master_unique_sequences": len(target_set),
            "recovered_from_current_raw_qc": len(recovered),
            "missing_from_current_raw_qc": len(missing),
            "note": "A current NCBI re-download may differ from the frozen 2025 release because source records can change. The frozen master FASTA remains the release source of truth."
        }, indent=2), encoding="utf-8")
        if missing:
            pd.Series(missing, name="sequence").to_csv(out / "master_sequences_not_recovered.csv", index=False)
        print(f"Master recovery: {len(recovered):,}/{len(target_set):,}")

    print("Wrote", out / "candidate_pool_prebalance.csv.gz")
    print("Wrote", out / "stage_counts.csv")


if __name__ == "__main__":
    main()
