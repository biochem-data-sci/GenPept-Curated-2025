#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import io
import os
import re
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd
import requests

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")
BIN_EDGES = [10, 50, 100, 150, 200]
BIN_LABELS = ["10-50", "51-100", "101-150", "151-200"]
LOW_QUALITY_PAT = re.compile(r"(fragment|hypothetical|putative|possible|predicted|LOW\s*QUALITY\s*PROTEIN)", re.I)
PRECURSOR_PAT = re.compile(r"\b(precursor|propeptide|preprotein|preproprotein|pre-proprotein|proprotein)\b", re.I)
FALSE_AMP_PAT = re.compile(r"\b(ampC|ampD|ampR|ampS)\b", re.I)
AMP_PAT = re.compile(
    "|".join(
        [
            r"\bantimicrobial peptide\b",
            r"\bantimicrobial\b",
            r"\bantibacterial peptide\b",
            r"\banti[ -]?bacterial\b",
            r"\bantifungal peptide\b",
            r"\banti[ -]?fungal\b",
            r"\bantiparasitic peptide\b",
            r"\banti[ -]?parasitic\b",
            r"\bbacteriocin\b",
            r"\blantibiotic\b",
            r"\bmicrocin\b",
            r"\bthiopeptide\b",
            r"\bsactipeptide\b",
            r"\bhost defense peptide\b",
            r"\bbactericidal\b",
            r"(?<![A-Za-z])AMP(?![A-Za-z])",
        ]
    ),
    re.I,
)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Apply QC, precursor flagging, annotation-based labeling, deduplication, and balancing."
    )
    ap.add_argument("--input-csv", required=True, help="01_raw_records.csv from step 01")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--target-total", type=int, default=11000)
    ap.add_argument("--target-per-label", type=int, default=5500)
    ap.add_argument("--min-len", type=int, default=10)
    ap.add_argument("--max-len", type=int, default=200)
    ap.add_argument("--dedup-mode", choices=["ipg", "exact", "none"], default="ipg")
    ap.add_argument("--fallback-to-exact", action="store_true", help="Fallback to exact-sequence dedup if IPG lookup fails")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"), help="NCBI email for IPG lookups")
    ap.add_argument("--sleep", type=float, default=0.34)
    ap.add_argument("--ipg-cache", default=None, help="Optional TSV cache for IPG lookups")
    return ap.parse_args()


def clean_seq(seq: str) -> str:
    return re.sub(r"[\s-]", "", str(seq).strip().upper())


def score_description(desc: str) -> tuple[int, int, str]:
    text = re.sub(r"\s+", " ", str(desc or "").strip())
    informative = 0 if not text else len(re.findall(r"[A-Za-z0-9]+", text))
    penalty = 1 if LOW_QUALITY_PAT.search(text) else 0
    return (penalty, -informative, text.lower())


def classify_qc(row: pd.Series, min_len: int, max_len: int) -> str:
    seq = row["sequence"]
    desc = str(row.get("description", "") or "")
    if not (min_len <= len(seq) <= max_len):
        return "fail_length"
    if not seq or any(ch not in ALLOWED_AA for ch in seq):
        return "fail_noncanonical"
    if LOW_QUALITY_PAT.search(desc):
        return "fail_low_quality_annotation"
    if PRECURSOR_PAT.search(desc):
        return "flag_precursor"
    return "pass"


def label_record(row: pd.Series) -> tuple[str, str]:
    desc = " ".join(
        [
            str(row.get("description", "") or ""),
            str(row.get("cds_products", "") or ""),
        ]
    )
    genes = str(row.get("cds_genes", "") or "")
    if FALSE_AMP_PAT.search(desc) or FALSE_AMP_PAT.search(genes):
        return "non-AMP", "false_amp_name"
    if AMP_PAT.search(desc):
        return "AMP", "annotation_keyword"
    return "non-AMP", "no_amp_keyword"


def read_cache(path: Path) -> dict[str, dict[str, str]]:
    cache: dict[str, dict[str, str]] = {}
    if not path.exists():
        return cache
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="	")
        for row in reader:
            cache[row["accession_version"]] = row
    return cache


def write_cache(path: Path, cache: dict[str, dict[str, str]]) -> None:
    fieldnames = ["accession_version", "ipg_id", "refseq_accession", "lookup_status"]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="	")
        writer.writeheader()
        for key in sorted(cache):
            writer.writerow({k: cache[key].get(k, "") for k in fieldnames})


def fetch_ipg_table(accession_version: str, email: str | None, sleep: float) -> pd.DataFrame:
    params_list = [
        {"db": "protein", "report": "ipg", "retmode": "text", "val": accession_version},
        {"db": "protein", "report": "ipg", "retmode": "text", "id": accession_version},
    ]
    headers = {"User-Agent": f"GenPept-Curated-2025-release ({email or 'no-email'})"}
    last_error: Exception | None = None
    for params in params_list:
        try:
            resp = requests.get(
                "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi",
                params=params,
                headers=headers,
                timeout=60,
            )
            resp.raise_for_status()
            text = resp.text.strip()
            if not text:
                continue
            sep = "	" if "	" in text.splitlines()[0] else ","
            df = pd.read_csv(io.StringIO(text), sep=sep)
            if df.empty:
                continue
            time.sleep(sleep)
            return df
        except Exception as exc:  # pragma: no cover - network path
            last_error = exc
            continue
    if last_error is not None:
        raise last_error
    return pd.DataFrame()


def normalize_ipg_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(c).strip().lower().replace(" ", "_") for c in out.columns]
    return out


def choose_refseq_accession(ipg_df: pd.DataFrame, taxon_name: str, existing: str) -> str:
    if existing:
        return existing
    if ipg_df.empty or "protein" not in ipg_df.columns:
        return existing
    proteins = [str(x).strip() for x in ipg_df["protein"].dropna().tolist()]
    if taxon_name in {"Bacteria", "Archaea"}:
        for acc in proteins:
            if acc.startswith("WP_"):
                return acc
    if taxon_name == "Fungi":
        for acc in proteins:
            if acc.startswith("NP_") or acc.startswith("XP_"):
                return acc
    return existing


def lookup_ipg(accession_version: str, taxon_name: str, existing_refseq: str, cache: dict[str, dict[str, str]], email: str | None, sleep: float) -> tuple[str, str, str]:
    cached = cache.get(accession_version)
    if cached is not None:
        return cached.get("ipg_id", ""), cached.get("refseq_accession", ""), cached.get("lookup_status", "")

    ipg_id = ""
    refseq_accession = existing_refseq or ""
    lookup_status = "not_found"
    ipg_df = fetch_ipg_table(accession_version, email=email, sleep=sleep)
    if not ipg_df.empty:
        ipg_df = normalize_ipg_columns(ipg_df)
        if "id" in ipg_df.columns:
            ipg_id = str(ipg_df["id"].iloc[0]).strip()
        refseq_accession = choose_refseq_accession(ipg_df, taxon_name=taxon_name, existing=refseq_accession)
        lookup_status = "ok"
    cache[accession_version] = {
        "accession_version": accession_version,
        "ipg_id": ipg_id,
        "refseq_accession": refseq_accession,
        "lookup_status": lookup_status,
    }
    return ipg_id, refseq_accession, lookup_status


def allocate_balanced_counts(df: pd.DataFrame, target_total: int) -> dict[str, int]:
    max_bal: dict[str, int] = {}
    total_bal = 0
    for b in BIN_LABELS:
        c = df[df["length_bin"] == b]["label"].value_counts()
        m = int(min(c.get("AMP", 0), c.get("non-AMP", 0)))
        max_bal[b] = 2 * m
        total_bal += 2 * m
    if total_bal < target_total:
        raise RuntimeError(
            f"Balanced capacity is {total_bal:,}, below target_total={target_total:,}."
        )

    raw_alloc = {b: target_total * (max_bal[b] / total_bal) for b in BIN_LABELS}
    alloc = {b: 2 * int(round(raw_alloc[b] / 2)) for b in BIN_LABELS}
    diff = target_total - sum(alloc.values())
    room = {b: max_bal[b] - alloc[b] for b in BIN_LABELS}
    step = 2 if diff > 0 else -2
    while diff != 0:
        moved = False
        for b in sorted(BIN_LABELS, key=lambda x: room[x], reverse=(diff > 0)):
            if (diff > 0 and room[b] >= 2) or (diff < 0 and alloc[b] >= 2):
                alloc[b] += step
                diff -= step
                moved = True
                if diff == 0:
                    break
        if not moved:
            break
    return alloc


def write_fasta(path: Path, df: pd.DataFrame) -> None:
    with path.open("w", encoding="utf-8") as f:
        for row in df.itertuples(index=False):
            f.write(f">{row.accession_version}|{row.label}|len:{row.length}\n{row.sequence}\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.target_total != 2 * args.target_per_label:
        print("ERROR: target_total must equal 2 * target_per_label for a balanced release.", file=sys.stderr)
        return 2

    cache_path = Path(args.ipg_cache) if args.ipg_cache else outdir / "ipg_lookup_cache.tsv"
    cache = read_cache(cache_path)

    df = pd.read_csv(args.input_csv)
    required = {"accession_version", "description", "sequence"}
    missing = required - set(df.columns)
    if missing:
        print(f"ERROR: input CSV missing columns: {sorted(missing)}", file=sys.stderr)
        return 2

    for col in ["cds_products", "cds_genes", "tax_id", "taxon_name", "refseq_accession"]:
        if col not in df.columns:
            df[col] = ""

    df["sequence"] = df["sequence"].map(clean_seq)
    df["length"] = df["sequence"].map(len)
    df["qc_status"] = df.apply(lambda r: classify_qc(r, min_len=args.min_len, max_len=args.max_len), axis=1)
    df[["label", "label_reason"]] = df.apply(lambda r: pd.Series(label_record(r)), axis=1)
    df["precursor_flag"] = df["qc_status"].eq("flag_precursor")

    precursor_df = df[df["precursor_flag"]].copy()
    main_df = df[df["qc_status"].eq("pass")].copy()

    ipg_ids: list[str] = []
    refseq_accessions: list[str] = []
    lookup_statuses: list[str] = []
    dedup_keys: list[str] = []
    dedup_mode_applied: list[str] = []

    for row in main_df.itertuples(index=False):
        if args.dedup_mode == "ipg":
            try:
                ipg_id, refseq_accession, lookup_status = lookup_ipg(
                    row.accession_version,
                    taxon_name=str(row.taxon_name or ""),
                    existing_refseq=str(getattr(row, "refseq_accession", "") or ""),
                    cache=cache,
                    email=args.email,
                    sleep=args.sleep,
                )
                if ipg_id:
                    dedup_key = f"ipg:{ipg_id}"
                    mode_used = "ipg"
                elif args.fallback_to_exact:
                    dedup_key = f"seq:{row.sequence}"
                    mode_used = "exact_fallback"
                else:
                    dedup_key = f"acc:{row.accession_version}"
                    mode_used = "accession_only"
            except Exception:
                if args.fallback_to_exact:
                    ipg_id, refseq_accession, lookup_status = "", str(getattr(row, "refseq_accession", "") or ""), "fallback_exact"
                    dedup_key = f"seq:{row.sequence}"
                    mode_used = "exact_fallback"
                else:
                    raise
        elif args.dedup_mode == "exact":
            ipg_id, refseq_accession, lookup_status = "", str(getattr(row, "refseq_accession", "") or ""), "not_requested"
            dedup_key = f"seq:{row.sequence}"
            mode_used = "exact"
        else:
            ipg_id, refseq_accession, lookup_status = "", str(getattr(row, "refseq_accession", "") or ""), "not_requested"
            dedup_key = f"acc:{row.accession_version}"
            mode_used = "none"

        ipg_ids.append(ipg_id)
        refseq_accessions.append(refseq_accession)
        lookup_statuses.append(lookup_status)
        dedup_keys.append(dedup_key)
        dedup_mode_applied.append(mode_used)

    main_df["ipg_id"] = ipg_ids
    main_df["refseq_accession"] = refseq_accessions
    main_df["ipg_lookup_status"] = lookup_statuses
    main_df["dedup_group"] = dedup_keys
    main_df["dedup_mode_applied"] = dedup_mode_applied

    write_cache(cache_path, cache)

    rep_idx = (
        main_df.assign(_desc_score=main_df["description"].map(score_description))
        .sort_values(by=["dedup_group", "_desc_score", "accession_version"])
        .groupby("dedup_group", sort=False)
        .head(1)
        .index
    )
    curated_df = main_df.copy()
    curated_df["is_representative"] = False
    curated_df.loc[rep_idx, "is_representative"] = True
    dedup_df = curated_df[curated_df["is_representative"]].copy()

    dedup_df["length_bin"] = pd.cut(
        dedup_df["length"].astype(int),
        bins=BIN_EDGES,
        labels=BIN_LABELS,
        include_lowest=True,
    )
    dedup_df = dedup_df[dedup_df["length_bin"].notna()].copy()

    alloc = allocate_balanced_counts(dedup_df, target_total=args.target_total)
    parts: list[pd.DataFrame] = []
    for b in BIN_LABELS:
        sub = dedup_df[dedup_df["length_bin"] == b]
        per_label = alloc[b] // 2
        parts.append(sub[sub["label"] == "AMP"].sample(n=per_label, random_state=args.seed))
        parts.append(sub[sub["label"] == "non-AMP"].sample(n=per_label, random_state=args.seed))
    dataset_df = pd.concat(parts, ignore_index=True)
    dataset_df = dataset_df.sort_values(["label", "length_bin", "accession_version"]).reset_index(drop=True)

    export_cols = [
        "accession_version",
        "sequence",
        "length",
        "label",
        "length_bin",
        "tax_id",
        "taxon_name",
        "description",
        "cds_products",
        "cds_genes",
        "ipg_id",
        "refseq_accession",
        "precursor_flag",
    ]
    metadata_cols = [
        "accession_version",
        "description",
        "sequence",
        "length",
        "tax_id",
        "taxon_name",
        "cds_products",
        "cds_genes",
        "qc_status",
        "precursor_flag",
        "label",
        "label_reason",
        "ipg_id",
        "refseq_accession",
        "ipg_lookup_status",
        "dedup_group",
        "dedup_mode_applied",
        "is_representative",
    ]

    dataset_path = outdir / "dataset.csv"
    metadata_path = outdir / "metadata.csv"
    precursor_path = outdir / "precursor_flagged.csv"
    manifest_path = outdir / "accession_manifest.tsv"
    ipg_map_path = outdir / "ipg_mapping.tsv"
    summary_path = outdir / "dataset_summary_by_bin.csv"
    attrition_path = outdir / "attrition_summary.csv"
    all_fasta = outdir / "dataset.fasta"
    amp_fasta = outdir / "amp_sequences.fasta"
    nonamp_fasta = outdir / "nonamp_sequences.fasta"

    dataset_df[export_cols].to_csv(dataset_path, index=False)
    pd.concat([curated_df[metadata_cols], precursor_df.assign(ipg_id="", refseq_accession=precursor_df.get("refseq_accession", ""), ipg_lookup_status="", dedup_group="", dedup_mode_applied="", is_representative=False)[metadata_cols]], ignore_index=True).to_csv(metadata_path, index=False)
    precursor_df.assign(length=precursor_df["sequence"].map(len)).to_csv(precursor_path, index=False)

    manifest_df = curated_df[["accession_version", "label", "dedup_group", "ipg_id", "refseq_accession", "is_representative"]].copy()
    manifest_df["selected_in_dataset"] = manifest_df["accession_version"].isin(dataset_df["accession_version"])
    manifest_df.to_csv(manifest_path, sep="	", index=False)
    curated_df[["accession_version", "ipg_id", "refseq_accession", "dedup_group", "is_representative"]].drop_duplicates().to_csv(ipg_map_path, sep="	", index=False)

    summary = (
        dataset_df.groupby(["length_bin", "label"]).size().rename("count").reset_index().pivot(index="length_bin", columns="label", values="count").fillna(0).astype(int).reset_index()
    )
    summary["total"] = summary.get("AMP", 0) + summary.get("non-AMP", 0)
    summary.to_csv(summary_path, index=False)

    attr = [
        {"step": "raw_input_records", "n": int(len(df))},
        {"step": "removed_length", "n": int((df["qc_status"] == "fail_length").sum())},
        {"step": "removed_noncanonical", "n": int((df["qc_status"] == "fail_noncanonical").sum())},
        {"step": "removed_low_quality_annotation", "n": int((df["qc_status"] == "fail_low_quality_annotation").sum())},
        {"step": "precursor_flagged", "n": int(len(precursor_df))},
        {"step": "main_pool_after_qc", "n": int(len(main_df))},
        {"step": "representatives_after_dedup", "n": int(len(dedup_df))},
        {"step": "final_dataset_total", "n": int(len(dataset_df))},
        {"step": "final_dataset_amp", "n": int((dataset_df["label"] == "AMP").sum())},
        {"step": "final_dataset_nonamp", "n": int((dataset_df["label"] == "non-AMP").sum())},
    ]
    pd.DataFrame(attr).to_csv(attrition_path, index=False)

    write_fasta(all_fasta, dataset_df[export_cols])
    write_fasta(amp_fasta, dataset_df[dataset_df["label"] == "AMP"][export_cols])
    write_fasta(nonamp_fasta, dataset_df[dataset_df["label"] == "non-AMP"][export_cols])

    print(f"Wrote {dataset_path}")
    print(f"Wrote {metadata_path}")
    print(f"Wrote {precursor_path}")
    print(f"Wrote {manifest_path}")
    print(f"Wrote {ipg_map_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {attrition_path}")
    print(f"Wrote {amp_fasta}")
    print(f"Wrote {nonamp_fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
