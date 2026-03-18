#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
from pathlib import Path
from typing import Iterator

from Bio import Entrez, SeqIO

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")
TAXA = {
    "bacteria": {"query": "txid2[orgn]", "taxid": 2, "name": "Bacteria"},
    "archaea": {"query": "txid2157[orgn]", "taxid": 2157, "name": "Archaea"},
    "fungi": {"query": "txid4751[orgn]", "taxid": 4751, "name": "Fungi"},
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Retrieve raw GenPept records for the GenPept-10-200-2025 dataset."
    )
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"), help="NCBI email")
    ap.add_argument("--date-start", default="2025/04/01", help="Start date (YYYY/MM/DD)")
    ap.add_argument("--date-end", default="2025/07/31", help="End date (YYYY/MM/DD)")
    ap.add_argument("--min-len", type=int, default=10)
    ap.add_argument("--max-len", type=int, default=200)
    ap.add_argument("--page-retmax", type=int, default=20000)
    ap.add_argument("--batch-efetch", type=int, default=400)
    ap.add_argument("--sleep", type=float, default=0.34)
    ap.add_argument(
        "--taxa",
        nargs="+",
        default=["bacteria", "archaea", "fungi"],
        help="Subset of: bacteria archaea fungi",
    )
    return ap.parse_args()


def make_query(taxon_term: str, date_start: str, date_end: str, min_len: int, max_len: int) -> str:
    src_term = "srcdb_genbank[PROP]"
    len_term = f"{min_len}:{max_len}[SLEN]"
    date_term = f'("{date_start}"[PDAT] : "{date_end}"[PDAT])'
    return f"({taxon_term}) AND {date_term} AND {src_term} AND {len_term}"


def esearch_count(term: str) -> int:
    with Entrez.esearch(db="protein", term=term, retmax=0) as handle:
        rec = Entrez.read(handle)
    return int(rec.get("Count", "0"))


def esearch_all_ids(term: str, page_retmax: int, sleep: float) -> list[str]:
    total = esearch_count(term)
    ids: list[str] = []
    for start in range(0, total, page_retmax):
        with Entrez.esearch(db="protein", term=term, retstart=start, retmax=page_retmax) as handle:
            rec = Entrez.read(handle)
        ids.extend(rec.get("IdList", []))
        time.sleep(sleep)
    return ids


def iter_fetch_records(id_list: list[str], batch_size: int, sleep: float) -> Iterator:
    for start in range(0, len(id_list), batch_size):
        batch = id_list[start : start + batch_size]
        if not batch:
            continue
        with Entrez.efetch(db="protein", id=",".join(batch), rettype="gb", retmode="text") as handle:
            for rec in SeqIO.parse(handle, "genbank"):
                yield rec
        time.sleep(sleep)


def clean_seq(seq: str) -> str:
    return re.sub(r"[\s-]", "", str(seq).strip().upper())


def extract_tax_id(rec) -> int | None:
    for feat in rec.features:
        if feat.type != "source":
            continue
        for x in feat.qualifiers.get("db_xref", []):
            if x.startswith("taxon:"):
                try:
                    return int(x.split(":", 1)[1])
                except ValueError:
                    return None
        break
    return None


def record_to_row(rec) -> dict:
    accession = rec.annotations.get("accessions", [rec.id])[0]
    version = rec.annotations.get("sequence_version")
    accession_version = f"{accession}.{version}" if version and not accession.endswith(f".{version}") else accession

    cds_products: list[str] = []
    cds_genes: list[str] = []
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        cds_products.extend(feat.qualifiers.get("product", []))
        cds_genes.extend(feat.qualifiers.get("gene", []))

    tax_id = extract_tax_id(rec)
    taxon_name = next((v["name"] for v in TAXA.values() if v["taxid"] == tax_id), f"taxid:{tax_id}" if tax_id else "")

    refseq_accession = ""
    for x in rec.annotations.get("dbxrefs", []):
        if x.startswith("RefSeq:"):
            refseq_accession = x.split(":", 1)[1]
            break

    sequence = clean_seq(str(rec.seq))
    return {
        "accession_version": accession_version,
        "description": rec.description or "",
        "sequence": sequence,
        "length": len(sequence),
        "contains_only_20aa": all(ch in ALLOWED_AA for ch in sequence),
        "tax_id": tax_id if tax_id is not None else "",
        "taxon_name": taxon_name,
        "refseq_accession": refseq_accession,
        "cds_products": ";".join(cds_products),
        "cds_genes": ";".join(cds_genes),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = [
        "accession_version",
        "description",
        "sequence",
        "length",
        "contains_only_20aa",
        "tax_id",
        "taxon_name",
        "refseq_accession",
        "cds_products",
        "cds_genes",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8") as f:
        for row in rows:
            f.write(f">{row['accession_version']} taxid={row['tax_id']}\n{row['sequence']}\n")


def main() -> int:
    args = parse_args()
    if not args.email:
        print("ERROR: provide --email or set NCBI_EMAIL", file=sys.stderr)
        return 2

    Entrez.email = args.email
    Entrez.tool = "GenPept-Curated-2025-release"

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw_taxon_rows: list[dict] = []
    raw_ids: list[str] = []
    seen_ids: set[str] = set()

    for taxon in args.taxa:
        taxon_key = taxon.lower()
        if taxon_key not in TAXA:
            print(f"ERROR: unsupported taxon preset: {taxon}", file=sys.stderr)
            return 2
        meta = TAXA[taxon_key]
        term = make_query(meta["query"], args.date_start, args.date_end, args.min_len, args.max_len)
        count = esearch_count(term)
        raw_taxon_rows.append({"taxon": meta["name"], "raw_records": count})
        ids = esearch_all_ids(term, page_retmax=args.page_retmax, sleep=args.sleep)
        for uid in ids:
            if uid not in seen_ids:
                raw_ids.append(uid)
                seen_ids.add(uid)

    rows = [record_to_row(rec) for rec in iter_fetch_records(raw_ids, batch_size=args.batch_efetch, sleep=args.sleep)]

    raw_taxon_csv = outdir / "01_raw_taxon_summary.csv"
    raw_records_csv = outdir / "01_raw_records.csv"
    raw_records_fasta = outdir / "01_raw_records.fasta"
    toc_csv = outdir / "01_toc_left_values.csv"

    total_raw = sum(int(r["raw_records"]) for r in raw_taxon_rows)
    for row in raw_taxon_rows:
        row["pct_raw"] = round((row["raw_records"] / total_raw * 100.0), 2) if total_raw else 0.0

    with raw_taxon_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["taxon", "raw_records", "pct_raw"])
        writer.writeheader()
        writer.writerows(raw_taxon_rows)

    with toc_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["raw_records_total", "Bacteria_pct", "Archaea_pct", "Fungi_pct", "date_window", "length_window"],
        )
        writer.writeheader()
        pct = {row["taxon"]: row["pct_raw"] for row in raw_taxon_rows}
        writer.writerow(
            {
                "raw_records_total": total_raw,
                "Bacteria_pct": pct.get("Bacteria", 0.0),
                "Archaea_pct": pct.get("Archaea", 0.0),
                "Fungi_pct": pct.get("Fungi", 0.0),
                "date_window": f"{args.date_start} to {args.date_end}",
                "length_window": f"{args.min_len}-{args.max_len} aa",
            }
        )

    write_csv(raw_records_csv, rows)
    write_fasta(raw_records_fasta, rows)

    print(f"Wrote {raw_records_csv}")
    print(f"Wrote {raw_records_fasta}")
    print(f"Wrote {raw_taxon_csv}")
    print(f"Wrote {toc_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
