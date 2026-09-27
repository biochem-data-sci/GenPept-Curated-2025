#!/usr/bin/env python3
"""Execute the preserved ten-branch NCBI Protein query specification.

The query CSV is authoritative. This helper records live counts and can fetch a
bounded set of GenBank records per branch. Live counts are not frozen results
because NCBI Protein changes after the original retrieval date.
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from pathlib import Path

from Bio import Entrez, SeqIO


REQUIRED_COLUMNS = {"cohort", "label_branch", "length_bin", "precursor_constraint", "query"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--query-spec",
        default="reproducibility/query_specification_from_preserved_rules.csv",
        help="CSV containing the ten preserved executable queries",
    )
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    parser.add_argument(
        "--fetch-records",
        action="store_true",
        help="Fetch GenBank records after counting; omitted by default because the search space is large",
    )
    parser.add_argument(
        "--retmax-per-query",
        type=int,
        default=0,
        help="Maximum records fetched per query; 0 is valid only for count-only mode",
    )
    parser.add_argument("--batch-size", type=int, default=300)
    parser.add_argument("--sleep", type=float, default=0.34)
    return parser.parse_args()


def load_spec(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("Query specification is empty")
    missing = REQUIRED_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f"Query specification is missing columns: {sorted(missing)}")
    if len(rows) != 10:
        raise ValueError(f"Expected exactly 10 preserved query specifications, found {len(rows)}")
    if any(not row["query"].strip() for row in rows):
        raise ValueError("Every query row must contain an executable query string")
    return rows


def esearch(query: str, retmax: int) -> tuple[int, list[str]]:
    with Entrez.esearch(db="protein", term=query, retmax=retmax) as handle:
        record = Entrez.read(handle)
    return int(record["Count"]), list(record.get("IdList", []))


def fetch_records(ids: list[str], path: Path, batch_size: int, pause: float) -> int:
    written = 0
    with path.open("w", encoding="utf-8") as out:
        for start in range(0, len(ids), batch_size):
            batch = ids[start : start + batch_size]
            with Entrez.efetch(db="protein", id=",".join(batch), rettype="gb", retmode="text") as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    SeqIO.write(record, out, "genbank")
                    written += 1
            time.sleep(pause)
    return written


def main() -> int:
    args = parse_args()
    if not args.email:
        raise SystemExit("Provide --email or set NCBI_EMAIL")
    if args.retmax_per_query < 0:
        raise SystemExit("--retmax-per-query must be >= 0")
    if args.fetch_records and args.retmax_per_query == 0:
        raise SystemExit("Full retrieval is intentionally explicit: set a positive --retmax-per-query")

    Entrez.email = args.email
    Entrez.api_key = args.api_key
    Entrez.tool = "GenPept-Curated-2025-v1.1"

    rows = load_spec(Path(args.query_spec))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    count_rows: list[dict[str, object]] = []

    for index, row in enumerate(rows, start=1):
        fetch_limit = args.retmax_per_query if args.fetch_records else 0
        count, ids = esearch(row["query"], fetch_limit)
        branch = f"{index:02d}_{row['cohort']}_{row['label_branch']}_{row['length_bin']}".replace("/", "-")
        fetched = 0
        if args.fetch_records:
            fetched = fetch_records(ids, outdir / f"{branch}.gb", args.batch_size, args.sleep)
        count_rows.append(
            {
                "query_index": index,
                "cohort": row["cohort"],
                "label_branch": row["label_branch"],
                "length_bin": row["length_bin"],
                "precursor_constraint": row["precursor_constraint"],
                "live_count": count,
                "requested_fetch_cap": fetch_limit,
                "fetched_records": fetched,
                "query": row["query"],
            }
        )
        time.sleep(args.sleep)

    output = outdir / "query_counts.csv"
    with output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(count_rows[0]))
        writer.writeheader()
        writer.writerows(count_rows)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
