from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CANONICAL = ROOT / "data" / "GenPept_Curated_2025_primary_split_v1.1.csv"
EXPECTED_CANONICAL_SHA256 = "421e55265e1462052f633c961f8d8cc20ce5e510284ffc36afc6e29a40b79b6c"
EXPECTED_HEADER = [
    "label", "sequence", "length", "row_id", "sample_id", "sequence_clean",
    "len_clean", "cluster_id", "comp_id", "split", "length_bin"
]
ALLOWED = set("ACDEFGHIKLMNPQRSTVWY")
EXPECTED_SPLIT_HASHES = {
    "train": "0df35d882928d02d174141a8919c18ca838b3d390c612e4599f3819e49c69192",
    "val": "088a3efa9d9a31a7219aa4edc4c97f182d72744999fe5b51068b54fff5f9101f",
    "test": "97815940ed08b079ac93443b6e16dcd34b2b583fae3d2189d16201c60bed5c71",
}
SPLIT_FILES = {
    "train": ROOT / "data" / "splits" / "train.csv",
    "val": ROOT / "data" / "splits" / "validation.csv",
    "test": ROOT / "data" / "splits" / "test.csv",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def main() -> None:
    assert sha256(CANONICAL) == EXPECTED_CANONICAL_SHA256
    header, rows = load_csv(CANONICAL)
    assert header == EXPECTED_HEADER
    assert len(rows) == 11000

    sample_ids = [r["sample_id"] for r in rows]
    sequences = [r["sequence"] for r in rows]
    row_ids = [r["row_id"] for r in rows]
    assert len(set(sample_ids)) == 11000
    assert len(set(sequences)) == 11000
    assert len(set(row_ids)) == 11000

    labels = Counter(r["label"] for r in rows)
    splits = Counter(r["split"] for r in rows)
    split_labels = Counter((r["split"], r["label"]) for r in rows)
    assert labels == Counter({"AMP": 5500, "non-AMP": 5500})
    assert splits == Counter({"train": 7700, "test": 2310, "val": 990})
    assert split_labels == Counter({
        ("train", "AMP"): 3850, ("train", "non-AMP"): 3850,
        ("val", "AMP"): 495, ("val", "non-AMP"): 495,
        ("test", "AMP"): 1155, ("test", "non-AMP"): 1155,
    })

    comp_splits: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        seq = r["sequence"]
        n = int(r["length"])
        assert 10 <= n <= 200
        assert len(seq) == n
        assert set(seq) <= ALLOWED
        assert r["sequence_clean"] == seq
        assert int(r["len_clean"]) == n
        comp_splits[r["comp_id"]].add(r["split"])
    assert len(comp_splits) == 10666
    assert max(len(v) for v in comp_splits.values()) == 1

    for split, path in SPLIT_FILES.items():
        assert sha256(path) == EXPECTED_SPLIT_HASHES[split]
        split_header, exported = load_csv(path)
        assert split_header == EXPECTED_HEADER
        expected = [r for r in rows if r["split"] == split]
        assert exported == expected

    print("PASS: GenPept-Curated-2025 Zenodo dataset v1.1 validated.")


if __name__ == "__main__":
    main()
