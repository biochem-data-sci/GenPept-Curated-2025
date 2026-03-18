#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Reproduce summary tables and figures from AMP/non-AMP FASTA files or dataset CSV.")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--dataset-csv", help="Dataset CSV with label and sequence columns")
    group.add_argument("--amp", help="AMP FASTA")
    ap.add_argument("--nonamp", help="non-AMP FASTA (required when --amp is used)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--fingerprint-accession", default=None)
    ap.add_argument("--fingerprint-index", type=int, default=0)
    return ap.parse_args()


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    seqs: list[tuple[str, str]] = []
    header = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(chunks).strip()))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        seqs.append((header, "".join(chunks).strip()))
    return seqs


def normalize_label(value: object) -> str:
    s = str(value).strip().replace("_", "-").replace(" ", "").lower()
    if s == "amp":
        return "AMP"
    if s in {"non-amp", "nonamp", "namp"}:
        return "non-AMP"
    return str(value)


def load_from_dataset_csv(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    if "label" not in cols or "sequence" not in cols:
        raise ValueError("dataset CSV must contain 'label' and 'sequence' columns")
    label_col = cols["label"]
    seq_col = cols["sequence"]
    id_col = cols.get("accession_version") or cols.get("accession") or cols.get("id")
    amp, nonamp = [], []
    for i, row in df.iterrows():
        label = normalize_label(row[label_col])
        seq = str(row[seq_col]).strip().upper()
        if not seq:
            continue
        if id_col:
            header = str(row[id_col])
        else:
            header = f"{label}_{i+1:06d}"
        if label == "AMP":
            amp.append((header, seq))
        elif label == "non-AMP":
            nonamp.append((header, seq))
    return amp, nonamp


def first_token_id(header: str) -> str:
    return header.split()[0] if header else ""


def alphabet_checks(seq: str) -> dict[str, bool]:
    non_std = set("UOBJZX")
    illegal = set("*- .\t\r\n")
    seq_u = seq.upper()
    return {
        "has_nonstd": any(c in non_std for c in seq_u),
        "has_illegal": any(c in illegal for c in seq_u),
        "has_outside20": any((c not in AA20) and (c not in non_std) and (c not in illegal) for c in seq_u),
    }


def build_qc_tables(amp: list[tuple[str, str]], nonamp: list[tuple[str, str]], outdir: Path) -> None:
    rows = []
    for label, seqs in (("AMP", amp), ("non-AMP", nonamp)):
        lengths = [len(s) for _, s in seqs]
        ids = [first_token_id(h) for h, _ in seqs]
        rows.append(
            {
                "set": label,
                "n_sequences": len(seqs),
                "min_len": min(lengths) if lengths else 0,
                "median_len": float(np.median(lengths)) if lengths else 0.0,
                "mean_len": round(float(np.mean(lengths)), 2) if lengths else 0.0,
                "max_len": max(lengths) if lengths else 0,
                "dup_ids_within_set": sum(1 for c in Counter(ids).values() if c > 1),
                "any_nonstd_aa": sum(1 for _, s in seqs if alphabet_checks(s)["has_nonstd"]),
                "any_illegal_chars": sum(1 for _, s in seqs if alphabet_checks(s)["has_illegal"]),
                "any_outside20": sum(1 for _, s in seqs if alphabet_checks(s)["has_outside20"]),
            }
        )
    pd.DataFrame(rows).to_csv(outdir / "counts_summary.csv", index=False)


def plot_counts(amp: list[tuple[str, str]], nonamp: list[tuple[str, str]], outdir: Path) -> None:
    plt.figure()
    plt.bar(["AMP", "non-AMP"], [len(amp), len(nonamp)])
    plt.title("Sequence counts")
    plt.xlabel("Label")
    plt.ylabel("Count")
    plt.savefig(outdir / "amp_nonamp_bar.png", bbox_inches="tight")
    plt.close()


def plot_lengths(amp: list[tuple[str, str]], nonamp: list[tuple[str, str]], outdir: Path) -> None:
    amp_len = [len(s) for _, s in amp]
    non_len = [len(s) for _, s in nonamp]
    plt.figure()
    plt.hist(amp_len, bins=30, alpha=0.5, label="AMP")
    plt.hist(non_len, bins=30, alpha=0.5, label="non-AMP")
    plt.title("Length distribution")
    plt.xlabel("Length (aa)")
    plt.ylabel("Count")
    plt.legend()
    plt.savefig(outdir / "length_distribution.png", bbox_inches="tight")
    plt.close()

    plt.figure()
    plt.hist(amp_len, bins=list(range(1, 51)), alpha=0.5, label="AMP")
    plt.hist(non_len, bins=list(range(1, 51)), alpha=0.5, label="non-AMP")
    plt.title("Length distribution (10-50 aa)")
    plt.xlabel("Length (aa)")
    plt.ylabel("Count")
    plt.xlim(1, 50)
    plt.legend()
    plt.savefig(outdir / "length_distribution_10_50.png", bbox_inches="tight")
    plt.close()


def save_length_bins(amp: list[tuple[str, str]], nonamp: list[tuple[str, str]], outdir: Path) -> None:
    lo, hi = 10, 200

    def per_len(lengths: list[int]) -> list[int]:
        c = Counter(L for L in lengths if lo <= L <= hi)
        return [c.get(L, 0) for L in range(lo, hi + 1)]

    amp_len = [len(s) for _, s in amp]
    non_len = [len(s) for _, s in nonamp]
    df = pd.DataFrame({"length_aa": list(range(lo, hi + 1)), "AMP": per_len(amp_len), "non-AMP": per_len(non_len)})
    df["total"] = df["AMP"] + df["non-AMP"]
    df.to_csv(outdir / "length_bins_10_200.csv", index=False)


def aa_freqs(amp: list[tuple[str, str]], nonamp: list[tuple[str, str]], outdir: Path) -> None:
    aa_order = list("ACDEFGHIKLMNPQRSTVWY") + ["X"]
    valid = set(aa_order)
    all_sequences = [s.upper() for _, s in amp] + [s.upper() for _, s in nonamp]
    counts = Counter()
    for seq in all_sequences:
        for aa in seq:
            if aa in valid:
                counts[aa] += 1
    for aa in aa_order:
        counts.setdefault(aa, 0)
    df = pd.DataFrame({"Amino Acid": list(counts.keys()), "Count": list(counts.values())}).sort_values("Count", ascending=False)
    df.to_csv(outdir / "aa_frequency_all.csv", index=False)
    plt.figure()
    plt.bar(df["Amino Acid"], df["Count"])
    plt.title("Amino acid frequency")
    plt.xlabel("Amino acid")
    plt.ylabel("Count")
    plt.savefig(outdir / "amino_acid_frequency.png", bbox_inches="tight")
    plt.close()


def property_panels(outdir: Path) -> None:
    AA = list("ACDEFGHIKLMNPQRSTVWY")
    kd = {"A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3}
    vdw = {"A": 88, "R": 173, "N": 114, "D": 111, "C": 108, "Q": 143, "E": 138, "G": 60, "H": 153, "I": 166, "L": 166, "K": 168, "M": 162, "F": 189, "P": 112, "S": 89, "T": 116, "W": 227, "Y": 193, "V": 140}
    nonpolar, polar_uncharged, charged = set("AVILMFWGC"), set("STNQY"), set("DEKRH")
    pol_low, pol_med = set("ADGPST"), set("CEHNQV")
    helix, sheet, coil = set("AELMQKRH"), set("VIYFWT"), set("GNPSDCT")
    buried, exposed = set("ALFCGIVWYM"), set("DEKRQNHSTP")

    attribs = [
        ("Hydrophobicity", lambda a: "Hydrophobic" if kd[a] > 0 else ("Neutral" if -0.4 <= kd[a] <= 0.4 else "Polar"), ["Polar", "Neutral", "Hydrophobic"]),
        ("Vdw volume", lambda a: "Small" if vdw[a] < 115 else ("Medium" if vdw[a] <= 155 else "Large"), ["Small", "Medium", "Large"]),
        ("Polarity", lambda a: "Polar" if a in charged else ("Neutral" if a in polar_uncharged else "Nonpolar"), ["Polar", "Neutral", "Nonpolar"]),
        ("Polarizability", lambda a: "Low" if a in pol_low else ("Medium" if a in pol_med else "High"), ["Low", "Medium", "High"]),
        ("Charge", lambda a: "Negative" if a in ("D", "E") else ("Positive" if a in ("K", "R", "H") else "Neutral"), ["Negative", "Neutral", "Positive"]),
        ("Secondary structure", lambda a: "Helix" if a in helix else ("Sheet" if a in sheet else "Coil"), ["Helix", "Sheet", "Coil"]),
        ("Solvent accessibility", lambda a: "Buried" if a in buried else ("Exposed" if a in exposed else "Intermediate"), ["Buried", "Intermediate", "Exposed"]),
    ]

    color_triplet = ["yellow", "orange", "red"]
    panel_paths = []
    for title, fn, order in attribs:
        cats = [fn(a) for a in AA]
        color_map = {cat: color_triplet[i] for i, cat in enumerate(order)}
        colors = [color_map[c] for c in cats]
        plt.figure(figsize=(12, 1.4))
        plt.bar(range(len(AA)), np.ones(len(AA)), color=colors, tick_label=AA)
        plt.yticks([])
        plt.title(title, loc="left")
        out = outdir / f"panel_bar_{title.replace(' ', '_').lower()}.png"
        plt.savefig(out, bbox_inches="tight")
        plt.close()
        panel_paths.append(out)

    plt.figure(figsize=(12, 1.2))
    plt.axis("off")
    plt.text(0.5, 0.5, "Biochemical property panels", ha="center", va="center")
    title_path = outdir / "panel_bar_title.png"
    plt.savefig(title_path, bbox_inches="tight")
    plt.close()

    imgs = [Image.open(title_path)] + [Image.open(p) for p in panel_paths]
    width = max(im.width for im in imgs)
    pad = 8
    height = sum(im.height for im in imgs) + pad * (len(imgs) + 1)
    canvas = Image.new("RGB", (width + 20, height), "white")
    y = pad
    for im in imgs:
        x = (canvas.width - im.width) // 2
        canvas.paste(im, (x, y))
        y += im.height + pad
    canvas.save(outdir / "panel_bar_stacked.png")


def fingerprint(amp: list[tuple[str, str]], outdir: Path, accession: str | None, index: int) -> None:
    if not amp:
        return
    seq_idx = None
    if accession is not None:
        for i, (header, _) in enumerate(amp):
            tok = first_token_id(header)
            if accession in (tok, header) or accession in header:
                seq_idx = i
                break
    if seq_idx is None:
        seq_idx = max(0, min(index, len(amp) - 1))
    header, seq = amp[seq_idx]
    seq = seq.upper()

    properties = {
        "hydrophobicity": {"polar": "RKEDQN", "neutral": "GASTPHY", "hydrophobic": "CLVIMFW"},
        "vdw_volume": {"small": "GASTPDC", "medium": "NVEQIL", "large": "MHKFRYW"},
        "polarity": {"polar": "EDQNKR", "neutral": "GASTPHY", "nonpolar": "CLVIMFW"},
        "polarizability": {"low": "GASDT", "medium": "CPNVEQIL", "high": "MHKFRYW"},
        "charge": {"negative": "DE", "neutral": "ACFGHILMNPQSTVWY", "positive": "KR"},
        "secondary_structure": {"helix": "EALMQKRH", "sheet": "VIYCWT", "coil": "GNPSD"},
        "solvent_accessibility": {"buried": "ALFCGIVW", "intermediate": "MPSTHY", "exposed": "DEKNQR"},
    }
    orders = {
        "hydrophobicity": ["polar", "neutral", "hydrophobic"],
        "vdw_volume": ["small", "medium", "large"],
        "polarity": ["polar", "neutral", "nonpolar"],
        "polarizability": ["low", "medium", "high"],
        "charge": ["negative", "neutral", "positive"],
        "secondary_structure": ["helix", "sheet", "coil"],
        "solvent_accessibility": ["buried", "intermediate", "exposed"],
    }
    color_triplet = ["yellow", "orange", "red"]
    panel_paths = []
    rows = []
    for key in ["hydrophobicity", "vdw_volume", "polarity", "polarizability", "charge", "secondary_structure", "solvent_accessibility"]:
        prop = properties[key]
        order = orders[key]
        color_map = {cat: color_triplet[i] for i, cat in enumerate(order)}
        categories, colors = [], []
        for aa in seq:
            cat = next((c for c, members in prop.items() if aa in members), "unknown")
            categories.append(cat)
            colors.append(color_map.get(cat, "lightgray"))
        rows.append(pd.DataFrame({"position": list(range(1, len(seq) + 1)), "aa": list(seq), "attribute": key, "category": categories}))
        plt.figure(figsize=(max(6, len(seq) * 0.25), 1.6))
        plt.bar(range(len(seq)), np.ones(len(seq)), color=colors, tick_label=list(seq))
        plt.xticks(range(len(seq)), list(seq), fontsize=8)
        plt.yticks([])
        plt.title(key.replace("_", " ").capitalize(), loc="left")
        out = outdir / f"fingerprint_{key}.png"
        plt.savefig(out, bbox_inches="tight")
        plt.close()
        panel_paths.append(out)
    pd.concat(rows, ignore_index=True).to_csv(outdir / "fingerprint_categories_by_position.csv", index=False)

    title_path = outdir / "fingerprint_title.png"
    plt.figure(figsize=(12, 1.2))
    plt.axis("off")
    plt.text(0.5, 0.5, f"Peptide fingerprint: {header}", ha="center", va="center")
    plt.savefig(title_path, bbox_inches="tight")
    plt.close()

    imgs = [Image.open(title_path)] + [Image.open(p) for p in panel_paths]
    width = max(im.width for im in imgs)
    pad = 8
    height = sum(im.height for im in imgs) + pad * (len(imgs) + 1)
    canvas = Image.new("RGB", (width + 20, height), "white")
    y = pad
    for im in imgs:
        x = (canvas.width - im.width) // 2
        canvas.paste(im, (x, y))
        y += im.height + pad
    canvas.save(outdir / "fingerprint_stacked.png")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.dataset_csv:
        amp, nonamp = load_from_dataset_csv(Path(args.dataset_csv))
    else:
        if not args.amp or not args.nonamp:
            raise ValueError("--amp and --nonamp are both required when --dataset-csv is not used")
        amp = parse_fasta(Path(args.amp))
        nonamp = parse_fasta(Path(args.nonamp))

    build_qc_tables(amp, nonamp, outdir)
    plot_counts(amp, nonamp, outdir)
    plot_lengths(amp, nonamp, outdir)
    save_length_bins(amp, nonamp, outdir)
    aa_freqs(amp, nonamp, outdir)
    property_panels(outdir)
    fingerprint(amp, outdir, accession=args.fingerprint_accession, index=args.fingerprint_index)
    print(f"Wrote outputs to {outdir}")


if __name__ == "__main__":
    main()
