# GenPept-Curated-2025

Public dataset and reproducibility code for the manuscript **“GenPept-Curated-2025: A New Machine Learning Benchmark for Antimicrobial Peptide Prediction.”**

This repository update is aligned to the frozen public dataset **v1.1** and the current manuscript state supplied for review. The exact dataset used for the paper is the frozen canonical CSV below; live NCBI retrieval is provided only as a procedural audit path because upstream records can change over time.

## Frozen dataset used by the paper

- Version: **1.1**
- Zenodo DOI: **10.5281/zenodo.22994107**
- Canonical file: `data/GenPept_Curated_2025_primary_split_v1.1.csv`
- SHA-256: `421e55265e1462052f633c961f8d8cc20ce5e510284ffc36afc6e29a40b79b6c`
- Records: **11,000**
- AMP-labeled: **5,500**
- non-AMP-labeled: **5,500**
- Sequence length: **10–200 aa**
- Frozen split: **7,700 train / 990 validation / 2,310 test**
- Unique sequences: **11,000**
- High-identity components in the frozen release: **10,666**
- Components spanning more than one split: **0**

The `non-AMP` label is annotation-negative/unlabeled under the study's operational rule; it is **not** an experimentally confirmed inactive class.

## What to run to reproduce the paper

### 1. Validate the exact frozen dataset

```bash
python validate_dataset.py
```

### 2. Reproduce the length-stratified EDA

```bash
pip install -r requirements.txt
python 05_data_visualization.py \
  --dataset-csv data/GenPept_Curated_2025_primary_split_v1.1.csv \
  --outdir outputs/eda
```

This generates the length-bin summary and the current three-panel class-composition figure from the canonical CSV. No precomputed manuscript result table is required.

### 3. Prepare the final benchmark source

The authoritative final benchmark notebook supplied by the authors is preserved byte-for-byte at:

```text
benchmark/BenchMark17model_AUTHORITATIVE.ipynb
```

SHA-256:

```text
f9e8185de48cebba5c5eb3d0b7748b358d0dfb1ea9d078584bd7f97a58901141
```

The exact production trainer emitted by the final notebook is also preserved:

```text
benchmark/Cell32_Benchmark17_Production_Trainer_v9.py
```

SHA-256:

```text
fadd6d8f398b0c7fabcc4c73c4133484cdf95986705c3d5e8c54d6b805bf6a6b
```

Prepare a portable copy of the notebook:

```bash
python 06_benchmark_reproduction.py --prepare-only
```

After the external datasets, published-method source trees, and required environments are installed, execute it with:

```bash
python 06_benchmark_reproduction.py --execute
```

The current manuscript benchmark contains **13 models** and **10 seeds (100–109)** and performs reciprocal evaluation with **SSFGM-BD1**, **AMPlify-balanced**, and **SSFGM-BD3** after exact-sequence overlap removal. See `benchmark/README.md`, `benchmark/ENVIRONMENT.md`, `external_data/README.md`, and `external_methods/README.md`.

## Repository files

### Primary paper-reproduction path

- `data/GenPept_Curated_2025_primary_split_v1.1.csv` — exact frozen dataset used by the paper.
- `05_data_visualization.py` — current manuscript EDA/length-bin reproduction.
- `06_benchmark_reproduction.py` — portability/execution launcher for the authoritative final benchmark notebook.
- `benchmark/BenchMark17model_AUTHORITATIVE.ipynb` — original final benchmark source supplied by the authors.
- `benchmark/Cell32_Benchmark17_Production_Trainer_v9.py` — exact final production trainer emitted by the notebook.

### Dataset-construction / audit helpers

- `01_data_retrieval_genpept.py` — preserved GenPept/NCBI query runner using the stored query specification.
- `02_build_balanced_dataset.py` — preserved curation helper from the supplied data-pipeline archive.
- `03_cluster_split_cdhit.py` — supplied CD-HIT cluster-intact splitter, with the public-update default seed aligned to **42** and operational length bins recomputed from numeric length.
- `04_check_cross_split_leakage_mmseqs2.py` — supplied MMseqs2 post-split high-identity relation checker.

These four scripts are useful for procedural auditing, but **the exact paper reproduction must use the frozen v1.1 CSV rather than a new live NCBI download**. In particular, the supplied Step 02 source does not by itself preserve the complete IPG representative-selection provenance needed to claim byte-identical reconstruction of the frozen release. This limitation is stated explicitly rather than silently reconstructed.

## Main environment recorded by the final benchmark source

The preserved environment audit records:

- Python 3.11.16
- TensorFlow 2.21.0
- PyTorch 2.11.0+cu128 / CUDA 12.8
- scikit-learn 1.9.0
- LightGBM 4.7.0
- fair-esm 2.0.0
- Linux x86_64 under WSL2
- NVIDIA GeForce RTX 3060, 12 GB VRAM
- NVIDIA driver 610.88

Legacy published methods use separate environments. Exact versions that were present in the supplied source are documented under `benchmark/`.

## Third-party benchmark data

Third-party SSFGM and AMPlify data are not redistributed in this repository because redistribution permission was not established by the supplied project files. The exact expected local layout and filenames are documented in `external_data/README.md`.

## Data dictionary, provenance, and rights

- `docs/DATA_DICTIONARY.csv`
- `docs/SOURCE_AND_RIGHTS.md`
- `provenance/query_specification_from_preserved_rules.csv`

No repository-wide software license is asserted here because a code license was not established in the supplied source materials. The creator-owned dataset curation/documentation rights statement is in `docs/SOURCE_AND_RIGHTS.md`.

## Figures

- [Figure 1: Curation pipeline](figures/Fig1.pdf).
- [Figure 2: Length-dependent dataset structure and class distribution](figures/Fig2.pdf).

## Citation

Dataset DOI: **https://doi.org/10.5281/zenodo.22994107**

A machine-readable citation file is provided as `CITATION.cff`.
