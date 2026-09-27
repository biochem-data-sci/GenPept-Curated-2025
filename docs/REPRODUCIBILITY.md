# Reproducibility guide

## Exact paper reproduction

Use the frozen v1.1 canonical dataset, not a new NCBI download.

1. Validate the dataset:

```bash
python validate_dataset.py
```

2. Reproduce the current EDA:

```bash
python 05_data_visualization.py --dataset-csv data/GenPept_Curated_2025_primary_split_v1.1.csv --outdir outputs/eda
```

3. Place the external comparison datasets exactly as described in `external_data/README.md`.

4. Place the published-method source trees exactly as described in `external_methods/README.md`.

5. Recreate the main and legacy environments from `benchmark/ENVIRONMENT.md` and the audited requirement files under `benchmark/`.

6. Prepare the authoritative benchmark notebook:

```bash
python 06_benchmark_reproduction.py --prepare-only
```

7. Execute the prepared notebook only after the required external data/method environments are available:

```bash
python 06_benchmark_reproduction.py --execute
```

## Dataset-construction audit path

Scripts `01`–`04` are retained for procedural auditing of retrieval, curation, clustering, and post-split high-identity checks. They are **not** the authoritative paper-reproduction path because NCBI is mutable and the supplied Step 02 source does not contain the complete frozen IPG representative-selection provenance.

The exact released dataset identity is instead fixed by:

- version 1.1
- DOI `10.5281/zenodo.22994107`
- canonical SHA-256 `421e55265e1462052f633c961f8d8cc20ce5e510284ffc36afc6e29a40b79b6c`
- 11,000 records
- 5,500 AMP / 5,500 non-AMP
- 7,700 / 990 / 2,310 train/validation/test records
- 10,666 high-identity components
- zero components spanning multiple frozen splits
