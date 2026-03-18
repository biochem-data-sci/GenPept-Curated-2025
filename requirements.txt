# GenPept-Curated-2025 code package

This package is organized for journal submission and public release. It contains the processing scripts, the visualization script, and prepared input tables for the downstream reproducibility workflow.

## Package contents

- `01_data_retrieval_genpept.py` - retrieve raw GenPept records from NCBI Protein using the study scope.
- `02_build_balanced_dataset.py` - apply sequence QC, annotation-based labeling, precursor flagging, IPG-aware deduplication, and balanced dataset construction.
- `03_cluster_split_cdhit.py` - create the official cluster-intact train/validation/test split with CD-HIT.
- `04_check_cross_split_leakage_mmseqs2.py` - run the post-split all-vs-all leakage check with MMseqs2.
- `05_extract_ctd_features.py` - extract CTD-like composition and distribution features.
- `06_data_visualization.py` - reproduce summary tables and figures.
- `data/dataset.csv` - standardized downstream input table.
- `data/balanced_11000.csv` - compact convenience table for quick downstream runs.

## Standardized downstream input

`data/dataset.csv` is the standard downstream input table for steps 03-06.

Required core fields used by the downstream scripts:
- `accession_version`
- `sequence`
- `length`
- `label`
- `length_bin`
- `precursor_flag`

Optional retained metadata fields:
- `tax_id`
- `taxon_name`
- `description`
- `cds_products`
- `cds_genes`
- `ipg_id`
- `refseq_accession`

## Workflow order

### Step 1 - `01_data_retrieval_genpept.py`
Retrieve raw GenPept records from NCBI Protein using the study scope:
- source: GenPept / GenBank proteins
- taxa: Bacteria, Archaea, Fungi
- date window: 2025-04-01 to 2025-07-31
- sequence length window: 10-200 aa

Main outputs:
- `01_raw_records.csv`
- `01_raw_records.fasta`
- `01_raw_taxon_summary.csv`
- `01_toc_left_values.csv`

### Step 2 - `02_build_balanced_dataset.py`
Apply the curation workflow to the raw records:
- retain only sequences containing the 20 standard amino acids within the 10-200 aa window
- remove low-quality annotations
- flag precursor records and export them separately
- assign AMP / non-AMP labels from annotation rules
- deduplicate by IPG when available
- construct the balanced 11,000-sequence dataset

Main outputs:
- `dataset.csv`
- `dataset.fasta`
- `amp_sequences.fasta`
- `nonamp_sequences.fasta`
- `precursor_flagged.csv`
- `metadata.csv`
- `accession_manifest.tsv`
- `ipg_mapping.tsv`
- `dataset_summary_by_bin.csv`
- `attrition_summary.csv`

### Step 3 - `03_cluster_split_cdhit.py`
Create the official cluster-intact train/validation/test split.

Protocol encoded in the script:
- CD-HIT identity threshold: `0.90`
- coverage threshold on the shorter sequence: `0.80`
- split ratio: `70 : 9 : 21`
- stratification: `label x length_bin`
- random seed: `42`

For the balanced 11,000-sequence release, the target split sizes are:
- train: `7700` total (`3850 AMP`, `3850 non-AMP`)
- validation: `990` total (`495 AMP`, `495 non-AMP`)
- test: `2310` total (`1155 AMP`, `1155 non-AMP`)

Main outputs:
- `all_11000_for_cdhit.fasta`
- `cdhit_c0.9_n5.clstr`
- `release_cluster_split.csv`
- `cluster_assignment.csv`
- `train.csv`, `val.csv`, `test.csv`
- `train.fasta`, `val.fasta`, `test.fasta`
- `split_summary_by_bin.csv`

### Step 4 - `04_check_cross_split_leakage_mmseqs2.py`
Run the post-split all-vs-all leakage check.

Protocol encoded in the script:
- identity threshold: `0.90`
- coverage threshold on the shorter sequence: `0.80`
- leakage definition: any cross-split pair meeting both thresholds

Main outputs:
- `mmseqs_all_pairs.csv`
- `mmseqs_cross_split_leaks.csv`
- `mmseqs_leak_summary.txt`

### Step 5 - `05_extract_ctd_features.py`
Extract CTD-like composition and distribution features from either:
- a dataset CSV, or
- a FASTA file.

### Step 6 - `06_data_visualization.py`
Reproduce summary tables and figures from either:
- a dataset CSV, or
- paired AMP / non-AMP FASTA files.

## Recommended commands

### Full reconstruction mode
Use this path when regenerating the balanced dataset from raw GenPept retrieval.

#### Step 1 - retrieve raw GenPept records
```bash
python 01_data_retrieval_genpept.py   --email your_email@example.com   --outdir outputs/01_retrieval
```

#### Step 2 - curate, deduplicate, label, and build the balanced dataset
```bash
python 02_build_balanced_dataset.py   --input-csv outputs/01_retrieval/01_raw_records.csv   --outdir outputs/02_dataset   --dedup-mode ipg   --fallback-to-exact
```

#### Step 3 - create cluster-based train/validation/test splits
```bash
python 03_cluster_split_cdhit.py   --dataset-csv outputs/02_dataset/dataset.csv   --outdir outputs/03_split
```

#### Step 4 - verify cross-split leakage with MMseqs2
```bash
python 04_check_cross_split_leakage_mmseqs2.py   --release-csv outputs/03_split/release_cluster_split.csv   --outdir outputs/04_mmseqs
```

#### Step 5 - feature extraction
```bash
python 05_extract_ctd_features.py   --dataset-csv outputs/02_dataset/dataset.csv   --out outputs/05_features/dataset_ctd_features.csv
```

#### Step 6 - summary tables and figures
```bash
python 06_data_visualization.py   --dataset-csv outputs/02_dataset/dataset.csv   --outdir outputs/06_figures
```

### Quick-start mode from the prepared dataset
Use this path when reproducing the official downstream workflow from the bundled input table.

#### Step 3 - create train/validation/test splits
```bash
python 03_cluster_split_cdhit.py   --dataset-csv data/dataset.csv   --outdir outputs/03_split
```

#### Step 4 - leakage check
```bash
python 04_check_cross_split_leakage_mmseqs2.py   --release-csv outputs/03_split/release_cluster_split.csv   --outdir outputs/04_mmseqs
```

#### Step 5 - CTD-like features
```bash
python 05_extract_ctd_features.py   --dataset-csv data/dataset.csv   --out outputs/05_features/dataset_ctd_features.csv
```

#### Step 6 - summary tables and figures
```bash
python 06_data_visualization.py   --dataset-csv data/dataset.csv   --outdir outputs/06_figures
```

## Dependencies

Install Python packages:
```bash
pip install -r requirements.txt
```

External tools required:
- `cd-hit` for step 03
- `mmseqs` for step 04

On Windows, steps 03 and 04 can call the binaries through WSL if CD-HIT and MMseqs2 are installed inside WSL.
