# Output overview

Use `README.md` for the workflow order and command examples.

Bundled inputs:
- `data/dataset.csv` - standard downstream input table for steps 03-06
- `data/balanced_11000.csv` - compact convenience input for downstream runs

Main downstream products:
- Step 03: `release_cluster_split.csv`, `train.csv`, `val.csv`, `test.csv`, and the corresponding FASTA files
- Step 04: `mmseqs_all_pairs.csv`, `mmseqs_cross_split_leaks.csv`, and `mmseqs_leak_summary.txt`
- Step 05: a CTD-like feature table
- Step 06: summary tables and figure files
