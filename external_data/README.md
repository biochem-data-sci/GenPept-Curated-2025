# External benchmark datasets

The current manuscript performs reciprocal cross-dataset evaluation with:

1. SSFGM-BD1
2. AMPlify-balanced
3. SSFGM-BD3

These are third-party datasets. They are **not redistributed in this GitHub update** because the supplied project files do not establish redistribution permission for them.

## Expected local layout for the benchmark notebook

The portable launcher rewrites the author's machine-specific paths to the following repository-local locations:

```text
external_data/
├── SSFGM-Model-main/
│   └── Data/
│       ├── Benchmark dataset 1/
│       │   ├── amp_train2.fasta
│       │   ├── amp_eval2.fasta
│       │   ├── amp_test2.fasta
│       │   ├── non_amp_train2.fasta
│       │   ├── non_amp_eval2.fasta
│       │   └── non_amp_test2.fasta
│       ├── Benchmark dataset 2/        # historical discovery input; final paper does not retain BD2
│       └── Benchmark dataset 3/
│           ├── AMP/
│           │   ├── train.fasta
│           │   ├── valid.fasta
│           │   └── test.fasta
│           └── nonAMP/
│               ├── train.fasta
│               ├── valid.fasta
│               └── test.fasta
└── AMPlify/
    ├── AMPlify_AMP_train_common.fa
    ├── AMPlify_non_AMP_train_balanced.fa
    ├── AMPlify_AMP_test_common.fa
    └── AMPlify_non_AMP_test_balanced.fa
```

The supplied `Dataset.zip` contains the SSFGM source tree. Its original BD1 folder is named `Benchmark dataset1` (without a space), whereas the authoritative benchmark notebook discovers `Benchmark dataset 1`. For a clean reviewer run, rename/copy that folder to the exact expected name shown above; do not modify the FASTA contents.

The four AMPlify FASTA filenames above are hard-locked by the authoritative benchmark notebook. The underlying AMPlify FASTA files were not present in the supplied archives used to assemble this update, so this bundle does not invent checksums or download URLs for them.

The benchmark notebook creates its AMPlify validation partition deterministically from the original training data according to its locked protocol. Do not pre-create a different validation set.
