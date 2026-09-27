# Final benchmark source

The authoritative benchmark source for the current manuscript is:

`BenchMark17model_AUTHORITATIVE.ipynb`

SHA-256:

`f9e8185de48cebba5c5eb3d0b7748b358d0dfb1ea9d078584bd7f97a58901141`

The notebook contains the final 13-model scope, 10 random seeds (100–109), within-dataset evaluation, and reciprocal cross-dataset evaluation with SSFGM-BD1, AMPlify-balanced, and SSFGM-BD3.

`Cell32_Benchmark17_Production_Trainer_v9.py` is the exact final production trainer emitted from the notebook and has SHA-256:

`fadd6d8f398b0c7fabcc4c73c4133484cdf95986705c3d5e8c54d6b805bf6a6b`

Use `../06_benchmark_reproduction.py --prepare-only` to create a repository-portable copy. The launcher only removes superseded Cell 32 v4/v8 cells, clears stale outputs, and rewrites author-machine absolute paths. It does not redefine the scientific models, seeds, metrics, threshold rules, or evaluation directions.

A full rerun requires the external datasets and published-method source trees documented under `../external_data/` and `../external_methods/`, plus the separate environments described in `ENVIRONMENT.md`.
