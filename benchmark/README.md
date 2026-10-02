# Final benchmark source

The authoritative benchmark source for the current manuscript is:

`BenchMark17model_AUTHORITATIVE.ipynb`

SHA-256:

`97cbc1e680fca5118af5e758d09f3aaa25ea4f5c5c36390105b7ef908c28d4ac`

The notebook contains the final 13-model scope, 10 random seeds (100-109), within-dataset evaluation, and reciprocal cross-dataset evaluation with SSFGM-BD1, AMPlify-balanced, and SSFGM-BD3.

`Cell32_Benchmark17_Production_Trainer_v9.py` is the exact final production trainer emitted from the notebook and has SHA-256:

`bc2276c2377e0c0730892c78412e307949547b4771b1eca10073c6ff125e5b6d`

Use `../06_benchmark_reproduction.py --prepare-only` to create a repository-portable copy. The launcher only removes superseded Cell 32 v4/v8 cells, clears stale outputs, and rewrites author-machine absolute paths. It does not redefine the scientific models, seeds, metrics, threshold rules, or evaluation directions.

A full rerun requires the external datasets and published-method source trees documented under `../external_data/` and `../external_methods/`, plus the separate environments described in `ENVIRONMENT.md`.
