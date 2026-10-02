# Reviewer run status

The public bundle was validated for source identity and runnable preparation steps.

## Verified in the packaging environment

- frozen canonical dataset SHA-256 and row/class/split counts;
- 11,000 unique sequences and sample IDs;
- 10,666 components and zero components spanning frozen splits;
- Python compilation of all public `.py` files;
- `--help` execution for workflow scripts 01-04;
- execution of `validate_dataset.py`;
- execution of `05_data_visualization.py` on the canonical dataset;
- exact reproduction of the current manuscript length-bin counts;
- authoritative benchmark notebook SHA-256;
- exact final trainer SHA-256;
- `06_benchmark_reproduction.py --prepare-only`;
- Python compilation of code cells in the prepared notebook after IPython transformation.

## Not claimed

A fresh full 13-model x 10-seed training rerun was not completed during packaging. That run requires third-party comparison datasets, published-method source trees, incompatible legacy method environments, the main GPU environment, and substantial compute time. The repository therefore does not claim an unperformed end-to-end clean-machine run.
