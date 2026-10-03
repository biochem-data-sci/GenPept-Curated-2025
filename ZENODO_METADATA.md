# Zenodo record metadata - GenPept-Curated-2025 v1.1

- Resource type: Dataset
- Title: GenPept-Curated-2025: An annotation-derived, high-identity-controlled benchmark for antimicrobial peptide annotation prediction
- Version: 1.1
- Publication date: 2026-09-27
- DOI: 10.5281/zenodo.22994107
- Public record: https://zenodo.org/records/22994107
- Visibility: Public

## Creators

- Pham, Huynh Trong
- Nguyen-Vo, Thanh-Hoang
- Huynh, Bao

## Canonical dataset

The canonical file is `data/GenPept_Curated_2025_primary_split_v1.1.csv`, SHA256 `bdfe03fbabbd0c7d6ff95ff42905680bb86ad0056e24cf964a33598aefb407a2`. It contains 11,000 unique sequences: 5,500 AMP annotation-positive and 5,500 annotation-negative/unlabeled comparison records. The latter are not experimentally confirmed inactive peptides. The balanced composition is a controlled benchmark design, not a natural prevalence estimate.

The frozen partitions contain 7,700 train, 990 validation (`val`) and 2,310 test records. The 10,666 operational high-identity components do not cross these partitions under the preserved reported procedure; this does not establish universal remote-homology or biological-family independence.

## Release package

`GenPept-Curated-2025_ZENODO_DATASET_v1.1_Final.zip` is the official synchronized package for local distribution, GitHub Release v1.1 and this Zenodo version. GitHub main and tag v1.1 contain the same relative file paths and bytes. It includes the canonical CSV and direct splits, construction and benchmark reproduction sources, manuscript figures, environment documentation, source rights, preserved QA records, MANIFEST.csv and SHA256SUMS.txt.

Run `python validate_dataset.py` from the extracted directory. Preserved similarity-screen summaries and benchmark outputs are historical scientific evidence; packaging synchronization does not represent a new CD-HIT, MMseqs2 or model-training run. QA source hashes explicitly marked historical identify provenance and are not alternate current datasets.

The dictionary and source-rights notice are `docs/DATA_DICTIONARY.csv` and `docs/SOURCE_AND_RIGHTS.md`. QA and integrity records are in `validation/`. Cite this dataset DOI; any related manuscript DOI is a separate identifier.

## License and source rights

Creative Commons Attribution 4.0 International (CC BY 4.0; SPDX CC-BY-4.0) applies to creator-owned curation, operational labels, split assignments, documentation, validation metadata and compilation only to the extent the creators hold the relevant rights. It does not purport to relicense underlying third-party NCBI/GenBank rights. See `docs/SOURCE_AND_RIGHTS.md` for the complete source-rights notice and source policy links.
