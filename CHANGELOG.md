# Changelog

## v1.1 - 2026-09-27

- aligned the public repository to the frozen Zenodo v1.1 canonical dataset;
- replaced old convenience data with the canonical 11,000-record CSV;
- preserved the supplied dataset-construction/audit workflow sources;
- aligned Step 03 default split seed to 42 and made operational length-bin parsing robust to canonical display labels;
- replaced the previous CTD-only public Step 05 role with the current manuscript length-stratified EDA reproduction;
- added the authoritative final 13-model / 10-seed benchmark notebook and exact Cell 32 v9 trainer;
- added a portable benchmark launcher without rewriting model logic;
- documented required third-party datasets and published-method environments;
- linked the public dataset DOI `10.5281/zenodo.22994107`;
- intentionally did not bundle static Tables 2-4, old benchmark outputs, reviewer-response archives, or third-party benchmark sequences.

## Packaging audit corrections - 2026-10-02

- Removed ANSI color escapes from saved notebook tracebacks without changing any source cell, numerical output, error type or error message.
- Original notebook SHA-256: `62625fc0047a74083ad51ff82167cb0119c9d1679da333379c21fa27d985ca6d`. Output-format-cleaned notebook SHA-256: `97cbc1e680fca5118af5e758d09f3aaa25ea4f5c5c36390105b7ef908c28d4ac`.
- Recomputed the package manifest and SHA-256 checksums from the actual file bytes. The manifest inventories payload; SHA256SUMS.txt also covers MANIFEST.csv.

## v1.1 packaging synchronization - 2026-10-03

- Unified local, GitHub main/tag/release and Zenodo package contents under one official archive filename.
- Retained the public dataset QA and integrity provenance in validation/ and recorded current Zenodo metadata.
- Regenerated manifest and SHA256 checksums from final payload bytes.
- Kept the canonical CSV, direct split files, scientific source cells and meaningful saved outputs unchanged.
