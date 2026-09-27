# Pinned published-method source code

The authoritative benchmark notebook audits published implementations before using them. Do not substitute newer commits without re-auditing the benchmark.

## Methods reported in the current 13-model manuscript

### AMPScannerV2

- Repository recorded by the authoritative notebook: `dan-veltri/amp-scanner-v2`
- Pinned commit: `933052e2365631fe93098892120ee535e0ba381a`
- Expected local directory:
  `external_methods/amp-scanner-v2_933052e2365631fe93098892120ee535e0ba381a/`

### amPEPpy v1.0

- Repository recorded by the authoritative notebook: `tlawrence3/amPEPpy`
- Tag: `v1.0`
- Pinned commit: `aa1f694c6cb3d09b16bc9378bed77f59e4f1e780`
- Expected local directory:
  `external_methods/amPEPpy_v1.0_aa1f694c6cb3d09b16bc9378bed77f59e4f1e780/`

### AmPEP

The notebook uses the author-released AmPEP MATLAB source as the scientific reference and a faithful Python reimplementation of the 105 D_F representation plus the released RF recipe for the benchmark adapter.

- Official release channel recorded in the notebook: SourceForge project `axpep`, directory `AmPEP_MATLAB_code`
- Pinned public mirror recorded by the notebook: `lassebuur/bacteriocin_classifier`
- Pinned mirror commit: `1fedfe1c3203fdd5a417ba1c2e92123478f2cfc0`
- Mirror subdirectory: `code_modules/ampep/ampep-matlab-code`
- Expected local source directory: `external_methods/AmPEP_author_release/`

The authoritative notebook verifies required source files/hashes before opening the benchmark training gate.

## Candidate methods not reported in the final 13-model paper

The historical authoritative notebook also audits AmpGram, ampir, AI4AMP, and AMPlify as model candidates before the final scope lock. They are not part of the 13 models reported in the current manuscript. The original notebook is retained unmodified for provenance; the final scope lock is visible in that source.
