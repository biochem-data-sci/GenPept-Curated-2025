# Audited benchmark environments

Only versions actually recorded by the supplied authoritative notebook are stated here.

## Main benchmark environment

Recorded by the notebook environment audit:

- Python 3.11.16
- PyTorch 2.11.0+cu128
- Torch CUDA 12.8
- TensorFlow 2.21.0 (CUDA build; one GPU detected)
- scikit-learn 1.9.0
- LightGBM 4.7.0
- fair-esm 2.0.0
- NVIDIA driver 610.88
- NVIDIA GeForce RTX 3060, 12 GB VRAM
- Linux x86_64 under WSL2

The notebook says that a full `pip freeze` was written on the author machine, but that file was not among the supplied files used to assemble this GitHub update. Therefore this repository does not invent exact versions for packages not shown in the preserved environment audit.

## AMPScannerV2 legacy environment

The notebook's successful adapter self-test recorded:

- Python 3.6.13
- TensorFlow 1.2.1
- Keras 2.0.6
- NumPy 1.16.0
- h5py 2.6.0
- Biopython 1.69
- scikit-learn 0.20.0
- Theano 0.9.0 (installed dependency; TensorFlow is the runtime backend)
- backports.weakref 1.0rc1

## amPEPpy v1.0 legacy environment

The notebook's successful source/adapter self-test recorded:

- Python 3.8.3
- scikit-learn 0.23.1
- NumPy 1.18.5
- SciPy 1.5.0
- pandas 1.0.5
- Biopython 1.77

## AmPEP

The notebook audits the released MATLAB source as the reference implementation and validates a faithful Python adapter for the 105 D_F representation and RF recipe in the benchmark workflow. No separate MATLAB runtime is required for the final Python benchmark adapter, but the released source files are checked as provenance.
