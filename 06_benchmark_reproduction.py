#!/usr/bin/env python3
"""Prepare or execute a portable copy of the authoritative final benchmark notebook.

Scientific source of truth:
    benchmark/BenchMark17model_AUTHORITATIVE.ipynb

This launcher does NOT rewrite model architectures, hyperparameters, seeds, metrics,
threshold policies, or evaluation logic. It performs only two operational changes:
1) removes superseded Cell 32 v4/v8 notebook cells that were replaced by Cell 32 v9;
2) rewrites machine-specific absolute paths to paths under the cloned repository.

A full reviewer rerun still requires the external benchmark datasets and the pinned
published-method source trees/environments documented under external_data/,
external_methods/, and docs/REPRODUCIBILITY.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

AUTHORITATIVE_NOTEBOOK_SHA256 = "f9e8185de48cebba5c5eb3d0b7748b358d0dfb1ea9d078584bd7f97a58901141"
SUPERSEDED_CELL_INDICES = {31, 32}  # Cell 32 v4 and v8; Cell 32 v9 is index 33.
ORIGINAL_PROJECT_ROOT = "/mnt/d/SUABAI_GenPept-Curated-2025_ 12.8.2026"
ORIGINAL_WORK_ROOT = "/home/pc/genpept_benchmark17_work"
ORIGINAL_SOURCE_ROOT = "/home/pc/genpept_sources"
ORIGINAL_MASTER_PYTHON = "/home/pc/miniconda3/envs/genpept_benchmark10_master/bin/python"
ORIGINAL_CONDA = "/home/pc/miniconda3/bin/conda"
ORIGINAL_ENVS_ROOT = "/home/pc/miniconda3/envs/"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-only", action="store_true", help="Create the portable notebook but do not execute it.")
    parser.add_argument("--execute", action="store_true", help="Execute the prepared notebook with Jupyter nbconvert.")
    parser.add_argument("--project-root", default=None, help="Repository root. Default: directory containing this script.")
    parser.add_argument("--work-root", default=None, help="Heavy benchmark work directory.")
    parser.add_argument("--source-root", default=None, help="Pinned published-method source directory.")
    parser.add_argument("--conda-exe", default=None, help="Path to conda executable for legacy method environments.")
    parser.add_argument("--output-notebook", default="benchmark/BenchMark13model_PORTABLE_PREPARED.ipynb")
    return parser.parse_args()


def replace_exact_block(text: str, pattern: str, replacement: str, label: str) -> str:
    new_text, count = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"Expected exactly one {label} block; found {count}.")
    return new_text


def patch_source(
    source: str,
    project_root: Path,
    work_root: Path,
    source_root: Path,
    conda_exe: Path,
    envs_root: Path,
) -> str:
    # Pure machine-path substitutions.
    source = source.replace(ORIGINAL_PROJECT_ROOT, str(project_root))
    source = source.replace(ORIGINAL_WORK_ROOT, str(work_root))
    source = source.replace(ORIGINAL_SOURCE_ROOT, str(source_root))
    source = source.replace(ORIGINAL_MASTER_PYTHON, sys.executable)
    source = source.replace(ORIGINAL_CONDA, str(conda_exe))
    source = source.replace(ORIGINAL_ENVS_ROOT, str(envs_root) + os.sep)
    source = source.replace('"BenchMark17model.ipynb"', '"benchmark/BenchMark17model_AUTHORITATIVE.ipynb"')

    # Frozen v1.1 canonical dataset now lives directly under data/ in the repository.
    genpept_pattern = (
        r'GENPEPT_PRIMARY = \(\s*PROJECT_ROOT\s*/ "Bộ dữ liệu trả lời_FINAL_LOCKED"\s*'
        r'/ "01_CANONICAL_DATASET"\s*/ "GenPept_Curated_2025_primary_split_v1\.1\.csv"\s*\)'
    )
    genpept_replacement = (
        'GENPEPT_PRIMARY = (\n'
        '    PROJECT_ROOT\n'
        '    / "data"\n'
        '    / "GenPept_Curated_2025_primary_split_v1.1.csv"\n'
        ')'
    )
    if re.search(genpept_pattern, source, flags=re.MULTILINE):
        source = replace_exact_block(source, genpept_pattern, genpept_replacement, "GenPept canonical path")

    # Expected SSFGM extraction layout for the uploaded Dataset.zip source material.
    external_pattern = (
        r'EXTERNAL_DATA_ROOT = \(\s*PROJECT_ROOT\s*/ "Bộ dữ liệu tra lời"\s*'
        r'/ "02_BENCHMARK"\s*/ "FROM_ZIP"\s*/ "Dataset"\s*/ "Dataset"\s*'
        r'/ "Sota"\s*/ "SSFGM-Model-main"\s*/ "SSFGM-Model-main"\s*/ "Data"\s*\)'
    )
    external_replacement = (
        'EXTERNAL_DATA_ROOT = (\n'
        '    PROJECT_ROOT\n'
        '    / "external_data"\n'
        '    / "SSFGM-Model-main"\n'
        '    / "Data"\n'
        ')'
    )
    if re.search(external_pattern, source, flags=re.MULTILINE):
        source = replace_exact_block(source, external_pattern, external_replacement, "SSFGM external-data path")

    amplify_pattern = r'AMPLIFY_DIR = \(\s*PROJECT_ROOT\s*/ "AMPlify"\s*\)'
    amplify_replacement = (
        'AMPLIFY_DIR = (\n'
        '    PROJECT_ROOT\n'
        '    / "external_data"\n'
        '    / "AMPlify"\n'
        ')'
    )
    if re.search(amplify_pattern, source, flags=re.MULTILINE):
        source = replace_exact_block(source, amplify_pattern, amplify_replacement, "AMPlify path")

    return source


def prepare_notebook(args: argparse.Namespace) -> Path:
    script_path = Path(__file__).resolve()
    project_root = Path(args.project_root).resolve() if args.project_root else script_path.parent
    notebook_path = project_root / "benchmark" / "BenchMark17model_AUTHORITATIVE.ipynb"
    if not notebook_path.is_file():
        raise FileNotFoundError(notebook_path)

    observed = sha256_file(notebook_path)
    if observed != AUTHORITATIVE_NOTEBOOK_SHA256:
        raise RuntimeError(
            "Authoritative benchmark notebook hash mismatch.\n"
            f"Expected: {AUTHORITATIVE_NOTEBOOK_SHA256}\nObserved: {observed}"
        )

    work_root = Path(args.work_root).resolve() if args.work_root else project_root / "benchmark_work"
    source_root = Path(args.source_root).resolve() if args.source_root else project_root / "external_methods"

    conda_candidate = args.conda_exe or os.environ.get("CONDA_EXE") or shutil.which("conda")
    if conda_candidate:
        conda_exe = Path(conda_candidate).resolve()
        envs_root = conda_exe.parent.parent / "envs"
    else:
        # Preparation remains possible. Full execution of legacy-method cells will not be.
        conda_exe = project_root / "MISSING_CONDA_EXECUTABLE"
        envs_root = project_root / "MISSING_CONDA_ENVS"

    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    original_cell_count = len(notebook.get("cells", []))
    selected_cells = []
    for index, cell in enumerate(notebook.get("cells", [])):
        if index in SUPERSEDED_CELL_INDICES:
            continue
        cell = json.loads(json.dumps(cell))
        source = "".join(cell.get("source", []))
        if source:
            patched = patch_source(source, project_root, work_root, source_root, conda_exe, envs_root)
            cell["source"] = patched.splitlines(keepends=True)
        if cell.get("cell_type") == "code":
            # A reviewer rerun must not depend on stale outputs from the author machine.
            cell["outputs"] = []
            cell["execution_count"] = None
        selected_cells.append(cell)

    notebook["cells"] = selected_cells
    metadata = notebook.setdefault("metadata", {})
    metadata["genpept_portable_patch"] = {
        "source_sha256": AUTHORITATIVE_NOTEBOOK_SHA256,
        "source_file": "benchmark/BenchMark17model_AUTHORITATIVE.ipynb",
        "removed_superseded_cell_indices": sorted(SUPERSEDED_CELL_INDICES),
        "scientific_logic_modified": False,
        "changes": "machine-specific paths only; Cell 32 v4/v8 removed because v9 supersedes them",
        "original_cell_count": original_cell_count,
        "prepared_cell_count": len(selected_cells),
    }

    serialized = json.dumps(notebook, ensure_ascii=False, indent=1)
    forbidden = [ORIGINAL_PROJECT_ROOT, ORIGINAL_WORK_ROOT, ORIGINAL_SOURCE_ROOT, ORIGINAL_MASTER_PYTHON, ORIGINAL_CONDA]
    remaining = [item for item in forbidden if item in serialized]
    if remaining:
        raise RuntimeError(f"Portable patch incomplete; original absolute paths remain: {remaining}")

    output_path = project_root / args.output_notebook
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized + "\n", encoding="utf-8")

    print("AUTHORITATIVE_NOTEBOOK_SHA256:", observed)
    print("PREPARED_NOTEBOOK:", output_path)
    print("CELLS:", original_cell_count, "->", len(selected_cells))
    print("PROJECT_ROOT:", project_root)
    print("WORK_ROOT:", work_root)
    print("SOURCE_ROOT:", source_root)
    print("CONDA_EXE:", conda_exe)
    return output_path


def execute_notebook(path: Path) -> None:
    jupyter = shutil.which("jupyter")
    if not jupyter:
        raise RuntimeError("jupyter executable not found. Install Jupyter/nbconvert in the main benchmark environment.")
    cmd = [
        jupyter,
        "nbconvert",
        "--to",
        "notebook",
        "--execute",
        "--ExecutePreprocessor.timeout=-1",
        "--output",
        path.name.replace(".ipynb", "_EXECUTED.ipynb"),
        str(path),
    ]
    print("Executing:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=path.parent.parent)


def main() -> None:
    args = parse_args()
    if not args.prepare_only and not args.execute:
        args.prepare_only = True
    prepared = prepare_notebook(args)
    if args.execute:
        execute_notebook(prepared)


if __name__ == "__main__":
    main()
