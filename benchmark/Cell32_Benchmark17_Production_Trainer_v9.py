from __future__ import annotations

from pathlib import Path
from datetime import datetime
from types import SimpleNamespace
import gc
import hashlib
import json
import math
import os
import pickle
import random
import resource
import signal
import shutil
import subprocess
import tempfile
import time
import traceback

import numpy as np
import pandas as pd

from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC


# ============================================================
# CELL 32 - PRODUCTION TRAINER ONLY + RESUME-AWARE FINAL TRAINING LAUNCH v9
# Benchmark: benchmark4_FINAL_v1
#
# Operational-only environment controls:
#   Cell 32 auto-policy:
#   - GENPEPT_CELL32_MAX_JOBS=0 means all remaining jobs
#   - positive integer caps only this invocation (operational recovery/checkpoint)
#   - first clean invocation still trains exactly 1 locked job when MAX_JOBS=0
#   - rerunning the IDENTICAL cell resumes remaining jobs
#
# v9 operational/audit correction (NO scientific protocol change):
#   - retains the v4 Keras mask dtype fix
#   - AMPScannerV2 helper no longer imports pandas; the locked Python 3.6/TF1.2.1/Keras2.0.6 env does not provide it, and stdlib csv preserves identical score semantics
#   - recovers stale RUNNING rows left by KeyboardInterrupt without deleting evidence
#   - records/prints per-run wall/training/validation timing and rolling ETA
#   - records separable feature/token-load and preprocessing timing where it can
#     be measured without changing the locked method implementation
#   - emits a checkpoint resource log from canonical resources.json artifacts
#   - catches KeyboardInterrupt explicitly so the ledger never remains RUNNING
#
# IMPORTANT: v9 uses adaptive resource-aware parallel orchestration to reduce total wall-clock time.
# Per-run timing is still recorded, and v9 marks the execution context so parallel-run
# timing is not misrepresented as uncontended cross-model runtime. Scientific
# model/data/seed/threshold contracts are unchanged.
# ============================================================

PROJECT_ROOT = Path("/mnt/d/SUABAI_GenPept-Curated-2025_ 12.8.2026")
OUTPUT_ROOT = PROJECT_ROOT / "benchmark17_outputs"
AUDIT_DIR = OUTPUT_ROOT / "audit"
MANIFEST_DIR = OUTPUT_ROOT / "manifests"
RESULTS_ROOT = OUTPUT_ROOT / "results"
WORK_ROOT = Path("/home/pc/genpept_benchmark17_work")
RUN_ROOT = WORK_ROOT / "runs" / "benchmark4_FINAL_v1"
HELPER_ROOT = WORK_ROOT / "production_helpers" / "benchmark4_FINAL_v1"
HELPER_ROOT.mkdir(parents=True, exist_ok=True)

FULL_RUN_PLAN_PATH = AUDIT_DIR / "benchmark17_FINAL_full_680_run_plan.csv"
EXTERNAL_EVAL_PLAN_PATH = AUDIT_DIR / "benchmark17_FINAL_external_2040_eval_plan.csv"
RUN_LEDGER_PATH = AUDIT_DIR / "benchmark17_FINAL_run_ledger.csv"
LAUNCHER_CONTRACT_PATH = MANIFEST_DIR / "benchmark17_FINAL_training_launcher_contract.json"
FULL_EXECUTION_MANIFEST_PATH = MANIFEST_DIR / "benchmark17_FINAL_full_execution_manifest.json"
TRAINING_GATE_PATH = MANIFEST_DIR / "benchmark17_FINAL_training_gate.json"
BENCHMARK_MANIFEST_PATH = MANIFEST_DIR / "benchmark4_FINAL_v1_manifest.json"
CONTROLLED_CONTRACT_PATH = MANIFEST_DIR / "benchmark17_FINAL_controlled_model_contract.json"
EVALUATION_CONTRACT_PATH = MANIFEST_DIR / "benchmark17_FINAL_evaluation_contract.json"
FEATURE_CACHE_MANIFEST_PATH = MANIFEST_DIR / "benchmark17_FINAL_feature_cache_manifest.json"
PUBLISHED_ROUTE_CONTRACT_PATH = MANIFEST_DIR / "benchmark17_FINAL_published_method_route_contract.json"
RUN_ARTIFACT_SCHEMA_PATH = MANIFEST_DIR / "benchmark17_FINAL_run_artifact_schema.json"
CUDA_PATHS_FILE = MANIFEST_DIR / "cuda_runtime_paths.json"
ENVIRONMENT_AUDIT_PATH = MANIFEST_DIR / "environment_audit.json"
CELL32_RESOURCE_LOG_PATH = AUDIT_DIR / "benchmark17_CELL32_training_resource_log.csv"
CELL32_RESOURCE_SUMMARY_PATH = AUDIT_DIR / "benchmark17_CELL32_training_resource_summary.json"
CELL32_V9_EVENTS_PATH = AUDIT_DIR / "benchmark17_CELL32_v9_executor_events.csv"
CELL32_V9_SUMMARY_PATH = AUDIT_DIR / "benchmark17_CELL32_v9_executor_summary.json"
CELL32_TRAINER_SHA256 = os.environ.get("GENPEPT_CELL32_TRAINER_SHA256", "").strip()
SELF_TRAINER_PATH = PROJECT_ROOT / "scripts_manual" / "Cell32_Benchmark17_Production_Trainer_v9.py"

EXPECTED_DATASETS = [
    "GenPept_Curated_2025",
    "External_Benchmark_1",
    "AMPlify_balanced",
    "External_Benchmark_3",
]
EXPECTED_SEEDS = list(range(100, 110))
CONTROLLED_IDS = [
    "logistic_regression",
    "linear_svm_calibrated",
    "random_forest",
    "extra_trees",
    "hist_gradient_boosting",
    "lightgbm",
    "cnn1d",
    "bilstm",
    "bilstm_attention",
    "transformer_encoder",
]
PUBLISHED_IDS = [
    "ampep",
    "ampscannerv2",
    "ampgram",
    "ampir",
    "ampeppy",
    "ai4amp",
    "amplify",
]
ALL_IDS = CONTROLLED_IDS + PUBLISHED_IDS

PREDICTION_COLUMNS = [
    "run_id", "protocol", "model_id", "source_dataset", "target_dataset",
    "seed", "sample_id", "y_true", "y_score", "y_pred", "threshold",
]


def now_iso() -> str:
    return datetime.now().astimezone().isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()

# Exact-worker-routing guard: every child must execute this committed v9 file,
# never an inherited v6/v7 path from notebook globals/environment.
def _verify_self_trainer_path() -> str:
    assert SELF_TRAINER_PATH.is_file(), f"Missing committed v9 trainer: {SELF_TRAINER_PATH}"
    actual = sha256_file(SELF_TRAINER_PATH)
    if CELL32_TRAINER_SHA256:
        assert actual == CELL32_TRAINER_SHA256, (
            f"v9 self trainer hash mismatch: expected={CELL32_TRAINER_SHA256} actual={actual}"
        )
    return actual

if os.environ.get("GENPEPT_CELL32_SELF_CHECK_ONLY", "0") == "1":
    _self_sha = _verify_self_trainer_path()
    print(f"CELL32_V9_SELF_CHECK_PASS\t{SELF_TRAINER_PATH}\t{_self_sha}", flush=True)
    raise SystemExit(0)


def load_json(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def read_csv_exact(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False, keep_default_na=False, na_filter=False)


def atomic_write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=str(path.parent))
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False, sort_keys=True, allow_nan=False)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def atomic_write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as f:
            df.to_csv(f, index=False)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".tmp.", dir=str(path.parent))
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        tmp.write_text(text, encoding="utf-8", newline="\n")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def lock_helper(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        existing = path.read_text(encoding="utf-8")
        assert existing == text, f"Existing production helper differs; refusing overwrite: {path}"
    else:
        write_text_atomic(path, text)
    return sha256_file(path)


def run_cmd(cmd, *, env=None, timeout=None, cwd=None, label="command"):
    cmd = [str(x) for x in cmd]
    result = subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        timeout=timeout,
        cwd=None if cwd is None else str(cwd),
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"{label} failed with code {result.returncode}\n"
            f"COMMAND: {cmd}\nSTDOUT:\n{result.stdout[-12000:]}\n"
            f"STDERR:\n{result.stderr[-12000:]}"
        )
    return result


def bool_env(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "y"}


def artifact_entry(path: Path):
    path = Path(path)
    assert path.is_file(), f"Missing artifact: {path}"
    return {
        "path": str(path),
        "size_bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
    }


def process_peak_ram_mb() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    # Linux ru_maxrss is KiB.
    return value / 1024.0


def format_duration(seconds) -> str:
    if seconds is None or not np.isfinite(float(seconds)):
        return "NA"
    sec = max(0, int(round(float(seconds))))
    d, rem = divmod(sec, 86400)
    h, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    if d:
        return f"{d}d {h:02d}:{m:02d}:{s:02d}"
    return f"{h:02d}:{m:02d}:{s:02d}"


def existing_run_resource_rows(plan_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in plan_df.itertuples(index=False):
        run_dir = Path(str(row.run_dir))
        manifest_path = run_dir / "run_manifest.json"
        resources_path = run_dir / "resources.json"
        if not manifest_path.is_file() or not resources_path.is_file():
            continue
        try:
            manifest = load_json(manifest_path)
            resources = load_json(resources_path)
        except Exception:
            continue
        if str(manifest.get("status", "")) not in {"TRAINED", "EVALUATED", "COMPLETE"}:
            continue
        validator = globals().get("trained_manifest_valid")
        if callable(validator) and not validator(row):
            continue
        rec = {
            "run_id": str(row.run_id),
            "model_id": str(row.model_id),
            "source_dataset": str(row.source_dataset),
            "seed": int(row.seed),
            "trainer_version": str(manifest.get("cell32_trainer_version", "unknown")),
            "training_seconds": resources.get("training_seconds"),
            "validation_seconds": resources.get("validation_seconds"),
            "feature_or_token_load_seconds": resources.get("feature_or_token_load_seconds"),
            "preprocessing_seconds": resources.get("preprocessing_seconds"),
            "cell32_run_wall_seconds": resources.get("cell32_run_wall_seconds"),
            "peak_process_ram_mb": resources.get("peak_process_ram_mb"),
            "peak_gpu_vram_mb": resources.get("peak_gpu_vram_mb"),
            "model_parameter_count": resources.get("model_parameter_count"),
            "checkpoint_size_mb": resources.get("checkpoint_size_mb"),
            "device": resources.get("device"),
            "resources_path": str(resources_path),
            "resources_sha256": sha256_file(resources_path),
        }
        rows.append(rec)
    cols = [
        "run_id", "model_id", "source_dataset", "seed", "trainer_version",
        "training_seconds", "validation_seconds", "feature_or_token_load_seconds",
        "preprocessing_seconds", "cell32_run_wall_seconds", "peak_process_ram_mb",
        "peak_gpu_vram_mb", "model_parameter_count", "checkpoint_size_mb",
        "device", "resources_path", "resources_sha256",
    ]
    return pd.DataFrame(rows, columns=cols)


def write_resource_checkpoint(plan_df: pd.DataFrame, *, complete: bool) -> None:
    rdf = existing_run_resource_rows(plan_df)
    atomic_write_csv(CELL32_RESOURCE_LOG_PATH, rdf)
    numeric = {}
    for col in [
        "training_seconds", "validation_seconds", "feature_or_token_load_seconds",
        "preprocessing_seconds", "cell32_run_wall_seconds",
        "peak_process_ram_mb", "peak_gpu_vram_mb", "checkpoint_size_mb",
    ]:
        if col not in rdf.columns:
            continue
        vals = pd.to_numeric(rdf[col], errors="coerce").dropna()
        numeric[col] = {
            "n_measured": int(len(vals)),
            "sum": None if len(vals) == 0 else float(vals.sum()),
            "mean": None if len(vals) == 0 else float(vals.mean()),
            "median": None if len(vals) == 0 else float(vals.median()),
            "min": None if len(vals) == 0 else float(vals.min()),
            "max": None if len(vals) == 0 else float(vals.max()),
        }
    payload = {
        "schema_version": "1.0",
        "created_at": now_iso(),
        "benchmark_version": "benchmark4_FINAL_v1",
        "cell32_complete": bool(complete),
        "trained_resource_rows": int(len(rdf)),
        "timing_note": (
            "Existing pre-v6 runs retain their original resources.json. v6 does not fabricate "
            "missing load/preprocessing timings. Source-test inference timing remains forbidden "
            "in Cell 32 and will be measured in Cell 33."
        ),
        "resource_log": {
            "path": str(CELL32_RESOURCE_LOG_PATH),
            "sha256": sha256_file(CELL32_RESOURCE_LOG_PATH),
        },
        "numeric_summary": numeric,
    }
    atomic_write_json(CELL32_RESOURCE_SUMMARY_PATH, payload)


def gpu_peak_mb_best_effort():
    # Do NOT import TensorFlow in the notebook/master process.
    # Cell 2 established TensorFlow GPU access only in a subprocess with
    # the locked CUDA runtime search path. Deep controlled models use that
    # same subprocess strategy below.
    return None


def locked_tensorflow_subprocess_env(seed=None):
    assert CUDA_PATHS_FILE.is_file(), f"Missing locked CUDA path manifest: {CUDA_PATHS_FILE}"
    assert ENVIRONMENT_AUDIT_PATH.is_file(), f"Missing environment audit: {ENVIRONMENT_AUDIT_PATH}"
    cuda_payload = load_json(CUDA_PATHS_FILE)
    env_audit = load_json(ENVIRONMENT_AUDIT_PATH)
    dirs = [str(x) for x in cuda_payload.get("cuda_library_dirs", [])]
    assert len(dirs) > 1
    for d in dirs:
        assert Path(d).is_dir(), f"Locked CUDA library directory missing: {d}"
    assert int(env_audit["tensorflow"]["gpu_count"]) >= 1
    assert env_audit["tensorflow"]["is_cuda_build"] is True
    env = os.environ.copy()
    old = env.get("LD_LIBRARY_PATH", "")
    parts = list(dirs)
    if old:
        parts.append(old)
    env["LD_LIBRARY_PATH"] = ":".join(dict.fromkeys(parts))
    env["CUDA_VISIBLE_DEVICES"] = "0"
    env["TF_CPP_MIN_LOG_LEVEL"] = "2"
    env["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
    if seed is not None:
        env["PYTHONHASHSEED"] = str(int(seed))
    return env


def verify_tensorflow_gpu_subprocess():
    code = r"""
import json
import tensorflow as tf
gpus = tf.config.list_physical_devices('GPU')
assert len(gpus) >= 1, 'TensorFlow GPU unavailable under locked Cell-2 CUDA runtime path'
for gpu in gpus:
    try:
        tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError:
        pass
with tf.device('/GPU:0'):
    x = tf.constant([1.0, 2.0, 3.0], dtype=tf.float32)
    y = tf.reduce_sum(x * 2.0)
assert float(y.numpy()) == 12.0
print(json.dumps({'tensorflow': tf.__version__, 'gpu_count': len(gpus), 'gpu': str(gpus[0])}))
"""
    r = subprocess.run(
        [str(Path(os.sys.executable).resolve()), "-c", code],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=locked_tensorflow_subprocess_env(),
        check=False,
    )
    if r.returncode != 0:
        print("TensorFlow GPU preflight STDOUT:\n" + r.stdout)
        print("TensorFlow GPU preflight STDERR:\n" + r.stderr)
    assert r.returncode == 0, "Locked TensorFlow GPU preflight failed; refusing production training."
    lines = [x for x in r.stdout.splitlines() if x.strip()]
    assert lines
    info = json.loads(lines[-1])
    assert int(info["gpu_count"]) >= 1
    return info


DEEP_WORKER_CODE = r"""
from pathlib import Path
import json, sys, time
import numpy as np
import pandas as pd
import tensorflow as tf

MODEL_ID=sys.argv[1]
SEED=int(sys.argv[2])
TOKEN_PATH=Path(sys.argv[3])
LABEL_PATH=Path(sys.argv[4])
INDEX_PATH=Path(sys.argv[5])
TRAIN_KEY=sys.argv[6]
VAL_KEY=sys.argv[7]
CONTROLLED_CONTRACT=Path(sys.argv[8])
OUT_DIR=Path(sys.argv[9])
OUT_DIR.mkdir(parents=True, exist_ok=True)

gpus=tf.config.list_physical_devices('GPU')
assert gpus, 'GPU missing in deep production worker'
for gpu in gpus:
    try:
        tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError:
        pass
tf.keras.utils.set_random_seed(SEED)
try:
    tf.config.experimental.enable_op_determinism()
except Exception:
    pass

payload=json.load(open(CONTROLLED_CONTRACT, encoding='utf-8'))
contract=payload['contract']
configs=contract['model_configs']
deep=contract['deep_shared_training']
cfg=configs[MODEL_ID]['architecture']
load_t0=time.perf_counter()
X=np.load(TOKEN_PATH, mmap_mode='r', allow_pickle=False)
y=np.load(LABEL_PATH, mmap_mode='r', allow_pickle=False)
idx=np.load(INDEX_PATH, allow_pickle=False)
tri=np.asarray(idx[TRAIN_KEY], dtype=np.int64)
vai=np.asarray(idx[VAL_KEY], dtype=np.int64)
load_sec=time.perf_counter()-load_t0
prep_t0=time.perf_counter()
Xtr=np.asarray(X[tri], dtype=np.int32)
ytr=np.asarray(y[tri], dtype=np.float32)
Xv=np.asarray(X[vai], dtype=np.int32)
yv=np.asarray(y[vai], dtype=np.float32)
prep_sec=time.perf_counter()-prep_t0
K=tf.keras

class MaskedGlobalMaxPooling1D(K.layers.Layer):
    def call(self, inputs):
        x, mask=inputs
        mask=tf.cast(mask,x.dtype)
        masked=tf.where(tf.expand_dims(mask>0,-1),x,tf.cast(-1e9,x.dtype))
        return tf.reduce_max(masked,axis=1)

class MaskedGlobalAveragePooling1D(K.layers.Layer):
    def call(self, inputs):
        x,mask=inputs
        mask=tf.cast(mask,x.dtype)
        num=tf.reduce_sum(x*tf.expand_dims(mask,-1),axis=1)
        den=tf.maximum(tf.reduce_sum(mask,axis=1,keepdims=True),tf.cast(1.0,x.dtype))
        return num/den

class MaskedAdditiveAttention(K.layers.Layer):
    def __init__(self,units,**kwargs):
        super().__init__(**kwargs)
        self.proj=K.layers.Dense(units,activation='tanh')
        self.score=K.layers.Dense(1,use_bias=False)
    def call(self,inputs):
        x,mask=inputs
        mask=tf.cast(mask,tf.bool)
        s=tf.squeeze(self.score(self.proj(x)),axis=-1)
        s=tf.where(mask,s,tf.cast(-1e9,s.dtype))
        a=tf.nn.softmax(s,axis=1)
        return tf.reduce_sum(x*tf.expand_dims(a,-1),axis=1)

class AddLearnedPositionEmbedding(K.layers.Layer):
    def __init__(self,max_positions,model_dim,**kwargs):
        super().__init__(**kwargs)
        self.pos_embedding=K.layers.Embedding(max_positions,model_dim)
    def call(self,x):
        positions=tf.range(start=0,limit=tf.shape(x)[1],delta=1)
        return x+self.pos_embedding(positions)

tokens=K.Input(shape=(200,),dtype='int32',name='tokens')
mask=K.layers.Lambda(lambda t: tf.not_equal(t,0),name='padding_mask')(tokens)

if MODEL_ID=='cnn1d':
    z=K.layers.Embedding(21,int(cfg['embedding_dim']),mask_zero=True)(tokens)
    z=K.layers.Conv1D(int(cfg['conv1_filters']),int(cfg['conv1_kernel_size']),padding='same',activation='relu')(z)
    z=K.layers.Conv1D(int(cfg['conv2_filters']),int(cfg['conv2_kernel_size']),padding='same',activation='relu')(z)
    z=MaskedGlobalMaxPooling1D()([z,mask])
    z=K.layers.Dense(int(cfg['dense_units']),activation='relu')(z)
    z=K.layers.Dropout(float(cfg['dropout']))(z)
    out=K.layers.Dense(1,activation='sigmoid')(z)
elif MODEL_ID=='bilstm':
    z=K.layers.Embedding(21,int(cfg['embedding_dim']),mask_zero=True)(tokens)
    z=K.layers.Bidirectional(K.layers.LSTM(int(cfg['bilstm_units_per_direction']),return_sequences=False,dropout=float(cfg['lstm_dropout']),recurrent_dropout=float(cfg['recurrent_dropout'])))(z)
    z=K.layers.Dense(int(cfg['dense_units']),activation='relu')(z)
    z=K.layers.Dropout(float(cfg['dense_dropout']))(z)
    out=K.layers.Dense(1,activation='sigmoid')(z)
elif MODEL_ID=='bilstm_attention':
    z=K.layers.Embedding(21,int(cfg['embedding_dim']),mask_zero=True)(tokens)
    z=K.layers.Bidirectional(K.layers.LSTM(int(cfg['bilstm_units_per_direction']),return_sequences=True,dropout=float(cfg['lstm_dropout']),recurrent_dropout=float(cfg['recurrent_dropout'])))(z)
    z=MaskedAdditiveAttention(int(cfg['attention_units']))([z,mask])
    z=K.layers.Dense(int(cfg['dense_units']),activation='relu')(z)
    z=K.layers.Dropout(float(cfg['dense_dropout']))(z)
    out=K.layers.Dense(1,activation='sigmoid')(z)
elif MODEL_ID=='transformer_encoder':
    d=int(cfg['model_dimension'])
    tok=K.layers.Embedding(21,d,mask_zero=True)(tokens)
    z=AddLearnedPositionEmbedding(200,d,name='position_embedding')(tok)
    am=K.layers.Lambda(lambda m: tf.expand_dims(m,axis=1),name='attention_mask')(mask)
    for block in range(int(cfg['n_encoder_blocks'])):
        a=K.layers.MultiHeadAttention(num_heads=int(cfg['num_attention_heads']),key_dim=int(cfg['key_dim_per_head']),dropout=float(cfg['attention_dropout']),name='mha_%d'%(block+1))(z,z,attention_mask=am)
        z=K.layers.LayerNormalization(name='attn_ln_%d'%(block+1))(z+a)
        ff=K.layers.Dense(int(cfg['feed_forward_dimension']),activation='relu')(z)
        ff=K.layers.Dropout(float(cfg['feed_forward_dropout']))(ff)
        ff=K.layers.Dense(d)(ff)
        z=K.layers.LayerNormalization(name='ff_ln_%d'%(block+1))(z+ff)
    z=MaskedGlobalAveragePooling1D()([z,mask])
    z=K.layers.Dense(int(cfg['dense_units']),activation='relu')(z)
    z=K.layers.Dropout(float(cfg['dense_dropout']))(z)
    out=K.layers.Dense(1,activation='sigmoid')(z)
else:
    raise KeyError(MODEL_ID)

model=K.Model(tokens,out,name=MODEL_ID)
model.compile(
    optimizer=K.optimizers.Adam(
        learning_rate=float(deep['learning_rate']),
        beta_1=float(deep['beta_1']),
        beta_2=float(deep['beta_2']),
        epsilon=float(deep['epsilon']),
        clipnorm=float(deep['clipnorm']),
    ),
    loss='binary_crossentropy',
)
callbacks=[
    K.callbacks.EarlyStopping(
        monitor='val_loss', mode='min',
        patience=int(deep['early_stopping']['patience']),
        restore_best_weights=True,
        min_delta=float(deep['early_stopping']['min_delta']),
    ),
    K.callbacks.ReduceLROnPlateau(
        monitor='val_loss', mode='min',
        factor=float(deep['reduce_lr_on_plateau']['factor']),
        patience=int(deep['reduce_lr_on_plateau']['patience']),
        min_lr=float(deep['reduce_lr_on_plateau']['minimum_learning_rate']),
    ),
]
t0=time.perf_counter()
hist=model.fit(
    Xtr,ytr,
    validation_data=(Xv,yv),
    epochs=int(deep['maximum_epochs']),
    batch_size=int(deep['batch_size']),
    shuffle=True,
    callbacks=callbacks,
    verbose=0,
)
train_sec=time.perf_counter()-t0
pd.DataFrame({
    'epoch':np.arange(1,len(hist.history['loss'])+1),
    'loss':hist.history['loss'],
    'val_loss':hist.history['val_loss'],
}).to_csv(OUT_DIR/'history.csv',index=False)
cp=OUT_DIR/'checkpoint_controlled.weights.h5'
model.save_weights(cp)
t1=time.perf_counter()
score=np.asarray(model.predict(Xv,batch_size=256,verbose=0),dtype=float).reshape(-1)
val_sec=time.perf_counter()-t1
assert np.isfinite(score).all() and ((score>=0)&(score<=1)).all()
np.save(OUT_DIR/'validation_scores.npy',score,allow_pickle=False)
peak=None
try:
    peak=float(tf.config.experimental.get_memory_info('GPU:0').get('peak',0))/(1024.0**2)
except Exception:
    pass
summary={
    'tensorflow':tf.__version__,
    'gpu_count':len(tf.config.list_physical_devices('GPU')),
    'seed':SEED,
    'feature_or_token_load_seconds':load_sec,
    'preprocessing_seconds':prep_sec,
    'training_seconds':train_sec,
    'validation_seconds':val_sec,
    'model_parameter_count':int(model.count_params()),
    'peak_gpu_vram_mb':peak,
    'checkpoint':'checkpoint_controlled.weights.h5',
}
json.dump(summary,open(OUT_DIR/'deep_worker_summary.json','w'),indent=2)
"""


def run_controlled_deep_subprocess(model_id, seed, run_row, staging_dir):
    worker_path = staging_dir / "controlled_deep_worker.py"
    write_text_atomic(worker_path, DEEP_WORKER_CODE)
    env = locked_tensorflow_subprocess_env(seed)
    cmd = [
        str(Path(os.sys.executable).resolve()),
        str(worker_path),
        str(model_id),
        str(int(seed)),
        str(TOKEN_CACHE_PATH),
        str(LABEL_CACHE_PATH),
        str(index_path),
        str(run_row.train_index_key),
        str(run_row.validation_index_key),
        str(CONTROLLED_CONTRACT_PATH),
        str(staging_dir),
    ]
    r = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        check=False,
    )
    if r.returncode != 0:
        print("Controlled deep worker STDOUT:\n" + r.stdout)
        print("Controlled deep worker STDERR:\n" + r.stderr)
    assert r.returncode == 0, f"Controlled deep GPU worker failed for {model_id}"
    summary = load_json(staging_dir / "deep_worker_summary.json")
    assert int(summary["gpu_count"]) >= 1
    score = np.load(staging_dir / "validation_scores.npy", allow_pickle=False)
    assert score.ndim == 1
    assert np.isfinite(score).all()
    assert ((score >= 0.0) & (score <= 1.0)).all()
    checkpoint = staging_dir / str(summary["checkpoint"])
    assert checkpoint.is_file()
    (staging_dir / "validation_scores.npy").unlink()
    (staging_dir / "deep_worker_summary.json").unlink()
    (staging_dir / "controlled_deep_worker.py").unlink()
    return (
        score,
        checkpoint,
        float(summary["training_seconds"]),
        float(summary["validation_seconds"]),
        int(summary["model_parameter_count"]),
        summary.get("peak_gpu_vram_mb"),
        float(summary.get("feature_or_token_load_seconds", 0.0)),
        float(summary.get("preprocessing_seconds", 0.0)),
    )


def binary_ece(y_true, y_score, n_bins=10):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    assert len(y_true) == len(y_score) and len(y_true) > 0
    assert np.isfinite(y_score).all()
    assert (((y_score >= 0.0) & (y_score <= 1.0))).all()
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    n_total = len(y_true)
    for i in range(n_bins):
        left, right = edges[i], edges[i + 1]
        if i < n_bins - 1:
            mask = (y_score >= left) & (y_score < right)
        else:
            mask = (y_score >= left) & (y_score <= right)
        n_bin = int(mask.sum())
        if n_bin == 0:
            continue
        observed_fraction = float(y_true[mask].mean())
        mean_probability = float(y_score[mask].mean())
        ece += (n_bin / n_total) * abs(observed_fraction - mean_probability)
    return float(ece)


def select_validation_threshold(y_true, y_score):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    assert len(y_true) == len(y_score) and len(y_true) > 0
    assert set(np.unique(y_true).tolist()) == {0, 1}
    assert np.isfinite(y_score).all()
    assert (((y_score >= 0.0) & (y_score <= 1.0))).all()
    unique_scores = np.unique(y_score)
    if len(unique_scores) >= 2:
        midpoints = (unique_scores[:-1] + unique_scores[1:]) / 2.0
    else:
        midpoints = np.array([], dtype=float)
    candidates = np.unique(np.concatenate([
        np.array([0.0, 0.5, 1.0], dtype=float), unique_scores, midpoints
    ]))
    candidates = candidates[(candidates >= 0.0) & (candidates <= 1.0)]
    rows = []
    for tau in candidates:
        y_pred = (y_score >= tau).astype(int)
        rows.append({
            "tau": float(tau),
            "MCC": float(matthews_corrcoef(y_true, y_pred)),
            "BACC": float(balanced_accuracy_score(y_true, y_pred)),
            "F1": float(f1_score(y_true, y_pred, zero_division=0)),
            "distance_to_0_5": abs(float(tau) - 0.5),
        })
    table = pd.DataFrame(rows).sort_values(
        ["MCC", "BACC", "F1", "distance_to_0_5", "tau"],
        ascending=[False, False, False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    best = table.iloc[0]
    return {
        "tau": float(best["tau"]),
        "validation_MCC": float(best["MCC"]),
        "validation_BACC": float(best["BACC"]),
        "validation_F1": float(best["F1"]),
        "n_threshold_candidates": int(len(table)),
    }


def compute_binary_metrics(y_true, y_score, tau):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    assert len(y_true) == len(y_score) and len(y_true) > 0
    assert set(np.unique(y_true).tolist()) == {0, 1}
    assert np.isfinite(y_score).all()
    assert (((y_score >= 0.0) & (y_score <= 1.0))).all()
    assert 0.0 <= float(tau) <= 1.0
    y_pred = (y_score >= float(tau)).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    return {
        "AUROC": float(roc_auc_score(y_true, y_score)),
        "AUPRC": float(average_precision_score(y_true, y_score)),
        "MCC": float(matthews_corrcoef(y_true, y_pred)),
        "BACC": float(balanced_accuracy_score(y_true, y_pred)),
        "F1": float(f1_score(y_true, y_pred, zero_division=0)),
        "ACC": float(accuracy_score(y_true, y_pred)),
        "Sensitivity": float(sensitivity),
        "Specificity": float(specificity),
        "Brier": float(brier_score_loss(y_true, y_score)),
        "ECE10": float(binary_ece(y_true, y_score, n_bins=10)),
        "threshold": float(tau),
        "TP": int(tp), "TN": int(tn), "FP": int(fp), "FN": int(fn),
        "n": int(len(y_true)),
        "n_positive": int((y_true == 1).sum()),
        "n_negative": int((y_true == 0).sum()),
    }


def prediction_frame(*, run_id, protocol, model_id, source_dataset, target_dataset,
                     seed, sample_ids, y_true, y_score, tau):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    y_pred = (y_score >= float(tau)).astype(int)
    out = pd.DataFrame({
        "run_id": run_id,
        "protocol": protocol,
        "model_id": model_id,
        "source_dataset": source_dataset,
        "target_dataset": target_dataset,
        "seed": int(seed),
        "sample_id": list(map(str, sample_ids)),
        "y_true": y_true,
        "y_score": y_score,
        "y_pred": y_pred,
        "threshold": float(tau),
    })
    assert list(out.columns) == PREDICTION_COLUMNS
    return out


def write_na_history(path: Path, note: str):
    atomic_write_csv(path, pd.DataFrame({
        "epoch": [0], "loss": [np.nan], "val_loss": [np.nan], "note": [note]
    }))


def write_fasta(df: pd.DataFrame, path: Path):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for row in df.itertuples(index=False):
            f.write(f">{row.sample_id}\n{row.sequence}\n")


def write_class_fastas(df: pd.DataFrame, amp_path: Path, nonamp_path: Path):
    write_fasta(df.loc[df["label"].astype(int).eq(1)].copy(), amp_path)
    write_fasta(df.loc[df["label"].astype(int).eq(0)].copy(), nonamp_path)


# ============================================================
# PRE-FLIGHT: LOCKED CONTRACTS AND PLANS
# ============================================================

for p in [
    PROJECT_ROOT, OUTPUT_ROOT, AUDIT_DIR, MANIFEST_DIR, RESULTS_ROOT,
    WORK_ROOT, RUN_ROOT,
]:
    assert p.is_dir(), f"Required directory missing: {p}"

for p in [
    FULL_RUN_PLAN_PATH, EXTERNAL_EVAL_PLAN_PATH, RUN_LEDGER_PATH,
    LAUNCHER_CONTRACT_PATH, FULL_EXECUTION_MANIFEST_PATH, TRAINING_GATE_PATH,
    BENCHMARK_MANIFEST_PATH, CONTROLLED_CONTRACT_PATH, EVALUATION_CONTRACT_PATH,
    FEATURE_CACHE_MANIFEST_PATH, PUBLISHED_ROUTE_CONTRACT_PATH,
    RUN_ARTIFACT_SCHEMA_PATH,
]:
    assert p.is_file(), f"Required locked artifact missing: {p}"

launcher_wrapper = load_json(LAUNCHER_CONTRACT_PATH)
full_exec_wrapper = load_json(FULL_EXECUTION_MANIFEST_PATH)
gate = load_json(TRAINING_GATE_PATH)
benchmark_manifest = load_json(BENCHMARK_MANIFEST_PATH)
controlled_wrapper = load_json(CONTROLLED_CONTRACT_PATH)
evaluation_contract = load_json(EVALUATION_CONTRACT_PATH)
feature_manifest = load_json(FEATURE_CACHE_MANIFEST_PATH)
route_wrapper = load_json(PUBLISHED_ROUTE_CONTRACT_PATH)
artifact_schema_wrapper = load_json(RUN_ARTIFACT_SCHEMA_PATH)

assert launcher_wrapper["status"] == "TRAINING_LAUNCHER_CONTRACT_LOCKED"
launcher = launcher_wrapper["contract"]
assert launcher["status"] == "READY_TO_LAUNCH_FINAL_17_MODEL_BENCHMARK"
assert launcher["training_gate_open"] is True
assert launcher["total_models"] == 17
assert launcher["protocol_1_training_jobs"] == 680
assert launcher["protocol_2_external_evaluations"] == 2040
assert launcher["datasets"] == EXPECTED_DATASETS
assert launcher["model_seeds"] == EXPECTED_SEEDS
assert launcher["execution_order"]["model_order"] == ALL_IDS
assert launcher["execution_order"]["source_dataset_order"] == EXPECTED_DATASETS
assert launcher["execution_order"]["seed_order"] == EXPECTED_SEEDS
assert gate["open"] is True and gate["benchmark_training_allowed"] is True
assert gate["published_methods_implementation_validated"] == 7
assert sorted(gate["validated_methods"]) == sorted(PUBLISHED_IDS)
assert route_wrapper["status"] == "PUBLISHED_METHOD_ROUTES_LOCKED_BEFORE_TRAINING"
assert route_wrapper["contract"]["published_model_ids"] == PUBLISHED_IDS

# Re-hash every immutable input recorded by Cell 31.
for key, entry in launcher["immutable_inputs"].items():
    path = Path(entry["path"])
    assert path.exists(), f"Immutable input missing [{key}]: {path}"
    if path.is_file():
        assert sha256_file(path) == entry["sha256"], f"Immutable input changed [{key}]"

full_plan = read_csv_exact(FULL_RUN_PLAN_PATH)
external_plan = read_csv_exact(EXTERNAL_EVAL_PLAN_PATH)
ledger = read_csv_exact(RUN_LEDGER_PATH)
assert len(full_plan) == 680 and full_plan["run_id"].is_unique
assert len(external_plan) == 2040
assert len(ledger) == 680 and ledger["run_id"].is_unique
assert full_plan["model_id"].drop_duplicates().tolist() == ALL_IDS
assert pd.to_numeric(full_plan["seed"], errors="raise").drop_duplicates().astype(int).tolist() == EXPECTED_SEEDS
assert external_plan.groupby("run_id").size().eq(3).all()
assert pd.to_numeric(external_plan["residual_exact_overlap"], errors="raise").eq(0).all()

controlled_contract = controlled_wrapper["contract"]
assert controlled_contract["controlled_model_ids"] == CONTROLLED_IDS
assert controlled_contract["model_seeds"] == EXPECTED_SEEDS
MODEL_CONFIGS = controlled_contract["model_configs"]
DEEP_SHARED = controlled_contract["deep_shared_training"]

# Caches and global row mapping.
cache_files = feature_manifest["files"]
CLASSICAL_CACHE_PATH = Path(cache_files["X_classical"]["path"])
TOKEN_CACHE_PATH = Path(cache_files["X_tokens"]["path"])
LABEL_CACHE_PATH = Path(cache_files["y"]["path"])
ROW_METADATA_PATH = Path(cache_files["row_metadata"]["path"])
for key, path in [
    ("X_classical", CLASSICAL_CACHE_PATH), ("X_tokens", TOKEN_CACHE_PATH),
    ("y", LABEL_CACHE_PATH), ("row_metadata", ROW_METADATA_PATH),
]:
    assert path.is_file()
    assert sha256_file(path) == cache_files[key]["sha256"]

X_CLASSICAL = np.load(CLASSICAL_CACHE_PATH, mmap_mode="r", allow_pickle=False)
X_TOKENS = np.load(TOKEN_CACHE_PATH, mmap_mode="r", allow_pickle=False)
Y = np.load(LABEL_CACHE_PATH, mmap_mode="r", allow_pickle=False)
assert X_CLASSICAL.shape == (41059, 421)
assert X_TOKENS.shape == (41059, 200)
assert Y.shape == (41059,)

# Verify the canonical cohort artifact WITHOUT loading benchmark rows.
final_cohort_path = Path(benchmark_manifest["final_cohort"]["path"])
assert final_cohort_path.is_file()
assert sha256_file(final_cohort_path) == benchmark_manifest["final_cohort"]["sha256"]

index_path = Path(launcher["immutable_inputs"]["execution_index"]["path"])
assert index_path.is_file()
INDEX = np.load(index_path, allow_pickle=False)
assert len(INDEX.files) == 24

# Cell 32 may materialize SOURCE TRAIN + SOURCE VALIDATION metadata only.
# No source-test or cross-target index key is included here.
TRAINVAL_KEYS = []
for _dataset in EXPECTED_DATASETS:
    TRAINVAL_KEYS.extend([
        f"within__{_dataset}__train",
        f"within__{_dataset}__validation",
    ])
for _key in TRAINVAL_KEYS:
    assert _key in INDEX.files
TRAINVAL_INDICES = np.unique(np.concatenate([
    np.asarray(INDEX[_key], dtype=np.int64)
    for _key in TRAINVAL_KEYS
]))
TRAINVAL_INDEX_SET = set(TRAINVAL_INDICES.tolist())

ROW_META_TRAINVAL = pd.read_csv(
    ROW_METADATA_PATH,
    low_memory=False,
    skiprows=lambda line_no: (line_no > 0 and (line_no - 1) not in TRAINVAL_INDEX_SET),
)
required_meta_columns = {
    "row_index", "sample_id", "dataset", "split", "label",
    "sequence_length", "sequence_sha256",
}
assert required_meta_columns.issubset(ROW_META_TRAINVAL.columns)
assert len(ROW_META_TRAINVAL) == len(TRAINVAL_INDICES)
assert ROW_META_TRAINVAL["row_index"].astype(np.int64).is_unique
assert set(ROW_META_TRAINVAL["split"].astype(str).unique()) <= {"train", "validation"}
assert not ROW_META_TRAINVAL["split"].astype(str).eq("test").any()
ROW_META_BY_INDEX = ROW_META_TRAINVAL.set_index(
    ROW_META_TRAINVAL["row_index"].astype(np.int64),
    drop=False,
)

TOKEN_TO_ID = controlled_contract["deep_representation"]["token_to_id"]
assert TOKEN_TO_ID["PAD"] == 0
ID_TO_AA = {int(v): str(k) for k, v in TOKEN_TO_ID.items() if str(k) != "PAD"}
assert set(ID_TO_AA) == set(range(1, 21))
assert set(ID_TO_AA.values()) == set("ACDEFGHIKLMNPQRSTVWY")


def decode_cached_sequence(row_index: int, expected_length: int) -> str:
    tokens = np.asarray(X_TOKENS[int(row_index)], dtype=np.int64)
    nonpad = tokens[tokens != 0]
    assert len(nonpad) == int(expected_length)
    seq = "".join(ID_TO_AA[int(token)] for token in nonpad.tolist())
    assert len(seq) == int(expected_length)
    return seq


def df_from_index_key(key: str) -> pd.DataFrame:
    # Hard Cell-32 boundary: test/cross keys are forbidden.
    assert key in TRAINVAL_KEYS, f"Cell 32 forbidden index key: {key}"
    idx = np.asarray(INDEX[key], dtype=np.int64)
    assert set(idx.tolist()).issubset(TRAINVAL_INDEX_SET)
    meta = ROW_META_BY_INDEX.loc[idx].copy()
    assert meta["row_index"].astype(np.int64).tolist() == idx.tolist()

    sequences = []
    for row in meta.itertuples(index=False):
        seq = decode_cached_sequence(int(row.row_index), int(row.sequence_length))
        observed_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        assert observed_sha == str(row.sequence_sha256)
        sequences.append(seq)

    df = pd.DataFrame({
        "sample_id": meta["sample_id"].astype(str).tolist(),
        "sequence": sequences,
        "label": meta["label"].astype(np.int8).tolist(),
        "dataset": meta["dataset"].astype(str).tolist(),
        "split": meta["split"].astype(str).tolist(),
        "row_index": idx,
        "length": meta["sequence_length"].astype(int).tolist(),
    })
    assert np.array_equal(
        df["label"].to_numpy(np.int8),
        np.asarray(Y[idx], dtype=np.int8),
    )
    assert df["sequence"].map(
        lambda s: set(str(s)).issubset(set("ACDEFGHIKLMNPQRSTVWY"))
    ).all()
    assert df["length"].between(10, 200).all()
    return df.reset_index(drop=True)


# ============================================================
# CONTROLLED MODELS
# ============================================================

def set_master_seed(seed: int):
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    try:
        import tensorflow as tf
        tf.keras.utils.set_random_seed(seed)
        try:
            tf.config.experimental.enable_op_determinism()
        except Exception:
            pass
    except Exception:
        pass


def make_classical_model(model_id: str, seed: int):
    cfg = MODEL_CONFIGS[model_id]
    if model_id == "logistic_regression":
        e = cfg["estimator"]
        return Pipeline([
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            ("model", LogisticRegression(
                penalty=e["penalty"], C=float(e["C"]), solver=e["solver"],
                max_iter=int(e["max_iter"]), tol=float(e["tol"]),
                class_weight=None, random_state=None,
            )),
        ])
    if model_id == "linear_svm_calibrated":
        e = cfg["base_estimator"]
        base = LinearSVC(
            C=float(e["C"]), loss=e["loss"], penalty=e["penalty"],
            dual=e["dual"], tol=float(e["tol"]), max_iter=int(e["max_iter"]),
            class_weight=None, random_state=seed,
        )
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
        try:
            cal = CalibratedClassifierCV(estimator=base, method="sigmoid", cv=cv, ensemble=True)
        except TypeError:
            cal = CalibratedClassifierCV(base_estimator=base, method="sigmoid", cv=cv, ensemble=True)
        return Pipeline([("scaler", StandardScaler()), ("model", cal)])
    if model_id == "random_forest":
        e = cfg["estimator"]
        return RandomForestClassifier(
            n_estimators=int(e["n_estimators"]), criterion=e["criterion"],
            max_depth=None, min_samples_split=int(e["min_samples_split"]),
            min_samples_leaf=int(e["min_samples_leaf"]), max_features=e["max_features"],
            bootstrap=bool(e["bootstrap"]), class_weight=None, n_jobs=-1,
            random_state=seed,
        )
    if model_id == "extra_trees":
        e = cfg["estimator"]
        return ExtraTreesClassifier(
            n_estimators=int(e["n_estimators"]), criterion=e["criterion"],
            max_depth=None, min_samples_split=int(e["min_samples_split"]),
            min_samples_leaf=int(e["min_samples_leaf"]), max_features=e["max_features"],
            bootstrap=bool(e["bootstrap"]), class_weight=None, n_jobs=-1,
            random_state=seed,
        )
    if model_id == "hist_gradient_boosting":
        e = cfg["estimator"]
        return HistGradientBoostingClassifier(
            loss=e["loss"], learning_rate=float(e["learning_rate"]),
            max_iter=int(e["max_iter"]), max_leaf_nodes=int(e["max_leaf_nodes"]),
            max_depth=None, min_samples_leaf=int(e["min_samples_leaf"]),
            l2_regularization=float(e["l2_regularization"]), max_bins=int(e["max_bins"]),
            early_stopping=False, random_state=seed,
        )
    if model_id == "lightgbm":
        import lightgbm as lgb
        e = cfg["estimator"]
        return lgb.LGBMClassifier(
            objective="binary", n_estimators=int(e["n_estimators"]),
            learning_rate=float(e["learning_rate"]), num_leaves=int(e["num_leaves"]),
            max_depth=int(e["max_depth"]), min_child_samples=int(e["min_child_samples"]),
            subsample=float(e["subsample"]), subsample_freq=int(e["subsample_freq"]),
            colsample_bytree=float(e["colsample_bytree"]), reg_alpha=float(e["reg_alpha"]),
            reg_lambda=float(e["reg_lambda"]), class_weight=None, n_jobs=-1,
            random_state=seed, deterministic=True, force_col_wise=True, verbosity=-1,
            device_type="cpu",
        )
    raise KeyError(model_id)


def positive_probability(model, X):
    p = np.asarray(model.predict_proba(X), dtype=float)
    classes = np.asarray(model.classes_ if hasattr(model, "classes_") else model[-1].classes_)
    # Pipeline exposes classes_ in current sklearn; fallback handles older versions.
    if not hasattr(model, "classes_") and hasattr(model, "named_steps"):
        classes = np.asarray(model.named_steps["model"].classes_)
    pos = int(np.where(classes == 1)[0][0])
    out = p[:, pos]
    assert np.isfinite(out).all() and ((out >= 0.0) & (out <= 1.0)).all()
    return out


def build_deep_model(model_id: str):
    import tensorflow as tf
    K = tf.keras
    cfg = MODEL_CONFIGS[model_id]["architecture"]

    class MaskedGlobalMaxPooling1D(K.layers.Layer):
        def call(self, inputs):
            x, mask = inputs
            mask = tf.cast(mask, x.dtype)
            masked = tf.where(tf.expand_dims(mask > 0, -1), x, tf.cast(-1e9, x.dtype))
            return tf.reduce_max(masked, axis=1)

    class MaskedGlobalAveragePooling1D(K.layers.Layer):
        def call(self, inputs):
            x, mask = inputs
            mask = tf.cast(mask, x.dtype)
            num = tf.reduce_sum(x * tf.expand_dims(mask, -1), axis=1)
            den = tf.maximum(tf.reduce_sum(mask, axis=1, keepdims=True), tf.cast(1.0, x.dtype))
            return num / den

    class MaskedAdditiveAttention(K.layers.Layer):
        def __init__(self, units, **kwargs):
            super().__init__(**kwargs)
            self.proj = K.layers.Dense(units, activation="tanh")
            self.score = K.layers.Dense(1, use_bias=False)
        def call(self, inputs):
            x, mask = inputs
            mask = tf.cast(mask, tf.bool)
            s = tf.squeeze(self.score(self.proj(x)), axis=-1)
            s = tf.where(mask, s, tf.cast(-1e9, s.dtype))
            a = tf.nn.softmax(s, axis=1)
            return tf.reduce_sum(x * tf.expand_dims(a, -1), axis=1)

    class AddLearnedPositionEmbedding(K.layers.Layer):
        def __init__(self, max_positions, model_dim, **kwargs):
            super().__init__(**kwargs)
            self.pos_embedding = K.layers.Embedding(max_positions, model_dim)
        def call(self, x):
            positions = tf.range(start=0, limit=tf.shape(x)[1], delta=1)
            return x + self.pos_embedding(positions)

    tokens = K.Input(shape=(200,), dtype="int32", name="tokens")
    mask = K.layers.Lambda(lambda t: tf.not_equal(t, 0), name="padding_mask")(tokens)

    if model_id == "cnn1d":
        x = K.layers.Embedding(21, int(cfg["embedding_dim"]), mask_zero=True)(tokens)
        x = K.layers.Conv1D(int(cfg["conv1_filters"]), int(cfg["conv1_kernel_size"]),
                            padding="same", activation="relu")(x)
        x = K.layers.Conv1D(int(cfg["conv2_filters"]), int(cfg["conv2_kernel_size"]),
                            padding="same", activation="relu")(x)
        x = MaskedGlobalMaxPooling1D()([x, mask])
        x = K.layers.Dense(int(cfg["dense_units"]), activation="relu")(x)
        x = K.layers.Dropout(float(cfg["dropout"]))(x)
        out = K.layers.Dense(1, activation="sigmoid")(x)

    elif model_id == "bilstm":
        x = K.layers.Embedding(21, int(cfg["embedding_dim"]), mask_zero=True)(tokens)
        x = K.layers.Bidirectional(K.layers.LSTM(
            int(cfg["bilstm_units_per_direction"]), return_sequences=False,
            dropout=float(cfg["lstm_dropout"]), recurrent_dropout=float(cfg["recurrent_dropout"]),
        ))(x)
        x = K.layers.Dense(int(cfg["dense_units"]), activation="relu")(x)
        x = K.layers.Dropout(float(cfg["dense_dropout"]))(x)
        out = K.layers.Dense(1, activation="sigmoid")(x)

    elif model_id == "bilstm_attention":
        x = K.layers.Embedding(21, int(cfg["embedding_dim"]), mask_zero=True)(tokens)
        x = K.layers.Bidirectional(K.layers.LSTM(
            int(cfg["bilstm_units_per_direction"]), return_sequences=True,
            dropout=float(cfg["lstm_dropout"]), recurrent_dropout=float(cfg["recurrent_dropout"]),
        ))(x)
        x = MaskedAdditiveAttention(int(cfg["attention_units"]))([x, mask])
        x = K.layers.Dense(int(cfg["dense_units"]), activation="relu")(x)
        x = K.layers.Dropout(float(cfg["dense_dropout"]))(x)
        out = K.layers.Dense(1, activation="sigmoid")(x)

    elif model_id == "transformer_encoder":
        d = int(cfg["model_dimension"])
        tok = K.layers.Embedding(21, d, mask_zero=True)(tokens)
        x = AddLearnedPositionEmbedding(200, d, name="position_embedding")(tok)
        attn_mask = K.layers.Lambda(lambda m: tf.expand_dims(m, axis=1), name="attention_mask")(mask)
        for block in range(int(cfg["n_encoder_blocks"])):
            attn = K.layers.MultiHeadAttention(
                num_heads=int(cfg["num_attention_heads"]),
                key_dim=int(cfg["key_dim_per_head"]),
                dropout=float(cfg["attention_dropout"]),
                name=f"mha_{block+1}",
            )(x, x, attention_mask=attn_mask)
            x = K.layers.LayerNormalization(name=f"attn_ln_{block+1}")(x + attn)
            ff = K.layers.Dense(int(cfg["feed_forward_dimension"]), activation="relu")(x)
            ff = K.layers.Dropout(float(cfg["feed_forward_dropout"]))(ff)
            ff = K.layers.Dense(d)(ff)
            x = K.layers.LayerNormalization(name=f"ff_ln_{block+1}")(x + ff)
        x = MaskedGlobalAveragePooling1D()([x, mask])
        x = K.layers.Dense(int(cfg["dense_units"]), activation="relu")(x)
        x = K.layers.Dropout(float(cfg["dense_dropout"]))(x)
        out = K.layers.Dense(1, activation="sigmoid")(x)
    else:
        raise KeyError(model_id)

    model = K.Model(tokens, out, name=model_id)
    model.compile(
        optimizer=K.optimizers.Adam(
            learning_rate=float(DEEP_SHARED["learning_rate"]),
            beta_1=float(DEEP_SHARED["beta_1"]), beta_2=float(DEEP_SHARED["beta_2"]),
            epsilon=float(DEEP_SHARED["epsilon"]), clipnorm=float(DEEP_SHARED["clipnorm"]),
        ),
        loss="binary_crossentropy",
    )
    return model


def deep_predict(model, X):
    p = np.asarray(model.predict(X, batch_size=256, verbose=0), dtype=float).reshape(-1)
    assert np.isfinite(p).all() and ((p >= 0.0) & (p <= 1.0)).all()
    return p


# ============================================================
# PUBLISHED MODEL HELPERS
# ============================================================

AMPEP_PROPERTY_GROUPS = [
    ("charge", ["ACFGHILMNPQSTVWY", "DE", "KR"]),
    ("hydrophobicity", ["CFILMVW", "AGHPSTY", "DEKNQR"]),
    ("normalized_vander_waals", ["ACDGPST", "EILNQV", "FHKMRWY"]),
    ("polarity", ["CFILMVWY", "AGPST", "DEHKNQR"]),
    ("polarizability", ["ADGST", "CEILNPQV", "FHKMRWY"]),
    ("secondary_structure", ["DGNPS", "AEHKLMQR", "CFITVWY"]),
    ("solvent_accessibility", ["ACFGILVW", "HMPSTY", "DEKNRQ"]),
]


def ampep_distribution_positions(pos, n):
    if not pos:
        return [0.0] * 5
    count = len(pos)
    out = [100.0 * pos[0] / n]
    for frac in [0.25, 0.50, 0.75]:
        k = math.floor(count * frac)
        selected = pos[0] if k == 0 else pos[k - 1]
        out.append(100.0 * selected / n)
    out.append(100.0 * pos[-1] / n)
    return out


def ampep_feature_vector(seq: str):
    seq = str(seq).strip().upper()
    n = len(seq)
    vals = []
    for _, groups in AMPEP_PROPERTY_GROUPS:
        for group in groups:
            positions = [i for i, aa in enumerate(seq, start=1) if aa in group]
            vals.extend(ampep_distribution_positions(positions, n))
    out = np.asarray(vals, dtype=np.float64)
    assert out.shape == (105,)
    return out


def ampep_matrix(df):
    return np.vstack([ampep_feature_vector(s) for s in df["sequence"].astype(str)])


def run_ampep(train_df, eval_df, seed, run_dir):
    t0 = time.perf_counter()
    Xtr = ampep_matrix(train_df)
    ytr = train_df["label"].to_numpy(np.int64)
    Xev = ampep_matrix(eval_df)
    model = RandomForestClassifier(
        n_estimators=100, criterion="gini", max_depth=None,
        min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0,
        max_features=11, max_leaf_nodes=None, min_impurity_decrease=0.0,
        bootstrap=True, oob_score=False, n_jobs=1, random_state=int(seed),
        verbose=0, warm_start=False, class_weight=None, ccp_alpha=0.0,
        max_samples=None,
    )
    model.fit(Xtr, ytr)
    train_seconds = time.perf_counter() - t0
    checkpoint = run_dir / "checkpoint_ampep.pkl"
    with open(checkpoint, "wb") as f:
        pickle.dump(model, f, protocol=pickle.HIGHEST_PROTOCOL)
    score = model.predict_proba(Xev)[:, int(np.where(model.classes_ == 1)[0][0])]
    write_na_history(run_dir / "history.csv", "AmPEP non-iterative RF; history not applicable")
    return score, checkpoint, train_seconds


AMPIR_R = r'''args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 5L)
train_path <- args[1]; eval_path <- args[2]; model_path <- args[3]; pred_path <- args[4]; run_seed <- as.integer(args[5])
suppressPackageStartupMessages({library(ampir); library(caret); library(kernlab)})
predictors <- c("Amphiphilicity","Hydrophobicity","pI","Mw","Charge","Xc1.A","Xc1.R","Xc1.N","Xc1.D","Xc1.C","Xc1.E","Xc1.Q","Xc1.G","Xc1.H","Xc1.I","Xc1.L","Xc1.K","Xc1.M","Xc1.F","Xc1.P","Xc1.S","Xc1.T","Xc1.W","Xc1.Y","Xc1.V","Xc2.lambda.1","Xc2.lambda.2")
tr <- read.csv(train_path, stringsAsFactors=FALSE, check.names=FALSE)
ev <- read.csv(eval_path, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("seq_name","seq_aa","Label") %in% colnames(tr)))
tr$Label <- factor(tr$Label, levels=c("Bg","Tg")); ev$Label <- factor(ev$Label, levels=c("Bg","Tg"))
trf <- calculate_features(tr[,c("seq_name","seq_aa")], min_len=10)
evf <- calculate_features(ev[,c("seq_name","seq_aa")], min_len=10)
stopifnot(all(predictors %in% colnames(trf)), all(predictors %in% colnames(evf)))
frame <- trf[,predictors,drop=FALSE]; frame$Label <- tr$Label
tab <- table(tr$Label)
weights <- ifelse(tr$Label == "Tg", (1/as.numeric(tab[1]))*0.5, (1/as.numeric(tab[2]))*0.5)
ctrl <- trainControl(method="none", classProbs=TRUE)
set.seed(run_seed)
model <- train(Label ~ ., data=frame, method="svmRadial", trControl=ctrl, preProcess=c("center","scale"), weights=weights, tuneGrid=data.frame(sigma=0.07,C=1))
prob <- predict(model, newdata=evf[,predictors,drop=FALSE], type="prob")
stopifnot(identical(colnames(prob), c("Bg","Tg")))
saveRDS(model, model_path)
write.csv(data.frame(sample_id=ev$seq_name, y_score=as.numeric(prob[,"Tg"]), stringsAsFactors=FALSE), pred_path, row.names=FALSE)
'''


def run_ampir(train_df, eval_df, seed, run_dir):
    helper = HELPER_ROOT / "ampir_production.R"
    lock_helper(helper, AMPIR_R)
    train_path = run_dir / "ampir_train.csv"
    eval_path = run_dir / "ampir_eval.csv"
    model_path = run_dir / "checkpoint_ampir.rds"
    pred_path = run_dir / "ampir_all_predictions.csv"
    tr = pd.DataFrame({
        "seq_name": train_df["sample_id"].astype(str),
        "seq_aa": train_df["sequence"].astype(str),
        "Label": np.where(train_df["label"].astype(int).eq(1), "Tg", "Bg"),
    })
    ev = pd.DataFrame({
        "seq_name": eval_df["sample_id"].astype(str),
        "seq_aa": eval_df["sequence"].astype(str),
        "Label": np.where(eval_df["label"].astype(int).eq(1), "Tg", "Bg"),
    })
    tr.to_csv(train_path, index=False); ev.to_csv(eval_path, index=False)
    conda = shutil.which("conda") or str(Path.home() / "miniconda3/bin/conda")
    assert Path(conda).is_file()
    t0 = time.perf_counter()
    run_cmd([conda, "run", "-n", "genpept_r_models", "Rscript", "--vanilla", helper,
             train_path, eval_path, model_path, pred_path, str(int(seed))],
            timeout=86400, label="ampir production fit")
    seconds = time.perf_counter() - t0
    pred = pd.read_csv(pred_path)
    score_map = pred.set_index(pred["sample_id"].astype(str))["y_score"]
    scores = score_map.loc[eval_df["sample_id"].astype(str)].to_numpy(float)
    write_na_history(run_dir / "history.csv", "ampir fixed svmRadial fit; epoch history not applicable")
    return scores, model_path, seconds


AMPGRAM_R = r'''args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 7L)
train_path <- args[1]; eval_path <- args[2]; model_path <- args[3]; pred_path <- args[4]; selected_path <- args[5]; registry_path <- args[6]; run_seed <- as.integer(args[7])
suppressPackageStartupMessages({library(slam); library(biogram); library(ranger)})
stopifnot(as.character(getRversion()) == "3.6.2")
AA20 <- strsplit("ARNDCEQGHILKMFPSTWYV", "", fixed=TRUE)[[1]]
STAT_NAMES <- c("fraction_true","pred_mean","pred_median","n_peptide","n_pos","pred_min","pred_max","longest_pos","n_pos_10","frac_0_0.2","frac_0.2_0.4","frac_0.4_0.6","frac_0.6_0.8","frac_0.8_1")
make_mer <- function(input_df, training_mode=FALSE) {
  pieces <- lapply(seq_len(nrow(input_df)), function(i) {
    chars <- strsplit(input_df$sequence[i], "", fixed=TRUE)[[1]]; n_mers <- length(chars)-9L; stopifnot(n_mers>=1L)
    mer_matrix <- t(vapply(seq_len(n_mers), function(start) chars[start:(start+9L)], character(10L))); colnames(mer_matrix) <- paste0("X",seq_len(10L))
    g <- NA_character_; if(training_mode){if(input_df$length[i]>=11L && input_df$length[i]<=19L) g <- "[11,19]" else if(input_df$length[i]>19L && input_df$length[i]<=26L) g <- "(19,26]"}
    data.frame(mer_matrix, source_peptide=rep(input_df$sample_id[i],n_mers), mer_id=paste0(input_df$sample_id[i],"m",seq_len(n_mers)), target=rep(as.logical(input_df$label[i]==1L),n_mers), peptide_length=rep(as.integer(input_df$length[i]),n_mers), group=rep(g,n_mers), stringsAsFactors=FALSE, check.names=FALSE)
  }); out <- do.call(rbind,pieces); rownames(out)<-NULL; out
}
count_multi <- function(mer_df,ns,ds){m<-as.matrix(mer_df[,grep("^X",colnames(mer_df)),drop=FALSE]); biogram::binarize(biogram::count_multigrams(ns=ns,ds=ds,seq=m,u=AA20))}
select_features <- function(mer_df, bin, cutoff=0.05){keep<-mer_df$group %in% c("[11,19]","(19,26]"); test<-biogram::test_features(target=mer_df$target[keep],features=bin[keep,,drop=FALSE],criterion="ig",adjust="BH",threshold=1,quick=TRUE,times=1e5); cut(test,breaks=c(0,cutoff,1))[[1]]}
count_spec <- function(mer_df, imp){m<-as.matrix(mer_df[,grep("^X",colnames(mer_df)),drop=FALSE]); biogram::binarize(biogram::count_specified(m,imp))}
longest <- function(x){sp<-strsplit(paste0(as.numeric(x>0.5),collapse=""),"0")[[1]]; len<-unname(sapply(sp,nchar)); if(length(len[len>0])==0) 0 else len[len>0]}
stats <- function(pm){ids<-unique(as.character(pm$source_peptide)); rows<-lapply(ids,function(id){d<-pm[pm$source_peptide==id,,drop=FALSE]; p<-as.numeric(d$pred); r<-longest(p); n<-length(p); data.frame(source_peptide=id,target=as.logical(d$target[1]),fraction_true=mean(p>0.5),pred_mean=mean(p),pred_median=median(p),n_peptide=n,n_pos=sum(p>0.5),pred_min=min(p),pred_max=max(p),longest_pos=max(r),n_pos_10=sum(r>=10),frac_0_0.2=sum(p<=0.2)/n,frac_0.2_0.4=sum(p>0.2&p<=0.4)/n,frac_0.4_0.6=sum(p>0.4&p<=0.6)/n,frac_0.6_0.8=sum(p>0.6&p<=0.8)/n,frac_0.8_1=sum(p>0.8&p<=1)/n,stringsAsFactors=FALSE,check.names=FALSE)}); out<-do.call(rbind,rows); rownames(out)<-NULL; out$target<-factor(out$target); out}
tr <- read.csv(train_path,stringsAsFactors=FALSE,check.names=FALSE); ev <- read.csv(eval_path,stringsAsFactors=FALSE,check.names=FALSE); reg <- read.csv(registry_path,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(all(c("sample_id","sequence","label","length") %in% colnames(tr)), all(c("sample_id","sequence","label","length") %in% colnames(ev)), nrow(reg)==33620L)
set.seed(run_seed)
trm <- make_mer(tr,TRUE); evm <- make_mer(ev,FALSE)
a <- count_multi(trm,c(1,rep(2,4)),list(0,0,1,2,3)); b<-count_multi(trm,c(3,3),list(c(0,0),c(0,1))); c<-count_multi(trm,c(3,3),list(c(1,0),c(1,1))); bin<-cbind(a,b,c)
stopifnot(ncol(bin)==33620L, setequal(colnames(bin),reg$feature_name))
imp <- select_features(trm,bin,0.05); stopifnot(length(imp)>=2L); write.csv(data.frame(feature_name=imp),selected_path,row.names=FALSE)
keep<-trm$group %in% c("[11,19]","(19,26]"); rf1data<-data.frame(as.matrix(bin[keep,imp,drop=FALSE]),tar=as.factor(trm$target[keep]))
rf1<-ranger::ranger(dependent.variable.name="tar",data=rf1data,write.forest=TRUE,probability=TRUE,num.trees=2000,verbose=FALSE,seed=run_seed)
ptr<-predict(rf1,data.frame(as.matrix(bin[,imp,drop=FALSE])))[["predictions"]][,"TRUE"]; tmp<-trm; tmp$pred<-ptr; st<-stats(tmp)
rf2<-ranger::ranger(dependent.variable.name="target",data=st[,c("target",STAT_NAMES),drop=FALSE],write.forest=TRUE,probability=TRUE,num.trees=500,verbose=FALSE,classification=TRUE,seed=run_seed)
model<-list(rf_mers=rf1,rf_peptides=rf2,imp_features=imp); class(model)<-"ag_model"; saveRDS(model,model_path)
evbin<-count_spec(evm,imp); pev<-predict(rf1,data.frame(as.matrix(evbin)))[["predictions"]][,"TRUE"]; tmp2<-evm; tmp2$pred<-pev; sev<-stats(tmp2); score<-predict(rf2,data.frame(sev[,STAT_NAMES,drop=FALSE]))[["predictions"]][,"TRUE"]
write.csv(data.frame(sample_id=as.character(sev$source_peptide),y_score=as.numeric(score),stringsAsFactors=FALSE),pred_path,row.names=FALSE)
'''


def run_ampgram(train_df, eval_df, seed, run_dir):
    helper = HELPER_ROOT / "ampgram_production.R"
    lock_helper(helper, AMPGRAM_R)
    train_path = run_dir / "ampgram_train.csv"; eval_path = run_dir / "ampgram_eval.csv"
    model_path = run_dir / "checkpoint_ampgram.rds"; pred_path = run_dir / "ampgram_all_predictions.csv"
    selected = run_dir / "ampgram_selected_features.csv"
    cols = ["sample_id", "sequence", "label", "length"]
    train_df[cols].to_csv(train_path, index=False); eval_df[cols].to_csv(eval_path, index=False)
    cell23 = load_json(MANIFEST_DIR / "benchmark17_ampgram_FINAL_exact_R_environment_runtime_v5.json")["contract"]
    rscript = Path(cell23["runtime"]["Rscript"]); env_root = Path(cell23["runtime"]["environment"]); rlib = Path(cell23["runtime"]["library"])
    reg = AUDIT_DIR / "benchmark17_ampgram_33620_feature_registry.csv"
    env = os.environ.copy()
    for k in ["R_HOME","R_LIBS","R_LIBS_USER","R_ENVIRON","R_ENVIRON_USER","R_PROFILE","R_PROFILE_USER","R_MAKEVARS_USER"]:
        env.pop(k, None)
    env["PATH"] = str(env_root / "bin") + os.pathsep + env.get("PATH", "")
    env["R_LIBS_USER"] = str(rlib); env["R_LIBS_SITE"] = str(rlib)
    t0 = time.perf_counter()
    run_cmd([rscript, "--vanilla", helper, train_path, eval_path, model_path, pred_path, selected, reg, str(int(seed))], env=env, timeout=172800, label="AmpGram production fit")
    seconds = time.perf_counter() - t0
    pred = pd.read_csv(pred_path); score_map = pred.set_index(pred["sample_id"].astype(str))["y_score"]
    scores = score_map.loc[eval_df["sample_id"].astype(str)].to_numpy(float)
    write_na_history(run_dir / "history.csv", "AmpGram two-stage ranger recipe; epoch history not applicable")
    return scores, model_path, seconds


AMPSCANNER_PY = r'''from pathlib import Path
import csv, importlib.util, os, random, sys
os.environ["KERAS_BACKEND"]="tensorflow"; os.environ["CUDA_VISIBLE_DEVICES"]=""; os.environ["OMP_NUM_THREADS"]="1"; os.environ["MKL_NUM_THREADS"]="1"; os.environ["NUMEXPR_NUM_THREADS"]="1"; os.environ["TF_CPP_MIN_LOG_LEVEL"]="2"
import numpy as np, tensorflow as tf, keras
from Bio import SeqIO
from keras import backend as K
from keras.preprocessing import sequence
SOURCE=Path(sys.argv[1]); AMP=Path(sys.argv[2]); NON=Path(sys.argv[3]); VAMP=Path(sys.argv[4]); VNON=Path(sys.argv[5]); EVAL=Path(sys.argv[6]); OUT=Path(sys.argv[7]); SEED=int(sys.argv[8]); OUT.mkdir(parents=True,exist_ok=True)
assert sys.version_info[:2]==(3,6); assert tf.__version__=="1.2.1"; assert keras.__version__=="2.0.6"
K.set_session(tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=1,inter_op_parallelism_threads=1,device_count={"GPU":0})))
tf.set_random_seed(SEED); np.random.seed(SEED); random.seed(SEED)
spec=importlib.util.spec_from_file_location("amp_scanner_prod",str(SOURCE)); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
def seqs(path): return [(r.id,str(r.seq)) for r in SeqIO.parse(str(path),"fasta")]
def enc(records): return sequence.pad_sequences([[mod.aa2int[a] for a in s] for _,s in records],maxlen=200)
atr=seqs(AMP); ntr=seqs(NON); ava=seqs(VAMP); nva=seqs(VNON); ev=seqs(EVAL)
tr=atr+ntr; y=np.asarray([1]*len(atr)+[0]*len(ntr),dtype=np.int64); p=np.random.RandomState(SEED).permutation(len(y)); X=enc(tr)[p]; y=y[p]
vr=ava+nva; yv=np.asarray([1]*len(ava)+[0]*len(nva),dtype=np.int64); XV=enc(vr)
model_path=OUT/"checkpoint_ampscannerv2.h5"
model=mod.compile_model(X,y,XV,yv,saved_model_name=str(model_path),merge_train_and_val=False)
XE=enc(ev); score=np.asarray(model.predict(XE,batch_size=32,verbose=0)).reshape(-1)
with open(str(OUT/"all_predictions.csv"),"w",newline="") as fh:
 writer=csv.writer(fh); writer.writerow(["sample_id","y_score"])
 for (sample_id,_),value in zip(ev,score): writer.writerow([sample_id,repr(float(value))])
'''


def run_ampscannerv2(train_df, validation_df, eval_df, seed, run_dir):
    helper = HELPER_ROOT / "ampscannerv2_production_v2.py"
    lock_helper(helper, AMPSCANNER_PY)
    source = Path("/home/pc/genpept_sources/amp-scanner-v2_933052e2365631fe93098892120ee535e0ba381a") / "amp_scanner_v2_train_tf1.py"
    py = Path.home() / "miniconda3/envs/genpept_ampscannerv2_orig/bin/python"
    amp=run_dir/"train_amp.fa"; non=run_dir/"train_nonamp.fa"; va=run_dir/"validation_amp.fa"; vn=run_dir/"validation_nonamp.fa"; ev=run_dir/"all_eval.fa"
    write_class_fastas(train_df,amp,non); write_class_fastas(validation_df,va,vn); write_fasta(eval_df,ev)
    t0=time.perf_counter(); run_cmd([py,helper,source,amp,non,va,vn,ev,run_dir,str(int(seed))],timeout=172800,label="AMPScannerV2 production fit"); seconds=time.perf_counter()-t0
    pred=pd.read_csv(run_dir/"all_predictions.csv"); score_map=pred.set_index(pred["sample_id"].astype(str))["y_score"]; scores=score_map.loc[eval_df["sample_id"].astype(str)].to_numpy(float)
    write_na_history(run_dir/"history.csv","AMPScannerV2 official fixed 10-epoch source fit; history object not exposed")
    return scores, run_dir/"checkpoint_ampscannerv2.h5", seconds


AMPEPPY_PY = r'''from pathlib import Path
from types import SimpleNamespace
import os,pickle,sys
import pandas as pd, numpy as np, sklearn
from Bio import SeqIO
SOURCE=Path(sys.argv[1]); AMP=Path(sys.argv[2]); NON=Path(sys.argv[3]); EVAL=Path(sys.argv[4]); OUT=Path(sys.argv[5]); SEED=int(sys.argv[6]); OUT.mkdir(parents=True,exist_ok=True)
assert ".".join(str(x) for x in sys.version_info[:3])=="3.8.3"; assert sklearn.__version__=="0.23.1"
sys.path.insert(0,str(SOURCE)); from amPEPpy import amPEP as mod; from amPEPpy._version import __version__; assert __version__=="1.0"
cwd=Path.cwd(); os.chdir(str(OUT))
try:
 mod.train(SimpleNamespace(positive=str(AMP),negative=str(NON),drop_feature=None,max_tree=175,min_tree=23,seed=SEED,num_processes=1,num_trees=160,tree_test=False,feature_importance=False))
finally: os.chdir(str(cwd))
model=OUT/"amPEP.model"; assert model.is_file()
pred_raw=OUT/"source_predictions.tsv"
mod.predict(SimpleNamespace(model=str(model),seq_file=str(EVAL),out_file=str(pred_raw),drop_feature=None))
p=pd.read_csv(pred_raw,sep="\t")
# Official output names verified in Cell 28.
assert "probability_AMP" in p.columns and "seq_id" in p.columns
pd.DataFrame({"sample_id":p["seq_id"].astype(str),"y_score":p["probability_AMP"].astype(float)}).to_csv(OUT/"all_predictions.csv",index=False)
'''


def run_ampeppy(train_df, eval_df, seed, run_dir):
    helper=HELPER_ROOT/"ampeppy_production.py"
    lock_helper(helper, AMPEPPY_PY)
    source=Path("/home/pc/genpept_sources/amPEPpy_v1.0_aa1f694c6cb3d09b16bc9378bed77f59e4f1e780")
    py=Path.home()/"miniconda3/envs/genpept_ampeppy_v1/bin/python"
    amp=run_dir/"train_amp.fa"; non=run_dir/"train_nonamp.fa"; ev=run_dir/"all_eval.fa"; write_class_fastas(train_df,amp,non); write_fasta(eval_df,ev)
    t0=time.perf_counter(); run_cmd([py,helper,source,amp,non,ev,run_dir,str(int(seed))],timeout=172800,label="amPEPpy production fit"); seconds=time.perf_counter()-t0
    pred=pd.read_csv(run_dir/"all_predictions.csv"); score_map=pred.set_index(pred["sample_id"].astype(str))["y_score"]; scores=score_map.loc[eval_df["sample_id"].astype(str)].to_numpy(float)
    write_na_history(run_dir/"history.csv","amPEPpy v1.0 RandomForest fit; epoch history not applicable")
    return scores, run_dir/"amPEP.model", seconds


AI4AMP_PY = r'''from pathlib import Path
import importlib.util,os,random,sys
os.environ["KERAS_BACKEND"]="tensorflow"; os.environ["CUDA_VISIBLE_DEVICES"]=""; os.environ["OMP_NUM_THREADS"]="1"; os.environ["MKL_NUM_THREADS"]="1"; os.environ["OPENBLAS_NUM_THREADS"]="1"; os.environ["TF_CPP_MIN_LOG_LEVEL"]="2"
import numpy as np,pandas as pd,tensorflow as tf,keras
from keras import backend as K, optimizers
from keras.callbacks import ModelCheckpoint
from keras.layers import Input,Conv1D,LSTM,Dense
from keras.models import Model,load_model
PC6=Path(sys.argv[1]); TABLE=Path(sys.argv[2]); AMP=Path(sys.argv[3]); NON=Path(sys.argv[4]); VAMP=Path(sys.argv[5]); VNON=Path(sys.argv[6]); EVAL=Path(sys.argv[7]); OUT=Path(sys.argv[8]); SEED=int(sys.argv[9]); OUT.mkdir(parents=True,exist_ok=True)
assert ".".join(str(x) for x in sys.version_info[:3])=="3.6.9"; assert tf.__version__=="1.9.0"; assert keras.__version__=="2.2.4"
random.seed(SEED); np.random.seed(SEED); tf.set_random_seed(SEED); K.set_session(tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=1,inter_op_parallelism_threads=1,device_count={"GPU":0})))
spec=importlib.util.spec_from_file_location("pc6",str(PC6)); pc6=importlib.util.module_from_spec(spec); spec.loader.exec_module(pc6)
def enc(path):
 d=pc6.PC_6(str(path),length=200); return list(d.keys()),np.asarray(list(d.values()),dtype=np.float32)
ia,xa=enc(AMP); inn,xn=enc(NON); iva,xva=enc(VAMP); ivn,xvn=enc(VNON); ide,xe=enc(EVAL)
X=np.concatenate([xa,xn]); y=np.concatenate([np.ones(len(xa)),np.zeros(len(xn))]).astype(np.float32); p=np.random.RandomState(SEED).permutation(len(y)); X=X[p]; y=y[p]
XV=np.concatenate([xva,xvn]); yv=np.concatenate([np.ones(len(xva)),np.zeros(len(xvn))]).astype(np.float32)
inp=Input(shape=(200,6)); x=Conv1D(64,16,activation="relu",padding="same")(inp); x=LSTM(100,return_sequences=False)(x); out=Dense(1,activation="sigmoid")(x); model=Model(inp,out)
model.compile(optimizer=optimizers.Adam(lr=0.0003),loss="binary_crossentropy",metrics=["accuracy"])
batch=max(1,int(0.5*len(X))); best=OUT/"checkpoint_ai4amp.h5"; hist=model.fit(X,y,validation_data=(XV,yv),shuffle=True,epochs=200,batch_size=batch,callbacks=[ModelCheckpoint(str(best),monitor="val_loss",save_best_only=True,mode="min",verbose=0)],verbose=0)
pd.DataFrame({"epoch":np.arange(1,len(hist.history["loss"])+1),"loss":hist.history["loss"],"val_loss":hist.history["val_loss"]}).to_csv(OUT/"history.csv",index=False)
bm=load_model(str(best)); score=np.asarray(bm.predict(xe,batch_size=batch,verbose=0)).reshape(-1); pd.DataFrame({"sample_id":ide,"y_score":score}).to_csv(OUT/"all_predictions.csv",index=False)
'''


def run_ai4amp(train_df, validation_df, eval_df, seed, run_dir):
    helper=HELPER_ROOT/"ai4amp_production.py"
    lock_helper(helper, AI4AMP_PY)
    pc6dir=Path("/home/pc/genpept_sources/PC6-protein-encoding-method_15f9e9997d26f1dbb6863d9c19881b5cf49f0533")
    py=Path.home()/"miniconda3/envs/genpept_ai4amp_legacy/bin/python"
    amp=run_dir/"train_amp.fa"; non=run_dir/"train_nonamp.fa"; va=run_dir/"validation_amp.fa"; vn=run_dir/"validation_nonamp.fa"; ev=run_dir/"all_eval.fa"; write_class_fastas(train_df,amp,non); write_class_fastas(validation_df,va,vn); write_fasta(eval_df,ev)
    t0=time.perf_counter(); run_cmd([py,helper,pc6dir/"Protein_Encoding.py",pc6dir/"6-pc",amp,non,va,vn,ev,run_dir,str(int(seed))],timeout=259200,label="AI4AMP production fit"); seconds=time.perf_counter()-t0
    pred=pd.read_csv(run_dir/"all_predictions.csv"); score_map=pred.set_index(pred["sample_id"].astype(str))["y_score"]; scores=score_map.loc[eval_df["sample_id"].astype(str)].to_numpy(float)
    return scores, run_dir/"checkpoint_ai4amp.h5", seconds


AMPLIFY_PY = r'''from pathlib import Path
import importlib.util,os,random,sys,tarfile
os.environ["KERAS_BACKEND"]="tensorflow"; os.environ["CUDA_VISIBLE_DEVICES"]=""; os.environ["OMP_NUM_THREADS"]="1"; os.environ["MKL_NUM_THREADS"]="1"; os.environ["OPENBLAS_NUM_THREADS"]="1"; os.environ["NUMEXPR_NUM_THREADS"]="1"; os.environ["TF_CPP_MIN_LOG_LEVEL"]="2"
import numpy as np,pandas as pd,tensorflow as tf,keras
from Bio import SeqIO
from keras import backend as K
from keras.callbacks import EarlyStopping
from sklearn.model_selection import StratifiedKFold
SOURCE=Path(sys.argv[1]); AMP=Path(sys.argv[2]); NON=Path(sys.argv[3]); EVAL=Path(sys.argv[4]); OUT=Path(sys.argv[5]); SEED=int(sys.argv[6]); OUT.mkdir(parents=True,exist_ok=True); SRC=SOURCE/"src"; sys.path.insert(0,str(SRC))
assert ".".join(str(x) for x in sys.version_info[:3])=="3.6.7"; assert tf.__version__=="1.12.0"; assert keras.__version__=="2.2.4"
random.seed(SEED); np.random.seed(SEED); tf.set_random_seed(SEED); K.set_session(tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=1,inter_op_parallelism_threads=1,device_count={"GPU":0})))
def imp(name,path):
 s=importlib.util.spec_from_file_location(name,str(path)); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
trm=imp("amplify_train_prod",SRC/"train_amplify.py"); prm=imp("amplify_pred_prod",SRC/"AMPlify.py")
def read(path): return [(r.id,str(r.seq)) for r in SeqIO.parse(str(path),"fasta")]
a=read(AMP); n=read(NON); e=read(EVAL); pairs=[(s,1) for _,s in a]+[(s,0) for _,s in n]; random.Random(123).shuffle(pairs); seqs,ys=zip(*pairs); y=np.asarray(ys,dtype=np.int64); X=trm.one_hot_padding(list(seqs),200).astype(np.float32); XE=prm.one_hot_padding([s for _,s in e],200).astype(np.float32)
folds=list(StratifiedKFold(n_splits=5,shuffle=True,random_state=50).split(X,y)); preds=[]; hrows=[]; weights=[]
for fold,(ti,vi) in enumerate(folds,1):
 model=trm.build_model(); cb=EarlyStopping(monitor="val_acc",min_delta=0.001,patience=50,restore_best_weights=True); h=model.fit(X[ti],y[ti],epochs=1000,batch_size=32,validation_data=(X[vi],y[vi]),verbose=0,callbacks=[cb]); wp=OUT/("amplify_fold_%d.h5"%fold); model.save_weights(str(wp)); weights.append(wp); preds.append(np.asarray(model.predict(XE,batch_size=32,verbose=0)).reshape(-1));
 for i,loss in enumerate(h.history.get("loss",[])): hrows.append({"fold":fold,"epoch":i+1,"loss":loss,"val_loss":h.history.get("val_loss",[float("nan")]*len(h.history.get("loss",[])))[i]})
score=np.mean(np.vstack(preds),axis=0); pd.DataFrame({"sample_id":[i for i,_ in e],"y_score":score}).to_csv(OUT/"all_predictions.csv",index=False); pd.DataFrame(hrows).to_csv(OUT/"amplify_fold_history.csv",index=False)
with tarfile.open(str(OUT/"checkpoint_amplify_5fold_weights.tar.gz"),"w:gz") as tar:
 for p in weights: tar.add(str(p),arcname=p.name)
'''


def run_amplify(train_df, eval_df, seed, run_dir):
    helper=HELPER_ROOT/"amplify_production.py"
    lock_helper(helper, AMPLIFY_PY)
    source=Path("/home/pc/genpept_sources/AMPlify_v1.0.0_ba1f6d196873876cd53d2afb8a71938f6d1687d1")
    py=Path.home()/"miniconda3/envs/genpept_amplify_v1/bin/python"
    amp=run_dir/"train_amp.fa"; non=run_dir/"train_nonamp.fa"; ev=run_dir/"all_eval.fa"; write_class_fastas(train_df,amp,non); write_fasta(eval_df,ev)
    t0=time.perf_counter(); run_cmd([py,helper,source,amp,non,ev,run_dir,str(int(seed))],timeout=604800,label="AMPlify production fit"); seconds=time.perf_counter()-t0
    pred=pd.read_csv(run_dir/"all_predictions.csv"); score_map=pred.set_index(pred["sample_id"].astype(str))["y_score"]; scores=score_map.loc[eval_df["sample_id"].astype(str)].to_numpy(float)
    # Canonical run-level history keeps required fields; detailed fold history is preserved separately.
    foldh=pd.read_csv(run_dir/"amplify_fold_history.csv")
    hist=foldh.groupby("epoch",as_index=False)[["loss","val_loss"]].mean(numeric_only=True)
    atomic_write_csv(run_dir/"history.csv",hist[["epoch","loss","val_loss"]])
    return scores, run_dir/"checkpoint_amplify_5fold_weights.tar.gz", seconds


# ============================================================
# LOCK PRODUCTION HELPER SOURCES BEFORE ANY BENCHMARK TRAINING
# ============================================================

HELPER_SOURCES = {
    "ampir": (HELPER_ROOT / "ampir_production.R", AMPIR_R),
    "ampgram": (HELPER_ROOT / "ampgram_production.R", AMPGRAM_R),
    "ampscannerv2": (HELPER_ROOT / "ampscannerv2_production_v2.py", AMPSCANNER_PY),
    "ampeppy": (HELPER_ROOT / "ampeppy_production.py", AMPEPPY_PY),
    "ai4amp": (HELPER_ROOT / "ai4amp_production.py", AI4AMP_PY),
    "amplify": (HELPER_ROOT / "amplify_production.py", AMPLIFY_PY),
}
HELPER_HASHES = {}
for _method, (_path, _text) in HELPER_SOURCES.items():
    HELPER_HASHES[_method] = lock_helper(_path, _text)
    if _path.suffix == ".py":
        compile(_text, str(_path), "exec")


# ============================================================
# RUN EXECUTION
# ============================================================
# ============================================================
# RUN EXECUTION — TRAIN + SOURCE-VALIDATION ONLY
#
# SCIENTIFIC BOUNDARY
# -------------------
# Cell 32 MUST NOT read source TEST rows and MUST NOT read
# external target rows. It creates frozen source checkpoints,
# validation predictions, and source-validation thresholds only.
# Evaluation is a later Cell.
# ============================================================


def trained_manifest_valid(row) -> bool:
    run_dir = Path(str(row.run_dir))
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        return False
    try:
        manifest = load_json(manifest_path)
        if manifest.get("status") not in {"TRAINED", "EVALUATED", "COMPLETE"}:
            return False
        if manifest.get("run_id") != str(row.run_id):
            return False
        if manifest.get("model_id") != str(row.model_id):
            return False
        if manifest.get("source_dataset") != str(row.source_dataset):
            return False
        if int(manifest.get("seed")) != int(row.seed):
            return False
        if str(manifest.get("train_index_key")) != str(row.train_index_key):
            return False
        if str(manifest.get("validation_index_key")) != str(row.validation_index_key):
            return False
        if str(manifest.get("source_test_index_key")) != str(row.source_test_index_key):
            return False
        if "training_artifacts" not in manifest:
            return False
        if (run_dir / "source_test_predictions.csv").exists():
            return False
        if (run_dir / "source_test_metrics.json").exists():
            return False
        if (run_dir / "external").exists():
            return False
        required = manifest.get("training_artifacts", {})
        required_names = {
            "config.json",
            "validation_predictions.csv",
            "validation_threshold.json",
            "history.csv",
            "resources.json",
        }
        assert required_names.issubset(set(required))
        for relpath, entry in required.items():
            p = run_dir / relpath
            if not p.is_file():
                return False
            if int(p.stat().st_size) != int(entry["size_bytes"]):
                return False
            if sha256_file(p) != str(entry["sha256"]):
                return False
        checkpoint_rel = str(manifest.get("checkpoint_relative_path", ""))
        if not checkpoint_rel:
            return False
        checkpoint = run_dir / checkpoint_rel
        if not checkpoint.is_file():
            return False
        cp_entry = manifest.get("checkpoint_artifact", {})
        if int(checkpoint.stat().st_size) != int(cp_entry.get("size_bytes", -1)):
            return False
        if sha256_file(checkpoint) != str(cp_entry.get("sha256", "")):
            return False
        return True
    except Exception:
        return False


def update_ledger(run_id: str, **updates):
    ledger_now = read_csv_exact(RUN_LEDGER_PATH)
    mask = ledger_now["run_id"].astype(str).eq(str(run_id))
    assert int(mask.sum()) == 1
    identity_before = ledger_now.loc[mask, [
        "job_no", "run_id", "model_id", "model_group", "source_dataset", "seed"
    ]].copy()
    for key, value in updates.items():
        assert key in ledger_now.columns, f"Unknown ledger field: {key}"
        ledger_now.loc[mask, key] = value
    identity_after = ledger_now.loc[mask, [
        "job_no", "run_id", "model_id", "model_group", "source_dataset", "seed"
    ]].copy()
    pd.testing.assert_frame_equal(
        identity_before.reset_index(drop=True),
        identity_after.reset_index(drop=True),
        check_dtype=False,
    )
    atomic_write_csv(RUN_LEDGER_PATH, ledger_now)
    reread = read_csv_exact(RUN_LEDGER_PATH)
    assert len(reread) == 680
    assert reread["run_id"].is_unique


def training_artifact_entry(path: Path):
    path = Path(path)
    assert path.is_file()
    return {
        "size_bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
    }


def prepare_training_frames(run_row):
    # CRITICAL: only source TRAIN and source VALIDATION are materialized.
    train_df = df_from_index_key(str(run_row.train_index_key))
    validation_df = df_from_index_key(str(run_row.validation_index_key))

    assert len(train_df) == int(run_row.train_n)
    assert len(validation_df) == int(run_row.validation_n)
    assert set(train_df["dataset"].astype(str)) == {str(run_row.source_dataset)}
    assert set(validation_df["dataset"].astype(str)) == {str(run_row.source_dataset)}
    assert set(train_df["split"].astype(str)) == {"train"}
    assert set(validation_df["split"].astype(str)) == {"validation"}
    assert set(train_df["sample_id"].astype(str)).isdisjoint(
        set(validation_df["sample_id"].astype(str))
    )
    return train_df, validation_df


def execute_training_run(run_row):
    run_wall_start = time.perf_counter()
    run_id = str(run_row.run_id)
    model_id = str(run_row.model_id)
    source_dataset = str(run_row.source_dataset)
    seed = int(run_row.seed)
    final_run_dir = Path(str(run_row.run_dir))

    assert model_id in ALL_IDS
    assert source_dataset in EXPECTED_DATASETS
    assert seed in EXPECTED_SEEDS
    assert not final_run_dir.exists(), (
        "Final run directory already exists without a validated TRAINED/EVALUATED/COMPLETE manifest.\n"
        "Refusing silent overwrite:\n"
        f"{final_run_dir}"
    )

    final_run_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(tempfile.mkdtemp(
        prefix=final_run_dir.name + ".cell32_staging.",
        dir=str(final_run_dir.parent),
    ))

    try:
        frame_load_start = time.perf_counter()
        train_df, validation_df = prepare_training_frames(run_row)
        frame_load_seconds = time.perf_counter() - frame_load_start
        set_master_seed(seed)

        # Verify closure binding for this exact row before fitting.
        closure_path = Path(str(run_row.closure_manifest))
        assert closure_path.is_file()
        if model_id in CONTROLLED_IDS:
            assert sha256_file(closure_path) == str(run_row.closure_sha256)
        else:
            selftest_path = Path(str(run_row.selftest_manifest))
            assert selftest_path.is_file()
            observed_selftest_sha = sha256_file(selftest_path)
            assert observed_selftest_sha == str(run_row.selftest_sha256)
            expected_bundle_sha = sha256_text(
                sha256_file(PUBLISHED_ROUTE_CONTRACT_PATH)
                + "\n"
                + observed_selftest_sha
                + "\n"
            )
            assert expected_bundle_sha == str(run_row.closure_sha256)

        config_payload = {
            "run_id": run_id,
            "model_id": model_id,
            "source_dataset": source_dataset,
            "seed": seed,
            "implementation_route": str(run_row.implementation_route),
            "implementation_environment": str(run_row.implementation_environment),
            "train_index_key": str(run_row.train_index_key),
            "validation_index_key": str(run_row.validation_index_key),
            "source_test_index_key": str(run_row.source_test_index_key),
            "cell32_trainer_version": "v9_adaptive_throughput_timing_exact_route_recovery",
            "cell32_trainer_sha256": CELL32_TRAINER_SHA256 or None,
            "technical_fix": "v4 mask dtype correction retained; v9 retains exact scientific training and adds adaptive throughput scheduling/timing/recovery and removes an unnecessary pandas dependency from AMPScannerV2 CSV emission by using Python stdlib csv; no architecture/hyperparameter/data/protocol change",
            "cell32_scope": {
                "fit_rows": "SOURCE TRAIN ONLY",
                "selection_rows": "SOURCE VALIDATION ONLY",
                "source_test_read": False,
                "external_target_read": False,
                "source_test_evaluation": False,
                "protocol2_evaluation": False,
            },
            "closure_manifest": str(closure_path),
            "closure_sha256": str(run_row.closure_sha256),
        }
        if model_id in CONTROLLED_IDS:
            config_payload["locked_model_config"] = MODEL_CONFIGS[model_id]
        else:
            config_payload["locked_published_route"] = route_wrapper["contract"]["methods"][model_id]
            config_payload["selftest_manifest"] = str(run_row.selftest_manifest)
            config_payload["selftest_sha256"] = str(run_row.selftest_sha256)
            if model_id in HELPER_HASHES:
                config_payload["production_helper_sha256"] = HELPER_HASHES[model_id]

        atomic_write_json(staging_dir / "config.json", config_payload)

        fit_start = time.perf_counter()
        checkpoint = None
        model_parameter_count = None
        gpu_peak_override = None
        feature_or_token_load_seconds = None
        preprocessing_seconds = None
        timing_provenance = "not_separately_measured"

        if model_id in CONTROLLED_IDS[:6]:
            model = make_classical_model(model_id, seed)
            feature_load_start = time.perf_counter()
            train_idx = np.asarray(INDEX[str(run_row.train_index_key)], dtype=np.int64)
            val_idx = np.asarray(INDEX[str(run_row.validation_index_key)], dtype=np.int64)
            X_train = np.asarray(X_CLASSICAL[train_idx], dtype=np.float32)
            y_train = np.asarray(Y[train_idx], dtype=np.int64)
            X_val = np.asarray(X_CLASSICAL[val_idx], dtype=np.float32)
            feature_or_token_load_seconds = time.perf_counter() - feature_load_start
            preprocessing_seconds = 0.0
            timing_provenance = "cached_421D_features_loaded_per_run; no run-specific preprocessing"

            assert X_train.shape[0] == len(train_df)
            assert X_val.shape[0] == len(validation_df)
            model.fit(X_train, y_train)
            training_seconds = time.perf_counter() - fit_start
            validation_start = time.perf_counter()
            val_score = positive_probability(model, X_val)
            validation_seconds = time.perf_counter() - validation_start

            checkpoint = staging_dir / "checkpoint_controlled.pkl"
            with open(checkpoint, "wb") as f:
                pickle.dump(model, f, protocol=pickle.HIGHEST_PROTOCOL)
            write_na_history(
                staging_dir / "history.csv",
                f"{model_id} non-iterative estimator; epoch history not applicable",
            )

        elif model_id in CONTROLLED_IDS[6:]:
            (
                val_score,
                checkpoint,
                training_seconds,
                validation_seconds,
                model_parameter_count,
                gpu_peak_override,
                feature_or_token_load_seconds,
                preprocessing_seconds,
            ) = run_controlled_deep_subprocess(
                model_id, seed, run_row, staging_dir
            )
            timing_provenance = "deep worker measured mmap/index load and train/validation array materialization separately"

        else:
            # Published adapters retain their source-faithful wrappers. Several perform
            # feature extraction inside the same legacy subprocess/function as fitting;
            # v6 does not fabricate a separate preprocessing time for those routes.
            timing_provenance = "published_route_total_fit_wrapper_time; preprocessing not safely separable without changing the locked published route"
            validation_start = time.perf_counter()
            if model_id == "ampep":
                val_score, checkpoint, training_seconds = run_ampep(
                    train_df, validation_df, seed, staging_dir
                )
            elif model_id == "ampir":
                val_score, checkpoint, training_seconds = run_ampir(
                    train_df, validation_df, seed, staging_dir
                )
            elif model_id == "ampgram":
                val_score, checkpoint, training_seconds = run_ampgram(
                    train_df, validation_df, seed, staging_dir
                )
            elif model_id == "ampscannerv2":
                val_score, checkpoint, training_seconds = run_ampscannerv2(
                    train_df, validation_df, validation_df, seed, staging_dir
                )
            elif model_id == "ampeppy":
                val_score, checkpoint, training_seconds = run_ampeppy(
                    train_df, validation_df, seed, staging_dir
                )
            elif model_id == "ai4amp":
                val_score, checkpoint, training_seconds = run_ai4amp(
                    train_df, validation_df, validation_df, seed, staging_dir
                )
            elif model_id == "amplify":
                val_score, checkpoint, training_seconds = run_amplify(
                    train_df, validation_df, seed, staging_dir
                )
            else:
                raise KeyError(model_id)
            validation_seconds = max(
                0.0,
                time.perf_counter() - validation_start - float(training_seconds),
            )

        val_score = np.asarray(val_score, dtype=float).reshape(-1)
        assert val_score.shape == (len(validation_df),)
        assert np.isfinite(val_score).all()
        assert ((val_score >= 0.0) & (val_score <= 1.0)).all()

        threshold_info = select_validation_threshold(
            validation_df["label"].to_numpy(dtype=int),
            val_score,
        )
        tau = float(threshold_info["tau"])

        threshold_payload = {
            "run_id": run_id,
            "model_id": model_id,
            "source_dataset": source_dataset,
            "seed": seed,
            "threshold_name": "tau_star",
            "selection_dataset": "source_validation_only",
            "objective": "maximize_MCC",
            "prediction_rule": "AMP_if_score_greater_than_or_equal_to_tau_star",
            "tie_breaking_order": [
                "higher_validation_MCC",
                "higher_validation_BACC",
                "higher_validation_F1",
                "closest_to_0.5",
                "lower_threshold",
            ],
            "source_test_used": False,
            "external_target_used": False,
            **threshold_info,
        }
        atomic_write_json(staging_dir / "validation_threshold.json", threshold_payload)

        validation_predictions = prediction_frame(
            run_id=run_id,
            protocol="validation_threshold_selection",
            model_id=model_id,
            source_dataset=source_dataset,
            target_dataset=source_dataset,
            seed=seed,
            sample_ids=validation_df["sample_id"],
            y_true=validation_df["label"],
            y_score=val_score,
            tau=tau,
        )
        assert list(validation_predictions.columns) == PREDICTION_COLUMNS
        assert len(validation_predictions) == int(run_row.validation_n)
        atomic_write_csv(staging_dir / "validation_predictions.csv", validation_predictions)

        checkpoint = Path(checkpoint)
        assert checkpoint.is_file()
        assert checkpoint.is_relative_to(staging_dir)

        cell32_run_wall_seconds = time.perf_counter() - run_wall_start
        resources_payload = {
            "run_id": run_id,
            "model_id": model_id,
            "source_dataset": source_dataset,
            "seed": seed,
            "device": (
                "method_specific"
                if model_id in PUBLISHED_IDS
                else ("gpu_tensorflow_locked_cell2_subprocess" if model_id in CONTROLLED_IDS[6:] else "cpu")
            ),
            "feature_or_token_load_seconds": feature_or_token_load_seconds,
            "preprocessing_seconds": preprocessing_seconds,
            "training_frame_load_seconds": float(frame_load_seconds),
            "training_seconds": float(training_seconds),
            "validation_seconds": float(validation_seconds),
            "source_test_inference_seconds": None,
            "source_test_ms_per_sequence": None,
            "source_test_sequences_per_second": None,
            "peak_process_ram_mb": float(process_peak_ram_mb()),
            "peak_gpu_vram_mb": gpu_peak_override if model_id in CONTROLLED_IDS[6:] else gpu_peak_mb_best_effort(),
            "model_parameter_count": model_parameter_count,
            "checkpoint_size_mb": float(checkpoint.stat().st_size / (1024 ** 2)),
            "cell32_run_wall_seconds": float(cell32_run_wall_seconds),
            "timing_provenance": timing_provenance,
            "executor_version": "v9_adaptive_throughput_timing_exact_route_recovery",
            "executor_lane": os.environ.get("GENPEPT_CELL32_EXECUTOR_LANE", "serial_or_parent"),
            "executor_mode": os.environ.get("GENPEPT_CELL32_EXECUTOR_MODE", "resource_aware_parallel"),
            "executor_parallel_context": os.environ.get("GENPEPT_CELL32_PARALLEL_CONTEXT", "unknown"),
            "executor_cpu_capacity_units": int(os.environ.get("GENPEPT_CELL32_CPU_CAPACITY_UNITS", "0") or 0),
        }
        atomic_write_json(staging_dir / "resources.json", resources_payload)

        # Explicit anti-leakage assertions for Cell 32.
        assert not (staging_dir / "source_test_predictions.csv").exists()
        assert not (staging_dir / "source_test_metrics.json").exists()
        assert not (staging_dir / "external").exists()

        training_files = [
            staging_dir / "config.json",
            staging_dir / "validation_predictions.csv",
            staging_dir / "validation_threshold.json",
            staging_dir / "history.csv",
            staging_dir / "resources.json",
        ]
        for p in training_files:
            assert p.is_file(), f"Missing Cell-32 training artifact: {p}"

        manifest = {
            "run_id": run_id,
            "model_id": model_id,
            "source_dataset": source_dataset,
            "seed": seed,
            "benchmark_version": "benchmark4_FINAL_v1",
            "cell32_trainer_version": "v9_adaptive_throughput_timing_exact_route_recovery",
            "cell32_trainer_sha256": CELL32_TRAINER_SHA256 or None,
            "technical_fix": "v4 mask dtype correction retained; v9 retains exact scientific training and adds adaptive throughput scheduling/timing/recovery and replaces AMPScannerV2 helper pandas-only CSV output with stdlib csv; scientific contract unchanged",
            "controlled_contract_sha256": sha256_file(CONTROLLED_CONTRACT_PATH),
            "evaluation_contract_sha256": sha256_file(EVALUATION_CONTRACT_PATH),
            "feature_cache_manifest_sha256": sha256_file(FEATURE_CACHE_MANIFEST_PATH),
            "published_route_contract_sha256": sha256_file(PUBLISHED_ROUTE_CONTRACT_PATH),
            "train_index_key": str(run_row.train_index_key),
            "validation_index_key": str(run_row.validation_index_key),
            "source_test_index_key": str(run_row.source_test_index_key),
            "validation_threshold": tau,
            "status": "TRAINED",
            "trained_at": now_iso(),
            "cell32_source_test_read": False,
            "cell32_external_target_read": False,
            "checkpoint_relative_path": checkpoint.relative_to(staging_dir).as_posix(),
            "checkpoint_artifact": training_artifact_entry(checkpoint),
            "training_artifacts": {
                p.relative_to(staging_dir).as_posix(): training_artifact_entry(p)
                for p in training_files
            },
        }
        if model_id in PUBLISHED_IDS:
            manifest["selftest_manifest"] = str(run_row.selftest_manifest)
            manifest["selftest_sha256"] = str(run_row.selftest_sha256)

        atomic_write_json(staging_dir / "run_manifest.json", manifest)

        # Validate bundle before atomic directory promotion.
        manifest_check = load_json(staging_dir / "run_manifest.json")
        assert manifest_check["status"] == "TRAINED"
        required_fields = artifact_schema_wrapper["schema"]["run_manifest_required_fields"] \
            if "schema" in artifact_schema_wrapper else artifact_schema_wrapper["contract"]["run_manifest_required_fields"] \
            if "contract" in artifact_schema_wrapper else artifact_schema_wrapper["run_manifest_required_fields"]
        for field in required_fields:
            assert field in manifest_check, f"run_manifest missing required field: {field}"
        assert manifest_check["cell32_source_test_read"] is False
        assert manifest_check["cell32_external_target_read"] is False

        # Atomic commit of the complete TRAINED bundle.
        os.replace(staging_dir, final_run_dir)
        staging_dir = None

        # Final disk verification after promotion.
        assert trained_manifest_valid(run_row)
        return final_run_dir / "run_manifest.json"

    finally:
        if staging_dir is not None and Path(staging_dir).exists():
            # Failed staging is preserved under a diagnostic name rather than
            # being silently deleted or promoted as a valid run.
            failed_dir = Path(str(staging_dir) + ".FAILED")
            if not failed_dir.exists():
                os.replace(staging_dir, failed_dir)


# ============================================================
# LOCKED CELL-2 GPU PREFLIGHT / SINGLE-RUN WORKER MODE
# ============================================================
WORKER_RUN_ID = os.environ.get("GENPEPT_CELL32_WORKER_RUN_ID", "").strip()
WORKER_ROW = None
if WORKER_RUN_ID:
    _worker_matches = [
        row for row in full_plan.itertuples(index=False)
        if str(row.run_id) == WORKER_RUN_ID
    ]
    assert len(_worker_matches) == 1, f"Unknown/duplicate worker run_id: {WORKER_RUN_ID}"
    WORKER_ROW = _worker_matches[0]

if WORKER_ROW is None:
    TF_GPU_PREFLIGHT = verify_tensorflow_gpu_subprocess()
    print("TensorFlow GPU preflight: PASS -", TF_GPU_PREFLIGHT)
elif str(WORKER_ROW.model_id) in CONTROLLED_IDS[6:]:
    TF_GPU_PREFLIGHT = verify_tensorflow_gpu_subprocess()
    print("Worker TensorFlow GPU preflight: PASS -", TF_GPU_PREFLIGHT)
else:
    TF_GPU_PREFLIGHT = {"skipped": True, "reason": "CPU worker does not use the locked master TensorFlow GPU route"}

# A worker executes exactly one already-selected run and never mutates the shared ledger.
# Ledger mutation is parent-only to avoid concurrent CSV races.
if WORKER_ROW is not None:
    assert not trained_manifest_valid(WORKER_ROW), (
        "Worker was asked to retrain an already valid canonical run: " + WORKER_RUN_ID
    )
    _worker_manifest = execute_training_run(WORKER_ROW)
    assert trained_manifest_valid(WORKER_ROW)
    print(f"CELL32_V9_WORKER_SUCCESS\t{WORKER_RUN_ID}\t{_worker_manifest}", flush=True)
    raise SystemExit(0)


# ============================================================
# PARENT: RECONCILE LEDGER AGAINST VALID TRAINED ARTIFACTS
# ============================================================
ledger = read_csv_exact(RUN_LEDGER_PATH)
assert len(ledger) == 680
assert ledger["run_id"].is_unique

# Recover stale RUNNING rows left by earlier notebook interrupts. A stale row is
# retryable only when no canonical validated training manifest exists.
stale_running = []
for row in full_plan.itertuples(index=False):
    run_id = str(row.run_id)
    mask = ledger["run_id"].astype(str).eq(run_id)
    status = str(ledger.loc[mask, "status"].iloc[0])
    if status == "RUNNING" and not trained_manifest_valid(row):
        stale_running.append(run_id)

if stale_running:
    evidence = {
        "schema_version": "1.0",
        "created_at": now_iso(),
        "status": "STALE_RUNNING_RECOVERED_FOR_RETRY",
        "run_ids": stale_running,
        "reason": "Prior notebook execution was interrupted before canonical run promotion.",
        "scientific_change": False,
        "source_test_read": False,
        "external_target_read": False,
        "recovery_executor": "v9_adaptive_throughput_timing_exact_route_recovery",
    }
    stale_key = sha256_text("\n".join(sorted(stale_running)))[:12]
    interrupt_audit_path = AUDIT_DIR / f"benchmark17_CELL32_v9_interrupt_recovery_{stale_key}.json"
    if interrupt_audit_path.is_file():
        prior = load_json(interrupt_audit_path)
        assert sorted(prior["run_ids"]) == sorted(stale_running)
    else:
        atomic_write_json(interrupt_audit_path, evidence)
    for run_id in stale_running:
        update_ledger(
            run_id,
            status="FAILED",
            finished_at=now_iso(),
            run_manifest="",
            last_error="Recovered stale RUNNING; no valid canonical run manifest; scheduled for exact retry by v9.",
        )
    ledger = read_csv_exact(RUN_LEDGER_PATH)
    print("Stale RUNNING recovery   :", stale_running)
else:
    print("Stale RUNNING recovery   : NONE")

# Rebind any valid on-disk run that the ledger does not currently acknowledge.
for row in full_plan.itertuples(index=False):
    if trained_manifest_valid(row):
        mask = ledger["run_id"].astype(str).eq(str(row.run_id))
        current_status = str(ledger.loc[mask, "status"].iloc[0])
        if current_status not in {"TRAINED", "EVALUATED", "COMPLETE"}:
            update_ledger(
                str(row.run_id),
                status="TRAINED",
                finished_at=now_iso(),
                run_manifest=str(Path(str(row.run_dir)) / "run_manifest.json"),
                last_error="",
            )
            ledger = read_csv_exact(RUN_LEDGER_PATH)

ledger = read_csv_exact(RUN_LEDGER_PATH)
status_series = ledger["status"].astype(str)
allowed_status = {"PLANNED", "RUNNING", "TRAINED", "EVALUATED", "COMPLETE", "FAILED"}
assert set(status_series.unique()).issubset(allowed_status)

trained_like = {
    str(row.run_id)
    for row in full_plan.itertuples(index=False)
    if trained_manifest_valid(row)
}
claimed_trained = set(
    ledger.loc[
        status_series.isin(["TRAINED", "EVALUATED", "COMPLETE"]),
        "run_id",
    ].astype(str)
)
assert claimed_trained == trained_like, (
    "Ledger/artifact mismatch before v9 launch.\n"
    f"claimed_only={sorted(claimed_trained - trained_like)}\n"
    f"manifest_only={sorted(trained_like - claimed_trained)}"
)

pending_rows = [
    row for row in full_plan.itertuples(index=False)
    if str(row.run_id) not in trained_like
]

# Optional operational cap only; 0 means all remaining jobs.
MAX_JOBS = int(os.environ.get("GENPEPT_CELL32_MAX_JOBS", "0"))
assert MAX_JOBS >= 0
if MAX_JOBS > 0:
    pending_rows = pending_rows[:MAX_JOBS]

# ------------------------------------------------------------
# ADAPTIVE THROUGHPUT EXECUTOR POLICY (v9)
# ------------------------------------------------------------
# IMPORTANT:
# - Scientific run identity/order, data, model source, seed, hyperparameters,
#   threshold policy and artifact validation are unchanged.
# - v9 changes ONLY operational scheduling: independent workers may overlap.
# - Controlled deep workers keep TF memory-growth enabled and deterministic ops.
# - Published legacy helpers keep their locked environments/backends unchanged.
# - Timing is logged with explicit parallel-contention provenance and must not be
#   interpreted as uncontended cross-model speed unless separately measured.

def _mem_available_gb():
    try:
        with open('/proc/meminfo', 'r', encoding='utf-8') as f:
            kv = {}
            for line in f:
                key, val = line.split(':', 1)
                kv[key] = val.strip()
        kb = float(kv.get('MemAvailable', '0 kB').split()[0])
        return kb / (1024.0 ** 2)
    except Exception:
        return 0.0


def _gpu_inventory():
    info = {'name': 'UNKNOWN', 'total_mb': 0.0, 'query_ok': False}
    try:
        r = subprocess.run(
            ['nvidia-smi', '--query-gpu=name,memory.total', '--format=csv,noheader,nounits'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False,
        )
        if r.returncode == 0 and r.stdout.strip():
            first = r.stdout.strip().splitlines()[0]
            name, mem = first.rsplit(',', 1)
            info = {'name': name.strip(), 'total_mb': float(mem.strip()), 'query_ok': True}
    except Exception:
        pass
    return info


LOGICAL_CPUS = int(os.cpu_count() or 4)
AVAILABLE_RAM_GB = float(_mem_available_gb())
GPU_INFO = _gpu_inventory()

# Existing canonical timing/resource data are used only for operational scheduling
# and ETA. They do not influence model fitting, thresholds, predictions or metrics.
existing_resources = existing_run_resource_rows(full_plan)
_historical_model_seconds = {}
if len(existing_resources):
    _tmp = existing_resources.copy()
    if 'cell32_run_wall_seconds' in _tmp.columns:
        _tmp['_wall'] = pd.to_numeric(_tmp['cell32_run_wall_seconds'], errors='coerce')
    else:
        _tmp['_wall'] = np.nan
    _tmp['_train'] = pd.to_numeric(_tmp.get('training_seconds'), errors='coerce')
    _tmp['_val'] = pd.to_numeric(_tmp.get('validation_seconds'), errors='coerce')
    _tmp['_fallback'] = _tmp['_train'].fillna(0) + _tmp['_val'].fillna(0)
    _tmp['_est'] = _tmp['_wall'].where(_tmp['_wall'].notna() & (_tmp['_wall'] > 0), _tmp['_fallback'])
    for mid, g in _tmp.groupby('model_id'):
        vals = pd.to_numeric(g['_est'], errors='coerce')
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if len(vals):
            _historical_model_seconds[str(mid)] = float(vals.median())

_gpu_peak_reference_mb = 0.0
if len(existing_resources) and 'peak_gpu_vram_mb' in existing_resources.columns:
    _g = existing_resources.loc[
        existing_resources['model_id'].astype(str).isin(CONTROLLED_IDS[6:]),
        'peak_gpu_vram_mb',
    ]
    _g = pd.to_numeric(_g, errors='coerce')
    _g = _g[np.isfinite(_g) & (_g > 0)]
    if len(_g):
        _gpu_peak_reference_mb = float(_g.max())

# GPU concurrency is evidence-gated. Two workers are enabled only when prior
# canonical deep runs show enough VRAM headroom on the detected GPU. Otherwise 1.
AUTO_GPU_MAX_WORKERS = 1
if (
    GPU_INFO['query_ok']
    and GPU_INFO['total_mb'] >= 10000.0
    and _gpu_peak_reference_mb > 0.0
    and (2.0 * _gpu_peak_reference_mb) <= (0.72 * GPU_INFO['total_mb'])
    and AVAILABLE_RAM_GB >= 8.0
    and LOGICAL_CPUS >= 8
):
    AUTO_GPU_MAX_WORKERS = 2
GPU_MAX_WORKERS = int(os.environ.get('GENPEPT_CELL32_GPU_MAX_WORKERS', str(AUTO_GPU_MAX_WORKERS)))
GPU_MAX_WORKERS = max(1, min(2, GPU_MAX_WORKERS))

# On the observed 8-logical-CPU / ~10-11 GiB-available class of machine, two
# independent single-thread/legacy CPU workers are allowed. AmpGram stays
# exclusive below 16 GiB available because its 33,620-feature construction and
# two ranger forests are the materially highest-memory CPU route.
AUTO_CPU_MAX_WORKERS = 2 if (LOGICAL_CPUS >= 8 and AVAILABLE_RAM_GB >= 8.0) else 1
CPU_MAX_WORKERS = int(os.environ.get('GENPEPT_CELL32_CPU_MAX_WORKERS', str(AUTO_CPU_MAX_WORKERS)))
CPU_MAX_WORKERS = max(1, min(3, CPU_MAX_WORKERS))
CPU_EXCLUSIVE_MODELS = {'ampgram'} if AVAILABLE_RAM_GB < 16.0 else set()

# Bound total concurrent worker processes to avoid turning parallelism into swap/OOM.
if AVAILABLE_RAM_GB < 8.0:
    AUTO_TOTAL_ACTIVE_LIMIT = 2
elif AVAILABLE_RAM_GB < 12.0:
    AUTO_TOTAL_ACTIVE_LIMIT = 3
elif AVAILABLE_RAM_GB < 20.0:
    AUTO_TOTAL_ACTIVE_LIMIT = 4
else:
    AUTO_TOTAL_ACTIVE_LIMIT = 5
TOTAL_ACTIVE_LIMIT = int(os.environ.get('GENPEPT_CELL32_TOTAL_ACTIVE_LIMIT', str(AUTO_TOTAL_ACTIVE_LIMIT)))
TOTAL_ACTIVE_LIMIT = max(1, min(6, TOTAL_ACTIVE_LIMIT))
MIN_MEM_AVAILABLE_GB = float(os.environ.get('GENPEPT_CELL32_MIN_MEM_AVAILABLE_GB', '2.5'))
POLL_SECONDS = float(os.environ.get('GENPEPT_CELL32_POLL_SECONDS', '0.20'))
assert POLL_SECONDS > 0

# Compatibility/logging weight only. Scheduling itself uses explicit worker-count,
# exclusivity, total-active and live-memory guards instead of the flawed v8
# overweight fallback.
CPU_WEIGHTS = {
    'logistic_regression': 1,
    'linear_svm_calibrated': 1,
    'random_forest': 1,
    'extra_trees': 1,
    'hist_gradient_boosting': 1,
    'lightgbm': 1,
    'ampep': 1,
    'ampir': 1,
    'ampgram': 3,
    'ampscannerv2': 1,
    'ampeppy': 1,
    'ai4amp': 2,
    'amplify': 2,
}


def _lane(row):
    return 'GPU' if str(row.model_id) in CONTROLLED_IDS[6:] else 'CPU'


def _weight(row):
    return int(CPU_WEIGHTS.get(str(row.model_id), 1))


pending_gpu = [row for row in pending_rows if _lane(row) == 'GPU']
pending_cpu = [row for row in pending_rows if _lane(row) == 'CPU']
TOTAL_THIS_INVOCATION = len(pending_rows)

EVENT_COLUMNS = [
    'event_no', 'invocation_id', 'run_id', 'model_id', 'source_dataset', 'seed',
    'lane', 'cpu_weight_units', 'started_at', 'finished_at', 'outer_wall_seconds',
    'returncode', 'status', 'worker_log', 'training_seconds', 'validation_seconds',
    'feature_or_token_load_seconds', 'preprocessing_seconds', 'cell32_run_wall_seconds',
    'parallel_context', 'cpu_capacity_units', 'logical_cpus', 'available_ram_gb',
    'gpu_max_workers', 'cpu_max_workers', 'total_active_limit', 'gpu_name',
    'gpu_total_mb', 'gpu_peak_reference_mb', 'resource_retry_no',
]
if CELL32_V9_EVENTS_PATH.is_file():
    events_df = read_csv_exact(CELL32_V9_EVENTS_PATH)
    for col in EVENT_COLUMNS:
        if col not in events_df.columns:
            events_df[col] = ''
    events_df = events_df[EVENT_COLUMNS]
else:
    events_df = pd.DataFrame(columns=EVENT_COLUMNS)

INVOCATION_ID = datetime.now().astimezone().strftime('%Y%m%dT%H%M%S%z')
WORKER_LOG_DIR = WORK_ROOT / 'worker_logs' / 'benchmark4_FINAL_v1' / 'cell32_v9'
WORKER_LOG_DIR.mkdir(parents=True, exist_ok=True)

print('=' * 128)
print('CELL 32 v9 - ADAPTIVE FAST TRAINING + TIMING + EXACT ROUTING + SAFE RECOVERY')
print('=' * 128)
print('Locked training jobs     : 680')
print('Already trained          :', len(trained_like))
print('Pending this invocation  :', TOTAL_THIS_INVOCATION)
print('GPU-lane pending         :', len(pending_gpu), f'(max concurrent={GPU_MAX_WORKERS})')
print('CPU-lane pending         :', len(pending_cpu), f'(max concurrent={CPU_MAX_WORKERS})')
print('Total active worker cap  :', TOTAL_ACTIVE_LIMIT)
print('Logical CPUs             :', LOGICAL_CPUS)
print('Available RAM (GiB)      :', f'{AVAILABLE_RAM_GB:.1f}')
print('Live RAM launch floor    :', f'{MIN_MEM_AVAILABLE_GB:.1f} GiB')
print('GPU                      :', GPU_INFO['name'])
print('GPU total VRAM (MiB)     :', f"{GPU_INFO['total_mb']:.0f}" if GPU_INFO['query_ok'] else 'UNKNOWN')
print('Observed deep peak (MiB) :', f'{_gpu_peak_reference_mb:.1f}' if _gpu_peak_reference_mb > 0 else 'UNAVAILABLE')
print('AmpGram CPU exclusivity  :', 'TRUE' if 'ampgram' in CPU_EXCLUSIVE_MODELS else 'FALSE')
print('Execution mode           : adaptive GPU concurrency + bounded CPU parallel lane')
print('Worker trainer path      :', SELF_TRAINER_PATH)
print('Worker trainer SHA256    :', _verify_self_trainer_path())
print('AMPScannerV2 helper      :', HELPER_ROOT / 'ampscannerv2_production_v2.py')
print('Timing context           : parallel resource-aware; contention provenance logged')
print('Source TEST rows read    : FALSE')
print('External target rows read: FALSE')
print('Scientific protocol      : UNCHANGED')
print('=' * 128)

batch_wall_start = time.perf_counter()
completed_this_invocation = 0
active = {}
failure_records = []
observed_by_model = {}
observed_by_lane = {'GPU': [], 'CPU': []}
stop_launching = False
GPU_MAX_WORKERS_CURRENT = GPU_MAX_WORKERS
CPU_MAX_WORKERS_CURRENT = CPU_MAX_WORKERS
TOTAL_ACTIVE_LIMIT_CURRENT = TOTAL_ACTIVE_LIMIT
resource_retry_counts = {}
MAX_RESOURCE_RETRIES = 1


def _estimate_seconds(row):
    mid = str(row.model_id)
    vals = observed_by_model.get(mid, [])
    if vals:
        return float(np.median(vals))
    if mid in _historical_model_seconds:
        return float(_historical_model_seconds[mid])
    lane_vals = observed_by_lane.get(_lane(row), [])
    if lane_vals:
        return float(np.median(lane_vals))
    return 60.0


def _eta_parallel_seconds():
    rem_gpu = list(pending_gpu) + [r['row'] for r in active.values() if r['lane'] == 'GPU']
    rem_cpu = list(pending_cpu) + [r['row'] for r in active.values() if r['lane'] == 'CPU']
    gpu_sec = sum(_estimate_seconds(r) for r in rem_gpu) / max(1, GPU_MAX_WORKERS_CURRENT)
    exclusive = [r for r in rem_cpu if str(r.model_id) in CPU_EXCLUSIVE_MODELS]
    normal = [r for r in rem_cpu if str(r.model_id) not in CPU_EXCLUSIVE_MODELS]
    # Exclusive work is serialized within CPU lane; normal work is divided by CPU workers.
    cpu_sec = sum(_estimate_seconds(r) for r in exclusive)
    cpu_sec += sum(_estimate_seconds(r) for r in normal) / max(1, CPU_MAX_WORKERS_CURRENT)
    return max(gpu_sec, cpu_sec)


def _append_event(rec, status, returncode, outer_wall, rsrc=None, retry_no=0):
    global events_df
    row = rec['row']
    rsrc = rsrc or {}
    event = {
        'event_no': int(len(events_df) + 1),
        'invocation_id': INVOCATION_ID,
        'run_id': str(row.run_id),
        'model_id': str(row.model_id),
        'source_dataset': str(row.source_dataset),
        'seed': int(row.seed),
        'lane': rec['lane'],
        'cpu_weight_units': int(rec['weight']),
        'started_at': rec['started_at'],
        'finished_at': now_iso(),
        'outer_wall_seconds': float(outer_wall),
        'returncode': int(returncode),
        'status': status,
        'worker_log': str(rec['log_path']),
        'training_seconds': rsrc.get('training_seconds'),
        'validation_seconds': rsrc.get('validation_seconds'),
        'feature_or_token_load_seconds': rsrc.get('feature_or_token_load_seconds'),
        'preprocessing_seconds': rsrc.get('preprocessing_seconds'),
        'cell32_run_wall_seconds': rsrc.get('cell32_run_wall_seconds'),
        'parallel_context': 'adaptive_resource_aware_parallel',
        'cpu_capacity_units': CPU_MAX_WORKERS_CURRENT,
        'logical_cpus': LOGICAL_CPUS,
        'available_ram_gb': float(_mem_available_gb()),
        'gpu_max_workers': GPU_MAX_WORKERS_CURRENT,
        'cpu_max_workers': CPU_MAX_WORKERS_CURRENT,
        'total_active_limit': TOTAL_ACTIVE_LIMIT_CURRENT,
        'gpu_name': GPU_INFO['name'],
        'gpu_total_mb': GPU_INFO['total_mb'],
        'gpu_peak_reference_mb': _gpu_peak_reference_mb,
        'resource_retry_no': int(retry_no),
    }
    events_df = pd.concat([events_df, pd.DataFrame([event])], ignore_index=True)
    atomic_write_csv(CELL32_V9_EVENTS_PATH, events_df[EVENT_COLUMNS])


def _tail_text(path, max_chars=12000):
    try:
        text = Path(path).read_text(encoding='utf-8', errors='replace')
        return text[-max_chars:]
    except Exception as exc:
        return f'<unable to read worker log: {exc}>'


def _gpu_active_count():
    return sum(1 for rec in active.values() if rec['lane'] == 'GPU')


def _cpu_active_records():
    return [rec for rec in active.values() if rec['lane'] == 'CPU']


def _cpu_active_count():
    return len(_cpu_active_records())


def _exclusive_cpu_active():
    return any(str(rec['row'].model_id) in CPU_EXCLUSIVE_MODELS for rec in _cpu_active_records())


def _live_memory_allows_launch():
    current = float(_mem_available_gb())
    return (current <= 0.0) or (current >= MIN_MEM_AVAILABLE_GB)


def _select_cpu_index():
    if not pending_cpu:
        return None
    if len(active) >= TOTAL_ACTIVE_LIMIT_CURRENT:
        return None
    if _cpu_active_count() >= CPU_MAX_WORKERS_CURRENT:
        return None
    if not _live_memory_allows_launch():
        return None
    if _exclusive_cpu_active():
        return None
    cpu_active = _cpu_active_count()
    for idx, row in enumerate(pending_cpu):
        is_exclusive = str(row.model_id) in CPU_EXCLUSIVE_MODELS
        if is_exclusive and cpu_active > 0:
            continue
        return idx
    return None


def _launch(row, lane):
    run_id = str(row.run_id)
    current = read_csv_exact(RUN_LEDGER_PATH)
    mask = current['run_id'].astype(str).eq(run_id)
    assert int(mask.sum()) == 1
    attempts = int(pd.to_numeric(current.loc[mask, 'attempt_count'], errors='raise').iloc[0]) + 1
    update_ledger(
        run_id, status='RUNNING', attempt_count=attempts, started_at=now_iso(),
        finished_at='', run_manifest='', last_error='',
    )
    weight = 0 if lane == 'GPU' else _weight(row)
    log_path = WORKER_LOG_DIR / f'{run_id}__attempt_{attempts}.log'
    log_fh = open(log_path, 'w', encoding='utf-8')
    env = os.environ.copy()
    env['GENPEPT_CELL32_WORKER_RUN_ID'] = run_id
    env['GENPEPT_CELL32_EXECUTOR_LANE'] = lane
    env['GENPEPT_CELL32_EXECUTOR_MODE'] = 'adaptive_resource_aware_parallel'
    env['GENPEPT_CELL32_PARALLEL_CONTEXT'] = 'parallel_resource_contended_not_for_unqualified_cross_model_runtime_claims'
    env['GENPEPT_CELL32_CPU_CAPACITY_UNITS'] = str(CPU_MAX_WORKERS_CURRENT)
    worker_trainer_sha = _verify_self_trainer_path()
    env['GENPEPT_CELL32_TRAINER_PATH'] = str(SELF_TRAINER_PATH)
    env['GENPEPT_CELL32_TRAINER_SHA256'] = worker_trainer_sha
    cmd = [str(Path(os.sys.executable).resolve()), str(SELF_TRAINER_PATH)]
    proc = subprocess.Popen(
        cmd, stdout=log_fh, stderr=subprocess.STDOUT, text=True, env=env,
        start_new_session=True,
    )
    rec = {
        'proc': proc, 'row': row, 'lane': lane, 'weight': weight,
        'log_path': log_path, 'log_fh': log_fh, 'start_perf': time.perf_counter(),
        'started_at': now_iso(),
    }
    active[proc.pid] = rec
    print(
        f'[LAUNCH {completed_this_invocation + len(active)}/{TOTAL_THIS_INVOCATION}] '
        f'{lane:<3} {run_id}' + (f' | weight={weight}' if lane == 'CPU' else '') +
        f' | active_gpu={_gpu_active_count()} active_cpu={_cpu_active_count()} '
        f'ram_now={_mem_available_gb():.1f}GiB',
        flush=True,
    )


def _resource_failure(lane, rc, log_tail):
    text = (log_tail or '').lower()
    gpu_markers = [
        'resourceexhaustederror', 'cuda_error_out_of_memory', 'out of memory',
        'failed to allocate memory', 'cudnn_status_alloc_failed',
    ]
    if lane == 'GPU' and any(x in text for x in gpu_markers):
        return True
    # SIGKILL/137 under parallel execution is treated as operational resource
    # contention only once; any repeat becomes a hard failure for audit.
    if int(rc) in {-9, 137}:
        return True
    return False


def _terminate_active(reason):
    for rec in list(active.values()):
        proc = rec['proc']
        if proc.poll() is None:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except Exception:
                try:
                    proc.terminate()
                except Exception:
                    pass
    deadline = time.time() + 10.0
    for rec in list(active.values()):
        proc = rec['proc']
        while proc.poll() is None and time.time() < deadline:
            time.sleep(0.1)
        if proc.poll() is None:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except Exception:
                try:
                    proc.kill()
                except Exception:
                    pass
        try:
            rec['log_fh'].close()
        except Exception:
            pass
        row = rec['row']
        if trained_manifest_valid(row):
            update_ledger(
                str(row.run_id), status='TRAINED', finished_at=now_iso(),
                run_manifest=str(Path(str(row.run_dir)) / 'run_manifest.json'), last_error='',
            )
        else:
            update_ledger(
                str(row.run_id), status='FAILED', finished_at=now_iso(), run_manifest='',
                last_error=reason[:2000],
            )


try:
    while pending_gpu or pending_cpu or active:
        # GPU: launch up to the evidence-gated concurrency cap, subject to total
        # active-worker and live-RAM guards.
        if not stop_launching:
            while (
                pending_gpu
                and _gpu_active_count() < GPU_MAX_WORKERS_CURRENT
                and len(active) < TOTAL_ACTIVE_LIMIT_CURRENT
                and _live_memory_allows_launch()
            ):
                _launch(pending_gpu.pop(0), 'GPU')

        # CPU: bounded parallel workers. AmpGram is exclusive on low-RAM hosts.
        if not stop_launching:
            while pending_cpu:
                idx = _select_cpu_index()
                if idx is None:
                    break
                _launch(pending_cpu.pop(idx), 'CPU')

        finished_pids = []
        for pid, rec in list(active.items()):
            rc = rec['proc'].poll()
            if rc is None:
                continue
            finished_pids.append(pid)
            try:
                rec['log_fh'].flush()
                rec['log_fh'].close()
            except Exception:
                pass
            row = rec['row']
            run_id = str(row.run_id)
            outer_wall = time.perf_counter() - rec['start_perf']
            valid = (rc == 0) and trained_manifest_valid(row)
            if valid:
                manifest_path = Path(str(row.run_dir)) / 'run_manifest.json'
                update_ledger(
                    run_id, status='TRAINED', finished_at=now_iso(),
                    run_manifest=str(manifest_path), last_error='',
                )
                rsrc = load_json(Path(str(row.run_dir)) / 'resources.json')
                _append_event(rec, 'TRAINED', rc, outer_wall, rsrc, retry_no=resource_retry_counts.get(run_id, 0))
                completed_this_invocation += 1
                observed_by_model.setdefault(str(row.model_id), []).append(float(outer_wall))
                observed_by_lane[rec['lane']].append(float(outer_wall))
                elapsed = time.perf_counter() - batch_wall_start
                eta = _eta_parallel_seconds()
                print(
                    f'[{completed_this_invocation}/{TOTAL_THIS_INVOCATION}] TRAINED {run_id}\n'
                    '          '
                    f"lane={rec['lane']} | train={float(rsrc.get('training_seconds', float('nan'))):.2f}s | "
                    f"val={float(rsrc.get('validation_seconds', float('nan'))):.2f}s | "
                    f'worker_wall={outer_wall:.2f}s | elapsed={format_duration(elapsed)} | '
                    f'ETA_parallel_est={format_duration(eta)} | active_after={len(active)-1}',
                    flush=True,
                )
            else:
                log_tail = _tail_text(rec['log_path'])
                retry_no = int(resource_retry_counts.get(run_id, 0))
                if _resource_failure(rec['lane'], rc, log_tail) and retry_no < MAX_RESOURCE_RETRIES:
                    retry_no += 1
                    resource_retry_counts[run_id] = retry_no
                    _append_event(rec, 'RESOURCE_RETRY_SCHEDULED', rc, outer_wall, {}, retry_no=retry_no)
                    update_ledger(
                        run_id, status='PLANNED', finished_at=now_iso(), run_manifest='',
                        last_error=(f'v9 resource-contention retry scheduled after returncode={rc}; retry={retry_no}')[:2000],
                    )
                    if rec['lane'] == 'GPU':
                        GPU_MAX_WORKERS_CURRENT = 1
                        pending_gpu.insert(0, row)
                    else:
                        CPU_MAX_WORKERS_CURRENT = 1
                        pending_cpu.insert(0, row)
                    TOTAL_ACTIVE_LIMIT_CURRENT = max(1, min(TOTAL_ACTIVE_LIMIT_CURRENT, GPU_MAX_WORKERS_CURRENT + CPU_MAX_WORKERS_CURRENT))
                    print('-' * 128)
                    print(f'CELL 32 v9 RESOURCE BACKOFF — {run_id}')
                    print(f'returncode             : {rc}')
                    print(f'retry                  : {retry_no}/{MAX_RESOURCE_RETRIES}')
                    print(f'GPU max workers now    : {GPU_MAX_WORKERS_CURRENT}')
                    print(f'CPU max workers now    : {CPU_MAX_WORKERS_CURRENT}')
                    print(f'Total active cap now   : {TOTAL_ACTIVE_LIMIT_CURRENT}')
                    print('Scientific protocol    : UNCHANGED; operational retry only')
                    print('-' * 128)
                else:
                    update_ledger(
                        run_id, status='FAILED', finished_at=now_iso(), run_manifest='',
                        last_error=(f'v9 worker returncode={rc}; canonical TRAINED manifest valid={trained_manifest_valid(row)}')[:2000],
                    )
                    _append_event(rec, 'FAILED', rc, outer_wall, {}, retry_no=retry_no)
                    failure_records.append((run_id, rc, str(rec['log_path']), log_tail))
                    stop_launching = True
                    print('-' * 128)
                    print(f'CELL 32 v9 WORKER FAILURE — {run_id}')
                    print(f'returncode             : {rc}')
                    print(f'worker log             : {rec["log_path"]}')
                    print(log_tail)
                    print('New launches stopped    : TRUE')
                    print('Already-active workers  : allowed to finish and commit if valid')
                    print('-' * 128)

        for pid in finished_pids:
            active.pop(pid, None)

        if stop_launching and not active:
            break
        if active:
            time.sleep(POLL_SECONDS)
        elif (pending_gpu or pending_cpu) and not stop_launching:
            # If only a live-memory guard blocks new work, wait rather than spin.
            time.sleep(max(POLL_SECONDS, 0.5))

except KeyboardInterrupt:
    print('\nCELL 32 v9 INTERRUPT RECEIVED - terminating active worker process groups safely...', flush=True)
    _terminate_active('KeyboardInterrupt in v9 parent; active worker terminated before canonical completion; exact retry required.')
    write_resource_checkpoint(full_plan, complete=False)
    summary_interrupt = {
        'schema_version': '1.0',
        'created_at': now_iso(),
        'status': 'INTERRUPTED',
        'invocation_id': INVOCATION_ID,
        'completed_this_invocation': completed_this_invocation,
        'elapsed_seconds': time.perf_counter() - batch_wall_start,
        'gpu_max_workers_initial': GPU_MAX_WORKERS,
        'gpu_max_workers_final': GPU_MAX_WORKERS_CURRENT,
        'cpu_max_workers_initial': CPU_MAX_WORKERS,
        'cpu_max_workers_final': CPU_MAX_WORKERS_CURRENT,
        'total_active_limit_initial': TOTAL_ACTIVE_LIMIT,
        'total_active_limit_final': TOTAL_ACTIVE_LIMIT_CURRENT,
        'logical_cpus': LOGICAL_CPUS,
        'available_ram_gb_start': AVAILABLE_RAM_GB,
        'gpu_info': GPU_INFO,
        'gpu_peak_reference_mb': _gpu_peak_reference_mb,
        'source_test_read': False,
        'external_target_read': False,
    }
    atomic_write_json(CELL32_V9_SUMMARY_PATH, summary_interrupt)
    raise

# Final checkpoint and summary.
final_ledger = read_csv_exact(RUN_LEDGER_PATH)
final_counts = final_ledger["status"].astype(str).value_counts().to_dict()
trained_count = int(final_ledger["status"].astype(str).isin(["TRAINED", "EVALUATED", "COMPLETE"]).sum())
failed_count = int((final_ledger["status"].astype(str) == "FAILED").sum())
write_resource_checkpoint(full_plan, complete=(trained_count == 680 and failed_count == 0))

executor_summary = {
    "schema_version": "1.0",
    "created_at": now_iso(),
    "invocation_id": INVOCATION_ID,
    "executor_version": "v9_adaptive_throughput_timing_exact_route_recovery",
    "scientific_protocol_changed": False,
    "execution_mode": "adaptive_gpu_concurrency_plus_bounded_cpu_parallel_lane",
    "cpu_capacity_units": CPU_MAX_WORKERS_CURRENT,
    "gpu_cpu_reserve_units": 0,
    "gpu_max_workers_initial": GPU_MAX_WORKERS,
    "gpu_max_workers_final": GPU_MAX_WORKERS_CURRENT,
    "cpu_max_workers_initial": CPU_MAX_WORKERS,
    "cpu_max_workers_final": CPU_MAX_WORKERS_CURRENT,
    "total_active_limit_initial": TOTAL_ACTIVE_LIMIT,
    "total_active_limit_final": TOTAL_ACTIVE_LIMIT_CURRENT,
    "gpu_info": GPU_INFO,
    "gpu_peak_reference_mb": _gpu_peak_reference_mb,
    "logical_cpus": LOGICAL_CPUS,
    "available_ram_gb": AVAILABLE_RAM_GB,
    "completed_this_invocation": completed_this_invocation,
    "invocation_elapsed_seconds": float(time.perf_counter() - batch_wall_start),
    "trained_total": trained_count,
    "failed_total": failed_count,
    "source_test_read": False,
    "external_target_read": False,
    "event_log": str(CELL32_V9_EVENTS_PATH),
    "resource_log": str(CELL32_RESOURCE_LOG_PATH),
    "resource_summary": str(CELL32_RESOURCE_SUMMARY_PATH),
    "timing_interpretation": "Per-run times are measured under the recorded resource-aware parallel context; do not make unqualified cross-model runtime claims without accounting for contention.",
}
atomic_write_json(CELL32_V9_SUMMARY_PATH, executor_summary)

print("=" * 128)
if failure_records:
    print("CELL 32 v9 FAILURE CHECKPOINT - AT LEAST ONE WORKER FAILED")
elif trained_count == 680 and failed_count == 0:
    print("CELL 32 PASS - ALL 680 LOCKED SOURCE TRAININGS COMPLETED")
else:
    print("CELL 32 CHECKPOINT - TRAINING NOT YET COMPLETE")
print("=" * 128)
print("Ledger status            :", final_counts)
print("TRAINED/EVALUATED/COMPLETE:", trained_count)
print("FAILED                   :", failed_count)
print("Completed this invocation:", completed_this_invocation)
print("Invocation elapsed       :", format_duration(time.perf_counter() - batch_wall_start))
print("Executor event log       :", CELL32_V9_EVENTS_PATH)
print("Executor summary         :", CELL32_V9_SUMMARY_PATH)
print("Resource log             :", CELL32_RESOURCE_LOG_PATH)
print("Resource summary         :", CELL32_RESOURCE_SUMMARY_PATH)
print("Source TEST rows read    : FALSE")
print("External target rows read: FALSE")
print("Source-test evaluation   : FALSE")
print("Protocol-2 evaluation    : FALSE")
print("Scientific metrics shown : FALSE")
if trained_count == 680 and failed_count == 0:
    print("NEXT                      : CELL 33 - SOURCE TEST + 2040 EXTERNAL EVALUATIONS")
elif failure_records:
    print("NEXT                      : AUDIT THE FIRST FAILED WORKER LOG; DO NOT CONTINUE TO CELL 33")
else:
    print("NEXT                      : RERUN THE IDENTICAL CELL 32 v9 TO RESUME")
print("=" * 128)

if failure_records:
    run_id, rc, log_path, _ = failure_records[0]
    raise AssertionError(f"Cell 32 v9 worker failed: {run_id}; returncode={rc}; log={log_path}")
