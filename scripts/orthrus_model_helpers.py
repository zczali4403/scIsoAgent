# ------------------------------------------------------------------------------
# Script: orthrus_model_helpers.py
# Purpose: Helper CLI for Orthrus embedding inference, task prediction, and result aggregation.
# Note: Documentation/layout update only; analysis logic was not changed.
# ------------------------------------------------------------------------------

"""Helpers for Orthrus-based gene-specific isoform functional inference.

This helper supports four stages:
1. Extract one gene's transcript FASTA from a full transcript FASTA.
2. Infer Orthrus embeddings from the gene-specific FASTA.
3. Apply pretrained task heads to the embedding matrix.
4. Aggregate task outputs and render the aggregated plot as a PDF figure.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
from matplotlib.colors import LinearSegmentedColormap
from tqdm import tqdm

TASK_DISPLAY_NAMES = {
    "mrl-sugimoto": "Mean Ribosome Load",
    "rnahl-human": "RNA Half-life",
    "translation-efficiency-human": "Translation Efficiency",
    "prot-loc": "Protein Localization",
    "rna-lifecycle-ietswaart": "RNA Lifecycle",
    "rna-loc-fazal": "RNA Localization",
}

DEFAULT_PROBE_DIR = Path("/work/scIsoAgent/orthrus_trained_probes")
DEFAULT_LABEL_MAP_DIR = Path("/work/scIsoAgent/orthrus_label_maps")
DEFAULT_TASK_SPECS = [
    {
        "task_name": "mrl-sugimoto",
        "probe_filename": "mrl-sugimoto_probe.joblib",
        "task": "regression",
        "label_map_filename": "",
    },
    {
        "task_name": "rnahl-human",
        "probe_filename": "rnahl-human_probe.joblib",
        "task": "regression",
        "label_map_filename": "",
    },
    {
        "task_name": "translation-efficiency-human",
        "probe_filename": "translation-efficiency-human_probe.joblib",
        "task": "regression",
        "label_map_filename": "",
    },
    {
        "task_name": "prot-loc",
        "probe_filename": "prot-loc_probe.joblib",
        "task": "multilabel",
        "label_map_filename": "prot-loc_label_map.json",
    },
    {
        "task_name": "rna-lifecycle-ietswaart",
        "probe_filename": "rna-lifecycle-ietswaart_probe.joblib",
        "task": "multilabel",
        "label_map_filename": "rna-lifecycle-ietswaart_label_map.json",
    },
    {
        "task_name": "rna-loc-fazal",
        "probe_filename": "rna-loc-fazal_probe.joblib",
        "task": "multilabel",
        "label_map_filename": "rna-loc-fazal_label_map.json",
    },
]

# Brief: Create the parent directory for a Python output path when it is missing.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Ensures the output file parent directory exists; no analysis object is returned.
def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


# Brief: Read FASTA records as transcript ID and sequence pairs.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Loaded fasta object/table from disk or package resources.
def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_id: str | None = None
    current_seq: list[str] = []

    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq).upper()))
                current_id = line[1:].split()[0]
                if not current_id:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                current_seq = []
            else:
                if current_id is None:
                    raise ValueError(
                        f"Found sequence before first FASTA header at line {line_number}"
                    )
                current_seq.append(line)

    if current_id is not None:
        records.append((current_id, "".join(current_seq).upper()))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


# Brief: Write transcript ID and sequence pairs to a FASTA file.
# Args:
#   path: Filesystem path to create, read, or write.
#   records: FASTA records represented as transcript ID and sequence pairs.
# Returns: Writes the FASTA file; no object is returned.
def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    ensure_parent(path)
    with path.open("w") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


# Brief: Read and validate the transcript-to-gene mapping table used for gene lookup.
# Args:
#   path: Filesystem path to create, read, or write.
#   transcript_id_col: Transcript/isoform identifier or transcript-level data used by the analysis.
#   gene_id_col: Column name containing gene IDs in the mapping table.
#   gene_symbol_col: Column name containing gene symbols in the mapping table.
# Returns: Loaded mapping table object/table from disk or package resources.
def read_mapping_table(
    path: Path,
    transcript_id_col: str = "transcript_id",
    gene_id_col: str = "gene_id",
    gene_symbol_col: str = "gene_symbol",
) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = [transcript_id_col, gene_id_col, gene_symbol_col]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            "Mapping table is missing required columns: " + ", ".join(missing)
        )
    out = df[[transcript_id_col, gene_id_col, gene_symbol_col]].copy()
    out.columns = ["transcript_id", "gene_id", "gene_symbol"]
    out = out.drop_duplicates(subset=["transcript_id"])
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["gene_id"] = out["gene_id"].astype(str)
    out["gene_symbol"] = out["gene_symbol"].astype(str)
    return out


# Brief: Find all transcript records matching a gene symbol or gene ID.
# Args:
#   gene_query: Gene symbol or gene ID requested by the user.
#   mapping_df: Loaded transcript-to-gene mapping data frame.
# Returns: Resolved gene transcripts value or matching records.
def resolve_gene_transcripts(gene_query: str, mapping_df: pd.DataFrame) -> pd.DataFrame:
    matches = mapping_df[
        (mapping_df["gene_symbol"] == gene_query)
        | (mapping_df["gene_id"] == gene_query)
    ].copy()
    if matches.empty:
        raise ValueError(f"No transcripts found for gene query: {gene_query}")
    return matches.sort_values(["gene_symbol", "gene_id", "transcript_id"])


# Brief: Extract FASTA records for all transcripts belonging to one queried gene.
# Args:
#   mapping_csv: Transcript-to-gene mapping CSV used for gene-specific FASTA extraction.
#   input_fasta: Input transcript FASTA file.
#   gene_query: Gene symbol or gene ID requested by the user.
#   output_fasta: Output FASTA file to create.
#   transcript_id_col: Transcript/isoform identifier or transcript-level data used by the analysis.
#   gene_id_col: Column name containing gene IDs in the mapping table.
#   gene_symbol_col: Column name containing gene symbols in the mapping table.
# Returns: Path to the written gene-specific FASTA file.
def extract_gene_fasta(
    mapping_csv: Path,
    input_fasta: Path,
    gene_query: str,
    output_fasta: Path,
    transcript_id_col: str = "transcript_id",
    gene_id_col: str = "gene_id",
    gene_symbol_col: str = "gene_symbol",
) -> Path:
    mapping_df = read_mapping_table(
        mapping_csv,
        transcript_id_col=transcript_id_col,
        gene_id_col=gene_id_col,
        gene_symbol_col=gene_symbol_col,
    )
    gene_tx = resolve_gene_transcripts(gene_query, mapping_df)
    wanted = set(gene_tx["transcript_id"].astype(str))
    records = read_fasta(input_fasta)
    subset = [(record_id, seq) for record_id, seq in records if record_id in wanted]
    if not subset:
        raise ValueError(
            f"No FASTA records from {input_fasta} matched the transcripts for gene: {gene_query}"
        )
    write_fasta(output_fasta, subset)
    return output_fasta


class _MambaGenerationOutput(dict):
    # Brief: Initialize a compatibility dictionary/object used by patched model-generation code.
    # Args:
    #   kwargs: Compatibility keyword payload stored on the dict-like output object.
    # Returns: Initialized dict-like compatibility object with keyword values exposed as attributes.
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.__dict__.update(kwargs)


try:
    import transformers.generation as transformers_generation

    generation_output = getattr(
        transformers_generation,
        "GenerateDecoderOnlyOutput",
        _MambaGenerationOutput,
    )
    if not hasattr(transformers_generation, "GreedySearchDecoderOnlyOutput"):
        transformers_generation.GreedySearchDecoderOnlyOutput = generation_output
    if not hasattr(transformers_generation, "SampleDecoderOnlyOutput"):
        transformers_generation.SampleDecoderOnlyOutput = generation_output
    if not hasattr(transformers_generation, "TextStreamer"):
        transformers_generation.TextStreamer = object
except ImportError:
    pass

# Brief: Run the Orthrus/mRNA-bench model to generate one embedding per transcript sequence.
# Args:
#   records: FASTA records represented as transcript ID and sequence pairs.
#   model_name: mRNA-bench model family name used for embedding inference.
#   model_version: mRNA-bench model version used for embedding inference.
#   device: Torch device used for embedding inference.
#   mrna_bench_data_path: Optional local data directory for the mRNA-bench package.
# Returns: Inferred embeddings values.
def infer_embeddings(
    records: list[tuple[str, str]],
    model_name: str,
    model_version: str,
    device: str,
    mrna_bench_data_path: Path | None = None,
) -> np.ndarray:
    try:
        import mrna_bench as mb
    except ImportError as exc:
        raise ImportError(
            "mrna_bench is required for Orthrus embedding inference."
        ) from exc

    if mrna_bench_data_path is not None:
        mb.update_data_path(str(mrna_bench_data_path))

    mb_model = mb.load_model(
        model_name=model_name,
        model_version=model_version,
        device=torch.device(device),
    )
    embeddings = []
    for _, seq in tqdm(records, desc="Embedding"):
        seq = seq.replace("U", "T")
        emb = mb_model.embed_sequence(seq).detach().cpu().numpy()[0]
        embeddings.append(emb)
    return np.stack(embeddings)


# Brief: Write an embedding matrix to CSV with transcript IDs.
# Args:
#   path: Filesystem path to create, read, or write.
#   transcript_ids: Transcript identifiers aligned with sequence records or embeddings.
#   embeddings: Numeric transcript embedding matrix to write.
# Returns: Writes embedding CSV; no object is returned.
def write_embedding_csv(
    path: Path, transcript_ids: list[str], embeddings: np.ndarray
) -> None:
    ensure_parent(path)
    header = ["transcript_id"] + [f"embedding_{i}" for i in range(embeddings.shape[1])]
    with path.open("w") as handle:
        handle.write(",".join(header) + "\n")
        for transcript_id, embedding in zip(transcript_ids, embeddings):
            values = [transcript_id] + [f"{value:.8g}" for value in embedding]
            handle.write(",".join(values) + "\n")


# Brief: Read FASTA, infer Orthrus embeddings, and save NPZ/optional CSV outputs.
# Args:
#   fasta: Transcript FASTA file used for sequence-based analysis.
#   output: Output path for the current CLI subcommand.
#   csv_output: Optional CSV path for embedding export.
#   model_name: mRNA-bench model family name used for embedding inference.
#   model_version: mRNA-bench model version used for embedding inference.
#   device: Torch device used for embedding inference.
#   mrna_bench_data_path: Optional local data directory for the mRNA-bench package.
# Returns: Result object/table from the embedding inference step.
def run_embedding_inference(
    fasta: Path,
    output: Path,
    csv_output: Path | None,
    model_name: str,
    model_version: str,
    device: str,
    mrna_bench_data_path: Path | None,
) -> Path:
    records = read_fasta(fasta)
    transcript_ids = [record_id for record_id, _ in records]
    lengths = np.array([len(seq) for _, seq in records], dtype=np.int64)
    embeddings = infer_embeddings(
        records,
        model_name=model_name,
        model_version=model_version,
        device=device,
        mrna_bench_data_path=mrna_bench_data_path,
    )
    ensure_parent(output)
    np.savez_compressed(
        output,
        transcript_ids=np.array(transcript_ids),
        embeddings=embeddings,
        lengths=lengths,
    )
    if csv_output is not None:
        write_embedding_csv(csv_output, transcript_ids, embeddings)
    return output


# Brief: Load a pretrained Orthrus downstream probe and optional metadata.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Tuple of (model, metadata) loaded from a joblib probe file.
def _load_probe(path: Path) -> tuple[object, dict]:
    obj = joblib.load(path)
    if isinstance(obj, dict) and "model" in obj:
        return obj["model"], obj.get("metadata", {})
    return obj, {}


# Brief: Load embeddings and transcript IDs from an NPZ file.
# Args:
#   npz_path: NPZ path containing saved embeddings.
#   embedding_key: NPZ key containing the embedding matrix.
#   id_key: NPZ key containing transcript IDs.
# Returns: Tuple of (embedding_matrix, transcript_id_array).
def _load_embeddings(
    npz_path: Path, embedding_key: str | None, id_key: str
) -> tuple[np.ndarray, np.ndarray]:
    with np.load(npz_path, allow_pickle=False) as data:
        if embedding_key is not None:
            if embedding_key not in data:
                raise ValueError(
                    f"Embedding key '{embedding_key}' not found in {npz_path}."
                )
            embeddings = data[embedding_key]
        elif "embeddings" in data:
            embeddings = data["embeddings"]
        elif "embedding" in data:
            embeddings = data["embedding"]
        else:
            raise ValueError(
                f"No embedding key found in {npz_path}. Expected 'embeddings' or 'embedding'."
            )

        if id_key in data:
            transcript_ids = data[id_key].astype(str)
        elif "transcript_ids" in data:
            transcript_ids = data["transcript_ids"].astype(str)
        else:
            transcript_ids = np.array([str(i) for i in range(len(embeddings))])

    if embeddings.ndim != 2:
        raise ValueError(f"Expected 2D embedding matrix, got shape {embeddings.shape}.")
    if len(transcript_ids) != len(embeddings):
        raise ValueError(
            f"ID count {len(transcript_ids)} does not match embedding rows {len(embeddings)}."
        )
    return embeddings, transcript_ids


# Brief: Infer the prediction task type from probe metadata or model capabilities.
# Args:
#   model: Loaded prediction model/probe object.
#   metadata: Probe metadata dictionary.
#   task_override: Explicit task type overriding inferred probe metadata.
# Returns: Detected/inferred task label used by the Orthrus prediction workflow.
def _infer_task(model: object, metadata: dict, task_override: str | None) -> str:
    if task_override:
        return task_override
    if "task" in metadata:
        return metadata["task"]
    class_name = model.__class__.__name__.lower()
    if "multioutput" in class_name:
        return "multilabel"
    if hasattr(model, "predict_proba"):
        return "classification"
    return "regression"


# Brief: Load a label-index-to-name mapping JSON file.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Dictionary mapping class indices to labels, or None when no label map is provided.
def _load_label_map(path: Path | None) -> dict[int, str] | None:
    if path is None:
        return None
    with path.open("r") as handle:
        raw = json.load(handle)
    return {int(k): str(v) for k, v in raw.items()}


# Brief: Apply one pretrained probe to Orthrus embeddings and write prediction CSV output.
# Args:
#   probe_path: Path to a trained Orthrus downstream probe.
#   embedding_npz: NPZ file containing transcript embeddings.
#   output_csv: Prediction CSV path to write.
#   task: Prediction task type for a probe.
#   label_map_json: JSON file mapping label indices to class names.
#   embedding_key: NPZ key containing the embedding matrix.
#   id_key: NPZ key containing transcript IDs.
#   threshold: Probability threshold for multilabel positive calls.
#   save_metadata: Whether to write probe metadata alongside predictions.
# Returns: Path to the written probe prediction CSV file.
def predict_with_probe(
    probe_path: Path,
    embedding_npz: Path,
    output_csv: Path,
    task: str | None = None,
    label_map_json: Path | None = None,
    embedding_key: str | None = None,
    id_key: str = "transcript_ids",
    threshold: float = 0.5,
    save_metadata: bool = False,
) -> Path:
    model, metadata = _load_probe(probe_path)
    label_map = _load_label_map(label_map_json)
    X, transcript_ids = _load_embeddings(
        embedding_npz, embedding_key=embedding_key, id_key=id_key
    )

    if hasattr(model, "n_features_in_") and model.n_features_in_ != X.shape[1]:
        raise ValueError(
            "Embedding dimension mismatch between probe and input: "
            f"{model.n_features_in_} vs {X.shape[1]}"
        )

    resolved_task = _infer_task(model, metadata, task)
    out_df = pd.DataFrame({"transcript_id": transcript_ids})

    if resolved_task in ["regression", "reg_lin", "reg_ridge"]:
        y_pred = model.predict(X)
        if np.ndim(y_pred) == 1:
            out_df["prediction"] = y_pred
        else:
            n_targets = y_pred.shape[1]
            target_labels = metadata.get("target_labels") or [
                f"target_{i}" for i in range(n_targets)
            ]
            for i in range(n_targets):
                out_df[f"prediction_{target_labels[i]}"] = y_pred[:, i]

    elif resolved_task == "classification":
        pred = model.predict(X)
        out_df["prediction"] = pred
        if label_map is not None:
            out_df["prediction_label"] = [
                label_map.get(int(x), str(x)) for x in np.ravel(pred)
            ]
        if hasattr(model, "predict_proba"):
            probs = model.predict_proba(X)
            if np.ndim(probs) == 2 and probs.shape[1] > 1:
                out_df["score"] = probs[:, 1]
            else:
                out_df["score"] = np.ravel(probs)
        elif hasattr(model, "decision_function"):
            out_df["score"] = model.decision_function(X)

    elif resolved_task == "multilabel":
        probs = model.predict_proba(X)
        if isinstance(probs, list):
            score_mat = np.stack(
                [
                    p[:, 1] if p.ndim == 2 and p.shape[1] > 1 else np.ravel(p)
                    for p in probs
                ],
                axis=1,
            )
        else:
            score_mat = np.asarray(probs)
            if score_mat.ndim == 3 and score_mat.shape[2] > 1:
                score_mat = score_mat[:, :, 1]

        pred_mat = (score_mat >= threshold).astype(int)
        n_labels = score_mat.shape[1]
        target_labels = metadata.get("target_labels")
        if not isinstance(target_labels, list) or len(target_labels) != n_labels:
            target_labels = [f"label_{i}" for i in range(n_labels)]
        if label_map is not None:
            target_labels = [
                label_map.get(i, target_labels[i]) for i in range(n_labels)
            ]

        prediction_strings = []
        for row in pred_mat:
            positive_labels = [
                target_labels[i] for i, value in enumerate(row) if value >= 1
            ]
            prediction_strings.append(
                "|".join(positive_labels) if positive_labels else "None"
            )
        out_df["prediction"] = prediction_strings

        for i, label in enumerate(target_labels):
            out_df[f"pred_{label}"] = pred_mat[:, i]
            out_df[f"score_{label}"] = score_mat[:, i]
    else:
        raise ValueError(f"Unsupported task type: {resolved_task}")

    ensure_parent(output_csv)
    out_df.to_csv(output_csv, index=False)

    if save_metadata:
        meta_out = output_csv.with_suffix(output_csv.suffix + ".metadata.json")
        ensure_parent(meta_out)
        meta_out.write_text(json.dumps(metadata, indent=2, ensure_ascii=False))

    return output_csv


# Brief: Create the default Orthrus task manifest from packaged probes and label maps.
# Args:
#   output: Output path for the current CLI subcommand.
#   probe_dir: Directory containing packaged Orthrus probe files.
#   label_map_dir: Directory containing packaged Orthrus label-map files.
# Returns: Constructed default task manifest object/table.
def build_default_task_manifest(
    output: Path,
    probe_dir: Path = DEFAULT_PROBE_DIR,
    label_map_dir: Path = DEFAULT_LABEL_MAP_DIR,
) -> Path:
    rows = []
    for spec in DEFAULT_TASK_SPECS:
        label_map_path = ""
        if spec["label_map_filename"]:
            label_map_path = str(label_map_dir / spec["label_map_filename"])
        rows.append(
            {
                "task_name": spec["task_name"],
                "probe_path": str(probe_dir / spec["probe_filename"]),
                "task": spec["task"],
                "label_map_json": label_map_path,
                "threshold": 0.5,
                "save_metadata": False,
            }
        )
    ensure_parent(output)
    pd.DataFrame(rows).to_csv(output, index=False)
    return output


# Brief: Read and validate an Orthrus task manifest table.
# Args:
#   path: Filesystem path to create, read, or write.
# Returns: Loaded task manifest object/table from disk or package resources.
def load_task_manifest(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    required = ["task_name", "probe_path"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            "Task manifest is missing required columns: " + ", ".join(missing)
        )
    return df.fillna("")


# Brief: Run all Orthrus tasks listed in a manifest against one embedding file.
# Args:
#   task_manifest: CSV/TSV manifest describing Orthrus prediction tasks.
#   embedding_npz: NPZ file containing transcript embeddings.
#   output_dir: Directory where output files should be written.
#   embedding_key: NPZ key containing the embedding matrix.
#   id_key: NPZ key containing transcript IDs.
# Returns: List of prediction CSV paths produced by all manifest tasks.
def run_task_manifest(
    task_manifest: Path,
    embedding_npz: Path,
    output_dir: Path,
    embedding_key: str | None = None,
    id_key: str = "transcript_ids",
) -> list[Path]:
    manifest_df = load_task_manifest(task_manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    produced: list[Path] = []
    resolved_rows = []

    for row in manifest_df.to_dict(orient="records"):
        task_name = str(row["task_name"])
        probe_path = Path(str(row["probe_path"]))
        task_type = str(row.get("task", "")).strip() or None
        label_map_json = str(row.get("label_map_json", "")).strip()
        threshold = float(row.get("threshold", 0.5) or 0.5)
        save_metadata = str(row.get("save_metadata", "")).strip().lower() in {
            "1",
            "true",
            "yes",
        }
        output_csv = output_dir / f"{task_name}_predictions.csv"
        predict_with_probe(
            probe_path=probe_path,
            embedding_npz=embedding_npz,
            output_csv=output_csv,
            task=task_type,
            label_map_json=Path(label_map_json) if label_map_json else None,
            embedding_key=embedding_key,
            id_key=id_key,
            threshold=threshold,
            save_metadata=save_metadata,
        )
        produced.append(output_csv)
        resolved_rows.append(
            {
                "task_name": task_name,
                "probe_path": str(probe_path),
                "task": task_type or "",
                "label_map_json": label_map_json,
                "threshold": threshold,
                "output_csv": str(output_csv),
            }
        )

    pd.DataFrame(resolved_rows).to_csv(
        output_dir / "resolved_task_manifest.csv", index=False
    )
    return produced


# Brief: Detect whether a prediction table contains regression, classification, or multilabel output.
# Args:
#   df: Prediction table inspected for task-type detection.
# Returns: Detected/inferred task label used by the Orthrus prediction workflow.
def _detect_task_type(df: pd.DataFrame) -> str:
    has_pred_cols = any(c.startswith("pred_") for c in df.columns)
    has_score_cols = any(c.startswith("score_") for c in df.columns)
    if has_pred_cols and has_score_cols:
        return "multilabel"
    if "prediction_label" in df.columns:
        return "classification"
    if "prediction" not in df.columns:
        return "unknown"
    numeric = pd.to_numeric(df["prediction"], errors="coerce")
    if numeric.notna().all():
        return "regression"
    return "classification"


# Brief: Min-max scale numeric predictions for heatmap visualization.
# Args:
#   col: Numeric pandas Series to scale into the [0, 1] interval.
# Returns: Numeric values scaled into the [0, 1] interval.
def _minmax_scale(col: pd.Series) -> pd.Series:
    mn = col.min()
    mx = col.max()
    if pd.isna(mn) or pd.isna(mx):
        return pd.Series([0.0] * len(col), index=col.index)
    if mx == mn:
        return pd.Series([0.5] * len(col), index=col.index)
    return (col - mn) / (mx - mn)


# Brief: Sort Orthrus tasks in a preferred display order.
# Args:
#   task_name: Orthrus task identifier used for preferred report ordering.
# Returns: Sort key tuple used to order Orthrus tasks in reports.
def _task_sort_key(task_name: str) -> tuple[int, str]:
    preferred = ["prot-loc", "rna-lifecycle-ietswaart", "rna-loc-fazal"]
    if task_name in preferred:
        return (preferred.index(task_name), task_name)
    return (len(preferred), task_name)


# Brief: Build a transcript-by-class score matrix from multilabel prediction output.
# Args:
#   task_df: Per-task Orthrus prediction table containing transcript IDs and class columns.
#   transcript_order: Transcript/isoform identifier or transcript-level data used by the analysis.
# Returns: DataFrame matrix assembled for Orthrus report plotting.
def _build_multilabel_score_matrix(
    task_df: pd.DataFrame, transcript_order: list[str]
) -> pd.DataFrame:
    score_cols = [c for c in task_df.columns if c.startswith("score_")]
    pred_cols = [c for c in task_df.columns if c.startswith("pred_")]
    if not score_cols and not pred_cols:
        return pd.DataFrame(index=transcript_order)
    if score_cols:
        value_cols = score_cols
        rename = {c: c[len("score_") :] for c in score_cols}
    else:
        value_cols = pred_cols
        rename = {c: c[len("pred_") :] for c in pred_cols}
    score_df = task_df[["transcript_id", *value_cols]].copy()
    score_df["transcript_id"] = score_df["transcript_id"].astype(str)
    score_df = score_df.set_index("transcript_id")
    score_df = score_df.rename(columns=rename)
    score_df = score_df.apply(pd.to_numeric, errors="coerce")
    score_df = score_df.reindex(transcript_order).fillna(0.0)
    return score_df


# Brief: Build a transcript-by-class binary prediction matrix from multilabel output.
# Args:
#   task_df: Per-task Orthrus prediction table containing transcript IDs and class columns.
#   transcript_order: Transcript/isoform identifier or transcript-level data used by the analysis.
# Returns: DataFrame matrix assembled for Orthrus report plotting.
def _build_multilabel_pred_matrix(
    task_df: pd.DataFrame, transcript_order: list[str]
) -> pd.DataFrame:
    pred_cols = [c for c in task_df.columns if c.startswith("pred_")]
    if not pred_cols:
        return pd.DataFrame(index=transcript_order)
    pred_df = task_df[["transcript_id", *pred_cols]].copy()
    pred_df["transcript_id"] = pred_df["transcript_id"].astype(str)
    pred_df = pred_df.set_index("transcript_id")
    pred_df = pred_df.rename(columns={c: c[len("pred_") :] for c in pred_cols})
    pred_df = pred_df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    pred_df = pred_df.reindex(transcript_order).fillna(0.0)
    return pred_df


# Brief: Merge Orthrus task outputs, write summary tables, and render the PDF report.
# Args:
#   result_dir: Directory containing per-task Orthrus prediction CSV files.
#   output_table: CSV path for the aggregated Orthrus prediction table.
#   output_category_map: JSON path for categorical-code mappings.
#   output_pdf: PDF path for the aggregated Orthrus report.
#   task_files: Optional explicit task prediction filenames to aggregate.
#   join: Join mode used when merging per-task prediction tables.
# Returns: Tuple of output paths: aggregated table, category map JSON, and PDF figure.
def aggregate_results_to_pdf(
    result_dir: Path,
    output_table: Path,
    output_category_map: Path,
    output_pdf: Path,
    task_files: list[str] | None = None,
    join: str = "inner",
) -> tuple[Path, Path, Path]:
    if not result_dir.exists():
        raise FileNotFoundError(f"Result directory not found: {result_dir}")

    if task_files:
        files = [result_dir / name for name in task_files]
    else:
        files = sorted(result_dir.glob("*_predictions.csv"))
    files = [f for f in files if f.is_file() and f.name != output_table.name]
    if not files:
        raise FileNotFoundError(f"No *_predictions.csv found in {result_dir}")

    merged_df: pd.DataFrame | None = None
    task_types: dict[str, str] = {}
    category_code_map: dict[str, dict[str, int]] = {}
    task_raw_frames: dict[str, pd.DataFrame] = {}

    for file_path in files:
        task_name = file_path.name.replace("_predictions.csv", "")
        df = pd.read_csv(file_path, keep_default_na=False)
        if "transcript_id" not in df.columns:
            continue
        task_type = _detect_task_type(df)
        task_types[task_name] = task_type
        task_raw_frames[task_name] = df.copy()

        if task_type == "regression":
            task_col = pd.to_numeric(df["prediction"], errors="coerce")
        elif task_type == "classification":
            if "prediction_label" in df.columns:
                task_col = df["prediction_label"].astype(str)
            else:
                task_col = df["prediction"].astype(str)
        elif task_type == "multilabel":
            task_col = df["prediction"].astype(str)
        else:
            if "prediction" not in df.columns:
                continue
            task_col = df["prediction"].astype(str)

        slim = pd.DataFrame(
            {"transcript_id": df["transcript_id"].astype(str), task_name: task_col}
        )
        if merged_df is None:
            merged_df = slim
        else:
            merged_df = merged_df.merge(slim, on="transcript_id", how=join)

    if merged_df is None:
        raise RuntimeError("No valid task prediction files to aggregate.")

    merged_df = merged_df.sort_values("transcript_id").reset_index(drop=True)
    ensure_parent(output_table)
    merged_df.to_csv(output_table, index=False)

    task_cols = [c for c in merged_df.columns if c != "transcript_id"]
    for task_name in task_cols:
        task_type = task_types.get(task_name, "unknown")
        if task_type == "regression":
            continue
        values = merged_df[task_name].fillna("NA").astype(str)
        unique_vals = sorted(values.unique())
        category_code_map[task_name] = {val: i for i, val in enumerate(unique_vals)}

    ensure_parent(output_category_map)
    output_category_map.write_text(
        json.dumps(category_code_map, indent=2, ensure_ascii=False)
    )

    transcript_order = merged_df["transcript_id"].astype(str).tolist()
    reg_tasks = [t for t in task_cols if task_types.get(t) == "regression"]
    reg_df = pd.DataFrame(index=transcript_order)
    for task_name in reg_tasks:
        numeric = pd.to_numeric(merged_df[task_name], errors="coerce")
        reg_df[TASK_DISPLAY_NAMES.get(task_name, task_name)] = _minmax_scale(
            numeric
        ).to_numpy()

    multilabel_tasks = sorted(
        [t for t in task_cols if task_types.get(t) == "multilabel"],
        key=_task_sort_key,
    )[:3]

    sns.set_theme(style="white", context="talk")
    reg_cmap = LinearSegmentedColormap.from_list("reg_cmap", ["#d6e8f5", "#1f4e79"])
    cls_cmap = LinearSegmentedColormap.from_list("cls_cmap", ["#f4d9de", "#8c4257"])
    fig, axes = plt.subplots(2, 2, figsize=(22, 16))
    axes = axes.flatten()

    ax0 = axes[0]
    if reg_df.shape[1] > 0:
        sns.heatmap(
            reg_df,
            ax=ax0,
            cmap=reg_cmap,
            vmin=0.0,
            vmax=1.0,
            linewidths=0.35,
            linecolor="white",
            cbar_kws={"label": "Scaled Score"},
        )
        ax0.set_title("Regression Tasks (Min-Max Scaled)")
        ax0.set_xlabel("Task")
        ax0.set_ylabel("Transcript")
    else:
        ax0.axis("off")
        ax0.text(0.5, 0.5, "No regression tasks", ha="center", va="center")

    for i in range(1, 4):
        ax = axes[i]
        if i - 1 < len(multilabel_tasks):
            task_name = multilabel_tasks[i - 1]
            task_df = task_raw_frames[task_name]
            score_df = _build_multilabel_score_matrix(task_df, transcript_order)
            pred_df = _build_multilabel_pred_matrix(task_df, transcript_order)
            if score_df.shape[1] == 0:
                ax.axis("off")
                ax.text(
                    0.5,
                    0.5,
                    f"{task_name}\n(no score columns)",
                    ha="center",
                    va="center",
                )
                continue
            vmax = float(np.nanmax(score_df.to_numpy())) if score_df.size else 1.0
            if not np.isfinite(vmax) or vmax <= 0:
                vmax = 1.0
            sns.heatmap(
                score_df,
                ax=ax,
                cmap=cls_cmap,
                vmin=0.0,
                vmax=vmax,
                linewidths=0.35,
                linecolor="white",
                cbar_kws={"label": "Class Score"},
            )
            if pred_df.shape == score_df.shape:
                pred_arr = pred_df.to_numpy()
                for row_idx in range(pred_arr.shape[0]):
                    for col_idx in range(pred_arr.shape[1]):
                        if pred_arr[row_idx, col_idx] >= 0.5:
                            ax.text(
                                col_idx + 0.5,
                                row_idx + 0.5,
                                "*",
                                ha="center",
                                va="center",
                                color="black",
                                fontsize=11,
                                fontweight="bold",
                            )
            ax.set_title(f"{TASK_DISPLAY_NAMES.get(task_name, task_name)} (Multilabel)")
            ax.set_xlabel("Class")
            ax.set_ylabel("Transcript")
        else:
            ax.axis("off")

    for ax in axes:
        if ax.has_data():
            ax.tick_params(axis="x", labelrotation=35)
            ax.tick_params(axis="y", labelrotation=0, labelsize=10)

    plt.suptitle("Orthrus Multi-Task Isoform Prediction Summary", fontsize=22, y=0.995)
    plt.tight_layout(rect=(0, 0, 1, 0.985))
    ensure_parent(output_pdf)
    plt.savefig(output_pdf, format="pdf")
    plt.close(fig)
    return output_table, output_category_map, output_pdf


# Brief: Define the Orthrus helper command-line interface and subcommands.
# Args: None.
# Returns: Constructed parser object/table.
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Orthrus model helper utilities.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_extract = subparsers.add_parser(
        "extract-gene-fasta", help="Extract one gene's transcripts into a FASTA."
    )
    p_extract.add_argument("--mapping-csv", type=Path, required=True)
    p_extract.add_argument("--input-fasta", type=Path, required=True)
    p_extract.add_argument("--gene-query", required=True)
    p_extract.add_argument("--output-fasta", type=Path, required=True)
    p_extract.add_argument("--transcript-id-col", default="transcript_id")
    p_extract.add_argument("--gene-id-col", default="gene_id")
    p_extract.add_argument("--gene-symbol-col", default="gene_symbol")

    p_embed = subparsers.add_parser(
        "infer-embeddings", help="Infer Orthrus embeddings from a FASTA."
    )
    p_embed.add_argument("--fasta", type=Path, required=True)
    p_embed.add_argument("--output", type=Path, required=True)
    p_embed.add_argument("--csv-output", type=Path, default=None)
    p_embed.add_argument("--mb-model-name", default="Orthrus")
    p_embed.add_argument("--mb-model-version", default="orthrus-base-4-track")
    p_embed.add_argument("--mrna-bench-data-path", type=Path, default=None)
    p_embed.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Torch device, for example cuda, cuda:0, or cpu.",
    )

    p_manifest = subparsers.add_parser(
        "build-default-manifest",
        help="Build a task manifest from the packaged Orthrus probe and label-map directories.",
    )
    p_manifest.add_argument("--output", type=Path, required=True)
    p_manifest.add_argument("--probe-dir", type=Path, default=DEFAULT_PROBE_DIR)
    p_manifest.add_argument("--label-map-dir", type=Path, default=DEFAULT_LABEL_MAP_DIR)

    p_predict = subparsers.add_parser(
        "predict-tasks", help="Run all tasks from a manifest on one embedding NPZ."
    )
    p_predict.add_argument("--embedding-npz", type=Path, required=True)
    p_predict.add_argument("--task-manifest", type=Path, required=True)
    p_predict.add_argument("--output-dir", type=Path, required=True)
    p_predict.add_argument("--embedding-key", default=None)
    p_predict.add_argument("--id-key", default="transcript_ids")

    p_agg = subparsers.add_parser(
        "aggregate-results",
        help="Aggregate task outputs and render the aggregated plot as a PDF figure.",
    )
    p_agg.add_argument("--result-dir", type=Path, required=True)
    p_agg.add_argument("--output-table", type=Path, required=True)
    p_agg.add_argument("--output-category-map", type=Path, required=True)
    p_agg.add_argument("--output-pdf", type=Path, required=True)
    p_agg.add_argument("--task-files", nargs="*", default=None)
    p_agg.add_argument("--join", choices=["inner", "outer"], default="inner")

    return parser


# Brief: Run the script-specific command-line workflow.
# Args: None.
# Returns: No return value; executes the selected workflow and writes configured outputs.
def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "extract-gene-fasta":
        output = extract_gene_fasta(
            mapping_csv=args.mapping_csv,
            input_fasta=args.input_fasta,
            gene_query=args.gene_query,
            output_fasta=args.output_fasta,
            transcript_id_col=args.transcript_id_col,
            gene_id_col=args.gene_id_col,
            gene_symbol_col=args.gene_symbol_col,
        )
        print(f"Saved gene-specific FASTA: {output}")
        return

    if args.command == "infer-embeddings":
        output = run_embedding_inference(
            fasta=args.fasta,
            output=args.output,
            csv_output=args.csv_output,
            model_name=args.mb_model_name,
            model_version=args.mb_model_version,
            device=args.device,
            mrna_bench_data_path=args.mrna_bench_data_path,
        )
        print(f"Saved embeddings: {output}")
        if args.csv_output is not None:
            print(f"Saved embedding CSV: {args.csv_output}")
        return

    if args.command == "build-default-manifest":
        output = build_default_task_manifest(
            output=args.output,
            probe_dir=args.probe_dir,
            label_map_dir=args.label_map_dir,
        )
        print(f"Saved task manifest: {output}")
        return

    if args.command == "predict-tasks":
        produced = run_task_manifest(
            task_manifest=args.task_manifest,
            embedding_npz=args.embedding_npz,
            output_dir=args.output_dir,
            embedding_key=args.embedding_key,
            id_key=args.id_key,
        )
        print("Saved task prediction CSVs:")
        for path in produced:
            print(path)
        return

    if args.command == "aggregate-results":
        output_table, output_category_map, output_pdf = aggregate_results_to_pdf(
            result_dir=args.result_dir,
            output_table=args.output_table,
            output_category_map=args.output_category_map,
            output_pdf=args.output_pdf,
            task_files=args.task_files,
            join=args.join,
        )
        print(f"Saved aggregated table: {output_table}")
        print(f"Saved category map: {output_category_map}")
        print(f"Saved PDF figure: {output_pdf}")
        return

    raise ValueError(f"Unsupported command: {args.command}")


# Main CLI flow.
if __name__ == "__main__":
    main()
