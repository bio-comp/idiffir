"""Utilities for migrating legacy classifier archives to sklearn artifacts."""

from __future__ import annotations

import configparser
import json
import tempfile
import zipfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Protocol, TypedDict, TypeVar

import numpy as np

from iDiffIR.SpliceGrapher.predict.model_schema import (
    ArtifactFormat,
    ClassifierTask,
    FeatureMapping,
    ModelMetadata,
    save_model_metadata,
)


class LegacyCompareStatus(str, Enum):
    """Report status for legacy fixture comparison."""

    COMPLETED = "completed"
    MISSING_FIXTURE = "missing_fixture"
    INVALID_FIXTURE = "invalid_fixture"


class LabelBalance(TypedDict):
    """Binary class-count report for migrated models."""

    zeros: int
    ones: int


class ModelMigrationReport(TypedDict, total=False):
    """Per-model migration report schema."""

    name: str
    records: int
    label_balance: LabelBalance
    artifact_path: str
    metadata_path: str
    scaled_accuracy: float
    unscaled_accuracy: float
    accuracy_delta_scaled_minus_unscaled: float
    legacy_compare_status: str
    legacy_compare_reason: str
    legacy_fixture_records: int
    legacy_fixture_agreement: float


class ArchiveMigrationReport(TypedDict):
    """Per-archive migration report schema."""

    archive: str
    species: str
    models: list[ModelMigrationReport]


class MigrationSummary(TypedDict):
    """Aggregate migration summary."""

    archive_count: int
    model_count: int


class AggregateMigrationReport(TypedDict):
    """Aggregate report schema covering all migrated archives."""

    archives: list[ArchiveMigrationReport]
    summary: MigrationSummary


PipelineT = TypeVar("PipelineT", bound="SequenceClassifier")


class SequenceClassifier(Protocol):
    """Minimal classifier contract used by migration training and reports."""

    def fit(self: PipelineT, sequences: list[str], labels: np.ndarray) -> PipelineT:
        """Train the classifier from sequence features and labels."""

    def predict(self, sequences: list[str]) -> np.ndarray:
        """Predict labels for sequence features."""


@dataclass(frozen=True)
class ClassifierConfig:
    """Normalized subset of legacy classifier config values needed for retraining."""

    cfg_name: str
    fasta_path: str
    svm_path: str
    dimer: str
    c_value: float
    threshold: float
    mink: int
    maxk: int
    degree: int
    gamma: float
    acceptor: bool
    normalize: bool


@dataclass(frozen=True)
class TrainingResult:
    """Output values from training scaled and unscaled migration pipelines."""

    scaled_pipeline: SequenceClassifier
    scaled_accuracy: float
    unscaled_accuracy: float

    @property
    def accuracy_delta_scaled_minus_unscaled(self) -> float:
        """Return scaled-minus-unscaled training accuracy delta."""
        return self.scaled_accuracy - self.unscaled_accuracy


@dataclass(frozen=True)
class LegacyFixture:
    """Reference fixture payload for old/new output comparison."""

    sequences: list[str]
    legacy_predictions: np.ndarray


@dataclass(frozen=True)
class LegacyFixtureComparison:
    """Comparison metrics for migrated model vs legacy fixture predictions."""

    status: LegacyCompareStatus
    records: int
    agreement: float
    reason: str = ""


def _require_ml_stack() -> None:
    """Validate required ML dependencies and raise actionable import guidance."""
    try:
        import joblib  # noqa: F401
        import sklearn  # noqa: F401
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "scikit-learn and joblib are required for classifier migration. "
            "Install them with `uv sync --group dev`."
        ) from exc


def load_classifier_config(path: Path) -> ClassifierConfig:
    """Load one legacy classifier cfg into a typed config object."""
    parser = configparser.ConfigParser()
    parser.read(path)
    section = parser["SpliceGrapherConfig"]

    return ClassifierConfig(
        cfg_name=path.name,
        fasta_path=section.get("fasta_path", ""),
        svm_path=section.get("svm_path", ""),
        dimer=section.get("dimer", ""),
        c_value=section.getfloat("c", fallback=1.0),
        threshold=section.getfloat("threshold", fallback=0.0),
        mink=section.getint("mink", fallback=1),
        maxk=section.getint("maxk", fallback=1),
        degree=section.getint("degree", fallback=1),
        gamma=section.getfloat("gamma", fallback=-1.0),
        acceptor=section.getboolean("acceptor", fallback=False),
        normalize=section.getboolean("normalize", fallback=True),
    )


def load_labeled_fasta(path: Path) -> tuple[list[str], np.ndarray]:
    """Load labeled FASTA records where headers include a `label=<0|1>` token."""
    sequences: list[str] = []
    labels: list[int] = []

    header = ""
    sequence_parts: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header and sequence_parts:
                sequences.append("".join(sequence_parts).upper())
                labels.append(_extract_label(header))
            header = line[1:]
            sequence_parts = []
            continue
        sequence_parts.append(line)

    if header and sequence_parts:
        sequences.append("".join(sequence_parts).upper())
        labels.append(_extract_label(header))

    if not sequences:
        raise ValueError(f"No labeled FASTA records found in {path}")
    if len(set(labels)) < 2:
        observed_labels = sorted(set(labels))
        raise ValueError(
            f"Expected at least two classes in {path}, got labels {observed_labels}"
        )

    return sequences, np.asarray(labels, dtype=np.int64)


def _extract_label(header: str) -> int:
    """Extract integer label value from a FASTA header."""
    for token in header.split():
        if token.startswith("label="):
            value = token.split("=", 1)[1]
            if value not in {"0", "1"}:
                raise ValueError(f"Unsupported label token `{token}` in `{header}`")
            return int(value)
    raise ValueError(f"Missing `label=` token in FASTA header: {header}")


def _kernel_settings(config: ClassifierConfig) -> dict[str, str | float | int]:
    """Translate legacy config hints to sklearn SVC kernel parameters."""
    if config.degree > 1:
        return {
            "kernel": "poly",
            "degree": config.degree,
            "gamma": config.gamma if config.gamma > 0 else "scale",
        }
    if config.gamma > 0:
        return {"kernel": "rbf", "gamma": config.gamma}
    return {"kernel": "linear"}


def build_pipeline(config: ClassifierConfig, include_scaler: bool) -> SequenceClassifier:
    """Build a sklearn pipeline for legacy sequence classification."""
    _require_ml_stack()
    from sklearn.feature_extraction.text import CountVectorizer
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC

    steps = [
        (
            "vectorizer",
            CountVectorizer(
                analyzer="char",
                ngram_range=(config.mink, config.maxk),
                lowercase=False,
            ),
        )
    ]
    if include_scaler:
        # Sparse-safe scaling preserves sparse matrix shape during normalization.
        steps.append(("scaler", StandardScaler(with_mean=False)))
    steps.append(
        (
            "svc",
            SVC(
                C=config.c_value,
                **_kernel_settings(config),
            ),
        )
    )
    return Pipeline(steps)


def train_and_score(
    config: ClassifierConfig,
    sequences: list[str],
    labels: np.ndarray,
) -> TrainingResult:
    """Train scaled and unscaled pipelines and return key migration metrics."""
    _require_ml_stack()
    from sklearn.metrics import accuracy_score

    scaled_pipeline = build_pipeline(config, include_scaler=True)
    scaled_pipeline.fit(sequences, labels)
    scaled_predictions = scaled_pipeline.predict(sequences)
    scaled_accuracy = float(accuracy_score(labels, scaled_predictions))

    unscaled_pipeline = build_pipeline(config, include_scaler=False)
    unscaled_pipeline.fit(sequences, labels)
    unscaled_predictions = unscaled_pipeline.predict(sequences)
    unscaled_accuracy = float(accuracy_score(labels, unscaled_predictions))

    return TrainingResult(
        scaled_pipeline=scaled_pipeline,
        scaled_accuracy=scaled_accuracy,
        unscaled_accuracy=unscaled_accuracy,
    )


def _fixture_path_for_cfg(cfg_path: Path) -> Path:
    """Return the expected legacy fixture path for a cfg file."""
    return cfg_path.with_name(f"{cfg_path.stem}.legacy_fixture.json")


def _load_legacy_fixture(path: Path) -> LegacyFixture:
    """Load and validate a JSON legacy fixture for model-output comparison."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    schema_version = str(payload.get("schema_version", ""))
    if schema_version != "1":
        raise ValueError(f"Unsupported fixture schema version `{schema_version}` in {path}")

    sequences = payload.get("sequences")
    legacy_predictions = payload.get("legacy_predictions")
    if not isinstance(sequences, list) or not isinstance(legacy_predictions, list):
        raise ValueError(f"Invalid fixture payload structure in {path}")
    if len(sequences) != len(legacy_predictions):
        raise ValueError(f"Mismatched sequence/prediction counts in {path}")

    normalized_sequences: list[str] = []
    normalized_predictions: list[int] = []
    for sequence in sequences:
        if not isinstance(sequence, str):
            raise ValueError(f"Fixture sequence must be text in {path}")
        normalized_sequences.append(sequence.upper())

    for prediction in legacy_predictions:
        if prediction not in {0, 1}:
            raise ValueError(f"Fixture legacy prediction must be 0/1 in {path}")
        normalized_predictions.append(int(prediction))

    return LegacyFixture(
        sequences=normalized_sequences,
        legacy_predictions=np.asarray(normalized_predictions, dtype=np.int64),
    )


def _compare_with_fixture(
    pipeline: SequenceClassifier,
    fixture: LegacyFixture,
) -> LegacyFixtureComparison:
    """Compare migrated model predictions with legacy fixture predictions."""
    _require_ml_stack()
    from sklearn.metrics import accuracy_score

    predicted = pipeline.predict(fixture.sequences)
    agreement = float(accuracy_score(fixture.legacy_predictions, predicted))
    return LegacyFixtureComparison(
        status=LegacyCompareStatus.COMPLETED,
        records=int(fixture.legacy_predictions.shape[0]),
        agreement=agreement,
    )


def _compare_legacy_fixture(
    extracted_root: Path,
    cfg_path: Path,
    pipeline: SequenceClassifier,
) -> LegacyFixtureComparison:
    """Evaluate legacy fixture agreement when fixture exists for a cfg model."""
    fixture_path = extracted_root / _fixture_path_for_cfg(cfg_path).name
    if not fixture_path.exists():
        return LegacyFixtureComparison(
            status=LegacyCompareStatus.MISSING_FIXTURE,
            records=0,
            agreement=0.0,
            reason="No legacy fixture provided in archive.",
        )

    try:
        fixture = _load_legacy_fixture(fixture_path)
    except ValueError as exc:
        return LegacyFixtureComparison(
            status=LegacyCompareStatus.INVALID_FIXTURE,
            records=0,
            agreement=0.0,
            reason=str(exc),
        )
    return _compare_with_fixture(pipeline, fixture)


def migrate_archive(archive_path: Path, output_root: Path) -> ArchiveMigrationReport:
    """Migrate one legacy zip archive into joblib models plus metadata."""
    _require_ml_stack()
    from joblib import dump
    from sklearn import __version__ as sklearn_version

    if not archive_path.exists():
        raise FileNotFoundError(archive_path)
    if archive_path.suffix != ".zip":
        raise ValueError(f"Expected .zip archive, got {archive_path}")

    report_models: list[ModelMigrationReport] = []
    species_dir = output_root / archive_path.stem
    species_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"{archive_path.stem}_") as tmp_dir:
        with zipfile.ZipFile(archive_path, "r") as bundle:
            bundle.extractall(tmp_dir)
        extracted_root = Path(tmp_dir)
        cfg_files = sorted(extracted_root.glob("*.cfg"))
        if not cfg_files:
            raise ValueError(f"No .cfg files found in {archive_path}")

        for cfg_path in cfg_files:
            config = load_classifier_config(cfg_path)
            fasta_path = extracted_root / config.fasta_path
            if not fasta_path.exists():
                raise FileNotFoundError(
                    "Missing FASTA "
                    f"`{config.fasta_path}` referenced by `{cfg_path.name}` in {archive_path}"
                )

            sequences, labels = load_labeled_fasta(fasta_path)
            training_result = train_and_score(config, sequences, labels)

            stem = cfg_path.stem
            model_path = species_dir / f"{stem}.joblib"
            dump(training_result.scaled_pipeline, model_path)

            metadata = ModelMetadata(
                classifier_task=_classifier_task(config),
                artifact_format=ArtifactFormat.JOBLIB,
                artifact_path=model_path.name,
                sklearn_version=sklearn_version,
                feature_extractor="char_ngram_count_v1",
                feature_mappings=(
                    FeatureMapping(
                        legacy_name="sequence_windows",
                        sklearn_name=f"char_ngram[{config.mink}:{config.maxk}]",
                        notes="Legacy sequence windows encoded as character n-gram counts.",
                    ),
                    FeatureMapping(
                        legacy_name="normalize",
                        sklearn_name="standard_scaler_with_mean_false",
                        notes="Scaling is explicitly persisted in the sklearn pipeline.",
                    ),
                ),
                threshold=config.threshold,
                legacy_config_path=config.cfg_name,
                legacy_svm_path=config.svm_path or None,
            )
            metadata_path = species_dir / f"{stem}.metadata.json"
            save_model_metadata(metadata, metadata_path)

            legacy_compare = _compare_legacy_fixture(
                extracted_root=extracted_root,
                cfg_path=cfg_path,
                pipeline=training_result.scaled_pipeline,
            )
            model_report: ModelMigrationReport = {
                "name": stem,
                "records": int(labels.shape[0]),
                "label_balance": {
                    "zeros": int(np.sum(labels == 0)),
                    "ones": int(np.sum(labels == 1)),
                },
                "artifact_path": str(model_path),
                "metadata_path": str(metadata_path),
                "scaled_accuracy": training_result.scaled_accuracy,
                "unscaled_accuracy": training_result.unscaled_accuracy,
                "accuracy_delta_scaled_minus_unscaled": (
                    training_result.accuracy_delta_scaled_minus_unscaled
                ),
                "legacy_compare_status": legacy_compare.status.value,
            }

            if legacy_compare.status == LegacyCompareStatus.COMPLETED:
                model_report["legacy_fixture_records"] = legacy_compare.records
                model_report["legacy_fixture_agreement"] = legacy_compare.agreement
            else:
                model_report["legacy_compare_reason"] = legacy_compare.reason

            report_models.append(model_report)

    return {
        "archive": str(archive_path),
        "species": archive_path.stem,
        "models": report_models,
    }


def migrate_archives(
    archive_paths: list[Path],
    output_root: Path,
    report_path: Path,
) -> AggregateMigrationReport:
    """Migrate multiple archives and persist an aggregate report."""
    if not archive_paths:
        raise ValueError("No classifier archives were provided.")

    output_root.mkdir(parents=True, exist_ok=True)
    report_entries = [migrate_archive(path, output_root=output_root) for path in archive_paths]
    report: AggregateMigrationReport = {
        "archives": report_entries,
        "summary": {
            "archive_count": len(report_entries),
            "model_count": sum(len(entry["models"]) for entry in report_entries),
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def _classifier_task(config: ClassifierConfig) -> ClassifierTask:
    """Map legacy classifier cfg to target task enum."""
    if "gap" in config.cfg_name.lower():
        return ClassifierTask.GAP
    if config.acceptor:
        return ClassifierTask.SITE_ACCEPTOR
    return ClassifierTask.SITE_DONOR
