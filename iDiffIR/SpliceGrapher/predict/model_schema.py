"""Schema helpers for PyML -> scikit-learn classifier migration artifacts."""

from __future__ import annotations

import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "1"


class ClassifierTask(str, Enum):
    """Supported task types for migrated classifier artifacts."""

    SITE_ACCEPTOR = "site_acceptor"
    SITE_DONOR = "site_donor"
    GAP = "gap"


class ArtifactFormat(str, Enum):
    """Supported serialization formats for migrated classifier assets."""

    JOBLIB = "joblib"


@dataclass(frozen=True)
class FeatureMapping:
    """Maps one legacy feature name to the new extractor output name."""

    legacy_name: str
    sklearn_name: str
    notes: str = ""

    def to_dict(self) -> dict[str, str]:
        """Serializes this feature mapping to a dict."""
        return {
            "legacy_name": self.legacy_name,
            "sklearn_name": self.sklearn_name,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> FeatureMapping:
        """Builds a feature mapping from a dict payload."""
        return cls(
            legacy_name=str(value["legacy_name"]),
            sklearn_name=str(value["sklearn_name"]),
            notes=str(value.get("notes", "")),
        )


@dataclass(frozen=True)
class ModelMetadata:
    """Normalized metadata stored next to migrated model artifacts."""

    classifier_task: ClassifierTask
    artifact_format: ArtifactFormat
    artifact_path: str
    sklearn_version: str
    feature_extractor: str
    feature_mappings: tuple[FeatureMapping, ...]
    threshold: float | None = None
    legacy_config_path: str | None = None
    legacy_svm_path: str | None = None
    schema_version: str = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Serializes model metadata to a JSON-compatible dict."""
        return {
            "schema_version": self.schema_version,
            "classifier_task": self.classifier_task.value,
            "artifact_format": self.artifact_format.value,
            "artifact_path": self.artifact_path,
            "sklearn_version": self.sklearn_version,
            "feature_extractor": self.feature_extractor,
            "feature_mappings": [mapping.to_dict() for mapping in self.feature_mappings],
            "threshold": self.threshold,
            "legacy_config_path": self.legacy_config_path,
            "legacy_svm_path": self.legacy_svm_path,
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> ModelMetadata:
        """Builds model metadata from a dict payload."""
        schema_version = str(value.get("schema_version", SCHEMA_VERSION))
        if schema_version != SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported model schema version '{schema_version}'; expected '{SCHEMA_VERSION}'."
            )

        mappings = tuple(
            FeatureMapping.from_dict(item) for item in value.get("feature_mappings", [])
        )
        return cls(
            classifier_task=ClassifierTask(value["classifier_task"]),
            artifact_format=ArtifactFormat(value["artifact_format"]),
            artifact_path=str(value["artifact_path"]),
            sklearn_version=str(value["sklearn_version"]),
            feature_extractor=str(value["feature_extractor"]),
            feature_mappings=mappings,
            threshold=_as_optional_float(value.get("threshold")),
            legacy_config_path=_as_optional_str(value.get("legacy_config_path")),
            legacy_svm_path=_as_optional_str(value.get("legacy_svm_path")),
            schema_version=schema_version,
        )


def save_model_metadata(metadata: ModelMetadata, output_path: Path) -> None:
    """Saves model metadata as JSON."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(metadata.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_model_metadata(path: Path) -> ModelMetadata:
    """Loads model metadata from JSON."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected top-level object in {path}, got {type(payload).__name__}")
    return ModelMetadata.from_dict(payload)


def _as_optional_float(value: Any) -> float | None:
    """Converts a value to float or None."""
    if value is None:
        return None
    return float(value)


def _as_optional_str(value: Any) -> str | None:
    """Converts a value to str or None."""
    if value is None:
        return None
    return str(value)
