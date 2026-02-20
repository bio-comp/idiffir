from pathlib import Path

import pytest

from iDiffIR.SpliceGrapher.predict.model_schema import (
    ArtifactFormat,
    ClassifierTask,
    FeatureMapping,
    ModelMetadata,
    load_model_metadata,
    save_model_metadata,
)


def test_model_metadata_json_round_trip(tmp_path: Path) -> None:
    metadata = ModelMetadata(
        classifier_task=ClassifierTask.SITE_DONOR,
        artifact_format=ArtifactFormat.JOBLIB,
        artifact_path="models/gt_don.joblib",
        sklearn_version="1.5.2",
        feature_extractor="kmer_window_v1",
        feature_mappings=(
            FeatureMapping(legacy_name="mink", sklearn_name="kmer_min"),
            FeatureMapping(legacy_name="maxk", sklearn_name="kmer_max"),
        ),
        threshold=0.47,
        legacy_config_path="gt_don.cfg",
        legacy_svm_path="gt_don.svm",
    )

    metadata_path = tmp_path / "migrated" / "gt_don.metadata.json"
    save_model_metadata(metadata, metadata_path)

    loaded = load_model_metadata(metadata_path)
    assert loaded == metadata


def test_model_metadata_rejects_unsupported_schema_version() -> None:
    payload = {
        "schema_version": "99",
        "classifier_task": "gap",
        "artifact_format": "joblib",
        "artifact_path": "models/gap.joblib",
        "sklearn_version": "1.5.2",
        "feature_extractor": "gap_features_v1",
        "feature_mappings": [],
    }

    with pytest.raises(ValueError, match="Unsupported model schema version"):
        ModelMetadata.from_dict(payload)


def test_model_metadata_rejects_unknown_classifier_task() -> None:
    payload = {
        "schema_version": "1",
        "classifier_task": "unknown_task",
        "artifact_format": "joblib",
        "artifact_path": "models/model.joblib",
        "sklearn_version": "1.5.2",
        "feature_extractor": "features_v1",
        "feature_mappings": [],
    }

    with pytest.raises(ValueError, match="unknown_task"):
        ModelMetadata.from_dict(payload)
