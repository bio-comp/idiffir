from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from iDiffIR.splice_core.filtering import (
    SequenceSiteLabel,
    SpliceSiteSvmModel,
    TrainableSequenceRecord,
    load_splice_site_model,
)


def _training_records() -> list[TrainableSequenceRecord]:
    """Build a tiny deterministic labeled dataset for classifier tests."""
    return [
        TrainableSequenceRecord(sequence="GTAAAAAA", label=SequenceSiteLabel.POSITIVE),
        TrainableSequenceRecord(sequence="GTCCCCCC", label=SequenceSiteLabel.POSITIVE),
        TrainableSequenceRecord(sequence="AGAAAAAA", label=SequenceSiteLabel.NEGATIVE),
        TrainableSequenceRecord(sequence="AGCCCCCC", label=SequenceSiteLabel.NEGATIVE),
    ]


def test_splice_site_svm_model_trains_predicts_and_roundtrips(tmp_path: Path) -> None:
    pytest.importorskip("skops.io")
    records = _training_records()
    model = SpliceSiteSvmModel.train(records)

    predictions = model.predict_labels([record.sequence for record in records])
    expected = np.asarray([record.label.value for record in records], dtype=np.int64)
    accuracy = float((predictions == expected).mean())
    assert accuracy >= 0.75

    artifact_path = tmp_path / "donor_filter.skops"
    metadata_path = tmp_path / "donor_filter.metadata.json"
    model.save(artifact_path=artifact_path, metadata_path=metadata_path)
    restored = SpliceSiteSvmModel.load(artifact_path=artifact_path, metadata_path=metadata_path)
    restored_from_spec = load_splice_site_model(metadata_path)

    restored_predictions = restored.predict_labels([record.sequence for record in records])
    from_spec_predictions = restored_from_spec.predict_labels(
        [record.sequence for record in records]
    )
    assert np.array_equal(predictions, restored_predictions)
    assert np.array_equal(predictions, from_spec_predictions)


def test_training_cli_writes_skops_bundle(tmp_path: Path) -> None:
    pytest.importorskip("skops.io")
    dataset_path = tmp_path / "training.tsv"
    artifact_path = tmp_path / "donor_model.skops"
    metadata_path = tmp_path / "donor_model.metadata.json"
    dataset_path.write_text(
        "\n".join(
            [
                "sequence\tlabel",
                "GTAAAAAA\t1",
                "GTCCCCCC\t1",
                "AGAAAAAA\t0",
                "AGCCCCCC\t0",
                "",
            ]
        ),
        encoding="utf-8",
    )

    completed = subprocess.run(
        [
            sys.executable,
            "scripts/train_splice_site_svm.py",
            "--dataset",
            str(dataset_path),
            "--artifact",
            str(artifact_path),
            "--metadata",
            str(metadata_path),
            "--site-type",
            "donor",
            "--dimer",
            "GT",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert artifact_path.exists()
    assert metadata_path.exists()
