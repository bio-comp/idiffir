from __future__ import annotations

import json
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

from iDiffIR.SpliceGrapher.predict.model_schema import load_model_metadata
from iDiffIR.SpliceGrapher.predict.sklearn_migration import migrate_archive, migrate_archives

ROOT = Path(__file__).resolve().parents[1]


def _write_legacy_archive(
    archive_path: Path,
    *,
    with_fixture: bool,
) -> None:
    """Create a tiny synthetic legacy classifier bundle for migration tests."""
    cfg_text = "\n".join(
        [
            "[SpliceGrapherConfig]",
            "fasta_path = gt_don.fa",
            "svm_path = gt_don.svm",
            "dimer = gt",
            "c = 1.0",
            "threshold = 0.4",
            "mink = 1",
            "maxk = 1",
            "degree = 1",
            "gamma = -1.0",
            "acceptor = False",
            "normalize = True",
            "",
        ]
    )
    fasta_text = "\n".join(
        [
            ">site_1 label=1",
            "AAAAAA",
            ">site_2 label=1",
            "AAAATA",
            ">site_3 label=0",
            "CCCCCC",
            ">site_4 label=0",
            "CCCGCC",
            "",
        ]
    )
    fixture_payload = {
        "schema_version": "1",
        "sequences": ["AAAAAA", "AAAATA", "CCCCCC", "CCCGCC"],
        "legacy_predictions": [1, 1, 0, 0],
    }

    with zipfile.ZipFile(archive_path, "w") as bundle:
        bundle.writestr("gt_don.cfg", cfg_text)
        bundle.writestr("gt_don.fa", fasta_text)
        bundle.writestr("gt_don.svm", "legacy-pyml-bytes")
        if with_fixture:
            bundle.writestr("gt_don.legacy_fixture.json", json.dumps(fixture_payload))


def test_migrate_archive_writes_joblib_and_metadata(tmp_path: Path) -> None:
    archive_path = tmp_path / "Synthetic_species.zip"
    output_root = tmp_path / "out"
    _write_legacy_archive(archive_path, with_fixture=False)

    report = migrate_archive(archive_path=archive_path, output_root=output_root)

    assert report["species"] == "Synthetic_species"
    assert len(report["models"]) == 1
    model_report = report["models"][0]
    assert model_report["records"] == 4
    assert model_report["legacy_compare_status"] == "missing_fixture"

    model_path = Path(model_report["artifact_path"])
    metadata_path = Path(model_report["metadata_path"])
    assert model_path.exists()
    assert metadata_path.exists()

    metadata = load_model_metadata(metadata_path)
    assert metadata.artifact_path == model_path.name
    assert metadata.threshold == pytest.approx(0.4)


def test_migrate_archive_reports_legacy_fixture_agreement(tmp_path: Path) -> None:
    archive_path = tmp_path / "Fixture_species.zip"
    output_root = tmp_path / "out"
    _write_legacy_archive(archive_path, with_fixture=True)

    report = migrate_archive(archive_path=archive_path, output_root=output_root)

    model_report = report["models"][0]
    assert model_report["legacy_compare_status"] == "completed"
    assert model_report["legacy_fixture_agreement"] == pytest.approx(1.0)
    assert model_report["legacy_fixture_records"] == 4


def test_migrate_archives_writes_aggregate_report(tmp_path: Path) -> None:
    archive_path = tmp_path / "Aggregate_species.zip"
    output_root = tmp_path / "out"
    report_path = tmp_path / "reports" / "migration_report.json"
    _write_legacy_archive(archive_path, with_fixture=False)

    report = migrate_archives(
        archive_paths=[archive_path],
        output_root=output_root,
        report_path=report_path,
    )

    assert report["summary"]["archive_count"] == 1
    assert report["summary"]["model_count"] == 1
    assert report_path.exists()


def test_migration_cli_runs_and_writes_report(tmp_path: Path) -> None:
    archive_path = tmp_path / "Cli_species.zip"
    output_root = tmp_path / "out"
    report_path = tmp_path / "reports" / "migration_report.json"
    _write_legacy_archive(archive_path, with_fixture=True)

    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "migrate_classifier_archives.py"),
            "--output-root",
            str(output_root),
            "--report-path",
            str(report_path),
            "--archives",
            str(archive_path),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert report_path.exists()
