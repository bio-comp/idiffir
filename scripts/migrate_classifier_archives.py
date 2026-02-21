"""CLI for migrating legacy SpliceGrapher classifier archives to sklearn artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

from iDiffIR.SpliceGrapher.predict.sklearn_migration import migrate_archives


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for classifier archive migration."""
    parser = argparse.ArgumentParser(
        description=(
            "Migrate legacy classifier zip archives (.cfg/.fa/.svm) to "
            "sklearn .joblib artifacts and metadata reports."
        )
    )
    parser.add_argument(
        "--archives",
        nargs="+",
        type=Path,
        default=[],
        help="Explicit archive paths to migrate.",
    )
    parser.add_argument(
        "--archive-dir",
        type=Path,
        default=None,
        help="Optional directory containing classifier *.zip archives.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Output directory where migrated artifacts are written.",
    )
    parser.add_argument(
        "--report-path",
        type=Path,
        required=True,
        help="Path where aggregate migration report JSON is written.",
    )
    return parser.parse_args(argv)


def _collect_archive_paths(args: argparse.Namespace) -> list[Path]:
    """Collect unique archive paths from CLI flags."""
    archive_paths = [path.resolve() for path in args.archives]
    if args.archive_dir:
        archive_paths.extend(sorted(path.resolve() for path in args.archive_dir.glob("*.zip")))

    unique_paths = sorted(set(archive_paths))
    return unique_paths


def main(argv: list[str] | None = None) -> int:
    """Run classifier archive migration and write aggregate report."""
    args = _parse_args(argv)
    archive_paths = _collect_archive_paths(args)
    if not archive_paths:
        raise ValueError("No archives provided. Use --archives and/or --archive-dir.")

    report = migrate_archives(
        archive_paths=archive_paths,
        output_root=args.output_root,
        report_path=args.report_path,
    )
    summary = report["summary"]
    print(
        f"Migrated {summary['archive_count']} archive(s), "
        f"{summary['model_count']} model(s) -> {args.report_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
