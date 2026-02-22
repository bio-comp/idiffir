"""Train a modern splice-site SVM and persist it with skops metadata."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from iDiffIR.splice_core.filtering import (
    SequenceSiteLabel,
    SpliceSiteSvmModel,
    SpliceSiteType,
    TrainableSequenceRecord,
)


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser for splice-site model training."""
    parser = argparse.ArgumentParser(
        description="Train a splice-site SVM from tabular sequence labels.",
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Path to TSV with columns: sequence,label",
    )
    parser.add_argument(
        "--artifact",
        required=True,
        help="Output .skops model path.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Output metadata JSON path.",
    )
    parser.add_argument(
        "--site-type",
        choices=[SpliceSiteType.DONOR.value, SpliceSiteType.ACCEPTOR.value],
        default=SpliceSiteType.DONOR.value,
        help="Classifier site type.",
    )
    parser.add_argument(
        "--dimer",
        default="GT",
        help="Expected splice-site dimer for this classifier.",
    )
    parser.add_argument(
        "--exon-window",
        type=int,
        default=8,
        help="Exonic context window length.",
    )
    parser.add_argument(
        "--intron-window",
        type=int,
        default=15,
        help="Intronic context window length.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Decision-function threshold for positive calls.",
    )
    parser.add_argument(
        "--ngram-min",
        type=int,
        default=1,
        help="Minimum character n-gram size.",
    )
    parser.add_argument(
        "--ngram-max",
        type=int,
        default=3,
        help="Maximum character n-gram size.",
    )
    parser.add_argument(
        "--c-value",
        type=float,
        default=1.0,
        help="SVC regularization parameter C.",
    )
    return parser


def load_records(dataset_path: Path) -> list[TrainableSequenceRecord]:
    """Load labeled training records from a tab-separated file."""
    records: list[TrainableSequenceRecord] = []
    with dataset_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence", "label"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError("Dataset must include tab-separated columns: sequence,label")
        for row in reader:
            sequence = str(row["sequence"]).strip().upper()
            label = int(str(row["label"]).strip())
            if label not in (0, 1):
                raise ValueError(f"Unsupported label value '{label}'; expected 0 or 1.")
            records.append(
                TrainableSequenceRecord(
                    sequence=sequence,
                    label=SequenceSiteLabel(label),
                )
            )
    if not records:
        raise ValueError("Dataset is empty.")
    return records


def main(argv: list[str] | None = None) -> int:
    """Train and persist one splice-site model bundle."""
    parser = build_parser()
    args = parser.parse_args(argv)

    dataset_path = Path(args.dataset)
    artifact_path = Path(args.artifact)
    metadata_path = Path(args.metadata)
    records = load_records(dataset_path)
    model = SpliceSiteSvmModel.train(
        records,
        site_type=SpliceSiteType(args.site_type),
        dimer=args.dimer,
        exon_window=args.exon_window,
        intron_window=args.intron_window,
        threshold=args.threshold,
        ngram_min=args.ngram_min,
        ngram_max=args.ngram_max,
        c_value=args.c_value,
    )
    model.save(artifact_path=artifact_path, metadata_path=metadata_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
