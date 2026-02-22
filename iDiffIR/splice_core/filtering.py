"""Modern sklearn/skops splice-site filtering primitives."""

from __future__ import annotations

import json
from dataclasses import dataclass
from enum import Enum, IntEnum
from pathlib import Path
from typing import Protocol, cast

import numpy as np
import pysam


class SpliceSiteType(str, Enum):
    """Splice-site classes consumed by context extraction and classification."""

    DONOR = "donor"
    ACCEPTOR = "acceptor"


class SequenceSiteLabel(IntEnum):
    """Binary labels used for training splice-site SVM classifiers."""

    NEGATIVE = 0
    POSITIVE = 1


@dataclass(frozen=True)
class TrainableSequenceRecord:
    """One labeled sequence row for SVM training."""

    sequence: str
    label: SequenceSiteLabel


@dataclass(frozen=True)
class SequenceSiteContext:
    """Extracted sequence context around one genomic splice-site position."""

    chrom: str
    position: int
    strand: str
    site_type: SpliceSiteType
    exon_sequence: str
    intron_sequence: str
    model_sequence: str
    dimer: str


@dataclass(frozen=True)
class SpliceSiteModelMetadata:
    """Metadata persisted with secure sklearn model artifacts."""

    site_type: SpliceSiteType
    dimer: str
    exon_window: int
    intron_window: int
    threshold: float
    ngram_min: int
    ngram_max: int
    schema_version: str = "1"

    def to_dict(self) -> dict[str, str | int | float]:
        """Serialize metadata to JSON-compatible dict."""
        return {
            "schema_version": self.schema_version,
            "site_type": self.site_type.value,
            "dimer": self.dimer,
            "exon_window": self.exon_window,
            "intron_window": self.intron_window,
            "threshold": self.threshold,
            "ngram_min": self.ngram_min,
            "ngram_max": self.ngram_max,
        }

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> SpliceSiteModelMetadata:
        """Parse metadata from a JSON dict payload."""
        return cls(
            site_type=SpliceSiteType(str(payload["site_type"])),
            dimer=str(payload["dimer"]),
            exon_window=int(payload["exon_window"]),
            intron_window=int(payload["intron_window"]),
            threshold=float(payload["threshold"]),
            ngram_min=int(payload["ngram_min"]),
            ngram_max=int(payload["ngram_max"]),
            schema_version=str(payload.get("schema_version", "1")),
        )


class _SequenceModel(Protocol):
    """Minimal sklearn-like contract needed by this module."""

    def fit(self, sequences: list[str], labels: np.ndarray) -> _SequenceModel:
        """Train the model on sequence strings."""

    def predict(self, sequences: list[str]) -> np.ndarray:
        """Predict class labels for sequence strings."""

    def decision_function(self, sequences: list[str]) -> np.ndarray:
        """Return decision-function scores for sequence strings."""


class SpliceSiteSvmModel:
    """Securely persisted char-kmer SVM classifier for splice-site filtering."""

    def __init__(self, *, pipeline: _SequenceModel, metadata: SpliceSiteModelMetadata) -> None:
        """Store trained pipeline and model metadata."""
        self._pipeline = pipeline
        self.metadata = metadata

    @classmethod
    def train(
        cls,
        records: list[TrainableSequenceRecord],
        *,
        site_type: SpliceSiteType = SpliceSiteType.DONOR,
        dimer: str = "GT",
        exon_window: int = 8,
        intron_window: int = 15,
        threshold: float = 0.0,
        ngram_min: int = 1,
        ngram_max: int = 3,
        c_value: float = 1.0,
    ) -> SpliceSiteSvmModel:
        """Train one SVC model with explicit StandardScaler normalization."""
        if not records:
            raise ValueError("At least one training record is required.")

        from sklearn.feature_extraction.text import CountVectorizer
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC

        sequences = [record.sequence.upper() for record in records]
        labels = np.asarray([int(record.label) for record in records], dtype=np.int64)
        pipeline = cast(
            _SequenceModel,
            Pipeline(
                [
                    (
                        "vectorizer",
                        CountVectorizer(
                            analyzer="char",
                            ngram_range=(ngram_min, ngram_max),
                            lowercase=False,
                        ),
                    ),
                    ("scaler", StandardScaler(with_mean=False)),
                    ("svc", SVC(C=c_value, kernel="linear", class_weight="balanced")),
                ]
            ),
        )
        pipeline.fit(sequences, labels)
        metadata = SpliceSiteModelMetadata(
            site_type=site_type,
            dimer=dimer.upper(),
            exon_window=exon_window,
            intron_window=intron_window,
            threshold=threshold,
            ngram_min=ngram_min,
            ngram_max=ngram_max,
        )
        return cls(pipeline=pipeline, metadata=metadata)

    def predict_labels(self, sequences: list[str]) -> np.ndarray:
        """Predict integer class labels for sequence strings."""
        normalized_sequences = [sequence.upper() for sequence in sequences]
        return np.asarray(self._pipeline.predict(normalized_sequences), dtype=np.int64)

    def predict_scores(self, sequences: list[str]) -> np.ndarray:
        """Predict decision-function scores for sequence strings."""
        scores = self._pipeline.decision_function([sequence.upper() for sequence in sequences])
        return np.atleast_1d(np.asarray(scores, dtype=np.float64))

    def classify_context(self, context: SequenceSiteContext) -> tuple[int, float]:
        """Classify one extracted sequence context and return `(label, score)`."""
        score = float(self.predict_scores([context.model_sequence])[0])
        label = int(score > self.metadata.threshold)
        return label, score

    def save(self, *, artifact_path: Path, metadata_path: Path) -> None:
        """Persist model pipeline with skops plus sidecar JSON metadata."""
        import skops.io as skops_io

        artifact_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        skops_io.dump(self._pipeline, str(artifact_path))
        metadata_path.write_text(
            json.dumps(self.metadata.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load(cls, *, artifact_path: Path, metadata_path: Path) -> SpliceSiteSvmModel:
        """Load a persisted model and metadata pair from disk."""
        import skops.io as skops_io

        payload = json.loads(metadata_path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError(f"Expected JSON object in {metadata_path}")
        metadata = SpliceSiteModelMetadata.from_dict(payload)
        untrusted_types = skops_io.get_untrusted_types(file=str(artifact_path))
        pipeline = cast(_SequenceModel, skops_io.load(str(artifact_path), trusted=untrusted_types))
        return cls(pipeline=pipeline, metadata=metadata)


def load_splice_site_model(spec_path: str | Path) -> SpliceSiteSvmModel:
    """Load one model bundle from either `.skops` or `.metadata.json` path."""
    artifact_path, metadata_path = _resolve_bundle_paths(Path(spec_path))
    return SpliceSiteSvmModel.load(artifact_path=artifact_path, metadata_path=metadata_path)


def extract_sequence_context(
    *,
    fasta_path: str | Path,
    chrom: str,
    position: int,
    strand: str,
    site_type: SpliceSiteType,
    exon_window: int,
    intron_window: int,
) -> SequenceSiteContext:
    """Extract exon/intron sequence context with deterministic orientation handling."""
    if strand not in {"+", "-"}:
        raise ValueError("strand must be '+' or '-'")
    if position < 1:
        raise ValueError("position must be a 1-based genomic coordinate")
    if exon_window < 1 or intron_window < 2:
        raise ValueError("exon_window must be >=1 and intron_window must be >=2")

    with pysam.FastaFile(str(fasta_path)) as reference:
        exon_start, exon_end, intron_start, intron_end = _window_coordinates(
            position=position,
            strand=strand,
            site_type=site_type,
            exon_window=exon_window,
            intron_window=intron_window,
        )
        exon_sequence = _fetch_one_based(
            reference=reference,
            chrom=chrom,
            start=exon_start,
            end=exon_end,
        )
        intron_sequence = _fetch_one_based(
            reference=reference,
            chrom=chrom,
            start=intron_start,
            end=intron_end,
        )

    if strand == "-":
        exon_sequence = _reverse_complement(exon_sequence)
        intron_sequence = _reverse_complement(intron_sequence)

    if site_type == SpliceSiteType.DONOR:
        model_sequence = exon_sequence + intron_sequence
        dimer = intron_sequence[:2]
    else:
        model_sequence = intron_sequence + exon_sequence
        dimer = intron_sequence[-2:]

    return SequenceSiteContext(
        chrom=chrom,
        position=position,
        strand=strand,
        site_type=site_type,
        exon_sequence=exon_sequence.upper(),
        intron_sequence=intron_sequence.upper(),
        model_sequence=model_sequence.upper(),
        dimer=dimer.upper(),
    )


def _window_coordinates(
    *,
    position: int,
    strand: str,
    site_type: SpliceSiteType,
    exon_window: int,
    intron_window: int,
) -> tuple[int, int, int, int]:
    """Compute exon/intron windows around one splice-site position."""
    if site_type == SpliceSiteType.DONOR and strand == "+":
        return position - exon_window + 1, position, position + 1, position + intron_window
    if site_type == SpliceSiteType.DONOR and strand == "-":
        return position, position + exon_window - 1, position - intron_window, position - 1
    if site_type == SpliceSiteType.ACCEPTOR and strand == "+":
        return position, position + exon_window - 1, position - intron_window, position - 1
    return position - exon_window + 1, position, position + 1, position + intron_window


def _fetch_one_based(*, reference: pysam.FastaFile, chrom: str, start: int, end: int) -> str:
    """Fetch an inclusive 1-based reference interval from FASTA."""
    if start < 1 or end < start:
        raise ValueError(f"Invalid reference window [{start}, {end}]")
    return reference.fetch(chrom, start - 1, end)


def _reverse_complement(sequence: str) -> str:
    """Return DNA reverse complement for uppercase/lowercase ASCII sequence."""
    translation = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(translation)[::-1]


def _resolve_bundle_paths(spec_path: Path) -> tuple[Path, Path]:
    """Resolve artifact and metadata paths from one user-provided model spec path."""
    if spec_path.suffix == ".skops":
        artifact_path = spec_path
        metadata_path = spec_path.with_suffix(".metadata.json")
        return artifact_path, metadata_path
    if spec_path.name.endswith(".metadata.json"):
        artifact_name = spec_path.name[: -len(".metadata.json")] + ".skops"
        artifact_path = spec_path.with_name(artifact_name)
        metadata_path = spec_path
        return artifact_path, metadata_path
    if spec_path.suffix == ".json":
        artifact_path = spec_path.with_suffix(".skops")
        metadata_path = spec_path
        return artifact_path, metadata_path
    raise ValueError(
        f"Unsupported model specification path '{spec_path}'. "
        "Use either <name>.skops or <name>.metadata.json."
    )
