"""Modern splice-core APIs shared across iDiffIR and TAPIS."""

from iDiffIR.splice_core.events import (
    AlternativeSpliceSiteEvent,
    AlternativeSpliceSiteType,
    ExonSkippingEvent,
    IntronRetentionEvent,
    detect_alternative_splice_sites,
    detect_exon_skipping,
    detect_intron_retention,
)
from iDiffIR.splice_core.filtering import (
    SequenceSiteContext,
    SequenceSiteLabel,
    SpliceSiteSvmModel,
    SpliceSiteType,
    TrainableSequenceRecord,
    extract_sequence_context,
    load_splice_site_model,
)
from iDiffIR.splice_core.graph import IntronSupport, SpliceGraphContext, build_splice_graph

__all__ = [
    "AlternativeSpliceSiteEvent",
    "AlternativeSpliceSiteType",
    "ExonSkippingEvent",
    "IntronRetentionEvent",
    "IntronSupport",
    "SequenceSiteContext",
    "SequenceSiteLabel",
    "SpliceGraphContext",
    "SpliceSiteSvmModel",
    "SpliceSiteType",
    "TrainableSequenceRecord",
    "build_splice_graph",
    "detect_alternative_splice_sites",
    "detect_exon_skipping",
    "detect_intron_retention",
    "extract_sequence_context",
    "load_splice_site_model",
]
