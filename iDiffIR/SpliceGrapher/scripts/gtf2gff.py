#! /usr/bin/env python
"""Convert an ENSEMBL-style GTF annotation file into GFF3 output."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import tempfile
import sys

from iDiffIR.SpliceGrapher.formats.loader import loadGeneModels
from iDiffIR.SpliceGrapher.shared.utils import ezopen


USAGE = """%(prog)s GTF-file [options]

Converts an ENSEMBL GTF gene model annotation file into its GFF3 equivalent.
Note that it only accepts records with the protein_coding tag by default."""

ALL_ENSEMBL_SOURCES = [
    "3prime_overlapping_ncrna",
    "ambiguous_orf",
    "antisense",
    "disrupted_domain",
    "IG_C_gene",
    "IG_C_pseudogene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_J_pseudogene",
    "IG_M_gene",
    "IG_V_gene",
    "IG_V_pseudogene",
    "IG_Z_gene",
    "lincRNA",
    "miRNA",
    "miRNA_pseudogene",
    "misc_RNA",
    "misc_RNA_pseudogene",
    "Mt_rRNA",
    "Mt_tRNA",
    "Mt_tRNA_pseudogene",
    "ncRNA",
    "ncrna_host",
    "non_coding",
    "nonsense_mediated_decay",
    "non_stop_decay",
    "polymorphic_pseudogene",
    "processed_pseudogene",
    "processed_transcript",
    "protein_coding",
    "pseudogene",
    "retained_intron",
    "retrotransposed",
    "rRNA",
    "rRNA_pseudogene",
    "scRNA_pseudogene",
    "sense_intronic",
    "sense_overlapping",
    "snlRNA",
    "snoRNA",
    "snoRNA_pseudogene",
    "snRNA",
    "snRNA_pseudogene",
    "TEC",
    "transcribed_processed_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "TR_C_gene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_J_pseudogene",
    "tRNA",
    "tRNA_pseudogene",
    "TR_V_gene",
    "TR_V_pseudogene",
    "unitary_pseudogene",
    "unprocessed_pseudogene",
]
PROTEIN_CODING = "protein_coding"


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for the conversion script."""
    parser = argparse.ArgumentParser(usage=USAGE)
    parser.add_argument("gtf_file", help="Input GTF file")
    parser.add_argument(
        "-A",
        dest="alltypes",
        default=False,
        help="Accept all source types (overrides -E and -S) [default: %(default)s]",
        action="store_true",
    )
    parser.add_argument("-o", dest="output", default=None, help="Output file [default: stdout]")
    parser.add_argument(
        "-E",
        dest="ensembl",
        default=False,
        help="Accept all ENSEMBL source types (see --show-types) [default: %(default)s]",
        action="store_true",
    )
    parser.add_argument(
        "-S",
        dest="sources",
        default=PROTEIN_CODING,
        help="Comma-separated list of GTF source types to accept [default: %(default)s]",
    )
    parser.add_argument("-v", dest="verbose", default=False, help="Verbose mode [default: %(default)s]", action="store_true")
    parser.add_argument("-z", dest="gzip_output", default=False, help="Use gzip compression on output [default: %(default)s]", action="store_true")
    parser.add_argument("--show-types", dest="showtypes", default=False, help="Outputs ENSEMBLE source types and exits. [default: %(default)s]", action="store_true")
    return parser


def _filtered_gtf(input_path: Path, allowed_sources: set[str]) -> Path:
    """Write a temporary GTF filtered to source types listed in ``allowed_sources``."""
    tmp = tempfile.NamedTemporaryFile(mode="w", prefix="idiffir_gtf_filter_", suffix=".gtf", delete=False)
    tmp_path = Path(tmp.name)
    with tmp:
        for line in ezopen(str(input_path)):
            if not line or line.startswith("#"):
                tmp.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[1] in allowed_sources:
                tmp.write(line)
    return tmp_path


def _get_output_stream(output: str | None, gzip_output: bool):
    """Return an output stream for stdout/plain text/gzip output paths."""
    if output is None:
        return sys.stdout
    if gzip_output:
        return gzip.open(output, "wt", encoding="utf-8")
    return open(output, "w", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    """Convert a GTF file to GFF3 output using the shared loader backend."""
    parser = build_parser()
    opts = parser.parse_args(argv)

    if opts.showtypes:
        print("Known ENSEMBL source types:")
        for source_type in ALL_ENSEMBL_SOURCES:
            print("  ", source_type)
        return 0

    gtf_path = Path(opts.gtf_file)
    if not gtf_path.exists():
        parser.error(f"Input GTF file not found: {gtf_path}")

    temp_gtf_path: Path | None = None
    model_input = gtf_path
    if not opts.alltypes:
        known_sources = opts.sources.split(",") if not opts.ensembl else ALL_ENSEMBL_SOURCES
        model_input = _filtered_gtf(gtf_path, set(known_sources))
        temp_gtf_path = model_input

    try:
        gene_model = loadGeneModels(str(model_input), verbose=opts.verbose, alltypes=True)
        out_stream = _get_output_stream(opts.output, opts.gzip_output)
        try:
            gene_model.writeGFF(out_stream, verbose=opts.verbose)
        finally:
            if out_stream is not sys.stdout:
                out_stream.close()
    finally:
        if temp_gtf_path and temp_gtf_path.exists():
            temp_gtf_path.unlink()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
