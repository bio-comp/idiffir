# SpliceGrapher (Legacy Submodule)

SpliceGrapher is the legacy graph/prediction subsystem bundled inside iDiffIR.
It predicts splice graphs from gene models and sequencing evidence.

Original author:
- Mark F. Rogers (`rogersma@cs.colostate.edu`)

## Status in iDiffIR

This directory is being modernized incrementally for Python 3 compatibility.
Current work keeps behavior parity first, then removes or replaces legacy dependencies.

## Supported Inputs (legacy design)

- RNA-seq alignments: SAM, BED, WIG
- EST alignments: PSL
- Transcript descriptions: GTF
- Gene models: GFF3

## Usage in this Repository

- Preferred root workflow is documented in top-level `README.md`.
- Most entry points are run through `uv`, for example:

```bash
uv run python iDiffIR/SpliceGrapher/scripts/predict_splicegraph.py --help
```

## Legacy Notes

- Historical references in older SpliceGrapher docs mention `setup.py` and PyML-era tooling.
- Current migration issues track modern replacements for packaging, parser internals, and classifier tooling.

## Citation

Rogers MF, Thomas J, Reddy AS, Ben-Hur A. SpliceGrapher: detecting patterns of
alternative splicing from RNA-Seq data in the context of gene models and EST data.
Genome Biol. 2012 Jan 31;13(1):R4. doi: 10.1186/gb-2012-13-1-r4.
PMID: 22293517; PMCID: PMC3334585.

## Documentation and Examples

- User guide artifact: `doc/userguide.pdf`
- Example data/scripts: `iDiffIR/SpliceGrapher/examples`
