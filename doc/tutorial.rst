========
Tutorial
========
This tutorial is meant as a complete walk-through for identifying
differential IR from RNA-seq data.  The process is as follows:

1. Align reads
2. Filter alignments
3. Get read depths
4. Analyze depths for differential IR

Align reads
-----------
**iDiffIR** accepts any sorted, indexed BAM file for single- or 
paired-end reads.  It may be helpful to use an aligner that
provides an "NH" tag (such as STAR_ or TopHat_)
to identify unique alignments; however, there are other ways of
removing multi-reads (reads mapping to more than 
one genomic location).

Filtering alignments
--------------------
It is good practice to remove multi-reads and filter spliced 
alignments for false-positive junctions  before beginning 
transcriptome analysis of RNA-Seq.  Assuming your alignments 
are in a file named aligned.bam, run the following to extract 
only uniquely alinged reads:

.. code-block:: none
   
   $samtools view -h aligned.bam | grep -we 'NH:i:1' -e '^@*' > unique.sam
   
You may wish to filter out potential false-positive splice junctions 
using SpliceGrapher\ [SG]_ :program:`sam_filter.py` 
from uniquely aligned reads.  While the presence
of false-positive splice junctions may not adversely affect
identifying differential IR analysis, it is a critical step for
detecting all other forms of AS (exon skipping, alternative 
acceptor/donor events).

.. program-output:: sam_filter.py --help



Getting read depths
-------------------

Creating indexed, sorted BAM files
..................................
If your reads are in the SAM format, you can convert it to a indexed,
sorted BAM file using :program:`convertSam.py`.  

.. program-output:: convertSam.py --help

If your alignments are in a sorted BAM file, an index can be 
generated using:

.. code-block:: none

   $samtools index <BAMFILE>

Both methods will produce .bai file, indicating that the BAM file 
has been indexed.


Getting chromosomal read depths
...............................
Now we are prepared to generate chromsomal read depths from the
alignments using :program:`getDepths.py`.

.. program-output:: getDepths.py --help

.. note::

   It is *strongly* recommended to use a clear naming system for 
   the output directories (**-o** option) generated from 
   :program:`getDepths.py`.  For instance if you have RNA-seq 
   between an wildtype and mutant with 2 replicates, 
   you may want to use a naming convention such as:
   
   .. code-block:: none

      $getDepths.py -o wt_rep1_depths wt_rep1/filtered.bam
      $getDepths.py -o wt_rep2_depths wt_rep2/filtered.bam
      $getDepths.py -o mt_rep1_depths mt_rep1/filtered.bam
      $getDepths.py -o mt_rep2_depths mt_rep2/filtered.bam

Running **iDiffIR**
-------------------

.. command-output:: idiffir.py --help 

.. _STAR: https://code.google.com/p/rna-star/
.. _TopHat: http://ccb.jhu.edu/software/tophat/index.shtml

.. [SG] Rogers, MF, Thomas, J, Reddy, AS, Ben-Hur, A (2012). 
	SpliceGrapher: detecting patterns of alternative splicing 
	from RNA-Seq data in the context of gene models and 
	EST data. *Genome Biol*., 13, 1:R4.
