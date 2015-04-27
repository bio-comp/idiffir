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
only uniquely aligned reads:

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


Running **iDiffIR**
-------------------

.. command-output:: idiffir.py --help 

While **iDiffIR** offers many options, default values should work
for most cases.  Using the above naming convention, a way to run
**iDiffIR** would be:

.. code-block:: none

   $idiffir.py -l Mutant Wildtype -o my_result my_genome.gff mt_rep1.bam:mt_rep2.bam:mt_rep3.bam wt_rep1.bam:wt_rep2.bam:wt_rep3.bam

The **-l** option provides titles for the reads and will be displayed
on the splice graphs of genes containing significant differential IR
events.    You may wish to use the **-d** option to lower the
gene expression fold-change threshold. Introns in genes 
with a gene expression fold-change above this 
threshold will not be tested.  While **iDiffIR** does 
adjust for differential gene expression, in some extreme cases it may 
be necessary.

Interpreting iDiffIR output
---------------------------
**iDiffIR** builds an output directory as follows:

::

   $ tree my_result
     my_result
     |-- figures
     |   |-- gene1.pdf
     |   |-- gene2.pdf
     |   |-- mva.pdf
     |   |-- pvalues.pdf
     |-- figuresLog
     |   |-- gene1.pdf
     |   |-- gene2.pdf
     |-- lists
         |-- allDIRGenes.txt
	 |-- allIntrons.txt

     3 directories, 8 files

Two figures directories **figures** and **figuresLog** are created
and contain splice graphs figures of genes with at least 
one significant
differential IR event.  The **figuresLog** directory contains
graphs in :math:`\log` scale which can be helpful in cases where
intronic expression is orders of magnitude less than the exonic 
expression within the gene.  The **lists** directory contain
two important files.  **allDIRGenes.txt** contain the
gene IDs for all genes containing a significant differential IR
event.  This is should allow convenient downstream analyses, 
such as GO term enrichment.  Finally, **allIntrons.txt** is a
tab-delimited file containing intronic coordinates and statistics
for all tested introns.  The fields of each row are:

1. **geneID** 
     the gene's identifier

2. **lowExonCoords** 
     the coordinates of the lowest (W.R.T. genomic 
     position) flanking exon

3. **intronCoords** 
     the coordinates of the intronic region tested

4. **highExonCoords** 
     the coordinates of the highest (W.R.T. genomic 
     position) flanking exon

5. **pValue**
     the :math:`p`\ -value of the tested intron (using a 2-sided
     :math:`Z`\ -score test).

6. **adjPValue**
     the multiple testing adjusted `p`\ -value

7. **logFoldChange**
     the :math:`\log`\ -fold change of the tested intron (w.r.t. the
     first-given condition)

8. **intronExp**
      the expression of the intron, computed as 
      :math:`\displaystyle\frac{1}{2}\log\left( x_1 + x_2\right)`,
      where :math:`x_1, x_2` are the average read depth of the intron
      in condition 1 and condition 2, respectively.
	    
9. **statistic**
      the test statistic (before z-score conversion)

10. **bestA**
      the pseudo-count value (:math:`a`) that minimizes the 
      :math:`p`\ -value

11. **known**
      whether this intron is known as a retained intron


.. todo::

   Add differential exon skipping

.. todo::

   Add MISO script usage for testing Alt 5', 3'

.. todo::

   Add simulation
.. _STAR: https://code.google.com/p/rna-star/
.. _TopHat: http://ccb.jhu.edu/software/tophat/index.shtml

.. [SG] Rogers, MF, Thomas, J, Reddy, AS, Ben-Hur, A (2012). 
	SpliceGrapher: detecting patterns of alternative splicing 
	from RNA-Seq data in the context of gene models and 
	EST data. *Genome Biol*., 13, 1:R4.
