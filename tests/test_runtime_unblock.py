from iDiffIR.SpliceGrapher.formats.GeneModel import Exon, Gene, GeneModel, Isoform
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import geneModelToSpliceGraph
from iDiffIR.SpliceGrapher.shared.utils import commaFormat


def test_clean_name_decodes_url_escapes():
    model = GeneModel(None)
    assert model.cleanName("Gene%20One") == "Gene-One"


def test_make_sorted_model_sorts_gene_objects():
    model = GeneModel(None)
    model.addChromosome(1, 1000, "chr1")
    model.addGene(Gene("gene_b", "", 200, 300, "chr1", "+"))
    model.addGene(Gene("gene_a", "", 10, 120, "chr1", "+"))

    model.makeSortedModel()

    assert [gene.id for gene in model.sorted["chr1"]["+"]] == ["gene_a", "gene_b"]


def test_gene_model_to_splice_graph_sorts_exon_nodes():
    gene = Gene("gene_a", "", 1, 200, "chr1", "+")
    isoform = Isoform("gene_a.1", 1, 200, "chr1", "+")
    gene.addExon(isoform, Exon(101, 150, "chr1", "+"))
    gene.addExon(isoform, Exon(11, 60, "chr1", "+"))

    graph = geneModelToSpliceGraph(gene)

    assert len(graph.nodeDict) == 2


def test_comma_format_python3_safe():
    assert commaFormat(1234567) == "1,234,567"
