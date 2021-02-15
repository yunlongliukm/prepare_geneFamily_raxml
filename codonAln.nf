#!/usr/bin/env nextflow

ids = Channel.fromPath('*.id')

process codonAln {
    input:
    file geneID from ids
    output:

    """
    seqtk subseq cds.fasta geneID > gene.cds
    seqtk subseq pep.fasta geneID > gene.pep
    mafft --anysymbol gene.pep > pep.aln
    pal2nal.pl pep.aln gene.cds -output fasta  -nogap -nomismatch >gene.p2n

    """
}
