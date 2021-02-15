#!/usr/bin/env nextflow

params.cds = "cds.fasta"
params.pep = "pep.fasta"

cds = file(params.cds)
pep = file(params.pep)

ids = Channel.fromPath('*.id')

process codonAln {
    input:
    file geneID from ids
    output:

    """
    seqtk subseq $cds $geneID > gene.cds
    seqtk subseq $pep $geneID > gene.pep
    mafft --anysymbol gene.pep > pep.aln
    pal2nal.pl pep.aln gene.cds -output fasta  -nogap -nomismatch >gene.p2n

    """
}
