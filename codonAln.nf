#!/usr/bin/env nextflow

params.cds = "cds.fasta"
params.pep = "pep.fasta"

cds = file(params.cds)
pep = file(params.pep)

ids = Channel.fromPath('gene*.id')

/*
process faIDX {
    output:
    ${cds}.fai ${pep}.fai
    """
    [ ! -f ${cds}.fai ] && samtools faidx $cds
    [ ! -f ${pep}.fai ] && samtools faidx $pep
    """
}
*/
process codonAln {
    input:
    file geneID from ids
    output:
    file "${ID}.p2n"  into pal2nal_file
    script:
    ID= geneID.getSimpleName()
    """
    samtools faidx $cds \$(cat $geneID) >gene.cds
    samtools faidx $pep \$(cat $geneID) >gene.pep
    mafft --anysymbol gene.pep > pep.aln
    pal2nal.pl pep.aln gene.cds -output fasta  -nogap -nomismatch >gene.p2n
    seqtk rename gene.p2n ${ID}_ >${ID}.p2n
    """
}

process mergePal2nal {
    input:
    file p2n from pal2nal_file
    output:
    file "concat.fa" into concatAln
    script:
    p2nComb = p2n.collect{ "-V $it" }.join(' ')
    """
    printf $p2nComb
    AMAS.py concat -i $p2nComb -f fasta -d dna -t concat.fa
    """
}
