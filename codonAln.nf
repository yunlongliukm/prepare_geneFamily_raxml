
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
    file "${ID}.p2n.fa"  into pal2nal_file
    script:
    ID= geneID.getSimpleName()
    """
    faops order $cds $geneID  gene.cds
    faops order $pep $geneID  gene.pep
    mafft --anysymbol gene.pep > pep.aln
    pal2nal.pl pep.aln gene.cds -output fasta  -nogap -nomismatch >gene.p2n
    cat $geneID |awk '{print \$1"\tseq"NR}' >switch.id
    faops replace -s gene.p2n switch.id ${ID}.p2n.fa
    """
}

process mergePal2nal {
    input:
    file p2n from pal2nal_file
    output:
    file "concat.fa" into concatAln
    script:
    p2nComb = p2n.collect{ "$it" }.join(' ')

    """
    AMAS.py concat -i $p2nComb -f fasta -d dna -t concat.fa
    """
}
