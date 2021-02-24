#!/usr/bin/env nextflow

params.cds = "/ds3200_1/users_root/mapengfei/10bg/perfect_gene_seq/cds.fasta"
params.pep = "/ds3200_1/users_root/mapengfei/10bg/perfect_gene_seq/pep.fasta"

cds = file(params.cds)
pep = file(params.pep)


ids = Channel.fromPath('gene*.id')

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
    cat $geneID |sed -e 's/\\(LOC_Os[0-9]*\\)g.*/\\1/' -e 's/\\(Chr.*_Ol.*\\)g.*/\\1/' -e 's/\\(LG.*_[a-z]*\\).*/\\1/' -e 's/evm/rhi/' >id.short
    paste $geneID id.short >switch.id
    python -m jcvi.formats.fasta format gene.p2n ${ID}.p2n --switch=switch.id
    """
}

process mergePal2nal {
    input:
    file p2n from pal2nal_file.filter{ it.size()>0 }.collect()
    output:
    file "concat.fa" into concatAln
    """
    AMAS.py concat -i $p2n -f fasta -d dna -t concat.fa
    """
}



process raxmlBlock {
    cpus 10
    executor 'lsf'
    queue 'Q88C6T_X1'

    publishDir './raxml', mode: 'copy', overwrite: false

    input:
    file concat from concatAln
    output:
    file "RAxML*" into raxmlTree
    """
    raxmlHPC-PTHREADS-AVX2 -f a -T 10 -m GTRGAMMA -n tre -s concat.fa -p 54321 -N 200 -x 12345
    """
    

}

