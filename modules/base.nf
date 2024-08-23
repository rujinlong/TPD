process extract_seqs_from_gbk {
    label "tpd_base"
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p01_extract_seqs"
    publishDir "$params.report/$sampleID", pattern: "${sampleID}_*"

    input:
    tuple val(sampleID), path(gbk)

    output:
    tuple val(sampleID), path("${sampleID}_genome.fasta"), emit: genome_ch
    tuple val(sampleID), path("${sampleID}_gene.fna"), emit: gene_ch
    tuple val(sampleID), path("${sampleID}_protein.faa"), emit: protein_ch

    when:
    task.ext.when == null || task.ext.when

    """
    extract_CDS_from_gbk.py -i $gbk -o $sampleID
    """
}