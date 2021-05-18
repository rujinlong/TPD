process predict_CRISPR {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p03_CRISPR"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(draft_genome_fna)

    output:
    path("${sampleID}_CRISPR_*")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    crt crt $draft_genome_fna ${sampleID}_CRISPR_CRT.txt
    cctyper $draft_genome_fna ${sampleID}_CRISPR_cctyper
    """
}