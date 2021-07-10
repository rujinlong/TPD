process DFAST {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p02_dfast"
    publishDir "$params.report/$sampleID", pattern: "statistics.txt"

    input:
    tuple val(sampleID), path(genome)

    output:
    path("dfast_output/*")
    tuple val(sampleID), path("dfast_output/genome.gbk"), emit: draft_gbk
    tuple val(sampleID), path("dfast_output/genome.fna"), emit: draft_genome_fna
    tuple val(sampleID), path("dfast_output/cds.fna"), emit: cds_ch
    tuple val(sampleID), path("protein_LOCUS.faa"), emit: protein_locus
    path("statistics.txt")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    dfast --genome ${genome} -o dfast_output --use_original_name t --cpu $task.cpus --minimum_length 30
    sed 's/>.*|LOCUS_/>LOCUS_/' dfast_output/protein.faa > protein_LOCUS.faa
    ln -s dfast_output/statistics.txt .
    """
}