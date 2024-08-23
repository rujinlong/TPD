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


process run_RagTag {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p01_RagTag"
    publishDir "$params.report/$sampleID", pattern: "*.stats"

    input:
    tuple val(sampleID), path(ref_fna)

    output:
    tuple val(sampleID), path('ragtag_output/*') 
    tuple val(sampleID), path("${sampleID}_1k.fna"), emit: genome_ch
    tuple val(sampleID), path('*.stats') 

    when:
    params.mode == "all"|| params.mode == 'genome'

    """
    rawseqs="${baseDir}/data/${sampleID}_scaffolds.fna"
    ragtag.py scaffold $ref_fna \$rawseqs -t $task.cpus -o ragtag_output
    ln -s ragtag_output/*.stats .
    sed "s/^>/>${sampleID}_/" ragtag_output/ragtag.scaffolds.fasta | sed 's/ /_/g' > ${sampleID}.fna
    seqkit seq -m 1000 ${sampleID}.fna > ${sampleID}_1k.fna
    """
}