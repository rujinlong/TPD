#!/usr/bin/env nextflow
// Usage: nextflow run predict_prophages.nf -profile hpc_slurm --datadir "$(pwd)/genomes" -resume

nextflow.enable.dsl=2

log.info """\
NF - TPD PIPELINE
=========================
result: ${params.outdir}
report: ${params.report}
"""

include { DFAST } from './modules/anno_dfast.nf'
include { predict_prophage_phispy } from './modules/prophage_phispy'
include { predict_prophage_phigaro } from './modules/prophage_phigaro'
include { predict_prophage_phageboost } from './modules/prophage_phageboost'
include { update_gbk_prophage } from './modules/update_prophage'
include { extract_prophages } from './modules/extract_prophage'


workflow {
    genomes = "${params.datadir}/*.fna"
    genomes_ch = channel.fromPath(genomes).map { [it.toString().split("/")[-1].replaceAll(/.fna$/, ""), it] }.unique()
    println "genomes: ${genomes}"

    DFAST(genomes_ch)
    predict_prophage_phispy(DFAST.out.draft_gbk)
    predict_prophage_phigaro(DFAST.out.draft_genome_fna)
    predict_prophage_phageboost(DFAST.out.draft_genome_fna)
    update_gbk_prophage(DFAST.out.draft_gbk.join(predict_prophage_phispy.out.phispy_tsv, by:0).join(predict_prophage_phigaro.out.phigaro_tsv, by:0).join(predict_prophage_phageboost.out.phageboost_tsv, by:0))
    extract_prophages(update_gbk_prophage.out.prophages_ch)
}



