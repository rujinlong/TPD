#!/usr/bin/env nextflow
// Usage: nextflow run tpd.nf -profile hpc_slurm --mode "all" --datadir data -resume

nextflow.enable.dsl=2

log.info """\
NF - TPD PIPELINE
=========================
result: ${params.outdir}
report: ${params.report}
"""

include { extract_seqs_from_gbk } from './modules/extract_gbk_seqs'
include { viral_annotation_VOGDB } from './modules/anno_vogdb'
include { viral_annotation_pVOG } from './modules/anno_pvog'
include { viral_annotation_PFAM } from './modules/anno_pfam'
include { viral_annotation_KEGG } from './modules/anno_kegg'
include { cds_anno_ARG } from './modules/anno_arg'
include { predict_prophage_phispy } from './modules/prophage_phispy'
include { predict_prophage_phigaro } from './modules/prophage_phigaro'
include { update_gbk_cds_annotation } from './modules/update_anno_cds'
include { update_gbk_prophage } from './modules/update_prophage'
include { extract_prophages } from './modules/extract_prophage'
include { predict_CRISPR } from './modules/predict_crispr'


workflow {
    data = "${params.datadir}/*.gbk"
    samples = channel.fromPath(data).map { [it.toString().split("/")[-1].replaceAll(/.gbk$/, ""), it] }.unique()
    extract_seqs_from_gbk(samples)
    viral_annotation_VOGDB(extract_seqs_from_gbk.out.protein_ch)
    viral_annotation_pVOG(extract_seqs_from_gbk.out.protein_ch)
    viral_annotation_PFAM(extract_seqs_from_gbk.out.protein_ch)
    viral_annotation_KEGG(extract_seqs_from_gbk.out.protein_ch)
    cds_anno_ARG(extract_seqs_from_gbk.out.gene_ch)

    predict_prophage_phispy(samples)
    predict_prophage_phigaro(extract_seqs_from_gbk.out.genome_ch)

    update_gbk_cds_annotation(samples.join(viral_annotation_pVOG.out.vanno_pvog, by:0).join(viral_annotation_VOGDB.out.vanno_vogdb, by:0).join(viral_annotation_PFAM.out.vanno_pfam, by:0).join(viral_annotation_KEGG.out.vanno_kegg, by:0).join(cds_anno_ARG.out.arg2update, by:0))
    update_gbk_prophage(update_gbk_cds_annotation.out.update_prophage.join(predict_prophage_phispy.out.phispy_tsv, by:0).join(predict_prophage_phigaro.out.phigaro_tsv, by:0))

    extract_prophages(update_gbk_prophage.out.prophages_ch)
    predict_CRISPR(extract_seqs_from_gbk.out.genome_ch)
}



