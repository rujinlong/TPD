/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { viral_annotation_VOGDB; viral_annotation_pVOG; viral_annotation_PFAM; viral_annotation_KEGG; cds_anno_ARG } from '../modules/anno_prot'
include { run_dfast; run_RagTag; predict_CRISPR } from '../modules/anno_genome'
include { predict_prophage_phispy; predict_prophage_phigaro; predict_prophage_phageboost; extract_prophages } from '../modules/prophage_miner'
include { update_gbk_cds_annotation; update_gbk_prophage } from '../modules/update'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TPD {
    data = "${params.datadir}/*_*"
    sampleIDs = channel.fromPath(data).map { it.toString().split("/")[-1].split("_")[0] }.unique()
    ref_fna = channel.fromPath(params.ref_fna)
    ref_gbk = channel.fromPath(params.ref_gbk)
    run_RagTag(sampleIDs.combine(ref_fna))
    run_dfast(run_RagTag.out.genome_ch.combine(ref_gbk))

    viral_annotation_VOGDB(run_dfast.out.protein_locus)
    viral_annotation_PFAM(run_dfast.out.protein_locus)
    viral_annotation_pVOG(run_dfast.out.protein_locus)
    viral_annotation_KEGG(run_dfast.out.protein_locus)
    cds_anno_ARG(run_dfast.out.cds_ch)

    predict_prophage_phispy(run_dfast.out.draft_gbk)
    predict_prophage_phigaro(run_dfast.out.draft_genome_fna)
    predict_prophage_phageboost(run_dfast.out.draft_genome_fna)

    update_gbk_cds_annotation(run_dfast.out.draft_gbk.join(viral_annotation_pVOG.out.vanno_pvog, by:0).join(viral_annotation_VOGDB.out.vanno_vogdb, by:0).join(viral_annotation_PFAM.out.vanno_pfam, by:0).join(viral_annotation_KEGG.out.vanno_kegg, by:0).join(cds_anno_ARG.out.arg2update, by:0))
    update_gbk_prophage(update_gbk_cds_annotation.out.update_prophage.join(predict_prophage_phispy.out.phispy_tsv, by:0).join(predict_prophage_phigaro.out.phigaro_tsv, by:0).join(predict_prophage_phageboost.out.phageboost_tsv, by:0))

    extract_prophages(update_gbk_prophage.out.prophages_ch)
    predict_CRISPR(run_dfast.out.draft_genome_fna)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete { 
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/