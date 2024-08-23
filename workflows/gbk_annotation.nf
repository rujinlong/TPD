/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { viral_annotation_VOGDB; viral_annotation_pVOG; viral_annotation_PFAM; viral_annotation_KEGG; cds_anno_ARG } from '../modules/anno_prot'
include { predict_CRISPR } from '../modules/anno_genome'
include { predict_prophage_phispy; predict_prophage_phigaro; predict_prophage_phageboost; extract_prophages } from '../modules/prophage_miner'
include { update_gbk_cds_annotation; update_gbk_prophage } from '../modules/update'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GBKANNOTATION {
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
    predict_prophage_phageboost(extract_seqs_from_gbk.out.genome_ch)

    update_gbk_cds_annotation(samples.join(viral_annotation_pVOG.out.vanno_pvog, by:0).join(viral_annotation_VOGDB.out.vanno_vogdb, by:0).join(viral_annotation_PFAM.out.vanno_pfam, by:0).join(viral_annotation_KEGG.out.vanno_kegg, by:0).join(cds_anno_ARG.out.arg2update, by:0))
    update_gbk_prophage(update_gbk_cds_annotation.out.update_prophage.join(predict_prophage_phispy.out.phispy_tsv, by:0).join(predict_prophage_phigaro.out.phigaro_tsv, by:0).join(predict_prophage_phageboost.out.phageboost_tsv, by:0))

    extract_prophages(update_gbk_prophage.out.prophages_ch)
    predict_CRISPR(extract_seqs_from_gbk.out.genome_ch)
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
