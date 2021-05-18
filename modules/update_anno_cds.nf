process update_gbk_cds_annotation {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p04_update_gbk_cds_annotation"

    input:
    tuple val(sampleID), path(original_gbk), path(pvog), path(vogdb), path(pfam), path(kegg), path(abricate) 

    output:
    tuple val(sampleID), path("${sampleID}_cds_anno.gbk"), emit: update_prophage

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    update_gbk_cds_annotation.py -g $original_gbk -a $pvog -o t1.gbk -d pVOG
    update_gbk_cds_annotation.py -g t1.gbk -a $vogdb -o t2.gbk -d VOGDB
    update_gbk_cds_annotation.py -g t2.gbk -a $pfam -o t3.gbk -d PFAM
    update_gbk_cds_annotation.py -g t3.gbk -a $kegg -o t4.gbk -d KEGG
    update_gbk_cds_annotation.py -g t4.gbk -a $abricate -o ${sampleID}_cds_anno.gbk -d abricate
    """
}