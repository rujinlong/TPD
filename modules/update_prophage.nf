process update_gbk_prophage {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p06_update_gbk_prophage"
    publishDir "$params.report/$sampleID", pattern: "*.gbk"

    input:
    tuple val(sampleID), path(gbk), path(tsv_phispy), path(tsv_phigaro), path(tsv_phageboost)

    output:
    tuple val(sampleID), path("${sampleID}_long.gbk"), emit: prophages_ch
    tuple val(sampleID), path("${sampleID}.gbk")

    when:
    params.mode == "all"

    """
    awk 'FNR==1{}{print}' $tsv_phispy $tsv_phigaro $tsv_phageboost > predicted_prophages.tsv
    merge_prophage_overlapping.py -i predicted_prophages.tsv -o merged_prophages.tsv
    cat predicted_prophages.tsv merged_prophages.tsv > prophages.tsv
    old_prophage="${params.datadir}/${sampleID}${params.old_prophage_ext}"
    if [ -e \${old_prophage} ];then
        if [[ "$params.old_prophage_ext" = "_old.gbk" ]];then
            extract_prophages.py -s $sampleID -g ${params.datadir}/${sampleID}_old.gbk -f 0 -o ${sampleID}_prophage -p $params.old_prophage_feature_name
	    elif [[ "$params.old_prophage_ext" = "_old.fna" ]];then
            ln -s ${params.datadir}/${sampleID}_old.fna ${sampleID}_prophage.fasta
        fi
        mapping_prophage_to_new_ref.py -r $gbk -p ${sampleID}_prophage.fasta -o ${sampleID}_prophage.tsv
        awk 'FNR==1{}{print}' ${sampleID}_prophage.tsv >> prophages.tsv
    fi
    add_prophage_to_gbk.py -g $gbk -p prophages.tsv -o ${sampleID}.gbk
    filter_contig_length.py -g ${sampleID}.gbk -m $params.contig_minlen -o ${sampleID}_long
    """
}