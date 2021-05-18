process extract_prophages {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p06_prophages"
    publishDir "$params.report/$sampleID", pattern: "${sampleID}_prophages.*"

    input:
    tuple val(sampleID), path(prophages)

    output:
    path("${sampleID}_prophages*")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    extract_prophages.py -g $prophages -o ${sampleID}_prophages_all_flank -p misc_feature -t total -l $params.flank_len
    extract_prophages.py -g $prophages -o ${sampleID}_prophages_merged_flank -p misc_feature -t merged -l $params.flank_len
    extract_prophages.py -g $prophages -o ${sampleID}_prophages -p misc_feature -t merged -l 0
    """
}