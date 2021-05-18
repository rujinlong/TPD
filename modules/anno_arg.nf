process cds_anno_ARG {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p02_anno_ARG"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(cds)

    output:
    tuple val(sampleID), path("anno_abricate.tsv"), emit: arg2update

    when:
    params.mode == "all"

    """
    for abrdb in argannot card ecoh ncbi plasmidfinder resfinder vfdb;do
        abricate --db \$abrdb $cds > ARG_\${abrdb}.tsv
    done

    head -n1 ARG_argannot.tsv | cut -f2- > anno_abricate.tsv
    cat ARG_* | grep -v "^#FILE" | cut -f2- | sort -k1,2 >> anno_abricate.tsv
    """
}