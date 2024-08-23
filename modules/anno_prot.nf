process viral_annotation_VOGDB {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p02_vanno_VOGDB"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(protein)

    output:
    tuple val(sampleID), path("anno_VOGDB.tsv"), emit: vanno_vogdb

    when:
    params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_VOGDB $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,dbid > anno_temp.tsv
    csvtk join -t -f dbid anno_temp.tsv $params.desc_VOGDB > anno_VOGDB.tsv
    """
}


process viral_annotation_pVOG {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p02_vanno_pVOG"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(protein)

    output:
    tuple val(sampleID), path("anno_pVOG.tsv"), emit: vanno_pvog

    when:
    params.mode == "all"

    """
    db=pVOG
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_pVOG $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,dbid > anno_temp.tsv
    csvtk join -t -f dbid anno_temp.tsv $params.desc_pVOG > anno_pVOG.tsv
    """
}


process viral_annotation_PFAM {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p02_vanno_PFAM"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(protein)

    output:
    tuple val(sampleID), path("anno_PFAM.tsv"), emit: vanno_pfam

    when:
    params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_PFAM $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,dbid > anno_temp.tsv
    csvtk join -t -f dbid anno_temp.tsv $params.desc_PFAM > anno_PFAM.tsv
    """
}


process viral_annotation_KEGG {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p02_vanno_KEGG"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(protein)

    output:
    tuple val(sampleID), path("anno_KEGG.tsv"), emit: vanno_kegg

    when:
    params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_KEGG $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,koid > anno_temp.tsv
    csvtk join -t -f koid anno_temp.tsv $params.desc_KEGG | csvtk rename -t -f 2 -n dbid > anno_KEGG.tsv
    """
}


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