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