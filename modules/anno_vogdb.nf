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