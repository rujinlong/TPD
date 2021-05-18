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