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