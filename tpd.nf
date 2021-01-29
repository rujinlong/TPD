#!/usr/bin/env nextflow
// Usage: nextflow run pipa.nf --fids fids.csv --mode "all|genome" -resume
// Usage: nextflow run tpd.nf -profile hpc_slurm --mode "all" -resume

nextflow.enable.dsl=2

process run_RagTag {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p01_RagTag"
    publishDir "$params.report/$sampleID", pattern: "*.stats"

    input:
    val(sampleID) 
    path(ref_fna)
    path(ref_gbk)

    output:
    tuple val(sampleID), path('ragtag_output/*') 
    tuple val(sampleID), path("${sampleID}_1k.fna"), emit: genome_ch
    tuple val(sampleID), path('*.stats') 

    when:
    params.mode == "all"

    """
    rawseqs="${baseDir}/data/${sampleID}_scaffolds.fna"
    ragtag.py scaffold $ref_fna \$rawseqs -t $task.cpus -o ragtag_output
    ln -s ragtag_output/*.stats .
    sed "s/^>/>${sampleID}_/" ragtag_output/ragtag.scaffolds.fasta | sed 's/ /_/g' > ${sampleID}.fna
    seqkit seq -m 1000 ${sampleID}.fna > ${sampleID}_1k.fna
    """
}

process run_dfast {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p02_dfast"
    publishDir "$params.report/$sampleID", pattern: "statistics.txt"

    input:
    tuple val(sampleID), path(genome)
    path(ref_gbk)

    output:
    path("dfast_output/*")
    tuple val(sampleID), path("dfast_output/genome.gbk"), emit: draft_gbk
    tuple val(sampleID), path("dfast_output/genome.fna"), emit: draft_genome_fna
    tuple val(sampleID), path("dfast_output/cds.fna"), emit: cds_ch
    tuple val(sampleID), path("protein_LOCUS.faa"), emit: protein_locus
    path("statistics.txt")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    dfast --genome ${genome} -o dfast_output --use_original_name t --cpu $task.cpus --minimum_length 30 --references ${ref_gbk}
    sed 's/>.*|LOCUS_/>LOCUS_/' dfast_output/protein.faa > protein_LOCUS.faa
    ln -s dfast_output/statistics.txt .
    """
}

process viral_annotation_VOGDB {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p03_vanno_VOGDB"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(protein)

    output:
    tuple val(sampleID), path("anno_VOGDB.tsv"), emit: vanno_vogdb

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_VOGDB $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,dbid > anno_temp.tsv
    csvtk join -t -f dbid anno_temp.tsv $params.desc_VOGDB > anno_VOGDB.tsv
    """
}

process viral_annotation_PFAM {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p03_vanno_PFAM"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(protein)

    output:
    tuple val(sampleID), path("anno_PFAM.tsv"), emit: vanno_pfam

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_PFAM $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,dbid > anno_temp.tsv
    csvtk join -t -f dbid anno_temp.tsv $params.desc_PFAM > anno_PFAM.tsv
    """
}

process viral_annotation_pVOG {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p03_vanno_pVOG"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(protein)

    output:
    tuple val(sampleID), path("anno_pVOG.tsv"), emit: vanno_pvog

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    db=pVOG
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_pVOG $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,dbid > anno_temp.tsv
    csvtk join -t -f dbid anno_temp.tsv $params.desc_pVOG > anno_pVOG.tsv
    """
}

process viral_annotation_KEGG {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p03_vanno_KEGG"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(protein)

    output:
    tuple val(sampleID), path("anno_KEGG.tsv"), emit: vanno_kegg

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_KEGG $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,koid > anno_temp.tsv
    csvtk join -t -f koid anno_temp.tsv $params.desc_KEGG | csvtk rename -t -f 2 -n dbid > anno_KEGG.tsv
    """
}

process cds_anno_ARG {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p03_anno_ARG"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(cds)

    output:
    tuple val(sampleID), path("anno_abricate.tsv"), emit: arg2update

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    for abrdb in argannot card ecoh ncbi plasmidfinder resfinder vfdb;do
        abricate --db \$abrdb $cds > ARG_\${abrdb}.tsv
    done

    head -n1 ARG_argannot.tsv | cut -f2- > anno_abricate.tsv
    cat ARG_* | grep -v "^#FILE" | cut -f2- | sort -k1,2 >> anno_abricate.tsv
    """
}

process predict_prophage_phispy {
    tag "$sampleID"
    label "big"
    publishDir "$params.outdir/$sampleID/p04_prophage_phispy"

    input:
    tuple val(sampleID), path(genome)

    output:
    path("out_phispy/*")
    tuple val(sampleID), path("prophage_phispy.tsv"), emit: phispy_tsv
    tuple val(sampleID), path("phispy.gbk"), emit: phispy_gbk

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    PhiSpy.py $genome -o out_phispy --phmms $params.db_phispy --threads $task.cpus --color --output_choice 31
    merge_two_gbks.py -r $genome -a out_phispy/genome.gbk -o phispy.gbk
    sed 1d out_phispy/prophage.tsv | cut -f2-4 > tmp2.tsv
    add_column.py -i tmp2.tsv -m phispy -o prophage_phispy.tsv
    """
}

process predict_prophage_phigaro {
    tag "$sampleID"
    label "big"
    publishDir "$params.outdir/$sampleID/p04_prophage_phigaro"

    input:
    tuple val(sampleID), path(genome)

    output:
    path("phigaro*")
    tuple val(sampleID), path("prophage_phigaro.tsv"), emit: phigaro_tsv

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    phigaro -f $genome -c $params.cfg_phigaro -e html tsv gff bed -o phigaro --not-open -t $task.cpus -m basic -d
    sed 1d phigaro/genome.phigaro.tsv | cut -f1-3 > tmp2.tsv
    add_column.py -i tmp2.tsv -m phigaro -o prophage_phigaro.tsv
    """
}

process update_gbk_cds_annotation {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p05_update_gbk_cds_annotation"

    input:
    tuple val(sampleID), path(dfast_gbk)
    tuple val(sampleID), path(pvog)
    tuple val(sampleID), path(vogdb)
    tuple val(sampleID), path(pfam)
    tuple val(sampleID), path(kegg)
    tuple val(sampleID), path(abricate) 

    output:
    tuple val(sampleID), path("${sampleID}_cds_anno.gbk"), emit: update_prophage

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    update_gbk_cds_annotation.py -g $dfast_gbk -a $pvog -o t1.gbk -d pVOG
    update_gbk_cds_annotation.py -g t1.gbk -a $vogdb -o t2.gbk -d VOGDB
    update_gbk_cds_annotation.py -g t2.gbk -a $pfam -o t3.gbk -d PFAM
    update_gbk_cds_annotation.py -g t3.gbk -a $kegg -o t4.gbk -d KEGG
    update_gbk_cds_annotation.py -g t4.gbk -a $abricate -o ${sampleID}_cds_anno.gbk -d abricate
    """
}

process update_gbk_prophage {
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p06_update_gbk_prophage"
    publishDir "$params.report/$sampleID", pattern: "*.gbk"

    input:
    val(sampleID)
    tuple val(sampleID), path(gbk)
    tuple val(sampleID), path(tsv_phispy)
    tuple val(sampleID), path(tsv_phigaro)

    output:
    tuple val(sampleID), path("${sampleID}_long.gbk"), emit: prophages_ch
    tuple val(sampleID), path("${sampleID}.gbk")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    cat $tsv_phispy $tsv_phigaro > prophages.tsv
    old_prophage="${baseDir}/data/${sampleID}_prophage.tsv"
    if [ -e \${old_prophage} ];then
        cat \${old_prophage} >> prophages.tsv
    fi
    add_prophage_to_gbk.py -g $gbk -p prophages.tsv -o ${sampleID}.gbk
    filter_contig_length.py -g ${sampleID}.gbk -m $params.contig_minlen -o ${sampleID}_long
    """
}

process extract_prophages {
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p07_prophages"
    publishDir "$params.report/$sampleID", pattern: "${sampleID}_prophages.*"

    input:
    tuple val(sampleID), path(prophages)

    output:
    path("${sampleID}_prophages.*")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    extract_prophages.py -s $sampleID -g $prophages -f $params.flank_len -o ${sampleID}_prophages
    """
}

process predict_CRISPR {
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p03_CRISPR"
    publishDir "$params.report/$sampleID"

    input:
    tuple val(sampleID), path(draft_genome_fna)

    output:
    path("${sampleID}_CRISPR_*")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    crt crt $draft_genome_fna ${sampleID}_CRISPR_CRT.txt
    cctyper $draft_genome_fna ${sampleID}_CRISPR_cctyper
    """
}


workflow {
    sampleIDs = channel.fromPath(params.data).map { it.toString().split("/")[-1].split("_")[0] }.unique()
    sampleIDs.view()
    ref_fna = params.ref_fna
    ref_gbk = params.ref_gbk
    run_RagTag(sampleIDs, ref_fna, ref_gbk)
    run_dfast(run_RagTag.out.genome_ch, ref_gbk)

    viral_annotation_VOGDB(run_dfast.out.protein_locus)
    viral_annotation_PFAM(run_dfast.out.protein_locus)
    viral_annotation_pVOG(run_dfast.out.protein_locus)
    viral_annotation_KEGG(run_dfast.out.protein_locus)
    cds_anno_ARG(run_dfast.out.cds_ch)

    predict_prophage_phispy(run_dfast.out.draft_gbk)
    predict_prophage_phigaro(run_dfast.out.draft_genome_fna)

    update_gbk_cds_annotation(run_dfast.out.draft_gbk, viral_annotation_pVOG.out.vanno_pvog, viral_annotation_VOGDB.out.vanno_vogdb, viral_annotation_PFAM.out.vanno_pfam, viral_annotation_KEGG.out.vanno_kegg, cds_anno_ARG.out.arg2update)
    update_gbk_prophage(sampleIDs, update_gbk_cds_annotation.out.update_prophage, predict_prophage_phispy.out.phispy_tsv, predict_prophage_phigaro.out.phigaro_tsv)

    extract_prophages(update_gbk_prophage.out.prophages_ch)
    predict_CRISPR(run_dfast.out.draft_genome_fna)
}



