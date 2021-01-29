#!/usr/bin/env nextflow
// Usage: nextflow run pipa.nf --fids fids.csv --mode "all|genome" -resume

nextflow.enable.dsl=2

process run_RagTag {
    label "small"
    publishDir "$params.outdir/$params.prefix/p01_RagTag"
    publishDir "$params.report/$params.prefix", pattern: "*.stats"

    input:
    path(rawseqs) 
    path(ref_fna)
    path(ref_gbk)

    output:
    path('ragtag_output/*') 
    path('ragtag_output/ragtag.scaffolds.fasta'), emit: genome_ch
    path('*.stats') 

    when:
    params.mode == "all"

    """
    ragtag.py scaffold $ref_fna $rawseqs -t $task.cpus -o ragtag_output
    ln -s ragtag_output/*.stats .
    """
}


process run_dfast {
    label "medium"
    publishDir "$params.outdir/$params.prefix/p02_dfast"
    publishDir "$params.report/$params.prefix", pattern: "statistics.txt"

    input:
    path(genome)
    path(ref_gbk)

    output:
    file("dfast_output/*")
    path("dfast_output/genome.gbk"), emit: draft_gbk
    path("dfast_output/genome.fna"), emit: draft_genome_fna
    path("dfast_output/cds.fna"), emit: cds_ch
    path("protein_LOCUS.faa"), emit: protein_locus
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
    label "medium"
    publishDir "$params.outdir/$params.prefix/p03_vanno_VOGDB"
    publishDir "$params.report/$params.prefix"

    input:
    path(protein)

    output:
    path("anno_VOGDB.tsv"), emit: vanno_vogdb

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_VOGDB $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,dbid > anno_temp.tsv
    csvtk join -t -f dbid anno_temp.tsv $params.desc_VOGDB > anno_VOGDB.tsv
    """
}

process viral_annotation_PFAM {
    label "medium"
    publishDir "$params.outdir/$params.prefix/p03_vanno_PFAM"
    publishDir "$params.report/$params.prefix"

    input:
    path(protein)

    output:
    path("anno_PFAM.tsv"), emit: vanno_pfam

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_PFAM $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,dbid > anno_temp.tsv
    csvtk join -t -f dbid anno_temp.tsv $params.desc_PFAM > anno_PFAM.tsv
    """
}

process viral_annotation_pVOG {
    label "medium"
    publishDir "$params.outdir/$params.prefix/p03_vanno_pVOG"
    publishDir "$params.report/$params.prefix"

    input:
    path(protein)

    output:
    path("anno_pVOG.tsv"), emit: vanno_pvog

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
    label "medium"
    publishDir "$params.outdir/$params.prefix/p03_vanno_KEGG"
    publishDir "$params.report/$params.prefix"

    input:
    path(protein)

    output:
    path("anno_KEGG.tsv"), emit: vanno_kegg

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    hmmsearch --tblout hit.hmmtbl --noali -T 40 --cpu $task.cpus -o hit_temp.txt $params.db_KEGG $protein
    grep -v '^#' hit.hmmtbl | awk '{print \$1,\$3}' | csvtk space2tab | csvtk add-header -t -n protid,koid > anno_temp.tsv
    csvtk join -t -f koid anno_temp.tsv $params.desc_KEGG | csvtk rename -t -f 2 -n dbid > anno_KEGG.tsv
    """
}

process cds_anno_ARG {
    label "small"
    publishDir "$params.outdir/$params.prefix/p03_anno_ARG"
    publishDir "$params.report/$params.prefix"

    input:
    path(cds)

    output:
    path("anno_abricate.tsv"), emit: arg2update

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
    label "big"
    publishDir "$params.outdir/$params.prefix/p04_prophage_phispy"

    input:
    path(genome)

    output:
    path("out_phispy/*")
    path("prophage_phispy.tsv"), emit: phispy_tsv
    path("phispy.gbk"), emit: phispy_gbk

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    PhiSpy.py $genome -o out_phispy --phmms $params.db_phispy --threads $task.cpus --color --output_choice 31
    merge_two_gbks.py -r $genome -a out_phispy/genome.gbk -o phispy.gbk
    sed 1d prophage.tsv | cut -f2-4 > tmp2.tsv
    add_column.py -i tmp2.tsv -m phispy -o prophage_phispy.tsv
    """
}

process predict_prophage_phigaro {
    label "big"
    publishDir "$params.outdir/$params.prefix/p04_prophage_phigaro"

    input:
    path(genome)

    output:
    path("phigaro*")
    path("prophage_phigaro.tsv"), emit: phigaro_tsv

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    phigaro -f $genome -c $params.cfg_phigaro -e html tsv gff bed -o phigaro --not-open -t $task.cpus -m basic -d
    sed 1d phigaro/genome.phigaro.tsv | cut -f1-3 > tmp2.tsv
    add_column.py -i tmp2.tsv -m phigaro -o prophage_phigaro.tsv
    """
}

process update_gbk_cds_annotation {
    label "small"
    publishDir "$params.outdir/$params.prefix/p05_update_gbk_cds_annotation"

    input:
    path(phispy_gbk)
    path(dfast_gbk)
    path(pvog)
    path(vogdb)
    path(pfam)
    path(kegg)
    path(abricate) 

    output:
    path("${params.prefix}_cds_anno.gbk"), emit: update_prophage

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    update_gbk_cds_annotation.py -g $phispy_gbk -a $pvog -o t1.gbk -d pVOG
    update_gbk_cds_annotation.py -g t1.gbk -a $vogdb -o t2.gbk -d VOGDB
    update_gbk_cds_annotation.py -g t2.gbk -a $pfam -o t3.gbk -d PFAM
    update_gbk_cds_annotation.py -g t3.gbk -a $kegg -o t4.gbk -d KEGG
    update_gbk_cds_annotation.py -g t4.gbk -a $abricate -o ${params.prefix}_cds_anno.gbk -d abricate
    """
}

process update_gbk_prophage {
    publishDir "$params.outdir/$params.prefix/p06_update_gbk_prophage"
    publishDir "$params.report/$params.prefix", pattern: "*.gbk"

    input:
    path(gbk)
    path(tsvs)

    output:
    path("${params.prefix}_long.gbk"), emit: prophages_ch
    path("${params.prefix}.gbk")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    cat $tsvs > prophages.tsv
    add_prophage_to_gbk.py -g $gbk -p prophages.tsv -o ${params.prefix}.gbk
    filter_contig_length.py -g ${params.prefix}.gbk -m $params.contig_minlen -o ${params.prefix}_long
    """
}

process extract_prophages {
    publishDir "$params.outdir/$params.prefix/p07_prophages"
    publishDir "$params.report/$params.prefix", pattern: "${params.prefix}_prophages.*"

    input:
    path(prophages)

    output:
    path("${params.prefix}_prophages.*")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    extract_prophages.py -s $params.prefix -g $prophages -f $params.flank_len -o ${params.prefix}_prophages
    """
}

process predict_CRISPR {
    publishDir "$params.outdir/$params.prefix/p03_CRISPR"
    publishDir "$params.report/$params.prefix"

    input:
    path(draft_genome_fna)

    output:
    path("${params.prefix}_CRISPR_*")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    crt crt $draft_genome_fna ${params.prefix}_CRISPR_CRT.txt
    cctyper $draft_genome_fna ${params.prefix}_CRISPR_cctyper.txt
    """
}


workflow {
    contigs = channel.fromPath(params.contigs)
    ref_fna = channel.fromPath(params.ref_fna)
    ref_gbk = channel.fromPath(params.ref_gbk)
    run_RagTag(contigs, ref_fna, ref_gbk)
    run_dfast(run_RagTag.out.genome_ch, ref_gbk)

    viral_annotation_VOGDB(run_dfast.out.protein_locus)
    viral_annotation_PFAM(run_dfast.out.protein_locus)
    viral_annotation_pVOG(run_dfast.out.protein_locus)
    viral_annotation_KEGG(run_dfast.out.protein_locus)
    cds_anno_ARG(run_dfast.out.cds_ch)

    predict_prophage_phispy(run_dfast.out.draft_gbk)
    predict_prophage_phigaro(run_dfast.out.draft_genome_fna)

    update_gbk_cds_annotation(predict_prophage_phispy.out.phispy_gbk, run_dfast.out.draft_gbk, viral_annotation_pVOG.out.vanno_pvog, viral_annotation_VOGDB.out.vanno_vogdb, viral_annotation_PFAM.out.vanno_pfam, viral_annotation_KEGG.out.vanno_kegg, cds_anno_ARG.out.arg2update)
    
    if( params.extra_prophage_coord == "false" )
        tsvs = predict_prophage_phispy.out.phispy_tsv.mix(predict_prophage_phigaro.out.phigaro_tsv)
    else
        tsvs = channel.fromPath(params.extra_prophage_coord).mix(predict_prophage_phispy.out.phispy_tsv, predict_prophage_phigaro.out.phigaro_tsv)

    update_gbk_prophage(update_gbk_cds_annotation.out.update_prophage, tsvs)

    extract_prophages(update_gbk_prophage.out.prophages_ch)
    predict_CRISPR(run_dfast.out.draft_genome_fna)
}



