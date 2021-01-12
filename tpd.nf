#!/usr/bin/env nextflow
// Usage: nextflow run pipa.nf --fids fids.csv --mode "all|genome" -resume

if (params.mode == "all") {
    Channel
        .fromPath(params.fids)
        .splitCsv(header:true)
        .map { row -> tuple(row.sampleID, file(row.rawseqs), file(row.ref_fna), file(row.ref_gbk)) }
        .set { contigs_ch }

    process run_RagTag {
        tag "$sampleID"
        label "small"
        publishDir "$params.outdir/$sampleID/p01_RagTag"
        publishDir "$params.report/$sampleID", pattern: "*.stats"
        conda '/home/viro/jinlong.ru/conda3/envs/ragtag'

        input:
        set val(sampleID), file(rawseqs), file(ref_fna), file(ref_gbk) from contigs_ch

        output:
        tuple val(sampleID), file('ragtag_output/*') 
        tuple val(sampleID), file('ragtag_output/ragtag.scaffolds.fasta'),  file(ref_gbk) into genome_ch
        tuple val(sampleID), file('*.stats') 

        when:
        params.mode == "all"

        """
        ragtag.py scaffold $ref_fna $rawseqs -t $task.cpus -o ragtag_output
        ln -s ragtag_output/*.stats .
        """
    }
} else {
    Channel
        .fromPath(params.fids)
        .splitCsv(header:true)
        .map { row -> tuple(row.sampleID, file(row.rawseqs), file(row.ref_gbk)) }
        .set { genome_ch }
}


process run_dfast {
    tag "$sampleID"
    label "medium"
    publishDir "$params.outdir/$sampleID/p02_dfast"
    publishDir "$params.report/$sampleID", pattern: "statistics.txt"
    conda '/home/viro/jinlong.ru/conda3/envs/dfast'

    input:
    set val(sampleID), file(genome), file(ref_gbk) from genome_ch

    output:
    tuple val(sampleID), file("dfast_output/*")
    tuple val(sampleID), file("dfast_output/genome.gbk") into dfast2phispy
    tuple val(sampleID), file("dfast_output/genome.gbk") into dfast2update
    tuple val(sampleID), file("dfast_output/genome.fna") into dfast2phigaro
    tuple val(sampleID), file("dfast_output/genome.fna") into dfast2CRT
    tuple val(sampleID), file("dfast_output/cds.fna") into cds2abricate
    tuple val(sampleID), file("protein_LOCUS.faa") into protein2vogdb
    tuple val(sampleID), file("protein_LOCUS.faa") into protein2pfam
    tuple val(sampleID), file("protein_LOCUS.faa") into protein2pvog
    tuple val(sampleID), file("protein_LOCUS.faa") into protein2kegg
    tuple val(sampleID), file("statistics.txt")

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
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    set val(sampleID), file(protein) from protein2vogdb

    output:
    tuple val(sampleID), file("anno_VOGDB.tsv") into vanno_vogdb

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
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    set val(sampleID), file(protein) from protein2pfam

    output:
    tuple val(sampleID), file("anno_PFAM.tsv") into vanno_pfam

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
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    set val(sampleID), file(protein) from protein2pvog

    output:
    tuple val(sampleID), file("anno_pVOG.tsv") into vanno_pvog

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
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    set val(sampleID), file(protein) from protein2kegg

    output:
    tuple val(sampleID), file("anno_KEGG.tsv") into vanno_kegg

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
    conda '/home/viro/jinlong.ru/conda3/envs/binfo-ng'

    input:
    set val(sampleID), file(cds) from cds2abricate

    output:
    tuple val(sampleID), file("anno_abricate.tsv") into arg2update

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
    conda '/home/viro/jinlong.ru/conda3/envs/phispy'

    input:
    set val(sampleID), file(genome) from dfast2phispy

    output:
    tuple val(sampleID), file("out_phispy/*")
    tuple val(sampleID), file("out_phispy/prophage.tsv") into coord_phispy
    tuple val(sampleID), file("phispy.gbk") into phispy2update

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    PhiSpy.py $genome -o out_phispy --phmms $params.db_phispy --threads $task.cpus --color --output_choice 31
    merge_two_gbks.py -r $genome -a out_phispy/genome.gbk -o phispy.gbk
    """
}


process predict_prophage_phigaro {
    tag "$sampleID"
    label "big"
    publishDir "$params.outdir/$sampleID/p04_prophage_phigaro"
    conda '/home/viro/jinlong.ru/conda3/envs/phigaro_env'

    input:
    set val(sampleID), file(genome) from dfast2phigaro

    output:
    tuple val(sampleID), file("phigaro*")
    tuple val(sampleID), file("phigaro.tsv") into coord_phigaro

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    phigaro -f $genome -c $params.cfg_phigaro -e html tsv gff bed -o phigaro --not-open -t $task.cpus -m basic -d
    """
}


process update_gbk_cds_annotation {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/$sampleID/p05_update_gbk_cds_annotation"
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    set val(sampleID), file(phispy_gbk), file(dfast_gbk), file(pvog), file(vogdb), file(pfam), file(kegg), file(abricate) from phispy2update.join(dfast2update).join(vanno_pvog).join(vanno_vogdb).join(vanno_pfam).join(vanno_kegg).join(arg2update)

    output:
    tuple val(sampleID), file("${sampleID}_cds_anno.gbk") into update_prophage

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    update_gbk_cds_annotation.py -g $phispy_gbk -a $pvog -o t1.gbk -d pVOG
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
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    set val(sampleID), file(prophage_phigaro), file(gbk) from coord_phigaro.join(update_prophage)

    output:
    tuple val(sampleID), file("${sampleID}_long.gbk") into prophages_ch
    tuple val(sampleID), file("${sampleID}.gbk")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    update_gbk_prophage.py -g $gbk -i $prophage_phigaro -o ${sampleID}.gbk
    filter_contig_length.py -g ${sampleID}.gbk -m $params.contig_minlen -o ${sampleID}_long
    """
}

process extract_prophages {
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p07_prophages"
    publishDir "$params.report/$sampleID", pattern: "${sampleID}_prophages.*"
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    set val(sampleID), file(prophages) from prophages_ch

    output:
    tuple val(sampleID), file("${sampleID}_prophages.*")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    extract_prophages.py -s $sampleID -g $prophages -f $params.flank_len -o ${sampleID}_prophages
    """
}


process prophage_clustering {
    publishDir "$params.outdir/$sampleID/p08_prophage_clusters"
    publishDir "$params.report/$sampleID", pattern: "clusters_reps.*"
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    file(prophages_all) from prophages_ch

    output:
    file("clusters_reps.fna")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    cat ${prophages_all} > prophages_all.fna
    sortgenome.pl --genomes-file prophages_all.fna --sortedgenomes-file sorted.fna
    gclust -both -nuc -threads $task.cpus -memiden 4 sorted.fna > clusters.txt
    pretty_cdhit.py -i clusters.txt -o clusters.map
    cut -f1 clusters.map | sed '1d' | sort -u > clusters.list
    seqkit grep -f clusters.list prophages_all.fna > clusters_reps.fna
    # todo: jpnb
    """
}


process predict_CRISPR {
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p03_CRISPR"
    publishDir "$params.report/$sampleID", pattern: "${sampleID}_CRISPR*"
    conda '/home/viro/jinlong.ru/conda3/envs/tpd'

    input:
    set val(sampleID), file(draft_genome) from dfast2CRT

    output:
    tuple val(sampleID), file("${sampleID}_CRISPR_*")

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    crt crt $draft_genome ${sampleID}_CRISPR_CRT.txt
    cctyper $draft_genome ${sampleID}_CRISPR_cctyper.txt
    """
}


workflow.onComplete { 
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}


