process predict_prophage_phageboost {
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p04_prophage_phageboost"

    input:
    tuple val(sampleID), path(genome)

    output:
    path("phageboost_results/*")
    tuple val(sampleID), path("prophage_phageboost.tsv"), emit: phageboost_tsv

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    seqkit seq -m 20000 $genome > genome20k.fna
    PhageBoost -f genome20k.fna -o phageboost_results
    grep -v '^#' phageboost_results/phages_*.gff | cut -f1,4,5 | sed 's/n_genes.*=//' > tmp2.tsv
    add_column.py -i tmp2.tsv -m phageboost -o prophage_phageboost.tsv
    """
}


process predict_prophage_phigaro {
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p03_prophage_phigaro"

    input:
    tuple val(sampleID), path(genome)

    output:
    path("phigaro*")
    tuple val(sampleID), path("prophage_phigaro.tsv"), emit: phigaro_tsv

    when:
    params.mode == 'genome' || params.mode == "all"

    """
    ln -s $genome genome.fasta
    phigaro -f genome.fasta -c $params.cfg_phigaro -e html tsv gff bed -o phigaro --not-open -t $task.cpus -m basic -d
    sed 1d phigaro/genome.phigaro.tsv | cut -f1-3 > tmp2.tsv
    add_column.py -i tmp2.tsv -m phigaro -o prophage_phigaro.tsv
    """
}


process predict_prophage_phispy {
    tag "$sampleID"
    publishDir "$params.outdir/$sampleID/p03_prophage_phispy"

    input:
    tuple val(sampleID), path(gbk)

    output:
    path("out_phispy/*")
    tuple val(sampleID), path("prophage_phispy.tsv"), emit: phispy_tsv
    tuple val(sampleID), path("phispy.gbk"), emit: phispy_gbk

    when:
    params.mode == "all"

    """
    PhiSpy.py $gbk -o out_phispy --phmms $params.db_phispy --threads $task.cpus --color --output_choice 31
    merge_two_gbks.py -r $gbk -a out_phispy/genome.gbk -o phispy.gbk
    sed 1d out_phispy/prophage.tsv | cut -f2-4 > tmp2.tsv
    add_column.py -i tmp2.tsv -m phispy -o prophage_phispy.tsv
    """
}


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