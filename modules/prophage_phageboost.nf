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