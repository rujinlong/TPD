process predict_prophage_phigaro {
    tag "$sampleID"
    label "big"
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