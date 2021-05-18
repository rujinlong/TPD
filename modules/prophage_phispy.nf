process predict_prophage_phispy {
    tag "$sampleID"
    label "big"
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
    ln -s $gbk genome.gbk
    PhiSpy.py genome.gbk -o out_phispy --phmms $params.db_phispy --threads $task.cpus --color --output_choice 31
    merge_two_gbks.py -r genome.gbk -a out_phispy/genome.gbk -o phispy.gbk
    sed 1d out_phispy/prophage.tsv | cut -f2-4 > tmp2.tsv
    add_column.py -i tmp2.tsv -m phispy -o prophage_phispy.tsv
    """
}