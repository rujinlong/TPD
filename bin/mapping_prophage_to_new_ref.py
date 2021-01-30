#!/usr/bin/env python

import click
from Bio import SeqIO
import pandas as pd
from pymummer import coords_file, alignment, nucmer


def run_nucmer(qry, ref):
    with open('tmp_qry.fna', 'w') as fh:
        SeqIO.write(qry, fh, "fasta")
    with open('tmp_ref.fna', 'w') as fh:
        SeqIO.write(ref, fh, "fasta")

    runner = nucmer.Runner("tmp_ref.fna", "tmp_qry.fna", "tmp_mm.out", maxgap=200, mincluster=1000) 
    runner.run()
    file_reader = coords_file.reader("tmp_mm.out")
    alignments = [coord for coord in file_reader if not coord.is_self_hit()]
    return alignments


@click.command()
@click.option("--fref", '-r', help="Host reference genome")
@click.option("--fprophage", '-p', help="Predicted prophage FASTA")
@click.option("--fout", '-o', help="output file name")
def main(fref, fprophage, fout):
    refs = list(SeqIO.parse(fref, "genbank"))
    prophages = list(SeqIO.parse(fprophage, "fasta"))
    
    alns = []
    for qry in prophages:
        for ref in refs:
            aln = run_nucmer(qry, ref)
            if len(aln) > 0:
                ref_start = min([x.ref_start for x in aln])
                ref_end = max([x.ref_end for x in aln])
                alns.append([ref.id, ref_start, ref_end, qry.id])
    
    df = pd.DataFrame(alns)
    df.to_csv(fout, index=False, header=False, sep='\t')

if __name__ == '__main__':
    main()