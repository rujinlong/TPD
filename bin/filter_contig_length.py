#!/usr/bin/env python

import click
from Bio import SeqIO


@click.command()
@click.option("--gbk", '-g', help="Input file name")
@click.option("--min_len", '-m', default=1000, type=int, required=False, help="Minimum length of contig")
@click.option('--out', '-o', help="Prefix of output file name")
def main(gbk, min_len, out):
    records = list(SeqIO.parse(gbk, format="genbank"))
    rst = [x for x in records if len(x)>min_len]
    
    with open(out + ".gbk", "w") as fh:
        SeqIO.write(rst, fh, "genbank")
    
    with open(out + ".fasta", "w") as fh:
        SeqIO.write(rst, fh, "fasta")


if __name__ == '__main__':
    main()
