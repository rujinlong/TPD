#!/usr/bin/env python

import click
from Bio import SeqIO


@click.command()
@click.option("--gbk", '-g', help="Input Genbank file")
@click.option("--fna", '-f', help="Output fasta file")
def main(gbk, fna):
    genbank_file = open (gbk, "r")
    output_file = open(fna, "w")
    records = SeqIO.parse(genbank_file, "genbank")
    SeqIO.write(records, output_file, "fasta")
    output_file.close()
    genbank_file.close()


if __name__ == '__main__':
    main()

