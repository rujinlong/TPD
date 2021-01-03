#!/usr/bin/env python

import click
from Bio import SeqIO


@click.command()
@click.option("--gbk_all", '-r', help="Raw genbank file")
@click.option("--gbk_less", '-a', help="Genbank file which removed some contigs")
@click.option('--gbk_out', '-o', help="Output file name")
def main(gbk_all, gbk_less, gbk_out):
    recs_all = list(SeqIO.parse(gbk_all, format="genbank"))
    recs_less = list(SeqIO.parse(gbk_less, format="genbank"))
    removed_ids = set([x.id for x in recs_all]) - set([x.id for x in recs_less])
    removed_recs = [x for x in recs_all if x.id in removed_ids]
    records = recs_less + removed_recs
    
    with open(gbk_out, "w") as fh:
        SeqIO.write(records, fh, "genbank")


if __name__ == '__main__':
    main()