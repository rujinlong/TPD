#!/usr/bin/env python

import click
from Bio import SeqIO


def extract_prophage(bacteria_rec, misc_rec, flank_length):
    start = misc_rec.location.start - flank_length
    end = misc_rec.location.end + flank_length
    prophage = bacteria_rec[start:end]
    prophage.name = "Prophage_{}:{}_{}_{}".format( start, end, misc_rec.strand, prophage.name)
    prophage.id = prophage.name
    prophage.type = "Prophage"
    return prophage


@click.command()
@click.option("--gbk", '-g', help="input file name")
@click.option('--flank_length', '-f', default=1000, type=int, help="flanking length")
@click.option("--fout_prefix", '-o', help="Prefix of output file name")
def main(gbk, flank_length, fout_prefix):
    contigs = list(SeqIO.parse(gbk, format="genbank"))
    prophages = []
    for contig in contigs:
        misc_features = [x for x in contig.features if x.type=="misc_feature"]
        for misc in misc_features:
            prophage = extract_prophage(contig, misc, flank_length)
            prophages.append(prophage)
    with open("{}.gbk".format(fout_prefix), "w") as fh:
        SeqIO.write(prophages, fh, "genbank")
    
    with open("{}.fasta".format(fout_prefix), "w") as fh:
        SeqIO.write(prophages, fh, "fasta")
    

if __name__ == '__main__':
    main()