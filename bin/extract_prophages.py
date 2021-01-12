#!/usr/bin/env python

import click
from Bio import SeqIO


def extract_prophage(contig, misc_rec, flank_length):
    start = misc_rec.location.start - flank_length
    end = misc_rec.location.end + flank_length
    prophage = contig[start:end]
    prophage.name = "{}_{}:{}_{}_prophage".format(contig.id, start, end, misc_rec.strand)
    prophage.id = prophage.name
    prophage.type = "Prophage"
    return prophage


@click.command()
@click.option("--sample_id", '-s', help="Given an ID for your sample")
@click.option("--gbk", '-g', help="Input file name")
@click.option('--flank_length', '-f', default=1000, type=int, help="Flanking length")
@click.option("--fout_prefix", '-o', help="Prefix of output file name")
def main(sample_id, gbk, flank_length, fout_prefix):
    contigs = list(SeqIO.parse(gbk, format="genbank"))
    prophages = []
    for contig in contigs:
        if not contig.id.startswith(sample_id):
            contig.id = "{}".format(sample_id)
            contig.name = contig.id
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