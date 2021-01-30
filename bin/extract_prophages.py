#!/usr/bin/env python

import click
from Bio import SeqIO


def extract_prophage(sample_id, rec, prophage_feature, flank_length):
    start = prophage_feature.location.start - flank_length
    end = prophage_feature.location.end + flank_length
    length = abs(start-end)
    prophage = rec[start:end]
    if not prophage.id.startswith(sample_id):
        prophage.id = prophage_feature.qualifiers['ID'][0]
    else:
        prophage.id = "{}_{}_len{}_{}-{}".format(prophage.id, prophage_feature.qualifiers['ID'][0], str(length), start, end)
    prophage.name = prophage.id
    prophage.type = "prophage"
    return prophage


@click.command()
@click.option("--sample_id", '-s', help="Given an ID for your sample")
@click.option("--gbk", '-g', help="Input file name")
@click.option('--flank_length', '-f', default=1000, type=int, help="Flanking length")
@click.option("--fout_prefix", '-o', help="Prefix of output file name")
@click.option('--prophage_type_name', '-p', default="prophage", help="prophage type name in genbank file. [prophage | mics_feature]")
def main(sample_id, gbk, flank_length, fout_prefix, prophage_type_name):
    gbk_recs = list(SeqIO.parse(gbk, format="genbank"))
    prophages = []
    for rec in gbk_recs:
        # if not rec.id.startswith(sample_id):
        #     rec.id = "{}_{}".format(sample_id, rec.id)
            # rec.name = rec.id
        prophage_features = [x for x in rec.features if x.type==prophage_type_name]
        for pf in prophage_features:
            prophage = extract_prophage(sample_id, rec, pf, flank_length)
            prophages.append(prophage)
    
    with open("{}.gbk".format(fout_prefix), "w") as fh:
        SeqIO.write(prophages, fh, "genbank")
    
    with open("{}.fasta".format(fout_prefix), "w") as fh:
        SeqIO.write(prophages, fh, "fasta")
    

if __name__ == '__main__':
    main()