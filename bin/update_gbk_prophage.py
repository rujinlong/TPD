#!/usr/bin/env python

import click
from Bio import SeqFeature, SeqIO
import pandas as pd


def add_phigaro_prophage(record, prophage):
    start_pos = SeqFeature.ExactPosition(prophage['begin'])
    end_pos = SeqFeature.ExactPosition(prophage['end'])
    feature_location = SeqFeature.FeatureLocation(start_pos,end_pos, strand=1)
    qualifiers = {"transposable": prophage['transposable'], "taxonomy": prophage['taxonomy'], "note": "prophage region identified with Phigaro"}

    new_feature = SeqFeature.SeqFeature(feature_location,type="misc_feature", qualifiers=qualifiers)
    record.features.append(new_feature)
    return record


@click.command()
@click.option("--gbk", '-g', help="Genbank file contains all contigs")
@click.option("--prophage", '-i', help="Predicted prophage location")
@click.option('--fout', '-o', help="Output file name")
@click.option("--predictor", '-p', default="phigaro", required=False, help="Prediction method")
def main(gbk, prophage, predictor, fout):
    records = list(SeqIO.parse(gbk, format="genbank"))
    if predictor == "phigaro":
        prophages = pd.read_csv(prophage, sep='\t').to_dict(orient="records")
    
    for rec in records:
        for pg in prophages:
            if rec.id == pg['scaffold']:
                rec = add_phigaro_prophage(rec, pg)
    
    with open(fout, "w") as fh:
        SeqIO.write(records, fh, "genbank")



if __name__ == '__main__':
    main()