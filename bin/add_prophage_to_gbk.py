#!/usr/bin/env python

import click
from Bio import SeqFeature, SeqIO
import pandas as pd


def add_prophage(record, prophage_loc, note):
    start_pos = SeqFeature.ExactPosition(prophage_loc['start'])
    end_pos = SeqFeature.ExactPosition(prophage_loc['end'])
    feature_location = SeqFeature.FeatureLocation(start_pos-1, end_pos, strand=1)
    qualifiers = {"ID": prophage_loc['prophage_id'], "note": note}

    new_feature = SeqFeature.SeqFeature(feature_location, type="misc_feature", qualifiers=qualifiers)
    record.features.append(new_feature)
    return record


@click.command()
@click.option("--gbk", '-g', help="input file")
@click.option("--prophage_coords", '-p', help="Prophage coordinates file")
@click.option("--note", '-n', default="Prophage", required=False, help="Prophage|Contig")
@click.option("--fout", '-o', help="output file name")
def main(gbk, prophage_coords, fout, note):
    records = list(SeqIO.parse(gbk, format="genbank"))
    prophages = pd.read_csv(prophage_coords, sep='\t', names = ['ref_id', 'start', 'end', 'prophage_id']).to_dict(orient="records")
    
    for rec in records:
        print(rec.id)
        for pg_loc in prophages:
            if rec.id == pg_loc['ref_id']:
                rec = add_prophage(rec, pg_loc, note)
    
    with open(fout, "w") as fh:
        SeqIO.write(records, fh, "genbank")


if __name__ == '__main__':
    main()