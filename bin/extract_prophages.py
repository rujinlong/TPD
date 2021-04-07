#!/usr/bin/env python

import click
from Bio import SeqIO


def extract_prophage(rec, prophage_feature, flank_len):
    start = prophage_feature.location.start - flank_len
    end = prophage_feature.location.end + flank_len
    if start <= 1:
        start = 1
    if end > len(rec):
        end = len(rec)
    length = abs(start-end)
    prophage = rec[start-1:end]
    if not prophage.id.startswith(rec.id):
        # prophage.id = prophage_feature.qualifiers['ID'][0]
        prophage.id = "{}_{}_len{}_{}-{}".format(prophage.id, prophage_feature.qualifiers['ID'][0], str(length), start, end)
    else:
        prophage.id = "{}_{}_len{}_{}-{}".format(prophage.id, prophage_feature.qualifiers['ID'][0], str(length), start, end)
    prophage.name = prophage.id
    prophage.type = "prophage"
    print(prophage.id)
    return prophage


@click.command()
@click.option("--gbk", '-g', help="Input file name")
@click.option("--fout_prefix", '-o', help="Prefix of output file name")
@click.option('--prophage_type_name', '-p', default="prophage", help="prophage type name in genbank file. [prophage | mics_feature]")
@click.option('--output_type', '-t', default='merged', help="Output all prophages or only merged [merged | total]")
@click.option("--flank_length", '-l', default=1000, help="Extend sides of each prophage region with this length.")
def main(gbk, fout_prefix, prophage_type_name, output_type, flank_length):
    gbk_recs = list(SeqIO.parse(gbk, format="genbank"))
    prophages = []
    for rec in gbk_recs:
        # if not rec.id.startswith(sample_id):
        #     rec.id = "{}_{}".format(sample_id, rec.id)
            # rec.name = rec.id
        prophage_features = [x for x in rec.features if x.type==prophage_type_name]
        for pf in prophage_features:
            try:
                pf.qualifiers.get('ID')
                prophage = extract_prophage(rec, pf, flank_length)
                prophages.append(prophage)
            except:
                continue

    # Biopython 1.78 requires 'molecule_type' key in annotations
    for pp in prophages:
        if not pp.annotations.get('molecule_type'):
            pp.annotations['molecule_type'] = 'DNA'
    
    
    # Extract only merged prophages
    if output_type == 'merged':
        prophages = [x for x in prophages if 'merged_prophage' in x.id]
    
    for pp in prophages:
        print(pp.id)

    with open("{}.gbk".format(fout_prefix), "w") as fh:
        SeqIO.write(prophages, fh, "genbank")
    
    with open("{}.fasta".format(fout_prefix), "w") as fh:
        SeqIO.write(prophages, fh, "fasta")
    

if __name__ == '__main__':
    main()