#!/usr/bin/env python

import click
from Bio import SeqIO
import pandas as pd
import re

def match_region_on_ref(query_rec, ref_rec):
    query_new = str(query_rec.seq).replace('N', '')
    refseq_old = str(ref_rec.seq)
    refseq_new = refseq_old.replace('N', '')
    match = re.search(query_new, refseq_new)
    
    if match:
        ref_idx_map = {}
        new_idx = 0
        for i in range(len(refseq_old)):
            if refseq_old[i] != "N":
                ref_idx_map[new_idx] = i
                new_idx += 1

        match_start_in_ref_old_idx = ref_idx_map[match.start()]
        match_end_in_ref_old_idx = ref_idx_map[match.end()]
        return [match_start_in_ref_old_idx, match_end_in_ref_old_idx, query_rec.id]


def extract_prophage_seq_and_id(contig_rec, prophage_type_name="prophage"):
    prophage_regions = [x for x in contig_rec.features if x.type==prophage_type_name]
    prophages = []
    for rg in prophage_regions:
        prophage_rec = contig_rec[rg.location.start:rg.location.end]
        prophage_rec.id = rg.qualifiers['ID'][0]
        prophages.append(prophage_rec)
    return prophages


def map_old_prophage_to_new_ref(gbk_old, gbk_new):
    old_gbk = list(SeqIO.parse(gbk_old, "genbank"))
    new_gbk = list(SeqIO.parse(gbk_new, "genbank"))
    
    rst = []
    for rec_old in old_gbk:
        prophages = extract_prophage_seq_and_id(rec_old)
        for prophage_rec in prophages:
            for ref in new_gbk:
                idx = match_region_on_ref(prophage_rec, ref)
                if idx != None:
                    rst.append([ref.id] + idx)
    return rst




@click.command()
@click.option("--old_gbk", '-i', help="old genbank file")
@click.option("--new_gbk", '-n', help="new genbank file")
@click.option("--fout", '-o', help="old_prophage.tsv")
def main(old_gbk, new_gbk, fout):
    mapping = map_old_prophage_to_new_ref(old_gbk, new_gbk)
    df = pd.DataFrame(mapping)
    df.to_csv(fout, index=False, header=False, sep='\t')


if __name__ == '__main__':
    main()