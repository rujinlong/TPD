#!/usr/bin/env python

# Usage: pretty_cdhit.py -i cluster.clstr -o cluster.tsv

import pandas as pd
import click
import re

@click.command()
@click.option("--fn_clstr", '-i', help="cd-hit clstr output")
@click.option("--out", '-o', help="virsorter.tsv")
def pretty_cdhit(fn_clstr, out):
    with open(fn_clstr, 'r') as fh:
        f = fh.read().splitlines()

    p = re.compile(r'>(.*)\.\.\.')
    entries_all = []
    entries_rep = []
    for row in f:
        if row.startswith('>'):
            clst = row
        elif row[-1] == '*':
            ctgid = p.findall(row)[0]
            rep_seq = ctgid
            entries_rep.append([clst, rep_seq])
            entries_all.append([clst, ctgid])
        else:
            ctgid = p.findall(row)[0]
            entries_all.append([clst, ctgid])
            
    df_rep = pd.DataFrame(entries_rep, columns=['clstr', 'repid'])
    df_all = pd.DataFrame(entries_all, columns=['clstr', 'ctgid'])
    df_final = df_all.merge(df_rep, on='clstr', how='left')
    df_final[['repid', 'ctgid']].to_csv(out, index=False, sep='\t')
    

if __name__ == '__main__':
    pretty_cdhit()
