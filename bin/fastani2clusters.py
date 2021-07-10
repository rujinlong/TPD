#!/usr/bin/env python

import click
import re
import pandas as pd

def clustering(lsts):
    """https://github.com/rikpg/IntersectionMerge/blob/master/core.py"""
    sets = [set(lst) for lst in lsts if lst]
    merged = 1
    while merged:
        merged = 0
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = 1
                    common |= x
            results.append(common)
        sets = results
    return sets


@click.command()
@click.option("--fastani", '-i', help="FastANI output file")
@click.option("--similarity", '-s', type=int, help="Minimum similarity in a cluster")
@click.option("--fout", '-o', help="output file name")
def main(fastani, similarity, fout):
    edges = pd.read_csv(fastani, sep='\t', usecols=range(3), names=['source', 'target', 'similarity'])
    edges_selected = edges[edges.similarity >= similarity]
    nodes_pairs = edges_selected[['source', 'target']].values
    nodes_pairs = [list(x) for x in nodes_pairs]
    
    # create clusters
    clusters = clustering(nodes_pairs)
    
    # add cluster id
    tmp = []
    for cluster_id in range(len(clusters)):
        for phage in clusters[cluster_id]:
            tmp.append([cluster_id, phage])
    df_clusters = pd.DataFrame(tmp, columns=['clusterID', 'genome'])
    
    # add genome length to table
    df_clusters['genome_length'] = df_clusters.apply(lambda x:int(re.search('len(.+?)_', x['genome']).group(1)), axis=1)
    
    # select longest genome
    replen = df_clusters.groupby('clusterID').max('genome_length').reset_index()
    df_rep = df_clusters.merge(replen, on=['clusterID', 'genome_length'], how='inner')
    df_rep.drop_duplicates(['clusterID', 'genome_length'], inplace=True)
    df_rep.reset_index(drop=True, inplace=True)
    df_rep['representative'] = "yes"
    
    # add rep genome to cluster table
    tbl = df_clusters.merge(df_rep[['genome', 'representative']], on=['genome'], how='left')
    tbl.fillna("no", inplace=True)
    
    tbl.to_csv(fout, sep='\t', index=False)
    


if __name__ == '__main__':
    main()