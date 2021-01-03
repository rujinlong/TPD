#!/usr/bin/env python

import click
import pandas as pd
from Bio import SeqIO

def update_others(genome_region, df, anno_db):
    for feature in genome_region.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers.get('locus_tag')
            if locus_tag:
                locus_tag = locus_tag[0]
            df2 = df[df.protid==locus_tag]
            if df2.shape[0]>0:
                feature.qualifiers[anno_db] = list(df2.anno.unique())
    return genome_region

def reformat_abricate(df_abricate):
    arg_anno = df_abricate.to_dict(orient="records")
    
    rst = {}
    for anno in arg_anno:
        locus_tag = anno['SEQUENCE']
        anno = [anno['DATABASE'], anno['ACCESSION'], anno['PRODUCT']]
        if not rst.get(locus_tag):
            rst[locus_tag] = [anno]
        if anno not in rst.get(locus_tag):
            rst[locus_tag].append(anno)
            
    return rst


def update_abricate(genome_region, annotation_dict):
    for feature in genome_region.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers.get('locus_tag')
            if locus_tag:
                locus_tag = locus_tag[0]
                annotations = annotation_dict.get(locus_tag)
                if annotations:
                    for anno in annotations:
                        feature.qualifiers["annotation_abricate_" + anno[0]] = str(anno[1]) + "~~" + str(anno[2])
    return genome_region


@click.command()
@click.option("--gbk", '-g', help="Genbank file")
@click.option("--annotation", '-a', help="Annotation file")
@click.option("--db", '-d', help="VOGDB|pVOG|PFAM|KEGG|abricate")
@click.option("--fout", '-o', default='out.txt', help="Output file name")
def main(gbk, annotation, db, fout):
    records = list(SeqIO.parse(gbk, format="genbank"))
    anno_db = "annotation_" + db

    df = pd.read_csv(annotation, sep='\t')
    
    if df.shape[0] > 0:
        if db == "abricate":
            annotation_dict = reformat_abricate(df)
            for genome_region in records:
                genome_region = update_abricate(genome_region, annotation_dict)
        else:  # "VOGDB|pVOG|PFAM|KEGG"
            clms = list(df.columns)
            clms.remove('protid')
            df['anno'] = df.apply(lambda x:'~~'.join(x[clms]), axis=1)
            for genome_region in records:
                genome_region = update_others(genome_region, df, anno_db)

    with open(fout, "w") as fh:
        SeqIO.write(records, fh, "genbank")


if __name__ == '__main__':
    main()
