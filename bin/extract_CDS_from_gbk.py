#!/usr/bin/env python

import click
from Bio import SeqIO, SeqRecord, Seq


def extract_CDSs(gbks, fout_DNA, fout_AA):
    cds_DNAs = []
    cds_AAs = []
    for gbk in gbks:
        cds_features = [x for x in gbk.features if x.type=='CDS']
        for cds_feature in cds_features:
            locus_tag = cds_feature.qualifiers.get('locus_tag')[0]
            
            protein_seq = cds_feature.qualifiers.get('translation')
            if protein_seq:
                protein_seq = protein_seq[0]
                cds_AA = SeqRecord.SeqRecord(seq=Seq.Seq(protein_seq), id=locus_tag, description="")
                cds_AAs.append(cds_AA)
            
                cds_DNA =  cds_feature.location.extract(gbk)
                cds_DNA.id=locus_tag
                cds_DNA.description = ""
                cds_DNAs.append(cds_DNA)
        
    with open(fout_DNA, "w") as fh:
        SeqIO.write(cds_DNAs, fh, format="fasta")
    
    with open(fout_AA, "w") as fh:
        SeqIO.write(cds_AAs, fh, format="fasta")


def extract_genome(gbks, fout):
    with open(fout, "w") as fh:
        SeqIO.write(gbks, fh, format="fasta")


@click.command()
@click.option("--fgbk", '-i', help="Genbank file")
@click.option("--fout_prefix", '-o', help="Genbank file")
def main(fgbk, fout_prefix):
    gbks = list(SeqIO.parse(fgbk, format="genbank"))

    # Output
    fout_genome = "{}_genome.fasta".format(fout_prefix)
    fout_cds_DNA = "{}_gene.fna".format(fout_prefix)
    fout_cds_AA = "{}_protein.faa".format(fout_prefix)

    extract_genome(gbks, fout_genome)
    extract_CDSs(gbks, fout_cds_DNA, fout_cds_AA)


if __name__ == '__main__':
    main()