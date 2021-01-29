#!/usr/bin/env python

import click


@click.command()
@click.option("--fin", '-i', help="input file name")
@click.option("--fout", '-o', default='out.txt', help="output file name")
@click.option("--pred_method", '-m', default='phigaro', help="Prophage prediction method")
def main(fin, fout, pred_method):
    with open(fin, "r") as fh:
        f = fh.read().splitlines()
        prophages = [x.split('\t') for x in f]
    
    rst = []
    for i in range(len(prophages)):
        new_rec = prophages[i] + ["{}_{}".format(pred_method, str(i+1))]
        rst.append('\t'.join(new_rec))

    with open(fout, "w") as fh:
        fh.write("\n".join(rst))


if __name__ == '__main__':
    main()