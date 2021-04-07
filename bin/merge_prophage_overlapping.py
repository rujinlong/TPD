#!/usr/bin/env python

import click
import pandas as pd


def read_prophage_intervals(fname):
    df = pd.read_csv(fname, sep='\t', names=['sampleID', 'nstart', 'nend', 'prophage_id'])
    intervals = {x:df[df.sampleID==x][['nstart', 'nend']].values.tolist() for x in df.sampleID.unique()}
    return intervals


def merge_overlapping(intervals):
    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged


def prophage_regions(merged_intervals):
    df = pd.DataFrame(merged_intervals, columns=['sampleID', 'nstart', 'nend'])
    df['prophage_id'] = ['merged_prophage_{}'.format(str(x+1)) for x in df.index.to_list()]
    return df[['sampleID', 'nstart', 'nend', 'prophage_id']]


@click.command()
@click.option("--fin", '-i', help="Input file name")
@click.option("--fout", '-o', default='out.txt', help="Output file name")
def main(fin, fout):
    all_intervals = read_prophage_intervals(fin)
    all_merged = []
    for sample_id in all_intervals.keys():
        intervals = all_intervals[sample_id]
        merged = merge_overlapping(intervals)
        all_merged += [[sample_id] + x for x in merged]
    prophages = prophage_regions(all_merged)
    prophages.to_csv(fout, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()