import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict
import threading
import copy
import pickle

def get_family(annotation):
    if "LINE" in annotation:
        return "LINE"
    elif "SINE" in annotation:
        return "SINE"
    elif "ERV" in annotation:
        return "ERV"
    elif "SVA" in annotation:
        return "SVA"
    else:
        return "ambiguous"

parser = argparse.ArgumentParser( description='See which inserts are shared between samples')
parser.add_argument('--sample', required=True)
parser.add_argument('--suffix', default=".insertions.repbase_annotated.tsv")
parser.add_argument('--folder', default="assembly_analysis")
parser.add_argument('--max-distance', default=1000)
args = parser.parse_args()


found_in_samples = defaultdict(IntervalTree)
output_positions = defaultdict(str)

for sample in args.sample.split(','):
    sample_tsv = "./"+sample+"/"+args.folder+"/"+sample+args.suffix
    interval_ret = {}
    print(sample_tsv,file=sys.stderr)
    insert_ret = defaultdict(IntervalTree)
    with open(sample_tsv) as csvfile:
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            #if "PASS" not in row_args[8]:
            #    continue
            chrom = row_args[0]
            start = int(row_args[1]) - args.max_distance
            end = int(row_args[2]) + args.max_distance
            family = get_family(row_args[9])
            #if "ambiguous" in family:
            #    continue
            found_in_samples[chrom][start:end] = sample+"_"+family

for sample in args.sample.split(','):
    sample_tsv = "./"+sample+"/"+args.folder+"/"+sample+args.suffix
    interval_ret = {}
    print(sample_tsv,file=sys.stderr)
    insert_ret = defaultdict(IntervalTree)
    with open(sample_tsv) as csvfile:
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            #if "PASS" not in row_args[8]:
            #    continue
            family = get_family(row_args[9])
            #if "ambiguous" in family:
            #    continue
            chrom = row_args[0]
            start = int(row_args[1])
            end = int(row_args[2])
            nearby = found_in_samples[chrom][start:end]
            seen = {}
            if chrom+":"+row_args[1]+"-"+row_args[2] in output_positions:
                for item in output_positions[chrom+":"+row_args[1]+"-"+row_args[2]].split(','):
                    seen[item] = 1
            for item in nearby:
                seen[item.data] = 1
            ret = []
            for item in seen:
                ret.append(item)
            output_positions[chrom+":"+row_args[1]+"-"+row_args[2]] = ",".join(ret)

for item in output_positions:
    annotation = "FAIL"
    if get_family(output_positions[item]) != "ambiguous":
        annotation = "PASS"
    print(item+"\t"+output_positions[item]+"\t"+annotation)
