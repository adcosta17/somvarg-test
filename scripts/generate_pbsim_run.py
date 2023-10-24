import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import random
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser( description='Add insertion sequences to contigs')
parser.add_argument('--input', required=True)
parser.add_argument('--output-folder', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-tsv', required=True)
parser.add_argument('--pbsim-model', required=True)
parser.add_argument('--pbsim-path', required=True)
args = parser.parse_args()

with open(args.input, 'r') as in_tsv:
    with open(args.output_tsv, 'w') as out_tsv:
        for line in in_tsv:
            row = line.strip().split('\t')
            fa = args.output_folder+"/"+args.output_prefix+"."+str(row[0])+".fa"
            print(args.pbsim_path+" "+fa+" --prefix "+args.output_folder+"/"+args.output_prefix+"."+row[0]+".pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 100 --hmm_model "+args.pbsim_model+" --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+row[0]+".pbsim_0001.maf")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+row[0]+".pbsim_0001.ref")
            out_tsv.write(row[0]+"\t"+args.output_folder+"/"+args.output_prefix+"."+row[0]+".pbsim_0001.fastq\t"+fa+"\n")

