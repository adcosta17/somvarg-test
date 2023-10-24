import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re
import random

def get_name(n):
    tmp = hex(n)[2:].zfill(32)
    return tmp[0:8]+"-"+tmp[8:12]+"-"+tmp[12:16]+"-"+tmp[16:20]+"-"+tmp[20:]


parser = argparse.ArgumentParser( description='Spike in insertion supporting reads into fastq')
parser.add_argument('--tsv', required=True)
parser.add_argument('--inserts', required=True)
parser.add_argument('--mat', required=True)
parser.add_argument('--pat', required=True)
parser.add_argument('--output-fastq', required=True)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--position-fraction', type=float, default=0.25)
parser.add_argument('--position-number', type=int, default=500)
parser.add_argument('--max-reads', type=int, default=4)
args = parser.parse_args()

random.seed(args.rep)
insert_dict = defaultdict(list)
max_count = 0

with pysam.FastaFile(args.inserts) as fh, open(args.output_fastq, 'w') as fout, pysam.FastxFile(args.mat) as mat_fq, pysam.FastxFile(args.pat) as pat_fq, open(args.tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split("\t")
        insert_dict[row[0]].append(row[2])
        max_count = max(max_count, int(row[3]))
    max_count += 1
    for entry in mat_fq:
        name = get_name(max_count)
        fout.write("@"+name+"\n")
        fout.write(entry.sequence+"\n")
        fout.write("+\n")
        fout.write("="*len(entry.sequence)+"\n")
        max_count += 1
    for entry in pat_fq:
        name = get_name(max_count)
        fout.write("@"+name+"\n")
        fout.write(entry.sequence+"\n")
        fout.write("+\n")
        fout.write("="*len(entry.sequence)+"\n")
        max_count += 1
    if args.rep > 0:
        # 0 is a base case, don't add any insertions for 0 only otherwise 
        # Now randomly select insertions to insert
        for pos in insert_dict:
            n = random.random()
            if n < args.position_fraction:
                # Use this position, select n reads
                max_reads = min(len(insert_dict[pos]), args.max_reads)
                reads = random.sample(insert_dict[pos], random.randint(1,max_reads))
                for read in reads:
                    seq = fh.fetch(read)
                    fout.write("@"+read+"\n")
                    fout.write(seq+"\n")
                    fout.write("+\n")
                    fout.write("="*len(seq)+"\n")
                    print(pos+"\t"+read)



