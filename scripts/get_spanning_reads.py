import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import random
from intervaltree import Interval, IntervalTree
import mappy as mp

def get_name(n):
    tmp = hex(n)[2:].zfill(32)
    return tmp[0:8]+"-"+tmp[8:12]+"-"+tmp[12:16]+"-"+tmp[16:20]+"-"+tmp[20:]

parser = argparse.ArgumentParser( description='Get reads that span the insertion sequence')
parser.add_argument('--input', required=True)
parser.add_argument('--inserts', required=True)
parser.add_argument('--tsv', required=True)
parser.add_argument('--fastq', required=True)
parser.add_argument('--flank', type=int, default=2000)
args = parser.parse_args()

# Open up the input tsv
# Map the reads in the fastq to the ref for each fastq/ref pair
# If the read spans the expected insertion position then output it to the fastq
# Output the info about the read spanning to the tsv

insert_dict = {}
with open(args.inserts, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        insert_dict[row[0]] = row

seen = {}
with open(args.input, 'r') as in_tsv:
    with open(args.fastq, 'w') as fout, open(args.tsv, 'w') as out_tsv:
        count = 1
        for line in in_tsv:
            row = line.strip().split('\t')
            # Map the reads for this sequence to the target
            a = mp.Aligner(row[2])  # load or build index
            if not a: raise Exception("ERROR: failed to load/build index")
            for name, seq, qual in mp.fastx_read(row[1]): # read a fasta/q sequence
                for hit in a.map(seq): # traverse alignments
                    # Get the start and end positions and compare them to the expected insert positions. If overlap exists then write the read to the output
                    start = int(insert_dict[row[0]][3])
                    end = start + len(insert_dict[row[0]][8]) + len(insert_dict[row[0]][10])
                    #print(str(start - len(insert_dict[row[0]][10]) - 500)+"\t"+str(end + len(insert_dict[row[0]][10]) + 500)+"\t"+str(hit))
                    if hit.mapq == 60 and hit.r_st < (start - len(insert_dict[row[0]][10]) - args.flank) and hit.r_en > (end + len(insert_dict[row[0]][10]) + args.flank) and row[0]+"_"+name not in seen:
                        # Read spans insert with sufficient flank
                        seen[row[0]+"_"+name] = 1
                        n = get_name(count)
                        fout.write("@"+n+"\n")
                        fout.write(seq+"\n")
                        fout.write("+\n")
                        fout.write(qual+"\n")
                        out_tsv.write(row[0]+"\t"+name+"\t"+n+"\t"+str(count)+"\n")
                        count += 1