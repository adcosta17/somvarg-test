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


parser = argparse.ArgumentParser( description='Rename fastq')
parser.add_argument('--input', required=True)
parser.add_argument('--output-fastq', required=True)
args = parser.parse_args()

insert_dict = defaultdict(list)
max_count = 0

with open(args.output_fastq, 'w') as fout, pysam.FastxFile(args.input) as in_fq:
    for entry in in_fq:
        name = get_name(max_count)
        fout.write("@"+name+"\n")
        fout.write(entry.sequence+"\n")
        fout.write("+\n")
        fout.write("="*len(entry.sequence)+"\n")
        max_count += 1




