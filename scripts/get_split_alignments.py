import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re
import random

parser = argparse.ArgumentParser( description='Get split alignment numbers')
parser.add_argument('--bam', required=True)
parser.add_argument('--mat', required=True)
parser.add_argument('--pat', required=True)
args = parser.parse_args()


split_count_bam = defaultdict(int)
split_count_mat = defaultdict(int)
split_count_pat = defaultdict(int)
split_count_combined = defaultdict(int)
total_mapped_bam = {}
total_mapped_mat = {}
total_mapped_pat = {}
total_mapped_combined = 0

count = 0
sam_reader = pysam.AlignmentFile(args.bam)
print("BAM")
for record in sam_reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count)
    if record.mapping_quality == 60:
        split_count_bam[record.query_name] += 1
        total_mapped_bam[record.query_name] = 1

total_split_bam = 0
for read in split_count_bam:
    if split_count_bam[read] >= 2:
        total_split_bam +=1

sam_reader = pysam.AlignmentFile(args.mat)
print("MAT")
count = 0
for record in sam_reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count)
    if record.mapping_quality == 60:
        split_count_mat[record.query_name] += 1
        total_mapped_mat[record.query_name] = 1

total_split_mat = 0
for read in split_count_mat:
    if split_count_mat[read] >= 2:
        total_split_mat +=1

sam_reader = pysam.AlignmentFile(args.pat)
print("PAT")
count = 0
for record in sam_reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count)
    if record.mapping_quality == 60:
        split_count_pat[record.query_name] += 1
        total_mapped_pat[record.query_name] = 1

total_split_pat = 0
for read in split_count_pat:
    if split_count_pat[read] >= 2:
        total_split_pat +=1

total_combined = 0
for read in split_count_mat:
    if read in split_count_pat:
        total_mapped_combined += 1
        if split_count_mat[read] == 1 or split_count_pat[read] == 1:
            continue
        total_combined += 1
    else:
        total_mapped_combined += 1
        if split_count_mat[read] == 1:
            continue
        total_combined += 1
for read in split_count_pat:
    if read not in split_count_mat:
        total_mapped_combined += 1
        if split_count_pat[read] == 1:
            continue
        total_combined += 1

print("GRCh38:\t"+str(total_split_bam)+"\t"+str(len(total_mapped_bam))+"\t"+str(total_split_bam/len(total_mapped_bam)))
print("Mat:\t"+str(total_split_mat)+"\t"+str(len(total_mapped_mat))+"\t"+str(total_split_mat/len(total_mapped_mat)))
print("Pat:\t"+str(total_split_pat)+"\t"+str(len(total_mapped_pat))+"\t"+str(total_split_pat/len(total_mapped_pat)))
print("Combined:\t"+str(total_combined)+"\t"+str(total_mapped_combined)+"\t"+str((total_combined)/(total_mapped_combined)))

