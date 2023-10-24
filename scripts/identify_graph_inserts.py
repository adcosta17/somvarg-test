import argparse
import sys
import os
import csv
import re
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from multiprocessing import Lock, Queue
import threading
import time
import queue
import pysam

def get_family(annotation):
    if "L1" in annotation or "LINE" in annotation:
        return "LINE"
    elif "Alu" in annotation:
        return "Alu"
    elif "SVA" in annotation:
        return "SVA"
    else:
        return "Ambiguous"


def get_neb(avg_size, read_lengths, min_read_len, read_positions):
    total = 0
    for read in read_lengths:
        if read in read_positions:
            total += max(0, read_lengths[read]-1000-avg_size)
    return total/1000000000

def mapped_end_to_end_hap(cigarstring):
    # Count up the position on the read until we get to a deletion
    count = 0
    read_start = 0
    read_end = 0
    match_count = 0
    large_indel = []
    clipped_bases = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
            match_count += int(cg[:cg.find("M")])
        if cg.endswith('='):
            count += int(cg[:cg.find("=")])
            match_count += int(cg[:cg.find("=")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
            match_count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
            if int(cg[:cg.find("I")]) > 50:
                large_indel.append([count, "I", int(cg[:cg.find("I")])])
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) > 50:
                large_indel.append([count, "D", int(cg[:cg.find("D")])])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            if match_count > 0 and read_end == 0:
                read_end = count
            count += int(cg[:cg.find("S")])
            clipped_bases += int(cg[:cg.find("S")])
            if match_count == 0:
                read_start = count
        elif cg.endswith('H'):
            if match_count > 0 and read_end == 0:
                read_end = count
            count += int(cg[:cg.find("H")])
            clipped_bases += int(cg[:cg.find("H")])
            if match_count == 0:
                read_start = count
    mapped = True
    if len(large_indel) > 0 or clipped_bases/count > 0.1:
        mapped = False
    ret = []
    pos = read_start
    for i in range(len(large_indel)):
        ret.append([pos, large_indel[i][0]])
        if large_indel[i][1] == "I":
            pos = large_indel[i][2]+large_indel[i][0]
        else:
            pos = large_indel[i][0]
    ret.append([pos,read_end])
    return (mapped, ret)
    

def mapped_end_to_end(cigarstring, read_len, read_aln_start, read_aln_end):
    # Count up the position on the read until we get to a deletion
    count = 0
    read_start = read_aln_start
    read_end = read_aln_end
    match_count = 0
    large_indel = []
    clipped_bases = read_aln_start + (read_len - read_aln_end)
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
            match_count += int(cg[:cg.find("M")])
        if cg.endswith('='):
            count += int(cg[:cg.find("=")])
            match_count += int(cg[:cg.find("=")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
            match_count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
            if int(cg[:cg.find("I")]) > 50:
                large_indel.append([count, "I", int(cg[:cg.find("I")])])
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) > 50:
                large_indel.append([count, "D", int(cg[:cg.find("D")])])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
    mapped = True
    if len(large_indel) > 0 or clipped_bases/read_len > 0.1:
        mapped = False
    ret = []
    pos = read_start
    for i in range(len(large_indel)):
        ret.append([pos, large_indel[i][0]])
        if large_indel[i][1] == "I":
            pos = large_indel[i][2]+large_indel[i][0]
        else:
            pos = large_indel[i][0]
    ret.append([pos,read_end])
    return (mapped, ret)

parser = argparse.ArgumentParser( description='Parse Minigraph alignments of reads to an augmented graph. Identify Translocations')
parser.add_argument('--gaf-folder', required=True)
parser.add_argument('--gfa', required=True)
parser.add_argument('--tsv-folder', required=True)
parser.add_argument('--max-distance', type=int, default=250)
parser.add_argument('--samples', required=True)
args = parser.parse_args()
contig_nodes = {}

node_lengths = defaultdict(int)
with open(args.gfa, 'r') as in_graph:
    for line in in_graph:
        row = line.strip().split('\t')
        if row[0] == "S":
            # Have a segment line
            node_lengths[row[1]] = len(row[2])
            chrom = row[4].split(':')[2]
            pos = int(row[5].split(':')[2])
            contig_nodes[row[0]] = chrom+":"+row[2]+"_"#+row[3]
            length = node_lengths[row[0]]


chrs_to_use = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
reads_to_use = {}

insert_count = 0
for sample in args.samples.split(','):
    with open(args.tsv_folder+"/"+sample+"_inserts.tsv" ,'r') as in_tsv:
        count = 0
        for line in in_tsv:
            row = line.strip().split('\t')
            if count == 0:
                row.append("AssemblyAligned\tSamples")
                count = 1
                #print("Sample\t"+'\t'.join(row))
                continue
            #if (int(row[13])-int(row[12]))/int(len(row[4])) < 0.5:
            #    continue
            if row[7] == "PASS" and "Polymorphic" not in line:
                insert_count += 1
                reads = row[6].split(',')
                read_samples = {}
                for read in reads:
                    read_info = read.split(':')
                    read_name = read_info[1]
                    reads_to_use[read_name] = 1


read_positions = defaultdict(list)
read_lengths = defaultdict(int)
polymorphic_reads = {}
for sample in args.samples.split(','):
    with open(sample+"/"+args.gaf_folder+"/"+sample+".gaf", 'r') as in_alignment:
        for line in in_alignment:
            row = line.strip().split("\t")
            read_lengths[row[0]] = int(row[1])
            if row[0] not in reads_to_use:
                continue
            mapped, positions = mapped_end_to_end(row[18].split(':')[2], int(row[1]), int(row[2]), int(row[3]))
            if mapped:
                polymorphic_reads[row[0]] = 1
            read_positions[row[0]].extend(positions)


for sample in args.samples.split(','):
    with open(args.tsv_folder+"/"+sample+"_inserts.tsv", 'r') as in_tsv, open(args.tsv_folder+"/"+sample+"_updated_inserts.tsv", 'w') as out_tsv:
        count = 0
        polymorphic_count = 0
        novel_count = 0
        na_count = 0
        for line in in_tsv:
            row = line.strip().split('\t')
            if count == 0:
                row.append("AssemblyAligned\tSamples")
                count = 1
                #print("Sample\t"+'\t'.join(row))
                continue
            if "PASS" in line and "Assembly" in line:
                polymorphic = "PossiblyNovel_NotInAssembly"
                reads = row[6].split(',')
                read_samples = {}
                found_count = 0
                #print(line, file=sys.stderr)
                for read in reads:
                    read_info = read.split(':')
                    read_name = read_info[1]
                    orientation = read_info[2]
                    read_insert_start = int(read_info[3].split('-')[0])
                    read_insert_end = int(read_info[3].split('-')[1])
                    if orientation == '-':
                        tmp = read_insert_start
                        read_insert_start = read_lengths[read_name] - read_insert_end
                        read_insert_end = read_lengths[read_name] - read_insert_start
                    read_samples[read_info[0]] = 1
                    if "Polymorphic" in polymorphic:
                        found_count += 1
                        continue
                    if read_name in polymorphic_reads:
                        polymorphic = "AssemblyPolymorphic_E2E"
                        found_count += 1
                    elif read_name in read_positions:
                        found_count += 1
                        for position in read_positions[read_name]:
                            if (read_insert_start - position[0]) > args.max_distance and  (position[1] - read_insert_end) > args.max_distance:
                                polymorphic = "AssemblyPolymorphic_Insert"
                if found_count == 0:
                    polymorphic = "NotInAssembly"
                row[len(row)-2] = polymorphic
                out_tsv.write("\t".join(row)+"\n")
        print(sample+"\t"+str(polymorphic_count)+"\t"+str(novel_count)+"\t"+str(na_count))

