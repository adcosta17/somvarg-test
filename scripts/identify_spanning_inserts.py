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
parser.add_argument('--gaf', required=True)
parser.add_argument('--gfa', required=True)
parser.add_argument('--tsv', required=True)
parser.add_argument('--hap1', type=str, required=True)
parser.add_argument('--hap2', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=500)
parser.add_argument('--normalized', type=str, required=True)
parser.add_argument('--normalized-hap', type=str, required=True)
parser.add_argument('--sample', required=True)
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
with open(args.tsv, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        row = line.strip().split('\t')
        if count == 0:
            row.append("AssemblyAligned\tSamples")
            count = 1
            #print("Sample\t"+'\t'.join(row))
            continue
        if args.sample not in row[5]:
            continue
        #if (int(row[13])-int(row[12]))/int(len(row[4])) < 0.5:
        #    continue
        if row[6] == "PASS" and "Polymorphic" not in line:
            insert_count += 1
            reads = row[5].split(',')
            read_samples = {}
            for read in reads:
                read_info = read.split(':')
                read_name = read_info[1]
                reads_to_use[read_name] = 1


read_positions = defaultdict(list)
read_lengths = defaultdict(int)
polymorphic_reads = {}
with open(args.gaf, 'r') as in_alignment:
    for line in in_alignment:
        row = line.strip().split("\t")
        read_lengths[row[0]] = int(row[1])
        if row[0] not in reads_to_use or read_lengths[row[0]] < 7000:
            continue
        mapped, positions = mapped_end_to_end(row[18].split(':')[2], int(row[1]), int(row[2]), int(row[3]))
        if mapped:
            polymorphic_reads[row[0]] = 1
        read_positions[row[0]].extend(positions)

read_positions_hap = defaultdict(list)
polymorphic_reads_hap = {}
reader = pysam.AlignmentFile(args.hap1)
for record in reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count, flush=True, file=sys.stderr)
    if record.query_name not in reads_to_use or read_lengths[record.query_name] < 7000:
        continue
    mapped, positions = mapped_end_to_end_hap(record.cigarstring)
    if mapped:
        polymorphic_reads_hap[record.query_name] = 1
    read_positions_hap[record.query_name].extend(positions)

print("HAP1 Input", file=sys.stderr)

reader = pysam.AlignmentFile(args.hap2)
for record in reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count, flush=True, file=sys.stderr)
    if record.query_name not in reads_to_use or read_lengths[record.query_name] < 7000:
        continue
    mapped, positions = mapped_end_to_end_hap(record.cigarstring)
    if mapped:
        polymorphic_reads_hap[record.query_name] = 1
    read_positions_hap[record.query_name].extend(positions)

print("HAP2 Input", file=sys.stderr)

counts = {}
single_sample_counts = {}
sizes = defaultdict(list)
counts["LINE"] = 0
counts["All"] = 0
counts["Alu"] = 0
counts["SVA"] = 0
counts["Ambiguous"] = 0
single_sample_counts["LINE"] = 0
single_sample_counts["All"] = 0
single_sample_counts["Alu"] = 0
single_sample_counts["SVA"] = 0
single_sample_counts["Ambiguous"] = 0

with open(args.tsv, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        row = line.strip().split('\t')
        if count == 0:
            row.append("AssemblyAligned\tSamples")
            count = 1
            #print("Sample\t"+'\t'.join(row))
            continue
        if args.sample not in row[5]:
            continue
        #if (int(row[13])-int(row[12]))/int(len(row[4])) < 0.5:
        #    continue
        if row[6] == "PASS" and "Polymorphic" not in line:
            polymorphic = "PossiblyNovel_NotInAssembly"
            reads = row[5].split(',')
            read_samples = {}
            for read in reads:
                read_info = read.split(':')
                read_name = read_info[1]
                read_samples[read_info[0]] = 1
                if "Polymorphic" in polymorphic:
                    continue
                if read_name in polymorphic_reads:
                    polymorphic = "AssemblyPolymorphic_E2E"
                elif read_name in read_positions:
                    read_insert_start = int(read_info[3].split('-')[0])
                    read_insert_end = int(read_info[3].split('-')[1])
                    for position in read_positions[read_name]:
                        if (read_insert_start - position[0]) > args.max_distance and  (position[1] - read_insert_end) > args.max_distance:
                            polymorphic = "AssemblyPolymorphic_Insert"
                else:
                    polymorphic = "NotInAssembly"
            row.append(polymorphic)
            single_sample = False
            if len(read_samples) == 1:
                row.append("SingleSample")
                single_sample = True
            else:
                row.append("MultiSample")
            if polymorphic == "PossiblyNovel_NotInAssembly":
                family = get_family(row[7])
                # Have a passing insert
                sizes[family].append(len(row[4]))
                sizes["All"].append(len(row[4]))
                counts[family] += 1
                counts["All"] += 1
                if single_sample:
                    single_sample_counts[family] += 1
            print("Graph\t"+"\t".join(row))
        else:
            row.append("NA\tNA")

avg_size = {} 
avg_size["All"] = 0
if counts["All"] > 0:
    avg_size["All"] = sum(sizes["All"])/len(sizes["All"])
avg_size["LINE"] = 0
if counts["LINE"] > 0:
    avg_size["LINE"] = sum(sizes["LINE"])/len(sizes["LINE"])
avg_size["Alu"] = 0
if counts["Alu"] > 0:
    avg_size["Alu"] = sum(sizes["Alu"])/len(sizes["Alu"])
avg_size["SVA"] = 0
if counts["SVA"] > 0:
    avg_size["SVA"] = sum(sizes["SVA"])/len(sizes["SVA"])
avg_size["Ambiguous"] = 0
if counts["Ambiguous"] > 0:
    avg_size["Ambiguous"] = sum(sizes["Ambiguous"])/len(sizes["Ambiguous"])


with open(args.normalized, 'w') as out_norm:
    out_norm.write("Sample\tAll_Count\tAll_Single_Sample\tLINE_Count\tLINE_Single_Sample\tAlu_Count\tAlu_Single_Sample\tSVA_Count\tSVA_Single_Sample\tAmbiguous_Count\tAmbiguous_Single_Sample\n")
    out_norm.write(args.sample)
    for item in avg_size:
        out_norm.write("\t"+str(counts[item])+"\t"+str(single_sample_counts[item]))
    out_norm.write("\n")
    out_norm.write(str(insert_count)+"\t"+str(len(read_positions))+'\n')


counts = {}
single_sample_counts = {}
sizes = defaultdict(list)
counts["LINE"] = 0
counts["All"] = 0
counts["Alu"] = 0
counts["SVA"] = 0
counts["Ambiguous"] = 0
single_sample_counts["LINE"] = 0
single_sample_counts["All"] = 0
single_sample_counts["Alu"] = 0
single_sample_counts["SVA"] = 0
single_sample_counts["Ambiguous"] = 0

with open(args.tsv, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        row = line.strip().split('\t')
        if count == 0:
            row.append("AssemblyAligned\tSamples")
            count = 1
            #print("Sample\t"+'\t'.join(row))
            continue
        if args.sample not in row[5]:
            continue
        #if (int(row[13])-int(row[12]))/int(len(row[4])) < 0.5:
        #    continue
        if row[6] == "PASS" and "Polymorphic" not in line:
            polymorphic = "PossiblyNovel_NotInAssembly"
            reads = row[5].split(',')
            read_samples = {}
            for read in reads:
                read_info = read.split(':')
                read_name = read_info[1]
                read_samples[read_info[0]] = 1
                if "Polymorphic" in polymorphic:
                    continue
                if read_name in polymorphic_reads_hap:
                    polymorphic = "AssemblyPolymorphic_E2E"
                elif read_name in read_positions_hap:
                    read_insert_start = int(read_info[3].split('-')[0])
                    read_insert_end = int(read_info[3].split('-')[1])
                    for position in read_positions_hap[read_name]:
                        if (read_insert_start - position[0]) > args.max_distance and  (position[1] - read_insert_end) > args.max_distance:
                            polymorphic = "AssemblyPolymorphic_Insert"
                else:
                    polymorphic = "NotInAssembly"
            row.append(polymorphic)
            single_sample = False
            if len(read_samples) == 1:
                row.append("SingleSample")
                single_sample = True
            else:
                row.append("MultiSample")
            if polymorphic == "PossiblyNovel_NotInAssembly":
                family = get_family(row[7])
                # Have a passing insert
                sizes[family].append(len(row[4]))
                sizes["All"].append(len(row[4]))
                counts[family] += 1
                counts["All"] += 1
                if single_sample:
                    single_sample_counts[family] += 1
            print("Hap\t"+"\t".join(row))
        else:
            row.append("NA\tNA")

avg_size = {} 
avg_size["All"] = 0
if counts["All"] > 0:
    avg_size["All"] = sum(sizes["All"])/len(sizes["All"])
avg_size["LINE"] = 0
if counts["LINE"] > 0:
    avg_size["LINE"] = sum(sizes["LINE"])/len(sizes["LINE"])
avg_size["Alu"] = 0
if counts["Alu"] > 0:
    avg_size["Alu"] = sum(sizes["Alu"])/len(sizes["Alu"])
avg_size["SVA"] = 0
if counts["SVA"] > 0:
    avg_size["SVA"] = sum(sizes["SVA"])/len(sizes["SVA"])
avg_size["Ambiguous"] = 0
if counts["Ambiguous"] > 0:
    avg_size["Ambiguous"] = sum(sizes["Ambiguous"])/len(sizes["Ambiguous"])


with open(args.normalized_hap, 'w') as out_norm:
    out_norm.write("Sample\tAll_Count\tAll_Single_Sample\tLINE_Count\tLINE_Single_Sample\tAlu_Count\tAlu_Single_Sample\tSVA_Count\tSVA_Single_Sample\tAmbiguous_Count\tAmbiguous_Single_Sample\n")
    out_norm.write(args.sample)
    for item in avg_size:
        out_norm.write("\t"+str(counts[item])+"\t"+str(single_sample_counts[item]))
    out_norm.write("\n")
    out_norm.write(str(insert_count)+"\t"+str(len(read_positions_hap))+'\n')
