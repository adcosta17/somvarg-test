import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
from intervaltree import Interval, IntervalTree


def get_read_insert_position(cigar, ref_start, min_size=100):
    ref_count = ref_start
    read_count = 0
    read_insert_positions = []
    for cg in re.findall('[0-9]*[A-Z=]', cigar):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
            ref_count += int(cg[:cg.find("M")])
        if cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
            ref_count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            read_start = read_count
            read_count += int(cg[:cg.find("I")])
            read_end = read_count
            if abs(read_end - read_start) >= min_size:
                read_insert_positions.append([ref_count, read_start,read_end])
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('P'):
            read_count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
    return read_insert_positions

def get_updated_insert_pos(offset_per_pos, pos, insert_pos):
    nearby = offset_per_pos[insert_pos:insert_pos+1]
    for item in nearby:
        return insert_pos + item.data
    return insert_pos


parser = argparse.ArgumentParser( description='Parse Minigraph alignments of reads to an augmented graph. ')
parser.add_argument('--control-bam', required=True)
parser.add_argument('--test-bam', required=True)
parser.add_argument('--tsv', required=True)
args = parser.parse_args()

# Read in the tsv of reads and the insertions they are supposed to support
# TSV based on read alignments to the assembly contigs, and alignments of assembly contigs to reference
test_data = defaultdict(list)
test_positions = defaultdict(int)
control_data = defaultdict(list)
control_positions = defaultdict(int)
missing_positions = defaultdict(int)
missing_count = 0
with open(args.tsv, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        # Each row:
        # Read  Chromosome    Pos    
        row = line.strip().split("\t")
        if int(row[5]) - int(row[4]) < 100:
            # Need to have at least 100 bp on the read to support the insertion
            continue
        if "Test" in row[0]:
            #print(row)
            test_data[row[1]].append((row[2], int(row[3])))
            #print(row[1]+"\t"+str(control_data[row[1]]))
            test_positions[row[2]+"_"+row[3]] += 1
        else:
            control_data[row[1]].append((row[2], int(row[3])))
            control_positions[row[2]+"_"+row[3]] += 1


chrs_to_use = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

# Identify if we have a TP, a path alignment to an SV within 500bp of the expected position
control_tp = {}
control_missed_fn = {}
control_seen_fn = {}
control_seen_fn_elsewhere = {}
samfile = pysam.AlignmentFile(args.control_bam)
for chrom in chrs_to_use:
    print(chrom)
    for record in samfile.fetch(contig=chrom):
        if record.query_name in control_data and record.mapq > 0:
            # Parse insertion and get insert positions 
            inserts = get_read_insert_position(record.cigarstring, record.reference_start)
            for data in control_data[record.query_name]:
                for insert in inserts:
                    if abs(data[1] - insert[0]) < 500:
                        # Matching TP
                        control_tp[record.query_name+"_"+data[0]+"_"+str(data[1])] = 1

for chrom in chrs_to_use:
    print(chrom)
    for record in samfile.fetch(contig=chrom):
        if record.query_name in control_data and record.mapq > 0:
            for data in control_data[record.query_name]:
                if record.query_name+"_"+data[0]+"_"+str(data[1]) in control_tp:
                    continue
                if record.reference_start < data[1] and record.reference_end > data[1]:
                    control_seen_fn[record.query_name+"_"+data[0]+"_"+str(data[1])] = 1

for chrom in chrs_to_use:
    print(chrom)
    for record in samfile.fetch(contig=chrom):
        if record.query_name in control_data and record.mapq > 0:
            for data in control_data[record.query_name]:
                if record.query_name+"_"+data[0]+"_"+str(data[1]) in control_tp or record.query_name+"_"+data[0]+"_"+str(data[1]) in control_seen_fn:
                    continue
                control_seen_fn_elsewhere[record.query_name+"_"+data[0]+"_"+str(data[1])] = 1

for read in control_data:
    for data in control_data[read]:
        if read+"_"+data[0]+"_"+str(data[1]) in control_tp or read+"_"+data[0]+"_"+str(data[1]) in control_seen_fn or read+"_"+data[0]+"_"+str(data[1]) in control_seen_fn_elsewhere:
            continue
        control_missed_fn[read+"_"+data[0]+"_"+str(data[1])] = 1


print("Control:")
print(len(control_tp))
print(len(control_seen_fn))
print(len(control_seen_fn_elsewhere))
print(len(control_missed_fn))


# Identify if we have a TP, a path alignment to an SV within 500bp of the expected position
test_tp = {}
test_missed_fn = {}
test_seen_fn = {}
test_seen_fn_elsewhere = {}
samfile = pysam.AlignmentFile(args.test_bam)
for chrom in chrs_to_use:
    for record in samfile.fetch(contig=chrom):
        if record.query_name in test_data and record.mapq > 0:
            # Parse insertion and get insert positions 
            inserts = get_read_insert_position(record.cigarstring, record.reference_start)
            for data in test_data[record.query_name]:
                for insert in inserts:
                    if abs(data[1] - insert[0]) < 500:
                        # Matching TP
                        test_tp[record.query_name+"_"+data[0]+"_"+str(data[1])] = 1

for chrom in chrs_to_use:
    for record in samfile.fetch(contig=chrom):
        if record.query_name in test_data and record.mapq > 0:
            for data in test_data[record.query_name]:
                if record.query_name+"_"+data[0]+"_"+str(data[1]) in test_tp:
                    continue
                if record.reference_start < data[1] and record.reference_end > data[1]:
                    test_seen_fn[record.query_name+"_"+data[0]+"_"+str(data[1])] = 1

for chrom in chrs_to_use:
    for record in samfile.fetch(contig=chrom):
        if record.query_name in test_data and record.mapq > 0:
            for data in test_data[record.query_name]:
                if record.query_name+"_"+data[0]+"_"+str(data[1]) in test_tp or record.query_name+"_"+data[0]+"_"+str(data[1]) in test_seen_fn:
                    continue
                test_seen_fn_elsewhere[record.query_name+"_"+data[0]+"_"+str(data[1])] = 1

for read in test_data:
    for data in test_data[read]:
        if read+"_"+data[0]+"_"+str(data[1]) in test_tp or read+"_"+data[0]+"_"+str(data[1]) in test_seen_fn or read+"_"+data[0]+"_"+str(data[1]) in test_seen_fn_elsewhere:
            continue
        test_missed_fn[read+"_"+data[0]+"_"+str(data[1])] = 1


print("Test:")
print(len(test_tp))
print(len(test_seen_fn))
print(len(test_seen_fn_elsewhere))
print(len(test_missed_fn))

