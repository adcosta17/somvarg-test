import pysam
import argparse
import sys
import os
import csv
import re
#import pyabpoa as pa
import mappy as mp
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from multiprocessing import Lock, Queue
import threading
import time
import queue

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def get_seen_nearby(chrom, start, end, seen):
    new_seen = {}
    new_seen = set()
    i = start
    count = 0
    if chrom not in seen:
        return new_seen
    while i < end:
        if i in seen[chrom]:
            new_seen = new_seen.union(seen[chrom][i])
            count += 1
        i += 1
    return new_seen

def get_bp_pos(alignment_i, i_read_pos, node_lengths, contig_nodes):
    i_pos = ""
    path = alignment_i[5].replace('>',' ').replace('<', ' ').split(' ')[1:]
    item = path[0]
    chrom = ""
    pos = ""
    #if alignment_i[5][0] == '>':
        # Forward
    if i_read_pos == int(alignment_i[2]):
        # Get the start chrom position
        item = path[0]
        length = node_lengths[item]
        #if "Insert" in item:
        #    length = 0
        #    if len(path) > 1:
        #        item = path[1]
        chrom = contig_nodes[item].split(':')[0]
        pos = str(int(contig_nodes[item].split(':')[1].split('_')[0])+int(alignment_i[7]))
    elif i_read_pos == int(alignment_i[3]):
        item = path[len(path) - 1]
        length = node_lengths[item]
        #if "Insert" in item:
        #    length = 0
        #    if len(path) > 1:
        #        item = path[len(path) - 2]
        #        length = node_lengths[item]
        chrom = contig_nodes[item].split(':')[0]
        pos = str(int(contig_nodes[item].split(':')[1].split('_')[0])+length-(int(alignment_i[6])-int(alignment_i[8])))
#    else:
#        # Reverse alignment
#        if i_read_pos == int(alignment_i[3]):
#            # Get the start chrom position
#            item = path[0]
#            length = node_lengths[item]
#            #if "Insert" in item:
#            #    length = 0
#            #    if len(path) > 1:
#            #        item = path[1]
#            #        length = node_lengths[item]
#            chrom = contig_nodes[item].split(':')[0]
#            pos = str(int(contig_nodes[item].split(':')[1].split('_')[0])+int(alignment_i[7]))
#        elif i_read_pos == int(alignment_i[2]):
#            item = path[len(path) - 1]
#            length = node_lengths[item]
#            #if "Insert" in item:
#            #    length = 0
#            #    if len(path) > 1:
#            #        item = path[len(path) - 2]
#            #        length = node_lengths[item]
#            chrom = contig_nodes[item].split(':')[0]
#            pos = str(int(contig_nodes[item].split(':')[1].split('_')[0])+length-(int(alignment_i[6])-int(alignment_i[8])))
    return chrom+":"+pos

class myThread2(threading.Thread):
    def __init__(self, i, q, t_lock, positions, read_seqs, nearby_read_alignments, multiple_aligned_reads, read_path_to_positions, aligner):
        threading.Thread.__init__(self)
        self.i = i
        self.q = q
        self.t_lock = t_lock
        self.positions = positions
        self.read_seqs = read_seqs
        self.nearby_read_alignments = nearby_read_alignments
        self.multiple_aligned_reads = multiple_aligned_reads
        self.read_path_to_positions = read_path_to_positions
        self.aligner = aligner
    def run(self):
        print("Running "+str(self.i), file=sys.stderr)
        while(True):
            if self.q.empty():
                print("Queue is Empty", file=sys.stderr)
                return
            try:
                position = self.q.get(True, 100)
                #print(position)
            except queue.Empty:
                print("Timeout Queue is Empty", file=sys.stderr)
                return
            chrom1 = position.split('-')[0].split(':')[0]
            pos1 = int(position.split('-')[0].split(':')[1])
            chrom2 = position.split('-')[1].split(':')[0]
            pos2 = int(position.split('-')[1].split(':')[1])
            read_seqs = []
            read_list = ""
            chrom1_len = 0
            chrom1_read = ""
            chrom2_len = 0
            chrom2_read = ""
            chrom1_read_name = ""
            chrom2_read_name = ""
            for read in self.positions[position]:
                seq = self.read_seqs[read]
                read_seqs.append(seq)
                read_list += ","+read
                for i in range(len(self.nearby_read_alignments[read])):
                    i_pos = self.nearby_read_alignments[read][i][0]
                    j_pos = self.nearby_read_alignments[read][i][1]
                    i_read_pos = self.nearby_read_alignments[read][i][4]
                    j_read_pos = self.nearby_read_alignments[read][i][5]
                    alignment_i = list(self.multiple_aligned_reads[read].values())[i_pos]
                    alignment_j = list(self.multiple_aligned_reads[read].values())[j_pos]
                    read_len = int(alignment_i[1])
                    read_i_start = int(alignment_i[2])
                    read_i_end = int(alignment_i[3])
                    ref_i_chrom = self.read_path_to_positions[read][alignment_i[5]][0]
                    ref_i_start = self.read_path_to_positions[read][alignment_i[5]][1]
                    ref_i_end = self.read_path_to_positions[read][alignment_i[5]][2]
                    read_j_start = int(alignment_j[2])
                    read_j_end = int(alignment_j[3])
                    ref_j_chrom = self.read_path_to_positions[read][alignment_j[5]][0]
                    ref_j_start = self.read_path_to_positions[read][alignment_j[5]][1]
                    ref_j_end = self.read_path_to_positions[read][alignment_j[5]][2]
                    if ref_i_chrom == chrom1 and (abs(ref_i_end - pos1) < 100  or abs(ref_i_start - pos1) < 100):
                        if read_i_start == i_read_pos:
                            # start to end
                            if (read_len - read_i_start) > chrom1_len:
                                chrom1_len = read_len - read_i_start
                                chrom1_read = seq
                                chrom1_read_name = read
                        else:
                            if read_i_end > chrom1_len:
                                chrom1_len = read_i_end
                                chrom1_read = seq
                                chrom1_read_name = read
                    if ref_i_chrom == chrom2 and (abs(ref_i_end - pos2) < 100  or abs(ref_i_start - pos2) < 100):
                        if read_i_start == i_read_pos:
                            # start to end
                            if (read_len - read_i_start) > chrom2_len:
                                chrom2_len = read_len - read_i_start
                                chrom2_read = seq
                                chrom2_read_name = read
                        else:
                            if read_i_end > chrom2_len:
                                chrom2_len = read_i_end
                                chrom2_read = seq
                                chrom2_read_name = read
                    if ref_j_chrom == chrom1 and (abs(ref_j_end - pos1) < 100  or abs(ref_j_start - pos1) < 100):
                        if read_j_start == j_read_pos:
                            # start to end
                            if (read_len - read_j_start) > chrom1_len:
                                chrom1_len = read_len - read_j_start
                                chrom1_read = seq
                                chrom1_read_name = read
                        else:
                            if read_j_end > chrom1_len:
                                chrom1_len = read_j_end
                                chrom1_read = seq
                                chrom1_read_name = read
                    if ref_j_chrom == chrom2 and (abs(ref_j_end - pos2) < 100  or abs(ref_j_start - pos2) < 100):
                        if read_j_start == j_read_pos:
                            # start to end
                            if (read_len - read_j_start) > chrom2_len:
                                chrom2_len = read_len - read_j_start
                                chrom2_read = seq
                                chrom2_read_name = read
                        else:
                            if read_j_end > chrom2_len:
                                chrom2_len = read_j_end
                                chrom2_read = seq
                                chrom2_read_name = read
            sorted_read_seqs = sorted(read_seqs, key=len)
            if chrom1_read != chrom2_read:
                # Compute a target sequence to use for the concencus calling
                a = mp.Aligner(seq=chrom1_read)
                target = ""
                for hit in a.map(chrom2_read):
                    #print(hit, flush=True)
                    if hit.q_st < hit.r_st:
                        # read1 is first then read2
                        target = chrom1_read[0:hit.r_st]
                        if hit.strand == -1:
                            target += reverse_complement(chrom2_read)[hit.q_st:]
                        else:
                            target += chrom2_read[hit.q_st:]
                    else:
                        # read2 then read 1
                        target = chrom2_read[0:hit.q_st]
                        if hit.strand == -1:
                            target += reverse_complement(chrom1_read)[hit.r_st:]
                        else:
                            target += chrom1_read[hit.r_st:]
                if target == "":
                    self.t_lock.acquire()
                    print(str(self.i)+"\tTarget Empty\t"+position+"\t"+chrom1_read_name+"\t"+str(chrom1_len)+"\t"+chrom2_read_name+"\t"+str(chrom2_len), file=sys.stderr)
                    self.t_lock.release()
                #print("Target\t"+target, flush=True)
                sorted_read_seqs.insert(0,target)
            else:
                sorted_read_seqs.insert(0, chrom1_read)
            con_seq = ""
#            if self.aligner == "abpoa":
#                # Compute concencus with abPOA
#                a = pa.msa_aligner()
#                res=a.msa(sorted_read_seqs, out_cons=True, out_msa=False)
#                for seq in res.cons_seq:
#                    con_seq = seq
#                    break
            if self.aligner == "racon":
                # Use Racon, Target is the first seq in sorted_read_seqs
                with open("target_"+str(self.i)+".fa", 'w') as out_fa:
                    out_fa.write(">target\n"+sorted_read_seqs[0]+"\n")
                with open("sequences_"+str(self.i)+".fa", 'w') as out_fa:
                    for i in range(len(sorted_read_seqs)):
                        if i == 0:
                            continue
                        out_fa.write(">seq_"+str(i)+"\n"+sorted_read_seqs[i]+"\n")
                # Compute the overlaps file
                hit_count = 0
                with open("overlaps_"+str(self.i)+".paf", 'w') as out_paf:
                    a = mp.Aligner("target_"+str(self.i)+".fa") 
                    if not a: raise Exception("ERROR: failed to load/build index")
                    for name, seq, qual in mp.fastx_read("sequences_"+str(self.i)+".fa"): # read a fasta/q sequence
                        for hit in a.map(seq): # traverse alignments
                            out_paf.write(name+"\t"+str(len(seq))+"\t"+str(hit)+"\n")
                            hit_count += 1
                # Run racon
                #print("Running Racon", flush=True)
                if hit_count > 0:
                    os.system(args.racon+" sequences_"+str(self.i)+".fa overlaps_"+str(self.i)+".paf target_"+str(self.i)+".fa > concencus_"+str(self.i)+".fa")
                else:
                    os.system("rm sequences_"+str(self.i)+".fa overlaps_"+str(self.i)+".paf target_"+str(self.i)+".fa")
                    self.t_lock.acquire()
                    print(str(self.i)+"\tNo overlap found:\t"+position+"\t"+str(len(sorted_read_seqs))+"\t"+str(len(target)), file=sys.stderr)
                    print(chrom1+"\t"+str(pos1)+"\t"+chrom2+"\t"+str(pos2)+"\t"+str(len(sorted_read_seqs) - 1)+"\tFalse\tNA\tNA")
                    self.t_lock.release()
                    continue
                #print("Finished Racon", flush=True)
                # Read in concecncus
                for name, seq, qual in mp.fastx_read("concencus_"+str(self.i)+".fa"):
                    con_seq = seq
                    break
                os.system("rm sequences_"+str(self.i)+".fa overlaps_"+str(self.i)+".paf target_"+str(self.i)+".fa concencus_"+str(self.i)+".fa")
            else:
                print("Error: Must select one of [abpoa, racon] as the aligner for concencus calling")
            #print(len(con_seq))
            # Align the seq to the reference and see where it maps to 
            a = mp.Aligner(args.ref)
            found1 = "False"
            found2 = "False"
            ret = "Expected\t"+chrom1+":"+str(pos1)+"\t"+chrom2+":"+str(pos2)+"\n"
            for hit in a.map(con_seq):
                #print(hit)
                if hit.ctg == chrom1 and (abs(hit.r_st-pos1) < 50 or abs(hit.r_en-pos1) < 50):
                    found1 = "True"
                elif hit.ctg == chrom2 and (abs(hit.r_st-pos2) < 50 or abs(hit.r_en-pos2) < 50):
                    found2 = "True"
                ret += "\t"+hit.ctg+":"+str(hit.r_st)+"-"+str(hit.r_en)+"\n"
            #if found1 and found2:
            #    print(chrom1+"\t"+str(pos1)+"\t"+chrom2+"\t"+str(pos2)+"\t"+read_list[1:])
            self.t_lock.acquire()
            print(ret, file=sys.stderr)
            print(chrom1+"\t"+str(pos1)+"\t"+chrom2+"\t"+str(pos2)+"\t"+str(len(sorted_read_seqs) - 1)+"\tTrue\t"+found1+"\t"+found2)
            self.t_lock.release()


parser = argparse.ArgumentParser( description='Parse Minigraph alignments of reads to an augmented graph. Identify Translocations')
parser.add_argument('--gaf', required=True)
parser.add_argument('--gfa', required=True)
parser.add_argument('--min-reads', required=False, default=3, type=int)
parser.add_argument('--read-gap', required=False, default=100, type=int)
parser.add_argument('--nearby-window', required=False, default=50, type=int)
parser.add_argument('--read-fraction', required=False, default=0.85, type=float)
parser.add_argument('--same-chrom', action='store_true')
parser.add_argument('--nodes', required=True)
parser.add_argument('--fastq', required=True)
parser.add_argument('--ref', required=True)
parser.add_argument('--aligner', required=False, default="racon")
parser.add_argument('--racon', required=False, default="racon")
parser.add_argument('--threads', required=False, default=1, type=int)

args = parser.parse_args()

gfa_positions = defaultdict(IntervalTree)
contig_nodes = {}

node_lengths = defaultdict(int)
with open(args.gfa, 'r') as in_graph:
    for line in in_graph:
        row = line.strip().split('\t')
        if row[0] == "S":
            # Have a segment line
            node_lengths[row[1]] = len(row[2])

with open(args.nodes, 'r') as in_graph:
    for line in in_graph:
        row = line.strip().split('\t')
        chrom = row[1].replace("alt_", '')
        pos = int(row[2])
        contig_nodes[row[0]] = chrom+":"+row[2]+"_"#+row[3]
        length = node_lengths[row[0]]
        if "Insert" in row[0]:
            gfa_positions[chrom][pos-500:pos+500] = row[0]



chrs_to_use = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

# In order to have a translocation call we need a split alignment between two regions of the graph
# First identify which reads have multiple alignments
alignment_count = defaultdict(int)
multiple_aligned_reads = {}
with open(args.gaf, 'r') as in_alignment:
    for line in in_alignment:
        row = line.strip().split("\t")
        if int(row[11]) < 60 or int(row[3])-int(row[2])<21:
            continue
        alignment_count[row[0]] += 1

with open(args.gaf, 'r') as in_alignment:
    for line in in_alignment:
        row = line.strip().split("\t")
        if alignment_count[row[0]] >= 2:
            if row[0] not in multiple_aligned_reads:
                multiple_aligned_reads[row[0]] = {}
            multiple_aligned_reads[row[0]][len(multiple_aligned_reads[row[0]])] = row[:12]

print(len(multiple_aligned_reads), file=sys.stderr)
non_overlapping_reads = {}
read_path_to_positions = {}
# Now iterate over reads with multiple alignments
# First filter out reads that have overlapping alignments
for read in multiple_aligned_reads:
    node_positions = defaultdict(IntervalTree)
    read_positions = IntervalTree()
    read_path_to_positions[read] = {}
    for row in list(multiple_aligned_reads[read].values()):
        path = row[5].replace('>',' ').replace('<', ' ').split(' ')[1:]
        start = -1
        end = 0
        length = int(row[8]) - int(row[7])
        for item in path:
            # get the position of the node
            chrom = contig_nodes[item].split(':')[0]
            if chrom not in chrs_to_use:
                continue
            pos = int(contig_nodes[item].split(':')[1].split('_')[0])
            if start < 0:
                start = pos
            if pos < start:
                start = pos
            if "Insert" in item:
                length -= node_lengths[item]
                continue
        start += int(row[7])
        end = start + length
        if end > start:
            node_positions[chrom][start:end] = 1
        read_path_to_positions[read][row[5]] = [chrom, start, end]
        read_positions[int(row[2])+10:int(row[3])-10] = 1
    overlap = False 
    for chrom in node_positions:
        for item in node_positions[chrom]:
            nearby = node_positions[chrom][item.begin:item.end]
            if len(nearby) > 1:
                # Have overlapping alignments on the reference
                overlap = True
    for item in read_positions:
        nearby = read_positions[item.begin:item.end]
        if len(nearby) > 1:
            overlap = True
    if not overlap:
        non_overlapping_reads[read] = True

print(len(non_overlapping_reads), file=sys.stderr)

#if "d7dc142d-f1f2-4f1a-afb5-ea1858e65e6d" in non_overlapping_reads:
#    print("d7dc142d-f1f2-4f1a-afb5-ea1858e65e6d")

# Check each pair of alignments for non overlapping reads
# If the alignments are adjacent to each other on the read and the reference, flag them

nearby_read_alignments = defaultdict(list)
possible_missed_insert_read = defaultdict(list)
location_map = {}
for read in non_overlapping_reads:
    vals = list(multiple_aligned_reads[read].values())
    for i in range(len(vals)):
        alignment_i = vals[i]
        read_len = int(alignment_i[1])
        read_i_start = int(alignment_i[2])
        read_i_end = int(alignment_i[3])
        ref_i_chrom = read_path_to_positions[read][alignment_i[5]][0]
        ref_i_start = read_path_to_positions[read][alignment_i[5]][1]
        ref_i_end = read_path_to_positions[read][alignment_i[5]][2]
        for j in range(len(vals)):
            if i <= j:
                continue
            alignment_j = vals[j]
            read_j_start = int(alignment_j[2])
            read_j_end = int(alignment_j[3])
            ref_j_chrom = read_path_to_positions[read][alignment_j[5]][0]
            ref_j_start = read_path_to_positions[read][alignment_j[5]][1]
            ref_j_end = read_path_to_positions[read][alignment_j[5]][2]
            # Compare the positions of i to j
            read_distance = min(abs(read_i_start - read_j_start), abs(read_i_start - read_j_end), abs(read_i_end - read_j_start), abs(read_i_end - read_j_end))
            if read_distance > args.read_gap:
                # Read positions are too far apart
                continue
            if ref_i_chrom == ref_j_chrom:
                if args.same_chrom:
                    ref_distance = min(abs(ref_i_start - ref_j_start), abs(ref_i_start - ref_j_end), abs(ref_i_end - ref_j_start), abs(ref_i_end - ref_j_end))
                    if ref_distance < read_len:
                        #both alignments are on the same chromsome within the read length of each other
                        continue
                else:
                    continue
            # Check that both alignments together represent and end-to-end or almost end-to-end alignment
            total = abs(read_i_end - read_i_start) + abs(read_j_end - read_j_start)
            #if read == "d7dc142d-f1f2-4f1a-afb5-ea1858e65e6d":
            #    print(read+"\t"+str(total)+"\t"+str(read_len))
            if total/read_len < args.read_fraction:
                continue
            # Have a candidate pair. Get the locations 
            posi = ""
            posj = ""
            read_i_pos = 0
            read_j_pos = 0
            if abs(read_i_start - read_j_end) == read_distance:
                # get chrom position of read_i_start and read_j_end
                posi = get_bp_pos(alignment_i, read_i_start, node_lengths, contig_nodes)
                posj = get_bp_pos(alignment_j, read_j_end, node_lengths, contig_nodes)
                read_i_pos = read_i_start
                read_j_pos = read_j_end
            elif abs(read_j_start - read_i_end) == read_distance:
                # get chrom position of read_i_start and read_j_end
                posi = get_bp_pos(alignment_i, read_i_end, node_lengths, contig_nodes)
                posj = get_bp_pos(alignment_j, read_j_start, node_lengths, contig_nodes)
                read_i_pos = read_i_end
                read_j_pos = read_j_start
            elif abs(read_i_start - read_j_start) == read_distance:
                # get chrom position of read_i_start and read_j_start
                posi = get_bp_pos(alignment_i, read_i_start, node_lengths, contig_nodes)
                posj = get_bp_pos(alignment_j, read_j_start, node_lengths, contig_nodes)
                read_i_pos = read_i_start
                read_j_pos = read_j_start
            elif abs(read_j_end - read_i_end) == read_distance:
                # get chrom position of read_i_end and read_j_end
                posi = get_bp_pos(alignment_i, read_i_end, node_lengths, contig_nodes)
                posj = get_bp_pos(alignment_j, read_j_end, node_lengths, contig_nodes)
                read_i_pos = read_i_end
                read_j_pos = read_j_end
            #if read == "d7dc142d-f1f2-4f1a-afb5-ea1858e65e6d":
            #    print(str(read_i_pos)+"\t"+str(posi)+"\t"+str(read_j_pos)+"\t"+str(posj))
            if ref_i_chrom not in location_map:
                location_map[ref_i_chrom] = defaultdict(set)
            if ref_j_chrom not in location_map:
                location_map[ref_j_chrom] = defaultdict(set)
            location_map[ref_i_chrom][int(posi.split(':')[1])].add(read)
            location_map[ref_j_chrom][int(posj.split(':')[1])].add(read)
            nearby_read_alignments[read].append([i, j, posi, posj, read_i_pos, read_j_pos])

print(len(nearby_read_alignments), file=sys.stderr)

candidate_list = defaultdict(list)
positions = defaultdict(set)
seen = {}
# Group pairs by genomic position, apply min read filter here
for read in nearby_read_alignments:
    # get the breakpoint positions
    for pos in nearby_read_alignments[read]:
        alignment_1 = list(multiple_aligned_reads[read].values())[pos[0]]
        alignment_2 = list(multiple_aligned_reads[read].values())[pos[1]]
        alignmentpos1 = pos[2]
        alignmentpos2 = pos[3]
        chrom1 = alignmentpos1.split(':')[0]
        pos1 = int(alignmentpos1.split(':')[1])
        nearby1 = get_seen_nearby(chrom1, pos1-args.nearby_window, pos1+args.nearby_window, location_map)
        chrom2 = alignmentpos2.split(':')[0]
        pos2 = int(alignmentpos2.split(':')[1])
        nearby2 = get_seen_nearby(chrom2, pos2-args.nearby_window, pos2+args.nearby_window, location_map)
        #if chrom1 == "chr10" and chrom2 == "chr1":
        #    print(read+"\t"+str(pos1)+"\t"+str(pos2))
        #elif chrom2 == "chr10" and chrom1 == "chr1":
        #    print(read+"\t"+str(pos2)+"\t"+str(pos1))
        #print(nearby1)
        #print(nearby2)
        if len(nearby2) < args.min_reads or len(nearby1) < args.min_reads:
            #if read == "d7dc142d-f1f2-4f1a-afb5-ea1858e65e6d":
            #    print("Not enough nearby")
            #    print(nearby1)
            #    print(nearby2)
            continue
        # Get the reads that are on each side and check that there are at least min read number shared
        a1_reads = set()
        a2_reads = set()
        for r in nearby1:
            a1_reads.add(r)
        for r in nearby2:
            a2_reads.add(r)
        combined = a1_reads.intersection(a2_reads)
        if len(combined) < args.min_reads:
            #if read == "d7dc142d-f1f2-4f1a-afb5-ea1858e65e6d":
            #    print("Not enough combined")
            #    print(combined)
            continue
        candidate_list[read].append(pos)
        found = False
        for seen_read in combined:
            if seen_read in seen:
                found = True
            seen[seen_read] = 1
        if found:
            continue
        if chrom1 < chrom2:
            positions[alignmentpos1+"-"+alignmentpos2] = combined
        elif chrom1 == chrom2:
            if pos1 < pos2:
                positions[alignmentpos1+"-"+alignmentpos2] = combined
            else:
                positions[alignmentpos2+"-"+alignmentpos1] = combined
        else:
            positions[alignmentpos2+"-"+alignmentpos1] = combined

print(len(candidate_list), file=sys.stderr)

print("Chrom1\tPos1\tChrom2\tPos2\tSupportingReadCount\tAssembled\tAtExpectedLocation1\tAtExpectedLocation2")
for position in positions:
    chrom1 = position.split('-')[0].split(':')[0]
    pos1 = int(position.split('-')[0].split(':')[1])
    chrom2 = position.split('-')[1].split(':')[0]
    pos2 = int(position.split('-')[1].split(':')[1])
    count = len(positions[position])
    print(chrom1+"\t"+str(pos1)+"\t"+chrom2+"\t"+str(pos2)+"\t"+str(count)+"\tFalse\tNA\tNA")
    #print(positions[position])


exit(0)

# compute a concencsus for the canadidates
fh_in = pysam.FastaFile(args.fastq)
thread_list = [None]*args.threads
t_lock = threading.Lock()
q = queue.Queue()
read_seqs = {}
for position in positions:
    q.put(position)
    for read in positions[position]:
        seq = fh_in.fetch(read)
        read_seqs[read] = seq
for i in range(args.threads):
    print("Submitted "+str(i), file=sys.stderr)
    thread_list[i] = myThread2(i, q, t_lock, positions, read_seqs, nearby_read_alignments, multiple_aligned_reads, read_path_to_positions, args.aligner)
    thread_list[i].start()
for i in range(args.threads):
    thread_list[i].join()
    print("Joined "+str(i), file=sys.stderr)




