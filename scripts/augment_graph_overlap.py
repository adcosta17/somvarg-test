import pysam
import argparse
import sys
import csv
from collections import defaultdict
from intervaltree import Interval, IntervalTree

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def get_rc(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))

def get_min_max_pos(name):
    positions = []
    for item in name.split(','):
        r = item.split('_')
        pos = int(r[1].split(':')[1])
        positions.append(pos)
    return min(positions), max(positions)


def get_position_string(name):
    positions = []
    for item in name.split(','):
        r = item.split('_')
        positions.append([int(r[0]), int(r[1]), int(r[2].split(':')[1])])
    return positions


parser = argparse.ArgumentParser( description='Extract SVs from a vcf or bedpe and apply them to the reference to get a set of alternative haplotypes to augment a graph with')
parser.add_argument('--tsv', required=True)
parser.add_argument('--overlap', default=51, type=int, required=False)
parser.add_argument('--ref', required=True)
parser.add_argument('--chrom', required=False, default="")
parser.add_argument('--nodes', required=True)
args = parser.parse_args()

overlap_size = args.overlap

ref_file = pysam.FastaFile(args.ref)

current_positions = defaultdict(int)
svs_to_write = defaultdict(list)
ref_seqs = {}


entries = defaultdict(IntervalTree)
with open(args.tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        contig = row[0]
        if args.chrom != "":
            if args.chrom != contig:
                continue
        window_pos = int(row[1])
        positions = get_position_string(row[2])
        alt_seq = row[3]
        if contig not in entries:
            seq = ref_file.fetch(contig)
            entries[contig].add(Interval(0, len(seq))) 
            ref_seqs[contig] = seq
        for pos in positions:
            nearby = entries[contig][pos[2]:pos[2]+1]
            # Split up the node and replace it
            for item in nearby:
                entries[contig].remove(Interval(item.begin, item.end))
                # Check to see if pos[2] is equal to either the item.begin or end, ie we have 2 inserts at the same position
                if item.begin != pos[2]:
                    entries[contig].add(Interval(item.begin, pos[2]))
                if item.end != pos[2]:
                    entries[contig].add(Interval(pos[2], item.end))

# Assign each node an ID
count = 0
nodes = {}
node_sequences = {}
node_start_end = {}
for contig in entries:
    nodes[contig] = {}
    node_sequences[contig] = {}
    node_start_end[contig] = {}
    for item in entries[contig]:
        nodes[contig][item.begin] = "s_"+contig+"_"+str(count)
        start = item.begin - overlap_size
        end = item.end + overlap_size
        if start < 0: 
            start = 0
        if end > len(ref_seqs[contig]):
            end = len(ref_seqs[contig])
        node_sequences[contig]["s_"+contig+"_"+str(count)] = ref_seqs[contig][start:end]
        node_start_end[contig]["s_"+contig+"_"+str(count)] = [start, end]
        count += 1

edges = {}
for contig in entries:
    for item in entries[contig]:
        start = item.begin - 1
        if start < 0 :
            start = 0
        end = item.end + 1
        nearby = entries[contig][start:end]
        for iv in nearby:
            if iv.begin == item.begin and iv.end == item.end:
                continue
            if iv.begin == item.end:
                edges["L\t"+nodes[contig][item.begin]+"\t+\t"+nodes[contig][iv.begin]+"\t+\t"+str(overlap_size)+"M"] = 1

# Once we have nodes computed, Parse the tsv again and output the insert seqs with the matching nodes
nodes["Inserts"] = defaultdict(list)
node_sequences["Inserts"] = {}
node_start_end["Inserts"] = {}
with open(args.tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        contig = row[0]
        if args.chrom != "":
            if args.chrom != contig:
                continue
        window_pos = int(row[1])
        positions = get_position_string(row[2])
        alt_seq = row[3]
        for pos in positions:
            nodes["Inserts"][contig+":"+str(pos[2])].append("s_"+contig+"_"+str(count)+"Insert")
            start = pos[2] - overlap_size
            end = pos[2] + overlap_size
            if start < 0:
                start = 0
            if end > len(ref_seqs[contig]):
                end = len(ref_seqs[contig])
            node_sequences["Inserts"]["s_"+contig+"_"+str(count)+"Insert"] = ref_seqs[contig][start:pos[2]]+alt_seq[pos[0]:pos[0]+pos[1]]+ref_seqs[contig][pos[2]:end]
            node_start_end["Inserts"]["s_"+contig+"_"+str(count)+"Insert"] = [start, end]
            # Create a node for the insert and then add a link between it and the contig nodes it matches with 
            nearby = entries[contig][pos[2]-1:pos[2]+1]
            for item in nearby:
                if item.begin == pos[2]:
                    # Link between insert and node
                    edges["L\ts_"+contig+"_"+str(count)+"Insert"+"\t+\t"+nodes[contig][item.begin]+"\t+\t"+str(overlap_size)+"M"] = 1
                elif item.end == pos[2]:
                    # Link between node and insert
                    edges["L\t"+nodes[contig][item.begin]+"\t+\ts_"+contig+"_"+str(count)+"Insert"+"\t+\t"+str(overlap_size)+"M"] = 1
            count += 1

print("H\tVN:Z:1.0")
# Print out the nodes and edges
with open(args.nodes, 'w') as out_nodes:
    for contig in nodes:
        if "Inserts" == contig:
            for pos in nodes[contig]:
                chrom = pos.split(':')[0]
                position = pos.split(":")[1]
                for item in nodes[contig][pos]:
                    print("S\t"+item+"\t"+node_sequences[contig][item]+"\tLN:i:"+str(len(node_sequences[contig][item])))
                    out_nodes.write(item+"\t"+chrom+"\t"+position+"\n")
        else:
            for pos in nodes[contig]:
                print("S\t"+nodes[contig][pos]+"\t"+node_sequences[contig][nodes[contig][pos]]+"\tLN:i:"+str(len(node_sequences[contig][nodes[contig][pos]])))
                out_nodes.write(nodes[contig][pos]+"\t"+contig+"\t"+str(pos)+"\n")

for item in edges:
    print(item)
