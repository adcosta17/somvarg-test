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
        positions.append(r[0]+"."+r[1].split(':')[1])
    return ",".join(positions)


parser = argparse.ArgumentParser( description='Extract SVs from a vcf or bedpe and apply them to the reference to get a set of alternative haplotypes to augment a graph with')
#parser.add_argument('--cute-vcf', required=False)
#parser.add_argument('--sniffles-vcf', required=False)
parser.add_argument('--tsv', required=True)
parser.add_argument('--flank', default=100000, type=int, required=False)
parser.add_argument('--ref', required=True)
parser.add_argument('--output-folder', required=True)
parser.add_argument('--min-insertion-length', type=int, default=50)
args = parser.parse_args()


ref_file = pysam.FastaFile(args.ref)

current_positions = defaultdict(int)
svs_to_write = defaultdict(list)
prev_contig = "chr10"

contig_lengths = defaultdict(int)

# Get shared SV positions between Sniffles and CuteSV
#cuteSV_positions = defaultdict(IntervalTree)
#bcf_in = pysam.VariantFile(args.cute_vcf)
#contigs = bcf_in.header.contigs.items()
#for item in contigs:
#    contig_lengths[item[0]] = int(item[1].length)
#for rec in bcf_in.fetch():
#    if rec.pos < args.flank:
#        continue 
#    if rec.info["SVTYPE"] == "INS":
#        alt_seq = rec.alts[0].upper()
#        if len(alt_seq) < args.min_insertion_length:
#            continue
#        if rec.pos <= args.flank or contig_lengths[rec.contig] - rec.pos <= args.flank:
#            continue
#        cuteSV_positions[rec.contig][rec.pos-50:rec.pos+50] = 1
#
#sniffles_positions = defaultdict(IntervalTree)
#bcf_in = pysam.VariantFile(args.sniffles_vcf)
#for rec in bcf_in.fetch():
#    if rec.pos < args.flank:
#        continue 
#    if rec.info["SVTYPE"] == "INS":
#        alt_seq = rec.alts[0].upper()
#        if len(alt_seq) < args.min_insertion_length:
#            continue
#        if rec.pos <= args.flank or contig_lengths[rec.contig] - rec.pos <= args.flank:
#            continue
#        sniffles_positions[rec.contig][rec.pos-50:rec.pos+50] = 1

entries = defaultdict(IntervalTree)
with open(args.tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        contig = row[0]
        window_pos = int(row[1])
        positions = get_position_string(row[2])
        alt_seq = row[3]
        min_pos, max_pos = get_min_max_pos(row[2])
        entries[contig][min_pos:max_pos+1] = ">"+contig+"_"+str(window_pos)+"_"+str(min_pos)+"_"+str(max_pos)+"_"+positions+"_"+str(len(positions))+"\n"+alt_seq+"\n"
        print(contig+"_"+str(window_pos)+"_"+str(min_pos)+"_"+str(max_pos)+"_"+positions)

for contig in entries:
    current_positions = defaultdict(int)
    for item in entries[contig]:
        start = item.begin
        end = item.end
        i = 0
        while True:
            if (start - args.flank) - current_positions[i] > args.flank or current_positions[i] == 0:
                # No overlaps here, add to list
                current_positions[i] = end + args.flank
                svs_to_write[i].append(item.data)
                break
            else:
                i += 1
                #print(str(i) + " " + str(rec.pos - args.flank) + " " + str(current_positions[i]) , file=sys.stderr)

for i in svs_to_write:
    with open(args.output_folder+"/SVS_"+str(i)+".fa", 'w') as out_fa:
        for item in svs_to_write[i]:
            out_fa.write(item)
