import pysam
import argparse
from pysam import VariantFile
import sys
import csv
import re
from collections import defaultdict
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser( description='Parse Minigraph alignments of reads to an augmented graph. Identify Translocations')
parser.add_argument('--control', required=True)
parser.add_argument('--test', required=True)
parser.add_argument('--min-reads', required=False, default=2)
args = parser.parse_args()

break_points_control = defaultdict(IntervalTree)
bcf_in = VariantFile(args.control)
for rec in bcf_in.fetch():
     if rec.info["SVTYPE"] == "BND" and len(rec.info["RNAMES"]) >= int(args.min_reads):
        alt_pos = rec.alts[0].replace('[', '').replace(']', '').replace('N','')
        break_points_control[rec.chrom][rec.pos-50:rec.pos+50] = alt_pos
        alt_pos_chrom = alt_pos.split(':')[0]
        alt_pos_i = int(alt_pos.split(':')[1])
        break_points_control[alt_pos_chrom][alt_pos_i-50:alt_pos_i+50] = rec.chrom+":"+str(rec.pos)


fp_count = 0
bcf_in = VariantFile(args.test)
for rec in bcf_in.fetch():
     if rec.info["SVTYPE"] == "BND" and len(rec.info["RNAMES"]) >= int(args.min_reads):
        alt_pos = rec.alts[0].replace('[', '').replace(']', '').replace('N','')
        nearby = break_points_control[rec.chrom][rec.pos:rec.pos+1]
        if len(nearby) == 0:
            fp_count += 1
        else:
            alt_pos_chrom = alt_pos.split(':')[0]
            alt_pos_i = int(alt_pos.split(':')[1])
            found = False
            for item in nearby:
                nearby_pos = item.data
                nearby_chrom = nearby_pos.split(':')[0]
                nearby_pos_i = int(nearby_pos.split(":")[1])
                if nearby_chrom == alt_pos_chrom and abs(nearby_pos_i - alt_pos_i) <= 50:
                    found = True
            if not found:
                fp_count += 1

print(fp_count)

