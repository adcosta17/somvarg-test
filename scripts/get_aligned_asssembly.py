import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re

parser = argparse.ArgumentParser( description='Extract Contigs from assemblies aligned to chroms and generate subset assemblies')
parser.add_argument('--chrom-list', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--assembly-path', default="/.mounts/labs/ont/scratch/HPRC_human_data")
parser.add_argument('--assembly-folder', default="assembly")
parser.add_argument('--input-path', required=True)
parser.add_argument('--input-folder', default="assembly_mapped")
parser.add_argument('--output-folder', default="assembly_subset")
parser.add_argument('--min-mapping-qual', type=int, default=60)
parser.add_argument('--min-mapped-fraction', type=float, default=0.3)
parser.add_argument('--min-count', type=int, default=100000)
args = parser.parse_args()

sample = args.sample
contigs = {}
contigs[sample] = {}
with pysam.FastxFile(args.assembly_path+'/'+sample+"/"+args.assembly_folder+"/"+sample+".mat.fa") as in_mat:
	contigs[sample]["mat"] = defaultdict(str)
	for entry in in_mat:
		contigs[sample]["mat"][entry.name] = entry.sequence
with pysam.FastxFile(args.assembly_path+'/'+sample+"/"+args.assembly_folder+"/"+sample+".pat.fa") as in_pat:
	contigs[sample]["pat"] = defaultdict(str)
	for entry in in_pat:
		contigs[sample]["pat"][entry.name] = entry.sequence

contigs_to_output = {}
contigs_to_output[sample] = {}
for hp in ["mat", "pat"]:
	samfile = pysam.AlignmentFile(args.input_path+'/'+sample+"/"+args.input_folder+"/"+sample+"."+hp+".sorted.bam", "rb")
	contigs_to_output[sample][hp] = defaultdict(int)
	for chrom in args.chrom_list.split(','):
		contigs_per_sample = defaultdict(list)
		contigs_to_use = {}
		sam_iter = samfile.fetch(region=chrom)
		for record in sam_iter:
			contigs_per_sample[record.query_name].append((record.query_alignment_start, record.query_alignment_end-record.query_alignment_start))
			if record.mapping_quality < args.min_mapping_qual:
				continue
			contigs_to_use[record.query_name] = 1
		for contig in contigs_to_use:
			positions = sorted(contigs_per_sample[contig])
			count = 0
			pos = 0
			for item in positions:
				start = item[0]
				length = item[1]
				if start+length <= pos:
					# This position has already been aligned
					continue
				if start <= pos:
					# have an overlapping region 
					count = count+(start+length-pos)
					pos = start+length
				else:
					# have a gap
					count = count + length
					pos = start+length
			if len(contigs[sample][hp][contig]) > 0:
				contigs_to_output[sample][hp][contig] = 1

for hp in ["mat", "pat"]:
	with open(args.input_path+'/'+sample+"/"+args.output_folder+"/"+sample+"."+hp+".subset.fa", 'w') as out_fa:
		for contig in contigs_to_output[sample][hp]:
			out_fa.write(">"+contig+"\n"+contigs[sample][hp][contig]+"\n")






