import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import random
from intervaltree import Interval, IntervalTree


def get_polya_len(seq):
    i = -1
    n = 0 
    while(n < len(seq)):
        if seq[i] != "A":
            break
        i -= 1
        n += 1
    return n

def get_tsd(seq):
    n = random.randint(6,22)
    return seq[len(seq)-n-1:]

def get_insert(seqs):
    n = random.randint(1,3)
    rt = ""
    if n == 1:
        rt = "LINE"
    elif n == 2:
        rt = "Alu"
    else:
        rt = "SVA"
    if rt == "LINE":
        # Add a randomly long Poly-A Tail to the sequence beetween 10 and 40 bases
        name = list(seqs[rt].keys())[0]
        n = random.randint(10,40)
        polya = "A"*n
        sequence = seqs[rt][name] + polya
        # determine how much of the sequence we're going to take, truncate from end
        n = random.randint(0, 99)
        if n < 25:
            # ~25% chance of a full length LINE
            return[name, sequence, len(polya)]
        else:
            nbases = int(float(n)/100 * len(sequence))
            if nbases < 100:
                nbases = 100
            #print(str(len(sequence))+"\t"+str(nbases))
            return[name, sequence[(len(sequence)-nbases):], len(polya)]
    else:
        # SVA and Alu seqs are Poly-A tailed randomly select one
        n = random.randint(0, len(seqs[rt])-1)
        name = list(seqs[rt].keys())[n]
        sequence = seqs[rt][name]
        n = random.randint(0, 99)
        if n < 50:
            # ~50% chance of a full length SVA or ALu
            return[name, sequence, get_polya_len(sequence)]
        else:
            nbases = int(float(n)/100 * len(sequence))
            if nbases < 100:
                nbases = 100
            #print(str(len(sequence))+"\t"+str(nbases))
            return[name, sequence[(len(sequence)-nbases):], get_polya_len(sequence[(len(sequence)-nbases):])]


def get_position(cigar, ref_pos, ref_start):
    ref_count = ref_start
    read_count = 0
    prev_ref_count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigar):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith('D'):
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
        if ref_count >= ref_pos:
            # reached the position we need
            if cg.endswith('M') or cg.endswith('X') or cg.endswith("="):
                diff = ref_count - ref_pos
                return read_count - diff
            else:
                return read_count
    return -1


def get_contig_insertion_seqs(hap_positions, hap_contigs, chromosome_lengths, centromeres, seqs):
    chromosomes = list(hap_positions.keys())
    n = random.randint(0, len(chromosomes)-1)
    chrom = chromosomes[n]
    # Randomly select a position on this chromosome to make an insertion into and select the best contig at the position
    if chrom not in chromosome_lengths:
        return ""
    pos = random.randint(50000,chromosome_lengths[chrom]-50001)
    nearby = centromeres[chrom][pos:pos+1]
    positions = hap_positions[chrom][pos:pos+1]
    found = False
    i = 0
    while(i < 50000):
        j = 0
        while(len(nearby) > 0):
            pos = random.randint(50000,chromosome_lengths[chrom]-50001)
            nearby = centromeres[chrom][pos:pos+1]
            j += 1
            if j > 50000:
                break
        positions = hap_positions[chrom][pos:pos+1]
        if len(positions) > 0:
            # Have a contig that we can use and select a position from to insert
            found = True
            break
        i += 1
    if found:
        length = 0
        name = ""
        cigar = ""
        start = 0
        end = 0
        for item in positions:
            if item.data[2] - item.data[1] > length:
                length = item.data[2] - item.data[1]
                name = item.data[0]
                cigar = item.data[3]
                start = item.begin
                end = item.end
        # Find the position on the contig that lines up with our position n
        contig_pos = get_position(cigar, pos, start)
        #if contig_pos < 0:
        #    print("Could not get position for "+name+" "+str(start)+" "+str(pos)+" "+str(end))
        left_flank = 50000
        right_flank = 50000
        if contig_pos < 50000:
            left_flank = contig_pos
        if len(hap_contigs[name]) - contig_pos < 50000:
            right_flank = len(hap_contigs[name]) - contig_pos - 1
        insert = get_insert(seqs)
        tsd = get_tsd(hap_contigs[name][contig_pos - left_flank:contig_pos])
        out_seq = hap_contigs[name][contig_pos - left_flank:contig_pos] + insert[1] + tsd + hap_contigs[name][contig_pos:contig_pos+right_flank]
        #print(name+"\t"+str(contig_pos)+"\t"+chrom+"\t"+str(pos)+"\t"+insert[0]+"\t"+insert[1]+"\t"+str(insert[2]))
        return [out_seq, name+"\t"+str(contig_pos)+"\t"+str(left_flank)+"\t"+str(right_flank)+"\t"+chrom+"\t"+str(pos)+"\t"+insert[0]+"\t"+insert[1]+"\t"+str(insert[2])+"\t"+tsd]
    return ""

parser = argparse.ArgumentParser( description='Add insertion sequences to contigs')
parser.add_argument('--input', required=True)
parser.add_argument('--mat-bam', required=True)
parser.add_argument('--pat-bam', required=True)
parser.add_argument('--mat-fasta', required=True)
parser.add_argument('--pat-fasta', required=True)
parser.add_argument('--mat-out-fa', required=True)
parser.add_argument('--pat-out-fa', required=True)
parser.add_argument('--centromeres', required=True)
parser.add_argument('--total', required=True, type=int)
parser.add_argument('--chrom-lengths', required=True)
parser.add_argument('--output-folder', required=True)
parser.add_argument('--output-prefix', required=True)
args = parser.parse_args()

random.seed(10)

centromeres = defaultdict(IntervalTree)
with open(args.centromeres) as in_cf:
    count = 0
    for line in in_cf:
        if count == 0:
            count = 1
            continue
        line_args = line.strip().split('\t')
        chrom = line_args[1]
        start = int(line_args[2])
        end = int(line_args[3])
        key = chrom + ":" + str(start) + "-" + str(end)
        centromeres[chrom][start:end] = key

chromosome_lengths = defaultdict(int)
with open(args.chrom_lengths) as in_cf:
    for line in in_cf:
        line_args = line.strip().split('\t')
        chrom = line_args[0]
        start = int(line_args[1])
        chromosome_lengths[chrom] = start

# Read in and store base concencus seuqnces
seqs = {}
seqs["LINE"] = defaultdict(str)
seqs["Alu"] = defaultdict(str)
seqs["SVA"] = defaultdict(str)
with pysam.FastxFile(args.input) as fh:
    for entry in fh:
        if "LINE" in entry.name:
            seqs["LINE"][entry.name] = entry.sequence.upper()
        elif "Alu" in entry.name:
            seqs["Alu"][entry.name] = entry.sequence.upper()
        elif "SVA" in entry.name:
            seqs["SVA"][entry.name] = entry.sequence.upper()

# Read in contig sequences
mat_contigs = defaultdict(str)
pat_contigs = defaultdict(str)
with pysam.FastxFile(args.mat_fasta) as fh:
    for entry in fh:
        mat_contigs[entry.name] = entry.sequence
with pysam.FastxFile(args.pat_fasta) as fh:
    for entry in fh:
        pat_contigs[entry.name] = entry.sequence        

mapped_count = defaultdict(IntervalTree)
mat = pysam.AlignmentFile(args.mat_bam, "rb")
pat = pysam.AlignmentFile(args.pat_bam, "rb")
for record in mat.fetch():
    if record.mapping_quality > 0 and len(mat_contigs[record.query_name]) > 100000:
        mapped_count[record.query_name][record.query_alignment_start:record.query_alignment_end] = 1
for record in pat.fetch():
    if record.mapping_quality > 0 and len(pat_contigs[record.query_name]) > 100000:
        mapped_count[record.query_name][record.query_alignment_start:record.query_alignment_end] = 1

mat_positions = defaultdict(IntervalTree)
pat_positions = defaultdict(IntervalTree)
# Read in contig to ref bam, store contigs that map with qual > 20, and positions that are not centromere alignments
# Store alignment recods for each of these contigs and store by chromsome positions
# Will allow us to later on randomly select a chromsomal position, see if a contig aligns and add an insertion 
mat = pysam.AlignmentFile(args.mat_bam, "rb")
seen = {}
with open(args.mat_out_fa, 'w') as out_fa:
    for record in mat.fetch():
        nearby = mapped_count[record.query_name][record.query_alignment_start:record.query_alignment_end]
        if len(nearby) > 1:
            continue
        if record.mapping_quality >=50 and record.reference_length > 100000:
            #print(record.query_name+"\t"+str(record.reference_length)+"\t"+str(record.query_alignment_start)+"\t"+str(record.query_alignment_end)+"\t"+record.cigarstring)
            mat_positions[record.reference_name][record.reference_start:record.reference_end] = [record.query_name, record.query_alignment_start, record.query_alignment_end, record.cigarstring]
            if record.query_name not in seen:
                out_fa.write(">"+record.query_name+'\n'+mat_contigs[record.query_name]+"\n")
                seen[record.query_name] = 1
pat = pysam.AlignmentFile(args.pat_bam, "rb")
with open(args.pat_out_fa, 'w') as out_fa:
    for record in pat.fetch():
        nearby = mapped_count[record.query_name][record.query_alignment_start:record.query_alignment_end]
        if len(nearby) > 1:
            continue
        if record.mapping_quality >=50 and record.reference_length > 100000:
            #print(record.query_name+"\t"+str(record.reference_length)+"\t"+str(record.query_alignment_start)+"\t"+str(record.query_alignment_end)+"\t"+record.cigarstring)
            pat_positions[record.reference_name][record.reference_start:record.reference_end] = [record.query_name, record.query_alignment_start, record.query_alignment_end, record.cigarstring]
            if record.query_name not in seen:
                out_fa.write(">"+record.query_name+'\n'+pat_contigs[record.query_name]+"\n")
                seen[record.query_name] = 1


# Have a set of contigs aligned to the reference. Now randomly select mat or pat, a chromosome and a contig from the set and a position on the contig to insert into
# Once we have selected a position to insert into, have to select which insert to get
count = 0
while(count < args.total):
    ret = ""
    n = random.randint(0,1)
    if n == 1: 
        # Mat
        ret = get_contig_insertion_seqs(mat_positions, mat_contigs, chromosome_lengths, centromeres, seqs)
    else:
        # Pat
        ret = get_contig_insertion_seqs(pat_positions, pat_contigs, chromosome_lengths, centromeres, seqs)
    if len(ret) > 0:
        print(str(count)+"\t"+ret[1])
        with open(args.output_folder+"/"+args.output_prefix+"."+str(count)+".fa", 'w') as out_fa:
            out_fa.write(">"+str(count)+"\n"+ret[0]+"\n")
        count += 1 


