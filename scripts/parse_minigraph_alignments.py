import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
from intervaltree import Interval, IntervalTree

def get_insert_in_cigar(cigar, start, read, insert_pos, min_size=50):
    ref_count = start
    read_count = 0
    found_insert = False
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
            read_count += int(cg[:cg.find("I")])
            if abs(ref_count - insert_pos) < 500 and int(cg[:cg.find("I")]) >= min_size:
                found_insert = True
            #if int(cg[:cg.find("I")]) >= min_size and read == "6c4b9283-f68c-46a8-a6be-03ee0009779c":
            #    print("Insert: " + str(ref_count) + " " + cg, file=sys.stderr)
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('P'):
            read_count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
    #if read == "6c4b9283-f68c-46a8-a6be-03ee0009779c":
    #    print(ref_count, file=sys.stderr)
    #    print(read_count, file=sys.stderr)
    return found_insert

def get_read_insert_position(cigar, start, read, insert_pos, min_size=50):
    ref_count = 0
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
                read_insert_positions.append([start+read_start,start+read_end])
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('P'):
            read_count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
    return read_insert_positions

parser = argparse.ArgumentParser( description='Parse Minigraph alignments of reads to an augmented graph. ')
parser.add_argument('--gaf', required=True)
parser.add_argument('--gfa', required=True)
parser.add_argument('--tsv', required=True)
parser.add_argument('--vcf', required=True)
#parser.add_argument('--mat-ref-inserts', required=True)
#parser.add_argument('--pat-ref-inserts', required=True)
args = parser.parse_args()

# Read in the tsv of reads and the insertions they are supposed to support
# TSV based on read alignments to the assembly contigs, and alignments of assembly contigs to reference
truth_data = defaultdict(list)
truth_positions = defaultdict(int)
control_data = defaultdict(list)
control_positions = defaultdict(int)
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
        if "Control" in row[0]:
            control_data[row[1]].append((row[2], int(row[3])))
            #print(row[1]+"\t"+str(control_data[row[1]]))
            control_positions[row[2]+"_"+row[3]] += 1
        else:
            truth_data[row[1]].append((row[2], int(row[3])))
            truth_positions[row[2]+"_"+row[3]] += 1

# Read in the VCF of insertions that are used to augment the graph 
vcf_positions = defaultdict(IntervalTree)
#bcf_in = pysam.VariantFile(args.vcf)
#for rec in bcf_in.fetch():
#    if rec.info["SVTYPE"] == "INS":
#        vcf_positions[rec.contig][(rec.pos-500):(rec.pos+500)] = len(rec.alts[0])

with open(args.vcf, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        row = line.strip().split('\t')
        positions = row[2].split(',')
        for pos in positions:
            chrom = pos.split('_')[1].split(':')[0]
            start = int(pos.split('_')[1].split(':')[1])
            vcf_positions[chrom][start-500:start+500] = 1


# Read in GFA File of what SVs actually got integrated into the graph
# Some SV calls won't get added to the graph, need to check which ones did to evaluate FPs
gfa_positions = defaultdict(IntervalTree)
with open(args.gfa, 'r') as in_graph:
    for line in in_graph:
        row = line.strip().split('\t')
        if row[0] == "S":
            # Have a segment line
            name = row[4].split(':')[2]
            if '_' in name:
                chrom = name.split('_')[0]
                min_pos = int(name.split('_')[2])
                max_pos = int(name.split('_')[3])
                gfa_positions[chrom][min_pos-500:max_pos+500] = name


chrs_to_use = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

# Identify if we have a TP, a path alignment to an SV within 500bp of the expected position
control_tp = {}
control_tp_positions = defaultdict(list)
control_missed_fn = {}
control_missed_fn_positions = defaultdict(int)
control_seen_fn = {}
control_seen_fn_positions = defaultdict(list)
control_seen_fn_elsewhere = {}
control_seen_fn_elsewhere_positions = defaultdict(list)
control_seen_fn_inserts = {}
control_seen_fn_inserts_positions = defaultdict(list)
with open(args.gaf, 'r') as in_alignment:
    for line in in_alignment:
        row = line.strip().split("\t")
        if '>' in row[5] or '<' in row[5]:
            # Check the position of where the read aligns is where we expect or if it is an fp
            path = row[5].replace('>',' ').replace('<', ' ').split(' ')
            for item in path:
                if "_" in item:
                    if "random" in item:
                        continue
                    #print(item)
                    # Have an SV position we've augmented
                    chrom = item.split('_')[0]
                    if chrom not in chrs_to_use:
                        continue
                    positions = item.split('_')[4].split(',')
                    # compare the pos to the truth data, see if tp or fp
                    #print(control_data[row[0]])
                    for sv in control_data[row[0]]:
                        if chrom != sv[0]:
                            continue
                        #print(str(row)+"\t"+str(sv[1]))
                        for pos in positions:
                            p = pos.replace('[','').replace(']','')
                            contig_pos = int(p.split('.')[1])
                            if abs(contig_pos-sv[1]) < 500:
                                control_tp[row[0]+"_"+sv[0]+"_"+str(sv[1])] = 1
                                control_tp_positions[sv[0]+"_"+str(sv[1])].append(row[0]+"_"+chrom+"_"+str(contig_pos))

# Idenfity FNs Start with FNs that intersect the expected insertion
count = 0
with open(args.gaf, 'r') as in_alignment:
    for line in in_alignment:
        row = line.strip().split("\t")
        # Check if the read spans an insertion we augmented with
        if '>' in row[5] or '<' in row[5]:
            continue
        elif row[0] in control_data:
            chrom = row[5]
            start = int(row[7])
            end = int(row[8])
            for item in control_data[row[0]]:
                # Check if the position intersects the expected one, if not then don't count this as a FN as it could be a sup/secondary alignment
                if row[0]+"_"+item[0]+"_"+str(item[1]) in control_tp:
                    continue
                if start < item[1] and end > item[1] and chrom == item[0]:
                    # Have an alignment that intersects the insertion position
                    # Identify if the CIGAR string contains an insertion within 500bp of the expected position
                    # If so note the read as seem with insert in CIGAR
                    #if row[0] == "6c4b9283-f68c-46a8-a6be-03ee0009779c" and count == 0:
                    #    print(row, file=sys.stderr)
                    #    print(row[2], file=sys.stderr)
                    #    print(row[3], file=sys.stderr)
                    #    count += 1
                    read_start = int(row[2])
                    if row[4] == '-':
                        read_start = int(row[3])
                    positions = get_read_insert_position(row[18].split(':')[2], read_start, row[0], item[1])
                    if get_insert_in_cigar(row[18].split(':')[2], start, row[0], item[1]):
                        control_seen_fn_inserts[row[0]+"_"+item[0]+"_"+str(item[1])] = 1
                        control_seen_fn_inserts_positions[item[0]+"_"+str(item[1])].append(row[0]+"_"+chrom+":"+str(start)+"-"+str(end)+"_"+str(positions))
                    else:
                        control_seen_fn[row[0]+"_"+item[0]+"_"+str(item[1])] = 1
                        control_seen_fn_positions[item[0]+"_"+str(item[1])].append(row[0]+"_"+chrom+":"+str(start)+"-"+str(end)+"_"+str(positions))
                    #print(item[0]+"\t"+str(item[1])+"\t"+row[0]+"\t"+row[5]+"\t"+row[7]+"\t"+row[8])
        else:
            continue

# Next look at FNs that align somewhere else that we don't expect them too
with open(args.gaf, 'r') as in_alignment:
    for line in in_alignment:
        row = line.strip().split("\t")
        # Check if the read spans an insertion we augmented with
        if '>' in row[5] or '<' in row[5]:
            continue
        elif row[0] in control_data:
            chrom = row[5]
            start = int(row[7])
            end = int(row[8])
            for item in control_data[row[0]]:
                # Check if the position intersects the expected one, if not then don't count this as a FN as it could be a sup/secondary alignment
                if row[0]+"_"+item[0]+"_"+str(item[1]) in control_tp or row[0]+"_"+item[0]+"_"+str(item[1]) in control_seen_fn or row[0]+"_"+item[0]+"_"+str(item[1]) in control_seen_fn_inserts:
                    continue
                # If here we don't have any alignment for this read that intersects the expected position, add to elsewhere FN list        
                control_seen_fn_elsewhere[row[0]+"_"+item[0]+"_"+str(item[1])] = 1
                control_seen_fn_elsewhere_positions[item[0]+"_"+str(item[1])].append(row[0]+"_"+chrom+":"+str(start)+"-"+str(end))
        else:
            continue

# Look to see which FNs we missed completly, Reads that did not align to the graph at all
for read in control_data:
    # check if we've seen the SV for the read
    for item in control_data[read]:
        if read+"_"+item[0]+"_"+str(item[1]) in control_tp or read+"_"+item[0]+"_"+str(item[1]) in control_seen_fn or read+"_"+item[0]+"_"+str(item[1]) in control_seen_fn_inserts or read+"_"+item[0]+"_"+str(item[1]) in control_seen_fn_elsewhere:
            continue
        else:
            control_missed_fn[read+"_"+item[0]+"_"+str(item[1])] = 1
            control_missed_fn_positions[item[0]+"_"+str(item[1])] += 1


print(len(control_tp))
print(len(control_seen_fn_inserts))
print(len(control_seen_fn))
print(len(control_seen_fn_elsewhere))
print(len(control_missed_fn))
exit(0)

print("AssemblyChrom\tAssemblyPos\tInVCF_POS_SVLEN\tInGraph_POS\tTP_Count\tFN_WithInsert_Count\tFN_NoInsert_Count\tFN_Elsewhere_Count\tFN_Missed_Count\tFN_WithInsert_Reads\tFN_NoInsert_Reads\tFN_Elsehwere_Reads")
for pos in control_positions:
    tp_count = 0
    if pos in control_tp_positions:
        tp_count = len(control_tp_positions[pos])
    nearby = vcf_positions[pos.split('_')[0]][int(pos.split('_')[1]):(int(pos.split('_')[1])+1)]
    in_vcf = "False"
    if len(nearby) > 0:
        in_vcf = "True:_"
        for item in nearby:
            in_vcf += str(item.begin+500)+"_"+str(item.data)+","
    seen_with_insert_count = 0
    seen_with_insert = ""
    if pos in control_seen_fn_inserts_positions:
        seen_with_insert_count = len(control_seen_fn_inserts_positions[pos])
        seen_with_insert = ','.join(control_seen_fn_inserts_positions[pos])
    else:
        seen_with_insert = "-"
    seen_elsewhere_count = 0
    seen_elsewhere = ""
    if pos in control_seen_fn_elsewhere_positions:
        seen_elsewhere_count = len(control_seen_fn_elsewhere_positions[pos])
        seen_elsewhere = ','.join(control_seen_fn_elsewhere_positions[pos])
    else:
        seen_elsewhere = "-"
    seen_no_insert_count = 0
    seen_no_insert = ""
    if pos in control_seen_fn_positions:
        seen_no_insert_count = len(control_seen_fn_positions[pos])
        seen_no_insert = ','.join(control_seen_fn_positions[pos])
    else:
        seen_no_insert = "-"
    missed_count = 0
    if pos in control_missed_fn_positions:
        missed_count = control_missed_fn_positions[pos]
    in_graph = "False"
    nearby = gfa_positions[pos.split('_')[0]][int(pos.split('_')[1]):(int(pos.split('_')[1])+1)]
    if len(nearby) > 0:
        in_graph = "True:_"
        for item in nearby:
            in_graph += item.data+","
    print(pos.split('_')[0]+"\t"+pos.split('_')[1]+"\t"+in_vcf+"\t"+in_graph+"\t"+str(tp_count)+"\t"+str(seen_with_insert_count)+"\t"+str(seen_no_insert_count)+"\t"+str(seen_elsewhere_count)+"\t"+str(missed_count)+"\t"+seen_with_insert+"\t"+seen_no_insert+"\t"+seen_elsewhere)
