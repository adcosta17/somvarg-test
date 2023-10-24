import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import gzip
import random
from intervaltree import Interval, IntervalTree

random.seed(10)

def get_read_length(cigarstring):
    # Gets the read length based on the Cigar String
    count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            count += int(cg[:cg.find("H")])
    return count


def get_insert_position_for_read(record, ref_start, ref_end, control):
    ref_count = record.reference_start
    prev_ref_count = 0
    failed_start = False
    read_count = 0
    read_start = -1
    read_end = -1
    for cg in re.findall('[0-9]*[A-Z=]', record.cigarstring):
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
        if ref_count > ref_start and read_start < 0:
            read_start = read_count
            if cg.endswith('M') or cg.endswith('X') or cg.endswith('='):
                read_start = read_start - (ref_count - ref_start)
            if prev_ref_count < 500 and not control:
                failed_start = True
        if ref_count > ref_end and read_end < 0:
            read_end = read_count
            if cg.endswith('M') or cg.endswith('X') or cg.endswith('='):
                read_end = read_end - (ref_count - ref_end)
    if (ref_count - ref_end < 500 or failed_start) and not control:
        read_end = -1
    return (read_start, read_end)



def get_insert_position_for_control(record, ref_start, ref_end):
    ref_count = record.reference_start
    prev_ref_count = 0
    failed_start = False
    read_count = 0
    read_start = -1
    read_end = -1
    positions = {}
    count = 0
    for cg in re.findall('[0-9]*[A-Z=]', record.cigarstring):
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
            start = -1
            end = -1
            if(int(cg[:cg.find("I")])) >= 100 and ref_count > ref_start and ref_count < ref_end:
                start = read_count
            read_count += int(cg[:cg.find("I")])
            if(int(cg[:cg.find("I")])) >= 100 and ref_count > ref_start and ref_count < ref_end:
                end = read_count
                positions[count] = [start,end]
                count += 1
        elif cg.endswith('D'):
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
        if ref_count > ref_start and read_start < 0:
            read_start = read_count
            if cg.endswith('M') or cg.endswith('X') or cg.endswith('='):
                read_start = read_start - (ref_count - ref_start)
            if prev_ref_count < 500:
                failed_start = True
        if ref_count > ref_end and read_end < 0:
            read_end = read_count
            if cg.endswith('M') or cg.endswith('X') or cg.endswith('='):
                read_end = read_end - (ref_count - ref_end)
    if ref_count - ref_end < 500 or failed_start:
        read_end = -1
    return positions

def mapped_end_to_end(record, read_length):
    ref_count = record.reference_start
    read_count = 0
    read_start = -1
    read_end = -1
    large_indel = False
    for cg in re.findall('[0-9]*[A-Z=]', record.cigarstring):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            read_count += int(cg[:cg.find("I")])
            if int(cg[:cg.find("I")]) >= 100:
                large_indel = True
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
            if int(cg[:cg.find("D")]) >= 100:
                large_indel = True
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
            if read_start == -1:
                read_start = read_count
            elif read_end == -1:
                read_end = read_count - int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
            if read_start == -1:
                read_start = read_count
            elif read_end == -1:
                read_end = read_count - int(cg[:cg.find("H")])
        if read_start == -1:
            read_start = read_count
    if read_end == -1:
        read_end = read_count
    if (read_end - read_start)/read_length < 0.5:
        # Need at least half the read to map
        return False
    if large_indel:
        return False
    return True


def get_inserts_for_control(ref_bam, mat_and_pat_inserts, chr_list):
    ret = {}
    bam = pysam.AlignmentFile(ref_bam, "rb")
    for chrom in mat_and_pat_inserts:
        if chrom in chr_list:
            for item in mat_and_pat_inserts[chrom]:
                insert_dict = item.data
                read_positions = {}
                to_print = False
                if "mat" in insert_dict:
                    mat_pos = insert_dict["mat"][0]+":"+str(int(insert_dict["mat"][1])-500)+"-"+str(int(insert_dict["mat"][2])+500)
                    for record in bam.fetch(region=mat_pos):
                        # Check read mapped end to end with no inserts
                        positions = get_insert_position_for_control(record, int(insert_dict["mat"][1])-500, int(insert_dict["mat"][2])+500)
                        if len(positions) == 0:
                            continue
                        read_start = -1
                        read_end = -1
                        for p in positions:
                            if read_start == -1:
                                read_start = positions[p][0]
                            if read_end == -1:
                                read_end = positions[p][1]
                            if positions[p][0] < read_start:
                                read_start = positions[p][0]
                            if positions[p][1] > read_end:
                                read_end = positions[p][1]
                        read_length = get_read_length(record.cigarstring)
                        if read_start < 0 or read_end < 0:
                            continue
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        read_positions[record.query_name] = [read_start, read_end, record.cigarstring, record.is_reverse, ]
                if "pat" in insert_dict:
                    pat_pos = insert_dict["pat"][0]+":"+str(int(insert_dict["pat"][1])-500)+"-"+str(int(insert_dict["pat"][2])+500)
                    for record in bam.fetch(region=pat_pos):
                        positions = get_insert_position_for_control(record, int(insert_dict["pat"][1])-500, int(insert_dict["pat"][2])+500)
                        if len(positions) == 0:
                            continue
                        read_start = -1
                        read_end = -1
                        for p in positions:
                            if read_start == -1:
                                read_start = positions[p][0]
                            if read_end == -1:
                                read_end = positions[p][1]
                            if positions[p][0] < read_start:
                                read_start = positions[p][0]
                            if positions[p][1] > read_end:
                                read_end = positions[p][1]
                        read_length = get_read_length(record.cigarstring)
                        if read_start < 0 or read_end < 0:
                            continue
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        read_positions[record.query_name] = [read_start, read_end, record.cigarstring, record.is_reverse]
                ret[chrom+":"+str(item.begin)+"-"+str(item.end)] = read_positions
    return ret

def spanning_read_check(cigarstring, read_length):
    read_start =  -1
    read_end = -1
    read_count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            if read_start == -1:
                # First match
                read_start = read_count
            read_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            if read_start == -1:
                # First match
                read_start = read_count
            read_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            if read_start == -1:
                # First match
                read_start = read_count
            read_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith('S'):
            if read_start > -1 and read_end == -1:
                read_end = read_count
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            if read_start > -1 and read_end == -1:
                read_end = read_count
            read_count += int(cg[:cg.find("H")])
    if read_start < 250 and (read_length - read_end) < 250:
        return True
    return False 

def get_inserts_for_sample(mat_bam, pat_bam, mat_and_pat_inserts, chr_list, file, control):
    ret = {}
    depths = {}
    mat = pysam.AlignmentFile(mat_bam, "rb")
    pat = pysam.AlignmentFile(pat_bam, "rb")
    for chrom in mat_and_pat_inserts:
        if chrom in chr_list:
            for item in mat_and_pat_inserts[chrom]:
                insert_dict = item.data
                read_positions = {}
                to_print = False
                #if chrom == "chr7" and item.begin > 85262406 and item.end < 85282407:
                #    to_print = True
                depth = 0
                if "mat" in insert_dict and "pat" in insert_dict:
                    # Get the reads that aligned to the mat first, and then for pat
                    mat_pos = insert_dict["mat"][3]+":"+str(int(insert_dict["mat"][4])-100)+"-"+str(int(insert_dict["mat"][5])+100)
                    for record in mat.fetch(region=mat_pos):
                        if record.mapping_quality < 20:
                            continue
                        read_length = get_read_length(record.cigarstring)
                        if not spanning_read_check(record.cigarstring, read_length):
                            continue
                        read_start, read_end = get_insert_position_for_read(record, int(insert_dict["mat"][4]), int(insert_dict["mat"][5]), control)
                        if file is not None:
                            with open(file, "a") as myfile:
                                myfile.write(record.query_name+"\n")
                        if control:
                            # delete partial insertions as well
                            if read_start < 0 and read_end < 0:
                                continue
                            if read_end < 0:
                                # Entered but didn't get to the end of the insert on the read
                                read_end = read_length
                        elif read_start < 0 or read_end < 0:
                            continue
                        if read_end - read_start < 100:
                            # Require at least 100 bp of read to overlap the insert sequence
                            continue
                        depth += 1
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        #if record.query_name == "90f746b1-3df1-4d99-989c-4507e3a3ca28":
                        #    print(record.query_name+"\t"+str(read_start)+"\t"+str(read_end))
                        read_positions[record.query_name] = [read_start, read_end, insert_dict["mat"][3]+":"+str(int(insert_dict["mat"][4]))+"-"+str(int(insert_dict["mat"][5]))]
                    pat_pos = insert_dict["pat"][3]+":"+str(int(insert_dict["pat"][4])-100)+"-"+str(int(insert_dict["pat"][5])+100)
                    for record in pat.fetch(region=pat_pos):
                        if record.mapping_quality < 20:
                            continue
                        read_length = get_read_length(record.cigarstring)
                        if not spanning_read_check(record.cigarstring, read_length):
                            continue
                        read_start, read_end = get_insert_position_for_read(record, int(insert_dict["pat"][4]), int(insert_dict["pat"][5]), control)
                        if file is not None:
                            with open(file, "a") as myfile:
                                myfile.write(record.query_name+"\n")
                        if control:
                            # delete partial insertions as well
                            if read_start < 0 and read_end < 0:
                                continue
                            if read_end < 0:
                                # Entered but didn't get to the end of the insert on the read
                                read_end = read_length
                        if read_start < 0 or read_end < 0:
                            continue
                        if read_end - read_start < 100:
                            # Require at least 100 bp of read to overlap the insert sequence
                            continue
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        if record.query_name not in read_positions:
                            depth += 1
                            #if record.query_name == "90f746b1-3df1-4d99-989c-4507e3a3ca28":
                            #   print(record.query_name+"\t"+str(read_start)+"\t"+str(read_end))
                            read_positions[record.query_name] = [read_start, read_end, insert_dict["pat"][3]+":"+str(int(insert_dict["pat"][4]))+"-"+str(int(insert_dict["pat"][5]))]
                elif "mat" in insert_dict:
                    mat_pos = insert_dict["mat"][3]+":"+str(int(insert_dict["mat"][4])-100)+"-"+str(int(insert_dict["mat"][5])+100)
                    for record in mat.fetch(region=mat_pos):
                        if record.mapping_quality < 20:
                            continue
                        read_length = get_read_length(record.cigarstring)
                        if not spanning_read_check(record.cigarstring, read_length):
                            continue
                        # Check read mapped end to end with no inserts
                        if file is not None:
                            with open(file, "a") as myfile:
                                myfile.write(record.query_name+"\n")
                        #if not mapped_end_to_end(record, read_length):
                        #    continue
                        read_start, read_end = get_insert_position_for_read(record, int(insert_dict["mat"][4]), int(insert_dict["mat"][5]), control)
                        if control:
                            # delete partial insertions as well
                            if read_start < 0 and read_end < 0:
                                continue
                            if read_end < 0:
                                # Entered but didn't get to the end of the insert on the read
                                read_end = read_length
                        if read_start < 0 or read_end < 0:
                            continue
                        if read_end - read_start < 100:
                            # Require at least 100 bp of read to overlap the insert sequence
                            continue
                        depth += 1
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        #if record.query_name == "90f746b1-3df1-4d99-989c-4507e3a3ca28":
                        #    print(record.query_name+"\t"+str(read_start)+"\t"+str(read_end))
                        read_positions[record.query_name] = [read_start, read_end, insert_dict["mat"][3]+":"+str(int(insert_dict["mat"][4]))+"-"+str(int(insert_dict["mat"][5]))]
                elif "pat" in insert_dict:
                    pat_pos = insert_dict["pat"][3]+":"+str(int(insert_dict["pat"][4])-100)+"-"+str(int(insert_dict["pat"][5])+100)
                    for record in pat.fetch(region=pat_pos):
                        if record.mapping_quality < 20:
                            continue
                        read_length = get_read_length(record.cigarstring)
                        if not spanning_read_check(record.cigarstring, read_length):
                            continue
                        if file is not None:
                            with open(file, "a") as myfile:
                                myfile.write(record.query_name+"\n")
                        #if not mapped_end_to_end(record, read_length):
                        #    continue
                        read_start, read_end = get_insert_position_for_read(record, int(insert_dict["pat"][4]), int(insert_dict["pat"][5]), control)
                        if control:
                            # delete partial insertions as well
                            if read_start < 0 and read_end < 0:
                                continue
                            if read_end < 0:
                                # Entered but didn't get to the end of the insert on the read
                                read_end = read_length
                        if read_start < 0 or read_end < 0:
                            continue
                        if read_end - read_start < 100:
                            # Require at least 100 bp of read to overlap the insert sequence
                            continue
                        depth += 1
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        read_positions[record.query_name] = [read_start, read_end, insert_dict["pat"][3]+":"+str(int(insert_dict["pat"][4]))+"-"+str(int(insert_dict["pat"][5]))]
                ret[chrom+":"+str(item.begin)+"-"+str(item.end)] = read_positions
                depths[chrom+":"+str(item.begin)+"-"+str(item.end)] = depth
    return ret, depths


def select_reads_to_delete(insert_positions, sample):
    for pos in insert_positions:
        for read in insert_positions[pos]:
            #print(pos,file=sys.stderr)
            #print(pos.split(':')[0]+"\t"+pos.split(':')[0].split("-")[0],file=sys.stderr)
            #print(insert_positions[pos][read],file=sys.stderr)
            #print(insert_positions[pos][read],file=sys.stderr)
            print(sample+"\t"+read+"\t"+pos.split(':')[0]+"\t"+pos.split(':')[1].split("-")[0]+"\t"+str(insert_positions[pos][read][0])+"\t"+str(insert_positions[pos][read][1])+"\t"+str(insert_positions[pos][read][2]))

def get_fastq(sample, args):
    if sample == "Control":
        return args.input_fastq_control
    if sample == "Test":
        return args.input_fastq_test
    return None

def get_output_fastq(sample, args):
    if sample in ["Control", "Test"]:
        return args.output_fasta_folder+"/"+sample+".fastq"
    return None

parser = argparse.ArgumentParser( description='Generate New fastqs with deletions in reads')
parser.add_argument('--input-fastq-control', required=True)
parser.add_argument('--input-fastq-test', required=True)
parser.add_argument('--centromeres', required=True)
parser.add_argument('--chr-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
parser.add_argument('--mat-ref-inserts', required=True)
parser.add_argument('--pat-ref-inserts', required=True)
parser.add_argument('--mat-bam-control', required=True)
parser.add_argument('--mat-bam-test', required=True)
parser.add_argument('--pat-bam-control', required=True)
parser.add_argument('--pat-bam-test', required=True)
args = parser.parse_args()

#print("centromeres")
# Read in Centromeres list
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

#print("mat")
# Store positions of hap to ref inserts on the ref. Will then identify reads that align to these positions on the haplotypes end to end
# Will then delete the insert sequence on some fraction of the reads that do align. 
mat_ref_inserts = defaultdict(IntervalTree)
combined_mat_pat_ref_inserts = defaultdict(IntervalTree)
count = 0
with open(args.mat_ref_inserts, 'r') as in_tsv:
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        nearby = centromeres[chrom][start-500:end+500]
        if row[8] != "PASS" or len(nearby) > 0:
            continue
        mat_ref_inserts[chrom][start:end] = row

#print("pat")
pat_ref_inserts = defaultdict(IntervalTree)
count = 0
with open(args.pat_ref_inserts, 'r') as in_tsv:
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        nearby = centromeres[chrom][start-500:end+500]
        if row[8] != "PASS" or len(nearby) > 0:
            continue
        pat_ref_inserts[chrom][start:end] = row
        nearby = mat_ref_inserts[chrom][start-500:end+500]
        if len(nearby) > 0:
            # Add both the mat and pat inserts to the combined at this position, should only have one item given the polished assembly
            for item in nearby:
                mat_row = item.data
                combined_mat_pat_ref_inserts[chrom][start:end] = {"pat": row, "mat": mat_row}
        else:
            # Add just the pat row
            combined_mat_pat_ref_inserts[chrom][start:end] = {"pat": row}

#print("mat2")
count = 0
with open(args.mat_ref_inserts, 'r') as in_tsv:
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        nearby = centromeres[chrom][start-500:end+500]
        if row[8] != "PASS" or len(nearby) > 0:
            continue
        nearby = pat_ref_inserts[chrom][start-500:end+500]
        if len(nearby) > 0:
            # Already added this insert to the combined set
            continue
        else:
            combined_mat_pat_ref_inserts[chrom][start:end] = {"mat": row}


# Read the bams for each sample against the mat inserts and the pat inserts. Decide which reads map to each sample at the position. 

# Now read in each of the read fastq inserts
# If the insert passes and lines up within 500 bp of an insert in the paternal contig set, flag it
# In sample 1 instantly collapse these to generate an insert free set, while in sample 3 and 5 collapse some fraction of the reads at each position
reads_with_inserts_per_sample = {}
depth_per_sample = {}
reads_with_inserts_per_sample["Control"], depth_per_sample["Control"] = get_inserts_for_sample(args.mat_bam_control, args.pat_bam_control, combined_mat_pat_ref_inserts, args.chr_list.split(","),None, True)
reads_with_inserts_per_sample["Test"], depth_per_sample["Test"] = get_inserts_for_sample(args.mat_bam_test, args.pat_bam_test, combined_mat_pat_ref_inserts, args.chr_list.split(","),None, True)


# Print header
print("Sample\tRead\tReadStart\tReadEnd\tRefPosition")

# Once we have insertions on reads, we now have to decide at each position how many reads to delete sequence and how many to retain
# For Sample 1 delete all inserts at every position - clean control sample
# For Samples 3 and 5 decide randomly based on some fraction of which to retian
for sample in ["Control", "Test"]:
    select_reads_to_delete(reads_with_inserts_per_sample[sample], sample)

