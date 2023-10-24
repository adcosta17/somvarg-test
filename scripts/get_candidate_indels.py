import pysam
import argparse
import sys
import re
from collections import defaultdict
from intervaltree import Interval, IntervalTree

def calculate_sdust_score(seq):
    if seq == "":
        return 0
    triplet_counts = defaultdict(int)
    for i in range(0, len(seq) - 2):
        triplet_counts[seq[i:i+3]] += 1
    sum_score = 0
    for triplet in triplet_counts:
        count = triplet_counts[triplet]
        s = float(count * (count - 1)) / 2.0
        sum_score += s
    if len(seq) - 1 == 0:
        return 0
    sum_score /= (len(seq) - 1)
    return sum_score


def get_deletion_pos(cigarstring):
    # Count up the position on the read until we get to a deletion
    deletion_positions = []
    count = 0
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            count += int(cg[:cg.find("H")])
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) >= 100:
                deletion_positions.append((count,int(cg[:cg.find("D")])))
    return deletion_positions

def get_hard_start(cigarstring):
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
        if cg.endswith('H'):
            return int(cg[:cg.find("H")])
        else:
            break
    return 0

def update_annotation(annotation, update):
    if annotation == "PASS":
        annotation = update
    else:
        annotation = annotation+","+update
    return annotation


def get_inserts_for_bam(args, bam, read_seqs):
    mapped_count = defaultdict(IntervalTree)
    max_ref_gap_at_candidate = args.reference_gap_minimum
    records_to_output = defaultdict(list)
    sam_reader = pysam.AlignmentFile(bam)
    for record in sam_reader.fetch():
        if record.is_unmapped:
            continue
        mapped_count[record.query_name][record.query_alignment_start:record.query_alignment_end] = 1
        read_annotation = "PASS"
        if record.mapq < args.min_mapq:
            read_annotation = "mapq<20" 
        # check for any long insertions
        aligned_pairs = record.get_aligned_pairs(matches_only=True)
        hard_clip = get_hard_start(record.cigarstring)
        read_seq = read_seqs[record.query_name]
        for idx in range(0, len(aligned_pairs) - 1):
            read_gap = aligned_pairs[idx + 1][0] - aligned_pairs[idx][0]
            ref_gap = aligned_pairs[idx + 1][1] - aligned_pairs[idx][1]
            if read_gap >= args.min_detected_inclusion_length and ref_gap <= max_ref_gap_at_candidate:
                ref_start = aligned_pairs[idx][1]
                ref_end = aligned_pairs[idx+1][1]
                read_start = hard_clip+aligned_pairs[idx][0]
                read_end = hard_clip+aligned_pairs[idx+1][0]
                insertion_sequence = ""
                if read_seq is not None:
                    insertion_sequence = read_seq[read_start:read_end]
                orientation = 0
                if record.is_reverse:
                    orientation = 1
                annotation = read_annotation
                if not ((abs(ref_start - record.reference_start) > int(args.min_flank_size)) and 
                    (abs(record.reference_end - ref_end) > int(args.min_flank_size))):
                    annotation = update_annotation(annotation, "flank_size")
                if not read_gap >= args.min_insertion_length:
                    annotation = update_annotation(annotation, "min_insertion_length")
                # convert read start and end positions if the alignment is rc
                if record.is_reverse and insertion_sequence != "":
                    tmp = len(read_seq) - read_start
                    read_start = len(read_seq) - read_end
                    read_end = tmp
                records_to_output[record.query_name].append(["%s\t%d\t%d\t%s\t%d\t%d\t%.1f\t%s" % (record.reference_name, ref_start, ref_end, record.query_name, read_start, read_end, orientation, insertion_sequence), annotation, read_start, read_end])
                #if record.query_name == "h1tg000024l":
                #    print("%s\t%d\t%d\t%s\t%d\t%d\t%.1f\t%d\t%d\t%d\t%d\t%s\t%s" % (record.reference_name, ref_start, ref_end, record.query_name, read_start, read_end, orientation, record.query_alignment_start, record.query_alignment_end, record.reference_start, record.reference_end, record.cigarstring[:10], record.cigarstring[len(record.cigarstring)-11:]), file=sys.stderr)
        if record.query_name in records_to_output:
            # Have at least one insertion for this read
            # Parse throught the record and find any deletions > 100 bp
            deletion_positions = get_deletion_pos(record.cigarstring)
            if len(deletion_positions) == 0:
                continue
            for i in range(len(records_to_output[record.query_name])):
                annotation = records_to_output[record.query_name][i][1]
                line_arr = records_to_output[record.query_name][i][0].split('\t')
                for del_pos in deletion_positions:
                    if read_seq is not None:
                        if record.is_reverse:
                            tmp = len(read_seq) - int(line_arr[4])
                            line_arr[4] = str(len(read_seq) - int(line_arr[5]))
                            line_arr[5] = str(tmp)
                        if ((abs(del_pos[0] - int(line_arr[4])) < 2*del_pos[1] or abs(del_pos[0] - int(line_arr[5])) < 2*del_pos[1]) and 
                            (int(line_arr[5]) - int(line_arr[4]))*0.75 < del_pos[1]):
                            # update the annotation to indicate a possible nearby deletion => mapping artifact
                            annotation = update_annotation(annotation, "deletion_possible_mapping_artifact")
                            records_to_output[record.query_name][i][1] = annotation
    for read in records_to_output:
        for item in records_to_output[read]:
            annotation = item[1]
            nearby = mapped_count[read][item[2]:item[3]]
            if len(nearby) > 1:
                annotation = update_annotation(annotation, "multiply_mapped")
            print(item[0]+"\t"+annotation)


parser = argparse.ArgumentParser( description='Extract reads with long insertions')
parser.add_argument('--fasta', type=str, required=True)
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--min-insertion-length', type=int, default=100)
parser.add_argument('--min-detected-inclusion-length', type=int, default=50)
parser.add_argument('--min-flank-size', required=False, default=100)
parser.add_argument('--min-mapq', required=False, type=int, default=20)
parser.add_argument('--reference-gap-minimum', type=int, default=100)
args = parser.parse_args()

print("\t".join(["chromosome", "reference_insertion_start", "reference_insertion_end", "read_name", "read_insertion_start", "read_insertion_end", "dust_score", "insertion_sequence","pass_fail"]))
read_seqs = {}
with pysam.FastxFile(args.fasta) as fh:
    for entry in fh:
        read_seqs[entry.name] = entry.sequence

#print("Sequences", file=sys.stderr)
get_inserts_for_bam(args, args.bam, read_seqs)
