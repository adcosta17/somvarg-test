import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pysam
from collections import defaultdict
import mappy as mp
import matplotlib.font_manager as font_manager
import matplotlib
import os.path

def get_diff_list(list1, list2):
    ret = []
    for i in range(len(list1)):
        ret.append(list2[i]-list1[i])
    return ret


parser = argparse.ArgumentParser( description='Generate plots of before and after realignment data')
parser.add_argument('--input-dir', required=True)
parser.add_argument('--samples-file', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-dir', required=True)
args = parser.parse_args()

data = {}
with open(args.samples_file, 'r') as in_data:
    count = 0
    for line in in_data:
        row = line.strip().split('\t')
        if count == 0:
            count = 1
            continue
        data[row[0]] = row 

sniffles_found = {}
sniffles_calls = defaultdict(int)
sniffles_time = defaultdict(int)
normal_graph_found = {}
normal_graph_calls = defaultdict(int)
normal_graph_time = defaultdict(int)
hprc_graph_found = {}
hprc_graph_calls = defaultdict(int)
hprc_graph_time = defaultdict(int)
max_sniffles = 0
max_graph = 0
for sample in data:
    sniffles_found[sample] = False
    normal_graph_found[sample] = False
    hprc_graph_found[sample] = False
    with open(args.input_dir+"/"+sample+"/sniffles/"+sample+".sniffles.vcf", 'r') as in_sniffles:
        for line in in_sniffles:
            row = line.strip().split('\t')
            if '#' in line:
                continue
            if "BND" in row[2]:
                sniffles_calls[sample] += 1
                if (row[0] == data[sample][2] and data[sample][3] in row[4]) or (row[0] == data[sample][3] and data[sample][2] in row[4]):
                    sniffles_found[sample] = True
    print(sample+"\t"+str(sniffles_calls[sample])+"\t"+str(sniffles_found[sample]))
    with open(args.input_dir+"/old_translocations/"+sample+".translocations.tsv", 'r') as in_graph:
        count = 0
        for line in in_graph:
            if count == 0:
                count = 1 
                continue
            row = line.strip().split('\t')
            normal_graph_calls[sample] += 1
            if (row[0] == data[sample][2] and data[sample][3] == row[2]) or (row[0] == data[sample][3] and data[sample][2] == row[2]):
                normal_graph_found[sample] = True
    print(sample+"\t"+str(normal_graph_calls[sample])+"\t"+str(normal_graph_found[sample]))
    with open(args.input_dir+"/"+sample+"/hprc_minigraph/"+sample+".translocations.tsv", 'r') as in_graph:
        count = 0
        for line in in_graph:
            if count == 0:
                count = 1 
                continue
            row = line.strip().split('\t')
            hprc_graph_calls[sample] += 1
            if (row[0] == data[sample][2] and data[sample][3] == row[2]) or (row[0] == data[sample][3] and data[sample][2] == row[2]):
                hprc_graph_found[sample] = True
    print(sample+"\t"+str(hprc_graph_calls[sample])+"\t"+str(hprc_graph_found[sample]))
    with open(args.input_dir+"/benchmarks_old/"+sample+"/sniffles/"+sample+".txt", 'r') as in_sniffles_bench:
        count = 0
        for line in in_sniffles_bench:
            if count == 0:
                count = 1 
                continue
            row = line.strip().split('\t')
            sniffles_time[sample] = int(float(row[0]))
            break
    with open(args.input_dir+"/benchmarks_old/"+sample+"/translocations/"+sample+".rgfa.txt", 'r') as in_graph_bench:
        count = 0
        for line in in_graph_bench:
            if count == 0:
                count = 1 
                continue
            row = line.strip().split('\t')
            normal_graph_time[sample] = int(float(row[0]))
            break
    with open(args.input_dir+"/benchmarks_old/"+sample+"/translocations/"+sample+".hprc.txt", 'r') as in_graph_bench:
        count = 0
        for line in in_graph_bench:
            if count == 0:
                count = 1 
                continue
            row = line.strip().split('\t')
            hprc_graph_time[sample] = int(float(row[0]))
            break

# Get the number of calls for each approach 

plot_data = {}
plot_data['X'] = []
plot_data['X1'] = []
plot_data['X2'] = []
plot_data['X3'] = []
plot_data['X4'] = []
plot_data['X5'] = []
plot_data["XString"] = []
plot_data["Coverage"] = []
plot_data["normal_graph_time"] = []
plot_data["sniffles_time"] = []
plot_data["hprc_graph_time"] = []
plot_data["normal_graph_calls"] = []
plot_data["sniffles_calls"] = []
plot_data["hprc_graph_calls"] = []
plot_data["hprc_offset_time"] = []
plot_data["hprc_offset_calls"] = []
plot_data["sniffles_offset_time"] = []
plot_data["sniffles_offset_calls"] = []
count = 1
tra_count = defaultdict(int)
for sample in data:
    if normal_graph_found[sample]:
        plot_data['X'].append(count)
        plot_data['X1'].append(count-0.25)
        plot_data['X2'].append(count+0.25)
        plot_data['X3'].append(float(count-0.375))
        plot_data['X4'].append(float(count-0.125))
        plot_data['X5'].append(float(count+0.125))
        plot_data["Coverage"].append(float(data[sample][4]))
        count += 1
        label = str(chr(63+count))#+"-"+data[sample][1]
        #if not normal_graph_found[sample]:
        #    label += '*'
        print(label)
        plot_data["XString"].append(label)
        plot_data["normal_graph_time"].append(normal_graph_time[sample])
        plot_data["hprc_offset_time"].append(normal_graph_time[sample])
        plot_data["sniffles_time"].append(sniffles_time[sample])
        plot_data["sniffles_offset_time"].append(normal_graph_time[sample]+hprc_graph_time[sample])
        plot_data["hprc_graph_time"].append(hprc_graph_time[sample])
        plot_data["normal_graph_calls"].append(normal_graph_calls[sample])
        plot_data["sniffles_calls"].append(sniffles_calls[sample])
        plot_data["hprc_graph_calls"].append(hprc_graph_calls[sample])
        plot_data["sniffles_offset_calls"].append(normal_graph_calls[sample]+hprc_graph_calls[sample])
        plot_data["hprc_offset_calls"].append(normal_graph_calls[sample])
        if sniffles_calls[sample] > max_sniffles:
            max_sniffles = sniffles_calls[sample]
        if normal_graph_calls[sample] > max_graph:
            max_graph = normal_graph_calls[sample]
        if hprc_graph_calls[sample] > max_graph:
            max_graph = hprc_graph_calls[sample]

print(max_sniffles)
print(max_graph)
print(count-1)

plt.bar(plot_data['X'],plot_data["normal_graph_time"], color='r', label="Normal Graph Time")
plt.bar(plot_data['X'],plot_data["hprc_graph_time"], bottom=plot_data["hprc_offset_time"], color='g', label="HPRC Graph Time")
plt.bar(plot_data['X'],plot_data["sniffles_time"], color='b', bottom=plot_data["sniffles_offset_time"], label="Sniffles2 Time")
plt.xticks(plot_data['X'], plot_data['XString'], fontsize=12, rotation=45)
plt.yticks(fontsize=15)
plt.title('Translocation Detection Time', fontsize=20)
plt.xlabel('Samples', fontsize=17)
plt.ylabel('Time(s)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_time.png")
plt.clf()

plt.bar(plot_data['X'],plot_data["normal_graph_calls"], color='r', label="Normal Graph Translocation Calls")
plt.bar(plot_data['X'],plot_data["hprc_graph_calls"], bottom=plot_data["hprc_offset_calls"], color='g', label="HPRC Graph Translocation Calls")
plt.bar(plot_data['X'],plot_data["sniffles_calls"], color='b', bottom=plot_data["sniffles_offset_calls"], label="Sniffles2 Translocation Calls")
plt.xticks(plot_data['X'], plot_data['XString'], fontsize=12, rotation=45)
plt.yticks(fontsize=15)
plt.title('Translocations Detected', fontsize=20)
plt.xlabel('Samples', fontsize=17)
plt.ylabel('Translocations Detected', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_calls.png")
plt.clf()


plt.bar(plot_data['X'],plot_data["normal_graph_time"], color='r', label="Normal Graph Time")
plt.bar(plot_data['X'],plot_data["sniffles_time"], color='b', bottom=plot_data["normal_graph_time"], label="Sniffles2 Time")
plt.xticks(plot_data['X'], plot_data['XString'], fontsize=12, rotation=45)
plt.yticks(fontsize=15)
plt.title('Translocation Detection Time', fontsize=20)
plt.xlabel('Samples', fontsize=17)
plt.ylabel('Time(s)', fontsize=17)
plt.yscale("symlog")
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_time_normal_only.png")
plt.clf()

plt.bar(plot_data['X'],plot_data["normal_graph_calls"], color='r', label="Normal Graph Translocation Calls")
plt.bar(plot_data['X'],plot_data["sniffles_calls"], color='b', bottom=plot_data["normal_graph_calls"], label="Sniffles2 Translocation Calls")
plt.xticks(plot_data['X'], plot_data['XString'], fontsize=12, rotation=45)
plt.yticks(fontsize=15)
plt.title('Translocations Detected', fontsize=20)
plt.xlabel('Samples', fontsize=17)
plt.ylabel('Translocations Detected', fontsize=17)
plt.yscale("symlog")
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_calls_normal_only.png")
plt.clf()


plt.bar(plot_data['X1'],plot_data["normal_graph_time"], 0.5, color='r', label="Normal Graph Time")
plt.bar(plot_data['X2'],plot_data["sniffles_time"], 0.5, color='b', label="Sniffles2 Time")
plt.xticks(plot_data['X'], plot_data['XString'], fontsize=12, rotation=45)
plt.yticks(fontsize=15)
plt.title('Translocation Detection Time', fontsize=20)
plt.xlabel('Samples', fontsize=17)
plt.ylabel('Time(s)', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_time_normal_only_group.png")
plt.clf()

plt.bar(plot_data['X1'],plot_data["normal_graph_calls"], 0.5, color='r', label="Normal Graph Translocation Calls")
plt.bar(plot_data['X2'],plot_data["sniffles_calls"], 0.5, color='b', label="Sniffles2 Translocation Calls")
plt.xticks(plot_data['X'], plot_data['XString'], fontsize=12, rotation=45)
plt.yticks(fontsize=15)
plt.title('Translocations Detected', fontsize=20)
plt.xlabel('Samples', fontsize=17)
plt.ylabel('Translocations Detected', fontsize=17)
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_calls_normal_only_group.png")
plt.clf()

plt.bar(plot_data['X1'],plot_data["normal_graph_calls"], 0.5, color='r', label="Normal Graph Translocation Calls")
plt.bar(plot_data['X2'],plot_data["sniffles_calls"], 0.5, color='b', label="Sniffles2 Translocation Calls")
plt.xticks(plot_data['X'], plot_data['XString'], fontsize=12, rotation=45)
plt.yticks(fontsize=15)
plt.title('Translocations Detected', fontsize=20)
plt.xlabel('Samples', fontsize=17)
plt.ylabel('Translocations Detected', fontsize=17)
plt.yscale("symlog")
plt.legend(loc='upper left', fontsize=15)
fig = plt.gcf()
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_calls_normal_only_group_symlog.png")
plt.clf()

####################################

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.bar(plot_data['X1'],plot_data["normal_graph_calls"], 0.5, color='r', label="Normal Graph Translocation Calls")
ax1.bar(plot_data['X2'],plot_data["sniffles_calls"], 0.5, color='b', label="Sniffles2 Translocation Calls")
ax2.plot(plot_data['X'],plot_data["Coverage"], color='g')
ax1.set_xticks(plot_data['X'])
ax1.set_xticklabels(plot_data['XString'], fontsize=12, rotation=35, ha='right')
#ax1.set_yticks(fontsize=15)
ax1.set_title('Translocations Detected', fontsize=20)
ax1.set_xlabel('Samples', fontsize=17)
ax1.set_ylabel('Translocations Detected', fontsize=17)
ax2.set_ylabel('Coverage', fontsize=17)
ax1.set_yscale("symlog")
ax1.legend(loc='upper left', fontsize=15)
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_calls_normal_only_group_symlog_cov.png")
plt.clf()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.bar(plot_data['X1'],plot_data["normal_graph_calls"], 0.5, color='r', label="Normal Graph Translocation Calls")
ax1.bar(plot_data['X2'],plot_data["sniffles_calls"], 0.5, color='b', label="Sniffles2 Translocation Calls")
ax2.plot(plot_data['X'],plot_data["Coverage"], color='g')
ax1.set_xticks(plot_data['X'])
ax1.set_xticklabels(plot_data['XString'], fontsize=20, rotation=35, ha='right')
#ax1.set_yticks(fontsize=15)
ax1.set_title('Translocations Detected', fontsize=30)
ax1.set_xlabel('Samples', fontsize=17)
ax1.set_ylabel('Translocations Detected', fontsize=20)
ax2.set_ylabel('Coverage', fontsize=17)
ax1.legend(loc='upper left', fontsize=15)
fig.set_size_inches(18.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_calls_normal_only_group_cov.png")
plt.clf()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.bar(plot_data['X3'],plot_data["normal_graph_time"], 0.25, color='green', label="SomvarG 10 Sample Graph Time", alpha=0.5)
ax1.bar(plot_data['X4'],plot_data["hprc_graph_time"], 0.25, color='orange', label="SomvarG HPRC Graph Time", alpha=0.5)
ax1.bar(plot_data['X5'],plot_data["sniffles_time"], 0.25, color='purple', label="Sniffles2 Time", alpha=0.5)
ax2.plot(plot_data['X4'],plot_data["Coverage"], 'bo')
ax2.tick_params(axis='y', colors='b',labelsize=20)
ax1.set_xticks(plot_data['X4'])
ax1.set_xticklabels(plot_data['XString'], fontsize=20, rotation=35, ha='right')
ax1.tick_params(axis='y', labelsize=20)
ax1.set_title('Translocation Detection Elapsed Time', fontsize=35)
ax1.set_xlabel('Samples', fontsize=25)
ax1.set_ylabel('Time(s)', fontsize=25)
ax2.set_ylabel('Coverage', fontsize=25)
ax2.yaxis.label.set_color('b')
ax2.set_ylim(top=45)
ax1.legend(loc='upper left', fontsize=17)
fig.set_size_inches(23.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_time_cov.png")
plt.clf()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.bar(plot_data['X3'],plot_data["normal_graph_calls"], 0.25, color='green', label="SomvarG 10 Sample Translocation Calls", alpha=0.5)
ax1.bar(plot_data['X4'],plot_data["hprc_graph_calls"], 0.25, color='orange', label="SomvarG HPRC Graph Translocation Calls", alpha=0.5)
ax1.bar(plot_data['X5'],plot_data["sniffles_calls"], 0.25, color='purple', label="Sniffles2 Translocation Calls", alpha=0.5)
ax2.plot(plot_data['X4'],plot_data["Coverage"], 'bo')
ax2.tick_params(axis='y', colors='b',labelsize=20)
ax1.set_xticks(plot_data['X4'])
ax1.set_xticklabels(plot_data['XString'], fontsize=20, rotation=35, ha='right')
ax1.tick_params(axis='y', labelsize=20)
ax1.set_title('Translocations Detected', fontsize=35)
ax1.set_xlabel('Samples', fontsize=25)
ax1.set_ylabel('Translocations Detected', fontsize=25)
ax2.set_ylabel('Coverage', fontsize=25)
ax1.set_yscale("symlog")
ax2.yaxis.label.set_color('b')
ax2.set_ylim(top=45)
ax1.legend(loc='upper left', fontsize=17)
fig.set_size_inches(23.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_calls_symlog_cov.png")
plt.clf()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.bar(plot_data['X3'],plot_data["normal_graph_calls"], 0.25, color='green', label="SomvarG 10 Sample Translocation Calls", alpha=0.5)
ax1.bar(plot_data['X4'],plot_data["hprc_graph_calls"], 0.25, color='orange', label="SomvarG HPRC Graph Translocation Calls", alpha=0.5)
ax1.bar(plot_data['X5'],plot_data["sniffles_calls"], 0.25, color='purple', label="Sniffles2 Translocation Calls", alpha=0.5)
ax2.plot(plot_data['X4'],plot_data["Coverage"], 'bo')
ax2.tick_params(axis='y', colors='b',labelsize=20)
ax1.set_xticks(plot_data['X4'])
ax1.set_xticklabels(plot_data['XString'], fontsize=20, rotation=35, ha='right')
ax1.tick_params(axis='y', labelsize=20)
ax1.set_title('Translocations Detected', fontsize=35)
ax1.set_xlabel('Samples', fontsize=25)
ax1.set_ylabel('Translocations Detected', fontsize=25)
ax2.set_ylabel('Coverage', fontsize=25)
ax2.yaxis.label.set_color('b')
ax2.set_ylim(top=45)
ax1.legend(loc='upper left', fontsize=17)
fig.set_size_inches(23.5, 12.5)
fig.savefig(args.output_dir+"/"+args.output_prefix+"_translocation_calls_cov.png")
plt.clf()