import os
import glob


rule all:
    input:
        expand("{s}/read_annotation/{s}.annotation.{l}.tsv", s=config["samples"], l=config["length"]),
        expand("{s}/augmented_graph/{s}.{r}.{l}.gaf", s=config["samples"], l=config["length"], r=config["read_sets"])
        #expand("{s}/graph/{s}.{r}.{l}.gaf", s=config["samples"], l=config["length"], r=config["read_sets"])

rule all_annotation:
    input:
        expand("{s}/read_annotation/{s}.annotation.{l}.tsv", s=config["samples"], l=config["length"])

rule all_assembly:
    input:
        expand("{s}/assembly_mapped/{s}.{h}.sorted.bam", s=config["samples"], h=["mat", "pat"]),
        expand("{s}/assembly_mapped/{s}.{h}.sorted.bam.bai", s=config["samples"], h=["mat", "pat"])

rule all_split:
    input:
        expand("{s}/graph_split/{s}.{c}.gaf", s=config["samples"], c=["00","01","02","03","04","05","06","07","08","09"]),
        expand("{s}/graph_split/{s}.{c}.gaf", s=config["samples"], c=list(range(10,90))),
        expand("{s}/graph_split/{s}.{c}.gaf", s=config["samples"], c=list(range(9000,9400)))

rule all_haps:
    input:
        expand("{s}/combined_realign_haps_only_/{s}.{r}.{l}.tsv", s=config["samples"], l=config["length"], r=config["read_sets"])

rule all_extract_tsv:
    input:
        expand("{s}/reads_realign/{s}.{r}.{l}.tsv",s=config["samples"], l=config["length"], r=config["read_sets"])

rule all_assembly_subset:
    input:
        expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.{r}.{l}.bam", s=config["samples"], l=config["length"], r=config["read_sets"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.{r}.{l}.bam", s=config["samples"], l=config["length"], r=config["read_sets"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.{r}.{l}.bam.bai", s=config["samples"], l=config["length"], r=config["read_sets"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.{r}.{l}.bam.bai", s=config["samples"], l=config["length"], r=config["read_sets"])

rule all_base_translocations:
    input:
        expand("{s}/augmented_graph/{s}.test.{l}.translocations.tsv", s=config["samples"], l=config["length"], r=config["read_sets"]),
        expand("{s}/reads_mapped/{s}.test.{l}.cutesv_translocations.tsv", s=config["samples"], l=config["length"], r=config["read_sets"]),
        expand("{s}/reads_mapped/{s}.test.{l}.sniffles_translocations.tsv", s=config["samples"], l=config["length"], r=config["read_sets"])

include: "rules/mapping.smk"
include: "rules/graph.smk"
