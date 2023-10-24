import os
import glob


rule all:
    input:
        expand("{s}/augmented_graph_rgfa/{s}.gaf", s=config["samples"]),
        expand("{s}/hg38_graph/{s}.gaf", s=config["samples"]),
        #expand("{s}/graph_aligner/{s}.gaf", s=config["samples"]),
        #expand("{s}/sniffles/{s}.sniffles.vcf", s=config["samples"]),
        #"combined_normals/combined.sniffles.vcf"

rule all_translocations:
    input:
        expand("{s}/augmented_graph/{s}.translocations.tsv", s=config["samples"]),
        expand("{s}/hprc_minigraph/{s}.translocations.tsv", s=config["samples"]),
        expand("{s}/augmented_graph_rgfa/{s}.translocations.tsv", s=config["samples"])

rule all_norm:
    input:
        expand("{s}/augmented_graph_rgfa/{s}.normalized.tsv", s=config["samples"])

include: "rules/ddts.smk"
