import os
import glob


rule all:
    input:
        expand("{s}/augmented_graph/{s}.gaf", s=config["samples"]),
        #expand("{s}/graph_aligner/{s}.gaf", s=config["samples"]),
        expand("{s}/sniffles/{s}.sniffles.vcf", s=config["samples"]),
        #"combined_normals/combined.sniffles.vcf"

rule all_translocations:
    input:
        expand("{s}/sniffles/{s}.sniffles.vcf", s=config["samples"]),
        expand("{s}/hprc_minigraph/{s}.translocations.tsv", s=config["samples"]),
        expand("{s}/augmented_graph_rgfa/{s}.translocations.tsv", s=config["samples"])

include: "rules/colo.smk"
