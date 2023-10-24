import os
import glob


rule all:
    input:
        expand("{s}/augmented_graph/{s}.gaf", s=config["samples"]),
        expand("{s}/sniffles/{s}.sniffles.vcf", s=config["samples"]),
        #expand("{s}/cutesv/{s}.cutesv.vcf", s=config["samples"])

rule all_sniffles:
    input:
        expand("{s}/sniffles/{s}.sniffles.combined.vcf", s=config["samples"])

rule all_translocations:
    input:
        expand("{s}/hprc_minigraph/{s}.translocations.tsv", s=config["samples"]),
        expand("{s}/augmented_graph_rgfa/{s}.translocations.tsv", s=config["samples"]),
        expand("{s}/sniffles/{s}.sniffles.vcf", s=config["samples"])

include: "rules/mls.smk"
