import os
import glob


rule all:
    input:
        expand("{s}/sniffles/{s}.sniffles.vcf", s=config["samples"]),
        expand("{s}/hap1/{s}.sniffles.vcf", s=config["samples"]),
        expand("{s}/hap2/{s}.sniffles.vcf", s=config["samples"]),
        expand("{s}/augmented_graph_rgfa/{s}.translocations.tsv", s=config["samples"]),
        expand("{s}/hprc_minigraph/{s}.translocations.tsv", s=config["samples"]),
        expand("{s}/hg38_rgfa/{s}.translocations.tsv", s=config["samples"])

include: "rules/hprc_test.smk"
