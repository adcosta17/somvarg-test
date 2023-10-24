import os
import glob

rule all:
    input:
        expand("HG{s}/HG{s}.spike_in.{r}.gaf",s=config["samples"],r=config["replicates"]),
        expand("HG{s}/HG{s}.spike_in.{r}.bam",s=config["samples"],r=config["replicates"])


include: "rules/simulation.smk"

