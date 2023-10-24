#
# Mapping related rules
#

def get_sample(wildcards):
    return wildcards.sample

def get_chrom_list(wildcards):
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def get_base_dir(wildcards):
    return config["base_dir"]

def get_reference_base(wildcards):
    return config["reference_all"]

def get_graph_ref_base(wildcards):
    return config["reference_graph_base"]

def get_fastq(wildcards):
    return config["sample_data_folder"]+"nanopore/"+wildcards.sample+".fastq.gz"

def get_hap1(wildcards):
    return config["sample_data_folder"]+"assembly/"+config["full_sample"]+".mat.fa"

def get_hap2(wildcards):
    return config["sample_data_folder"]+"assembly/"+config["full_sample"]+".pat.fa"

#def get_truth_translocations(wildcards):
#    return config["sample_analysis_folder"]+"/"+wildcards.sample+"/fastq/"+wildcards.sample+".fastq.gz"

## Index a bam
rule make_bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "samtools index {input}"


rule map_reads_to_ref:
    output:
        bam="{sample}/reads_mapped/{sample}.bam"
    params:
        memory_per_thread="10G",
        fastq=get_fastq,
        ref_to_use=get_graph_ref_base
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax map-ont {params.ref_to_use} {params.fastq} | samtools sort -o {output.bam}
        """

rule map_hap1:
    output:
        bam="{sample}/hap1/{sample}.bam"
    params:
        memory_per_thread="10G",
        fastq=get_fastq,
        ref_to_use=get_hap1
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax map-ont {params.ref_to_use} {params.fastq} | samtools sort -o {output.bam}
        """


rule map_hap2:
    output:
        bam="{sample}/hap2/{sample}.bam"
    params:
        memory_per_thread="10G",
        fastq=get_fastq,
        ref_to_use=get_hap2
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax map-ont {params.ref_to_use} {params.fastq} | samtools sort -o {output.bam}
        """


rule extract_inserts:
    input:
        bam="{sample}/reads_mapped/{sample}.bam",
        bam_index="{sample}/reads_mapped/{sample}.bam.bai"
    output:
        tsv="{sample}/reads_realign/{sample}.tsv",
        merged_reads="{sample}/reads_realign/{sample}.merged_reads.txt"
    threads: 10
    params:
        memory_per_thread="15G",
        fastq=get_fastq,
        base_dir=get_base_dir
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-tsv {params.base_dir}/{output.tsv} --output-merged {params.base_dir}/{output.merged_reads} --fastq-file {params.fastq} --threads {threads} --min-read-len 1000
        cd {params.base_dir}
        """
 

rule realign_haplotypes_only:
    input:
        bam="{sample}/reads_mapped/{sample}.bam",
        tsv="{sample}/reads_realign/{sample}.tsv"
    output:
        tsv="{sample}/combined_realign_haps_only/{chrom}.realigned.tsv"
    threads: 10
    params:
        memory_per_thread="20G",
        fastq=get_fastq,
        base_dir=get_base_dir,
        ref=get_reference_base
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.base_dir}/{input.bam} --tsv-list {params.base_dir}/{input.tsv} --fastq-list {params.fastq} --output-dir {params.base_dir}/{wildcards.sample}/combined_realign_haps_only --tsv-prefix realigned --bam-prefix sorted --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 10000 --max-depth 500 --chromosome-list {wildcards.chrom} --max-haps 3 --only-haps
        cd {params.base_dir}
        """

rule merge_realign_haplotypes_only:
    input:
        tsv = expand("{sample}/combined_realign_haps_only/{chrom}.realigned.tsv",chrom=config["chroms"],sample=config["normals"])
    output:
        tsv="combined_normals_haps_only/combined.tsv"
    threads: 1
    params:
        memory_per_thread="96G",
    shell:
        """
        cat {input.tsv} >> {output.tsv}.tmp
        sort -h -k1,2n -t '\t' {output.tsv}.tmp > {output.tsv}
        rm {output.tsv}.tmp
        """
# Graph Generation Rules

rule generate_augmented_graph_rgfa:
    input:
        tsv="combined_normals_haps_only/combined.tsv"
    output:
        gfa="combined_normals/augmented_graph_rgfa/combined.gfa",
    threads: 1
    params:
        memory_per_thread="36G", 
        ref=get_graph_ref_base,
        script=srcdir("../scripts/augment_graph_rgfa.py")
    benchmark:
        "benchmarks/augmented_graph_rgfa.txt"
    shell:
        """
        python {params.script} --variants {input.tsv} --ref {params.ref} --format tsv > {output.gfa}
        """

rule generate_hg38_graph:
    output:
        "hg38_graph/hg38.gfa"
    threads: 10
    params:
        memory_per_thread="36G",
        ref=get_graph_ref_base
    shell:
        """
        {config[minigraph]} -cxggs -t10 {params.ref} {params.ref} > {output}
        """ 

rule align_hg38_rgfa:
    input:
        "hg38_graph/hg38.gfa"
    output:
        "{sample}/hg38_rgfa/{sample}.gaf"
    threads: 10
    params:
        memory_per_thread="6G",
        fastq=get_fastq
    shell:
        """
        {config[minigraph]} -c -x lr -t {threads} {input} {params.fastq} > {output}
        """


rule align_augmented_graph_rgfa:
    input:
        "combined_normals/augmented_graph_rgfa/combined.gfa",
    output:
        "{sample}/augmented_graph_rgfa/{sample}.gaf"
    threads: 10
    params:
        memory_per_thread="6G",
        fastq=get_fastq
    benchmark:
        "benchmarks/{sample}/augmented_graph_rgfa/{sample}.txt"
    shell:
        """
        {config[minigraph]} -c -x lr -t {threads} {input} {params.fastq} > {output}
        """

rule generate_sniffles_vcf:
    input:
        bam="{sample}/reads_mapped/{sample}.bam",
        bam_index="{sample}/reads_mapped/{sample}.bam.bai"
    output:
        vcf="{sample}/sniffles/{sample}.sniffles.vcf",
        snf="{sample}/sniffles/{sample}.sniffles.snf"
    threads: 10
    params:
        memory_per_thread="10G",
        ref=get_graph_ref_base
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --snf {output.snf} --output-rnames
        """

rule sniffles_hap1:
    input:
        bam="{sample}/hap1/{sample}.bam",
        bam_index="{sample}/hap1/{sample}.bam.bai"
    output:
        vcf="{sample}/hap1/{sample}.sniffles.vcf",
        snf="{sample}/hap1/{sample}.sniffles.snf"
    threads: 10
    params:
        memory_per_thread="10G",
        ref=get_graph_ref_base
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --snf {output.snf} --output-rnames
        """

rule sniffles_hap2:
    input:
        bam="{sample}/hap2/{sample}.bam",
        bam_index="{sample}/hap2/{sample}.bam.bai"
    output:
        vcf="{sample}/hap2/{sample}.sniffles.vcf",
        snf="{sample}/hap2/{sample}.sniffles.snf"
    threads: 10
    params:
        memory_per_thread="10G",
        ref=get_graph_ref_base
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --snf {output.snf} --output-rnames
        """


rule align_hprc_minigraph:
    output:
        "{sample}/hprc_minigraph/{sample}.gaf"
    threads: 10
    params:
        memory_per_thread="6G",
        fastq=get_fastq
    shell:
        """
        {config[minigraph]} -c -x lr -t {threads} {config[hprc_minigraph]} {params.fastq} > {output}
        """

rule call_translocations_hg38_rgfa_graph:
    input:
        gaf="{sample}/hg38_rgfa/{sample}.gaf",
        gfa="hg38_graph/hg38.gfa"
    output:
        "{sample}/hg38_rgfa/{sample}.translocations.tsv"
    threads: 1
    params:
        memory_per_thread="96G",
        script=srcdir("../scripts/identify_translocations_rgfa.py"),                                                                                                    ref=get_graph_ref_base,
    shell:
        """
        python {params.script} --gaf {input.gaf} --gfa {input.gfa} > {output}
        """


rule call_translocations_augmented_rgfa_graph:
    input:
        gaf="{sample}/augmented_graph_rgfa/{sample}.gaf",
        gfa="combined_normals/augmented_graph_rgfa/combined.gfa"
    output:
        "{sample}/augmented_graph_rgfa/{sample}.translocations.tsv"
    threads: 1
    params:
        memory_per_thread="96G",
        script=srcdir("../scripts/identify_translocations_rgfa.py"),                                                                                                    ref=get_graph_ref_base,
    shell:
        """
        python {params.script} --gaf {input.gaf} --gfa {input.gfa} > {output}
        """



rule call_translocations_hprc_graphs:
    input:
        hprc_minigraph="{sample}/hprc_minigraph/{sample}.gaf"
    output:
        hprc_minigraph="{sample}/hprc_minigraph/{sample}.translocations.tsv"
    threads: 1
    params:
        memory_per_thread="36G",
        script=srcdir("../scripts/identify_translocations_rgfa.py"),
        ref=get_graph_ref_base,
    shell:
        """
        python {params.script} --gaf {input.hprc_minigraph} --gfa {config[hprc_minigraph]} > {output.hprc_minigraph}
        """


