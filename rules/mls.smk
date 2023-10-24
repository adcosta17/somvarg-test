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
    return config["sample_data_folder"]+"/"+wildcards.sample+"/fastq/"+wildcards.sample+".fastq.gz"

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
    threads: 20
    benchmark:
        "benchmarks/{sample}/reads_mapped/{sample}.txt"
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
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-tsv {params.base_dir}/{output.tsv} --output-merged {params.base_dir}/{output.merged_reads} --fastq-file {params.fastq} --threads {threads} --min-read-len 500
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

rule generate_augmented_graph:
    input:
        tsv="combined_normals_haps_only/combined.tsv"
    output:
        gfa="combined_normals/augmented_graph/combined.gfa",
        nodes="combined_normals/augmented_graph/combined.nodes.txt"
    threads: 1
    params:
        memory_per_thread="36G", 
        ref=get_graph_ref_base,
        script=srcdir("../scripts/augment_graph.py")
    benchmark:
        "benchmarks/augmented_graph.txt"
    shell:
        """
        python {params.script} --tsv {input.tsv} --ref {params.ref} --nodes {output.nodes} > {output.gfa}
        """

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
        python {params.script} --variants {input.tsv} --format tsv --ref {params.ref} > {output.gfa}
        """

rule align_augmented_graph:
    input:
        "combined_normals/augmented_graph/combined.gfa",
    output:
        "{sample}/augmented_graph/{sample}.gaf"
    threads: 10
    params:
        memory_per_thread="6G",
        fastq=get_fastq
    benchmark:
        "benchmarks/{sample}/augmented_graph/{sample}.txt"
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

rule generate_cute_vcf:
    input:
        bam="{sample}/reads_mapped/{sample}.bam",
        bai="{sample}/reads_mapped/{sample}.bam.bai"
    output:
        vcf="{sample}/cutesv/{sample}.cutesv.vcf"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_graph_ref_base
    shell:
        """
        cuteSV {input.bam} {params.ref} {output.vcf} {wildcards.sample}/reads_mapped/cutesv --threads {threads} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --report_readid
        """


rule generate_sniffles_snf:
    input:
        bam="{sample}/reads_mapped/{sample}.bam",
        bai="{sample}/reads_mapped/{sample}.bam.bai"
    output:
        vcf="{sample}/sniffles/{sample}.sniffles.vcf",
        snf="{sample}/sniffles/{sample}.sniffles.snf"
    threads: 10
    params:
        memory_per_thread="8G",
        ref=get_graph_ref_base
    benchmark:
        "benchmarks/{sample}/sniffles/{sample}.txt"
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --minsupport 3 --output-rnames --snf {output.snf}
        """

rule generate_sniffles_pop_vcf:
    input:
        snf_sample="{sample}/sniffles/{sample}.sniffles.snf",
        snf_normals=expand("{sample}/sniffles/{sample}.sniffles.snf",sample=config["normals"])
    output:
        vcf="{sample}/sniffles/{sample}.sniffles.combined.vcf"
    threads: 10
    params:
        memory_per_thread="8G",
        ref=get_graph_ref_base
    benchmark:
        "benchmarks/{sample}/sniffles/{sample}.pop.txt"
    shell:
        """
        sniffles -t {threads} -i {input.snf_sample} {input.snf_normals}  -v {output.vcf} --minsupport 3 --output-rnames --combine-output-filtered
        """


#
#rule call_translocations_base_ref:
#    input:
#        cute_vcf="{sample}/cutesv/{sample}.cutesv.vcf",
#        sniffles_vcf="{sample}/sniffles/{sample}.sniffles.vcf"
#    output:
#        cute="{sample}/reads_mapped/{sample}.test.{length}.cutesv_translocations.tsv",
#        sniffles="{sample}/reads_mapped/{sample}.test.{length}.sniffles_translocations.tsv"
#    threads: 1
#    params:
#        memory_per_thread="12G",
#        script=srcdir("../scripts/identify_translocations_ref.py")
#    shell:
#        """
#        python {params.script} --control {input.sniffles_vcf} --test {input.sniffles_test_vcf} > {output.sniffles}
#        python {params.script} --control {input.cute_vcf} --test {input.cute_test_vcf} > {output.cute}
#        """

rule call_translocations_augmented_graph:
    #input:
    #    gaf="{sample}/augmented_graph/{sample}.gaf",
    #    gfa="combined_normals/augmented_graph/combined.gfa",
    #    nodes="combined_normals/augmented_graph/combined.nodes.txt"
    output:
        "{sample}/augmented_graph/{sample}.translocations.tsv"
    threads: 1
    params:
        memory_per_thread="96G",
        script=srcdir("../scripts/identify_translocations.py"),
        ref=get_graph_ref_base,
        gaf="{sample}/augmented_graph/{sample}.gaf",
        gfa="combined_normals/augmented_graph/combined.gfa",
        nodes="combined_normals/augmented_graph/combined.nodes.txt"
    benchmark:
        "benchmarks/{sample}/translocations/{sample}.txt"
    shell:
        """
        python {params.script} --gaf {params.gaf} --gfa {params.gfa} --nodes {params.nodes} --fastq {config[sample_data_folder]}/{wildcards.sample}/fastq/{wildcards.sample}.fastq.gz --ref {params.ref} --racon /u/adcosta/SimpsonLab/racon/build/bin/racon > {output}
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
        ref=get_graph_ref_base,
    benchmark:
        "benchmarks/{sample}/translocations/{sample}.rgfa.txt"
    shell:
        """
        python config[somvarg_dir]/somvarg.py --gaf {input.gaf} --gfa {input.gfa} > {output}
        """


rule call_translocations_hprc_graphs:
    input:
        hprc_minigraph="{sample}/hprc_minigraph/{sample}.gaf"
    output:
        hprc_minigraph="{sample}/hprc_minigraph/{sample}.translocations.tsv"
    threads: 1
    params:
        memory_per_thread="36G",
        ref=get_graph_ref_base,
    benchmark:
        "benchmarks/{sample}/translocations/{sample}.hprc.txt"
    shell:
        """
        python config[somvarg_dir]/somvarg.py --gaf {input.hprc_minigraph} --gfa {config[hprc_minigraph]} > {output.hprc_minigraph}
        """

