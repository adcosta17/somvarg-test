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
    if wildcards.length == "long":
        return config["assembly_folder"]+"/"+wildcards.sample+"/nanopore/"+wildcards.reads+".fastq.gz"
    else:
        return config["assembly_folder"]+"/"+wildcards.sample+"/nanopore/split/"+wildcards.reads+".fastq.gz"

def get_test_fastq(wildcards):
    if wildcards.length == "long":
        return config["assembly_folder"]+"/"+wildcards.sample+"/nanopore/test.fastq.gz"
    else:
        return config["assembly_folder"]+"/"+wildcards.sample+"/nanopore/split/test.fastq.gz"

def get_control_fastq(wildcards):
    if wildcards.length == "long":
        return config["assembly_folder"]+"/"+wildcards.sample+"/nanopore/control.fastq.gz"
    else:
        return config["assembly_folder"]+"/"+wildcards.sample+"/nanopore/split/control.fastq.gz"

def get_assembly(wildcards):
    return config["assembly_folder"]+"/"+wildcards.sample+"/assembly/"+wildcards.sample+"."+wildcards.hap+".fa"

def get_output_bam_list(wildcards):
    bam_list = config['base_dir']+wildcards.sample+"/reads_realign/"+wildcards.sample+"."+wildcards.reads+"."+wildcards.length+".tmp.bam"
    return bam_list

rule align_assembly:
    output:
        "{sample}/assembly_mapped/{sample}.{hap}.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_reference_base,
        fasta=get_assembly
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax asm5  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """


rule get_subset_assemblies:
    input:
        bams=expand("{{sample}}/assembly_mapped/{{sample}}.{hap}.sorted.bam", hap=config["haps"]),
        bam_bai=expand("{{sample}}/assembly_mapped/{{sample}}.{hap}.sorted.bam.bai", hap=config["haps"])
    output:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa"
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_aligned_asssembly.py"),
        chr_list=get_chrom_list,
        sample=get_sample,
        input_path=get_base_dir
    threads: 1
    shell:
        """
        python {params.script} --chrom-list {params.chr_list} --sample {params.sample} --input-path {params.input_path}
        """


rule map_subset_to_ref:
    input:
        fa="{sample}/assembly_subset/{sample}.{hap}.subset.fa"
    output:
        "{sample}/assembly_subset/{sample}.{hap}_subset_to_ref.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_reference_base
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm5 {params.ref_to_use} {input.fa} | samtools sort -o {output}
        """

rule map_reads_to_subset:
    input:
        fa="{sample}/assembly_subset/{sample}.{hap}.subset.fa",
    output:
        first="{sample}/reads_assembly_subset_mapped/{sample}.{hap}.sorted.{reads}.{length}.bam"
    params:
        memory_per_thread="10G",
        fastq=get_fastq
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.fa} {params.fastq} | samtools sort -o {output.first}
        """



rule get_pat_and_mat_ref_insertions:
    input:
        mat_bam="{sample}/assembly_subset/{sample}.mat_subset_to_ref.bam",
        mat_bam_index="{sample}/assembly_subset/{sample}.mat_subset_to_ref.bam.bai",
        mat_fa="{sample}/assembly_subset/{sample}.mat.subset.fa",
        pat_bam="{sample}/assembly_subset/{sample}.pat_subset_to_ref.bam",
        pat_bam_index="{sample}/assembly_subset/{sample}.pat_subset_to_ref.bam.bai",
        pat_fa="{sample}/assembly_subset/{sample}.pat.subset.fa",
    output:
        mat="{sample}/assembly_subset/{sample}.mat.indels.tsv",
        pat="{sample}/assembly_subset/{sample}.pat.indels.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_indels.py"),
        memory_per_thread="36G"
    shell:
        """
        python {params.candidate_insertion_script} --bam {input.mat_bam} --fasta {input.mat_fa} --min-insertion-length 100 --min-mapq 20 --min-detected-inclusion-length 100 > {output.mat}
        python {params.candidate_insertion_script} --bam {input.pat_bam} --fasta {input.pat_fa} --min-insertion-length 100 --min-mapq 20 --min-detected-inclusion-length 100 > {output.pat}
        """

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



rule get_read_annotation:
    input:
        mat_inserts="{sample}/assembly_subset/{sample}.mat.indels.tsv",
        pat_inserts="{sample}/assembly_subset/{sample}.pat.indels.tsv",
        reads_mat_control_bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.control.{length}.bam",
        reads_mat_test_bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.test.{length}.bam",
        reads_pat_control_bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.control.{length}.bam",
        reads_pat_test_bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.test.{length}.bam",
        reads_mat_control_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.control.{length}.bam.bai",
        reads_mat_test_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.test.{length}.bam.bai",
        reads_pat_control_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.control.{length}.bam.bai",
        reads_pat_test_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.test.{length}.bam.bai"
    output:
        tsv="{sample}/read_annotation/{sample}.annotation.{length}.tsv",
    threads: 1
    params:
        script=srcdir("../scripts/get_truth_data.py"),
        memory_per_thread="64G",
        control_fastq=get_control_fastq,
        test_fastq=get_test_fastq,
        sample=get_sample
    shell:
        """
        python {params.script} --input-fastq-control {params.control_fastq} --input-fastq-test {params.test_fastq} --mat-ref-inserts {input.mat_inserts} --pat-ref-inserts {input.pat_inserts} --mat-bam-control {input.reads_mat_control_bam} --pat-bam-control {input.reads_pat_control_bam} --mat-bam-test {input.reads_mat_test_bam} --pat-bam-test {input.reads_pat_test_bam} --centromeres {config[centromere_filter]} > {output.tsv}
        """


rule map_reads_to_ref:
    output:
        bam="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam"
    params:
        memory_per_thread="10G",
        fastq=get_fastq,
        ref_to_use=get_graph_ref_base
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {params.ref_to_use} {params.fastq} | samtools sort -o {output.bam}
        """

rule extract_inserts:
    input:
        bam="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam",
        bam_index="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam.bai"
    output:
        tsv="{sample}/reads_realign/{sample}.{reads}.{length}.tsv",
        merged_reads="{sample}/reads_realign/{sample}.{reads}.{length}.merged_reads.txt"
    threads: 20
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
 

rule realign_inserts:
    input:
        bam="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam",
        tsv="{sample}/reads_realign/{sample}.{reads}.{length}.tsv"
    output:
        bam="{sample}/combined_realign_winnow_{reads}.{length}/{chrom}.sorted.{reads}.{length}.bam",
        tsv="{sample}/combined_realign_winnow_{reads}.{length}/{chrom}.realigned.{reads}.{length}.tsv"
    threads: 10
    params:
        memory_per_thread="30G",
        fastq=get_fastq,
        base_dir=get_base_dir,
        ref=get_reference_base
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.base_dir}/{input.bam} --tsv-list {params.base_dir}/{input.tsv} --fastq-list {params.fastq} --output-dir {params.base_dir}/{wildcards.sample}/combined_realign_winnow_{wildcards.reads}.{wildcards.length} --tsv-prefix realigned.{wildcards.reads}.{wildcards.length} --bam-prefix sorted.{wildcards.reads}.{wildcards.length} --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 8000 --max-depth 50 --chromosome-list {wildcards.chrom} --only-realign
        cd {params.base_dir}
        """

rule merge_realign_inserts:
    input:
        bam= "{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam",
        bai= "{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam.bai",
        cbams = expand("{{sample}}/combined_realign_winnow_{{reads}}.{{length}}/{chrom}.sorted.{{reads}}.{{length}}.bam",chrom=config["chroms"]),
        cbams_index = expand("{{sample}}/combined_realign_winnow_{{reads}}.{{length}}/{chrom}.sorted.{{reads}}.{{length}}.bam.bai",chrom=config["chroms"]),
        tsv = expand("{{sample}}/combined_realign_winnow_{{reads}}.{{length}}/{chrom}.realigned.{{reads}}.{{length}}.tsv",chrom=config["chroms"])
    output:
        bam="{sample}/reads_realign/{sample}.{reads}.{length}.tmp.bam"
        #tsv="combined_realign_winnow_all/all.realigned.tsv"
    threads: 1
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        out_bam_list=get_output_bam_list,
        ref=get_reference_base
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py merge --output-dir {params.base_dir}/{wildcards.sample}/combined_realign_winnow_{wildcards.reads}.{wildcards.length} --bam-list {params.base_dir}/{input.bam} --tsv-prefix realigned.{wildcards.reads}.{wildcards.length} --bam-prefix sorted.{wildcards.reads}.{wildcards.length}.tmp --output-bam-list {params.out_bam_list}
        cd {params.base_dir}
        """

rule sort_realign_bams:
    input:
        bam="{sample}/reads_realign/{sample}.{reads}.{length}.tmp.bam"
    output:
        bam="{sample}/reads_realign_sorted/{sample}.{reads}.{length}.bam"
    threads: 1
    params:
        memory_per_thread="24G"
    shell:
        """
        samtools sort {input.bam} > {output.bam}
        rm {input.bam}
        """


rule realign_haplotypes_only:
    input:
        bam="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam",
        tsv="{sample}/reads_realign/{sample}.{reads}.{length}.tsv"
    output:
        tsv="{sample}/combined_realign_haps_only_{reads}.{length}/{chrom}.realigned.{reads}.{length}.tsv"
    threads: 10
    params:
        memory_per_thread="20G",
        fastq=get_fastq,
        base_dir=get_base_dir,
        ref=get_reference_base
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.base_dir}/{input.bam} --tsv-list {params.base_dir}/{input.tsv} --fastq-list {params.fastq} --output-dir {params.base_dir}/{wildcards.sample}/combined_realign_haps_only_{wildcards.reads}.{wildcards.length} --tsv-prefix realigned.{wildcards.reads}.{wildcards.length} --bam-prefix sorted.{wildcards.reads}.{wildcards.length} --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 10000 --max-depth 500 --chromosome-list {wildcards.chrom} --max-haps 2 --only-haps
        cd {params.base_dir}
        """

rule merge_realign_haplotypes_only:
    input:
        bam= "{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam",
        bai= "{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam.bai",
        tsv = expand("{{sample}}/combined_realign_haps_only_{{reads}}.{{length}}/{chrom}.realigned.{{reads}}.{{length}}.tsv",chrom=config["chroms"])
    output:
        tsv="{sample}/combined_realign_haps_only_/{sample}.{reads}.{length}.tsv"
    threads: 1
    params:
        memory_per_thread="24G",
    shell:
        """
        cat {input.tsv} >> {output.tsv}.tmp
        sort -h -k1,2n -t '\t' {output.tsv}.tmp > {output.tsv}
        rm {output.tsv}.tmp
        """