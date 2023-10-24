
# Delete sequence to simulate somatic insertions

def get_sample(wildcards):
    return config["full_sample"]

def get_assembly(wildcards):
    sample = get_sample(wildcards)
    return config["assembly_folder"]+"/"+sample+"/assembly/"+sample+"."+wildcards.hap+".fa"

def get_fastq(wildcards):
    return config["fastq"]

def get_ref(wildcards):
    return config["reference"]

def max_depth(wildcards):
    return 500

def get_base_dir(wildcards):
    return config["base_dir"]

def get_repbase(wildcards):
    return config["repbase"]

def get_centromeres(wildcards):
    return config["centromeres"]

def get_telomeres(wildcards):
    return config["telomeres"]

def get_maternal(wildcards):
    sample = get_sample(wildcards)
    return config["assembly_folder"]+"/"+sample+"/assembly/"+sample+".mat.fa"

def get_paternal(wildcards):
    sample = get_sample(wildcards)
    return config["assembly_folder"]+"/"+sample+"/assembly/"+sample+".pat.fa"

def get_chrom_list(wildcards):
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def get_chrom_lengths(wildcards):
    return config["chrom_lengths"]

def get_pbsim_model(wildcards):
    return config["pbsim_model"]

def get_graph_ref_base(wildcards):
    return config["reference_graph_base"]

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



rule align_maternal:
    output:
        "HG{sample}.mat.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_maternal
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule align_paternal:
    output:
        "HG{sample}.pat.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_paternal
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule get_insert_seqs:
    input:
        mat_bam="HG{sample}.mat.sorted.bam",
        mat_bai="HG{sample}.mat.sorted.bam.bai",
        pat_bam="HG{sample}.pat.sorted.bam",
        pat_bai="HG{sample}.pat.sorted.bam.bai"
    output:
        tsv="HG{sample}/HG{sample}.insert_added_sequences.tsv", 
        mat="HG{sample}/HG{sample}.mat.assembly_subset.fa",
        pat="HG{sample}/HG{sample}.pat.assembly_subset.fa"
    threads: 1
    params:
        memory_per_thread="24G",
        repbase_seqs=get_repbase,
        mat_fasta=get_maternal,
        pat_fasta=get_paternal,
        centromeres=get_centromeres,
        chrom_lengths=get_chrom_lengths,
        script=srcdir("../scripts/get_insert_seqs.py")
    shell:
        """
        python {params.script} --input {params.repbase_seqs} --mat-bam {input.mat_bam} --pat-bam {input.pat_bam} --mat-fasta {params.mat_fasta} --pat-fasta {params.pat_fasta} --centromeres {params.centromeres} --total 500 --chrom-lengths {params.chrom_lengths} --output-folder HG{wildcards.sample} --output-prefix HG{wildcards.sample}_inserted_sequence --mat-out-fa {output.mat} --pat-out-fa {output.pat} > {output.tsv}
        """


rule simulate_seqs:
    input:
        "HG{sample}/HG{sample}.insert_added_sequences.tsv"
    output:
        tsv="HG{sample}/HG{sample}.simulated_fastq_list.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        pbsim_model=get_pbsim_model,
        script=srcdir("../scripts/generate_pbsim_run.py")
    shell:
        """

        python {params.script} --input {input} --output-folder HG{wildcards.sample} --output-prefix HG{wildcards.sample}_inserted_sequence --output-tsv {output.tsv} --pbsim-model {params.pbsim_model} --pbsim-path {config[pbsim]} > run_pbsim.sh
        chmod +x run_pbsim.sh
        ./run_pbsim.sh
        """

rule find_spanning_reads:
    input:
        tsv="HG{sample}/HG{sample}.simulated_fastq_list.tsv",
        inserts="HG{sample}/HG{sample}.insert_added_sequences.tsv"
    output:
        fastq="HG{sample}/HG{sample}.all_spanning_reads.fastq",
        tsv="HG{sample}/HG{sample}.spanning_reads_list.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_spanning_reads.py")
    shell:
        """
        python {params.script} --input {input.tsv} --inserts {input.inserts} --fastq {output.fastq} --tsv {output.tsv}
        """

rule simulate_fastq:
    input: 
        mat="HG{sample}/HG{sample}.mat.assembly_subset.fa",
        pat="HG{sample}/HG{sample}.pat.assembly_subset.fa"
    output:
        mat="HG{sample}/HG{sample}.mat_full.fastq",
        pat="HG{sample}/HG{sample}.pat_full.fastq"
    threads: 1
    params:
        memory_per_thread="64G",
        pbsim_model=get_pbsim_model
    shell:
        """
        {config[pbsim]} {input.mat} --prefix HG{wildcards.sample}.mat_full.pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 20 --hmm_model {params.pbsim_model} --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98
        rm HG{wildcards.sample}.mat_full.pbsim*.maf
        rm HG{wildcards.sample}.mat_full.pbsim*.ref
        cat HG{wildcards.sample}.mat_full.pbsim* >> HG{wildcards.sample}/HG{wildcards.sample}.mat_full.fastq
        rm HG{wildcards.sample}.mat_full.pbsim*
        {config[pbsim]} {input.pat} --prefix HG{wildcards.sample}.pat_full.pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 20 --hmm_model {params.pbsim_model} --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98
        rm HG{wildcards.sample}.pat_full.pbsim*.maf
        rm HG{wildcards.sample}.pat_full.pbsim*.ref
        cat HG{wildcards.sample}.pat_full.pbsim* >> HG{wildcards.sample}/HG{wildcards.sample}.pat_full.fastq
        rm HG{wildcards.sample}.pat_full.pbsim*
        """

rule simulate_base_fastq:
    input: 
        mat="HG{sample}/HG{sample}.mat.assembly_subset.fa",
        pat="HG{sample}/HG{sample}.pat.assembly_subset.fa"
    output:
        "HG{sample}/HG{sample}.combined.fastq"
    threads: 1
    params:
        memory_per_thread="64G",
        pbsim_model=get_pbsim_model,
        script=srcdir("../scripts/rename_fastq.py")
    shell:
        """
        {config[pbsim]} {input.mat} --prefix HG{wildcards.sample}.mat_combined.pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 20 --hmm_model {params.pbsim_model} --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98
        rm HG{wildcards.sample}.mat_combined.pbsim*.maf
        rm HG{wildcards.sample}.mat_combined.pbsim*.ref
        cat HG{wildcards.sample}.mat_combined.pbsim* >> HG{wildcards.sample}/HG{wildcards.sample}.mat_combined.fastq
        rm HG{wildcards.sample}.mat_combined.pbsim*
        {config[pbsim]} {input.pat} --prefix HG{wildcards.sample}.pat_combined.pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 20 --hmm_model {params.pbsim_model} --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98
        rm HG{wildcards.sample}.pat_combined.pbsim*.maf
        rm HG{wildcards.sample}.pat_combined.pbsim*.ref
        cat HG{wildcards.sample}.pat_combined.pbsim* >> HG{wildcards.sample}/HG{wildcards.sample}.pat_combined.fastq
        rm HG{wildcards.sample}.pat_combined.pbsim*
        cat HG{wildcards.sample}/HG{wildcards.sample}.mat_combined.fastq >> HG{wildcards.sample}/HG{wildcards.sample}.combined.tmp.fastq
        cat HG{wildcards.sample}/HG{wildcards.sample}.pat_combined.fastq >> HG{wildcards.sample}/HG{wildcards.sample}.combined.tmp.fastq
        rm HG{wildcards.sample}/HG{wildcards.sample}.mat_combined.fastq
        rm HG{wildcards.sample}/HG{wildcards.sample}.pat_combined.fastq
        python {params.script} --input HG{wildcards.sample}/HG{wildcards.sample}.combined.tmp.fastq --output-fastq {output}
        rm HG{wildcards.sample}/HG{wildcards.sample}.combined.tmp.fastq
        """

rule zip_fastq:
    input:
        "{prefix}.fastq"
    output:
        "{prefix}.fastq.gz"
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "bgzip -i {input}"


rule spike_fastqs:
    input:
        mat="HG{sample}/HG{sample}.mat_full.fastq.gz",
        pat="HG{sample}/HG{sample}.pat_full.fastq.gz",
        fastq="HG{sample}/HG{sample}.all_spanning_reads.fastq",
        tsv="HG{sample}/HG{sample}.spanning_reads_list.tsv"
    output:
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq",
        tsv="HG{sample}/HG{sample}.spike_in.{rep}.insert_list.txt"
    params:
        memory_per_thread="36G",
        script=srcdir("../scripts/spike_in_fastqs.py")
    shell:
        """
        python {params.script} --tsv {input.tsv} --inserts {input.fastq} --mat {input.mat} --pat {input.pat} --output-fastq {output.fastq} --rep {wildcards.rep} > {output.tsv}
        """


rule map_fastqs:
    input:
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq.gz"
    output:
        bam="HG{sample}/HG{sample}.spike_in.{rep}.bam"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.ref} {input.fastq} | samtools sort > {output.bam}
        """

rule map_reads_to_ref:
    input:
        "HG{sample}/HG{sample}.combined.fastq.gz"
    output:
        bam="HG{sample}/HG{sample}.combined.bam"
    params:
        memory_per_thread="10G",
        fastq=get_fastq,
        ref_to_use=get_graph_ref_base
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {params.ref_to_use} {input} | samtools sort -o {output.bam}
        """

rule extract_inserts:
    input:
        bam="HG{sample}/HG{sample}.combined.bam",
        bam_index="HG{sample}/HG{sample}.combined.bam.bai"
    output:
        tsv="HG{sample}/HG{sample}.combined.tsv",
        merged_reads="HG{sample}/HG{sample}.combined.merged_reads.txt"
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
 

rule realign_haplotypes_only:
    input:
        bam="HG{sample}/HG{sample}.combined.bam",
        tsv="HG{sample}/HG{sample}.combined.tsv"
    output:
        tsv="HG{sample}/{chrom}.HG{sample}.combined.realigned.tsv"
    threads: 10
    params:
        memory_per_thread="20G",
        fastq=get_fastq,
        base_dir=get_base_dir,
        ref=get_ref
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.base_dir}/{input.bam} --tsv-list {params.base_dir}/{input.tsv} --fastq-list {params.fastq} --output-dir {params.base_dir}/HG{wildcards.sample}/ --tsv-prefix HG{wildcards.sample}.combined.realigned --bam-prefix HG{wildcards.sample}.combined.sorted --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 10000 --max-depth 500 --chromosome-list {wildcards.chrom} --max-haps 3 --only-haps
        cd {params.base_dir}
        """

rule merge_realign_haplotypes_only:
    input:
        tsv = expand("HG{{sample}}/{chrom}.HG{{sample}}.combined.realigned.tsv",chrom=config["chroms"])
    output:
        tsv="HG{sample}/HG{sample}.combined.realigned.tsv"
    threads: 1
    params:
        memory_per_thread="24G",
    shell:
        """
        cat {input.tsv} >> {output.tsv}.tmp
        sort -h -k1,2n -t '\t' {output.tsv}.tmp > {output.tsv}
        rm {output.tsv}.tmp
        """
# Graph Generation Rules

rule generate_augmented_graph:
    input:
        tsv="HG{sample}/HG{sample}.combined.realigned.tsv"
    output:
        gfa="HG{sample}/HG{sample}.combined.gfa",
        nodes="HG{sample}/HG{sample}.combined.nodes.txt"
    threads: 1
    params:
        memory_per_thread="36G", 
        ref=get_graph_ref_base,
        script=srcdir("../scripts/augment_graph.py")
    shell:
        """
        python {params.script} --tsv {input.tsv} --ref {params.ref} --nodes {output.nodes} > {output.gfa}
        """

rule align_augmented_graph:
    input:
        gfa="HG{sample}/HG{sample}.combined.gfa",
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq.gz"
    output:
        "HG{sample}/HG{sample}.spike_in.{rep}.gaf"
    threads: 10
    params:
        memory_per_thread="6G"
    shell:
        """
        {config[minigraph]} -c -x lr -t {threads} {input.gfa} {input.fastq} > {output}
        """