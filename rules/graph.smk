
# Graph Generation Rules

def get_ref(wildcards):
    return config["reference_all"]

def get_graph_ref_base(wildcards):
    return config["reference_graph_base"]

def get_fastq(wildcards):
    reads = "test"
    if "control" in wildcards.reads:
        reads = "control"
    if wildcards.length == "long":
        return config["assembly_folder"]+"/"+wildcards.sample+"/nanopore/"+reads+".fastq.gz"
    else:
        return config["assembly_folder"]+"/"+wildcards.sample+"/nanopore/split/"+reads+".fastq.gz"

def get_gfa(wildcards):
    return wildcards.sample+"/graph/"+wildcards.sample+".control.long.gfa"

rule generate_cute_vcf:
    input:
        bam="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam",
        bai="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam.bai",
    output:
        vcf="{sample}/cutesv.{reads}.{length}/{sample}.{reads}.{length}.cutesv.vcf"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_graph_ref_base
    shell:
        """
        cuteSV {input.bam} {params.ref} {output.vcf} {wildcards.sample}/cutesv.{wildcards.reads}.{wildcards.length} --threads {threads} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --report_readid
        """


rule generate_sniffles_vcf:
    input:
        bam="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam",
        bai="{sample}/reads_mapped/{sample}.sorted.{reads}.{length}.bam.bai",
    output:
        vcf="{sample}/sniffles.{reads}.{length}/{sample}.{reads}.{length}.sniffles.vcf"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_graph_ref_base
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --output-rnames
        """


rule alt_seq_generation:
    input:
        #cute_vcf="{sample}/reads_mapped/{sample}.control.{length}.cutesv.vcf",
        #sniffles_vcf="{sample}/reads_mapped/{sample}.control.{length}.sniffles.vcf",
        tsv="{sample}/combined_realign_haps_only_/{sample}.control.{length}.tsv"
    output:
        fa="{sample}/augmented_seqs_{length}/{sample}.control.{length}.alt_seqs.txt"
    threads: 1
    params:
        memory_per_thread="200G",
        ref=get_graph_ref_base, 
        script=srcdir("../scripts/get_svs_with_flank.py")
    shell:
        """
        python {params.script} --tsv {input.tsv} --ref {params.ref} --output-folder {wildcards.sample}/augmented_seqs_{wildcards.length} > {output.fa}
        """

rule minigraph_augmented_graph:
    input:
        fa="{sample}/augmented_seqs_{length}/{sample}.control.{length}.alt_seqs.txt"
    output:
        "{sample}/graph/{sample}.control.{length}.gfa"
    threads: 16
    params:
        memory_per_thread="5G", 
        ref=get_graph_ref_base
    shell:
        """
        {config[minigraph]} -c -x ggs -t {threads} {params.ref} {wildcards.sample}/augmented_seqs_{wildcards.length}/*.fa > {output}
        """

rule generate_augmented_graph:
    input:
        tsv="{sample}/combined_realign_haps_only_/{sample}.control.{length}.tsv"
    output:
        gfa="{sample}/augmented_graph/{sample}.control.{length}.gfa",
        nodes="{sample}/augmented_graph/{sample}.control.{length}.nodes.txt"
    threads: 1
    params:
        memory_per_thread="36G", 
        ref=get_graph_ref_base,
        script=srcdir("../scripts/augment_graph.py")
    shell:
        """
        python {params.script} --tsv {input.tsv} --ref {params.ref} --nodes {output.nodes} > {output.gfa}
        """

rule align_minigraph_graph:
    input:
        "{sample}/graph/{sample}.control.{length}.gfa"
    output:
        "{sample}/graph/{sample}.{reads}.{length}.gaf"
    threads: 10
    params:
        memory_per_thread="6G", 
        ref=get_ref,
        fastq=get_fastq
    shell:
        """
        {config[minigraph]} -c -x lr -t {threads} {input} {params.fastq} > {output}
        """

rule align_augmented_graph:
    input:
        "{sample}/augmented_graph/{sample}.control.{length}.gfa"
    output:
        "{sample}/augmented_graph/{sample}.{reads}.{length}.gaf"
    threads: 10
    params:
        memory_per_thread="6G", 
        ref=get_ref,
        fastq=get_fastq
    shell:
        """
        {config[minigraph]} -c -x lr -t {threads} {input} {params.fastq} > {output}
        """

rule call_translocations_base_ref:
    input:
        cute_vcf="{sample}/cutesv.control.long/{sample}.control.long.cutesv.vcf",
        sniffles_vcf="{sample}/sniffles.control.long/{sample}.control.long.sniffles.vcf",
        cute_test_vcf="{sample}/cutesv.test.long/{sample}.test.long.cutesv.vcf",
        sniffles_test_vcf="{sample}/sniffles.test.long/{sample}.test.long.sniffles.vcf"
    output:
        cute="{sample}/reads_mapped/{sample}.test.long.cutesv_translocations.tsv",
        sniffles="{sample}/reads_mapped/{sample}.test.long.sniffles_translocations.tsv"
    threads: 1
    params:
        memory_per_thread="12G",
        script=srcdir("../scripts/identify_translocations_ref.py")
    shell:
        """
        python {params.script} --control {input.sniffles_vcf} --test {input.sniffles_test_vcf} > {output.sniffles}
        python {params.script} --control {input.cute_vcf} --test {input.cute_test_vcf} > {output.cute}
        """

rule call_translocations_augmented_graph:
    input:
        gaf="{sample}/augmented_graph/{sample}.{reads}.{length}.gaf",
        gfa="{sample}/augmented_graph/{sample}.control.{length}.gfa",
        nodes="{sample}/augmented_graph/{sample}.control.{length}.nodes.txt"
    output:
        "{sample}/augmented_graph/{sample}.{reads}.{length}.translocations.tsv"
    threads: 1
    params:
        memory_per_thread="12G",
        script=srcdir("../scripts/identify_translocations.py")
    shell:
        """
        python {params.script} --gaf {input.gaf} --gfa {input.gfa} --nodes {input.nodes} > {output}
        """

rule zip_fastq:
    input:
        "{prefix}.fastq"
    output:
        "{prefix}.fastq.gz"
    params:
        memory_per_thread="6G"
    threads: 1
    shell:
        "bgzip -i {input}"


rule align_augmented_graph_subset:
    input:
        fastq="{sample}/fastq_split/{sample}.{count}.fastq.gz"
    output:
        "{sample}/graph_split/{sample}.{count}.gaf"
    threads: 5
    params:
        memory_per_thread="6G", 
        ref=get_ref,
        gfa=get_gfa
    shell:
        """
        {config[minigraph]} -c -x lr -t {threads} {params.gfa} {input.fastq} > {output}
        """

