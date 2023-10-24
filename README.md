# somvarg-test

This repository contains scripts and snakemake rules to evaluate SomvarG

The following tools are required. For some tools the path to the install folder must be specified in the ```project_config.yaml``` file. See the config file and/or the sections below to see which tools are required for which analysis. 

- [SomvarG](https://github.com/adcosta17/somvarg)
- [somrit](https://github.com/adcosta17/somrit)
- [sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- [cuteSV](https://github.com/tjiangHIT/cuteSV)
- python and matplotlib
- samtools
- bgzip
- [minigraph](https://github.com/lh3/minigraph)
- [minimap2](https://github.com/lh3/minimap2)

## Running Tests

For each evaluation there is an assosiated snakemake rule to run. 
Each evaluation expects that a fastq file be provided in a sample specific subfolder within the working directory, eg: <sample>/fastq/<sample>.fastq.gz, and be bgziped

### MLL Evaluation 

To run the evaluation of the MLL samples the following rule can be run.

```
snakemake -s run_mls.smk --configfile project_config.yaml all_translocations
```

### COLO829 Evaluation 

To run the evaluation of the COLO829 tumor and COLO829BL normal samples the following rule can be run.

```
snakemake -s run_colo.smk --configfile project_config.yaml all_translocations
```

### HPRC Baseline Evaluation 

To run the evaluation of the baseline HPRC evaluation the following rule can be run. This rule assumes that there are two samples, one control and one test both replicates from the same HPRC normal sample with no assumed geninune somatic events present.

```
snakemake -s run_hprc_graph_test.smk --configfile project_config.yaml all
```

### BT16 Graph vs Diploid Filtering Evaluation 

To run the evaluation of the 11 BT16 samples using an augmented graph for polymorphic filtering, the read to graph alignments can be generated as follows:
Note this requires the reads used to generate the assembly be noted as the normal while the 11 samples noted as the tumor in the config file. 

```
snakemake -s run_ddts.smk --configfile project_config.yaml all
```

Then to identify which insertions are polymorphic the following script can be used that takes in the calls made previously by the RTD-Diploid pipeline 
```
python scripts/identify_graph_inserts.py --gfa combined_normals/augmented_graph_rgfa/combined.gfa --gaf-folder augmented_graph_rgfa --samples <11_sample_names_csv> --tsv-folder <RTD-Diploid_per_sample_tsv_folder>
```
