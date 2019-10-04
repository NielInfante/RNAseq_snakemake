# Snakefile, config and scripts to run an RNA seq analysis


To use this pipeline, create a new repo using this one as a template (click the Use Template button).

### Set up Environment

Now, you need to confirm all the needed software is installed. The easiest way to do this is to install [miniconda](https://docs.conda.io/en/latest/miniconda.html), then create a conda environment to run in. I suggest:

```{sh}
conda env create -f envs/snake_no_conda
```

This will create a conda environment called snake_no_conda. If you want a different name, change it in the yaml file before creating.


### Set up Config and Metadata

You will need to create a metadata file for the samples.

You will need to create a configuration file.

The configuration file format is:
```
samples:
	Con_1:
		F: Con_1_1.fq.gz
		R: Con_1_2.fq.gz
	T_120d_1:
		F: T_120d_1_1.fq.gz
		R: T_120d_1_2.fq.gz

reads_folder: reads
input_file_suffix: .fq.gz

reference_base: Data/hg38_rna
salmon_bootstraps: 30
metadata_file: meta.txt
sample_column: Sample

# Differential Analysis
experiments:
	120d:
		prefix: 120d
		filter: Time == '120d' | Time == '0'
		PCA_Group: Time
		design: Time
		contrast: c('Time', '120d', '0')
		display_column: Sample

```
For each sample, list the sample name, as well as the forward and reverse reads file. If you have a lot of files, you may find the following shell command helpful:

```{sh}
# Assuming files look like {sample_name}_R1.fastq.gz

 for f in *R1*; do nm=${f%_R*}; echo "  ${nm}:"; echo "    F: $f"; echo "    R: ${f/R1/R2}"; done

```

Note the spacing is important to get the indent levels correct. Also note that you should avoid using an integer as the sample name, this will likely cause issues.

The reference base should be in the Data folder, it should be a fasta file with all the transcripts you wish to consider. It should have a suffix of .fa.

salmon_bootstraps gives how many bootstraps perform.

The metadata file should be tab delimited. There should be a column that matches the sample names listed in the config file; specify this with the sample_column parameter.

For each different contrast you want to do, list a separate experiment. Filter should be standard dplyr format, to select just the samples you want. PCA_Group specifies how samples will be colored in cerain plots. Design and contrast are input into DESeq. design can be an aribitray R formula. Display_column is if you wish to have samples labeled by something other than sample name in output figures.

Use the files meta.txt and config.yaml as guides.

Reference transcriptome should be put in a folder called Data, and its name referenced in the config file. Note the file should have a suffix of .fa, which is not included in the config file.

### To Run the Analysis

```{sh}
conda activate snake_no_conda

snakemake -n
snakemake --cores 32

conda deactivate
```

The first snakemake command is a dry run, which will build the dag and determine which jobs will need to be run. This is good practice, and helps trouble shoot issues with the config file or filenames. The second snakemake command will run the pipeline. Specify how many cores you want to use.

### The Results

A quality report will be in the QC folder.

For each contrast, a separate folder is created under deseq. In each folder, there are files for all results, significant results, and an html report with figures.


### The Pipeline

![dag](dag.svg)
