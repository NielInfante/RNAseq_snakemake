# Snakefile, config and scripts to run an RNA seq analysis

## Readme for noConda branch

This is the setup to run whe not using the --use-conda directive to Snakemake.


To clone this branch, use the command:

```{sh}

git clone -b noConda https://github.com/NielInfante/RNAseq_snakemake.git

```
Name the resulting folder what you like, cd into it, then:
```{sh}
rm .git
git init
```

This will initialize your git project.

Now, you need to confirm all the needed software is installed. The easiest way to do this is to install miniconda, then create a conda environment to run in. I suggest:

```{sh}
conda env create -f envs/snake_no_conda
```

This will create a conda environment called snake_no_conda. If you want a different name, change it in the yaml file before creating.


You will need to create a metadata file for the samples.

You will need to create a configuration file.

Use the files meta.txt and config_small.yaml as guides.

Reads should be put in a folder called reads.

Reference transcriptome should be put in a folder called Data, and its name referenced i the config file. Note the file should have a suffix of .fa, which is not included in the config file.

Finally:

```{sh}
conda activate snake_no_conda

snakemake -n
snakemake --cores 32

conda deactivate
```

The first snakemake command is a dry run, which will build the dag and determine which jobs will need to be run. This is good practice, and helps trouble shoot issues with the config file or filenames.


### The Pipeline

![dag](dag.svg)
