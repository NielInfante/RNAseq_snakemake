# Snakefile, config and scripts to run an RNA seq analysis

<table style="width:100%">
  <tr>
    <th>Firstname</th>
    <th>Lastname</th>
    <th>Age</th>
  </tr>
  <tr>
    <td>Jill</td>
    <td>Smith</td>
    <td>50</td>
  </tr>
  <tr>
    <td>Eve</td>
    <td>Jackson this line is very long. Will it wrap automatically, or will I have to put manual line breaks in? I'm not sure, so that's why I'm testing things, to see if it wlll wrap automattically or not..  WExtra characters are not a bad thing here. More is better.</td>
    <td>94</td>
  </tr>
</table>


To use this pipeline, create a new repo using this one as a template (click the Use Template button).

### Set up Environment

Now, you need to install all the needed software. The easiest way to do this is to install [miniconda](https://docs.conda.io/en/latest/miniconda.html), and run using the --use-conda option. To create a minimal environment, which is just enough to run Snakemake, use

```{sh}
conda env create -f envs/snake.yaml
```

This will create a conda environment called snake, if you want a different name, change it in the yaml file before creating. The you can activate it by:
```{sh}
conda activate snake
```
The environment snake_no_conda should have everything needed to execute all the steps, but it has not been tested extensively, and on some systems it will not install at all.


### Set up Config File and Metadata

#### You will need to create a metadata file for the samples.
The metadata file should be a tab delimited file with a line for each sample. There should be a column the gives the sample ID which matches the sample ID in the configuration file.


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

# Note: organism_db also needs to be added to rnaseq_GO.yaml in envs
organism_db: org.Hs.eg.db
kegg_db: hsa

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



| H1 | H2 | H3 |
| ---- | ---- | --- |
| Hey| `trey this`<br>`Too` | and this |

| Format   | Tag example |
| -------- | ----------- |
| Headings | =heading1=<br>==heading2==<br>===heading3=== |
| New paragraph | A blank line starts a new paragraph |
| Source code block |  // all on one line<br> {{{ if (foo) bar else   baz }}} |


| Tables        | Are           | Cool  |
| ------------- |:-------------:| -----:|
| col 3 is      | right-aligned | $1600 |
| col 2 is      | centered      |   $12 |
| zebra stripes | are neat      |    $1 |
| <ul><li>item1</li><li>item2</li></ul>| See the list | from the first column|


For each sample, list the sample name, as well as the forward and reverse reads file. If you have a lot of files, you may find the following shell command helpful:

```{sh}
# Assuming files look like {sample_name}_R1.fastq.gz

 for f in *R1*; do nm=${f%_R*}; echo "  ${nm}:"; echo "    F: $f"; echo "    R: ${f/R1/R2}"; done

```

Note the spacing is important to get the indent levels correct. Also note that you should avoid using an integer as the sample name, this will likely cause issues.

The reference base should be in the Data folder, it should be a fasta file with all the transcripts you wish to consider. It should have a suffix of .fa.

salmon_bootstraps gives how many bootstraps perform.

AnnotationDBI organism db,
prebuilt list is [here](http://bioconductor.org/packages/release/BiocViews.html#___OrgDb)
Additional resources are at [AnnotationHub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html), or build your own using [AnnotationForge](https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html)
Kegg IDs can be found [here](https://www.genome.jp/kegg/catalog/org_list.html)


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
