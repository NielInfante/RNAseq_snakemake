# Snakefile, config and scripts to run an RNA seq analysis

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



| Item | Example | Comment |
| ---- | ---- | --- |
| `samples:` | `Con_1:`<br>&nbsp`; Too` | and this |


<table style="width:100%">
  <tr>
    <th>Item</th>
    <th>Example</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td><pre>samples:<pre></td>
    <td><pre>Con1:
  F: Con_1_1.fq.gz
  R: Con_1_2.fq.gz</pre></td>
    <td>Sample ID, followed by the forward and reverse file names. If you have a lot of files, you may find the following shell command helpful, assuming the files look like {sample_name}_R1.fastq.gz:<br><pre>for f in *R1*; do nm=${f%_R*};
 echo "  ${nm}:";
 echo "    F: $f";
 echo "    R: ${f/R1/R2}";
done</pre><br>Note the spacing is important to get the indent levels correct. Also note that you should avoid using an integer as the sample name, this will likely cause issues.</td>
  </tr>
  <tr>
    <td><pre>reads_folder:</pre></td>
    <td><pre>reads</pre</td>
    <td>The folder where the reads are stored</td>
  </tr>
  <tr>
		<td><pre>input_file_suffix:</pre></td>
		<td><pre>.fq.gz</pre></td>
		<td>The suffix of the read files</td>
	</tr>
	<tr>
		<td><pre>reference_base:</pre></td>
		<td><pre>Data/hg38_rna</pre></td>
		<td>Location and name of the reference transcripts. The actual filename should end with .fa. I generally get this from <a href="http://useast.ensembl.org/info/data/ftp/index.htm"l>Ensembl</a>, getting both CDS and ncRNA, then cat them into a single file.</td>
	</tr>
	<tr>
		<td><pre>salmon_bootstraps:</pre></td>
		<td><pre>30</pre></td>
		<td>Number of bootstraps in the salmon quantification. This is mainly important if you are going to do transcript level differential expression (outside this pipeline) and don't want to re-run salmon for that.</td>
	</tr>
	<tr>
		<td><pre>organism_db:</pre></td>
		<td><pre>org.Hs.eg.db</pre></td>
		<td>What organism is this? A list of prebuild annotations is <a href="http://bioconductor.org/packages/release/BiocViews.html#___OrgDb">here</a>. Additional resources are at <a href="https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html">AnnotationHub</a>, or build your own using <a href=https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html">AnnotationForge</a>.<br>This organism also needs to be added to rnaseq_</td>
	</tr>
	<tr>
		<td><pre>kegg_db:</pre></td>
		<td><pre>hsa</pre></td>
		<td>Three letter KEGG organism ID. <a href="https://www.genome.jp/kegg/catalog/org_list.html">Here</a> is a list of available IDs.</td>
	</tr>
	<tr>
		<td><pre>metadata_file:</pre></td>
		<td><pre>meta.txt</pre></td>
		<td>The metadata file. See above for a description.</td>
	</tr>
	<tr>
		<td><pre>sample_column:</pre></td>
		<td><pre>Sample</pre></td>
		<td>The name of the column in the meta data file that corresponds to the sample ID.</td>
	</tr>
	<tr>
		<td><pre>experiments:</pre></td>
		<td><pre>120d:</pre></td>
		<td>A list of every contrast to make on the data</td>
	</tr>
	<tr>
		<th colspan="3"> Experiment Parameters</th>
	</tr>
	<tr>
		<td><pre></pre></td>
		<td><pre></pre></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>prefix:</pre></td>
		<td><pre>120d</pre></td>
		<td>A repeat of the experiment name. Output folders will have this name.</td>
	</tr>
	<tr>
		<td><pre>filter:</pre></td>
		<td><pre>Time == '120d' | Time == '0'</pre></td>
		<td>Select the samples you want using standard R logic. This will be executed inside a dplyr filter command. Any columns in the metadata file can be used. If you do not want to do any filtering, simple use TRUE.</td>
	</tr>
	<tr>
		<td><pre>PCA_Group:</pre></td>
		<td><pre>Time</pre></td>
		<td>A column from the metadata file, it will be used for color groups in the PCA plot.</td>
	</tr>
	<tr>
		<td><pre>design:</pre></td>
		<td><pre>Time</pre></td>
		<td>Design of the experiment. See the DESeq2 <a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html">documentation</a> for details.</td>
	</tr>
	<tr>
		<td><pre>contrast:</pre></td>
		<td><pre>c('Time', '120d', '0')</pre></td>
		<td>The contrast to be used; what do you want compared. See the DESeq2 <a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html">documentation</a> for details.</td>
	</tr>
	<tr>
		<td><pre>display_column:</pre></td>
		<td><pre>Sample</pre></td>
		<td>A column from the metadata file that you want samples labeled with in figures.</td>
	</tr>

</table>




The reference base should be in the Data folder, it should be a fasta file with all the transcripts you wish to consider. It should have a suffix of .fa.



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
