"""
Author: Niel Infante
Aim: Do RNA seq analysis
Date: 10 September 2019
"""


# read config info into this namespace
configfile: "config.yaml"

# Create file list
files = list()
fastqFiles = {}
for s in config['samples']:
    file = config['samples'][s]['F']

    fileName = re.sub(config['input_file_suffix'], "", file)
    files.append(file)
    fastqFiles[file] = fileName

    file = config['samples'][s]['R']
    fileName = re.sub(config['input_file_suffix'], "", file)
    files.append(file)
    fastqFiles[file] = fileName



#print(fastqFiles)

READ_FOLDER = config['reads_folder']


rule all:
    input:
        expand("QC/FastQC/{sample_file}_fastqc.zip", sample_file=fastqFiles.values()),       #FastQC output
        expand("salmon/{sample}/quant.sf", sample=config["samples"]),               # salmon quantification
        expand("results/{experiment}/report.html", experiment=config["experiments"]),  # DESeq reports
        "QC/report.html"



rule fastqc:
    input:
        lambda wildcards: f"{config['reads_folder']}/{wildcards.sample_file}{config['input_file_suffix']}"
#        f"{config['reads_folder']}/{sample_file}{config['input_file_suffix']}"
#        lambda wildcards: str(config['reads_folder'] / f"{wildcards.sample_file}")
#        expand("QC/FastQC/{sample_file}_fastqc.zip", sample_file=fastqFiles.values())
    output:
#        "QC/FastQC/{fastqFiles[sample_file]}_fastqc.zip"
        "QC/FastQC/{sample_file}_fastqc.zip"
    threads: 16
    conda:
        "envs/rnaseq_QC.yaml"
    shell:
        "fastqc -t {threads} {input} -o QC/FastQC"
#        "cat {input} | fastqc -t {threads} stdin:{sample_file} -o QC/FastQC"


rule multiQC:
    input:
#        expand("QC/FastQC/{sample_file}_fastqc.zip", sample_file=files),       #FastQC output
        expand("QC/FastQC/{sample_file}_fastqc.zip", sample_file=fastqFiles.values()),       #FastQC output
        expand("salmon/{sample}/quant.sf", sample=config["samples"])           # salmon quantification
    output:
        "QC/report.html"
    conda:
        "envs/rnaseq_QC.yaml"
    shell:
        "multiqc -f -n {output} QC/FastQC salmon"


rule salmon_build_index:
    input:
        ref=f"{config['reference_base']}.fa"
    output:
        directory(f"{config['reference_base']}")
    conda:
        "envs/rnaseq_salmon.yaml"
    shell:
        "salmon index -t {input.ref} -i {output} --type quasi -k 31"


rule salmon_quantification:
    input:
        r1 = lambda wildcards: f"{config['reads_folder']}/{config['samples'][wildcards.sample]['F']}",
        r2 = lambda wildcards: f"{config['reads_folder']}/{config['samples'][wildcards.sample]['R']}",
        index = f"{config['reference_base']}"
    output:
        quant = 'salmon/{sample}/quant.sf',
        lib = 'salmon/{sample}/lib_format_counts.json'
    log:
        'logs/salmon/{sample}.log'
    params:
        # optional parameters
        libtype ="A",
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra=f"--gcBias --validateMappings --numBootstraps {config['salmon_bootstraps']}"
    threads: 16
    conda:
        "envs/rnaseq_salmon.yaml"
    wrapper:
        "0.38.0/bio/salmon/quant"


# Parses transcript file and writes two output files with data in format for R to use
rule parseTranscripts:
    input:
        ref=f"{config['reference_base']}.fa"
    output:
        ID="Data/IDs",
        bio="Data/Biotype"
    shell:
        "perl scripts/parseFasta.pl {input.ref} {output.ID} {output.bio}"

# Rule to install packages.
# This should only have to run once.
# Actually, shouldn't need this at all switching to --use-conda
rule init_R:
    output:
        "envs/R_initialized"
    conda:
#        "envs/forRMD.yaml"
        "envs/minimal_R.yaml"
    log:
        "logs/R_init.log"
    params:
        test = "My param"
#    shell:
#        "Rscript scripts/init.R {config[organism]}"
    script:
        "scripts/init.R"


# This rule will create an R config file to be used by doDESeq.R
rule create_config_for_deseq:
    input:
#        "envs/R_initialized"      # Only include if we really want to initialize
    output:
        file="results/{experiment}/deseq/config.R"
    run:
        exp=output[0].split("/")[1]
        print(exp)

        file = open(output[0],'w')

        file.write("outPrefix <- '{}'\n".format(config['experiments'][exp]['prefix']) )
        file.write("PCA_Group <- '{}'\n".format(config['experiments'][exp]['PCA_Group']))
        file.write("design =~ {}\n".format(config['experiments'][exp]['design']))
        file.write("contrast <- {}\n\n".format(config['experiments'][exp]['contrast']))

        file.write("meta <- meta %>% filter({})\n".format(config['experiments'][exp]['filter']))
        file.write("samples <- meta${}\n".format(config['sample_column']))
        file.write("meta$Graph_Display <- meta${}\n".format(config['experiments'][exp]['display_column']))
        file.write("sample_column <- '{}'".format(config['sample_column']))

        file.close()


rule do_deseq:
    input:
        lambda wildcards: f"results/{wildcards.experiment}/deseq/config.R",
        expand("salmon/{sample}/quant.sf", sample=config['samples']),
        id="Data/IDs"
    output:
        "results/{experiment}/deseq/dds.rds",
        "results/{experiment}/deseq/results.txt"
    params:
        exp = lambda wildcards: f"{wildcards.experiment}"
    log:
        "logs/R_deseq_{experiment}.log"
    conda:
        "envs/rnaseq_deseq.yaml"
#    shell:
#        "Rscript scripts/doDESeq.R {params.exp} {config[metadata_file]} {config[tx2gene]}"
    script:
        "scripts/doDESeq.R"
        
        
rule do_GO:
    input:
        lambda wildcards: f"results/{wildcards.experiment}/deseq/dds.rds"
    output:
        "results/{experiment}/GO/all_genes_BP_results.txt",
    params:
        exp = lambda wildcards: f"{wildcards.experiment}"
    log:
        "logs/R_GO_{experiment}.log"
    conda:
        "envs/rnaseq_GO.yaml"
    script:
        "scripts/do_GO.R"
        

rule do_KEGG:
    input:
        lambda wildcards: f"results/{wildcards.experiment}/deseq/dds.rds"
    output:
        "results/{experiment}/KEGG/KEGG_results.txt",
    params:
        exp = lambda wildcards: f"{wildcards.experiment}"
    log:
        "logs/R_KEGG_{experiment}.log"
    conda:
        "envs/rnaseq_GO.yaml"
    script:
        "scripts/KEGG.R"
        

rule create_report:
    input:
        lambda wildcards: f"results/{wildcards.experiment}/deseq/dds.rds",
        lambda wildcards: f"results/{wildcards.experiment}/GO/all_genes_BP_results.txt",
        lambda wildcards: f"results/{wildcards.experiment}/KEGG/KEGG_results.txt"
    output:
        "results/{experiment}/report.html"
    params:
        exp = lambda wildcards: f"{wildcards.experiment}"
    log:
        "logs/R_report_{experiment}.log"
    conda:
        "envs/forRMD.yaml"
    script:
        "scripts/create_reports.Rmd"

rule report:
    output:
        "my_report.html"
    shell:
        "touch my_report.html"
