"""
Author: Niel Infante
Aim: Do RNA seq analysis
Date: 10 September 2019
"""


# read config info into this namespace
configfile: "config_small.yaml"
#configfile: "config.yaml"
#print (config['samples'])


READ_FOLDER = config['reads_folder']
READ_SUFFIX = config['input_file_suffix']



rule all:
    input:
        expand("FastQC/{sample}_{num}_fastqc.zip", sample=config["samples"], num=['1', '2']),

        "my_report.html",
#        expand("salmon/{sample}", sample=config["samples"]),
#        expand("deseq/{experiment}/config.R", experiment=config["experiments"]),
        expand("deseq/{experiment}/report.html", experiment=config["experiments"])


rule fastqc:
    input:
        lambda wildcards: f"{READ_FOLDER}/{config['samples'][wildcards.sample]}_{wildcards.num}{READ_SUFFIX}"
    output:
#        "FastQC/{sample}_{num,\d+}{suffix}_fastqc.zip"
        "FastQC/{sample}_{num,\d+}_fastqc.zip"
    threads: 16
    shell:
        "fastqc -t {threads} {input} -o FastQC"



rule salmon_build_index:
    input:
        ref=f"{config['reference_base']}.fa"
#    params:
#        outdir=f"{config['reference_base']}"
    output:
        directory(f"{config['reference_base']}")
    shell:
        "salmon index -t {input.ref} -i {output} --type quasi -k 31"


rule salmon_quant_reads:
    input:
        r1 = lambda wildcards: f"{READ_FOLDER}/{config['samples'][wildcards.sample]}_R1{READ_SUFFIX}",
        r2 =lambda wildcards: f"{READ_FOLDER}/{config['samples'][wildcards.sample]}_R2{READ_SUFFIX}",
        index = f"{config['reference_base']}"
    output:
        quant = 'salmon/{sample}/quant.sf',
        lib = 'salmon/{sample}/lib_format_counts.json'
    log:
        'logs/salmon3/{sample}.log'
    params:
        # optional parameters
        libtype ="A",
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra=f"--gcBias --validateMappings --numBootstraps {config['salmon_bootstraps']}"
    threads: 16
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
# This shold only have to run once.
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
rule create_cofig_for_deseq:
    input:
        "envs/R_initialized"
    output:
#        temp("deseq/{experiment}/config.R")
        file="deseq/{experiment}/config.R"
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

#        file.write("{}\n".format(config['experiments'][exp]['']))
#        file.write("{}\n".format(config['experiments'][exp]['']))
#        file.write("{}\n".format(config['experiments'][exp]['']))
#        file.write("{}\n".format(config['experiments'][exp]['']))

        file.close()


rule do_deseq:
    input:
        lambda wildcards: f"deseq/{wildcards.experiment}/config.R",
        expand("salmon/{sample}/quant.sf", sample=config["samples"]),
        id="Data/IDs"
    output:
        "deseq/{experiment}/dds.rds",
        "deseq/{experiment}/results.txt"
    params:
        exp = lambda wildcards: f"{wildcards.experiment}"
    log:
        "logs/R_deseq_{experiment}.log"
    conda:
#        "envs/forRMD.yaml"
        "envs/minimal_R.yaml"
#    shell:
#        "Rscript scripts/doDESeq.R {params.exp} {config[metadata_file]} {config[tx2gene]}"
    script:
        "scripts/doDESeq.R"

rule deseq_report:
    input:
        lambda wildcards: f"deseq/{wildcards.experiment}/dds.rds"
    output:
        "deseq/{experiment}/report.html"
    params:
        exp = lambda wildcards: f"{wildcards.experiment}"
    log:
        "logs/R_report_{experiment}.log"
    script:
        "scripts/create_reports.Rmd"

rule report:
    output:
        "my_report.html"
    shell:
        "touch my_report.html"
