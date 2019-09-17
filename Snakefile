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

        "report.html",
        expand("salmon/{sample}", sample=config["samples"]),
        expand("deseq/{experiment}/config.R", experiment=config["experiments"]),
#        expand("deseq/{experiment}/dds.rds", experiment=config["experiments"])


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


rule salmon_quantification:
    input:
        forward=lambda wildcards: f"{READ_FOLDER}/{config['samples'][wildcards.sample]}_1{READ_SUFFIX}",
        reverse=lambda wildcards: f"{READ_FOLDER}/{config['samples'][wildcards.sample]}_2{READ_SUFFIX}",
        ref=f"{config['reference_base']}"
#        confirm_ref=f"{config['reference_base']}/hash.bin"
    output:
        out=directory("salmon/{sample}")

    params:
        bootstraps = {config['salmon_bootstraps']}
    threads: 16
    shell:
        "salmon quant -i {input.ref} --libType A --gcBias --numBootstraps {params.bootstraps} -p {threads} -1 {input.forward} -2 {input.reverse} -o {output.out} --validateMappings"


# This rule will create an R config file to be used by doDESeq.R
rule create_cofig_for_deseq:
    output:
#        temp("deseq/{experiment}/config.R")
        file="deseq/{experiment}/config.R"
    run:
        exp=output[0].split("/")[1]
        print(exp)

        file = open(output[0],'w')

        file.write("library('{}')\n".format(config['experiments'][exp]['organism']))
        file.write("orgDB <- {})\n\n".format(config['experiments'][exp]['organism']))

        file.write("outPrefix <- '{}'\n".format(config['experiments'][exp]['prefix']) )
        file.write("PCA_Group <- '{}'\n".format(config['experiments'][exp]['PCA_Group']))
        file.write("design =~ {}\n".format(config['experiments'][exp]['design']))
        file.write("contrast <- {}\n\n".format(config['experiments'][exp]['contrast']))

        file.write("meta <- meta %>% filter({})\n".format(config['experiments'][exp]['filter']))
        file.write("tx2gene <- read_tsv({})\n".format(config['tx2gene']))
        file.write("samples <- meta${}\n".format(config['sample_column']))
        file.write("meta$Display <- meta${}\n".format(config['experiments'][exp]['display_column']))

#        file.write("{}\n".format(config['experiments'][exp]['']))
#        file.write("{}\n".format(config['experiments'][exp]['']))
#        file.write("{}\n".format(config['experiments'][exp]['']))
#        file.write("{}\n".format(config['experiments'][exp]['']))

        file.close()


rule do_deseq:
    input:
        lambda wildcards: f"deseq/{[wildcards.experiment]}/config.R"
    output:
        "deseq/{experiment}/dds.rds",
        "deseq/{experiment}/results.txt"
    params:
        exp = lambda wildcards: f"{wildcards.experiment}"
    shell:
        "Rscript scripts/doDESeq.R {params.exp} {config[metadata_file]} {config[tx2gene]}"



rule report:
    output:
        "report.html"
    shell:
        "touch report.html"
