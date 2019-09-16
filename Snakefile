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
#        expand("FastQC/{sample}_{num}{suffix}_fastqc.zip", sample=config["samples"], num=['1', '2'], suffix=config['input_file_suffix']),
        expand("FastQC/{sample}_{num}_fastqc.zip", sample=config["samples"], num=['1', '2']),

        "report.html",
#        expand("salmon/{sample}", sample=config["samples"])
        expand("salmon/{sample}", sample=config["samples"])


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


rule report:
    output:
        "report.html"
    shell:
        "touch report.html"
