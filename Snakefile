# Snakefile to run trimmomatic and qc ahead of dada2 amplicon pipeline
## Last updated 07-08-2020 SHu
configfile: "config.yaml"

import io 
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#

PROJ = config["proj_name"]
INPUTDIR = config["raw_data"]
SCRATCH = config["scratch"]
OUTPUTDIR = config["outputDIR"]
ADAPTERS = config["primers"]

SUF = config["suffix"]
R1_SUF = str(config["r1_suf"])
R2_SUF = str(config["r2_suf"])

# Use glob statement to find all samples in 'raw_data' directory
## Wildcard '{num}' must be equivalent to 'R1' or '1', meaning the read pair designation.
SAMPLE_LIST,NUMS = glob_wildcards(INPUTDIR + "/{sample}_{num}" + SUF)
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
SET_NUMS = set(NUMS)

#----DEFINE RULES----#

rule all:
  input:
    # fastqc output before trimming
    raw_html = expand("{scratch}/fastqc/{sample}_{num}_fastqc.html", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    raw_zip = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),

    # Cut-adapt output
    cutadapt_out = expand("{scratch}/trim-cutadapt/{sample}_{num}_trim.fastq.gz", scratch = SCRATCH, sample = SAMPLE_SET, num=SET_NUMS),
    cutadapt_qc = expand("{scratch}/trim-cutadapt/{sample}.qc.txt", scratch = SCRATCH, sample = SAMPLE_SET),
    
    # Trimmomatic output
    matic_out = expand("{scratch}/trim-matic/{sample}_{num}_trim.fastq.gz", scratch = SCRATCH, sample = SAMPLE_SET, num = SET_NUMS),

    # fastqc output after trimming primers
    trim_html = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.html", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    trim_zip = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),

	# Multi-qc output
    raw_multi_html = SCRATCH + "/fastqc/raw_multiqc.html",
    raw_multi_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
    trim_multi_html = SCRATCH + "/fastqc/trimmed_multiqc.html", #next change to include proj name
    trim_multi_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt"

rule fastqc:
  input:    
    INPUTDIR + "/{sample}_{num}.fastq.gz"
  output:
    html = SCRATCH + "/fastqc/{sample}_{num}_fastqc.html",
    zip = SCRATCH + "/fastqc/{sample}_{num}_fastqc.zip"
  params: ""
  log:
    SCRATCH + "/logs/fastqc/{sample}_{num}.log"
  wrapper:
    "0.35.2/bio/fastqc"

rule cutadapt:
    input:
        [INPUTDIR + "/{sample}_" + R1_SUF + SUF, INPUTDIR + "/{sample}_" + R2_SUF + SUF]
    output:
        fastq1 = SCRATCH + "/trim-cutadapt/{sample}_" + R1_SUF + "_trim.fastq.gz",
        fastq2 = SCRATCH + "/trim-cutadapt/{sample}_" + R2_SUF + "_trim.fastq.gz",
        qc = SCRATCH + "/trim-cutadapt/{sample}.qc.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters = "-a CCAGCASCYGCGGTAATTCC -A ACTTTCGTTCTTGATYRA",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        others = "--minimum-length 100 -q 10 --max-n 0 -e 0.4 --pair-filter=both"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 4 # set desired number of threads here
    wrapper:
        "0.35.2/bio/cutadapt/pe"

rule trimmomatic_pe:
  input:
    r1 = INPUTDIR + "/{sample}_" + R1_SUF + SUF,
    r2 = INPUTDIR + "/{sample}_" + R2_SUF + SUF
  output:
    r1 = SCRATCH + "/trim-matic/{sample}_" + R1_SUF + "_trim.fastq.gz",
    r2 = SCRATCH + "/trim-matic/{sample}_" + R2_SUF + "_trim.fastq.gz",
    r1_unpaired = SCRATCH + "/trim-matic/{sample}_1.unpaired.fastq.gz",
    r2_unpaired = SCRATCH + "/trim-matic/{sample}_2.unpaired.fastq.gz"
  log:
    SCRATCH + "/trim-matic/logs/trimmomatic/{sample}.log"
  params:
    trimmer = ["ILLUMINACLIP:{}:8:12:8".format(ADAPTERS), "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:10:2", "MINLEN:150"],
    extra = ""
  wrapper:
    "0.35.2/bio/trimmomatic/pe"

rule fastqc_trim:
  input:
    SCRATCH + "/trim-matic/{sample}_{num}_trim.fastq.gz"
  output:
    html = SCRATCH + "/fastqc/{sample}_{num}_trimmed_fastqc.html",
    zip = SCRATCH + "/fastqc/{sample}_{num}_trimmed_fastqc.zip"
  params: ""
  log:
    SCRATCH + "/logs/fastqc/{sample}_{num}_trimmed.log"
  wrapper:
    "0.35.2/bio/fastqc"

rule multiqc:
  input:
    raw_qc = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
    trim_qc = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS)
  output:
    raw_multi_html = SCRATCH + "/fastqc/raw_multiqc.html", 
    raw_multi_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
    trim_multi_html = SCRATCH + "/fastqc/trimmed_multiqc.html", 
    trim_multi_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt"
  conda:
   "envs/multiqc-env.yaml"
  shell: 
    """
    multiqc -n multiqc.html {input.raw_qc} #run multiqc
    mv multiqc.html {output.raw_multi_html} #rename html
    mv multiqc_data/multiqc_general_stats.txt {output.raw_multi_stats} #move and rename stats
    rm -rf multiqc_data #clean-up
    #repeat for trimmed data
    multiqc -n multiqc.html {input.trim_qc} #run multiqc
    mv multiqc.html {output.trim_multi_html} #rename html
    mv multiqc_data/multiqc_general_stats.txt {output.trim_multi_stats} #move and rename stats
    rm -rf multiqc_data	#clean-up
    """ 
