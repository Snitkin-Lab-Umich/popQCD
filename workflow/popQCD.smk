# Author: Dhatri Badri 
configfile: "config/config.yaml"

include: "popQCD_report.smk"

import pandas as pd
import os
import json

# Run workflow until coverage
samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])
PREFIX = config["prefix"]

if not os.path.exists(f"results/"):
    os.makedirs(f"results/")

rule all:
    input:
        QC_summary_report = expand("results/{prefix}/{prefix}_Report/{prefix}_QC_summary.csv", prefix=PREFIX),
        kraken = expand("results/{prefix}/{prefix}_Report/{prefix}_Final_Kraken_Report.csv", prefix=PREFIX),
        coverage = expand("results/{prefix}/{prefix}_Report/{prefix}_Final_Coverage.csv", prefix=PREFIX),
        relative_abundance = expand("results/{prefix}/{prefix}_Report/{prefix}_Final_Relative_Abundance.csv", prefix=PREFIX),
        fastqc_raw = expand("results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        fastqc_aftertrim = expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        # kraken_report = expand("results/{prefix}/kraken/{sample}/{sample}_kraken_report.tsv", sample=SAMPLE, prefix=PREFIX),
        # instrain_results = expand( "results/{prefix}/inStrain/{sample}/{sample}_dummy_completion_file.txt", sample=SAMPLE, prefix=PREFIX),
        instrain_results = expand("results/{prefix}/inStrain/{sample}/output/{sample}_genome_info.tsv", sample=SAMPLE, prefix=PREFIX),

rule quality_raw:
    input:
        r1 = lambda wildcards: expand(str(config["pop_data"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["pop_data"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
    output:
        raw_fastqc_report_fwd = f"results/{{prefix}}/quality_raw/{{sample}}/{{sample}}_Forward/{{sample}}_R1_fastqc.html",
        raw_fastqc_report_rev = f"results/{{prefix}}/quality_raw/{{sample}}/{{sample}}_Reverse/{{sample}}_R2_fastqc.html",
    log:
        "logs/{prefix}/quality_raw/{sample}/{sample}.log"
    params:
        outdir="results/{prefix}/quality_raw/{sample}/{sample}"
    #conda:
    #    "envs/fastqc.yaml"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    #envmodules:
    #    "Bioinformatics",
    #    "fastqc"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """

rule trimmomatic_pe:
    input:    
        r1 = lambda wildcards: expand(str(config["pop_data"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["pop_data"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
    output:
        r1 = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R1_trim_paired.fastq.gz",
        r2 = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R2_trim_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R2_trim_unpaired.fastq.gz",
    params:
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
        threads = config["ncores"],
    log:
        "logs/{prefix}/trimmomatic/{sample}/{sample}.log"
    #conda:
    #    "envs/trimmomatic.yaml"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    #envmodules:
    #    "Bioinformatics",
    #    "trimmomatic"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {params.threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}"

rule quality_aftertrim:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        aftertrim_fastqc_report_fwd = f"results/{{prefix}}/quality_aftertrim/{{sample}}/{{sample}}_Forward/{{sample}}_R1_trim_paired_fastqc.html",
        aftertrim_fastqc_report_rev = f"results/{{prefix}}/quality_aftertrim/{{sample}}/{{sample}}_Reverse/{{sample}}_R2_trim_paired_fastqc.html",
        # fastqc_stats = f"results/{{prefix}}/quality_aftertrim/{{sample}}/{{sample}}_Forward/{{sample}}_R1_trim_paired_fastqc/fastqc_data.txt",       
    log:
        "logs/{prefix}/{sample}/quality_aftertrim/{sample}.log"
    params:
        outdir="results/{prefix}/quality_aftertrim/{sample}/{sample}",
        sample = "{sample}"
    #conda:
    #    "envs/fastqc.yaml"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    #envmodules:
    #    "Bioinformatics",
    #    "fastqc"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        cd {params.outdir}_Forward/
        unzip {params.sample}_R1_trim_paired_fastqc.zip
        """

rule bowtie2:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),     
    output:
        sam_file = "results/{prefix}/bowtie2/{sample}/{sample}.sam"
    params:
        bowtie2_db = config["bowtie2_db"], 
    singularity:
        "docker://staphb/bowtie2:2.5.4"
    # envmodules:
    #    "Bioinformatics",
    #    "bowtie2/2.4.2-3pufpzz"
    shell:
        "bowtie2 -p 10 -x {params.bowtie2_db}/UHGG_reps.fasta.bt2 -1 {input.r1} -2 {input.r2} > {output.sam_file}"

# rule kraken:
#     input:
#         r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
#         r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
#     output:
#         kraken_out = f"results/{{prefix}}/kraken/{{sample}}/{{sample}}_kraken_out",
#         kraken_report = f"results/{{prefix}}/kraken/{{sample}}/{{sample}}_kraken_report.tsv",
#     params:
#         db = config["kraken_db"],
#         threads = 12
#         # threads = config["threads"]
#     #conda:
#     #    "envs/kraken.yaml"
#     singularity:
#         "docker://staphb/kraken2:2.1.3"
#     shell:
#         "kraken2 --memory-mapping  --db {params.db} --threads {params.threads} --paired --gzip-compressed --quick --output {output.kraken_out} --report {output.kraken_report} {input.r1} {input.r2}"


rule instrain_botwie2:
    input:
        sam_file = "results/{prefix}/bowtie2/{sample}/{sample}.sam"
    output:
        # instrain_profile = dir("results/{prefix}/inStrain/{sample}/"),
        instrain_results = "results/{prefix}/inStrain/{sample}/output/{sample}_genome_info.tsv",
        instrain_completion = temp("results/{prefix}/inStrain/{sample}/{sample}_dummy_completion_file.txt")
    params:
        bowtie2_db = config["bowtie2_db"], 
        out_dir = "results/{prefix}/inStrain",
        sample = "{sample}"
    conda:
       "envs/instrain.yaml"        
    shell:
        """
        inStrain profile {input.sam_file} {params.bowtie2_db}/UHGG_reps.fasta -o {params.out_dir}/{params.sample} -p 10 -g {params.bowtie2_db}/UHGG_reps.genes.fna -s {params.bowtie2_db}/UHGG.stb --database_mode
        touch {output.instrain_completion}
        """

"""
END OF PIPELINE
"""

