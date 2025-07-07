# Author: Dhatri Badri 
configfile: "config/config.yaml"

# include: "popQCD_report.smk"

import pandas as pd
import os
import json
import numpy as np

# Run workflow until coverage
samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])
PREFIX = config["prefix"]

if not os.path.exists(f"results/"):
    os.makedirs(f"results/")

def ra_cov_species(outdir, prefix):
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = f"{outdir}/{prefix}_Report"
    metadata_path = os.path.join(config["bowtie2_db"], "genomes-all_metadata.tsv")

    # Load genome metadata once
    genome_metadata = pd.read_csv(metadata_path, sep='\t')
    genome_metadata['Genome'] = genome_metadata['Genome'].astype(str)

    def extract_fastqc_metrics(fastqc_path):
        total_sequences = None
        total_bases = None
        with open(fastqc_path, 'r') as f:
            for line in f:
                if line.startswith("Total Sequences"):
                    total_sequences = int(line.strip().split('\t')[-1])
                elif line.startswith("Total Bases"):
                    bases_str = line.strip().split('\t')[-1].strip()
                    if 'Gbp' in bases_str:
                        total_bases = int(float(bases_str.replace('Gbp', '').strip()) * 1e9)
                    elif 'Mbp' in bases_str:
                        total_bases = int(float(bases_str.replace('Mbp', '').strip()) * 1e6)
                    elif 'Kbp' in bases_str:
                        total_bases = int(float(bases_str.replace('Kbp', '').strip()) * 1e3)
                    else:
                        total_bases = int(bases_str)
                if total_sequences is not None and total_bases is not None:
                    break
        if total_sequences is None:
            raise ValueError(f"'Total Sequences' not found in {fastqc_path}")
        if total_bases is None:
            raise ValueError(f"'Total Bases' not found in {fastqc_path}")
        return total_sequences, total_bases

    summary_records = []
    inStrain_dir = os.path.join(outdir, 'inStrain')
    
    for sample in os.listdir(inStrain_dir):
        genome_info_path = os.path.join(inStrain_dir, sample, 'output', f"{sample}_genome_info.tsv")
        fastqc_path = f"{outdir}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc/fastqc_data.txt"

        if not os.path.exists(genome_info_path):
            print(f"Skipping {sample}: genome info file not found.")
            continue
        if not os.path.exists(fastqc_path):
            print(f"Skipping {sample}: FastQC file not found.")
            continue

        genome_df = pd.read_csv(genome_info_path, sep='\t')
        total_sequences, total_bases = extract_fastqc_metrics(fastqc_path)

        genome_df_filtered = genome_df[genome_df['breadth'] > 0.5].copy()

        sample_record = {'sample': sample, 'Total_bases': total_bases}

        if genome_df_filtered.empty:
            # Fill NA for all species fields
            for rank in range(1, 4):
                sample_record[f"top_species_{rank}"] = "NA"
                sample_record[f"top_species_{rank}_RA"] = 0.0
                sample_record[f"top_species_{rank}_coverage"] = 0.0
            summary_records.append(sample_record)
            continue
        # Sort and take top 3 by read count
        top_hits = genome_df_filtered.sort_values(by='filtered_read_pair_count', ascending=False).head(3)

        total_sequences, total_bases = extract_fastqc_metrics(fastqc_path)

        sample_record = {'sample': sample, 'Total_bases': total_bases}
        
        for rank, (_, row) in enumerate(top_hits.iterrows(), start=1):
            genome_name = row['genome'].replace(".fna", "")
            filtered_reads = row['filtered_read_pair_count']
            coverage = row['coverage']

            # Relative abundance
            ra = filtered_reads / total_sequences if total_sequences > 0 else 0.0

            # Match species from metadata
            matched = genome_metadata[genome_metadata['Genome'] == genome_name]
            if not matched.empty:
                lineage = matched.iloc[0]['Lineage']
                species = lineage.split("s__")[-1].strip() or "NA"
            else:
                species = "NA"

            sample_record[f"top_species_{rank}"] = species
            sample_record[f"top_species_{rank}_RA"] = round(ra, 6)
            sample_record[f"top_species_{rank}_coverage"] = round(coverage, 4)
        
        # Fill NA for missing ranks
        for rank in range(1, 4):
            species_key = f"top_species_{rank}"
            ra_key = f"{species_key}_RA"
            cov_key = f"{species_key}_coverage"

            if species_key not in sample_record:
                sample_record[species_key] = "NA"
                sample_record[ra_key] = 0.0
                sample_record[cov_key] = 0.0

        summary_records.append(sample_record)


    # Final dataframe
    final_df = pd.DataFrame(summary_records)
    final_df['QC_Check'] = np.where(final_df['Total_bases'] > 10**9, 'PASS', 'FAIL')
    report_path = os.path.join(report_dir, f"{prefix}_QC_summary.csv")
    final_df.to_csv(report_path, index=False)
    print(f"Report saved to: {report_path}")

rule all:
    input:
        # summary_rule = expand("results/{prefix}/inStrain/{sample}_done.txt", prefix=PREFIX, sample=SAMPLE),
        QC_summary_report = expand("results/{prefix}/{prefix}_Report/{prefix}_QC_summary.csv", prefix=PREFIX),
        # kraken = expand("results/{prefix}/{prefix}_Report/{prefix}_Final_Kraken_Report.csv", prefix=PREFIX),
        # coverage = expand("results/{prefix}/{prefix}_Report/{prefix}_Final_Coverage.csv", prefix=PREFIX),
        # relative_abundance = expand("results/{prefix}/{prefix}_Report/{prefix}_Final_Relative_Abundance.csv", prefix=PREFIX),
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
        # instrain_completion = temp("results/{prefix}/inStrain/{sample}/{sample}_dummy_completion_file.txt")
    params:
        bowtie2_db = config["bowtie2_db"], 
        out_dir = "results/{prefix}/inStrain",
        sample = "{sample}"
    conda:
       "envs/instrain.yaml"        
    shell:
        """
        inStrain profile {input.sam_file} {params.bowtie2_db}/UHGG_reps.fasta -o {params.out_dir}/{params.sample} -p 10 -g {params.bowtie2_db}/UHGG_reps.genes.fna -s {params.bowtie2_db}/UHGG.stb --database_mode
        """


rule Summary:
    input:
        instrain_results = lambda wildcards: expand("results/{prefix}/inStrain/{sample}/output/{sample}_genome_info.tsv", sample=SAMPLE, prefix=PREFIX),
        fastqc_data = lambda wildcards: expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html",sample=SAMPLE,prefix=PREFIX)
    output:
        QC_summary_report = "results/{prefix}/{prefix}_Report/{prefix}_QC_summary.csv",
    params:
        prefix = "{prefix}",
        metadata_path = config["bowtie2_db"],
        outdir = "results/{prefix}/"
    run:
        ra_cov_species({params.outdir}, {params.prefix})

        
"""
END OF PIPELINE
"""

