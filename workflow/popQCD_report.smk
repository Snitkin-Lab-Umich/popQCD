# Author: Dhatri Badri
configfile: "config/config.yaml"

import pandas as pd
import os
import json
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
import csv
from pathlib import Path

PREFIX = config["prefix"]

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)

# Calculate relative abundance 
# Relative abundance = the number of reads aligned to a reference/the total number of reads 
# in the sample. 
# Filter anything that is less than 50% breath and add the rest of the abundances 
# The number of reads aligned to a reference can be found in the genome.tsv 
# (filtered_read_pair_count), while the total number of reads can be found in the fastqc report.
# After calculating their abundance's just add them up.
# def relative_abundance(outdir, prefix):
#     prefix = prefix.pop()
#     outdir = "results/%s" % prefix
#     report_dir = str(outdir) + "/%s_Report" % prefix

#     def extract_total_sequences(fastqc_path):
#         with open(fastqc_path, 'r') as f:
#             for line in f:
#                 if line.startswith("Total Sequences"):
#                     return int(line.strip().split('\t')[-1])
#         raise ValueError(f"'Total Sequences' not found in {fastqc_path}")

#     results = []
#     inStrain_dir = os.path.join(outdir, 'inStrain') 
#     for sample in os.listdir(inStrain_dir):
#         genome_info_path = f"results/{prefix}/inStrain/{sample}/output/{sample}_genome_info.tsv"
#         fastqc_path = f"results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc/fastqc_data.txt"

#         if not os.path.exists(genome_info_path):
#             print(f"Skipping {sample}: genome info file not found.")
#             continue
#         if not os.path.exists(fastqc_path):
#             print(f"Skipping {sample}: FastQC file not found.")
#             continue

#         # Read and filter genome info
#         genome_df = pd.read_csv(genome_info_path, sep='\t')
#         genome_df_filtered = genome_df[genome_df['breadth'] > 0.5]

#         # Sum of filtered read pair counts
#         filtered_reads = genome_df_filtered['filtered_read_pair_count'].sum()

#         # Get total sequences from FastQC
#         total_sequences = extract_total_sequences(fastqc_path)

#         if total_sequences == 0:
#             relative = 0.0
#         else:
#             relative = filtered_reads / total_sequences

#         results.append({'sample': sample, 'relative_abundance': relative})

#     # Save results
#     df_results = pd.DataFrame(results)
#     results_file_path = os.path.join(report_dir, f'{prefix}_Final_Relative_Abundance.csv')  # Save final result to CSV
    
#     df_results.to_csv(results_file_path, index=False)
    # print(f"Saved relative abundance report to {output_csv}")

# def coverage:
# Coverage = no. of reads/ total genome length

# def coverage(outdir, prefix):
#     prefix = prefix.pop()
#     outdir = "results/%s" % prefix
#     report_dir = str(outdir) + "/%s_Report" % prefix

#     results = []

#     quality_aftertrim_dir = os.path.join(outdir, 'quality_aftertrim') 

#     for sample in os.listdir(quality_aftertrim_dir):
#         fastqc_path = f"results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc/fastqc_data.txt"

#         if not os.path.exists(fastqc_path):
#             print(f"Skipping {sample}: fastqc file not found.")
#             continue

#         total_sequences = None
#         total_bases = None

#         with open(fastqc_path, 'r') as f:
#             for line in f:
#                 if line.startswith("Total Sequences"):
#                     total_sequences = int(line.strip().split('\t')[-1])
#                 elif line.startswith("Total Bases"):
#                     total_bases = int(line.strip().split('\t')[-1])

#         if total_sequences is None or total_bases is None:
#             print(f"Skipping {sample} could not extract Total Sequences or Total Bases.")
#             continue

#         if total_bases == 0:
#             cov = 0.0
#         else:
#             cov = total_sequences / total_bases # 

#         results.append({'sample': sample, 'coverage': cov})

#     df_results = pd.DataFrame(results)
#     results_file_path = os.path.join(report_dir, f'{prefix}_Final_Coverage.csv') 
    
#     df_results.to_csv(results_file_path, index=False)
#     # print(f"Saved coverage report to {output_csv}")


# Kraken
# Pick out top 3 kraken species

# def kraken_report(outdir, prefix):
#     """
#     Create a report file for the the
#     """
#     prefix = prefix.pop()
#     outdir = "results/%s" % prefix
#     report_dir = str(outdir) + "/%s_Report" % prefix

#     kraken_dir = os.path.join(outdir, 'kraken') 
#     # report_files = glob.glob(os.path.join(outdir, prefix, 'kraken', 'sample', '*_kraken_report.tsv'))

#     output_data = []
#     for sample in os.listdir(kraken_dir):
#     # for report_file in report_files:
#         # sample_name = os.path.basename(report_file).replace('_kraken_report.tsv', '')
#         report_file = f"results/{prefix}/kraken/{sample}/{sample}_kraken_report.tsv"
#         top_species = []

#         with open(report_file, 'r') as f:
#             species_lines = []
#             for line in f:
#                 cols = line.strip().split('\t')
#                 if len(cols) < 5:
#                     continue
#                 # Clean up extra whitespace from species names
#                 rank_code = cols[3].strip()
#                 name = cols[5].strip()
#                 if rank_code == 'S':  # Species
#                     try:
#                         percentage = float(cols[0].strip())
#                         species_lines.append((percentage, name))
#                     except ValueError:
#                         continue

#             # Sort species by percentage in descending order
#             species_lines.sort(key=lambda x: x[0], reverse=True)

#             # Get top 3 species names
#             top_species = [s[1] for s in species_lines[:3]]
#             # Pad with 'NA' if fewer than 3 species
#             top_species += ['NA'] * (3 - len(top_species))

#         output_data.append([sample] + top_species)

#     # Write to output CSV
#     output_csv = os.path.join(report_dir, f'{prefix}_Final_Kraken_Report.csv')
#     with open(output_csv, 'w', newline='') as f:
#         writer = csv.writer(f)
#         writer.writerow(['sample', 'top_species_1', 'top_species_2', 'top_species_3'])
#         writer.writerows(output_data)

# def summary(outdir, prefix):
#     """
#     Combines Kraken, relative abundance, and coverage reports into a single QC summary.

#     Parameters:
#     - outdir (str): Output directory where final CSVs are located and summary will be saved.
#     - prefix (str): Common prefix for the input CSV files (e.g., 'final').
#     """

#     prefix = prefix.pop()
#     outdir = "results/%s" % prefix
#     report_dir = str(outdir) + "/%s_Report" % prefix
    
#     # Define file paths
#     # kraken_path = os.path.join(report_dir, f"{prefix}_Final_Kraken_Report.csv")
#     abundance_path = os.path.join(report_dir, f"{prefix}_Final_Relative_Abundance.csv")
#     coverage_path = os.path.join(report_dir, f"{prefix}_Final_Coverage.csv")
    
#     # Read CSV files
#     # kraken_df = pd.read_csv(kraken_path)
#     abundance_df = pd.read_csv(abundance_path)
#     coverage_df = pd.read_csv(coverage_path)

#     # Merge dataframes on 'sample' column
#     qc_summary_df = kraken_df.merge(abundance_df, on="sample", how="outer") \
#                              .merge(coverage_df, on="sample", how="outer")

#     # Output file
#     summary_path = os.path.join(report_dir, f"{prefix}_QC_summary.csv")
#     qc_summary_df.to_csv(summary_path, index=False)

    # print(f"QC summary saved to: {summary_path}")

# rule relative_abundance_report:
#     input:
#         outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
#         genome_info_path = lambda wildcards: expand(f"results/{wildcards.prefix}/inStrain/{wildcards.sample}/output/{wildcards.sample}_genome_info.tsv"),
#         fastqc_path = lambda wildcards: expand(f"results/{wildcards.prefix}/quality_aftertrim/{wildcards.sample}/{wildcards.sample}_Forward/{wildcards.sample}_R1_trim_paired_fastqc/fastqc_data.txt"),
#     output:
#         relative_abundance = f"results/{{prefix}}/{{prefix}}_Report/{{prefix}}_Final_Relative_Abundance.csv",
#     params:
#         prefix = "{prefix}",
#     run:
#         coverage_report({params.prefix}, {input.outdir})

# rule coverage_report:
#     input:
#         outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
#         fastqc_path = lambda wildcards: expand(f"results/{wildcards.prefix}/quality_aftertrim/{wildcards.sample}/{wildcards.sample}_Forward/{wildcards.sample}_R1_trim_paired_fastqc/fastqc_data.txt")
#     output:
#         coverage = f"results/{{prefix}}/{{prefix}}_Report/{{prefix}}_Final_Coverage.csv",
#     params:
#         prefix = "{prefix}",
#     run:
#         coverage_report({params.prefix}, {input.outdir})

# rule kraken_report:
#     input:
#         outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
#         kraken_report = lambda wildcards: expand(f"results/{wildcards.prefix}/kraken/{wildcards.sample}/{wildcards.sample}_kraken_report.tsv")
#     output:
#         kraken = f"results/{{prefix}}/{{prefix}}_Report/{{prefix}}_Final_Kraken_Report.csv",
#     params:
#         prefix = "{prefix}",
#     run:
#         kraken_report({params.prefix}, {input.outdir})

def ra_cov_species(outdir, prefix, bowtie2_db):
    prefix = prefix.pop()
    outdir = f"results/{prefix}"
    report_dir = f"{outdir}/{prefix}_Report"
    metadata_path = os.path.join(bowtie2_db, "genomes-all_metadata.tsv")

    # Load genome metadata once
    genome_metadata = pd.read_csv(metadata_path, sep='\t')
    genome_metadata['Genome'] = genome_metadata['Genome'].astype(str)

    def extract_total_sequences(fastqc_path):
        with open(fastqc_path, 'r') as f:
            for line in f:
                if line.startswith("Total Sequences"):
                    return int(line.strip().split('\t')[-1])
        raise ValueError(f"'Total Sequences' not found in {fastqc_path}")

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
        genome_df_filtered = genome_df[genome_df['breadth'] > 0.5].copy()

        if genome_df_filtered.empty:
            continue

        # Sort and take top 3 by read count
        top_hits = genome_df_filtered.sort_values(by='filtered_read_pair_count', ascending=False).head(3)

        total_sequences = extract_total_sequences(fastqc_path)

        sample_record = {'sample': sample}
        
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
        
        for rank in range(1, 4):
            species_key = f"top_species_{rank}"
            ra_key = f"{species_key}_RA"
            cov_key = f"{species_key}_coverage"

            if species_key not in sample_record:
                sample_record[species_key] = "NA"
                sample_record[ra_key] = 0.0
                sample_record[cov_key] = 0.0


        summary_records.append(sample_record)

    # Final DataFrame
    final_df = pd.DataFrame(summary_records)
    report_path = os.path.join(report_dir, f"{prefix}_QC_summary.csv")
    final_df.to_csv(report_path, index=False)
    print(f"Report saved to: {report_path}")

rule Summary:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        # relative_abundance = f"results/{{prefix}}/{{prefix}}_Report/{{prefix}}_Final_Relative_Abundance.csv",
        # coverage = f"results/{{prefix}}/{{prefix}}_Report/{{prefix}}_Final_Coverage.csv",
        # kraken = f"results/{{prefix}}/{{prefix}}_Report/{{prefix}}_Final_Kraken_Report.csv",
    output:
        QC_summary_report = f"results/{{prefix}}/{{prefix}}_Report/{{prefix}}_QC_summary.csv",
    params:
        prefix = "{prefix}",
        genome_metadata = config["bowtie2_db"]
    run:
        summary({input.outdir}, {params.prefix}, {params.genome_metadata} )

    
        

