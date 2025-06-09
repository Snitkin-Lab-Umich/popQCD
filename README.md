# popQCD
QC pipeline for population data 


## Summary
This pipeline performs the following steps:

- Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on raw sequencing reads to assess their quality before any processing.

- Trims adapters and low-quality bases from raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

- Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) again on the trimmed reads to ensure quality improvement.

- Aligns the trimmed reads to a reference genome database using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Produces SAM files containing the alignment results for each sample.


- Profiles population microdiversity and coverage using [inStrain](https://instrain.readthedocs.io/) on the Bowtie2 alignments. Generates detailed genome info tables for each sample.

<!-- - The pipeline includes (but currently commented out) rules for running [Kraken2](https://ccb.jhu.edu/software/kraken2/) to classify reads taxonomically. -->

- The main output is the **QC summary file**:  
`results/{prefix}/{prefix}_Report/{prefix}_QC_summary.csv` . This file compiles quality control metrics for all samples, providing a comprehensive overview of data quality throughout the pipeline. 


The workflow generates all the output in the prefix folder set in `config/config.yaml`. Each workflow steps gets its own individual folder as shown below. This structure provides a general view of how outputs are organized, with each tool or workflow step having its own directory. **Note that this overview does not capture all possible outputs from each tool; it only highlights the primary directories and SOME of their contents.**


## Installation 


> If you are using Great Lakes HPC, ensure you are cloning the repository in your scratch directory. Change `your_uniqname` to your uniqname. 

```

cd /scratch/esnitkin_root/esnitkin1/your_uniqname/

```

> Clone the github directory onto your system. 

```

git clone https://github.com/Snitkin-Lab-Umich/popQCD.git

```

> Ensure you have successfully cloned QCD. Type `ls` and you should see the newly created directory **_popQCD_**. Move to the newly created directory.

```

cd popQCD

```

> Load Bioinformatics, snakemake and singularity modules from Great Lakes modules.

```

module load Bioinformatics snakemake singularity mamba

```
<!--
```

module load snakemake singularity

```
-->

This workflow makes use of singularity containers available through [State Public Health Bioinformatics group](https://github.com/StaPH-B/docker-builds). If you are working on Great Lakes (umich cluster)—you can load snakemake and singularity modules as shown above. However, if you are running it on your local or other computing platform, ensure you have snakemake and singularity installed.

## Setup config, samples and cluster files

**_If you are just testing this pipeline, the config and sample files are already loaded with test data, so you do not need to make any additional changes to them. However, it is a good idea to change the prefix (name of your output folder) in the config file to give you an idea of what variables need to be modified when running your own samples on popQCD._**

### Config
As an input, the snakemake file takes a config file where you can set the path to `samples.csv`, path to your raw sequencing reads, path to adapter fasta file etc. Instructions on how to modify `config/config.yaml` is found in `config.yaml`. 

### Samples
Add samples to `config/samples.csv` following the explanation provided below. `samples.csv` should be a comma seperated file consisting of two columns—`sample_id` and `illumina_r1`.

* `sample_id` is the prefix that should be extracted from your FASTQ reads. For example, in  your raw FASTQ files directory, if you have a file called `Rush_KPC_POP_110_R1.fastq.gz`, your sample_id would be `Rush_KPC_POP_110`.

* `illumina_r1` is the name of the entire raw FASTQ file. In the same directory,  if your file is called `Rush_KPC_POP_110_R1.fastq.gz`, your sample_id would be `Rush_KPC_POP_110_R1.fastq.gz`. **_Only include forward reads._**

You can create samples.csv file using the following for loop. Replace *path_to_your_raw_reads* below with the actual path to your raw sequencing reads.

```

echo "sample_id,illumina_r1" > config/samples.csv

for read1 in path_to_your_raw_reads/*_R1.fastq.gz; do
    sample_id=$(basename $read1 | sed 's/_R1.fastq.gz//g')
    read1_basename=$(basename $read1)
    echo $sample_id,$read1_basename
done >> config/samples.csv

```

### Cluster file

Increase/reduce the walltime in `config/cluster.json` to ensure the jobs are being submitted in a timely manner. 

## Quick start

### Run popQCD on a set of samples.

> Preview the steps in popQCD by performing a dryrun of the pipeline. 

```

snakemake -s workflow/popQCD.smk --dryrun -p

```

> Run popQCD locally

```

snakemake -s workflow/popQCD.smk -p --configfile config/config.yaml --cores all

```

>Run popQCD directly on terminal (**_note: if you close your computer/shut down terminal, the pipeline will stop running. Terminal window has to be open until pipeline runs to completion._**)

```

# Load necessary modules
module load Bioinformatics snakemake singularity mamba

# Run Snakemake
snakemake -s workflow/popQCD.smk -p --use-conda --use-singularity --conda-frontend mamba -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 


```
> Submit popQCD as a batch job (**reccommended**)

Change these `SBATCH` commands: `--job-name` to a more descriptive name like run_popQCD, `--mail-user` to your email address, `--time` depending on the number of samples you have (should be more than what you specified in `cluster.json`). Feel free to make changes to the other flags if you are comfortable doing so. Once you have made the necessary changes, save the below script as `run_popQCD.sbat`. Don't forget to submit to Slurm! `sbatch bash_script_to_run_popQCD.sbat`.

```
#!/bin/bash

#SBATCH --job-name=popQCD
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dhatrib@umich.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=12:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard

# Load necessary modules
module load Bioinformatics snakemake singularity mamba

# Run Snakemake
snakemake -s workflow/popQCD.smk -p --use-conda --use-singularity --conda-frontend mamba -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 

echo Finish time: `date`

```

## DAG

## Dependencies

### Near essentials
- [Snakemake](https://snakemake.readthedocs.io/) – Workflow management system 

### Tool stack used in the workflow
 
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) 
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
- [inStrain](https://instrain.readthedocs.io/) 
<!---  - [Kraken2](https://ccb.jhu.edu/software/kraken2/) – Taxonomic sequence classification system ---> 
- [Pandas](https://pandas.pydata.org/) 
- [Python](https://www.python.org/)  

