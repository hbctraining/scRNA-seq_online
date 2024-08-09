#!/bin/bash

#SBATCH --job-name=cr-count             # Job name
#SBATCH --partition=short               # Partition name
#SBATCH --time=0-06:00                  # Runtime in D-HH:MM format
#SBATCH --nodes=1                       # Number of nodes (keep at 1)
#SBATCH --ntasks=1                      # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=16              # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=64G                      # Memory needed per node (total)
#SBATCH --error=jobid_%j.err            # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out           # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                 # Type of email notification (BEGIN, END, FAIL, ALL)

module load gcc
module load cellranger/7.1.0

# Inputs for cellranger
project_name="03_ctrl"                 # Name of output
path_fastq="02_FASTQ_SRR5398238/"          # Path to folder with FASTQ files for one sample
path_ref="refdata-gex-GRCh38-2024-A"        # Path to cellranger compatible reference
local_cores=16
local_mem=64


cellranger count \
    --id=${project_name} \
    --fastqs=${path_fastq} \
    --transcriptome=${path_ref} \
    --localcores=${local_cores} \
    --localmem=${local_mem}

