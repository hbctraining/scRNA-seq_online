#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=64G
#SBATCH -o 02_bamtofastq_%j.out
#SBATCH -e 02_bamtofastq_%j.err  

bamtofastq --nthreads=8 01_BAM/SRR5398238.bam 02_FASTQ_SRR5398238/ --cr11
bamtofastq --nthreads=8 01_BAM/SRR5398239.bam 02_FASTQ_SRR5398239/ --cr11


