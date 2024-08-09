#!/bin/bash
#SBATCH -c 1
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH --mem=64G
#SBATCH -o SRR5398239_%j.out
#SBATCH -e SRR5398239_%j.err  


wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398239/2.2.bam.1
mv 2.2.bam.1 SRR5398239.bam
wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398238/2.1.bam.1
mv 2.1.bam.1 SRR5398238.bam
