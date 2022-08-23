#!/bin/bash
#SBATCH --workdir=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/PhantomPeak/slurm
#SBATCH --error=ENCFF678EGD.%J.err
#SBATCH --output=ENCFF678EGD.%J.out
#SBATCH --cpus-per-task=3
#SBATCH --mem=1000
Rscript /home/weser/phantompeakqualtools/run_spp.R -c=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/ENCFF678EGD.bam -out=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/PhantomPeak/ENCFF678EGD.txt -odir=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/PhantomPeak -savp
