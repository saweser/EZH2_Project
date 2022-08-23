#!/bin/bash
#SBATCH --workdir=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2/slurm
#SBATCH --error=ENCFF678EGD.%J.err
#SBATCH --output=ENCFF678EGD.%J.out
#SBATCH --cpus-per-task=3
#SBATCH --mem=1000
/opt/anaconda2/bin/macs2 callpeak -f BAM -t /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/ENCFF678EGD.bam -c /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/ENCFF392XRJ_control.bam --tempdir /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2 --outdir /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2 --broad --broad-cutoff 0.1 -g hs -n ENCFF678EGD -B -q 0.01 --keep-dup 1
