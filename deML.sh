#!/bin/bash
#SBATCH --error=deML.%J.err
#SBATCH --output=deML.%J.out
#SBATCH --workdir=/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/fastq_files/deML
#SBATCH --mem=100000





deML -i /data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/fastq_files/deML/index.txt \
-f /data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/fastq_files/lane1_R3.fastq.gz \
-if1 /data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/fastq_files/lane1_R1.fastq.gz \
-o /data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/fastq_files/deML \
-s summary.txt

# you should have more than 80% of reads successfully demultiplexed
