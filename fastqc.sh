#!/bin/bash
#SBATCH --error=fastqc.%J.err
#SBATCH --output=fastqc.%J.out
#SBATCH --workdir=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/fastqc
#SBATCH --mem=100000



fastqc --threads 10 \
--outdir /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/fastqc \
/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/ENCFF392XRJ_control.bam \
/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/ENCFF678EGD.bam \
/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/ENCFF912PNB.bam \
