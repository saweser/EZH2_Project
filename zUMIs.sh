#!/bin/bash
#SBATCH --error=zUMIs.%J.err
#SBATCH --output=zUMIs.%J.out
#SBATCH --cpus-per-task=20
#SBATCH --workdir=/data/htp/A07/EZH2_Project/RNA_Seq/run_zUMIs

srun /home/weser/bin/zUMIs/zUMIs-master.sh -c -y /data/htp/A07/EZH2_Project/RNA_Seq/EZH2.yaml
