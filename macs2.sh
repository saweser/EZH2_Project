#!/bin/bash

slurms="/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2/slurm"

for file in `cut -f1 /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2/samples.txt`;
do

# MAKING THE HEADER
echo '#!/bin/bash' > $slurms/$file.macs2.sh
echo '#SBATCH --workdir='$slurms >> $slurms/$file.macs2.sh
echo '#SBATCH --error='$file'.%J.err' >>$slurms/$file.macs2.sh
echo '#SBATCH --output='$file'.%J.out' >> $slurms/$file.macs2.sh
echo '#SBATCH --cpus-per-task=3' >> $slurms/$file.macs2.sh
echo '#SBATCH --mem=1000' >> $slurms/$file.macs2.sh

echo "/opt/anaconda2/bin/macs2 callpeak \
-f BAM \
-t /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/"$file".bam \
-c /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/ENCFF392XRJ_control.bam \
--tempdir /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2 \
--outdir /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2 \
--broad \
--broad-cutoff 0.1 \
-g hs \
-n "$file" \
-B \
-q 0.01 \
--keep-dup 1" >> $slurms/$file.macs2.sh

sbatch $slurms/"$file".macs2.sh
done

# run the _model.r script for every sample (generated while running macs2)
# Rscript "$file"_model.r
