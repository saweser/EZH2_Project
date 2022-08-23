#!/bin/bash

# run the phantom peakqual tool
slurms="/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/PhantomPeak/slurm"

for file in `cut -f1 /data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2/samples.txt`;do
# MAKING THE HEADER
echo '#!/bin/bash' > $slurms/$file.phantom.sh
echo '#SBATCH --workdir='$slurms >> $slurms/$file.phantom.sh
echo '#SBATCH --error='$file'.%J.err' >>$slurms/$file.phantom.sh
echo '#SBATCH --output='$file'.%J.out' >> $slurms/$file.phantom.sh
echo '#SBATCH --cpus-per-task=3' >> $slurms/$file.phantom.sh
echo '#SBATCH --mem=1000' >> $slurms/$file.phantom.sh
# run phantompeakqual
echo "Rscript /home/weser/phantompeakqualtools/run_spp.R \
-c=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/BamFiles_hg38/"$file".bam \
-out=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/PhantomPeak/"$file".txt \
-odir=/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/PhantomPeak \
-savp" >> $slurms/$file.phantom.sh
sbatch $slurms/$file.phantom.sh
done
