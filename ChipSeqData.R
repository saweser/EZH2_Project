##### retrieve data from ENCODE ####################################################

# bam files for analysis downloaded from Encode
# Samples: https://www.encodeproject.org/experiments/ENCSR000AQE/
# Control: https://www.encodeproject.org/experiments/ENCSR000AKY/
# Files:ENCFF912PNB.bam, ENCFF678EGD.bam, ENCFF392XRJ_control.bam
# ChipSeq, K562, EZH2
# sequenced:
  # ENCFF912PNB: Illumina Genome Analyzer IIx
  # ENCFF678EGD: Illumina Genome Analyzer
  # ENCFF392XRJ: Illumina genome analyzer IIe
# mapped to hg38





##### data is also available at GEO ################################################
rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/ChipSeq")


library(GEOquery)
gseC <- getGEO("GSE95865", GSEMatrix=TRUE)
gseS <- getGEO("GSM1003576", GSEMatrix=TRUE)

show(gseC)

# get the pheno data
dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 1:3])

#if separate data table is included use this:
#df1 <- getGSEDataTables("GSE95865")
#lapply(df1, head)

# download the count data (download into working diirectory)
# control:
filePaths = getGEOSuppFiles("GSE95865")
#samples:
filePaths = getGEOSuppFiles("GSM1003576")

# control sample that comes from illumina genome analyzer
filePaths = getGEOSuppFiles("GSM2527339")

# read in the broad peak files.
a<-read.table("/data/htp/A07/EZH2_Project/ChipSeq/GSM1003576/GSM1003576_hg19_wgEncodeBroadHistoneK562Ezh239875StdPk.broadPeak.gz")
colnames(a)<-c("chrom", "chromStart", "chromENd", "name", "score","strand","signalValue", "pValue","qValue")
