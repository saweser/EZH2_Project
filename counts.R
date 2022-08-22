##### count table ##############################################################################################################

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/RNA_Seq")
#library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)

AllCounts<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/zUMIs_output/expression/EZH2_rev2.dgecounts.rds")
anno<-read.table("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/Bria_EZH2_sample_info_bc.txt", header = T)

rownames(anno)<-anno$BC

###############################################################################################################################

##### counts: downsampled UMI counts #############################################################################################
names(AllCounts$umicount$exon$downsampling$downsampled_)

names(AllCounts$umicount$exon$all)
names(AllCounts$umicount$exon$downsampling)


counts<-AllCounts$umicount$exon$all%>%
  as.matrix() %>%
  as.data.frame()

# leave all samples in
# downsampling plot?

rownames(anno)==colnames(counts)
anno<-anno[colnames(counts),]
colnames(counts)<-anno$name

rownames(anno)<-anno$name
##########################################################################################################

##### save objects  ##################################################################################
saveRDS(anno, "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/anno.rds")
saveRDS(counts, "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/counts.rds")




