##### Overlap of peaks with RNA and Proteom

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/ChipSeq")

library(plyranges)
library(tidyverse)

peaks<-readRDS("/data/htp/A07/EZH2_Project/ChipSeq/peakAnno.rds")
ov<-read.csv("/data/htp/A07/EZH2_Project/Proteom/overlap_rna_prot.csv", row.names = 1)


peak.ov<-peaks %>% filter(gene_name%in%ov$Gene.Name)

##### peaks in promoter ############################################################################################################
peak.ov.prom <-peak.ov %>% filter(grepl("Promoter", annotation))
unique(peak.ov.prom$gene_name)# 21 have peaks in promoter
# "GNG12"   "CNN3"    "AKT3"    "CAMK1D"  "MPP7"    "PLEKHA7" "FERMT2"  "MYO5A"   "ALDH1A2" "AKAP13"  "SPECC1"  "RBMS1"   "MAGI1"  
# "PPIC"    "ABRACL"  "CALD1"   "TPD52"   "ASS1"    "AIF1L"   "PIR"     "FHL1"   

##### peaks in other regions ######################################################################################################
peak.ov.other<-peak.ov %>% filter(!grepl("Promoter", annotation))
unique(peak.ov.other$gene_name) #16

#### check for enhancer ##############################################################################################################
enhancer<-read_xlsx(sheet = 1, "/data/htp/A07/EZH2_Project/ChipSeq/GeneHancer_version_4-4.xlsx") %>%
  select(1,4,5,7,2,3,6,8,9) 
colnames(enhancer)[1]<-"seqnames"
enhancer$strand<-rep('*', nrow(enhancer))
# separate the attributes column
enhancer<-enhancer %>% separate(attributes, sep = ";", into = c("genehancer_id", "connected_gene1", "score1", "connected_gene2","score2", "connected_gene3", "score3", "h", "i","j", "k", "l", "m","n","o","p","q","r","s","t","u"))%>%
  as.data.frame() %>%
  mutate(genehancer_id=gsub("genehancer_id=","", genehancer_id)) %>%
  mutate(connected_gene1=gsub("connected_gene=","", connected_gene1)) %>%
  mutate(score1=gsub("score=","", score1)) %>%
  mutate(connected_gene2=gsub("connected_gene=","", connected_gene2)) %>%
  mutate(score2=gsub("score=","", score2)) %>%
  mutate(connected_gene3=gsub("connected_gene=","", connected_gene3)) %>%
  mutate(score3=gsub("score=","", score3))
# as granges
enhancer<-enhancer %>% as_granges()

# table for peaks not in promoter
colnames(peak.ov.other)[1]<-"seqnames"
colnames(peak.ov.other)[2]<-"start"
colnames(peak.ov.other)[3]<-"end"
colnames(peak.ov.other)[4]<-"width"
colnames(peak.ov.other)[5]<-"strand"
peak.ov.other.gr<-peak.ov.other %>% as_granges()

ov.gr<-find_overlaps(enhancer, peak.ov.other.gr) %>% as.data.frame() %>% 
  select_if(~sum(!is.na(.)) > 0) %>%
  select(-source.x, -feature.name, -frame)

saveRDS(ov.gr, "/data/htp/A07/EZH2_Project/ChipSeq/overlap_enhancer_peaks_41rnaprot.rds")

a<-ov.gr %>% select(genehancer_id, connected_gene1, connected_gene2, connected_gene3, gene_name)%>%
  arrange(gene_name)
# need to check the GH ids in gene cards for the cell type used
# new GH ids have J instead of F otherwise everything the same
# use the J id in the paper
  # 7:GH07F134833, CALD1, : K562
  # 10:GH10F012365, CAMK1D, common myeloid progenitor, CD34-positive
  # 15:GH05F160163, FABP6, common myeloid progenitor, CD34-positive
  # 17:GH0XF136156, FHL1, Common myeloid progenitor CD34+
  #18:GH15F052462, MYO5A, common myeloid progenitor, CD34-positive
  #27:GH11F016775, PLEKHA7, Common myeloid progenitor CD34+
  #28:GH11F016777, PLEKHA7, Common myeloid progenitor CD34+
  #29:GH08F080078, TPD52, CD14+ monocytes, Common myeloid progenitor CD34+
  #30:GH03F023746, UBE2E1	, CD14-positive monocyte

unique(ov.gr$gene_name) #11
# 11 enhancer did overlap with peaks
# of those 9 enhancer in the right cell type, 8 different genes, two enhancer for PLEKHA7

# table of the 9 enhancer
table<-ov.gr %>% filter(genehancer_id %in% c("GH07F134833","GH10F012365","GH05F160163","GH0XF136156",
                                             "GH15F052462","GH11F016775","GH11F016777","GH08F080078",
                                             "GH03F023746"))
# do the gene name and connected gene fit?
table.fit<-table[c(3,4,5,6,9),]
saveRDS(table.fit, "/data/htp/A07/EZH2_Project/ChipSeq/enhancer.final.rds")

table.un<-table[-c(3,4,5,6,9),]
# GH10F012365	not CAMK1D	
# GH15F052462	not MYO5A
# GH11F016775 not PLEKHA7


##### Most important genes overview
# FHL1
# promoter: yes
# enhaner: yes

# UBE2E
# promoter: no
# enhancer: yes

# CA2
# not in peak list

# CNN3
# promoter: yes
# enhancer: no

# AKAP13
# Promoter: yes
# Enhancer: no

# PDK3
# not in peak list

# TPD52
# Promoter: Yes
# Enhancer: Yes

# MYO5A:
# Pormoter: Yes
# Enhancer: No

# AKT3
# Promoter: Yes
# Enhancer: No

# SPECC1
# Promoter: Yes
# Enhancer: No







