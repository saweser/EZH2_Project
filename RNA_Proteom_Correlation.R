#### Proteom and RNASeq

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/Proteom")

library(UpSetR)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(limma)
library(tibble)


prot<-readRDS("/data/htp/A07/EZH2_Project/Proteom/diffProtTopTable.rds")
rna<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/DE_Genes_without_WT13_KO7.rds")

prot.sig<-prot %>% filter(adj.P.Val<0.05)
rna.sig<-rna %>% filter(adj.P.Val<0.05)

### Correlation  ##########################################################################
# lfc rna expression against lfc prot expression level per gene
cor<-prot %>% 
  dplyr::select(logFC, adj.P.Val, T..PG.Genes, dir) %>%
  full_join(rna[,c(2, 8, 20, 36)], by=c("T..PG.Genes"="gene_name")) %>%
  na.omit()

colnames(cor)<-c("Prot.logFC", "Prot.adj.P.Val", "Gene.Name", "Prot.dir", "RNA.logFC", "RNA.adj.P.Val", "RNA.dir")

library(ggpubr)
p<-ggplot(data=cor, aes(x=Prot.logFC,y=RNA.logFC))
p<-p+geom_point(col="grey",alpha=0.5)+
  geom_smooth(method="lm",se=F, color="black", size=0.7)+
  geom_point(data=filter(cor, Prot.dir!="same" & RNA.dir!="same" ),
             aes(x=Prot.logFC,y=RNA.logFC),alpha=0.5)+
 # scale_color_manual(values = wes_palette("BottleRocket2",n=3))+
  theme_bw()+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  xlab("logFC of mRNA expression")+ylab("logFC of protein expression")+
  theme(legend.position = "None")+
  stat_cor(method = "pearson") +
  geom_vline(xintercept=0, color="grey30")+
  geom_hline(yintercept=0, color="grey30")

p
ggsave("correlation.pdf", plot = p, width= 120, height= 80, units= "mm", device = "pdf"
       ,path = "/data/htp/A07/EZH2_Project/Proteom/")

#### numbers ##################################################################################
cor %>% filter(Prot.dir!="same" & RNA.dir!="same") #41
cor %>% filter(Prot.dir=="up" & RNA.dir=="up") #37
cor %>% filter(Prot.dir=="down" & RNA.dir=="down") #2
# 2 in wrong direction
length(intersect(prot$T..PG.Genes, rna$gene_name)) #5139
###############################################################################################

##### table of overlap #######################################################################
ov<-cor %>% filter(Prot.dir!="same" & RNA.dir!="same") #41
ov_up<-cor %>% filter(Prot.dir=="up" & RNA.dir=="up") #37
ov_down<-cor %>% filter(Prot.dir=="down" & RNA.dir=="down") #2
ov_wrong1<-cor %>% filter(Prot.dir=="down" & RNA.dir=="up") #1
ov_wrong2<-cor %>% filter(Prot.dir=="up" & RNA.dir=="down") #1

#save
write.csv(ov_up[,c(3,5,1)], "/data/htp/A07/EZH2_Project/Proteom/overlap_rna_prot_up.csv")
write.csv(ov[,c(3,5,1)], "/data/htp/A07/EZH2_Project/Proteom/overlap_rna_prot.csv")

##### heatmap of de rna and proteom ##########################################################
library(ComplexHeatmap)
raw<-read.table("/data/htp/A07/EZH2_Project/Proteom/output_tables/Log2 values only.txt",
                 header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=F, dec=",")
raw.mean<-raw %>%
  mutate(KO1=(KO1_1 + KO1_2 +KO1_3)/3) %>%
  mutate(KO2=(KO2_1 + KO2_2 +KO2_3)/3) %>%
  mutate(KO3=(KO3_1 + KO3_2 +KO3_3)/3) %>%
  mutate(KO4=(KO4_1 + KO4_2 +KO4_3)/3) %>%
  mutate(KO5=(KO5_1 + KO5_2 +KO5_3)/3) %>%
  mutate(KO6=(KO6_1 + KO6_2 +KO6_3)/3) %>%
  dplyr::select(KO1, KO2, KO3, KO4, KO5, KO6, WT_1, WT_2, WT_3, WT_4, WT_5, T..PG.Genes, T..PG.ProteinAccessions)
colnames(raw.mean)<-c("KO_1", "KO_3", "KO_7", "KO_10", "KO_18", "KO_23", "WT_1", "WT_2", "WT_3", "WT_5", "WT_13", "T..PG.Genes", "T..PG.ProteinAccessions")


raw.ov <- raw.mean %>% filter(T..PG.Genes %in% ov$Gene.Name)%>%
  filter(T..PG.ProteinAccessions!="H7C4S7")

temp<-t(apply(as.matrix(raw.ov[,c(1:11)]),1,scale))
colnames(temp)<-colnames(raw.ov[,c(1:11)])
rownames(temp)<-ov$Gene.Name
Heatmap(temp, name="zscore",show_column_names = T, show_row_dend = F, show_column_dend = F, clustering_distance_rows = "spearman")

#### rna counts
norm.counts<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/voom.rds")$E

#norm.counts<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/norm.counts.rds")

# add gene name 
library(biomartr)
gff<-plyranges::read_gff("/data/htp/A07/EZH2_Project/RNA_Seq/anno/Homo_sapiens.GRCh38.100.chr.chrNamesUCSC.gtf") %>% 
  as.data.frame  %>%
  filter(type=="gene") # subset the gff by gene

norm.counts2<-norm.counts %>% 
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(gff, by="gene_id") %>%
  filter(gene_name %in% ov$Gene.Name) %>%
  column_to_rownames("gene_name")

temp2<-t(apply(as.matrix(norm.counts2[,c(2:10)]),1,scale))
colnames(temp2)<-colnames(norm.counts2[,c(2:10)])

#reorder rownames 
temp2<-temp2[rownames(temp),]

Heatmap(t(temp2), name="zscore",show_column_names = T, show_row_dend = F, show_column_dend = F, 
        clustering_distance_rows = "pearson")

# plot both heatmaps together
ht1 = Heatmap(t(temp), name="zscore", show_column_names = T, show_row_dend = F, 
              show_column_dend = F, clustering_distance_rows = "spearman", 
              row_title = "Proteom")
ht2 = Heatmap(t(temp2), name="zscore2", show_heatmap_legend = F, show_column_names = T, show_row_dend = F, 
              show_column_dend = F, clustering_distance_rows = "spearman", 
              row_title = "RNA_Seq")

ht_list = ht1 %v% ht2
pdf("/data/htp/A07/EZH2_Project/Proteom/Heatmap_DeGenes_Overlap.pdf", width=10, height=5)
draw(ht_list, ht_gap = unit(0.5, "cm"), merge_legend = TRUE)
dev.off()




# heatmap_legend_param = list(direction = "horizontal")
# show_heatmap_legend = F
# ht_list = ht1 + ht2
# draw(ht_list, ht_gap = unit(1, "cm"))



# spearman for log transformed
# rna counts use pearson
# z score? 
# what i do is scale, scale is the zscore





##### cancer mutatios from cosmic data base ##################################################
library(data.table)
# would not read in with read.table use fread instead
hal<- fread("/data/htp/A07/CosmicData/Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv", select = c(1:7), na.strings=c("","NA"))
resist<- fread("/data/htp/A07/CosmicData/CosmicResistanceMutations-2.tsv", select = c(1:7), na.strings=c("","NA"))
cancerGenes<- fread("/data/htp/A07/CosmicData/cancer_gene_census-2.csv", na.strings=c("","NA"))


##### DE in prot and rna
intersect(hal$GENE_NAME, ov$Gene.Name) #0
intersect(resist$`Gene Name`, ov$Gene.Name) #0
intersect(cancerGenes$`Gene Symbol`, ov$Gene.Name) #3 # "AKT3"     "HSP90AA1" "MYO5A"    "SPECC1"

canGenSub<-cancerGenes %>% filter(`Gene Symbol` %in% intersect(cancerGenes$`Gene Symbol`, ov$Gene.Name))

# AKT3: apoptosis 
# HSP90AA1: 
# MYO5A: apoptosis
# SPECC1: apoptosis


##### DE in prot
intersect(hal$GENE_NAME, prot.sig$T..PG.Genes) #14
intersect(resist$`Gene Name`, prot.sig$T..PG.Genes) #0
intersect(cancerGenes$`Gene Symbol`, prot.sig$T..PG.Genes) #31

hal_prot<-hal %>% filter(GENE_NAME %in% intersect(hal$GENE_NAME, prot.sig$T..PG.Genes))
canGenProt<-cancerGenes %>% filter(`Gene Symbol` %in% intersect(cancerGenes$`Gene Symbol`, prot.sig$T..PG.Genes))


##### DE in rna
intersect(hal$GENE_NAME, rna.sig$gene_name) #3
intersect(resist$`Gene Name`, rna.sig$gene_name) #1

# SMO, resistance, not detected in proteom, 
# AFF3: upregulation mediates tamoxifen resistance in breast cancers
# EIF3E: down, found in many cancer
# NOTCH2: up, but silenced in AML











change<-c("VII.vs.III","VII.vs.III","VII.vs.I","VII.vs.I","III.vs.I","III.vs.I")
direction<-c("up","down","up","down","up","down")
value<-as.numeric(c(1,0,35,22,2,1))

library(wesanderson)
bla<-as.data.frame(cbind(change, direction, value))
bla$value<-as.numeric(as.character(bla$value))
bla1=ggplot(bla, aes(x=change,y=value, fill=change, alpha=direction)) + geom_bar(position = "dodge2",stat="identity") +
  coord_flip()+theme_light()+theme(legend.position = "right")+
  scale_fill_manual(values = wes_palette("BottleRocket2",n=3))+
  theme(axis.title.y = element_blank())+
  scale_alpha_discrete(range = c(1, 0.5))+
  labs(y="Number of Proteins")

bla1
ggsave("diffProteinExpression.pdf", plot = bla1, width= 150, height= 80, units= "mm", device = "pdf",path = "/data/htp/long_term_data/Prot/")





















