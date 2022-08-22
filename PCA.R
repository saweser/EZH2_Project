########## Normalization and PCA #######################################################################
# normalize count data with edgeR
# then calculcate UMIs per million (UPM)
# PCA of normalized upm

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/RNA_Seq")

source("/data/htp/A07/RNA_Seq/common_files/functions.R")
library(edgeR)
library(ggplot2)
library(tidyr)
library(cowplot)
library(gridExtra)
library(ggrepel)
library(ggpubr)
library(ggbeeswarm)


anno<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/anno.rds")
counts<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/counts.rds")

######################################################################################################################

##### leave out WT13 and KO7
anno <-anno %>% dplyr::filter(!name %in% c("WT_13", "KO_7"))
counts<-counts %>% dplyr::select(-WT_13, -KO_7)

colnames(counts)==rownames(anno)
#anno<-anno[colnames(counts),]

### Normalization: EdgeR ##########################################################################################
# calculate UPM of normalized counts
rownames(anno)==colnames(counts)

sums <- colSums(counts)

# normalization: RLE
nf <- edgeR::calcNormFactors(counts, method = "RLE")
upm<-as.data.frame(t((t(counts*nf)/sums)* 1e+06)) # transform, because R will calculate row wise

#### PCA ##############################################################################################################
PCA_cond <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "condition")+
  theme_classic()+
  geom_text_repel(label = anno$name, force = 1, point.padding = 0.1, show.legend = FALSE)+
  scale_colour_manual(values = alpha(c("KO" = "blue", "WT" = "black"),1))+
  theme(legend.position = "none")+
  theme(axis.line = element_line(colour = 'black', size = 0.8))+
  theme(axis.title=element_text(size=14,face="bold"))
  
PCA_cond
ggsave("PCA_condition.pdf", plot = PCA_cond, width= 120, height= 150, units= "mm", device = "pdf",path = "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/Without_WT13_KO7/")



#resistance phenotype
# KO samples with high IC50=YES, same IC50 than WT =NO, non KO samples =NA (fig4c)
#PCA_res <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "resistance")+
#  theme_classic()
#PCA_res
#ggsave("PCA_resistance_PC1_2_WTKO.pdf", plot = PCA_res, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/")
#######################################################################################################################

###### expression of EZH2 #############################################################################################
a<-upm %>% filter(rownames(upm) %in% c("ENSG00000106462")) %>% 
  t() %>%
  as.data.frame() %>%
  mutate(sample=rownames(.))%>%
  mutate(gene=rep("EZH2",9))
df_a<-a %>%
  left_join(anno, by=c("sample"="name"))

bla1=ggplot(df_a, aes(x=condition,y=ENSG00000106462, fill=condition)) + 
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  theme(axis.title.x = element_blank())+
  ylab("normalized counts (upm)")+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values=c("grey50", "grey28"))+
  theme(legend.position = "NONE")
bla1
ggsave("EZH2_expression_KOWT.pdf", plot = bla1, width= 140, height= 80, units= "mm", device = "pdf"
       ,path = "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/Without_WT13_KO7/")
########################################################################################################################

##### Target Genes of EZH2 #############################################################################################
#### HOXB7
# should be upregulated in KO
# color the high and low resitance phenotype KOs
a<-upm %>% filter(rownames(upm) %in% c("ENSG00000260027")) %>% 
  t() %>%
  as.data.frame() %>%
  mutate(sample=rownames(.))%>%
  mutate(gene=rep("HOXB7",9))
df_a<-a %>%
  left_join(anno, by=c("sample"="name"))

bla1=ggplot(df_a, aes(x=condition,y=ENSG00000260027, fill=condition)) + 
  geom_boxplot()+ 
  #geom_beeswarm(priority='density',cex=2.5, aes(color=resistance)) +
  theme_classic()+
  theme(axis.title.x = element_blank())+
  ylab("normalized counts (upm)")+
  ggtitle("HOXB7")
bla1

ggsave("HOXB7_expression_KOWT.pdf", plot = bla1, width= 140, height= 80, units= "mm", device = "pdf"
       ,path = "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/Without_WT13_KO7/")


#bla2=ggplot(df_a, aes(x=resistance,y=V1, fill=resistance)) + 
#  geom_boxplot()+ 
#  geom_point() +
#  theme_classic()+
#  theme(axis.title.x = element_blank())+
#  ylab("normalized counts (upm)")+
#  ggtitle("HOXB7")
#bla2

#ggsave("HOXB7_KOWT_resistance.pdf", plot = bla2, width= 140, height= 80, units= "mm", device = "pdf"
       ,path = "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2")


#### check if housekeeping genes are expressed
# ENSG00000075624 ACTB
# ENSG00000111640 GAPDH
a<-upm %>% filter(rownames(upm) %in% c("ENSG00000075624")) %>% 
  t() %>%
  as.data.frame() %>%
  mutate(sample=rownames(.))
df_a<-a %>%
  left_join(anno, by=c("sample"="name"))

bla1=ggplot(df_a, aes(x=condition,y=ENSG00000075624)) + 
  geom_boxplot()+ 
  theme_classic()+
  theme(axis.title.x = element_blank())+
  ylab("normalized counts (upm)")+
  ggtitle("GAPDH")
bla1

##### heatmap ###########################################################################################
library(gplots)
library(RColorBrewer)
library(matrixStats)
upm<-as.matrix(upm)
rowVar <-rowVars(upm)
mv100 <- order(rowVar, decreasing = T )[1:20] # take out the  most variable genes

heatmap_genes <- as.matrix(upm[mv100,]) #select expression values of the top500 DE genes

library(viridis)
library(scales)
hmcol<-viridis_pal(option = "magma",direction = -1)(9) # color for heatmap
# show_col(viridis_pal(option = "magma")(9))

pdf("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/Without_WT13_KO7/heatmap.pdf")
heatmap.2(heatmap_genes, trace="none", margins = c(9,9), 
          dendrogram = "column", labRow = rownames(heatmap_genes), keysize = 1.2, symkey= FALSE, density.info="none", 
          key.xlab= "Expression (UPM)",cexCol=1, key.par=list(), key.xtickfun = NULL, 
          key.ytickfun = NULL,col = hmcol, colsep=1:ncol(heatmap_genes), rowsep=1:nrow(heatmap_genes),sepwidth=c(0.01,0.01))
dev.off()
# cluster in 2 separate groups! but does not fit the annotation



####
#####


# where are those genes?
library(biomartr)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

attributes = listAttributes(ensembl)


getBM(attributes=c("chromosome_name", "ensembl_gene_id","external_gene_name"), 
      filters = "ensembl_gene_id", 
      values = rownames(heatmap_genes), 
      mart = ensembl)

# many are still mitochondria genes


####### removed MT genes ##################################################################################
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)

chromName<-getBM(attributes=c("chromosome_name", "ensembl_gene_id","external_gene_name"), 
                 filters = "ensembl_gene_id", 
                 values = rownames(upm), 
                 mart = ensembl)

a<-chromName %>% dplyr::group_by(chromosome_name) %>% dplyr::count() # 36 MT genes

chromName<-chromName %>%
  filter(!chromosome_name=="MT")
counts<-counts[chromName$ensembl_gene_id,]











