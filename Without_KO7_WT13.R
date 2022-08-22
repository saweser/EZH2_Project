########## Differential gene expression: Limma ########################################

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/RNA_Seq")

library(limma)
library(edgeR)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(ggrepel)
library(extrafont)
library(tibble)
loadfonts()
source("/data/htp/A07/RNA_Seq/common_files/functions.R")

anno<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/anno.rds")
counts<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/counts.rds")
#######################################################################################

##### leave out WT13 and KO7
anno <-anno %>% filter(!name %in% c("WT_13", "KO_7"))
counts<-counts %>% dplyr::select(-WT_13, -KO_7)

colnames(counts)==rownames(anno)
# anno<-anno[colnames(counts),]

### for GEO upload ###################################################################
# annoGEO<-anno[order(anno$condition, anno$number),]
# #for GEO upload: annotation
# annoGEO<-annoGEO[,c(1:3)]
# colnames(annoGEO)<-c("sample_name", "clone", "condition")
# write.table(annoGEO, "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/Without_WT13_KO7/GEO_Upload/annoGEO.txt")
# #for GEO upload: count table
# countsGEO<-counts[,rownames(annoGEO)]
# write.table(countsGEO, "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/Without_WT13_KO7/GEO_Upload/EZH2_KO_WT_counts.txt")
# # index file for deML
# index<-anno[,c("BC","name")]
# write.table(index, "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/fastq_files/deML/index.txt", 
#             quote = FALSE, sep = " ", row.names = FALSE)
#####################################################################################

##### Diff Exp
# make DGE object
dge <- DGEList(counts = counts, lib.size = colSums(counts), samples=anno, remove.zeros = TRUE)
dge <- edgeR::calcNormFactors(dge, method = "RLE")


# filter very little expressed genes
# find CPM value that corresponds to a read count of 10 
min.cpm <- cpm(10, mean(dge$samples$lib.size))
# retain genes with at least 1 sample above this threshold 
keep <- rowSums(cpm(dge) > as.vector(min.cpm)) >= 1
table(keep)
dge <- dge[keep, ]

# recalculate norm factors
dge <- edgeR::calcNormFactors(dge, method = "RLE")

norm.counts<-cpm(dge)
saveRDS(norm.counts, "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/norm.counts.rds")


#### model matrix and contrast:  ###############################################################
design<-stats::model.matrix(~0+condition,data=anno)
colnames(design)<-gsub("condition","",colnames(design))

cont.matrix <- makeContrasts(
  KO_WT = KO - WT,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
saveRDS(v, "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/voom.rds")

fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

res_KO_WT<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)%>%
  rownames_to_column("ensembl_gene_id")

#### add gene name 
library(biomartr)
gff<-plyranges::read_gff("/data/htp/A07/EZH2_Project/RNA_Seq/anno/Homo_sapiens.GRCh38.100.chr.chrNamesUCSC.gtf") %>% 
  as.data.frame  %>%
  filter(type=="gene") # subset the gff by gene

res_KO_WT<-res_KO_WT %>% left_join(gff, by=c("ensembl_gene_id"="gene_id"))

##### add direction (up and down)
res_KO_WT<-res_KO_WT %>%
  dplyr::mutate(dir=ifelse(adj.P.Val<0.05 & logFC>0 ,"up",ifelse(adj.P.Val<0.05 & logFC<0, "down","same"))) 

##### save DE table #################################################################################
saveRDS(res_KO_WT, "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/DE_Genes_without_WT13_KO7.rds")

##### volcano plot ##################################################################################
res<-res_KO_WT%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p <- ggplot(data = res, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " adjusted p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 16, colour = "black",family = "Arial"),
        axis.text = element_text(size = 12, colour = "black",family = "Arial"), strip.text = element_text(size = 10, family = "Arial")) 
volcano.p 
ggsave("volcanoDE.png", plot = volcano.p, width= 150, height= 120, units= "mm", device = "png"
       ,path = "/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2")
#######################################################################################################

##### genes with highest logFC ########################################################################
hiFC<-res_KO_WT_p %>%
  filter(logFC>4)
######################################################################################################

##### reactome pathway enrichment #####################################################################
res_KO_WT<-readRDS("/data/htp/A07/EZH2_Project/RNA_Seq/2020_09_EZH2/DE_Genes_without_WT13_KO7.rds")

library(ReactomePA)
keytypes(org.Hs.eg.db)

entrez<-AnnotationDbi::select(org.Hs.eg.db, keys= res_KO_WT$ensembl_gene_id,
                              keytype="ENSEMBL", columns=c("ENSEMBL","ENTREZID","GENENAME") )

df<-res_KO_WT %>% left_join(entrez, by=c("ensembl_gene_id"= "ENSEMBL"))
df.sig<-df %>% filter(adj.P.Val<0.05)

x <- enrichPathway(gene= unique(df.sig$ENTREZID), readable=T, pvalueCutoff = 0.05)
# set readable=TRUE to see gene symbols in the plots
dotplot(x, showCategory=20)

xx <- enrichPathway(gene= unique(df.sig$ENTREZID), universe=unique(entrez$ENTREZID), readable=T, pvalueCutoff = 0.05)
dotplot(xx, showCategory=20)

# dotplot
dotplot(x, showCategory=20)

#### Over representation analysis (ORA) ###########################################################################
library(clusterProfiler)
# plots from enrichplot
# enrichGO
gse<- enrichGO(unique(df.sig$ENTREZID),OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP",pAdjustMethod = "BH", 
               pvalueCutoff=0.05,qvalueCutoff= 0.05,readable = TRUE)
# BF: Biological Process
# MF: Molecular Function
# CC: 
gse.data<-as.data.frame(gse)
dotplot(gse, showCategory=20)


# simplify result, combine closely related terms
gse2 <- simplify(gse, cutoff=0.7, by="p.adjust", select_fun=min)
gse2_result<-gse2@result

# dotplot
dotplot(gse2, showCategory=20)







