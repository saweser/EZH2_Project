#### Proteom data Limma

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/Proteom")

library(UpSetR)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(limma)
library(tibble)

prot<-read.table("/data/htp/A07/EZH2_Project/Proteom/output_tables/Log2 values only.txt",
              header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=F, dec=",")

#### mean of the 3 technical replicates
prot<-prot %>%
  mutate(KO1=(KO1_1 + KO1_2 +KO1_3)/3) %>%
  mutate(KO2=(KO2_1 + KO2_2 +KO2_3)/3) %>%
  mutate(KO3=(KO3_1 + KO3_2 +KO3_3)/3) %>%
  mutate(KO4=(KO4_1 + KO4_2 +KO4_3)/3) %>%
  mutate(KO5=(KO5_1 + KO5_2 +KO5_3)/3) %>%
  mutate(KO6=(KO6_1 + KO6_2 +KO6_3)/3) %>%
  dplyr::select(KO1, KO2, KO3, KO4, KO5, KO6, WT_1, WT_2, WT_3, WT_4, WT_5, T..PG.Genes, T..PG.ProteinAccessions)

colnames(prot)<-c("KO_1", "KO_3", "KO_7", "KO_10", "KO_18", "KO_23", "WT_1", "WT_2", "WT_3", "WT_5", "WT_13", "T..PG.Genes", "T..PG.ProteinAccessions")

# sample names had to be changed back to the real ones
# Enes changed them
# WT1 à WT1
# WT2 à WT2
# WT3 à WT3
# WT4 à WT5
# WT5 à WT13
# KO1 à KO1
# KO2 à KO3
# KO3 à KO7
# KO4 à KO10
# KO5 à KO18
# KO6 à KO23

##### annotation 
anno<-colnames(prot)[1:11] %>% as.data.frame() %>%
  mutate(condition= c(rep("KO",6), rep("WT",5)))

colnames(anno)<-c("sample", "condition")
rownames(anno)<-anno$sample

##### log2 of expression
# values are already log transformed

##### add row names 
# some protein names are duplicated (some reviewed some unreviewd in UniProt)
dup<-prot[duplicated(prot$T..PG.Genes),]
# for now use the aceesion number as row names
rownames(prot)<-prot$T..PG.ProteinAccessions

library(RColorBrewer)
library(pheatmap)


##### Descriptive Statistics ##################################################################
plot.dat <- melt(as.matrix(prot[,c(1:11)]), value.name = "Expression", varnames = c("bla","Sample"), factorsAsStrings=FALSE)
plot.dat$Sample<-as.character(plot.dat$Sample)
#violin plots
p1 <- ggplot(plot.dat, aes(x=Sample, y=Expression)) +
  geom_violin(position=position_dodge(width=5))+
  geom_boxplot(width=0.2) + 
  coord_flip() + theme_bw() + 
  theme(axis.text= element_text(size=4), axis.title = element_text(size=6))+
  labs(x="log2 (Expression)")
p1
# denstiy plot
p2 <- ggplot(plot.dat, aes(x=Expression, col=Sample)) +
  geom_density() + ylab("Density") + 
  theme_bw() + theme(legend.position="right") + 
  theme(axis.text= element_text(size=5), axis.title = element_text(size=7))+
  labs(x="log2 (Expression)")
p2

#ggsave("Density.pdf",plot = cowplot::plot_grid(p1, p2, nrow=2), width= 200, height= 200, units= "mm", device = "pdf",path = "/home/weser/AML_LT/Proteome/plots")

# spearman rank correlation
cormat <- round(cor(prot[,c(1:11)], method = "spearman", use = "complete.obs"),2)
head(cormat)

ord=hclust(1-as.dist(cormat))$order
co=melt(cormat[ord,ord])

co$Var1<-factor(co$Var1, levels = rownames(cormat)[ord])
co$Var2<-factor(co$Var2, levels = rownames(cormat)[ord])

# plot
p<-ggplot(data = co, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile( colour = "white")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"))+
  scale_fill_gradient(low = "#132B43", high = "#56B1F7",
                      space = "Lab", na.value = "grey50", guide = "colourbar",
                      aesthetics = "colour")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())

p
ggsave("Spearman.pdf",plot = p, width= 150, height= 100, units= "mm", device = "pdf",path = "/data/htp/A07/EZH2_Project/Proteom/")

### PCA
library(wesanderson)
stage_col<-wes_palette("Cavalcanti1",n=5)[5:1]

source("/data/htp/A07/RNA_Seq/common_files/functions.R")

# rownames anno same oder than colnames counts
colnames(prot[,c(1:11)])==rownames(anno)
#anno<-anno[order(anno$name),]

PCA_cond <- pcaFunction(mat = prot[,c(1:11)], inf = anno, ngenes = 500, col = "condition")+ theme_classic()+
  scale_color_manual(values=c("#972D15", "#A2A475", "#D8B70A"))
PCA_cond

ggsave("PCA_cond.pdf", plot = PCA_cond, width= 150, height= 100, units= "mm", device = "pdf",path = "/data/htp/A07/EZH2_Project/Proteom/")

#################################################################################################

###### expression of EZH2 #############################################################################################
a<-prot %>% filter(T..PG.Genes %in% c("EZH2")) %>% 
  select(1:11) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sample=rownames(.))%>%
  mutate(gene=rep("EZH2",11))
df_a<-a %>%
  left_join(anno, by="sample")
colnames(df_a)[1]<-"EZH2"

bla1=ggplot(df_a, aes(x=condition,y=EZH2, fill=condition)) + 
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  theme(axis.title.x = element_blank())+
  ylab("log2 intensity")+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values=c("grey50", "grey28"))+
  theme(legend.position = "NONE")
bla1
ggsave("EZH2_expression_Prot.pdf", plot = bla1, width= 100, height= 90, units= "mm", device = "pdf"
       ,path = "/data/htp/A07/EZH2_Project/Proteom/")
########################################################################################################################

#### Differential Protein Expression ############################################################################################################
# log transformed values
# otherwise no normalization or filtering
# limma F test

# anno has same order ?
rownames(anno)==colnames(prot[,c(1:11)])
#anno<-anno[colnames(prot[,c(1:11)]),]

# leave out proteins with NA in all samples
keep <- apply(prot[,c(1:11)], 1, function(x) !all(is.na(x)))
prot <- prot[keep,]
# no preteins have been excluded!

#### model matrix and contrast 
design<-stats::model.matrix(~0+condition,data=anno)

cont.matrix <- makeContrasts(
  KO.vs.WT = conditionKO - conditionWT,
  levels = design)

# linear model fit
fit <- lmFit(object = prot[,c(1:11)],design = design, ndups = 1, spacing = 1, weights = NULL, method = "ls")
fit2 <- contrasts.fit(fit, contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.05, trend = T, robust = T)

# use limma F test so topTableF not decideTests
#res<-topTableF(fit2, number = Inf, adjust.method = "BH", sort.by = "none") %>%
#  rownames_to_column("Acession") %>%
#  left_join(prot[,c(12,13)], by=c("Acession"="T..PG.ProteinAccessions"))
#res.sig<-res %>% filter(adj.P.Val<0.05)
# same result topTable and topTableF only 2 conditions!

# TopTable
KO.vs.WT<-topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "none") %>%
  rownames_to_column("Acession") %>%
  left_join(prot[,c(12,13)], by=c("Acession"="T..PG.ProteinAccessions"))

# add direction(up and down)
KO.vs.WT<-KO.vs.WT %>%
  dplyr::mutate(dir=ifelse(adj.P.Val<0.05 & logFC>0 ,"up",ifelse(adj.P.Val<0.05 & logFC<0, "down","same"))) 

# significant genes only
KO.vs.WT.sig<-KO.vs.WT %>% filter(adj.P.Val<0.05)

# save as rds
saveRDS(KO.vs.WT, "/data/htp/A07/EZH2_Project/Proteom/diffProtTopTable.rds")
###################################################################################

##### volcano plot ##################################################################################
KO.vs.WT<-readRDS("/data/htp/A07/EZH2_Project/Proteom/diffProtTopTable.rds")

res<-KO.vs.WT%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p <- ggplot(data = res, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " adjusted p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 16, colour = "black",family = "Arial"),
        axis.text = element_text(size = 12, colour = "black",family = "Arial"), strip.text = element_text(size = 10, family = "Arial")) 
volcano.p 
ggsave("volcanoDE.png", plot = volcano.p, width= 150, height= 120, units= "mm", device = "png"
       ,path = "/data/htp/A07/EZH2_Project/Proteom/")
#######################################################################################################

##### reactome pathway enrichment #####################################################################
library(ReactomePA)

keytypes(org.Hs.eg.db)

entrez<-AnnotationDbi::select(org.Hs.eg.db, keys= KO.vs.WT$T..PG.Genes,
                                  keytype="SYMBOL", columns=c("ENSEMBL","ENTREZID","GENENAME") )

df<-KO.vs.WT %>% left_join(entrez, by=c("T..PG.Genes"="SYMBOL"))
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
gse<- enrichGO(unique(df.sig$ENSEMBL),OrgDb = "org.Hs.eg.db", keyType="ENSEMBL", ont = "BP",pAdjustMethod = "BH", 
               pvalueCutoff=0.05,qvalueCutoff= 0.05,readable = TRUE)
# BF: Biological Process
# MF: Molecular Function
# CC: 
gse.data<-as.data.frame(gse)

# simplify result, combine closely related terms
gse2 <- simplify(gse, cutoff=0.7, by="p.adjust", select_fun=min)
gse2_result<-gse2@result

# dotplot
dotplot(gse2, showCategory=20)

# emapplot
emapplot(gse)
dev.off()
# cnetplot(gse)
#goplot(gse)
#heatplot(gse)

##### KEGG pathways
# KEGG only takes EntezIDs
ekegg <- enrichKEGG(gene = df.sig$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
# dotplot
dotplot(ekegg, showCategory=20)

# interesting Go here is Ras signaling pathway
# plot this one as cnetplot to see the genes involved

# to show gene symbols instead entrezIDs

y <- setReadable(ekegg, 'org.Hs.eg.db', keyType="ENTREZID")

# subset the enrichResult for interesting pathway/s
y@result = y@result[y@result$Description %in% c("Ras signaling pathway"), ]

cnetplot(y)
dev.off()

















