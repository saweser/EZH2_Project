##### Peak Annotation ####################################################################
rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/ChipSeq")

library(plyranges)
library(tidyverse)

gapped_1<-read.table("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2/ENCFF678EGD_peaks.gappedPeak")
broad_1<-read.table("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2/ENCFF678EGD_peaks.broadPeak")
gapped_2<-read.table("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2/ENCFF912PNB_peaks.gappedPeak")
broad_2<-read.table("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/macs2/ENCFF912PNB_peaks.broadPeak")

colnames<-c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")


colnames(gapped_1)<-colnames
colnames(broad_1)<-colnames
colnames(gapped_2)<-colnames
colnames(broad_2)<-colnames

broad_1$strand<-rep("*",nrow(broad_1))
broad_2$strand<-rep("*",nrow(broad_2))

broad_1<-as_granges(broad_1)
broad_2<-as_granges(broad_2)

#### replikate1 only has 1000 peaks while replicate2 has 30,000
a<-find_overlaps(broad_1, broad_2) %>% as.data.frame()
a<-find_overlaps(broad_2, broad_1)
# almost all peaks of broad1 overlap wih peaks in broad2
# probably broad1 did not work so well
# much more likely that the 30,000 number is correct, as we have TF as well has histone modification function.
# proceed with broad2 only!

##### annotate peaks with nearest gene
txdb<-GenomicFeatures::makeTxDbFromGFF(file = "/data/htp/A07/EZH2_Project/RNA_Seq/anno/Homo_sapiens.GRCh38.100.chr.chrNamesUCSC.gtf", format = "gtf")

library(ChIPseeker)
peakAnno <- annotatePeak(broad_2,TxDb=txdb, sameStrand = TRUE, level = "transcript")
# level=transcript is the default
peakAnno.data<-as.data.frame(peakAnno)

# tssRegion=c(-1500, 1000)


#### add gene name
gff<-read_gff("/data/htp/A07/EZH2_Project/RNA_Seq/anno/Homo_sapiens.GRCh38.100.chr.chrNamesUCSC.gtf") %>%
  as.data.frame %>%
  filter(type=="gene") # subset the gff by gene

peakAnno.data<-peakAnno.data %>% left_join(gff, by=c("geneId"="gene_id"))
saveRDS(peakAnno.data, "/data/htp/A07/EZH2_Project/ChipSeq/peakAnno.rds")


#### where are peaks located?
p<-plotAnnoBar(peakAnno) +
  theme_classic() +
  theme(legend.position = "bottom")
p
ggsave("FeatureDistribution.pdf", plot = p, width= 170, height= 60, units= "mm", device = "pdf"
       ,path = "/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/")


#### peaks in promoter regions
# -1500, 1000
prom<-peakAnno.data %>%
  filter(annotation=="Promoter")

#### peaks in promoter regions protein coding
promCod<-prom %>%
  filter(gene_biotype=="protein_coding") %>%
  select(geneId, gene_name) %>%
  filter(!duplicated(.))
summary(duplicated(promCod))

#### Over representation analysis (ORA) ###########################################################################
library(clusterProfiler)
# enrichGO
gse<- enrichGO(promCod$geneId,OrgDb = "org.Hs.eg.db", keyType="ENSEMBL", ont = "BP",pAdjustMethod = "BH",
               pvalueCutoff=0.05,qvalueCutoff= 0.05,readable = TRUE)
# BF: Biological Process
# MF: Molecular Function
# CC:
gse.data<-as.data.frame(gse)

# dotplot
pdf("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/enrichGo_dotplot.pdf", width=8, height=5)
dotplot(gse, showCategory=20)
dev.off()

# emapplot
pdf("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/enrichGo_emapplot.pdf", width=8, height=5)
emapplot(gse)
dev.off()
# cnetplot(gse)
#goplot(gse)
#heatplot(gse)

##### KEGG pathways
# KEGG only takes EntezIDs
entrez<-AnnotationDbi::select(org.Hs.eg.db, keys =promCod$geneId,
                              keytype = "ENSEMBL",
                              columns = c("ENTREZID"))

ekegg <- enrichKEGG(gene = entrez$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
# dotplot
pdf("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/kegg_dotplot.pdf", width=8, height=5)
dotplot(ekegg, showCategory=20)
dev.off()

# interesting Go here is Ras signaling pathway
# plot this one as cnetplot to see the genes involved

# to show gene symbols instead entrezIDs
y <- setReadable(ekegg, 'org.Hs.eg.db', keyType="ENTREZID")

# subset the enrichResult for interesting pathway/s
y@result = y@result[y@result$Description %in% c("Ras signaling pathway"), ]

pdf("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/kegg_cnetplot_ras.pdf", width=8, height=5)
cnetplot(y)
dev.off()

#### reactome pathway enrichment
library(ReactomePA)
x <- enrichPathway(gene= entrez$ENTREZID, readable=T, pvalueCutoff = 0.05)

# dotplot
pdf("/data/htp/A07/EZH2_Project/ChipSeq/EncodeData/reactome_dotplot.pdf", width=8, height=5)
dotplot(x, showCategory=20)
dev.off()


# do the IDR
