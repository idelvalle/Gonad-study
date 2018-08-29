

rm(list=ls())
setwd("~/Desktop/RNA-Seq-Brain-June-2018/Gonads/")

library("DESeq2")
library(RColorBrewer)
library(gplots)
library(ggplot2)
library("ggbeeswarm")
library(reshape2)
library(plyr)
library("dplyr")
library("magrittr")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("clusterProfiler")
library(org.Hs.eg.db)
library("AnnotationDbi")
library("Homo.sapiens")
library("limma")
library("maSigPro")
library("MASS")
library(edgeR)
library(NOISeq)



countdata <- read.table(file="counts.all.csv", header=TRUE, sep=",", row.names=1)

head(countdata)

colnames(countdata)

phenodata <- read.table("phenodata.csv", sep=",", header=TRUE, row.names=1)

phenodata

ncol(countdata)
nrow(phenodata)

colnames(countdata) <- phenodata$Sequencing_ID

phenodata.gonads <- phenodata %>% filter(Area=="Gonad")
phenodata.gonads
countdata.gonads <- countdata[,phenodata.gonads$Sequencing_ID]

phenodata.ovary <- phenodata %>% filter(Tissue=="Ovary")
phenodata.ovary
countdata.ovary <- countdata[,phenodata.ovary$Sequencing_ID]


#DESeq2

countdata <- as.matrix(countdata)
head(countdata)

tissue <- factor(phenodata.gonads$Tissue)
tissue


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata.gonads), tissue)
coldata
dds <- DESeqDataSetFromMatrix(countData=countdata.gonads, colData=coldata, design=~tissue)

dds

dds <- dds[ rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

#BUILDING RESULTS TABLE
resultsNames(dds)

res <- results(dds, name = "tissue_Testis_vs_Ovary", cooksCutoff = FALSE)


mcols(res, use.names=TRUE)
summary(res)

columns(Homo.sapiens)

res.symbol <- res

res.symbol$symbol <- mapIds(Homo.sapiens,keys=row.names(res),column="SYMBOL",keytype="ENSEMBL",multiVals="first")

res.symbol <- res.symbol[complete.cases(res.symbol[, 6:7]),] #Remove NA values from Gene Symbol & padj
nrow(res.symbol)

res.dup <- res.symbol[duplicated(res.symbol$symbol),] #Visualize duplicates

nrow(res.dup)

write.csv(res.symbol, file="results.testis.vs.ovary.csv")
write.csv(res.dup, file="results.testis.vs.ovary.dup.csv")

#DESeq2

stage <- factor(phenodata.ovary$Stage)
stage


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata.ovary <- data.frame(row.names=colnames(countdata.ovary), stage)
coldata.ovary
dds.ovary <- DESeqDataSetFromMatrix(countData=countdata.ovary, colData=coldata.ovary, design=~stage)

dds.ovary

dds.ovary <- dds.ovary[ rowSums(counts(dds.ovary)) > 10, ]

dds.ovary <- DESeq(dds.ovary)

#BUILDING RESULTS TABLE
resultsNames(dds.ovary)

res.ovary <- results(dds.ovary, name = "stage_F3_vs_F2", cooksCutoff = FALSE)


mcols(res.ovary, use.names=TRUE)
summary(res.ovary)

columns(Homo.sapiens)

res.ovary.symbol <- res.ovary

res.ovary.symbol$symbol <- mapIds(Homo.sapiens,keys=row.names(res.ovary),column="SYMBOL",keytype="ENSEMBL",multiVals="first")

res.ovary.symbol <- res.ovary.symbol[complete.cases(res.ovary.symbol[, 6:7]),] #Remove NA values from Gene Symbol & padj
nrow(res.ovary.symbol)

res.ovary.dup <- res.ovary.symbol[duplicated(res.ovary.symbol$symbol),] #Visualize duplicates

nrow(res.ovary.dup)

write.csv(res.ovary.symbol, file="results.F3.vs.F2.ovary.csv")
write.csv(res.ovary.dup, file="results.F3.vs.F2.ovary.dup.csv")


### MA PLOT #


DESeq2::plotMA(res, alpha = 0.05, ylim=c(-10,10)) #Points colored in red if padj is lower than 0.05.Points which fall out of the window are plotted as open triangles pointing either up or down.

DESeq2::plotMA(res.ovary, alpha = 0.05, ylim=c(-10,10))

# Regularized log transformation for clustering/heatmaps, etc. rlog tends to work well on small datasets (n <30)
rld <- rlogTransformation(dds, blind = TRUE)
head(assay(rld),3)
hist(assay(rld))
names(colData(dds))

plotPCA(rld, intgroup="tissue")


pcaData <- plotPCA(rld, intgroup="tissue", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=tissue, shape=tissue, label = phenodata.gonads$Stage)) +
  geom_point(size=3) + geom_text(aes(label=phenodata.gonads$Stage),hjust=1.5, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

rld.ovary <- rlogTransformation(dds.ovary, blind = TRUE)

names(colData(dds.ovary))

plotPCA(rld.ovary, intgroup="stage")


pcaData.ovary <- plotPCA(rld.ovary, intgroup="stage", returnData=TRUE)
percentVar.ovary <- round(100 * attr(pcaData.ovary, "percentVar"))
ggplot(pcaData.ovary, aes(PC1, PC2, color=stage, shape=stage, label = phenodata.ovary$EmbryoID)) +
  geom_point(size=3) + geom_text(aes(label=phenodata.ovary$EmbryoID),hjust=0.5, vjust=1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#### SAMPLE DISTANCES #####

# Sample Distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Sample Distances
sampleDists.ovary <- dist(t(assay(rld.ovary)))
sampleDistMatrix.ovary <- as.matrix( sampleDists.ovary)
rownames(sampleDistMatrix.ovary) <- colnames(rld.ovary)
colnames(sampleDistMatrix.ovary) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix.ovary,
         clustering_distance_rows = sampleDists.ovary,
         clustering_distance_cols = sampleDists.ovary,
         col = colors)

#### TPM PLOTS #####

gene.length <- read.table("gene.length.csv", row.names=1, sep=",",header=T)
head(gene.length)

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

tpms <- apply(countdata.gonads, 2, function(x) tpm(x, gene.length$Length))
tpms.ovary <- apply(countdata.ovary, 2, function(x) tpm(x, gene.length$Length))

colSums(tpms)
colMeans(tpms)

colSums(tpms.ovary)
colMeans(tpms.ovary)

write.table(tpms, file="TPM.Gonads.csv", sep=",")
write.table(tpms.ovary, file="TPM.Ovary.csv", sep=",")

tpm <- read.table(file="TPM.Gonads.csv", sep=",", header=TRUE, row.names = 1)
head(tpm)

tpm.ovary <- read.table(file="TPM.Ovary.csv", sep=",", header=TRUE, row.names = 1)
head(tpm.ovary)

theme_set(theme_gray(base_size = 25))

phenodata.gonads
phenodata.ovary

gene <- as.data.frame(t(tpm["ENSG00000177508",]))
gene$Stage <- phenodata.gonads$Stage
gene$Tissue <- phenodata.gonads$Tissue
#gene$Sex <- phenodata.gonads$Sex
gene
write.table(gene, file="DHX37.tpm.csv", sep=",")
#melt_gene<-melt(gene)

gene.ovary <- as.data.frame(t(tpm.ovary["ENSG00000095596",]))
gene.ovary$Stage <- phenodata.ovary$Stage
#gene$Sex <- phenodata.gonads$Sex
gene.ovary
write.table(gene, file="DHX37.tpm.csv", sep=",")
#melt_gene<-melt(gene)

ggplot(gene.ovary, aes(x=Stage ,y=ENSG00000095596,fill=Stage))+geom_boxplot() + ylab("CYP26A1")


save.image(file="Gonads.RData")
