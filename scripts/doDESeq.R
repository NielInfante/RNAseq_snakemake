# Error logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(tximport)
#library(rjson)
library(DESeq2)
library('AnnotationDbi')
library("RColorBrewer")
library(ggrepel)
library(gplots)


#DelayedArray


organism <- snakemake@config$organism
library(organism)

# Assume this script will be invoked with the arguments (in order):
# prefix  
# metadata_file


# I'm not currently using this, but I don't want to loose the technique
#iris %>% filter(Sepal.Length > 5.3 & Species == 'setosa')

#x <- "Sepal.Length > 5.3 & Species == 'setosa'"
#fc <- rlang::parse_expr(x)
#iris %>% filter(!!fc)


#args <- commandArgs()

#meta <- read_tsv(args[2])

meta_file_name <- snakemake@config$metadata_file
meta <- read_tsv(meta_file_name)

exp <- snakemake@params$exp

snakemake@source(paste0('deseq/', exp, '/config.R'))


outDir <- "deseq"

files <- paste0('salmon/', samples, '/quant.sf')

txi <- tximport(files, type='salmon', tx2gene = tx2gene)

dds <- DESeqDataSetFromTximport(txi, meta, design)

# How many genes, out of those with at least a single count, have three samples with a count of 10 or more
dds <- dds[rowSums(counts(dds)) > 0,]
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,] # filter them out

dds <- DESeq(dds)

res<-results(dds, contrast=contrast)
res<-res[order(res$padj),]
res <- as.data.frame(res)

# Save dds
saveRDS(dds, paste0(outDir, '/', outPrefix, '/dds.rds'))


my_concat <- function(x){paste(x, sep="|", collapse="|")}

# Get gene names
res$Gene <- mapIds(orgDB, keys=row.names(res), column='SYMBOL', keytype='ENSEMBL', multiVals=my_concat)
res$ID <- row.names(res)

# Use Gene ID in place where there is no gene name
idx <- is.na(res$Gene)
res$Gene[idx] <- res$ID[idx]
idx <- which(res$Gene == 'NA')  # mapIDs with my_concat can return "NA", not NA
res$Gene[idx] <- res$ID[idx]


# Write Results
outResults <- data.frame(GeneID=res$ID, Gene=res$Gene, baseMean=res$baseMean, stat=res$stat, log2FoldChange=res$log2FoldChange, pvalue=res$pvalue, padj=res$padj)
name <- paste(outDir, '/', outPrefix, '/results.txt', sep="") 
write.table(outResults, file=name, sep="\t", quote=F, row.names=F)

# Significant genes
r2 <- res[!(is.na(res$padj)),]
resSig <- r2[ r2$padj < 0.05, ]
resTable <- data.frame(GeneID=row.names(resSig), Gene=resSig$Gene, baseMean=resSig$baseMean, stat=resSig$stat, log2FoldChange=resSig$log2FoldChange, pvalue=resSig$pvalue, padj=resSig$padj)
write.table(resTable,file=paste(outDir, "/", outPrefix, "_significant.txt", sep=""), sep="\t", quote=F, row.names=F)



##########  Sanity Check
# Plot counts of most significant, to check if fold change is right
png(paste0(outDir, '/',outPrefix,'/sanity.check.png'))
plotCounts(dds, gene=res[1,]$ID, intgroup = PCA_Group, main=res[1,]$Gene, pch=19)
dev.off()



