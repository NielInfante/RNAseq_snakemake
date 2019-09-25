# Error logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(tximport)
#library(rjson)
library(DESeq2)
library("RColorBrewer")
library(ggrepel)
library(gplots)



#outPrefix <- '120d'
#PCA_Group <- 'Time'
#design =~ Time
#contrast <- c('Time', '120d', '0')

#meta <- meta %>% filter(Time == '120d' | Time == '0')
#samples <- meta$Sample
#meta$Graph_Display <- meta$Sample




print("Loaded packages")
meta_file_name <- snakemake@config$metadata_file
meta <- read_tsv(meta_file_name)

tx2gene <- read_tsv(snakemake@input[['id']])


# This gets the prefix we need
exp <- snakemake@params$exp

# This sources the config file the previous rule created
snakemake@source(paste0('../deseq/', exp, '/config.R'))


# script runs in .snakemake/scripts
outDir <- "deseq"

files <- paste0('salmon/', samples, '/quant.sf')

getwd()
print("Files are:")
print(files)

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
res$GeneID <- row.names(res)

# Save dds
saveRDS(dds, paste0(outDir, '/', outPrefix, '/dds.rds'))


biotype <- read_tsv('Data/Biotype')
biotype <- left_join(biotype, tx2gene)

# Only need one transcript per gene
biotype <- biotype[!duplicated(biotype$GeneID),]

res <- left_join(res, biotype)

# Write Results
outResults <- data.frame(GeneID=res$GeneID, Gene=res$GeneName, baseMean=res$baseMean, stat=res$stat, log2FoldChange=res$log2FoldChange, pvalue=res$pvalue, padj=res$padj)
name <- paste(outDir, '/', outPrefix, '/results.txt', sep="") 
write.table(outResults, file=name, sep="\t", quote=F, row.names=F)

# Significant genes
r2 <- res[!(is.na(res$padj)),]
resSig <- r2[ r2$padj < 0.05, ]
resTable <- data.frame(GeneID=resSig$GeneID, Gene=resSig$GeneName, baseMean=resSig$baseMean, stat=resSig$stat, log2FoldChange=resSig$log2FoldChange, pvalue=resSig$pvalue, padj=resSig$padj)
write.table(resTable,file=paste(outDir, "/", outPrefix, "/significant.txt", sep=""), sep="\t", quote=F, row.names=F)



##########  Sanity Check
# Plot counts of most significant, to check if fold change is right
png(paste0(outDir, '/',outPrefix,'/sanity.check.png'))
title=paste(res[1,]$GeneName, "\nFold Change:",res[1,]$log2FoldChange)
plotCounts(dds, gene=res[1,]$GeneID, intgroup = PCA_Group, main=res[1,]$GeneName, 
					 sub=paste('FC:', format(res[1,]$log2FoldChange, digits=2)), pch=19)
dev.off()



