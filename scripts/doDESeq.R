# Error logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(tximport)
#library(rjson)
library(DESeq2)



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
#snakemake@source(paste0('../experiements/', exp, '/deseq/', 'config.R'))
source(paste0('results/', exp, '/deseq/config.R'))


# script runs in .snakemake/scripts
outDir <- paste0("results/", exp, "/deseq/")

files <- paste0('salmon/', samples, '/quant.sf')

print("Files are:")
print(files)

txi <- tximport(files, type='salmon', tx2gene = tx2gene)

tpm <- as.data.frame(txi$abundance)
names(tpm) <- paste0(samples, '_TPM')
tpm$meanTPM <- rowMeans(tpm)
tpm$GeneID <- row.names(tpm)

print('Did tximport')

# Import into DESeq object
dds <- DESeqDataSetFromTximport(txi, meta, design)

# How many genes, out of those with at least a single count, have three samples with a count of 10 or more
dds <- dds[rowSums(counts(dds)) > 0,]
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,] # filter them out

# Do the DESeq analysis
dds <- DESeq(dds)

# Save dds
saveRDS(dds, paste0(outDir, 'dds.rds'))


# Get results and extra data
res<-results(dds, contrast=contrast)
res<-res[order(res$padj),]
res <- as.data.frame(res)
res$GeneID <- row.names(res)

# Add count data
cnt <- as.data.frame(counts(dds, normalized=T))
names(cnt) <- colData(dds)[,sample_column]
cnt$GeneID <- row.names(cnt)
res <- inner_join(res, cnt)


# Add Biotypes
biotype <- read_tsv('Data/Biotype')
biotype <- left_join(biotype, tx2gene)

# Only need one transcript per gene
biotype <- biotype[!duplicated(biotype$GeneID),]
res <- left_join(res, biotype)
res$TransID <- NULL

# Add TPM
res <- left_join(res, tpm)

# Significant genes

## Get an appropriate base mean cutoff value
# Old way, now trying clustering

den <- density(log2(res$meanTPM + 1))
maxY <- max(den$y)
#maxAt <- which.max(den$y)

#for (i in (maxAt + 1):length(den$y)){
#	diff <- den$y[i] - den$y[i-1]
#	if (diff > 0){
#		cut_index <- i - 1
#		break
#	}
#}

#if (cut_index < 0){
#	cut_value <- median(log2(res$baseMean))
#} else {
#	cut_value <- den$x[cut_index]	
#}


# Use k means clustering on log2 mean TPM plus 1
# First cluster should be basal expression,
# Second cluster should be genes we're interested in.
km <- kmeans(log2(res$meanTPM + 1), 2)

res$Cluster <- km$cluster
# Want  to keep the cluster with the larger expression value
if (km$centers[1] > km$centers[2]){
	res <- res %>% mutate(Cluster=ifelse(Cluster==1, 'In','Out'))
} else {
	res <- res %>% mutate(Cluster=ifelse(Cluster==1, 'Out','In'))
}

inMin <- min(res[res$Cluster == 'In',]$meanTPM)
outMax <- max(res[res$Cluster == 'Out',]$meanTPM)
cut_value <- log2(mean(inMin, outMax) + 1)

# write out cutoff so I can use it in report
write.table(data.frame(CV=c(cut_value), MY=c(maxY)), file=paste(outDir, "basemean_cutoff.txt", sep=""), sep="\t", quote=F, row.names=F)



# Write Results
outResults <- data.frame(GeneID=res$GeneID, Gene=res$GeneName, baseMean=res$baseMean, stat=res$stat, log2FoldChange=res$log2FoldChange, pvalue=res$pvalue, padj=res$padj)
res <- res %>% select(GeneID, GeneName, everything())
name <- paste0(outDir, 'results.txt') 
write_tsv(res, name)



resSig <- res[!(is.na(res$padj)),] # Keep only genes that have a calculated adjusted p
resSig <- resSig[ resSig$padj < 0.05, ] # Keep adjusted p less than 0.05
#resSig <- resSig[resSig$baseMean > (2^cut_value), ] # Keep genes with enough expression
resSig <- resSig %>% filter(Cluster == 'In') # Keep only the second cluster
resSig <- resSig[abs(resSig$log2FoldChange) > 0.585, ] # Kep genes with a fold change of at least 1.5
	
write_tsv(resSig, paste0(outDir, "significant.txt"))



##########  Sanity Check
# Plot counts of most significant, to check if fold change is right
png(paste0(outDir, 'sanity.check.png'))
title=paste(res[1,]$GeneName, "\nFold Change:",res[1,]$log2FoldChange)
plotCounts(dds, gene=res[1,]$GeneID, intgroup = PCA_Group, main=res[1,]$GeneName, 
					 sub=paste('FC:', format(res[1,]$log2FoldChange, digits=2)), pch=19)
dev.off()



