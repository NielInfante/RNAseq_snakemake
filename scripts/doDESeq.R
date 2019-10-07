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

# Save dds
saveRDS(dds, paste0(outDir, '/', outPrefix, '/dds.rds'))


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

# Write Results
outResults <- data.frame(GeneID=res$GeneID, Gene=res$GeneName, baseMean=res$baseMean, stat=res$stat, log2FoldChange=res$log2FoldChange, pvalue=res$pvalue, padj=res$padj)
res <- res %>% select(GeneID, GeneName, everything())
name <- paste(outDir, '/', outPrefix, '/results.txt', sep="") 
write.table(res, file=name, sep="\t", quote=F, row.names=F)

# Significant genes

## Get an appropriate base mean cutoff value
den <- density(log2(res$baseMean))
maxY <- max(den$y)
maxAt <- which.max(den$y)

for (i in (maxAt + 1):length(den$y)){
	diff <- den$y[i] - den$y[i-1]
	if (diff > 0){
		cut_index <- i - 1
		break
	}
}

if (cut_index < 0){
	cut_value <- median(log2(res$baseMean))
} else {
	cut_value <- den$x[cut_index]	
}

# write out cutoff so I can use it in report
write.table(data.frame(CV=c(cut_value), MY=c(maxY)), file=paste(outDir, "/", outPrefix, "/basemean_cutoff.txt", sep=""), sep="\t", quote=F, row.names=F)

resSig <- res[!(is.na(res$padj)),] # Keep only genes that have a calculated adjusted p
resSig <- resSig[ resSig$padj < 0.05, ] # Keep adjusted p less than 0.05
resSig <- resSig[resSig$baseMean > (2^cut_value), ] # Keep genes with enough expression

write.table(resSig,file=paste(outDir, "/", outPrefix, "/significant.txt", sep=""), sep="\t", quote=F, row.names=F)



##########  Sanity Check
# Plot counts of most significant, to check if fold change is right
png(paste0(outDir, '/',outPrefix,'/sanity.check.png'))
title=paste(res[1,]$GeneName, "\nFold Change:",res[1,]$log2FoldChange)
plotCounts(dds, gene=res[1,]$GeneID, intgroup = PCA_Group, main=res[1,]$GeneName, 
					 sub=paste('FC:', format(res[1,]$log2FoldChange, digits=2)), pch=19)
dev.off()



