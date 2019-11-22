# Error logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(clusterProfiler)
library(pathview)

# Get organism DB and KEGG DB from config
kegg_db <- snakemake@config$kegg_db
organism_db <- snakemake@config$organism_db
library(organism_db, character.only = T)
orgDB <- organism_db


# This gets the prefix we need
exp <- snakemake@params$exp

outDir <- paste0("results/", exp, "/KEGG/")

# Get list of significant genes
sig <- read_tsv(paste0('results/', exp, '/deseq/significant.txt'))

# Check if there are any significant genes. If not, write empty result and return
if (dim(sig)[1] == 0){
	d <- tibble(no_genes = character())
	write_tsv(d,  paste0(outDir, 'KEGG_results.txt'))	
	write_tsv(d,  paste0(outDir, 'KEGG_GSEA_results.txt'))	
	return()
}

# get Entrez gene IDs
genes <- bitr(sig$GeneID, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = orgDB)
dim(genes)
genes <- left_join(genes, sig, by=c("ENSEMBL" = 'GeneID'))
geneList <- genes$ENTREZID

# Overrepresentation Test
kegg_over <- enrichKEGG(gene=geneList, organism=kegg_db, pvalueCutoff = 0.05)

print(paste("Got", nrow(kegg_over@result), 'sig stuff'))
kegg_res <- kegg_over@result
kegg_res <- filter(kegg_res, p.adjust < 0.05)

write_tsv(kegg_res, paste0(outDir, 'KEGG_results.txt'))


# GSEA

# Get list of all genes
res <- read_tsv(paste0('results/', exp, '/deseq/results.txt'))

# get Entrez gene IDs
genes <- bitr(res$GeneID, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = orgDB)
genes <- left_join(genes, res, by=c("ENSEMBL" = 'GeneID'))
genes <- na.omit(genes)

# Pick how to sort
#genes <- genes %>% arrange(desc(abs(log2FoldChange)))
#genes <- genes %>% arrange(desc(log2FoldChange))
#genes <- genes %>% arrange(desc(abs(stat)))
#genes <- genes %>% arrange(desc(stat))
genes <- genes %>% arrange(desc(-1*log10(padj)))
geneList <- -1 * log10(genes$padj)
#geneList <- abs(genes$log2FoldChange)
#geneList <- genes$log2FoldChange
#geneList <- abs(genes$stat)
#geneList <- genes$stat

names(geneList) <- genes$ENTREZID

kegg_gsea <- gseKEGG(geneList = geneList, keyType='ncbi-geneid',organism=kegg_db, nPerm=1000, 
										 minGSSize=120, pvalueCutoff=0.05, verbose=F)



kegg_res <- kegg_gsea@result

write_tsv(kegg_res, paste0(outDir, 'KEGG_GSEA_results.txt'))

# Do some plots
if (nrow(kegg_res) == 0){
	# No significant pathways
	return
}

# Don't want absolute FC for figure
geneList <- genes$log2FoldChange
names(geneList) <- genes$ENTREZID

# Need to be in output directory, as pathview automatically writes figures to working directory.
setwd(outDir)

# Loop through each significant pathway
for (idx in 1:nrow(kegg_res)){
	
	pathview(gene.data=geneList, pathway.id=kegg_res$ID, species=kegg_db,
					 limit=list(gene=max(abs(geneList)), cpd=1),
					 out.suffix='kegg_gsea')

	p <- enrichplot::gseaplot2(kegg_res, geneSetID = idx,  title = kegg_res$Description[idx])
	ggsave(p, file=paste('kegg_gsea_', kegg_res$ID, '_running.png'))
	
}



