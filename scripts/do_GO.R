# Error logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(clusterProfiler)


# Get Organism DB from config
organism_db <- snakemake@config$organism_db
library(organism_db, character.only = T)
orgDB <- organism_db


# This gets the prefix we need
exp <- snakemake@params$exp

outDir <- paste0("results/", exp, "/GO/")

# Get list of significant genes
sig <- read_tsv(paste0('results/', exp, '/deseq/significant.txt'))

# Check if there are any significant genes. If not, write empty result and return
if (dim(sig)[1] == 0){
	d <- tibble(no_genes = character())
	write_tsv(d,  paste0(outDir, 'BP_results.txt'))	
	return()
}

# get Entrez gene IDs
genes <- bitr(sig$GeneID, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = orgDB)
geneList <- genes$ENTREZID

# Get all genes for the universe
res <- read_tsv(paste0('results/', exp, '/deseq/results.txt'))
univ <- bitr(res$GeneID, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = orgDB)
universe <- univ$ENTREZID


# Define function to do GO analysis and make figures

myGoFunc <- function(geneList, ontology){
	
	ego <- enrichGO(gene        = geneList,
									universe      = universe,
									OrgDb         = orgDB,
									ont           = ontology,
									pAdjustMethod = "BH",
									pvalueCutoff  = 0.01,
									qvalueCutoff  = 0.05,
									readable      = TRUE)
	
	
	ego_s <- simplify(ego, by="p.adjust", select_fun=min)
	
	saveRDS(ego_s, file=paste0(outDir, ontology, '_ego_object.rds'))
	
	# Draw several plots
	p <- goplot(ego_s)
	ggsave(p, filename=paste0(outDir, ontology, '_GO_graph.png' ))
	
	p <- barplot(ego_s, showCategory = 40)
	ggsave(p, filename=paste0(outDir, ontology, '_bar.png' ))
	
	p <- dotplot(ego_s)
	ggsave(p, filename=paste0(outDir, ontology, '_dot.png' ))
	
	p <- cnetplot(ego_s)
	ggsave(p, filename=paste0(outDir, ontology, '_concept.png' ))
	
	# Write results
	res <- as_tibble(ego_s@result)
	write_tsv(res, paste0(outDir, ontology, '_results.txt'))
	
}



# Run the function using MF, CC, and BP
myGoFunc(geneList, 'MF')
myGoFunc(geneList, 'CC')
myGoFunc(geneList, 'BP')

