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
	for (p in c('sig_all_genes_', 'sig_up_genes_', 'sig_down_genes_')){
		for (ont in c('BP','MF','CC')){
#			print(paste0(outDir, p, ont, "_results.txt"))
			write_tsv(d, paste0(outDir, p, ont, "_results.txt"))
		}
	}
	#	write_tsv(d,  paste0(outDir, 'BP_results.txt'))	
#	write_tsv(d,  paste0(outDir, 'MF_results.txt'))	
#	write_tsv(d,  paste0(outDir, 'CC_results.txt'))	
	return()
}

# get Entrez gene IDs
genes <- bitr(sig$GeneID, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = orgDB)
sig <- left_join(genes, sig, by=c("ENSEMBL" = 'GeneID'))

# Get all genes for the universe
res <- read_tsv(paste0('results/', exp, '/deseq/results.txt'))
univ <- bitr(res$GeneID, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = orgDB)
universe <- univ$ENTREZID


# Define function to do GO analysis and make figures

myGoFunc <- function(geneList, prefix, ontology){
	
	ego <- enrichGO(gene          = geneList,
									universe      = universe,
									OrgDb         = orgDB,
									ont           = ontology,
									pAdjustMethod = "BH",
									pvalueCutoff  = 0.01,
									qvalueCutoff  = 0.05,
									readable      = TRUE)
	
	

	if(is.null(ego) | dim(as.data.frame(ego))[1] == 0) {
		write_tsv(tibble(no_genes = character()), paste0(outDir, prefix, ontology, '_results.txt'))
		return()
	} else {
		ego_s <- simplify(ego, by="p.adjust", select_fun=min)
		# Write results
		res <- as_tibble(ego_s@result)
		write_tsv(res, paste0(outDir, prefix, ontology, '_results.txt'))
	}
	
	saveRDS(ego_s, file=paste0(outDir, prefix, ontology, '_ego_object.rds'))
	
	# Draw several plots
	p <- goplot(ego_s)
	ggsave(p, filename=paste0(outDir, prefix, ontology, '_GO_graph.png' ))
	
	p <- barplot(ego_s, showCategory = 40)
	ggsave(p, filename=paste0(outDir, prefix, ontology, '_bar.png' ), height=5, width=10)
	
	p <- dotplot(ego_s)
	ggsave(p, filename=paste0(outDir, prefix, ontology, '_dot.png' ), height=5, width=10)
	
	p <- cnetplot(ego_s)
	ggsave(p, filename=paste0(outDir, prefix, ontology, '_concept.png' ))
	
	
}







# All
geneList <- genes$ENTREZID

for (ont in c('BP','MF','CC')){
	myGoFunc(geneList, 'all_genes_', ont)	
}


# Up
geneList <- sig %>% filter(padj < 0.05 & log2FoldChange > 0) %>% dplyr::select(ENTREZID)

for (ont in c('BP','MF','CC')){
	myGoFunc(geneList, 'all_up_genes_', ont)	
}

# Down
geneList <- sig %>% filter(padj < 0.05 & log2FoldChange < 0) %>% dplyr::select(ENTREZID)

for (ont in c('BP','MF','CC')){
	myGoFunc(geneList, 'all_down_genes_', ont)	
}




