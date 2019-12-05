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

p <- dotplot(kegg_over)
ggsave(p, filename=paste0(outDir, "over_overview_dot.png"))
p <- barplot(kegg_over)
ggsave(p, filename=paste0(outDir, "over_overview_bar.png"))

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
genes <- genes %>% arrange(desc(stat))
geneList <- genes$stat
#genes <- genes %>% arrange(desc(-1*log10(padj)))
#geneList <- -1 * log10(genes$padj)
#geneList <- abs(genes$log2FoldChange)
#geneList <- genes$log2FoldChange
#geneList <- abs(genes$stat)

names(geneList) <- genes$ENTREZID

kegg_gsea <- gseKEGG(geneList = geneList, keyType='ncbi-geneid',organism=kegg_db, nPerm=1000, 
										 minGSSize=50, pvalueCutoff=0.05, verbose=F)


kegg_res <- as_tibble(kegg_gsea@result)

write_tsv(kegg_res, paste0(outDir, 'KEGG_GSEA_results.txt'))

# Do some plots
if (nrow(kegg_res) == 0){
	# No significant pathways
	return
}

# Overview Plots
#kegg_res <- kegg_res %>% mutate(GeneRatio = str_extract(leading_edge, "[0-9]+"))
kegg_res <- kegg_res %>% mutate(GeneRatio = ( (1 + str_count(core_enrichment,"/")) / setSize))
kegg_res$Description <- as.factor(kegg_res$Description)
kegg_res <- kegg_res %>%  arrange(GeneRatio)
p <- ggplot(kegg_res, aes(x=GeneRatio, y= fct_reorder(Description, GeneRatio))) + 
	geom_point(aes(color=p.adjust, size=setSize)) +
	labs(y="", x="Ratio of Significant Genes", title='KEGG Pathways', color="p value", size="Genes in\nPathway") +
	theme_minimal()
ggsave(p, filename=paste0(outDir, "gsea_overview_ratio.png"))

p <- ggplot(kegg_res, aes(reorder(Description, NES), NES)) +
	geom_col() +
	coord_flip() +
	labs(x="", y="Normalized Enrichment Score",
			 title="KEGG Pathways") + 
	theme_minimal()
ggsave(p, filename=paste0(outDir, "gsea_overview_nes.png"))


# Pathway level plots
# Don't want absolute FC for figure
geneList <- genes$log2FoldChange
names(geneList) <- genes$ENTREZID

# Need to be in output directory, as pathview automatically writes figures to working directory.
setwd(outDir)


# Get genes in each pathway
path_list  <- KEGGREST::keggLink("pathway", kegg_db) %>% 
	tibble(pathway = ., ENTREZID = sub(paste0(kegg_db, ":"), "", names(.)))

path_list <- path_list %>% mutate(pathway=substr(pathway, 6, 99999))


# Loop through each significant pathway
for (idx in 1:nrow(kegg_res)){
	print(paste("Doing num", idx))
	
	# This will error for some pathways, so wrap it in try catch
	tmp <- tryCatch(
		expr = {
			pathview(gene.data=geneList, pathway.id=kegg_res$ID[idx], species=kegg_db,
						 limit=list(gene=max(abs(geneList)), cpd=1),
						 out.suffix='kegg_gsea')
		},
		error=function(e){
			print(paste("Pathview error:",e,"\nHappened with pathway", kegg_res$ID[idx]))
		}
	)

	p <- enrichplot::gseaplot2(kegg_gsea, geneSetID = kegg_res$ID[idx],  title = kegg_res$Description[idx])
	ggsave(p, filename=paste0('kegg_gsea_', kegg_res$ID[idx], '_running.png'))

	
	pg <- path_list %>% filter(pathway == kegg_res$ID[idx])
	pathway.genes <- genes %>% filter(ENTREZID %in% pg$ENTREZID)
	
#	print(paste("pathway has", dim(pathway.genes)[1], "rows, with", names(pathway.genes)))	
	p <- ggplot(pathway.genes, aes(GeneName, log2FoldChange)) +
		geom_col(aes(fill=padj<0.05)) +
		coord_flip() + 
		labs(x="Gene", y="Log 2 Fold Change", title=kegg_res$Description[idx], fill="Significant") + 
		theme_minimal() + scale_fill_manual(values=c('grey','red'))
	ggsave(p, filename=paste0("bar_", kegg_res$ID[idx], ".png"))
		
	
}






