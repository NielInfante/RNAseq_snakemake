# Error logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")



organism <- snakemake@config$organism

#str(snakemake)

# Pass organism ID in

args <- commandArgs()


if (!require(organism, quietly = TRUE))
  BiocManager::install(organism, update=F)


x <- '1'
write.table(x, file="envs/R_initialized")
