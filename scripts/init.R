# Error logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


# This allows me to install an R package that is 
# specified in the config.yaml file.
# This way, the user doesn't have to edit the environment files.
# Also, this script should run only once, and not try to install the
# same package multiple times.

organism_db <- snakemake@config$organism_db

#str(snakemake)


if (!require(organism_db, quietly = TRUE))
  BiocManager::install(organism_db, update=F)


# Indicator to say the script ran, and snakemake doesn't call it again.
x <- '1'
write.table(x, file="envs/R_initialized")
