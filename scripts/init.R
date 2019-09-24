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


install.packages('tidyverse', repos='http://cran.us.r-project.org')

#require(devtools)
#install_version("ggplot2", version = "0.9.1", repos = "http://cran.us.r-project.org")

#install.packages('RColorBrewer', repos='http://cran.us.r-project.org')
#install.packages('ggrepel', repos='http://cran.us.r-project.org')
#install.packages('gplots', repos='http://cran.us.r-project.org')





x <- '1'
write.table(x, file="envs/R_initialized")
