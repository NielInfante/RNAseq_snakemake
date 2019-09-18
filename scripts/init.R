# Pass organism ID in

args <- commandArgs()



install.packages('tidyverse', repos='http://cran.us.r-project.org')
install.packages('RColorBrewer', repos='http://cran.us.r-project.org')
install.packages('ggrepel', repos='http://cran.us.r-project.org')
install.packages('gplots', repos='http://cran.us.r-project.org')

install.packages("BiocManager", repos='http://cran.us.r-project.org')
BiocManager::install('DESeq2', update=F)
BiocManager::install('tximport', update=F)
BiocManager::install('AnnotationDbi', update=F)

if (!require(args[1], quietly = TRUE))
  BiocManager::install(args[1], update=F)


library(AnnotationDbi)
library("org.Hs.eg.db")

x <- '1'
write.table(x, file="envs/R_initialized")
