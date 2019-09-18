# Pass organism ID in

args <- commandArgs()



install.packages('tidyverse')
install.packages('RColorBrewer')
install.packages('ggrepel')
install.packages('gplots')
install.packages()

install.packages("BiocManager")
BiocManager::install('DESeq2', update=F)
BiocManager::install('tximport', update=F)
BiocManager::install('AnnotationDbi', update=F)

if (!require(args[1], quietly = TRUE))
  BiocManager::install(args[1], update=F)


library(AnnotationDbi)
library('org.Hs.eg.db')

x <- '1'
write.table(x, file="envs/R_initialized")
