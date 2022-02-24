# Installations
if(FALSE){
  
  # Install package for general analysis
  install.packages("pheatmap")
  install.packages("data.table")
  install.packages("dplyr")
  
}

library(dplyr)
library(data.table)
library(pheatmap)

tab <- fread("../data/GSE74912_ATACseq_All_Counts.txt.gz")
meta <- fread("../data/samplelookup.tsv")[1:80,]

# Some ugly code just to subset to the samples that we have in the meta data
counts_mat <- data.matrix(data.frame(tab[,meta$SampleName, with = FALSE]) )

# Make the names pretty
colnames(counts_mat) <- make.unique(meta$CellType)

pheatmap(cor(counts_mat))
