# Installations
if(FALSE){
  
  # Install package for support with gene names / coordinates
  install.packages("devtools")
  devtools::install_github("stephenturner/annotables")
  
  # install DEseq2
  install.packages("BiocManager")
  BiocManager::install("DESeq2") # for differential analysis
  BiocManager::install("GenomicRanges") # for gene annotation
  BiocManager::install("apeglm") # for volcano plot fold changes
  
  
}

library(annotables)
library(DESeq2)
library(GenomicRanges)
library(dplyr)
library(data.table)

# Build a simplified data frame of gene coordinates
data.frame(
  chr = paste0("chr", grch37$chr), 
  start = grch37$start,
  end = grch37$end,
  gene = grch37$symbol
) %>% makeGRangesFromDataFrame() -> gene_coordinates

# Import the ATAC-seq data
tab <- fread("../data/GSE74912_ATACseq_All_Counts.txt.gz")
meta <- fread("../data/samplelookup.tsv") %>%
  filter(CellType %in% c("CD8Tcell", "NKcell"))

# Some ugly code just to subset to the samples that we have in the meta data
counts_mat <- data.matrix(data.frame(tab[,meta$SampleName, with = FALSE]) )

# Now run DEseq2
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = data.frame(meta),
                              design= ~ CellType) 
dds <- DESeq(dds) # this takes a couple of minutes
resultsNames(dds) # lists the coefficients
res <- results(dds, name="CellType_NKcell_vs_CD8Tcell")

# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="CellType_NKcell_vs_CD8Tcell", type="apeglm")

# Now put all of the results into a dataframe
data.frame(
  res
)