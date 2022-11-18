# Differential Gene Expression - STAR
setwd("~/differential_expression") # Navigate to where you put your stuff
install.packages("tidyverse")
install.packages("DESeq2")
install.packages("pheatmap")
install.packages("ggrepel")

install.packages("BiocManager")
BiocManager::install("apeglm")
BiocManager:: install("biomaRt")






library(DESeq2) # Then load up DESeq2


# 1. IMPORTING THE DATA -----------------
# Importing the meta data table



# Create the DESeq2 object



# 2. FILTERING -----------------
### To start, we should filter out lowly expressing genes. These take up computing power and lower the speed of the analysis.





# 3. QUALITY CONTROL -----------------
## PCA




## Hierarchical clustering
library(pheatmap)




# 4. ANNOTATION -----------------
library(biomaRt)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://mar2017.archive.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")



annot <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "description", "external_gene_name"),
               filters = "ensembl_gene_id", values = gene_ids, mart = mart)



# 5. DIFFERENTIAL GENE EXPRESSION ANALYSIS -----------------
## Recreate DESeq2 object



## Remember to refilter




## Getting results
library(apeglm)




# 6. VISUALISATION -----------------






library(ggrepel)


