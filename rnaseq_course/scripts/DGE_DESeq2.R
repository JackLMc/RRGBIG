# Differential Gene Expression - STAR
setwd("./rnaseq_course/differential_expression") # So we'll first navigate to the differential gene expression folder that I gave you all.
library(DESeq2) # Then load up DESeq2

# 1. IMPORTING THE DATA -----------------
# Importing the meta data table
### To start with, we'll load in the meta data table, with header as true as this contains our column names
`sampletable <- read.table("sample_sheet_foxc1.txt", header = T, sep = "\t")`

### We'll next use the filename column to add the SRA codes as the rownames. This is just so they're there for safe keeping.
`rownames(sampletable) <- gsub("_counts.txt", "", sampletable$FileName)`

### A quick check the number of rows and columns
`nrow(sampletable)` # Should be 10  
`ncol(sampletable)` # Should be 4


### If we have a look at the top of the meta data, we see there are two parameters that could be of interest. 
### For this, lets say we're interested in the Condition.
`View(sampletable)`


### Our next move is to load in the count data from STAR into a DESeq2 object. These objects are like big lists, but instead of us having to make the separate bits, it has pre-made slots for us.
se_star <- DESeqDataSetFromHTSeqCount(sampleTable = sampletable,
                                      directory = "STAR_counts",
                                      design = ~ Condition) 
### The design is the column we're interested in and is evaluate gene expression changes with respect to those different levels


# 2. FILTERING --------
### To start, we should filter out lowly expressing genes. These take up computing power and lower the speed of the analysis.
### This isn't vital before, but it's good practice.

### Let's have a look at how many genes are there before filtering
nrow(se_star) # Quite a lot. Remember these are ensembl IDs and not symbols like we all know and love


### So to filter, we'll do it based on the individual genes. Usually, if a gene has a count lower than 10 across all samples, it's quite low. So we'll filter the matrix for those genes which have a count of above 10 across all the samples.
se_star <- se_star[rowSums(counts(se_star)) > 10, ]

nrow(se_star) ### now if we look at the number of genes, there were quite a lot lost!
# se_star[rowSums(counts(se_star)) == 0, ]


# 3. QUALITY CONTROL -------------
### Firstly, we need to normalise the counts like I went through in the seminar last week. Well, we don't actually need to do this for the quality control workflow, all we need is the Size Factors which will be used in the transformation step so we can visualise the data a bit better. 
### To do this, we first estimate the size factors. We can put this back into the dds object since there is a slot for it.
se_star <- estimateSizeFactors(se_star)



### The transformation step uses regularised logarithm of the counts. This step just helps with visualising the data a bit better. We can do this by using the `rlog()` function. 
### This function can be quite slow if doing it for > 20 samples, so if you are performing it on more than 20 then use the `vst()` function instead which performs similar stuff, but is a lot quicker.
### We also set an argument (blind) to TRUE here, which performs the transformation blind of the sample condition. 
### This is important to do in quality control since we don't want any transformations that increase the similarity between the sample conditions.
rld <- rlog(se_star, blind = T)

### Quality control is usually completed using two different measures. PCA and hierarchical clustering. To perform PCA, DESeq has an inbuilt function.
plotPCA(rld, intgroup = "Condition")
### When using the condition as the input group, we see that it doesn't separate the samples all that well. So perhaps there are more meaningful comparisons using other measures.
plotPCA(rld, intgroup = "Differentiation")
### Which there are if we use differentiation. We can use them together. And we see quite interesting things when we combine. 
plotPCA(rld, intgroup = c("Condition", "Differentiation"))



### The other method is by hierarchical clustering.
### DESeq2 doesn't have an inbuilt function for heatmaps, so we're going to use the `pheatmap()` function from the pheatmap package. 
### This package requires a matrix as it's input, so we have to retrieve that from the rlog transformed values first
### We can do this using a function from SummarisedExperiment, which loads when you load up DESeq2
rld_mat <- assay(rld)

### We then find the correlations between samples. This is how well they correlate with one another in terms of all genes. This can be achieved using the `cor()` function in base R
rld_cor <- cor(rld_mat)
library(pheatmap)
pheatmap(rld_cor)

### We can also prepare a metadata object to add to the heatmap
metadata <- sampletable[,c("Differentiation", "Condition")]
rownames(metadata) <- sampletable$SampleName
pheatmap(rld_cor, annotation_col = metadata)

### Based on these two metrics, we decide whether our data is good enough to go forward with. Ideally, we want our samples from the same group clustering together. 

# 4. ANNOTATION -------------
### Once we've done the clustering stuff to check how nice the data looks. i.e. that the conditions largely stick together. We can do differential gene expression analysis.
### A sort of prerequisite for this, so we know the genes a bit better is to annotate them with gene symbols.
### It doesn't have to be done, but it makes for more friendly viewing.
### We can use the biomaRt package to do so.
# BiocManager::install("biomaRt")
library(biomaRt)


### So this package is a collection of annotations for all the various versions 
### You can list the various genome version annotations using 
listEnsemblArchives() ### When we aligned, we used version 88
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://mar2017.archive.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")


### We then have to choose "filters", which correspond to the input WE want to retrieve more annotation for. In our case this is the ensembl IDs
# a list of available filters can be obtained with listFilters(mart)
head(listFilters(mart))


# We can see what is available in terms of ensembl ID
grep("ensembl", listFilters(mart)[,1], value = T)
# ensembl_gene_id is what we want!


### We then have to choose "attributes." These correspond to the kind of annotation you want to retrieve from it. 
# a list of available attributes can be obtained with listAttributes(mart)
head(listAttributes(mart))

# list of ENSEMBL IDs we want to annotate
gene_ids <- rownames(se_star)

# annotate, we'll take out the gene ID, chromosome, start position, end position, a description of the gene, and the gene name.
annot <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "description", "external_gene_name"),
               filters = "ensembl_gene_id", values = gene_ids, mart = mart)

# Now we have our annotation matrix for all genes which are in our experiment that passed the filter
head(annot)


# 5. DIFFERENTIAL GENE EXPRESSION ANALYSIS -----------------
### Firstly, recreate the matrix. Just in case you had outliers or anything. Make sure to remove the raw file from the directory.
se_star <- DESeqDataSetFromHTSeqCount(sampleTable = sampletable,
                                      directory = "STAR_counts",
                                      design = ~ Condition) 

### And if you do do this, don't forget to refilter the low genes
se_star <- se_star[rowSums(counts(se_star)) > 10, ]

### After finding out that Differentiation status impacted the variance, we could add in Differentiation by doing this:
### But we'll leave it be for now. 

### Then to run DESeq2, it's all one function
se_star <- DESeq(se_star) 
# Which does everything that I went over in the seminar last week, so I won't cover them again.

### You can grab some of the 'behind the scenes' data with bespoke functions
# Plot dispersion estimates
plotDispEsts(se_star)


### To get the output of the differential gene expression, we now have to build results tables for the various comparisons that we're interested in. 
### For this experiment we're only interested in WT vs KO
### DESeq2 automatically uses the levels as the input for the comparisons. 
levels(se_star$Condition)

### It doesn't really matter in this instance. But if there were a third situation, like a knock-in of the gene. We'd be comparing it to the KO in this situation, whereas we'd probably want to compare it to the WT
### So if we wanted to change it, we could use the relevel function
relevel(se_star$Condition, ref = "WT")



### To build a results table we use the `results()` function. 
### To tell DESeq2 which groups we want to compare, we have to supply the contrasts we would like using the `contrasts` argument. 
### The alpha argument here, is what P-value we determine as significant. By default, this is 0.1 but we'll change this to 0.05.
contrasts_ko <- c("Condition", "WT", "KO") ### So DESeq2 doesn't perform the shrinkage stuff I was on about last week by default anymore

### So we'll save it as unshrunken.
res_tableKO_unshrunken <- results(se_star, contrast = contrasts_ko, alpha = 0.05)

### We'll perform the shrinkage stuff using the `lfcShrink()` function. We first have to get what it's called in the table, though
resultsNames(se_star)
library(apeglm)
res_tableKO <- lfcShrink(se_star, coef = "Condition_WT_vs_KO" , res = res_tableKO_unshrunken)


### We can next have a quick look at how the results table looks
library(tidyverse)
res_tableKO %>% data.frame() %>% head() %>% View()
### So although we've done some filtering ourselves, DESeq adds in more filtering. And those with a low baseMean will have their adjusted P-value set to NA

### We can find out what the various columns mean by using the `mcols` function
mcols(res_tableKO, use.names = T)

### Which apparently doesn't fit the screen.
### baseMean = mean of the normalised counts across all samples
### log2FoldChange = change in expression level between the groups. In this case, because we have the WT first, a positive value will indicate enrichment in the WT group
### P value is the p value
### adjusted P value is adjusted with benjamini hochberg by default
### lfcSE is the standard error of the fold change
### stat is the Wald statistic, which is the log fold change divided by its standard error

## We're also viewing the ensembl genes, which perhaps aren't that informative. We can use the annotation table we created earlier to find out a bit more
de_symbols <- merge(data.frame(ID = rownames(res_tableKO), res_tableKO, check.names = F), annot, by.x = "ID", by.y = "ensembl_gene_id", all=F)
head(de_symbols)


### We can filter the object for the significant values, using both the adjusted P value. And if you want to be more stringent, the log fold change.
### As we're working in log2() transformed values, a 0.58 cutoff will mean a 1.5 fold change between them
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

sigOE <- de_symbols %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

### To order these in terms of p value, we can write
sigOE[order(sigOE$padj),] %>% head()

### Or to see if a certain gene is significant
sigOE[sigOE$external_gene_name == "FOXC1", ]

# 6. VISUALISATION ------------
### Now we have our list of differentially expressed genes, we can visualise these. Perhaps the most famous method of doing so is using a volcano plot
### These plots have the log-transformed adjusted P-values on the y axis and the log2 fold change values on the x axis
### To generate these, we need a column in our results data indicating whether or not the gene is considered differentially expressed based on adjusted p-values
de_symbols <- de_symbols[!is.na(de_symbols$padj),]
res_tableKO_tb <- de_symbols %>% 
  mutate(threshold_KO = padj < 0.05 & abs(log2FoldChange) >= 0.58)



### We can then create a ggplot object, demonstrating this
ggplot(res_tableKO_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_KO)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value")

### Which isn't super fancy. But it demonstrates the point
### You can then annotate the graph for the top 10 most differentially expressed genes 
## Create a column to indicate which genes to label
res_tableKO_tb <- res_tableKO_tb %>% arrange(padj) %>% mutate(genelabels = "") ### So we first arrange by adjusted P value and create a column for the gene names
res_tableKO_tb$genelabels[1:10] <- res_tableKO_tb$external_gene_name[1:10] ### We then populate the genelabel column with the external gene name


library(ggrepel)
ggplot(res_tableKO_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_KO)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(aes(label = genelabels))
### And add this on to the end of the ggplot object

### So this was for the Condition variable, what happens if we do it for the Differentiation status? 

