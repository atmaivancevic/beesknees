### Install and load DESeq2.
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
biocLite("pheatmap")
library("pheatmap")

### Read in the output from featureCounts.
countdata <- read.table("/Users/CL_Shared/data/atma/beenomics_rna/featureCounts/featureCounts.txt", header = TRUE, skip = 1)
head(countdata)

# Set first column to be the rownames
rownames(countdata) <- countdata[,1]
head(countdata)

### Remove the first six columns (geneid, chr, start, end, strand, length).
countdata <- countdata[ ,7:ncol(countdata)]
head(countdata)

# Edit column names (e.g. to remove the dir path)
colnames(countdata) <- gsub("X.Users.CL_Shared.data.atma.beenomics_rna.bam.", "", colnames(countdata))
colnames(countdata) <- gsub("_L004.uniq.bam", "", colnames(countdata))
colnames(countdata)
head(countdata)

# Reorder sample order
countdata <- countdata[,c(1,2,3,5,4,6,7,8)]
head(countdata)

### Convert countdata table into a matrix - necessary for running DESeq2. 
countdata <- as.matrix(countdata)
head(countdata)

# Assign control vs test samples (or different conditions)
condition <- factor(c(rep("fullbody", 2), rep("body", 2), rep("head", 2), rep("brain", 1), rep("queen", 1)))

# Create a "coldata" table containing the sample names with their appropriate condition (e.g. control versus cancer sample)
coldata <- data.frame(row.names=colnames(countdata), condition)
coldata
head(countdata)

# Construct a DESeqDataSet

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds

dds <- DESeq(dds)

# print the deseq normalized read count table
normalized_countdata <- as.data.frame(counts(dds, normalized = TRUE))
head(normalized_countdata,5)
head(countdata,5)

# make a pearsons correlation matrix using normalized counts
cormat <- cor(normalized_countdata, method="pearson")
cormat
symnum(cormat)
?cor

# plot the correlation matrix
library(reshape2)
melted_cormat <- melt(cormat)
melted_cormat

library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# make a function to reorder the correlation matrix by clustering
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# now apply the reorder function
cormat <- reorder_cormat(cormat)

# now replot it, with samples appearing in the clustered order
melted_cormat <- melt(cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# can we make a correlation matrix of genes, not samples?
# erm let's try
# need to transpose the normalized count table
cormat_genes <- cor(t(normalized_countdata), method="pearson")
# oh apparently the standard deviation is 0 lol. Guess this was a bad idea.

# how about a pca plot?
# first, recommended to do a variance stabilizing transformation (like another sort of normalizing)
# https://master.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#starting-from-count-matrices
vsd <- vst(dds)
class(dds)
class(vsd)
colData(vsd)

# make a PCA plot (note: this plotPCA is part of the deseq2 package)
plotPCA(vsd, returnData=FALSE, ntop=100000)

# can also use ggplot to do a pca plot
# by using returnData=TRUE to return the data points rather than plotting them
# and then using ggplot to plot
# this gives us more control for changing things like shape of points
data <- plotPCA(vsd, returnData=TRUE, ntop=500)
data
percentVar <- round(100 * attr(data, "percentVar"))
percentVar

library("ggplot2")
ggplot(data, aes(PC1, PC2, colour=condition, shape=condition)) + geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

### To get your differential expression results, you need to specify between what conditions you are comparing. 
#To do this, you use the ''contrast'' argument. 
#Always reference experimental group before control group for fold change to be correct. 

# heatmap of the most significant genes
# here we're choosing the top 30 genes, ranking by adjusted p value
# need to set the condition comparison (e.g. queen vs body seems to be the default one it chooses)

#### QUEEN VS BODY COMPARISON
queen_vs_body <- results(dds, contrast=c("condition", "queen", "body"))
nrow(queen_vs_body)
plotMA(queen_vs_body)
queen_vs_body_mat <- assay(vsd)[head(order(queen_vs_body$padj),30), ]
queen_vs_body_mat <- queen_vs_body_mat - rowMeans(queen_vs_body_mat)
df <- as.data.frame(colData(vsd))
pheatmap(queen_vs_body_mat, annotation_col=df)
# the top 2 genes (haven't checked the rest) were also found to be DE in this paper:
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0173438
# about infected honeybees
# does that mean that Orit's bees may be infected? (which may be causing the colony collapse)
# or are these genes really important?
# need to do an igv check of these genes to confirm that they are really DE in our samples

#### QUEEN VS FULLBODY COMPARISON
# make heatmap of the top 30 DE genes
queen_vs_fullbody <- results(dds, contrast=c("condition", "queen", "fullbody"))
nrow(queen_vs_fullbody)
plotMA(queen_vs_fullbody)
queen_vs_fullbody_mat <- assay(vsd)[head(order(queen_vs_fullbody$padj),30),]
queen_vs_fullbody_mat
queen_vs_fullbody_mat <- queen_vs_fullbody_mat - rowMeans(queen_vs_fullbody_mat)
df <- as.data.frame(colData(vsd))
pheatmap(queen_vs_fullbody_mat, annotation_col=df,show_rownames=T)

#export results table to csv
queen_vs_fullbody_ordered <- queen_vs_fullbody[order(queen_vs_fullbody$padj),]
head(queen_vs_fullbody_ordered)
queen_vs_fullbody_ordered_table <- as.data.frame(queen_vs_fullbody_ordered)
queen_vs_fullbody_ordered_table
write.csv(queen_vs_fullbody_ordered_table, file="/Users/ativ2716/repos/beesknees/project_files/queen_vs_fullbody_ordered_table.csv")

#### QUEEN VS BRAIN COMPARISON
# make heatmap of the top 30 DE genes
queen_vs_brain <- results(dds, contrast=c("condition", "queen", "brain"))
nrow(queen_vs_brain)
plotMA(queen_vs_brain)
queen_vs_brain_mat <- assay(vsd)[head(order(queen_vs_brain$padj),30), ]
queen_vs_brain_mat
queen_vs_brain_mat <- queen_vs_brain_mat - rowMeans(queen_vs_brain_mat)
df <- as.data.frame(colData(vsd))
pheatmap(queen_vs_brain_mat, annotation_col=df)

#export results table to csv
queen_vs_brain_ordered <- queen_vs_brain[order(queen_vs_brain$padj),]
head(queen_vs_brain_ordered)
queen_vs_brain_ordered_table <- as.data.frame(queen_vs_brain_ordered)
queen_vs_brain_ordered_table
write.csv(queen_vs_brain_ordered_table, file="/Users/ativ2716/repos/beesknees/project_files/queen_vs_brain_ordered_table.csv")

#### HEAD VS BODY COMPARISON
# make heatmap of the top 30 DE genes
head_vs_body <- results(dds, contrast=c("condition", "head", "body"))
nrow(head_vs_body)
plotMA(head_vs_body)
head_vs_body_mat <- assay(vsd)[head(order(head_vs_body$padj),1000), ]
head_vs_body_mat
head_vs_body_mat <- head_vs_body_mat - rowMeans(head_vs_body_mat)
df <- as.data.frame(colData(vsd))
pheatmap(head_vs_body_mat, annotation_col=df,show_rownames=F)
?pheatmap

#export results table to csv
head_vs_body_ordered <- head_vs_body[order(head_vs_body$padj),]
head(head_vs_body)
head_vs_body_ordered_table <- as.data.frame(head_vs_body_ordered)
head_vs_body_ordered_table
write.csv(head_vs_body_ordered_table, file="/Users/ativ2716/repos/beesknees/project_files/head_vs_body_ordered_table.csv")

# if you want to zoom into a certain gene, deseq also provides a plotCounts function
# e.g. for gene GB43248 (identified as DE in the head vs body comparison):
plotCounts(dds, gene="GB43248", intgroup="condition")

# next step for bees:
# gene ontology analysis (e.g. using CL_Wiki entry for biomart)
# pathway analysis (beebase?) or using drosophila ortholog genes 
# look at what other bee papers did for help

