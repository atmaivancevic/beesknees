### Install and load DESeq2.
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

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

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds

dds <- DESeq(dds)
resultsNames(dds)

### To get your differential expression results, you need to specify between what conditions you are comparing. 
#To do this, you use the ''contrast'' argument. 
#You should always reference your experimental group before your control group for your fold change to be correct. 
queen_vs_fullbody <- results(dds, contrast=c("condition", "queen", "fullbody"))

### Then, filter out all "N/A" entries using na.omit() and sort your res table by ascending p-value. 
# List the number of entries that have an adjusted p-value less than 0.05.
queen_vs_fullbody <- na.omit(queen_vs_fullbody)
queen_vs_fullbody <- queen_vs_fullbody[order(queen_vs_fullbody$padj), ]

# Separate data into three subsets

value1 = subset(queen_vs_fullbody, padj<0.05)
value2 = subset(queen_vs_fullbody, abs(log2FoldChange)>1)
value3 = subset(queen_vs_fullbody, padj<0.05 & abs(log2FoldChange)>1)

nrow(value1)
nrow(value2)
nrow(value3)

