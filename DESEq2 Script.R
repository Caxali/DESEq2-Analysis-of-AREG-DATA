#======================================
#   RNA-seq Analysis with DESEq2
#======================================
library(DESeq2)
library(tidyverse)

#Read your Data
data <- read.table('C260219002A-AWR.gene_counts.tsv', h = T)

#Process the data
filtered <- data[, c('gene_name', 'C1', 'C2','C3','C4','C5','C6','C7','C8')]
table(duplicated(filtered$gene_name))
filtered$gene_name <- make.unique(filtered$gene_name, sep = '.')
table(duplicated(filtered$gene_name))
rownames(filtered) <- filtered$gene_name


#Create a index to identify the transcript of duplicated genes. 
ref <- data.frame(Gene_Name = filtered$gene_name,
                  Ensembl = data$gene_id)
filtered <- filtered[,-1]

#Creating the Meta data
meta <- data.frame(Samples = colnames(filtered),
                   Group = as.factor(c(rep('AAV-GFP', 4), rep('AAV-AREG', 4))))
rownames(meta) <- meta$Samples

#Check if samples order is correct
all(colnames(filtered) %in% rownames(meta))
all(colnames(filtered) == rownames(meta))

#Both true, create a DESEq2 object
dds <- DESeqDataSetFromMatrix(countData = filtered,
                              colData = meta,
                              design = ~ Group)

dds


#Removing genes with low genes counts. Minimum 10 genes reades
keep <- rowSums(counts(dds)) >= 10
keep
dds <- dds[keep,]
dds

#Set factor level
dds$Group <- relevel(dds$Group, ref = 'AAV-GFP')
dds$Group

#Run DESEq
dds <- DESeq(dds)
res <- results(dds)
res
res <- as.data.frame(res)
res$FDR <- p.adjust(res$pvalue, method = 'BH')
res$log_10FDR <- -log10(res$FDR)

up <- res[which(res$log2FoldChange >= 0.5 & res$pvalue < 0.05),]
down <- res[which(res$log2FoldChange <= -0.5 & res$pvalue < 0.05),]

#Saving all tables
write.table(meta, 'Meta_Data.txt', sep = "\t")
write.table(ref, 'Index.txt', sep = '\t')
write.table(res, 'DESEq2 Results.txt', sep = '\t')
write.table(up , 'Up Genes.txt', sep = '\t')
write.table(down, 'Down Genes.txt', sep = '\t')



