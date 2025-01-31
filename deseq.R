# Before running, ensure that the libraries have already been downloaded and installed. For installation refer to Bioconductor

library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(AnnotationDbi)

# Set Input file
summary_feature_counts = read.table("summary_feature_counts_new.txt", header=TRUE, row.names = 1)

# Remove columns that were not going to be used
features_counts_res <- subset(summary_feature_counts, select = setdiff(names(summary_feature_counts), c("Chr", "Start", "End", "Strand", "Length")))

# Create dataframe conditions
colData <- data.frame(
  row.names = colnames(features_counts_res),
  condition = c("Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case",
                "Lung_WT_Case", "Lung_WT_Control", "Lung_WT_Control", "Lung_WT_Control",
                "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case",
                "Blood_WT_Control", "Blood_WT_Control", "Blood_WT_Control")
)

# Convert condition to categorical 
colData$condition = as.factor(colData$condition)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = features_counts_res,
                              colData = colData,
                              design = ~ condition)

# Run DESeq
dds <- DESeq(dds)

# Run variance stabilizing transformation
vsd = vst(dds, blind = TRUE)

# Plot the PCA 
plotPCA(vsd, intgroup="condition")

# Get results for differential expression
res_blood <- results(dds, contrast=c("condition","Blood_WT_Case","Blood_WT_Control"))
res_blood

## Convert ENSEMBL ID to Gene symbols ## 

#Extract the IDs from the rownames of DESeq object
ids <- rownames(res_blood)

#Map the IDs to respective gene symbols
symbols <- mapIds(org.Mm.eg.db, 
                  keys = ids, 
                  column = "SYMBOL", 
                  keytype = "ENSEMBL")

#Replace all NA values with the original IDs for genes that does not have symbols available
symbols[is.na(symbols)] <- rownames(res_blood)[is.na(symbols)]

#Reorder the symbols
symbols <- symbols[match(rownames(res_blood), names(symbols))]

#Replace ids (each rownames) of the DESeq object with gene symbols
rownames(res_blood) <- symbols



## Volcano Plot ##

# Remove NA values
res_clean <- na.omit(res_blood)

# Separate upregulated and downregulated genes
upregulated <- res_clean[res_clean$log2FoldChange > 0, ]
downregulated <- res_clean[res_clean$log2FoldChange < 0, ]

# Sort upregulated genes by p-value and select top 5
top_upregulated <- rownames(upregulated)[order(upregulated$pvalue)][1:5]

# Sort downregulated genes by p-value and select top 5
top_downregulated <- rownames(downregulated)[order(downregulated$pvalue)][1:5]

# Combine the top upregulated and downregulated genes
top_genes <- c(top_upregulated, top_downregulated)

# Plot the volcano plot
EnhancedVolcano(res_clean,
                lab = rownames(res_clean),
                title = "Volcano Plot of DEGs between Infected and Control Blood Samples",
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0,
                selectLab = top_genes,
                labSize = 5,
                max.overlaps = 20,
                drawConnectors = TRUE,
)



## Set the genes with high significant difference and find the downregulated and upregulated ##

# Set genes that are are differentially expressed 
res_blood_selected <- res_blood[which(res_blood$padj<0.05), ]

# Set genes that are upregulated
res_blood_upreg <- res_blood_selected[which(res_blood_selected$log2FoldChange>0), ]

# Set genes that are downregulated
res_blood_downreg <- res_blood_selected[which(res_blood_selected$log2FoldChange<0), ]

# Print how much total genes that are differentially expressed
nrow(res_blood_selected)

# Print how much genes that are upregulated
nrow(res_blood_upreg)

# Print how much genes that are downregulated
nrow(res_blood_downreg)



## Sanity check ##

# Selected genes from original publication
target_gene <- c("Cxcl9", "Gbp2", "Batf2")

# Iterate over each gene in the target_gene vector
for (gene in target_gene) {
  # Check if the gene is present in the rownames of the dataset
  if (gene %in% rownames(res_blood_selected)) {
    # Retrieve the log2FoldChange value for the gene
    log2fc <- res_blood_selected[gene, "log2FoldChange"]
    
    # Print the result
    print(paste("The gene", gene, "is present in the dataset with log2FoldChange of", log2fc))
  } else {
    # If the gene is not present
    print(paste("The gene", gene, "is NOT present in the dataset."))
  }
}



## Overrepresentation analysis ##

# Place the genes that are upregulated to a variable
blood_upreg_names = rownames(res_blood_upreg)

# Place all genes to a variable
all_names = rownames(features_counts_res)

# Run enrichment analysis for BP, CC, and MF gene ontologies
ego_bp_blood <- enrichGO(gene  = blood_upreg_names, # List of selected genes 
                universe      = all_names, # List of all genes
                keyType       = "SYMBOL",
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                readable      = TRUE)

ego_cc_blood <- enrichGO(gene  = blood_upreg_names, # List of selected genes 
                universe      = all_names, # List of all genes
                keyType       = "ENSEMBL",
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                readable      = TRUE)

ego_mf_blood <- enrichGO(gene  = blood_upreg_names, # List of selected genes 
                universe      = all_names, # List of all genes
                keyType       = "ENSEMBL",
                OrgDb         = org.Mm.eg.db,
                ont           = "MF",
                readable      = TRUE)

# Plot the tree plot of each of the ontologies
ego_bp_lung_sym <- setReadable(ego_bp_lung, 'org.Mm.eg.db', 'ENSEMBL')
ego_bp_lung_sym2 <- pairwise_termsim(ego_bp_lung_sym)
p1 <- treeplot(ego_bp_lung_sym2)

ego_cc_lung_sym <- setReadable(ego_cc_lung, 'org.Mm.eg.db', 'ENSEMBL')
ego_cc_lung_sym2 <- pairwise_termsim(ego_cc_lung_sym)
p2 <- treeplot(ego_cc_lung_sym2)

ego_mf_lung_sym <- setReadable(ego_mf_lung, 'org.Mm.eg.db', 'ENSEMBL')
ego_mf_lung_sym2 <- pairwise_termsim(ego_mf_lung_sym)
p3 <- treeplot(ego_mf_lung_sym2)




