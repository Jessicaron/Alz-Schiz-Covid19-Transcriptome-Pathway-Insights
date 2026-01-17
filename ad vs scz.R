# ============================================================
# Full workflow for GSE138082 (hippocampus) & GSE270454 (blood)
# ============================================================

# ----------------------------
# Step 0: Load libraries
# ----------------------------
suppressMessages({
  if(!require(DESeq2)) BiocManager::install("DESeq2")
  if(!require(dplyr)) install.packages("dplyr")
  if(!require(clusterProfiler)) BiocManager::install("clusterProfiler")
  if(!require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")
  if(!require(ggplot2)) install.packages("ggplot2")
  if(!require(reshape2)) install.packages("reshape2")
  library(DESeq2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(reshape2)
})

# ----------------------------
# Step 1: Load count matrices
# ----------------------------
counts1_file <- "C:/Users/Jessica/Downloads/GSE138082_Subfields_Counts.txt.gz"
counts2_file <- "C:/Users/Jessica/Downloads/GSE270454_RNAseq-combined-counts-matrix.csv.gz"

counts1 <- read.delim(counts1_file, header = TRUE, stringsAsFactors = FALSE)
rownames(counts1) <- counts1[,1]; counts1 <- counts1[,-1]

counts2 <- read.csv(counts2_file, header = TRUE, stringsAsFactors = FALSE)
rownames(counts2) <- counts2[,1]; counts2 <- counts2[,-1]

# Make sure columns are unique
colnames(counts1) <- make.names(colnames(counts1))
colnames(counts2) <- make.names(colnames(counts2))

# Keep only common genes for combined analysis
common_genes <- intersect(rownames(counts1), rownames(counts2))
counts1 <- counts1[common_genes,,drop=FALSE]
counts2 <- counts2[common_genes,,drop=FALSE]

# ----------------------------
# Step 2: Create phenotype tables
# ----------------------------
# GSE138082: hippocampus SZ vs CTL
n1 <- ncol(counts1)
pheno1 <- data.frame(
  sample = colnames(counts1),
  condition = factor(c(rep("CTL",12), rep("SZ", n1-12))), # adjust if needed
  dataset = "GSE138082"
)
rownames(pheno1) <- pheno1$sample

# GSE270454: blood RNA-seq
n2 <- ncol(counts2)
# Example: 11 ASM, 14 ASO, 10 MCI, remaining AD
pheno2 <- data.frame(
  sample = colnames(counts2),
  condition = factor(c(rep("ASM",11), rep("ASO",14), rep("MCI",10), rep("AD", n2-35))),
  dataset = "GSE270454"
)
rownames(pheno2) <- pheno2$sample

# ----------------------------
# Step 3: DESeq2 for hippocampus (SZ vs CTL)
# ----------------------------
dds1 <- DESeqDataSetFromMatrix(
  countData = counts1,
  colData = pheno1,
  design = ~ condition
)
dds1 <- dds1[rowSums(counts(dds1)) > 10, ]
dds1 <- DESeq(dds1)

res_SZ <- results(dds1, contrast = c("condition","SZ","CTL"))
res_SZ <- as.data.frame(res_SZ[order(res_SZ$padj),])
write.csv(res_SZ, "DEG_SZ_vs_CTL.csv")

# ----------------------------
# Step 4: DESeq2 for blood (AD vs ASM)
# ----------------------------
dds2 <- DESeqDataSetFromMatrix(
  countData = counts2,
  colData = pheno2,
  design = ~ condition
)
dds2 <- dds2[rowSums(counts(dds2)) > 10, ]
dds2 <- DESeq(dds2)

res_AD <- results(dds2, contrast = c("condition","AD","ASM"))
res_AD <- as.data.frame(res_AD[order(res_AD$padj),])
write.csv(res_AD, "DEG_AD_vs_ASM.csv")

# ----------------------------
# Step 5: KEGG enrichment function
# ----------------------------
run_kegg <- function(res, name){
  sig <- res %>% filter(!is.na(padj) & padj < 0.05)
  if(nrow(sig) == 0){
    cat("No significant genes for", name, "- skipping KEGG\n")
    return(NULL)
  }
  entrez <- mapIds(org.Hs.eg.db, keys = rownames(sig),
                   column = "ENTREZID", keytype = "SYMBOL", multiVals = "first") %>% na.omit()
  kegg <- enrichKEGG(gene = entrez, organism = "hsa", pvalueCutoff = 0.05)
  if(!is.null(kegg)){
    write.csv(kegg@result, paste0("KEGG_",name,".csv"), row.names = FALSE)
    cat(name, "KEGG pathways:", nrow(kegg@result), "\n")
  }
}

# Run KEGG
run_kegg(res_SZ, "SZ_vs_CTL")
run_kegg(res_AD, "AD_vs_ASM")

# ----------------------------
# Step 6: Plot top DEGs
# ----------------------------
plot_top_genes <- function(res, counts_matrix, filename, top_n=10){
  # Get top genes
  top_genes <- head(rownames(res[order(res$padj),]), top_n)
  genes <- intersect(top_genes, rownames(counts_matrix))
  
  if(length(genes)==0) {
    cat("No matching genes for", filename,"\n")
    return(NULL)
  }
  
  # Extract counts safely and force matrix format
  mat <- as.matrix(counts_matrix[genes,,drop=FALSE])
  mat <- log2(mat + 1)
  
  # Melt safely
  mat_long <- reshape2::melt(mat)
  colnames(mat_long) <- c("Gene","Sample","log2count")  # safe now
  
  # Plot
  p <- ggplot(mat_long, aes(x=Sample, y=log2count, fill=Gene)) +
    geom_bar(stat="identity", position="dodge") +
    theme(axis.text.x=element_text(angle=90,hjust=1)) +
    labs(title=paste("Top genes:", filename), y="log2(count+1)", x="Samples")
  
  print(p)
  ggsave(paste0(filename,".png"), plot=p, width=10, height=6)
  cat("Plot saved:", paste0(filename,".png"), "\n")
}

cat("✅ Workflow complete. DEGs, KEGG results, and plots saved.\n")

library(ggplot2)

plot_volcano <- function(res, title="Volcano Plot") {
  res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "yes", "no")
  
  p <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c("no"="grey", "yes"="red")) +
    theme_minimal() +
    labs(title=title, x="log2 Fold Change", y="-log10(adj p-value)")
  
  print(p)   # display in R
  ggsave(paste0(gsub(" ", "_", title), ".png"), plot=p, width=8, height=6)
}

# Example usage
plot_volcano(res_SZ, "SZ vs CTL Volcano")
plot_volcano(res_AD, "AD vs CTL Volcano")

library(pheatmap)

plot_top_heatmap <- function(res, counts, top_n=20, title="Top Genes Heatmap") {
  # Get top genes by adjusted p-value
  top_genes <- rownames(head(res[order(res$padj),], top_n))
  
  mat <- counts[top_genes, , drop=FALSE]
  mat <- log2(mat + 1)
  
  p <- pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, 
                main=title, scale="row", fontsize_row=8)
  
  ggsave(paste0(gsub(" ", "_", title), ".png"), width=10, height=8)
}

# Example usage
plot_top_heatmap(res_SZ, counts1, top_n=20, title="SZ Top Genes")
plot_top_heatmap(res_AD, counts2, top_n=20, title="AD Top Genes")


