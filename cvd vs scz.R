# =====================================================
# Full RNA-seq DE & KEGG workflow for AD, SZ, COVID
# =====================================================

suppressMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
})

# ----------------------------
# 1. Load counts
# ----------------------------
counts_AD    <- read.csv("C:/Users/Jessica/Downloads/GSE270454_RNAseq-combined-counts-matrix.csv.gz", row.names=1, check.names=FALSE)
counts_SZ    <- read.csv("C:/Users/Jessica/Downloads/GSE138082_Subfields_Counts.txt.gz", sep="\t", row.names=1, check.names=FALSE)
counts_COVID <- read.delim("C:/Users/Jessica/Downloads/GSE157103_genes.ec.tsv.gz", row.names=1, check.names=FALSE)

# Convert counts to integers
counts_AD    <- round(as.matrix(counts_AD))
counts_SZ    <- round(as.matrix(counts_SZ))
counts_COVID <- round(as.matrix(counts_COVID))

# ----------------------------
# 2. Create phenotype tables
# ----------------------------
make_pheno <- function(counts, pos_pattern, neg_pattern, dataset_name){
  ctl <- grep(neg_pattern, colnames(counts), value=TRUE)
  case <- grep(pos_pattern, colnames(counts), value=TRUE)
  if(length(ctl)==0 | length(case)==0){
    stop("No samples matched patterns for ", dataset_name)
  }
  pheno <- data.frame(
    sample = c(ctl, case),
    condition = factor(c(rep("CTL", length(ctl)), rep(pos_pattern, length(case))),
                       levels=c("CTL", pos_pattern)),
    dataset = dataset_name
  )
  rownames(pheno) <- pheno$sample
  counts <- counts[, c(ctl, case)]
  return(list(counts=counts, pheno=pheno))
}

# ----------------------------
# AD dataset
# ----------------------------
AD_list <- make_pheno(counts_AD, pos_pattern="AD", neg_pattern="MCI", "AD")
counts_AD <- AD_list$counts
pheno_AD  <- AD_list$pheno

# ----------------------------
# SZ dataset
# ----------------------------
SZ_list <- make_pheno(counts_SZ, pos_pattern="SZ", neg_pattern="CTL", "SZ")
counts_SZ <- SZ_list$counts
pheno_SZ  <- SZ_list$pheno

# ----------------------------
# COVID dataset (automatic selection)
# ----------------------------
# Here we assume first 20 samples are controls and the rest are COVID cases
ctl_COVID  <- colnames(counts_COVID)[1:20]  # adjust if you know which are controls
case_COVID <- setdiff(colnames(counts_COVID), ctl_COVID)

pheno_COVID <- data.frame(
  sample = c(ctl_COVID, case_COVID),
  condition = factor(c(rep("CTL", length(ctl_COVID)), rep("COVID", length(case_COVID))),
                     levels=c("CTL","COVID")),
  dataset="COVID"
)
rownames(pheno_COVID) <- pheno_COVID$sample
counts_COVID <- counts_COVID[, c(ctl_COVID, case_COVID)]

# ----------------------------
# 3. DESeq2
# ----------------------------
run_deseq <- function(counts, pheno){
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=pheno, design=~condition)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  return(list(dds=dds, res=res))
}

AD_DEG    <- run_deseq(counts_AD, pheno_AD)
SZ_DEG    <- run_deseq(counts_SZ, pheno_SZ)
COVID_DEG <- run_deseq(counts_COVID, pheno_COVID)

# ----------------------------
# 4. Save DEGs
# ----------------------------
write.csv(AD_DEG$res, "DEG_ADvsMCI.csv", row.names=FALSE)
write.csv(SZ_DEG$res, "DEG_SZvsCTL.csv", row.names=FALSE)
write.csv(COVID_DEG$res, "DEG_COVIDvsCTL.csv", row.names=FALSE)

# ----------------------------
# 5. KEGG enrichment
# ----------------------------
run_kegg <- function(res, name){
  sig_genes <- res %>% filter(!is.na(padj) & padj < 0.05)
  if(nrow(sig_genes)==0){
    cat("No significant genes for", name, "- skipping KEGG.\n")
    return(NULL)
  }
  entrez <- mapIds(org.Hs.eg.db,
                   keys=rownames(sig_genes),
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first") %>% na.omit()
  if(length(entrez)==0){
    cat("No valid ENTREZ IDs for", name, "- skipping KEGG.\n")
    return(NULL)
  }
  kegg <- enrichKEGG(gene=entrez, organism="hsa", pvalueCutoff=0.05)
  if(!is.null(kegg)){
    write.csv(kegg@result, paste0("KEGG_", name,".csv"), row.names=FALSE)
    cat(name,"KEGG pathways:", nrow(kegg@result), "\n")
  }
  return(kegg)
}

kegg_AD    <- run_kegg(AD_DEG$res, "ADvsMCI")
kegg_SZ    <- run_kegg(SZ_DEG$res, "SZvsCTL")
kegg_COVID <- run_kegg(COVID_DEG$res, "COVIDvsCTL")

# ----------------------------
# 6. Top gene plots
# ----------------------------
plot_top_genes <- function(res, counts_matrix, filename, n=10){
  sig <- res %>% filter(!is.na(padj)) %>% arrange(padj)
  top_genes <- head(sig$gene, n)
  genes <- intersect(top_genes, rownames(counts_matrix))
  if(length(genes)==0) return(cat("No matching genes for", filename, "\n"))
  
  mat <- counts_matrix[genes, , drop=FALSE]
  mat <- log2(mat + 1)
  mat_long <- melt(mat)
  colnames(mat_long) <- c("Gene","Sample","log2count")
  
  p <- ggplot(mat_long, aes(x=Sample, y=log2count, fill=Gene)) +
    geom_bar(stat="identity", position="dodge") +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    labs(title=paste("Top genes:", filename), y="log2(count+1)", x="Samples")
  
  print(p)
  ggsave(paste0(filename,".png"), plot=p, width=10, height=6)
  cat("Plot saved:", paste0(filename,".png"), "\n")
}

plot_top_genes(AD_DEG$res, counts_AD, "Top_AD_genes")
plot_top_genes(SZ_DEG$res, counts_SZ, "Top_SZ_genes")
plot_top_genes(COVID_DEG$res, counts_COVID, "Top_COVID_genes")

cat("✅ Workflow complete. DEGs, KEGG, and plots saved.\n")

# ----------------------------
# KEGG pathway visualization
# ----------------------------
suppressMessages({
  if(!requireNamespace("pathview", quietly = TRUE)) {
    BiocManager::install("pathview")
  }
  library(pathview)
})

# Function to plot top KEGG pathways
plot_kegg_pathways <- function(kegg_obj, res, name, top_n=3) {
  if(is.null(kegg_obj)) {
    cat("No KEGG pathways for", name, "- skipping visualization.\n")
    return(NULL)
  }
  
  # Select top pathways by p-value
  top_pathways <- head(kegg_obj@result$ID, top_n)
  
  # Prepare gene fold changes
  fc <- res$log2FoldChange
  names(fc) <- rownames(res)
  
  for(p in top_pathways){
    cat("Plotting KEGG pathway:", p, "for", name, "\n")
    pathview(gene.data = fc,
             pathway.id = p,
             species = "hsa",
             out.suffix = paste0(name,"_",p))
  }
}

# ----------------------------
# Plot KEGG maps
# ----------------------------
plot_kegg_pathways(kegg_AD, AD_DEG$res, "ADvsMCI")
plot_kegg_pathways(kegg_SZ, SZ_DEG$res, "SZvsCTL")
plot_kegg_pathways(kegg_COVID, COVID_DEG$res, "COVIDvsCTL")

cat("✅ KEGG pathway maps saved.\n")

