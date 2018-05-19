library(methods)
library(perturbLM)
library(ggplot2)

## Name of the sample to be analyzed
options <- commandArgs(trailingOnly = T)
file.handle <- options[[1]]
n.cores <- as.integer(options[[2]])

min.cl.pval <- 1e-12
min.diff.genes <- 50

n.rand <- 10000
min.set.genes <- 10
max.set.genes <- 500

## Load regression results
coefs.df <- read.table(paste(file.handle, "_regression.pvals.tsv", sep = ""), sep = "\t", 
                       header = T, stringsAsFactors = F)
top.coefs.df <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025)
n.diff.genes <- table(top.coefs.df$Group)
print(sort(n.diff.genes, decreasing = T))

genotype.summary <- read.table(paste(file.handle, "_summary.tsv", sep = ""), sep = "\t", header = T, row.names = 1)
genotypes.use <- rownames(subset(genotype.summary, abs(min_cluster_pval) < min.cl.pval | n_diff_genes >= min.diff.genes))
genotypes.use <- sort(unique(c(genotypes.use, "mCherry-int-ctrl"))); print(genotypes.use);

## Load genesets
genesets.file <- "hallmark_go_bp_genesets.gmt"
genesets <- LoadGenesets(genesets.file)
genesets <- FilterGenesets(genesets, unique(coefs.df$Gene), min.size = min.set.genes, max.size = max.set.genes)

## Run GSEA on TFs
coefs.list <- lapply(genotypes.use, function(g) {
  df <- subset(coefs.df, Group == g)
  cfs <- df$cf; names(cfs) <- df$Gene;
  cfs
})
names(coefs.list) <- genotypes.use

gsea.df <- MultipleGSEAEnrich(coefs.list, genesets, n.rand = n.rand, n.cores = n.cores)
write.table(gsea.df, paste(file.handle, "gsea.tsv", sep = "_"), sep = "\t")

## Run Fisher Enrichment on TFs
top.markers.list  <- lapply(genotypes.use, function(g) {
  df <- subset(top.coefs.df, Group == g)
  return(as.character(df$Gene))
})
names(top.markers.list) <- genotypes.use

fisher.df <- MultipleFisherEnrich(top.markers.list, genesets)
write.table(fisher.df, paste(file.handle, "geneset_fisher.tsv", sep = "_"), sep = "\t")
