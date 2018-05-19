#### Correlate stem cell media replicates ####
library(perturbLM)
library(swne)

ind.files <- c("Individual_Batch_Analysis/up-tf1_regression.pvals.tsv",
               "Individual_Batch_Analysis/up-tf2_regression.pvals.tsv",
               "Individual_Batch_Analysis/up-tf3_regression.pvals.tsv",
               "Individual_Batch_Analysis/up-tf4_regression.pvals.tsv",
               "Individual_Batch_Analysis/up-tf5_regression.pvals.tsv")
names(ind.files) <- paste("stem", seq(1,length(ind.files)), sep = "-")

cor.matrix <- matrix(0, length(ind.files), length(ind.files))
rownames(cor.matrix) <- colnames(cor.matrix) <- names(ind.files)

ind.coefs <- lapply(ind.files, function(fi) {
  coefs.df.batch <- read.table(fi, sep = "\t", header = T, stringsAsFactors = F)
  coefs.df.batch$Hit <- coefs.df.batch$FDR < 0.5
  rownames(coefs.df.batch) <- interaction(coefs.df.batch$Gene, coefs.df.batch$Group, sep = "_")
  coefs.df.batch
})
names(ind.coefs) <- names(ind.files)

file.comb <- combn(names(ind.files), m = 2)
for (i in 1:ncol(file.comb)) {
  fi_1 <- file.comb[1, i]
  fi_2 <- file.comb[2, i]
  
  df.1 <- ind.coefs[[fi_1]]
  df.2 <- ind.coefs[[fi_2]]
  
  common.idx <- intersect(rownames(df.1), rownames(df.2))
  sig.idx <- union(rownames(subset(df.1, Hit)), rownames(subset(df.2, Hit)))
  hits.idx <- intersect(common.idx, sig.idx)
  
  cfs_1 <- df.1$cf; names(cfs_1) <- rownames(df.1);
  cfs_2 <- df.2$cf; names(cfs_2) <- rownames(df.2);
  
  cfs.r <- cor(cfs_1[hits.idx], cfs_2[hits.idx])
  cor.matrix[fi_1, fi_2] <- cfs.r
  cor.matrix[fi_2, fi_1] <- cfs.r
  
  f.name <- paste("up-tf-stem_corr_", fi_1, "_vs_", fi_2, ".pdf", sep = "")
  pdf(f.name, width = 3, height = 3)
  print(PlotHexBin(cfs_1[hits.idx], cfs_2[hits.idx], n.bins = 25))
  # print(PlotCorrelation(cfs_1[hits.idx], cfs_2[hits.idx], show.corr = T, pt.size = 1.25, alpha = 0.5))
  dev.off()
}

diag(cor.matrix) <- NA
cor.matrix <- round(cor.matrix, 3)

col.palette <- colorRampPalette(c("white", "tomato"))

pdf("up-tf-stem_batch_corr_matrix.pdf", width = 4.5, height = 4.5)
gplots::heatmap.2(cor.matrix, Rowv = F, Colv = F, dendrogram = "none", cellnot = cor.matrix,
                  notecex = 1, col = col.palette, notecol = "black", trace = "none",
                  breaks = 10, cexRow = 1.25, cexCol = 1.25, density.info = "none",
                  margins = c(7,7), key = F)
dev.off()

#### Correlate multilineage media batches ####
library(perturbLM)
library(swne)

multi.coefs.df.1 <- read.table("Individual_Batch_Analysis/up-tf-multi-batch1_regression.pvals.tsv", 
                               sep = "\t", stringsAsFactors = F)
multi.cfs.1 <- multi.coefs.df.1$cf
rownames(multi.coefs.df.1) <- names(multi.cfs.1) <- interaction(multi.coefs.df.1$Gene, multi.coefs.df.1$Group, sep = "_")

multi.coefs.df.2 <- read.table("Individual_Batch_Analysis/up-tf-multi-batch2_regression.pvals.tsv", 
                         sep = "\t", stringsAsFactors = F)
multi.cfs.2 <- multi.coefs.df.2$cf
rownames(multi.coefs.df.2) <- names(multi.cfs.2) <- interaction(multi.coefs.df.2$Gene, multi.coefs.df.2$Group, sep = "_")

multi.common.idx <- intersect(names(multi.cfs.1), names(multi.cfs.2))
multi.sig.idx <- union(rownames(subset(multi.coefs.df.1, FDR < 0.5)), rownames(subset(multi.coefs.df.2, FDR < 0.5)))
multi.hits.idx <- intersect(multi.common.idx, multi.sig.idx)

pdf("up-tf-multi_corr_multi-1_vs_multi-2.pdf", width = 3, height = 3)
PlotHexBin(multi.cfs.1[multi.hits.idx], multi.cfs.2[multi.hits.idx], n.bins = 25)
dev.off()


#### Correlate multilineage media batches ####
library(perturbLM)
library(swne)

endo.coefs.df.1 <- read.table("Individual_Batch_Analysis/up-tf-endo-batch1_regression.pvals.tsv", 
                               sep = "\t", stringsAsFactors = F)
endo.cfs.1 <- endo.coefs.df.1$cf
rownames(endo.coefs.df.1) <- names(endo.cfs.1) <- interaction(endo.coefs.df.1$Gene, endo.coefs.df.1$Group, sep = "_")

endo.coefs.df.2 <- read.table("Individual_Batch_Analysis/up-tf-endo-batch2_regression.pvals.tsv", 
                               sep = "\t", stringsAsFactors = F)
endo.cfs.2 <- endo.coefs.df.2$cf
rownames(endo.coefs.df.2) <- names(endo.cfs.2) <- interaction(endo.coefs.df.2$Gene, endo.coefs.df.2$Group, sep = "_")

endo.common.idx <- intersect(names(endo.cfs.1), names(endo.cfs.2))
endo.sig.idx <- union(rownames(subset(endo.coefs.df.1, FDR < 0.5)), rownames(subset(endo.coefs.df.2, FDR < 0.5)))
endo.hits.idx <- intersect(endo.common.idx, endo.sig.idx)

pdf("up-tf-endo_corr_endo-1_vs_endo-2.pdf", width = 3, height = 3)
PlotHexBin(endo.cfs.1[endo.hits.idx], endo.cfs.2[endo.hits.idx], n.bins = 25)
dev.off()