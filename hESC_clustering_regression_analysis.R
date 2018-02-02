#### Create t-SNE plots, differential expression heatmaps, and cluster enrichment heatmaps
library(perturbLM)
library(Seurat)

# file.handle <- "up-tf-all" ## Name of the sample to be analyzed
file.handle <- "up-tf-endo"
# file.handle <- "up-tf-multi"
# file.handle <- "up-tf-klf"
# file.handle <- "up-tf-myc"

## Load regression results
coefs.df <- read.table(paste(file.handle, "_regression_trimmed.pvals.tsv", sep = ""), sep = "\t", 
                       header = T, stringsAsFactors = F)
top.hits <- subset(coefs.df, FDR < 0.3 & abs(cf) > 0.025 & !grepl(":", Group))
all.genes <- unique(coefs.df$Gene)
print(sort(table(top.hits$Group), decreasing = T))

## Load genotype cluster counts and calculate cluster enrichment
df <- read.table(paste(file.handle, "_cluster_genotype_counts.tsv", sep = ""), sep = "\t")
chisq.p <- GenotypeClusterChisq(df, "mCherry")
chisq.fdr <- p.adjust(chisq.p, method = "BH")
df.pvals <- GenotypeClusterPvals(df)
df.fdr <- matrix(p.adjust(df.pvals, method = "BH"), nrow(df.pvals), ncol(df.pvals), dimnames = list(rownames(df.pvals), colnames(df.pvals)))

## Summarize differential expression results
genotype.counts <- sort(table(top.hits$Group), decreasing = T)
genotype.counts <- genotype.counts[names(genotype.counts) %in% rownames(df.pvals)]
genotype.avg_diff <- tapply(top.hits$cf, top.hits$Group, function(x) mean(abs(x)))[names(genotype.counts)]

genotype.sizes <- sapply(ReadGenotypes(paste(file.handle, "_pheno_dict.csv", sep = "")), length)
genotype.trimmed.sizes <- sapply(ReadGenotypes(paste(file.handle, "_pheno_dict_trimmed.csv", sep = "")), length)

genotype.summary <- data.frame(n_diff_genes = as.integer(genotype.counts), avg_eff_size = genotype.avg_diff)
genotype.summary <- genotype.summary[!grepl(":", rownames(genotype.summary)), ]
genotype.summary$min.cluster.fdr <- apply(df.fdr[rownames(genotype.summary), ], 1, function(x) x[which.min(x)])
genotype.summary$chisq.fdr <- chisq.fdr[rownames(genotype.summary)]
genotype.summary$n_cells <- genotype.sizes[rownames(genotype.summary)]
genotype.summary$n_cells_trimmed <- genotype.trimmed.sizes[rownames(genotype.summary)]

# write.table(genotype.summary, file = paste(file.handle, "_genotype_summary.tsv", sep = ""), sep = "\t")

## Determine significant genotypes
# genotype.summary <- read.table(paste(file.handle, "_genotype_summary.tsv", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
sig.genotypes <- union(rownames(subset(genotype.summary, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                       rownames(subset(genotype.summary, n_diff_genes >= 100)))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry")))


## tSNE plots (Figure 1c-e)
se.obj <- readRDS(paste(file.handle, "_clustering.Robj", sep = ""))
se.obj@meta.data$tree.ident <- sapply(se.obj@meta.data$tree.ident, function(x) paste("C", x, sep = ""))

# pdf("stem_tsne_plot.pdf", width = 4, height = 4)
PlotDims(se.obj@dr$tsne@cell.embeddings[,1], se.obj@dr$tsne@cell.embeddings[,2], cluster = se.obj@ident, 
         do.label = T, alpha = 0.25, label.size = 8, pt.size = 0.8)
# dev.off()

# pdf("stem_tsne_plot_nolabel.pdf", width = 4, height = 4)
PlotDims(se.obj@dr$tsne@cell.embeddings[,1], se.obj@dr$tsne@cell.embeddings[,2], cluster = se.obj@ident,
         do.label = F, alpha = 0.25, label.size = 8, pt.size = 0.8)
# dev.off()


## Cluster enrichment heatmap (Figure 1c-e)
sig.clusters <- apply(df.fdr[sig.genotypes,], 1, function(x) names(x[x < 1e-2])); names(sig.clusters) <- sig.genotypes;

heat.mat <- -log10(df.pvals[sig.genotypes, unique(unlist(sig.clusters, F, F))])
heat.mat[heat.mat > 75] <- 75

# pdf("cluster_enrichment_heatmap.pdf", width = 3.25, height = 2.75)
ggheat(heat.mat, clustering = "none", x.lab.size = 15, y.lab.size = 15) + theme(text = element_text(family = "sans"))
# dev.off()


## Differential expression heatmap (Figure S3a-c)
top.hits.plot <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025 & Group %in% sig.genotypes)
top.hits.plot <- top.hits.plot[order(top.hits.plot$cf, decreasing = T),]
top.hits.plot <- top.hits.plot[order(top.hits.plot$Group),]
# top.hits.plot <- Reduce(rbind, by(top.hits.plot, top.hits.plot$Group, head, n = 5))

heat.df <- coefs.df[coefs.df$Gene %in% top.hits.plot$Gene & coefs.df$Group %in% top.hits.plot$Group,]
heat.df$Score <- mapply(function(p,cf) -log10(p) * sign(cf), heat.df$p_val, heat.df$cf)
heat.mat <- UnflattenDataframe(heat.df, output.name = "Score", row.col = 'Gene', col.col = 'Group')
heat.mat <- heat.mat[unique(top.hits.plot$Gene),]

# pdf("diff_expr_heatmap.pdf", width = 5, height = 3)
ggheat(t(heat.mat), clustering = "none", labRow = T, labCol = F, y.lab.size = 15, x.lab.size = 15)
# dev.off()
