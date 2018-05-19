#### Create t-SNE plots, differential expression heatmaps, and cluster enrichment heatmaps
library(perturbLM)
library(swne)
library(pagoda2)
library(ggplot2)
library(RColorBrewer)

## Name of the sample to be analyzed
file.handle <- "up-tf-stem"
# file.handle <- "up-tf-endo"
# file.handle <- "up-tf-multi"
# file.handle <- "up-tf-klf"
# file.handle <- "up-tf-myc"

## Load regression results
coefs.df <- read.table(paste(file.handle, "regression.pvals.tsv", sep = "_"), sep = "\t", 
                       header = T, stringsAsFactors = F)
top.coefs.df <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025)
print(sort(table(top.coefs.df$Group), decreasing = T))

# ## Write top regression results
# write.table(top.coefs.df, file = paste(file.handle, "hits.pvals.tsv", sep = "_"), sep = "\t")

## Load clustering
r <- readRDS(paste(file.handle, "clustering_p2.Robj", sep = "_"))
clusters <- r$clusters$PCA$multilevel
levels(clusters) <- paste("C", levels(clusters), sep = "")
clusters.list <- UnflattenGroups(clusters)

## Load genotypes
# genotypes.list <- ReadGenotypes("up-tf-klf-myc_pheno_dict.csv")
genotypes.list <- ReadGenotypes(paste(file.handle, "pheno_dict.csv", sep = "_"))
genotypes.list <- lapply(genotypes.list, function(x) x[x %in% names(clusters)])
genotypes.sizes <- sapply(genotypes.list, length)

##  Calculate genotype enrichment
genotypes.cl.tbl <- GenotypeClusterCounts(genotypes.list, clusters.list)
genotypes.cl.pvals <- GenotypeClusterPvals(genotypes.cl.tbl)

## Summarize analysis results
genotype.counts <- sort(table(top.coefs.df$Group), decreasing = T)
genotype.counts <- genotype.counts[names(genotype.counts) %in% rownames(genotypes.cl.pvals)]
genotype.avg_cf <- tapply(top.coefs.df$cf, top.coefs.df$Group, function(x) mean(abs(x)))[names(genotype.counts)]

genotype.summary <- data.frame(n_diff_genes = as.integer(genotype.counts), avg_cf = genotype.avg_cf)
genotype.summary$min_cluster_pval <- apply(genotypes.cl.pvals[rownames(genotype.summary), ], 1, function(x) x[which.min(abs(x))])
genotype.summary$n_cells <- genotypes.sizes[rownames(genotype.summary)]

# write.table(genotype.summary, file = paste(file.handle, "summary.tsv", sep = "_"), sep = "\t")

## Determine significant genotypes
min.cl.pval <- 1e-12
min.diff.genes <- 50
sig.genotypes <- rownames(subset(genotype.summary, abs(min_cluster_pval) < min.cl.pval | n_diff_genes >= min.diff.genes))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry-int-ctrl"))); print(sig.genotypes);

## Cluster enrichment heatmap (Figure 1c-e)
max.lp <- 50

genotypes.cl.lp <- apply(genotypes.cl.pvals[sig.genotypes,], 1:2, function(x) ifelse(x > 0, -log10(x), log10(-1 * x)))
genotypes.cl.lp[genotypes.cl.lp > max.lp] <- max.lp
genotypes.cl.lp[genotypes.cl.lp < -1 * max.lp] <- -1 * max.lp

pdf(paste(file.handle, "cl_enrich_heatmap.pdf", sep = "_"), width = 5, height = 7.5)
ggHeat(genotypes.cl.lp, clustering = "none", x.lab.size = 14, y.lab.size = 14)
dev.off()


## tSNE plots (Figure 1c-e)
tsne.emb <- r$embeddings$PCA$tSNE

pdf(paste(file.handle, "tsne_plot.pdf", sep = "_"), width = 4, height = 4)
PlotDims(tsne.emb, sample.groups = clusters, pt.size = 0.75, alpha.plot = 0.5, label.size = 6, do.label = T,
         show.legend = F)
dev.off()

pdf(paste(file.handle, "tsne_plot_nolabel.pdf", sep = "_"), width = 4, height = 4)
PlotDims(tsne.emb, sample.groups = clusters, pt.size = 0.75, alpha.plot = 0.5, label.size = 6, do.label = F,
         show.legend = F)
dev.off()

pdf(paste(file.handle, "tsne_plot_batch.pdf", sep = "_"), width = 4, height = 4)
PlotDims(tsne.emb, sample.groups = r$batch, pt.size = 1, alpha.plot = 0.5, label.size = 8, do.label = T,
         show.legend = F)
dev.off()


# ## Differential expression heatmap (Figure S3a-c)
# top.hits.plot <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025 & Group %in% sig.genotypes)
# top.hits.plot <- top.hits.plot[order(top.hits.plot$cf, decreasing = T),]
# top.hits.plot <- top.hits.plot[order(top.hits.plot$Group),]
# # top.hits.plot <- Reduce(rbind, by(top.hits.plot, top.hits.plot$Group, head, n = 5))
# 
# heat.df <- coefs.df[coefs.df$Gene %in% top.hits.plot$Gene & coefs.df$Group %in% top.hits.plot$Group,]
# heat.df$Score <- mapply(function(p,cf) -log10(p) * sign(cf), heat.df$p_val, heat.df$cf)
# heat.mat <- UnflattenDataframe(heat.df, output.name = "Score", row.col = 'Gene', col.col = 'Group')
# heat.mat <- heat.mat[unique(top.hits.plot$Gene),]
# 
# # pdf("diff_expr_heatmap.pdf", width = 5, height = 3)
# ggheat(t(heat.mat), clustering = "none", labRow = T, labCol = F, y.lab.size = 15, x.lab.size = 15)
# # dev.off()
