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
# file.handle <- "up-tf-neuron-nohygro"

## Load regression results
coefs.df <- read.table(paste(file.handle, "regression.pvals.tsv", sep = "_"), sep = "\t", 
                       header = T, stringsAsFactors = F)
top.coefs.df <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025)
print(sort(table(top.coefs.df$Group), decreasing = T))

## Load clustering
r <- readRDS(paste(file.handle, "clustering.Robj", sep = "_"))
clusters <- r@ident; names(clusters) <- r@cell.names;
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

genotypes.min.p <- apply(genotypes.cl.pvals, 1, function(x) min(abs(x)) * length(x))
metap::sumlog(genotypes.min.p)

## Summarize analysis results
genotype.counts <- sort(table(top.coefs.df$Group), decreasing = T)
genotype.counts <- genotype.counts[names(genotype.counts) %in% rownames(genotypes.cl.pvals)]
genotype.avg_cf <- tapply(top.coefs.df$cf, top.coefs.df$Group, function(x) mean(abs(x)))[names(genotype.counts)]

genotype.summary <- data.frame(n_diff_genes = as.integer(genotype.counts), avg_cf = genotype.avg_cf)
genotype.summary$min_cluster_pval <- apply(genotypes.cl.pvals[rownames(genotype.summary), ], 1, function(x) min(abs(x)))
genotype.summary$min_cluster <- apply(genotypes.cl.pvals[rownames(genotype.summary), ], 1, function(x) names(x[which.min(abs(x))]))
genotype.summary$n_cells <- genotypes.sizes[rownames(genotype.summary)]

# write.table(genotype.summary, file = paste(file.handle, "summary.tsv", sep = "_"), sep = "\t")

## Determine significant genotypes
min.cl.pval <- 1e-12
min.diff.genes <- 40
ctrl.cl <- genotype.summary["mCherry", "min_cluster"]
sig.genotypes <- rownames(subset(genotype.summary, (abs(min_cluster_pval) < min.cl.pval & min_cluster != ctrl.cl) 
                                 | n_diff_genes > min.diff.genes))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry-int-ctrl")), decreasing = T); print(sig.genotypes);


## Differential expression barplot
n.diff.genes <- genotype.summary[sig.genotypes, "n_diff_genes"]
names(n.diff.genes) <- sig.genotypes
gg.bar <- ggBarplot(n.diff.genes, fill.color = "purple")

# pdf(paste(file.handle, "diff_genes_barplot.pdf", sep = "_"), width = 3.5, height = 1.5)
# print(gg.bar)
# dev.off()

## Cluster enrichment heatmap (Figure 1c-e)
max.lp <- 50

genotypes.cl.lp <- apply(genotypes.cl.pvals[sig.genotypes,], 1:2, function(x) ifelse(x > 0, -log10(x), log10(-1 * x)))
genotypes.cl.lp[genotypes.cl.lp > max.lp] <- max.lp
genotypes.cl.lp[genotypes.cl.lp < -1 * max.lp] <- -1 * max.lp

# pdf(paste(file.handle, "cl_enrich_heatmap.pdf", sep = "_"), width = 4.5, height = 5)
# ggHeat(genotypes.cl.lp, clustering = "none", x.lab.size = 14, y.lab.size = 14)
# dev.off()

gg.heat <- ggHeat(t(genotypes.cl.lp), clustering = "none", x.lab.size = 14, y.lab.size = 14)
gg.bar <- gg.bar + theme(axis.text.x = element_blank())
gg.heat <- gg.heat + theme(legend.position = "none")

library(gtable)
gg.1 <- ggplotGrob(gg.heat)
gg.2 <- ggplotGrob(gg.bar)

gg.aligned <- rbind(gg.2, gg.1, size = "first")
gg.aligned$widths <- unit.pmax(gg.1$widths, gg.2$widths)

pdf(paste(file.handle, "stacked_diff_genes_cl_enrich_nolabels.pdf", sep = "_"), width = 5.5, height = 6)
grid.newpage()
grid.draw(gg.aligned)
dev.off()



## tSNE plots (Figure 1c-e)
plot.seed <- 32907098
tsne.emb <- GetCellEmbeddings(r, reduction.type = "tsne")

pdf(paste(file.handle, "tsne_plot.pdf", sep = "_"), width = 4, height = 4)
PlotDims(tsne.emb, sample.groups = clusters, pt.size = 0.75, alpha.plot = 0.5, label.size = 6, do.label = T,
         show.legend = F, seed = plot.seed)
dev.off()

pdf(paste(file.handle, "tsne_plot_nolabel.pdf", sep = "_"), width = 4, height = 4)
PlotDims(tsne.emb, sample.groups = clusters, pt.size = 0.75, alpha.plot = 0.5, label.size = 0, do.label = F,
         show.legend = F, seed = plot.seed)
dev.off()

pdf(paste(file.handle, "tsne_plot_batch.pdf", sep = "_"), width = 4, height = 4)
batch <- factor(r@meta.data$batch); names(batch) <- r@cell.names; batch <- batch[rownames(tsne.emb)];
PlotDims(tsne.emb, sample.groups = batch, pt.size = 1, alpha.plot = 0.5, label.size = 8, do.label = T,
         show.legend = F, seed = plot.seed)
dev.off()


## tSNE plot overlay
plot.seed <- 32907098
tsne.emb <- GetCellEmbeddings(r, reduction.type = "tsne")

TF <- "NEUROD1"

tf.groups <- as.character(clusters); names(tf.groups) <- names(clusters);
tf.groups[names(tf.groups) %in% genotypes.list[[TF]]] <- TF
tf.groups[!names(tf.groups) %in% genotypes.list[[TF]]] <- ""
tf.groups <- factor(tf.groups)

pdf(paste0(file.handle, "_tsne_plot_", TF, ".pdf"), width = 3, height = 3)
PlotDims(tsne.emb, sample.groups = tf.groups, pt.size = 0.35, alpha.plot = 0.4, label.size = 0, do.label = T,
         show.legend = F, seed = plot.seed) + scale_color_manual(values = c("grey", "red")) + ggtitle(TF)
dev.off()