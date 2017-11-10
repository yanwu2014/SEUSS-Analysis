#### Generate differential expression heatmaps, tSNE plots, cluster enrichment plots ####

source("hESC_functions.R")

file.handle <- "up-tf-all" ## Name of the sample to be analyzed
# file.handle <- "up-tf-endo"
# file.handle <- "up-tf-multi"
# file.handle <- "up-tf-klf"
# file.handle <- "up-tf-myc"

## Load regression results
coefs.df <- read.table(paste(file.handle, "_regression_trimmed.pvals.tsv", sep = ""), sep = "\t", header = T, stringsAsFactors = F)
# coefs.df <- read.table(paste(file.handle, "_regression.pvals.tsv", sep = ""), sep = "\t", header = T, stringsAsFactors = F)
top.hits <- subset(coefs.df, FDR < 0.3 & abs(cf) > 0.025 & !grepl(":", Group))
all.genes <- unique(coefs.df$Gene)
print(sort(table(top.hits$Group), decreasing = T))

## Load genotype cluster counts and calculate cluster enrichment
df <- read.table(paste(file.handle, "_cluster_genotype_counts.tsv", sep = ""), sep = "\t")
chisq.p <- genotype.cluster.chisq(df, "mCherry")
chisq.fdr <- p.adjust(chisq.p, method = "BH")
df.pvals <- genotype.cluster.pvals(df)
df.fdr <- matrix(p.adjust(df.pvals, method = "BH"), nrow(df.pvals), ncol(df.pvals), dimnames = list(rownames(df.pvals), colnames(df.pvals)))

se.obj <- readRDS(paste(file.handle, "_clustering.Robj", sep = ""))

## Summarize differential expression results
genotype.counts <- sort(table(top.hits$Group), decreasing = T)
genotype.counts <- genotype.counts[names(genotype.counts) %in% rownames(df.pvals)]
genotype.avg_diff <- tapply(top.hits$cf, top.hits$Group, function(x) mean(abs(x)))[names(genotype.counts)]

genotype.sizes <- sapply(read.genotypes(paste(file.handle, "_pheno_dict.csv", sep = "")), length)
genotype.trimmed.sizes <- sapply(read.genotypes(paste(file.handle, "_pheno_dict_trimmed.csv", sep = "")), length)

genotype.summary <- data.frame(n_diff_genes = as.integer(genotype.counts), avg_eff_size = genotype.avg_diff)
genotype.summary <- genotype.summary[!grepl(":", rownames(genotype.summary)), ]
genotype.summary$min.cluster.fdr <- apply(df.fdr[rownames(genotype.summary), ], 1, function(x) x[which.min(x)])
genotype.summary$chisq.fdr <- chisq.fdr[rownames(genotype.summary)]
genotype.summary$n_cells <- genotype.sizes[rownames(genotype.summary)]
genotype.summary$n_cells_trimmed <- genotype.trimmed.sizes[rownames(genotype.summary)]

write.table(genotype.summary, file = paste(file.handle, "_genotype_summary.tsv", sep = ""), sep = "\t")

## Determine significant genotypes
# genotype.summary <- read.table(paste(file.handle, "_genotype_summary.tsv", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
sig.genotypes <- union(rownames(subset(genotype.summary, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                       rownames(subset(genotype.summary, n_diff_genes >= 100)))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry")))


## tSNE plots
library(Seurat)
se.obj <- readRDS("up-tf-all_clustering.Robj")
se.obj@meta.data$tree.ident <- sapply(se.obj@meta.data$tree.ident, function(x) paste("C", x, sep = ""))

pdf("multi_tsne_plot.pdf", width = 4, height = 4)
tsne_plot(se.obj@dr$tsne@cell.embeddings[,1], se.obj@dr$tsne@cell.embeddings[,2], color.plot = se.obj@ident, do.label = T,
          alpha = 0.25, label.size = 8, pt.size = 0.8)
dev.off()

pdf("multi_tsne_plot_nolabel.pdf", width = 4, height = 4)
tsne_plot(se.obj@dr$tsne@cell.embeddings[,1], se.obj@dr$tsne@cell.embeddings[,2], color.plot = se.obj@ident, do.label = F,
          alpha = 0.25, label.size = 8, pt.size = 0.8)
dev.off()


## Cluster enrichment heatmap
sig.clusters <- apply(df.fdr[sig.genotypes,], 1, function(x) names(x[x < 1e-2])); names(sig.clusters) <- sig.genotypes;

heat.mat <- -log10(df.pvals[sig.genotypes, unique(unlist(sig.clusters, F, F))])
heat.mat[heat.mat > 75] <- 75

pdf("cluster_enrichment_heatmap.pdf", width = 3.25, height = 2.75)
ggheat(heat.mat, clustering = "none", x.lab.size = 15, y.lab.size = 15) + theme(text = element_text(family = "sans"))
dev.off()


## Trim genotypes
max.FDR <- 1e-6

genotypes.input.file <- paste(file.handle, "pheno_dict.csv", sep = "_")
trimmed.genotypes.output.file <- paste(file.handle, "pheno_dict_trimmed.csv", sep = "_")

genotypes.list <- read.genotypes(genotypes.input.file)
sig.cl.enrich <- intersect(names(chisq.fdr[chisq.fdr < max.FDR]), rownames(df.fdr[apply(df.fdr, 1, function(x) min(x) < max.FDR), ]))
genotype.hit.clusters <- apply(df.fdr[sig.cl.enrich,], 1, function(x) names(x[x < 1e-2]))

clusters <- sapply(as.character(se.obj@meta.data$tree.ident), function(i) paste("C", i, sep = ""))
names(clusters) <- se.obj@cell.names
clusters.list <- unflatten.cell.genotypes(clusters)

trim.genotypes.list <- lapply(names(genotypes.list), function(g) {
  g.cells <- genotypes.list[[g]]
  if (g %in% sig.cl.enrich) {
    cl <- genotype.hit.clusters[[g]]
    cl.cells <- unlist(clusters.list[cl], F, F)
    return(intersect(g.cells, cl.cells))
  }
  else {
    return(g.cells)
  }
})
names(trim.genotypes.list) <- names(genotypes.list)

print(sapply(genotypes.list[sig.cl.enrich], length))
print(sapply(trim.genotypes.list[sig.cl.enrich], length))
write.genotypes(trim.genotypes.list, trimmed.genotypes.output.file)


## Differential expression heatmap
top.hits.plot <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025 & Group %in% sig.genotypes)
top.hits.plot <- top.hits.plot[order(top.hits.plot$cf, decreasing = T),]
top.hits.plot <- top.hits.plot[order(top.hits.plot$Group),]
# top.hits.plot <- Reduce(rbind, by(top.hits.plot, top.hits.plot$Group, head, n = 5))

heat.df <- coefs.df[coefs.df$Gene %in% top.hits.plot$Gene & coefs.df$Group %in% top.hits.plot$Group,]
heat.df$Score <- mapply(function(p,cf) -log10(p) * sign(cf), heat.df$p_val, heat.df$cf)
heat.mat <- unflatten.dataframe(heat.df, output.name = "Score", row.col = 'Gene', col.col = 'Group')
heat.mat <- heat.mat[unique(top.hits.plot$Gene),]

pdf("diff_expr_heatmap.pdf", width = 5, height = 3)
ggheat(t(heat.mat), clustering = "none", labRow = T, labCol = F, y.lab.size = 15, x.lab.size = 15)
dev.off()

#### Orthologous mouse TF enrichment in hPSC media ####

library(liger)
source("hESC_functions.R")

## Create mouse TF targets genesets
mouse.targets.df <- read.table("mouse_tf_targets.txt", sep = "\t", header = T, stringsAsFactors = F)
genotypes.keep <- c("Ascl1", "Cdx2", "Klf4", "Otx2", "Myc", "Myod1", "Otx2")

mouse.targets.df <- subset(mouse.targets.df, TF %in% genotypes.keep & log.fc > 0.2 & FDR < 0.05)
mouse.targets.df <- mouse.targets.df[order(mouse.targets.df$FDR),]
mouse.targets.df <- Reduce(rbind, by(mouse.targets.df, mouse.targets.df$TF, head, n = 600))
print(table(mouse.targets.df$TF))

mouse.genesets <- lapply(unique(mouse.targets.df$TF), function(tf) subset(mouse.targets.df, TF == tf)$Genes)
names(mouse.genesets) <- sapply(unique(mouse.targets.df$TF), toupper)

human.genesets <- lapply(mouse.genesets, function(genes) unique(convertMouseGeneList(genes)))
human.genesets <- lapply(human.genesets, function(genes) genes[genes %in% unique(coefs.df$Gene)])
print(sapply(human.genesets, length))

## Trim mouse TF targets genesets
min.corr <- 0.05
max.genesets <- 2

se.obj <- readRDS("up-tf-all_clustering.Robj")
human.genesets <- lapply(human.genesets, function(genes) genes[genes %in% rownames(se.obj@data)])
sapply(human.genesets, length)

genesets_pca <- lapply(human.genesets, function(genes) {
  scores <- prcomp(t(as.matrix(se.obj@data[genes,])), rank = 1, center = T, scale = T)$x[,"PC1"]
  means <- colMeans(as.matrix(se.obj@data[genes,]))
  if (cor(scores, means) < 0) { scores <- scores * -1 }
  scores
})

genesets_corrs <- lapply(names(human.genesets), function(g) {
  pca_scores <- genesets_pca[[g]]
  genes <- human.genesets[[g]]
  apply(as.matrix(se.obj@data[genes,]), 1, function(x) cor(x, pca_scores))
})
names(genesets_corrs) <- names(human.genesets)

human.genesets.trimmed <- lapply(genesets_corrs, function(corr) names(corr[corr > min.corr]))

genesets.genes <- unique(unlist(human.genesets.trimmed, F, F))
genesets.genes.counts <- sapply(genesets.genes, function(gene) sum(sapply(human.genesets.trimmed, function(genes) gene %in% genes)))
genesets.genes.counts <- genesets.genes.counts[genesets.genes.counts <= max.genesets]
human.genesets.trimmed <- lapply(human.genesets.trimmed, function(genes) genes[genes %in% names(genesets.genes.counts)])

sapply(human.genesets, length)
sapply(human.genesets.trimmed, length)
write.genesets(human.genesets.trimmed, file.name = "hESC_bulk_mouse_ortholog_tf_genesets.gmt")


## GSEA enrichment for mouse ortholog genesets
coefs.df <- read.table("up-tf-all_regression_trimmed.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
genotype.summary <- read.table("up-tf-all_genotype_summary.tsv", header = T, sep = "\t", stringsAsFactors = F)
sig.genotypes <- union(rownames(subset(genotype.summary, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                       rownames(subset(genotype.summary, n_diff_genes >= 100)))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry")))
all.genes <- unique(coefs.df$Gene)
n.cores <- 16

genesets <- load.genesets("hESC_bulk_mouse_ortholog_tf_genesets.gmt")
genesets <- filter.genesets(genesets, all.genes, min.size = 10, max.size = 2000)
sapply(genesets, length)

scores.list <- lapply(sig.genotypes, function(g) {
  g.coefs.df <- subset(coefs.df, Group == g)
  g.coefs <- g.coefs.df$cf; names(g.coefs) <- g.coefs.df$Gene;
  g.coefs
}); names(scores.list) <- sig.genotypes;

genesets.enrich <- multiple.gsea.enrich(scores.list, genesets, n.rand = 10000, n.cores = n.cores, power = 1)

top.genesets.enrich <- subset(genesets.enrich, FDR < 5e-2)
heat.df <- subset(genesets.enrich, Group %in% sig.genotypes & genesets %in% top.genesets.enrich$genesets)
heat.df$Score <- mapply(function(p, cf) -log(p) * sign(cf), heat.df$p.val, heat.df$sscore)
heat.mat <- unflatten.dataframe(heat.df, output.name = "Score", row.col = 'genesets', col.col = 'Group')
heat.mat[is.na(heat.mat)] <- 0

pdf("bulk_ortholog_gsea_enrichment.pdf", width = 6, height = 2.5)
ggheat(heat.mat, clustering = "col", x.lab.size = 12, y.lab.size = 12)
dev.off()

write.table(genesets.enrich, file = "up-tf-all_TF_gene_module_enrichment.tsv", sep = "\t")


#### EMT transition analysis in hPSC media ####

## PCA on Hallmark EMT genes
library(Seurat)

hallmark.genesets <- load.genesets("hallmark_genesets.gmt")
emt.genes <- hallmark.genesets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

genotypes.trimmed <- flatten.genotype.list(read.genotypes(paste(file.handle, "_pheno_dict_trimmed.csv", sep = "")))
cells.use <- names(genotypes.trimmed[genotypes.trimmed %in% c("SNAI2", "KLF4")])

se.obj <- readRDS("up-tf-all_clustering.Robj")
se.obj <- SubsetData(se.obj, cells.use = cells.use)
se.obj <- ScaleData(se.obj, genes.use = emt.genes, model.use = "negbinom")
se.obj <- RunPCA(se.obj, pc.genes = emt.genes, pcs.compute = 10, do.print = F)
se.obj@ident <- factor(genotypes.trimmed[se.obj@cell.names])

se.obj.pc.scores <- PCAEmbed(se.obj, dims.use = 1:2)

pdf("EMT_transition_KLF4_SNAI2.pdf", width = 5, height = 5)
dim_plot(se.obj.pc.scores[,1], se.obj.pc.scores[,2], color.plot = factor(se.obj@ident), alpha = 1,
         label.size = 8, show.legend = F)
dev.off()

se.obj.pc.loadings <- sort(PCALoad(se.obj, dims.use = 1)[,1], decreasing = T)
print(head(se.obj.pc.loadings, n = 40))
print(tail(se.obj.pc.loadings, n = 40))

## Barplot of key markers
coefs.df <- read.table(paste(file.handle, "_regression_trimmed.pvals.tsv", sep = ""), sep = "\t", header = T, stringsAsFactors = F)
genotype.summary <- read.table(paste(file.handle, "_genotype_summary.tsv", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
sig.genotypes <- union(rownames(subset(genotype.summary, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                       rownames(subset(genotype.summary, n_diff_genes >= 100)))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry")))

coefs.matrix.stem <- unflatten.dataframe(coefs.df, "cf", row.col = "Gene", col.col = "Group")
coefs.matrix.stem <- coefs.matrix.stem[, sig.genotypes]

emt.genes.plot <- c("SPP1", "EPCAM", "LAMC1", "THY1", "VIM", "TPM2")

coefs.plot <- c(coefs.matrix.stem[emt.genes.plot, "KLF4"], coefs.matrix.stem[emt.genes.plot, "SNAI2"])
tfs.factor <- c(rep("KLF4", length(emt.genes.plot)), rep("SNAI2", length(emt.genes.plot)))
genes.factor <- rep(emt.genes.plot, 2)

barplot.df <- data.frame(coefs = coefs.plot, TF = tfs.factor, Gene = genes.factor)

pdf("EMT_markers_barplot.pdf", width = 6, height = 4.5)
ggplot(data = barplot.df, aes(x = Gene, y = coefs)) + geom_bar(aes(fill = TF), position = "dodge", stat = "identity") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 20), legend.text = element_text(size = 16))
dev.off()


#### Summarize media condition differential expression in stacked barplot ####

library(ggplot2)

## Load stem cell media results
genotype.summary <- read.table("up-tf-all_genotype_summary.tsv", header = T, sep = "\t", stringsAsFactors = F)
sig.genotypes <- union(rownames(subset(genotype.summary, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                       rownames(subset(genotype.summary, n_diff_genes >= 100)))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry")))

## Load endothelial results
genotype.summary.endo <- read.table("up-tf-endo_genotype_summary.tsv", sep = "\t")
sig.genotypes.endo <- union(rownames(subset(genotype.summary.endo, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                            rownames(subset(genotype.summary.endo, n_diff_genes >= 100)))
sig.genotypes.endo <- sort(unique(c(sig.genotypes.endo, "mCherry")))

## Load multilineage media results
genotype.summary.multi <- read.table("up-tf-multi_genotype_summary.tsv", sep = "\t")
sig.genotypes.multi <- union(rownames(subset(genotype.summary.multi, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                             rownames(subset(genotype.summary.multi, n_diff_genes >= 100)))
sig.genotypes.multi <- sort(unique(c(sig.genotypes.multi, "mCherry")))


## Stacked barplot summary
all.sig.genotypes <- unique(c(sig.genotypes, sig.genotypes.multi, sig.genotypes.endo))
n.sig <- length(all.sig.genotypes)

all.diff.genes <- c(genotype.summary[all.sig.genotypes, "n_diff_genes"],
                    genotype.summary.endo[all.sig.genotypes, "n_diff_genes"],
                    genotype.summary.multi[all.sig.genotypes, "n_diff_genes"])
all.diff.genes[is.na(all.diff.genes)] <- 0
media.types <- c(rep("Stem cell", n.sig), rep("Endothelial", n.sig), rep("Multilineage", n.sig))
barplot.df <- data.frame(Genotype = rep(all.sig.genotypes, 3), Diff_expr_genes = all.diff.genes, Media = media.types)

pdf("media_condition_summary_barplot.pdf", width = 6, height = 6)
ggplot(barplot.df, aes(x = Genotype, y = Diff_expr_genes, fill = Media)) + geom_bar(stat = "identity") + theme_classic() +
  theme(axis.title.x = element_blank()) + coord_flip()
dev.off()


#### Correlate screens to assess batch effects ####

source("hESC_functions.R")

## Correlate individual hPSC media samples with merged sample
coefs.df <- read.table("up-tf-all_regression_trimmed.pvals.tsv", sep = "\t", stringsAsFactors = F)
cfs <- coefs.df$cf
names(cfs) <- interaction(coefs.df$Gene, coefs.df$Group, sep = "_")
top.hits <- subset(coefs.df, FDR < 0.05, abs(cf) > 0.025)

individual.files <- c("hESC_Individual_Batch_Coefficients/up-tf1_regression.pvals.tsv",
                      "hESC_Individual_Batch_Coefficients/up-tf2_regression.pvals.tsv",
                      "hESC_Individual_Batch_Coefficients/up-tf3_regression.pvals.tsv",
                      "hESC_Individual_Batch_Coefficients/up-tf4_regression.pvals.tsv",
                      "hESC_Individual_Batch_Coefficients/up-tf5_regression.pvals.tsv")

for (i in 1:length(individual.files)) {
  fi <- individual.files[[i]]
  coefs.df.batch <- read.table(fi, sep = "\t", header = T, stringsAsFactors = F)
  cfs.batch <- coefs.df.batch$cf
  names(cfs.batch) <- interaction(coefs.df.batch$Gene, coefs.df.batch$Group, sep = "_")
  top.hits.batch <- subset(coefs.df.batch, FDR < 0.05, abs(cf) > 0.025)
  
  hits.use <- unique(union(interaction(top.hits$Gene, top.hits$Group, sep = "_"), 
                           interaction(top.hits.batch$Gene, top.hits.batch$Group, sep = "_")))
  hits.use <- hits.use[hits.use %in% names(cfs)]; hits.use <- hits.use[hits.use %in% names(cfs.batch)]
  f.name <- paste("up_tf", i, sep = "")
  pdf(paste(f.name, "correlation.pdf", sep = "_"), width = 3, height = 3.25)
  print(gghexbin(cfs[hits.use], cfs.batch[hits.use], n.bins = 100))
  dev.off()
}


## Correlate multilineage media batches

coefs.df.1 <- read.table("hESC_Individual_Batch_Coefficients/up-tf-multi-batch1_regression.pvals.tsv", sep = "\t", stringsAsFactors = F)
cfs.1 <- coefs.df.1$cf
names(cfs.1) <- interaction(coefs.df.1$Gene, coefs.df.1$Group, sep = "_")

coefs.df.2 <- read.table("hESC_Individual_Batch_Coefficients/up-tf-multi-batch2_regression.pvals.tsv", sep = "\t", stringsAsFactors = F)
cfs.2 <- coefs.df.2$cf
names(cfs.2) <- interaction(coefs.df.2$Gene, coefs.df.2$Group, sep = "_")

top.hits.1 <- subset(coefs.df.1, FDR < 0.05, abs(cf) > 0.025)
top.hits.2 <- subset(coefs.df.2, FDR < 0.05, abs(cf) > 0.025)
hits.use <- unique(union(interaction(top.hits.1$Gene, top.hits.1$Group, sep = "_"), 
                         interaction(top.hits.2$Gene, top.hits.2$Group, sep = "_")))
hits.use <- hits.use[hits.use %in% names(cfs.1)]; hits.use <- hits.use[hits.use %in% names(cfs.2)]

pdf("up-tf-multi_correlation.pdf", width = 3, height = 3.25)
gghexbin(cfs.1[hits.use], cfs.2[hits.use], n.bins = 100)
dev.off()

#### Split KLF family + MYC mutants counts matrix ####

library(Seurat)
source("hESC_functions.R")

counts <- Read10X("up-tf-klf-myc/")
genotypes.list <- read.genotypes("up-tf-klf-myc_pheno_dict.csv")
length(unique(unlist(genotypes.list, F, F)))

genotypes.matrix <- as.matrix(design.matrix.genotypes(genotypes.list, max.guides = 2, min.cells = 10, drop.cells = F))
dim(genotypes.matrix)

klf.genotypes <- grep("KLF", names(genotypes.list), value = T)
myc.genotypes <- grep("MYC|Myc", names(genotypes.list), value = T)

klf.genotypes.matrix <- genotypes.matrix[rowSums(genotypes.matrix[,myc.genotypes]) == 0,]
myc.genotypes.matrix <- genotypes.matrix[rowSums(genotypes.matrix[,klf.genotypes]) == 0,]
dim(klf.genotypes.matrix)
dim(myc.genotypes.matrix)

klf.cells <- rownames(klf.genotypes.matrix)
myc.cells <- rownames(myc.genotypes.matrix)
klf.counts <- as.matrix(counts[,klf.cells])
myc.counts <- as.matrix(counts[,myc.cells])

write.table(klf.counts, file = "up-tf-klf.counts.tsv", sep = "\t")
write.table(myc.counts, file = "up-tf-myc.counts.tsv", sep = "\t")


#### MYC/KLF Screen analysis ####

## Load data
source("hESC_functions.R")

coefs.df.myc <- read.table("up-tf-myc_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.df.myc <- subset(coefs.df.myc, Group != "metadata")
top.hits.myc <- subset(coefs.df.myc, FDR < 0.05 & abs(cf) > 0.01)

coefs.df.klf <- read.table("up-tf-klf_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.df.klf <- subset(coefs.df.klf, Group != "metadata")
top.hits.klf <- subset(coefs.df.klf, FDR < 0.05 & abs(cf) > 0.01)

## TF effects on gene modules
coefs.matrix.myc <- unflatten.dataframe(coefs.df.myc, "cf")
coefs.matrix.myc <- coefs.matrix.myc[apply(coefs.matrix.myc, 1, function(x) any(abs(x) > 0.025)),]

coefs.matrix.klf <- unflatten.dataframe(coefs.df.klf, "cf")
coefs.matrix.klf <- coefs.matrix.klf[apply(coefs.matrix.klf, 1, function(x) any(abs(x) > 0.025)),]

gene_module_mapping <- read.table("up-tf_gene_module_mapping.txt", sep = "\t", header = T)
gene_modules <- gene_module_mapping$Description; names(gene_modules) <- gene_module_mapping$Module;
gene_clusters_list <- load.genesets("up-tf_functional_gene_modules.gmt")
names(gene_clusters_list) <- sapply(names(gene_clusters_list), function(x) gene_modules[[x]])

tf.modules.matrix.myc <- calc.gene_modules.effect(coefs.matrix.myc, gene_clusters_list, gene_module_mapping, min.coef = 0.025)
tf.modules.matrix.klf <- calc.gene_modules.effect(coefs.matrix.klf, gene_clusters_list, gene_module_mapping, min.coef = 0.025)

pdf("myc-mutants_gene_module_effects.pdf", width = 8.0, height = 3.75)
ggheat(tf.modules.matrix.myc, clustering = "both", x.lab.size = 12, y.lab.size = 12)
dev.off()

pdf("klf-family_gene_module_effects.pdf", width = 8.0, height = 3.75)
ggheat(tf.modules.matrix.klf, clustering = "both", x.lab.size = 12, y.lab.size = 12)
dev.off()


## Load hPSC media screen results
coefs.df.stem <- read.table("up-tf-all_regression_trimmed.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)

## Correlate MYC between myc mutant screen and hPSC media screen
TF <- "MYC"
g.coefs.df <- subset(coefs.df.stem, Group == TF)
g.coefs.df.myc <- subset(coefs.df.myc, Group == TF)

cfs <- g.coefs.df$cf; names(cfs) <- g.coefs.df$Gene;
cfs.myc <- g.coefs.df.myc$cf
names(cfs.myc) <- g.coefs.df.myc$Gene
genes.use <- intersect(names(cfs), names(cfs.myc))

pdf("MYC_hPSC_myc-mutants_correlation.pdf", width = 4.5, height = 4)
ggcorrelation(cfs[genes.use], cfs.myc[genes.use], x.lab = "hPSC screen", y.lab = "MYC mutants", title = TF, 
              use.label = F, box = T, show.corr = T)
dev.off()

## Correlate KLF4 between KLF family screen and hPSC media screen
TF <- "KLF4"
g.coefs.df <- subset(coefs.df.stem, Group == TF)
g.coefs.df.klf <- subset(coefs.df.klf, Group == TF)

cfs <- g.coefs.df$cf; names(cfs) <- g.coefs.df$Gene;
cfs.klf <- g.coefs.df.klf$cf
names(cfs.klf) <- g.coefs.df.klf$Gene
genes.use <- intersect(names(cfs), names(cfs.klf))

pdf("KLF4_hPSC_klf_family_correlation.pdf", width = 4.5, height = 4)
ggcorrelation(cfs[genes.use], cfs.klf[genes.use], x.lab = "hPSC screen", y.lab = "klf mutants", title = TF, 
              use.label = F, box = T, show.corr = T)
dev.off()