#### Create t-SNE plots, differential expression heatmaps, and cluster enrichment heatmaps
library(perturbLM)
library(Seurat)

## Load regression results
coefs.df <- read.table("up-tf-neuron_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
top.hits <- subset(coefs.df, FDR < 0.3 & abs(cf) > 0.01)
all.genes <- unique(coefs.df$Gene)
print(sort(table(top.hits$Group), decreasing = T))

## Summarize differential expression results
genotype.counts <- sort(table(top.hits$Group), decreasing = T)
genotype.avg_diff <- tapply(top.hits$cf, top.hits$Group, function(x) mean(abs(x)))[names(genotype.counts)]

genotypes.matrix <- DesignMatrixGenotypes(ReadGenotypes("up-tf-neuron_pheno_dict.csv"), 
                                          max.genotypes = 3, min.cells = 15)
genotype.sizes <- Matrix::colSums(genotypes.matrix)
genotype.summary <- data.frame(n_diff_genes = as.integer(genotype.counts), avg_eff_size = genotype.avg_diff)
genotype.summary$n_cells <- genotype.sizes[rownames(genotype.summary)]
# write.table(genotype.summary, file = paste(file.handle, "_genotype_summary.tsv", sep = ""), sep = "\t")

## Determine significant genotypes
# genotype.summary <- read.table(paste(file.handle, "_genotype_summary.tsv", sep = ""), header = T, sep = "\t", stringsAsFactors = F)

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
