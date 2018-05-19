#### Create t-SNE plots, differential expression heatmaps, and cluster enrichment heatmaps
library(perturbLM)
library(Seurat)

## Load regression results
coefs.df <- read.table("up-tf-neuron_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
top.hits <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025)

top.genotypes <- sort(table(top.hits$Group), decreasing = T)
top.genotypes <- top.genotypes[grepl(":", names(top.genotypes))]
head(top.genotypes, n = 10)

## Summarize differential expression results
genotype.counts <- sort(table(top.hits$Group), decreasing = T)
genotype.avg_diff <- tapply(top.hits$cf, top.hits$Group, function(x) mean(abs(x)))[names(genotype.counts)]

genotypes.matrix <- DesignMatrixGenotypes(ReadGenotypes("up-tf-neuron_pheno_dict.csv"), 
                                          max.genotypes = 3, min.cells = 15)
genotype.sizes <- Matrix::colSums(genotypes.matrix)
genotype.summary <- data.frame(n_diff_genes = as.integer(genotype.counts), avg_eff_size = genotype.avg_diff)
genotype.summary$n_cells <- genotype.sizes[rownames(genotype.summary)]
write.table(genotype.summary, file = "up-tf-neuron_summary.tsv", sep = "\t")

## Combo Barplot
grp.1 <- "OTX2"
grp.2 <- "NEUROG3"
combo.grp <- paste(grp.1, grp.2, sep = ":")
combo.tfs.hits <- subset(top.hits, Group == combo.grp); rownames(combo.tfs.hits) <- combo.tfs.hits$Gene;
combo.tfs.hits <- subset(combo.tfs.hits, !grepl("LINC|AC0", Gene))

genes.use <- head(combo.tfs.hits[order(abs(combo.tfs.hits$cf), decreasing = T), "Gene"], n = 12)

cfs.grp.1 <- subset(coefs.df, Group == grp.1); rownames(cfs.grp.1) <- cfs.grp.1$Gene;
cfs.grp.2 <- subset(coefs.df, Group == grp.2); rownames(cfs.grp.2) <- cfs.grp.2$Gene;
cfs.grp.1 <- cfs.grp.1[genes.use, "cf"]; names(cfs.grp.1) <- genes.use;
cfs.grp.2 <- cfs.grp.2[genes.use, "cf"]; names(cfs.grp.2) <- genes.use;
cfs.expect <- cfs.grp.1 + cfs.grp.2

cfs.combo <- combo.tfs.hits[genes.use, "cf"]; names(cfs.combo) <- genes.use;
cfs.combo <- sort(cfs.combo, decreasing = T); cfs.expect <- cfs.expect[names(cfs.combo)]
cfs.combo <- cfs.combo + cfs.expect

cfs.plot <- c(cfs.combo, cfs.expect)
genes.plot <- rep(names(cfs.combo), 2)
cfs.type <- factor(c(rep("Observed", length(cfs.combo)), rep("Expected", length(cfs.combo))), 
                   levels = c("Observed", "Expected"))
diff.genes.barplot.df <- data.frame(Cf = cfs.plot, Coef = cfs.type, Gene = factor(genes.plot, levels = names(cfs.combo)))

pdf(paste("up-tf-neuron", grp.1, grp.2, "barplot.pdf", sep = "_"), width = 6.5, height = 4.5)
ggplot(data = diff.genes.barplot.df, aes(x = Gene, y = Cf, fill = Coef, width = 0.7)) +
  geom_bar(position = "dodge", stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(hjust = 1, angle = 90, size = 14), 
        axis.text.y = element_text(size = 12))
dev.off()

