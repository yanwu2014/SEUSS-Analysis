#### Create t-SNE plots, differential expression heatmaps, and cluster enrichment heatmaps
library(perturbLM)
library(Seurat)

## Load regression results
coefs.df <- read.table("up-tf-neuron-nohygro_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
top.hits <- subset(coefs.df, FDR < 0.05 & abs(cf) > 0.025)

top.genotypes <- sort(table(top.hits$Group), decreasing = T)
top.genotypes <- top.genotypes[grepl(":", names(top.genotypes))]
head(top.genotypes, n = 10)

## Summarize differential expression results
genotype.counts <- sort(table(top.hits$Group), decreasing = T)
genotype.avg_diff <- tapply(top.hits$cf, top.hits$Group, function(x) mean(abs(x)))[names(genotype.counts)]

genotypes.matrix <- DesignMatrixGenotypes(ReadGenotypes("up-tf-neuron-nohygro_pheno_dict.csv"), 
                                          max.genotypes = 2, min.cells = 1)
genotype.sizes <- Matrix::colSums(genotypes.matrix)
genotype.summary <- data.frame(n_diff_genes = as.integer(genotype.counts), avg_eff_size = genotype.avg_diff)
genotype.summary$n_cells <- genotype.sizes[rownames(genotype.summary)]
# write.table(genotype.summary, file = "up-tf-neuron-nohygro_summary.tsv", sep = "\t")

## Combo Barplot
head(top.genotypes, n = 10)
sort(genotype.sizes[grepl(":", names(genotype.sizes))], decreasing = T)

grp.1 <- "OTX2"
grp.2 <- "ASCL5"
combo.grp <- paste(grp.1, grp.2, sep = ":")
combo.tfs.hits <- subset(top.hits, Group == combo.grp); rownames(combo.tfs.hits) <- combo.tfs.hits$Gene;
combo.tfs.hits <- subset(combo.tfs.hits, !grepl("LINC|AC0|RP11", Gene))

# genes.use <- head(combo.tfs.hits[order(abs(combo.tfs.hits$cf), decreasing = T), "Gene"], n = 12)
genes.use <- combo.tfs.hits[,"Gene"]

cfs.grp.1 <- subset(coefs.df, Group == grp.1); rownames(cfs.grp.1) <- cfs.grp.1$Gene;
cfs.grp.2 <- subset(coefs.df, Group == grp.2); rownames(cfs.grp.2) <- cfs.grp.2$Gene;
cfs.grp.1 <- cfs.grp.1[genes.use, "cf"]; names(cfs.grp.1) <- genes.use;
cfs.grp.2 <- cfs.grp.2[genes.use, "cf"]; names(cfs.grp.2) <- genes.use;
cfs.expect <- cfs.grp.1 + cfs.grp.2

cfs.combo <- combo.tfs.hits[genes.use, "cf"]; names(cfs.combo) <- genes.use;
cfs.combo <- sort(cfs.combo, decreasing = T); cfs.expect <- cfs.expect[names(cfs.combo)]
cfs.combo <- cfs.combo + cfs.expect

cfs.ratio <- abs(cfs.combo/cfs.expect)
genes.plot <- names(head(sort(cfs.ratio, decreasing = T), n = 12))
cfs.combo <- cfs.combo[genes.plot]; cfs.expect <- cfs.expect[genes.plot];

cfs.plot <- c(cfs.combo, cfs.expect)
genes.plot <- rep(names(cfs.combo), 2)
cfs.type <- factor(c(rep("Observed", length(cfs.combo)), rep("Expected", length(cfs.combo))), 
                   levels = c("Observed", "Expected"))
diff.genes.barplot.df <- data.frame(Cf = cfs.plot, Coef = cfs.type, Gene = factor(genes.plot, levels = names(cfs.combo)))

pdf(paste("up-tf-neuron", grp.1, grp.2, "barplot_ratio.pdf", sep = "_"), width = 6.5, height = 4.5)
ggplot(data = diff.genes.barplot.df, aes(x = Gene, y = Cf, fill = Coef, width = 0.7)) +
  geom_bar(position = "dodge", stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(hjust = 1, angle = 90, size = 14), 
        axis.text.y = element_text(size = 12))
dev.off()


## Interaction diff gene barplot
top.genotypes.counts <- as.integer(top.genotypes)
names(top.genotypes.counts) <- names(top.genotypes)

pdf("up-tf-neuron_interaction_diff_genes.pdf", width = 6.5, height = 5)
ggBarplot(top.genotypes.counts, fill.color = "grey")
dev.off()


## Shuffling rate barplot
library(perturbLM)
shuff.rate <- c(0.6, 36.3, 0.8)
names(shuff.rate) <- c("Unpooled vector", "Pooled vector", "Pooled vector noHygro")

pdf("up-tf_shuffling_barplot.pdf", width = 1.5, height = 4.5)
ggBarplot(shuff.rate, y.lim = c(0,100), fill.color = "grey")
dev.off()
