#### Split KLF family + MYC mutants counts matrix ####
library(Seurat)
library(perturbLM)

counts <- ReadData("up-tf-klf-myc/")
genotypes.list <- ReadGenotypes("up-tf-klf-myc_pheno_dict.csv")
length(unique(unlist(genotypes.list, F, F)))

genotypes.matrix <- as.matrix(DesignMatrixGenotypes(genotypes.list, max.genotypes = 2, 
                                                    min.cells = 10, drop.cells = F))
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

dim(klf.counts)
dim(myc.counts)

# write.table(klf.counts, file = "up-tf-klf.counts.tsv", sep = "\t")
# write.table(myc.counts, file = "up-tf-myc.counts.tsv", sep = "\t")


#### MYC/KLF Differential Expr & Gene Module analysis ####
library(perturbLM)
library(swne)

file.handle <- "up-tf-klf"
# file.handle <- "up-tf-myc"

## Load data
coefs.df <- read.table(paste(file.handle, "regression.pvals.tsv", sep = "_"), sep = "\t", header = T, 
                           stringsAsFactors = F)
genotype.summary <- read.table(paste(file.handle, "summary.tsv", sep = "_"), sep = "\t", header = T, row.names = 1)

## Diff exp barplot
diff.genes <- genotype.summary$n_diff_genes
names(diff.genes) <- rownames(genotype.summary)

if (file.handle == "up-tf-klf") {
  klf.families <- c("Group1", "Group1", "Group1", "Group1", "Group1", "Group1", "Group2", "Group2", "Group2",
                    "Group3", "Group3", "Group3", "Group3", "Group3", "Group3", "Group4", "Group4")
  names(klf.families) <- c("KLF1", "KLF5", "KLF6", "KLF2", "KLF4", "KLF7", "KLF3", "KLF8", "KLF12", "KLF9",
                           "KLF10", "KLF11", "KLF13", "KLF14", "KLF16", "KLF15", "KLF17")
  klf.families <- factor(klf.families)
  klf.order.df <- data.frame(diff.genes, family = klf.families[names(diff.genes)])
  klf.order.df <- klf.order.df[order(klf.order.df$family),]
  diff.genes <- diff.genes[rownames(klf.order.df)]
}

gg.bar <- ggBarplot(diff.genes, fill.color = "lightgrey")
pdf(paste(file.handle, "diff_genes_barplot.pdf", sep = "_"), width = 5.25, height = 3.5)
gg.bar
dev.off()

## TF effects on gene modules
gene.module.mapping <- read.table("up-tf_module_names.txt", sep = "\t", header = T)
colnames(gene.module.mapping) <- c("Module", "Description")
# gene.module.mapping <- NULL
gene.modules.list <- LoadGenesets("up-tf_gene_modules.gmt")

coefs.matrix <- UnflattenDataframe(coefs.df, "cf")
tf.modules.matrix <- CalcGeneModuleEffect(coefs.matrix, gene.modules.list, 
                                          gene.module.mapping, min.coef = 0.01)
tf.modules.matrix <- tf.modules.matrix[apply(abs(tf.modules.matrix), 1, max) > 0.01,]
tf.modules.matrix <- tf.modules.matrix[,names(diff.genes)]

pdf(paste(file.handle, "heatmap_labels.pdf", sep = "_"), width = 6.0, height = 4.5)
ggHeat(tf.modules.matrix, clustering = "row", x.lab.size = 12, y.lab.size = 12)
dev.off()

rownames(tf.modules.matrix) <- colnames(tf.modules.matrix) <- NULL
gg.heat <- ggHeat(tf.modules.matrix, clustering = "row", x.lab.size = 12, y.lab.size = 12)

library(grid)
library(gridExtra)
library(ggplot2)

gg.bar <- gg.bar + theme(axis.text = element_blank())
gg.heat <- gg.heat + theme(legend.position = "none")

pdf(paste(file.handle, "stacked_barplot_heatmap_nolabels.pdf", sep = "_"), width = 6.0, height = 5)
grid.arrange(gg.bar, gg.heat, ncol = 1)
dev.off()


#### Endogenous MYC expression analysis ####
library(swne)

myc.counts <- ReadData("up-tf-myc.counts.tsv")
myc.counts <- FilterData(myc.counts, min.samples.frac = 0.001, trim = 0.01, min.nonzero.features = 0)
myc.counts <- ScaleCounts(myc.counts)

myc.expr <- c(myc.counts["MYC",], myc.counts["MYCL",], myc.counts["MYCN",])
gene.factor <- c(rep("MYC", ncol(myc.counts)), rep("MYCL", ncol(myc.counts)), rep("MYCN", ncol(myc.counts)))
myc.expr.df <- data.frame(Expr = myc.expr, Gene = factor(gene.factor))

pdf("MYC_endogenous_expression.pdf", width = 4.5, height = 3.5)
ggplot(myc.expr.df, aes(Gene, Expr)) + geom_violin(aes(fill = Gene))
dev.off()

