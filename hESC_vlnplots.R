library(perturbLM)
library(pagoda2)
library(ggplot2)
library(swne)

## Sample to analyze
file.handle <- "up-tf-stem"

## Load pagoda2 object
r <- readRDS(paste(file.handle, "clustering_p2.Robj", sep = "_"))

## Load genotypes
genotypes.list <- ReadGenotypes(paste(file.handle, "pheno_dict.csv", sep = "_"))
genotypes <- FlattenGroups(genotypes.list)
genotypes <- genotypes[names(genotypes) %in% rownames(r$counts)]

## Select significant genotypes
min.cl.pval <- 1e-12
min.diff.genes <- 50

genotype.summary <- read.table(paste(file.handle, "summary.tsv", sep = "_"), sep = "\t", header = T, row.names = 1)
sig.genotypes <- rownames(subset(genotype.summary, abs(min_cluster_pval) < min.cl.pval | n_diff_genes >= min.diff.genes))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry"))); print(sig.genotypes);

genotypes <- genotypes[genotypes %in% sig.genotypes]
norm.counts <- t(r$counts[names(genotypes),])

## Load regression results
coefs.df <- read.table(paste(file.handle, "regression.pvals.tsv", sep = "_"), sep = "\t", header = T)
top.pos.coefs.df <- subset(coefs.df, FDR < 0.05 & cf > 0.025)
top.pos.coefs.df <- top.pos.coefs.df[order(top.pos.coefs.df$cf, decreasing = T),]

## Load gene module results
gene.modules.list <- ReadGenesets("up-tf_gene_modules.gmt")
modules.names.df <- read.table("up-tf_module_names.txt", sep = "\t", header = T)
modules.names <- modules.names.df$Annotation; names(modules.names) <- modules.names.df$Module;
names(gene.modules.list) <- modules.names[names(gene.modules.list)]

genotype <- "CDX2"
module <- "Cell Proliferation"
genotype.module.coefs <- subset(top.pos.coefs.df, Group == genotype & Gene %in% gene.modules.list[[module]])
head(genotype.module.coefs, n = 10)

gene <- "KRT8"

pdf(paste(file.handle, gene, "vlnplot.pdf", sep = "_"), width = 7.5, height = 3.5)
PlotViolin(gene, norm.counts[gene,], factor(genotypes), x.lab.rot = T, y.log = F, x.title = "Genotype",
           remove.legend = T, size.x.use = 12, size.title.use = 12, point.size.use = -1)
dev.off()
