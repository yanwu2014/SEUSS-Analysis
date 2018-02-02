## Creates Figure 1f in Parekh, Wu et al

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