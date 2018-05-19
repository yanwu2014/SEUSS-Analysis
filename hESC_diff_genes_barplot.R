## Creates Figure 1f in Parekh, Wu et al
library(ggplot2)
library(swne)

min.cl.pval <- 1e-12
min.diff.genes <- 50

## Load stem cell media results
genotype.summary.stem <- read.table("up-tf-stem_summary.tsv", header = T, sep = "\t", stringsAsFactors = F)
genotypes.plot.stem <- union(rownames(subset(genotype.summary.stem, abs(min_cluster_pval) < min.cl.pval)),
                            rownames(subset(genotype.summary.stem, n_diff_genes >= min.diff.genes)))
genotypes.plot.stem <- sort(unique(c(genotypes.plot.stem, "mCherry-int-ctrl"))); print(genotypes.plot.stem);

## Load endothelial results
genotype.summary.endo <- read.table("up-tf-endo_summary.tsv", sep = "\t")
genotypes.plot.endo <- union(rownames(subset(genotype.summary.endo, abs(min_cluster_pval) < min.cl.pval)),
                            rownames(subset(genotype.summary.endo, n_diff_genes >= min.diff.genes)))
genotypes.plot.endo <- sort(unique(c(genotypes.plot.endo, "mCherry-int-ctrl"))); print(genotypes.plot.endo);

## Load multilineage media results
genotype.summary.multi <- read.table("up-tf-multi_summary.tsv", sep = "\t")
genotypes.plot.multi <- union(rownames(subset(genotype.summary.multi, abs(min_cluster_pval) < min.cl.pval)),
                             rownames(subset(genotype.summary.multi, n_diff_genes >= min.diff.genes)))
genotypes.plot.multi <- sort(unique(c(genotypes.plot.multi, "mCherry-int-ctrl"))); print(genotypes.plot.multi);


## Stacked barplot summary
genotypes.plot <- unique(c(genotypes.plot.stem, genotypes.plot.multi, genotypes.plot.endo))
n.sig <- length(genotypes.plot); print(n.sig);

all.diff.genes <- c(genotype.summary.stem[genotypes.plot, "n_diff_genes"],
                    genotype.summary.endo[genotypes.plot, "n_diff_genes"],
                    genotype.summary.multi[genotypes.plot, "n_diff_genes"])
all.diff.genes[is.na(all.diff.genes)] <- 0
media.types <- c(rep("Stem cell", n.sig), rep("Endothelial", n.sig), rep("Multilineage", n.sig))
barplot.df <- data.frame(Genotype = rep(genotypes.plot, 3), Diff_expr_genes = all.diff.genes, Media = media.types)

pdf("media_condition_summary_barplot.pdf", width = 6, height = 7)
ggplot(barplot.df, aes(x = Genotype, y = Diff_expr_genes, fill = Media)) + geom_bar(stat = "identity") + theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16)) + 
  coord_flip()
dev.off()