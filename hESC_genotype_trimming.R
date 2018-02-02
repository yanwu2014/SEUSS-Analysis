## Trim genotypes: For each genotype, only keep cells in clusters for which the
## genotype is highly enriched.
library(perturbLM)

file.handle <- "up-tf-all" ## Name of the sample to be analyzed
# file.handle <- "up-tf-endo"
# file.handle <- "up-tf-multi"
# file.handle <- "up-tf-klf"
# file.handle <- "up-tf-myc"

max.FDR <- 1e-6 ## Set FDR threshold

## Load genotype cluster counts and calculate enrichment
df <- read.table(paste(file.handle, "_cluster_genotype_counts.tsv", sep = ""), sep = "\t")
chisq.p <- GenotypeClusterChisq(df, "mCherry")
chisq.fdr <- p.adjust(chisq.p, method = "BH")
df.pvals <- GenotypeClusterPvals(df)
df.fdr <- matrix(p.adjust(df.pvals, method = "BH"), nrow(df.pvals), ncol(df.pvals), dimnames = list(rownames(df.pvals), colnames(df.pvals)))

## Load genotypes
genotypes.input.file <- paste(file.handle, "pheno_dict.csv", sep = "_")
trimmed.genotypes.output.file <- paste(file.handle, "pheno_dict_trimmed.csv", sep = "_")

## Load Seurat object
se.obj <- readRDS(paste(file.handle, "_clustering.Robj", sep = ""))
se.obj@meta.data$tree.ident <- sapply(se.obj@meta.data$tree.ident, function(x) paste("C", x, sep = ""))

genotypes.list <- ReadGenotypes(genotypes.input.file)
sig.cl.enrich <- intersect(names(chisq.fdr[chisq.fdr < max.FDR]), rownames(df.fdr[apply(df.fdr, 1, function(x) min(x) < max.FDR), ]))
genotype.hit.clusters <- apply(df.fdr[sig.cl.enrich,], 1, function(x) names(x[x < 1e-2]))

clusters <- as.character(se.obj@ident); names(clusters) <- se.obj@cell.names;
clusters.list <- UnflattenCellGenotypes(clusters)

trim.genotypes.list <- lapply(names(genotypes.list), function(g) {
  g.cells <- genotypes.list[[g]]
  if (g %in% sig.cl.enrich) {
    cl <- genotype.hit.clusters[[g]]
    cl.cells <- unlist(clusters.list[cl], F, F)
    print(g)
    print(cl)
    print(length(g.cells))
    return(intersect(g.cells, cl.cells))
  }
  else {
    return(g.cells)
  }
})
names(trim.genotypes.list) <- names(genotypes.list)

print(sapply(genotypes.list[sig.cl.enrich], length))
print(sapply(trim.genotypes.list[sig.cl.enrich], length))
# write.genotypes(trim.genotypes.list, trimmed.genotypes.output.file)