#### Orthologous mouse TF enrichment in hPSC media. Generates Figure S3d

library(liger)
library(perturbLM)

## Create mouse TF targets genesets
mouse.targets.df <- read.table("mouse_tf_targets.txt", sep = "\t", header = T, stringsAsFactors = F)
genotypes.keep <- c("Ascl1", "Cdx2", "Klf4", "Otx2", "Myc", "Myod1", "Otx2")

mouse.targets.df <- subset(mouse.targets.df, TF %in% genotypes.keep & log.fc > 0.2 & FDR < 0.05)
mouse.targets.df <- mouse.targets.df[order(mouse.targets.df$FDR),]
mouse.targets.df <- Reduce(rbind, by(mouse.targets.df, mouse.targets.df$TF, head, n = 600))
print(table(mouse.targets.df$TF))

mouse.genesets <- lapply(unique(mouse.targets.df$TF), function(tf) subset(mouse.targets.df, TF == tf)$Genes)
names(mouse.genesets) <- sapply(unique(mouse.targets.df$TF), toupper)

library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", verbose = F)
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", verbose = F)

human.genesets <- lapply(mouse.genesets, function(mousex) {
  human.genes.df = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mousex ,
                          mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  unique(human.genes.df[,2])
})

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
genesets.genes.counts <- sapply(genesets.genes, function(gene) sum(sapply(human.genesets.trimmed, 
                                                                          function(genes) gene %in% genes)))
genesets.genes.counts <- genesets.genes.counts[genesets.genes.counts <= max.genesets]
human.genesets.trimmed <- lapply(human.genesets.trimmed, function(genes) genes[genes %in% names(genesets.genes.counts)])

sapply(human.genesets, length)
sapply(human.genesets.trimmed, length)
# write.genesets(human.genesets.trimmed, file.name = "hESC_bulk_mouse_ortholog_tf_genesets.gmt")
rm(se.obj); invisible(gc());

## GSEA enrichment for mouse ortholog genesets
coefs.df <- read.table("up-tf-all_regression_trimmed.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
genotype.summary <- read.table("up-tf-all_genotype_summary.tsv", header = T, sep = "\t", stringsAsFactors = F)
sig.genotypes <- union(rownames(subset(genotype.summary, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                       rownames(subset(genotype.summary, n_diff_genes >= 100)))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry")))
all.genes <- unique(coefs.df$Gene)
n.cores <- 16

genesets <- LoadGenesets("hESC_bulk_mouse_ortholog_tf_genesets.gmt")
genesets <- FilterGenesets(genesets, all.genes, min.size = 10, max.size = 2000)
sapply(genesets, length)

scores.list <- lapply(sig.genotypes, function(g) {
  g.coefs.df <- subset(coefs.df, Group == g)
  g.coefs <- g.coefs.df$cf; names(g.coefs) <- g.coefs.df$Gene;
  g.coefs
}); names(scores.list) <- sig.genotypes;

genesets.enrich <- MultipleGSEAEnrich(scores.list, genesets, n.rand = 10000, n.cores = n.cores, power = 1)

top.genesets.enrich <- subset(genesets.enrich, FDR < 5e-2)
heat.df <- subset(genesets.enrich, Group %in% sig.genotypes & genesets %in% top.genesets.enrich$genesets)
heat.df$Score <- mapply(function(p, cf) -log(p) * sign(cf), heat.df$p.val, heat.df$sscore)
heat.mat <- UnflattenDataframe(heat.df, output.name = "Score", row.col = 'genesets', col.col = 'Group')
heat.mat[is.na(heat.mat)] <- 0

# pdf("bulk_ortholog_gsea_enrichment.pdf", width = 6, height = 2.5)
ggheat(heat.mat, clustering = "col", x.lab.size = 12, y.lab.size = 12)
# dev.off()

# write.table(genesets.enrich, file = "up-tf-all_TF_gene_module_enrichment.tsv", sep = "\t")