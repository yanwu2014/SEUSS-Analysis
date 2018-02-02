## Analyze SNAI2/KLF4 effects on genes related to Epithelial-Mesenchymal transition in hPSC medium

library(Seurat)
library(perturbLM)

## PCA on Hallmark EMT genes
hallmark.genesets <- LoadGenesets("hallmark_genesets.gmt")
emt.genes <- hallmark.genesets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

genotypes.trimmed <- FlattenGenotypeList(ReadGenotypes("up-tf-all_pheno_dict_trimmed.csv"))
cells.use <- names(genotypes.trimmed[genotypes.trimmed %in% c("SNAI2", "KLF4")])

se.obj <- readRDS("up-tf-all_clustering.Robj")
se.obj <- SubsetData(se.obj, cells.use = cells.use)
se.obj <- ScaleData(se.obj, genes.use = emt.genes, model.use = "negbinom")
se.obj <- RunPCA(se.obj, pc.genes = emt.genes, pcs.compute = 10, do.print = F)
se.obj@ident <- factor(genotypes.trimmed[se.obj@cell.names])

se.obj.pc.scores <- PCAEmbed(se.obj, dims.use = 1:2)

# pdf("EMT_transition_KLF4_SNAI2.pdf", width = 5, height = 5)
PlotDims(se.obj.pc.scores[,1], se.obj.pc.scores[,2], clusters = factor(se.obj@ident), alpha = 1,
         label.size = 8, show.legend = F)
# dev.off()

se.obj.pc.loadings <- sort(PCALoad(se.obj, dims.use = 1)[,1], decreasing = T)
print(head(se.obj.pc.loadings, n = 40))
print(tail(se.obj.pc.loadings, n = 40))

## Barplot of key markers
coefs.df <- read.table("up-tf-all_regression_trimmed.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
genotype.summary <- read.table("up-tf-all_genotype_summary.tsv", header = T, sep = "\t", stringsAsFactors = F)
sig.genotypes <- union(rownames(subset(genotype.summary, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                       rownames(subset(genotype.summary, n_diff_genes >= 100)))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry")))

coefs.matrix.stem <- UnflattenDataframe(coefs.df, "cf", row.col = "Gene", col.col = "Group")
coefs.matrix.stem <- coefs.matrix.stem[, sig.genotypes]

emt.genes.plot <- c("SPP1", "EPCAM", "LAMC1", "THY1", "VIM", "TPM2")

coefs.plot <- c(coefs.matrix.stem[emt.genes.plot, "KLF4"], coefs.matrix.stem[emt.genes.plot, "SNAI2"])
tfs.factor <- c(rep("KLF4", length(emt.genes.plot)), rep("SNAI2", length(emt.genes.plot)))
genes.factor <- rep(emt.genes.plot, 2)

barplot.df <- data.frame(coefs = coefs.plot, TF = tfs.factor, Gene = genes.factor)

# pdf("EMT_markers_barplot.pdf", width = 6, height = 4.5)
ggplot(data = barplot.df, aes(x = Gene, y = coefs)) + geom_bar(aes(fill = TF), position = "dodge", stat = "identity") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 20), legend.text = element_text(size = 16))
# dev.off()