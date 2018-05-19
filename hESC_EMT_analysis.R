## Analyze SNAI2/KLF4 effects on genes related to Epithelial-Mesenchymal transition in hPSC medium
library(pagoda2)
library(perturbLM)
library(swne)

## PCA on Hallmark EMT genes
genesets <- LoadGenesets("hallmark_go_bp_genesets.gmt")
emt.genes <- genesets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

genotypes.list <- ReadGenotypes("up-tf-stem_pheno_dict.csv")
genotypes <- FlattenGenotypeList(genotypes.list)

r <- readRDS("up-tf-stem_clustering_p2.Robj")

cells.use <- names(genotypes[genotypes %in% c("SNAI2", "KLF4")])
cells.use <- intersect(cells.use, rownames(r$counts))
emt.genes <- intersect(emt.genes, colnames(r$counts))

r$counts <- r$counts[cells.use,]
r$adjustVariance(gam.k = 10, verbose = F)
r$calculatePcaReduction(nPcs = 20, odgenes = emt.genes, verbose = F)

pc.emb <- r$reductions$PCA[,2:3]; rownames(pc.emb) <- rownames(r$counts);
clusters <- factor(genotypes[cells.use])

pdf("up-tf-stem_EMT_pca_plot.pdf", width = 5, height = 5)
PlotDims(pc.emb, sample.groups = clusters, alpha = 0.8, x.lab = "PC2", y.lab = "PC3", label.size = 6, 
         show.legend = F, pt.size = 1.5)
dev.off()

## Barplot of key EMT markers
library(ggplot2)
library(swne)

coefs.df <- read.table("up-tf-stem_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs <- UnflattenDataframe(coefs.df, "cf")

emt.genes.plot <- c("CDH1", "CDH2", "SPP1", "EPCAM", "LAMC1", "THY1", "VIM", "TPM2")

coefs.plot <- c(coefs[emt.genes.plot, "KLF4"], coefs[emt.genes.plot, "SNAI2"])
tfs.factor <- c(rep("KLF4", length(emt.genes.plot)), rep("SNAI2", length(emt.genes.plot)))
genes.factor <- rep(emt.genes.plot, 2)

barplot.df <- data.frame(coefs = coefs.plot, TF = tfs.factor, Gene = genes.factor)

pdf("up-tf-stem_EMT_barplot.pdf", width = 6, height = 4.5)
ggplot(data = barplot.df, aes(x = Gene, y = coefs)) + geom_bar(aes(fill = TF), position = "dodge", stat = "identity") + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 20), 
        legend.text = element_text(size = 16))
dev.off()
