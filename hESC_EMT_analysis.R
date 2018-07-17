## Analyze SNAI2/KLF4 effects on genes related to Epithelial-Mesenchymal transition in hPSC medium

## Barplot of key EMT markers
library(ggplot2)
library(swne)

coefs.df <- read.table("up-tf-stem_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs <- UnflattenDataframe(coefs.df, "cf")

emt.genes.plot <- c("CDH1", "EPCAM", "LAMC1", "SPP1", "THY1", "VIM", "TPM2")

coefs.plot <- c(coefs[emt.genes.plot, "KLF4"], coefs[emt.genes.plot, "SNAI2"])
tfs.factor <- c(rep("KLF4", length(emt.genes.plot)), rep("SNAI2", length(emt.genes.plot)))
genes.factor <- rep(emt.genes.plot, 2)

barplot.df <- data.frame(coefs = coefs.plot, TF = tfs.factor, Gene = genes.factor)

pdf("up-tf-stem_EMT_barplot.pdf", width = 5, height = 4.5)
ggplot(data = barplot.df, aes(x = Gene, y = coefs)) + geom_bar(aes(fill = TF), position = "dodge", stat = "identity") + 
  theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                          axis.text.x = element_text(angle = 90, hjust = 1, size = 20), 
                          legend.text = element_text(size = 16))
dev.off()


## PCA on Hallmark EMT genes
library(Seurat)
library(perturbLM)
library(swne)
library(ggplot2)

genesets <- LoadGenesets("hallmark_go_bp_genesets.gmt")
emt.genes <- union(genesets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, genesets$GO_EPITHELIAL_TO_MESENCHYMAL_TRANSITION)

genotypes.list <- ReadGenotypes("up-tf-stem_pheno_dict.csv")
genotypes <- FlattenGenotypeList(genotypes.list)

r <- readRDS("up-tf-stem_clustering.Robj")
r <- ScaleData(r, genes.use = emt.genes, vars.to.regress = c("nUMI", "batch", "percent.mt"), model.use = "negbinom",
               scale.max = 4)

cells.use <- names(genotypes[genotypes %in% c("SNAI2", "KLF4")])
cells.use <- intersect(cells.use, colnames(r@data))
emt.genes <- intersect(emt.genes, rownames(r@data))

r <- SubsetData(r, cells.use = cells.use)
r <- RunPCA(r, pc.genes = emt.genes, do.print = F)

pc.emb <- GetCellEmbeddings(r, reduction.type = "pca", dims.use = 1:3)
clusters <- factor(genotypes[cells.use], levels = c("SNAI2", "KLF4"))

pdf("up-tf-stem_EMT_pca_plot.pdf", width = 4.5, height = 4.5)
PlotDims(pc.emb[,c(1,2)], sample.groups = clusters, alpha = 0.75, x.lab = "Hallmark EMT: PC1", 
         y.lab = "Hallmark EMT: PC2", label.size = 6, show.legend = F, pt.size = 2, seed = 352356)
dev.off()