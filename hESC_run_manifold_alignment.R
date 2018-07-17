library(methods)
library(Seurat)
library(swne)
library(perturbLM)

## Command line parameters (can hardcode in script if desired)
options <- commandArgs(trailingOnly = T)
matrix.dir <- options[[1]]
output.file <- options[[2]]

## Other parameters
min.expr <- 0.1 ## Minimum average expression for variable genes
min.var.z <- 0.5 ## Minimum Z-score for variable genes

## Load the dataset
counts <- ReadData(matrix.dir)

## Pull out teratoma id
batch <- factor(sapply(colnames(counts), function(x) strsplit(x, split = "\\.")[[1]][[2]]))
levels(batch)

## Create and setup seurat objects for each teratoma
obj.list <- lapply(levels(batch), function(i) {
  batch.counts <- counts[,names(batch[batch == i])]

  mito.genes <- grep("^MT-", rownames(batch.counts), value = T)
  percent.mt <- Matrix::colSums(batch.counts[mito.genes, ])/Matrix::colSums(batch.counts)
  batch.counts <- batch.counts[!rownames(batch.counts) %in% mito.genes,]

  obj <- CreateSeuratObject(batch.counts, min.genes = 200, min.cells = round(0.01*ncol(batch.counts)),
                            normalization.method = "LogNormalize")
  obj <- AddMetaData(obj, percent.mt, col.name = "percent.mt")
  obj <- FindVariableGenes(obj, x.low.cutoff = min.expr, y.cutoff = min.var.z)
  obj <- ScaleData(obj, model.use = "negbinom", vars.to.regress = c("nUMI", "percent.mt"))
  obj <- RunPCA(obj, pc.genes = obj@var.genes, pcs.compute = 30)
  obj@meta.data$batch <- paste0("rep", i)
  obj
})

## Determine genes to use for CCA, must be highly variable in at least 2 datasets
var.genes <- c()
for (i in 1:length(obj.list)) {
  var.genes <- c(var.genes, obj.list[[i]]@var.genes)
}
var.genes.counts <- table(var.genes)
genes.use <- names(var.genes.counts[var.genes.counts >= 1])
common.genes <- Reduce(intersect, lapply(obj.list, function(obj) obj@var.genes))
genes.use <- intersect(genes.use, common.genes)
length(genes.use)

if (length(obj.list) >= 3) {
  obj.integrated <- RunMultiCCA(obj.list, genes.use = genes.use, num.ccs = 20)
} else {
  obj.integrated <- RunCCA(obj.list[[1]], obj.list[[2]], group.by = "batch", num.cc = 20,
                           genes.use = genes.use)
}

MetageneBicorPlot(obj.integrated, grouping.var = "batch", dims.eval = 1:20)

## Save output
print("Done with CCA")
save.image(output.file)

## Run rare non-overlapping filtering
ccs.use <- 20
obj.integrated <- CalcVarExpRatio(object = obj.integrated, reduction.type = "pca",
                                  grouping.var = "batch", dims.use = 1:ccs.use)
obj.integrated <- SubsetData(obj.integrated, subset.name = "var.ratio.pca", accept.low = 0.25)

## Alignment
obj.integrated <- AlignSubspace(obj.integrated, reduction.type = "cca", grouping.var = "batch",
                                dims.align = 1:ccs.use)

## Clustering
cluster.res <- 1.8
obj.integrated <- FindClusters(obj.integrated, reduction.type = "cca.aligned", dims.use = 1:ccs.use,
                               k.param = 30, save.SNN = T, resolution = cluster.res)

## Cluster evaluation and merging
min.conn <- 1e-4
conn <- Seurat:::CalcConnectivity(obj.integrated)
print(sum(conn > min.conn))

obj.integrated <- RunPCA(obj.integrated, pcs.compute = 20, do.print = F)
obj.integrated <- ValidateClusters(obj.integrated, pc.use = 1:20, top.genes = 15, min.connectivity = min.conn,
                                   acc.cutoff = 0.95, verbose = T)
obj.integrated <- BuildClusterTree(obj.integrated, do.reorder = T, reorder.numeric = T)

## Save output
save.image(output.file)

## Visualization
obj.integrated <- RunTSNE(obj.integrated, reduction.use = "cca.aligned",
                          dims.use = 1:ccs.use)
TSNEPlot(obj.integrated, pt.size = 0.5, do.label = T)
TSNEPlot(obj.integrated, group.by = "batch", pt.size = 0.5, do.label = F)

save.image(output.file)
saveRDS(obj.integrated, file = gsub(".RData", ".Robj", output.file))