library(swne)
library(perturbLM)

## Load bulk counts
files <- list.files(path = "bulk-validation/", pattern = "*ReadsPerGene.out.tab")
tfs <- sapply(files, ExtractField, field = 2)

bulk.counts <- do.call("cbind", lapply(files, function(fi) {
  fi.counts.df <- read.table(paste0("bulk-validation/", fi), sep = "\t", header = F, row.names = 1)
  fi.counts <- fi.counts.df[[1]]; names(fi.counts) <- rownames(fi.counts.df);
  fi.counts
}))
colnames(bulk.counts) <- tfs

## Convert ensembl IDs to gene names
ensembl.df <- read.table("bulk-validation/ensembl_ids.txt", header = T, sep = "\t")
ensembl.df <- ensembl.df[!duplicated(ensembl.df$Gene.stable.ID),]
rownames(ensembl.df) <- ensembl.df$Gene.stable.ID

rownames(bulk.counts) <- sapply(rownames(bulk.counts), ExtractField, field = 1, delim = "\\.")
bulk.counts <- bulk.counts[rownames(bulk.counts) %in% rownames(ensembl.df),]
rownames(bulk.counts) <- ensembl.df[rownames(bulk.counts), "Gene.name"]
dim(bulk.counts)

## Write bulk counts to file
write.table(bulk.counts, file = "bulk-validation/up-tf-klf_bulk_validations.tsv", sep = "\t")


## Correlate with single cell results
coefs.df <- read.table("up-tf-klf_regression.pvals.tsv", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)
coefs <- UnflattenDataframe(coefs.df, "cf")

bulk.norm.counts <- log(bulk.counts/colSums(bulk.counts) + 1)
bulk.logfc <- sapply(tfs[tfs != "mCh"], function(tf) {
  log((bulk.norm.counts[,tf] + 1e-3)/(bulk.norm.counts[,"mCh"] + 1e-3))
})
colnames(bulk.logfc) <- tfs[tfs != "mCh"]
dim(bulk.logfc)


genes.use <- intersect(rownames(coefs), rownames(bulk.logfc))

tf <- "KLF5"

pdf(paste0("up-tf-klf_bulk_corr_", tf, ".pdf"), width = 4, height = 4)
PlotHexBin(bulk.logfc[genes.use, tf], coefs[genes.use, tf], n.bins = 50)
dev.off()

