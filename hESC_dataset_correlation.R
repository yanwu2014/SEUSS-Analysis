#### Correlate stem cell media replicates ####
library(perturbLM)
library(swne)

stem.coefs.df.1 <- read.table("Individual_Batch_Analysis/up-tf2_regression.pvals.tsv", 
                               sep = "\t", stringsAsFactors = F)
stem.cfs.1 <- stem.coefs.df.1$cf
rownames(stem.coefs.df.1) <- names(stem.cfs.1) <- interaction(stem.coefs.df.1$Gene, stem.coefs.df.1$Group, sep = "_")

stem.coefs.df.2 <- read.table("Individual_Batch_Analysis/up-tf3_regression.pvals.tsv", 
                               sep = "\t", stringsAsFactors = F)
stem.cfs.2 <- stem.coefs.df.2$cf
rownames(stem.coefs.df.2) <- names(stem.cfs.2) <- interaction(stem.coefs.df.2$Gene, stem.coefs.df.2$Group, sep = "_")

stem.common.idx <- intersect(names(stem.cfs.1), names(stem.cfs.2))
stem.sig.idx <- union(rownames(subset(stem.coefs.df.1, FDR < 0.5)), rownames(subset(stem.coefs.df.2, FDR < 0.5)))
stem.hits.idx <- intersect(stem.common.idx, stem.sig.idx)

pdf("up-tf-stem_corr_stem-1_vs_stem-2.pdf", width = 3, height = 3)
PlotHexBin(stem.cfs.1[stem.hits.idx], stem.cfs.2[stem.hits.idx], n.bins = 25)
dev.off()



#### Correlate multilineage media batches ####
library(perturbLM)
library(swne)

multi.coefs.df.1 <- read.table("Individual_Batch_Analysis/up-tf-multi-batch1_regression.pvals.tsv", 
                               sep = "\t", stringsAsFactors = F)
multi.cfs.1 <- multi.coefs.df.1$cf
rownames(multi.coefs.df.1) <- names(multi.cfs.1) <- interaction(multi.coefs.df.1$Gene, multi.coefs.df.1$Group, sep = "_")

multi.coefs.df.2 <- read.table("Individual_Batch_Analysis/up-tf-multi-batch2_regression.pvals.tsv", 
                               sep = "\t", stringsAsFactors = F)
multi.cfs.2 <- multi.coefs.df.2$cf
rownames(multi.coefs.df.2) <- names(multi.cfs.2) <- interaction(multi.coefs.df.2$Gene, multi.coefs.df.2$Group, sep = "_")

multi.common.idx <- intersect(names(multi.cfs.1), names(multi.cfs.2))
multi.sig.idx <- union(rownames(subset(multi.coefs.df.1, FDR < 0.5)), rownames(subset(multi.coefs.df.2, FDR < 0.5)))
multi.hits.idx <- intersect(multi.common.idx, multi.sig.idx)

pdf("up-tf-multi_corr_multi-1_vs_multi-2.pdf", width = 3, height = 3)
PlotHexBin(multi.cfs.1[multi.hits.idx], multi.cfs.2[multi.hits.idx], n.bins = 25)
dev.off()


#### Correlate endothelial media batches ####
library(perturbLM)
library(swne)

endo.coefs.df.1 <- read.table("Individual_Batch_Analysis/up-tf-endo-batch1_regression.pvals.tsv", 
                              sep = "\t", stringsAsFactors = F)
endo.cfs.1 <- endo.coefs.df.1$cf
rownames(endo.coefs.df.1) <- names(endo.cfs.1) <- interaction(endo.coefs.df.1$Gene, endo.coefs.df.1$Group, sep = "_")

endo.coefs.df.2 <- read.table("Individual_Batch_Analysis/up-tf-endo-batch2_regression.pvals.tsv", 
                              sep = "\t", stringsAsFactors = F)
endo.cfs.2 <- endo.coefs.df.2$cf
rownames(endo.coefs.df.2) <- names(endo.cfs.2) <- interaction(endo.coefs.df.2$Gene, endo.coefs.df.2$Group, sep = "_")

endo.common.idx <- intersect(names(endo.cfs.1), names(endo.cfs.2))
endo.sig.idx <- union(rownames(subset(endo.coefs.df.1, FDR < 0.5)), rownames(subset(endo.coefs.df.2, FDR < 0.5)))
endo.hits.idx <- intersect(endo.common.idx, endo.sig.idx)

pdf("up-tf-endo_corr_endo-1_vs_endo-2.pdf", width = 3, height = 3)
PlotHexBin(endo.cfs.1[endo.hits.idx], endo.cfs.2[endo.hits.idx], n.bins = 25)
dev.off()



#### Correlate Neuronal Screen batches ####
library(perturbLM)
library(swne)

neuron.coefs.df.1 <- read.table("up-tf-neuron-original_regression.pvals.tsv", 
                                sep = "\t", stringsAsFactors = F)
neuron.coefs.df.1 <- subset(neuron.coefs.df.1, !grepl(":", Group))
neuron.cfs.1 <- neuron.coefs.df.1$cf
rownames(neuron.coefs.df.1) <- names(neuron.cfs.1) <- interaction(neuron.coefs.df.1$Gene, neuron.coefs.df.1$Group, sep = "_")

neuron.coefs.df.2 <- read.table("up-tf-neuron-nohygro_regression.pvals.tsv", 
                                sep = "\t", stringsAsFactors = F)
neuron.coefs.df.2 <- subset(neuron.coefs.df.2, !grepl(":", Group))
neuron.cfs.2 <- neuron.coefs.df.2$cf
rownames(neuron.coefs.df.2) <- names(neuron.cfs.2) <- interaction(neuron.coefs.df.2$Gene, neuron.coefs.df.2$Group, sep = "_")

max.fdr <- 0.5
neuron.common.idx <- intersect(names(neuron.cfs.1), names(neuron.cfs.2))
neuron.sig.idx <- union(rownames(subset(neuron.coefs.df.1, FDR < max.fdr)), 
                        rownames(subset(neuron.coefs.df.2, FDR < max.fdr)))
neuron.hits.idx <- intersect(neuron.common.idx, neuron.sig.idx)

pdf("neuron_shuffling_single_TF_correlation.pdf", width = 3.5, height = 3.5)
PlotHexBin(neuron.cfs.1[neuron.hits.idx], neuron.cfs.2[neuron.hits.idx], n.bins = 30)
dev.off()