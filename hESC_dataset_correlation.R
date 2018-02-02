#### Correlate screens to establish reproducibility ####

library(perturbLM)

## Correlate individual hPSC media samples with merged sample
coefs.df <- read.table("up-tf-all_regression_trimmed.pvals.tsv", sep = "\t", stringsAsFactors = F)
cfs <- coefs.df$cf; names(cfs) <- interaction(coefs.df$Gene, coefs.df$Group, sep = "_");
top.hits <- subset(coefs.df, FDR < 0.05, abs(cf) > 0.025)

individual.files <- c("hESC_Individual_Batch_Coefficients/up-tf1_regression.pvals.tsv",
                      "hESC_Individual_Batch_Coefficients/up-tf2_regression.pvals.tsv",
                      "hESC_Individual_Batch_Coefficients/up-tf3_regression.pvals.tsv",
                      "hESC_Individual_Batch_Coefficients/up-tf4_regression.pvals.tsv",
                      "hESC_Individual_Batch_Coefficients/up-tf5_regression.pvals.tsv")

for (i in 1:length(individual.files)) {
  fi <- individual.files[[i]]
  coefs.df.batch <- read.table(fi, sep = "\t", header = T, stringsAsFactors = F)
  cfs.batch <- coefs.df.batch$cf
  names(cfs.batch) <- interaction(coefs.df.batch$Gene, coefs.df.batch$Group, sep = "_")
  top.hits.batch <- subset(coefs.df.batch, FDR < 0.05, abs(cf) > 0.025)
  
  hits.use <- unique(union(interaction(top.hits$Gene, top.hits$Group, sep = "_"), 
                           interaction(top.hits.batch$Gene, top.hits.batch$Group, sep = "_")))
  hits.use <- hits.use[hits.use %in% names(cfs)]; hits.use <- hits.use[hits.use %in% names(cfs.batch)]
  f.name <- paste("up_tf", i, sep = "")
  # pdf(paste(f.name, "correlation.pdf", sep = "_"), width = 3, height = 3.25)
  print(gghexbin(cfs[hits.use], cfs.batch[hits.use], n.bins = 100))
  # dev.off()
}


## Correlate multilineage media batches

coefs.df.1 <- read.table("hESC_Individual_Batch_Coefficients/up-tf-multi-batch1_regression.pvals.tsv", sep = "\t", stringsAsFactors = F)
cfs.1 <- coefs.df.1$cf
names(cfs.1) <- interaction(coefs.df.1$Gene, coefs.df.1$Group, sep = "_")

coefs.df.2 <- read.table("hESC_Individual_Batch_Coefficients/up-tf-multi-batch2_regression.pvals.tsv", sep = "\t", stringsAsFactors = F)
cfs.2 <- coefs.df.2$cf
names(cfs.2) <- interaction(coefs.df.2$Gene, coefs.df.2$Group, sep = "_")

top.hits.1 <- subset(coefs.df.1, FDR < 0.05, abs(cf) > 0.025)
top.hits.2 <- subset(coefs.df.2, FDR < 0.05, abs(cf) > 0.025)
hits.use <- unique(union(interaction(top.hits.1$Gene, top.hits.1$Group, sep = "_"), 
                         interaction(top.hits.2$Gene, top.hits.2$Group, sep = "_")))
hits.use <- hits.use[hits.use %in% names(cfs.1)]; hits.use <- hits.use[hits.use %in% names(cfs.2)]

# pdf("up-tf-multi_correlation.pdf", width = 3, height = 3.25)
gghexbin(cfs.1[hits.use], cfs.2[hits.use], n.bins = 100)
# dev.off()