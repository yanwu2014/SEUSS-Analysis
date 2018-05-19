library(perturbLM)
library(swne)

orig.coefs.df <- read.table("up-tf-stem_regression_orig.pvals.tsv", sep = "\t", header = T)
orig.coefs <- UnflattenDataframe(orig.coefs.df, "cf")

coefs.df <- read.table("up-tf-stem_regression_lm.pvals.tsv")
coefs <- UnflattenDataframe(coefs.df, "cf")

genes.use <- intersect(rownames(orig.coefs), rownames(coefs))

## Compare coefficient/logFC magnitudes
TF <- "CDX2"
# pdf(paste("hESC", TF, "mdl_cor.pdf", sep = "_"), width = 4, height = 4)
PlotCorrelation(orig.coefs[genes.use, TF], coefs[genes.use, TF], x.lab = "Linear coefficients", y.lab = "LogFC", 
                box = F, title = TF, show.corr = T)
# dev.off()


## Compare statistical significance
top.hits.df <- subset(coefs.df, abs(cf) > 0.05 & FDR < 0.1)
sort(table(top.hits.df$Group), decreasing = T)

orig.top.hits.df <- subset(orig.coefs.df, abs(cf) > 0.05 & FDR < 0.05)
sort(table(orig.top.hits.df$Group), decreasing = T)

top.hits <- as.character(interaction(top.hits.df$Gene, top.hits.df$Group))
orig.top.hits <- as.character(interaction(orig.top.hits.df$Gene, orig.top.hits.df$Group))

n.int <- length(intersect(top.hits, orig.top.hits))
print(n.int/(length(top.hits) + length(orig.top.hits) - n.int))
