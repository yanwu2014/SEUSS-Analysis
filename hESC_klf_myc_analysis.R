#### Split KLF family + MYC mutants counts matrix ####
library(Seurat)
library(perturbLM)

counts <- ReadData("up-tf-klf-myc/")
genotypes.list <- ReadGenotypes("up-tf-klf-myc_pheno_dict.csv")
length(unique(unlist(genotypes.list, F, F)))

genotypes.matrix <- as.matrix(DesignMatrixGenotypes(genotypes.list, max.genotypes = 2, 
                                                    min.cells = 10, drop.cells = F))
dim(genotypes.matrix)

klf.genotypes <- grep("KLF", names(genotypes.list), value = T)
myc.genotypes <- grep("MYC|Myc", names(genotypes.list), value = T)

klf.genotypes.matrix <- genotypes.matrix[rowSums(genotypes.matrix[,myc.genotypes]) == 0,]
myc.genotypes.matrix <- genotypes.matrix[rowSums(genotypes.matrix[,klf.genotypes]) == 0,]
dim(klf.genotypes.matrix)
dim(myc.genotypes.matrix)

klf.cells <- rownames(klf.genotypes.matrix)
myc.cells <- rownames(myc.genotypes.matrix)
klf.counts <- as.matrix(counts[,klf.cells])
myc.counts <- as.matrix(counts[,myc.cells])

dim(klf.counts)
dim(myc.counts)

# write.table(klf.counts, file = "up-tf-klf.counts.tsv", sep = "\t")
# write.table(myc.counts, file = "up-tf-myc.counts.tsv", sep = "\t")


#### MYC/KLF Gene Module analysis ####
library(perturbLM)
library(swne)

## Load data
coefs.df.myc <- read.table("up-tf-myc_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.df.myc <- subset(coefs.df.myc, Group != "metadata")
top.hits.myc <- subset(coefs.df.myc, FDR < 0.05 & abs(cf) > 0.01)

coefs.df.klf <- read.table("up-tf-klf_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.df.klf <- subset(coefs.df.klf, Group != "metadata")
top.hits.klf <- subset(coefs.df.klf, FDR < 0.05 & abs(cf) > 0.01)

## TF effects on gene modules
coefs.matrix.myc <- UnflattenDataframe(coefs.df.myc, "cf")
coefs.matrix.klf <- UnflattenDataframe(coefs.df.klf, "cf")

gene.module.mapping <- read.table("up-tf_module_names.txt", sep = "\t", header = T)
colnames(gene.module.mapping) <- c("Module", "Description")
gene.modules.list <- LoadGenesets("up-tf_gene_modules.gmt")

tf.modules.matrix.myc <- CalcGeneModuleEffect(coefs.matrix.myc, gene.modules.list, 
                                              gene.module.mapping, min.coef = 0.025)
tf.modules.matrix.klf <- CalcGeneModuleEffect(coefs.matrix.klf, gene.modules.list, 
                                              gene.module.mapping, min.coef = 0.025)

pdf("up-tf-myc_gene_module_effects.pdf", width = 7.25, height = 4.0)
ggHeat(tf.modules.matrix.myc, clustering = "both", x.lab.size = 12, y.lab.size = 12)
dev.off()

pdf("up-tf-klf_gene_module_effects.pdf", width = 7.25, height = 4.0)
ggHeat(tf.modules.matrix.klf, clustering = "both", x.lab.size = 12, y.lab.size = 12)
dev.off()

#### Correlate MYC and KLF4 with hPSC media screen

## Load hPSC media screen results
coefs.df.stem <- read.table("up-tf-stem_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)

## Correlate MYC between myc mutant screen and hPSC media screen
TF <- "MYC"
g.coefs.df <- subset(coefs.df.stem, Group == TF)
g.coefs.df.myc <- subset(coefs.df.myc, Group == TF)

cfs <- g.coefs.df$cf; names(cfs) <- g.coefs.df$Gene;
cfs.myc <- g.coefs.df.myc$cf
names(cfs.myc) <- g.coefs.df.myc$Gene
genes.use <- intersect(names(cfs), names(cfs.myc))

pdf("up-tf-myc_hPSC-myc_correlation.pdf", width = 4.5, height = 4)
PlotCorrelation(cfs[genes.use], cfs.myc[genes.use], x.lab = "hPSC screen", y.lab = "MYC mutants", title = TF, 
                use.label = F, box = T, show.corr = T)
dev.off()

## Correlate KLF4 between KLF family screen and hPSC media screen
TF <- "KLF4"
g.coefs.df <- subset(coefs.df.stem, Group == TF)
g.coefs.df.klf <- subset(coefs.df.klf, Group == TF)

cfs <- g.coefs.df$cf; names(cfs) <- g.coefs.df$Gene;
cfs.klf <- g.coefs.df.klf$cf
names(cfs.klf) <- g.coefs.df.klf$Gene
genes.use <- intersect(names(cfs), names(cfs.klf))

pdf("up-tf-klf_hPSC-klf_correlation.pdf", width = 4.5, height = 4)
PlotCorrelation(cfs[genes.use], cfs.klf[genes.use], x.lab = "hPSC screen", y.lab = "klf mutants", title = TF, 
                use.label = F, box = T, show.corr = T)
dev.off()
