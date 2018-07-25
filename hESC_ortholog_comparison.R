#### Orthologous mouse TF enrichment in hPSC media. Generates Figure S3d
library(perturbLM)
library(swne)
library(cellMapper)

## Load regression results in stem cell media
coefs.df <- read.table("up-tf-stem_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)

#### Correlate coefficients with bulk microarray fold-change ####
bulk.tf.expr <- as.matrix(read.table("bulk_mouse_tf_screen_logfc.txt", sep = "\t", header = T, row.names = 1))

mouse2human <- MouseHumanMapping(rownames(bulk.tf.expr))
bulk.tf.expr <- bulk.tf.expr[names(mouse2human),]
rownames(bulk.tf.expr) <- mouse2human[rownames(bulk.tf.expr)]
dim(bulk.tf.expr)

coefs <- UnflattenDataframe(coefs.df, "cf")
genes.use <- intersect(rownames(coefs), rownames(bulk.tf.expr))
coefs <- coefs[genes.use,]; bulk.tf.expr <- bulk.tf.expr[genes.use,];

mouse.tf.mapping <- MouseHumanMapping(colnames(bulk.tf.expr))
bulk.tf.expr <- bulk.tf.expr[,names(mouse.tf.mapping)]
colnames(bulk.tf.expr) <- mouse.tf.mapping
# tfs.use <- intersect(colnames(bulk.tf.expr), colnames(coefs))
tfs.use <- c("CDX2", "KLF4", "MYC", "MYOD1", "ASCL1")

top.coefs.df <- subset(coefs.df, FDR < 0.1 & abs(cf) > 0.025)
sig.genes <- intersect(top.coefs.df$Gene, genes.use)

tfs.corr <- matrix(0, length(tfs.use), length(tfs.use)); rownames(tfs.corr) <- colnames(tfs.corr) <- tfs.use;
for(tf1 in tfs.use) {
  for(tf2 in tfs.use) {
    r <- cor(coefs[sig.genes,tf1], bulk.tf.expr[sig.genes,tf2], method = "spearman")
    tfs.corr[tf1, tf2] <- r
    tfs.corr[tf2, tf1] <- r
  }
}

# pdf("up-tf-stem_mouse_tf_cor.pdf", width = 4.5, height = 4)
ggHeat(tfs.corr, x.lab.size = 14, y.lab.size = 14)
# dev.off()



#### GSEA enrichment on ortholog target genesets ####
## Create mouse TF targets genesets
mouse.targets.df <- read.table("mouse_tf_targets.txt", sep = "\t", header = T, stringsAsFactors = F)
genotypes.keep <- c("Ascl1", "Cdx2", "Klf4", "Otx2", "Sox2", "Mef2c", "Myc", "Myod1", "Otx2")

mouse.targets.df <- subset(mouse.targets.df, TF %in% genotypes.keep & log.fc > 0.2 & FDR < 0.05)
mouse.targets.df <- mouse.targets.df[order(mouse.targets.df$FDR),]
mouse.targets.df <- Reduce(rbind, by(mouse.targets.df, mouse.targets.df$TF, head, n = 600))
print(table(mouse.targets.df$TF))

mouse.genesets <- lapply(unique(mouse.targets.df$TF), function(tf) subset(mouse.targets.df, TF == tf)$Genes)
names(mouse.genesets) <- sapply(unique(mouse.targets.df$TF), toupper)

all.mouse.genes <- unique(unlist(mouse.genesets, F, F))
mouse2human <- MouseHumanMapping(all.mouse.genes)

human.genesets <- lapply(mouse.genesets, function(mousex) {
  mousex <- intersect(mousex, names(mouse2human))
  humanx <- mouse2human[mousex]
  intersect(humanx, unique(coefs.df$Gene))
})
print(sapply(human.genesets, length))


## GSEA enrichment for mouse ortholog genesets
n.cores <- 24

coefs.list <- lapply(names(human.genesets), function(g) {
  df <- subset(coefs.df, Group == g)
  cf <- df$cf; names(cf) <- df$Gene;
  cf
})
names(coefs.list) <- names(human.genesets);

genesets.enrich <- MultipleGSEAEnrich(coefs.list, human.genesets, n.rand = 1e5, n.cores = n.cores, power = 1)

heat.df <- genesets.enrich
heat.df$lp <- mapply(function(p, cf) -log10(p) * sign(cf), heat.df$p.val, heat.df$sscore)
heat.lp <- UnflattenDataframe(heat.df, output.name = "lp", row.col = 'genesets', col.col = 'Group')
heat.lp[is.na(heat.lp)] <- 0

# pdf("up-tf-stem_bulk_ortho_gsea_enrich.pdf", width = 5.5, height = 2.75)
ggHeat(heat.lp, clustering = "both", x.lab.size = 12, y.lab.size = 12)
# dev.off()