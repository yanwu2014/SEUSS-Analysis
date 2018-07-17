#### Orthologous mouse TF enrichment in hPSC media. Generates Figure S3d
library(perturbLM)
library(swne)
library(cellMapper)

## Load regression results in stem cell media
coefs.df <- read.table("up-tf-stem_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
genotype.summary <- read.table("up-tf-stem_summary.tsv", header = T, sep = "\t", stringsAsFactors = F)

min.cl.pval <- 1e-12
min.diff.genes <- 40
ctrl.cl.stem <- genotype.summary["mCherry", "min_cluster"]
sig.genotypes <- rownames(subset(genotype.summary, (abs(min_cluster_pval) < min.cl.pval & min_cluster != ctrl.cl.stem) 
                                 | n_diff_genes >= min.diff.genes))

## Create mouse TF targets genesets
mouse.targets.df <- read.table("mouse_tf_targets.txt", sep = "\t", header = T, stringsAsFactors = F)
genotypes.keep <- c("Ascl1", "Cdx2", "Klf4", "Otx2", "Myc", "Myod1", "Otx2")

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

coefs.list <- lapply(sig.genotypes, function(g) {
  df <- subset(coefs.df, Group == g)
  cf <- df$cf; names(cf) <- df$Gene;
  cf
})
names(coefs.list) <- sig.genotypes;

genesets.enrich <- MultipleGSEAEnrich(coefs.list, human.genesets, n.rand = 50000, n.cores = n.cores, power = 1)
# genesets.enrich <- read.table("up-tf-stem_bulk_ortho_gsea_enrich.tsv", sep = "\t", header = T)

heat.df <- genesets.enrich
heat.df$lp <- mapply(function(p, cf) -log(p) * sign(cf), heat.df$p.val, heat.df$sscore)
heat.lp <- UnflattenDataframe(heat.df, output.name = "lp", row.col = 'genesets', col.col = 'Group')
heat.lp[is.na(heat.lp)] <- 0

pdf("up-tf-stem_bulk_ortho_gsea_enrich.pdf", width = 5.5, height = 2.75)
ggHeat(heat.lp, clustering = "col", x.lab.size = 12, y.lab.size = 12)
dev.off()

write.table(genesets.enrich, file = "up-tf-stem_bulk_ortho_gsea_enrich.tsv", sep = "\t")
