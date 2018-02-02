#### GATA4/MITF media condition specific analysis ####

## Load data
library(perturbLM)

coefs.stem <- read.table("up-tf-all_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.endo <- read.table("up-tf-endo_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.multi <- read.table("up-tf-multi_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)

gene.modules.list <- LoadGenesets("up-tf_functional_gene_modules.gmt")
gene.modules.mapping <- read.table("up-tf_gene_module_mapping.txt", sep = "\t", header = T, stringsAsFactors = F)
gene.modules.names <- gene.modules.mapping$Description; names(gene.modules.names) <- gene.modules.mapping$Module;
names(gene.modules.list) <- sapply(names(gene.modules.list), function(x) gene.modules.names[[x]])
gene.modules <- FlattenGenotypeList(gene.modules.list)

top.coefs.mitf <- subset(coefs.multi, Group == "MITF" & abs(cf) > 0.05 & FDR < 0.1)
top.coefs.mitf <- top.coefs.mitf[order(top.coefs.mitf$cf, decreasing = T),]
top.coefs.mitf$Module <- gene.modules[top.coefs.mitf$Gene]
table(top.coefs.mitf$Module)

top.coefs.gata4 <- subset(coefs.endo, Group == "GATA4" & abs(cf) > 0.05 & FDR < 0.1)
top.coefs.gata4 <- top.coefs.gata4[order(top.coefs.gata4$cf, decreasing = T),]
top.coefs.gata4$Module <- gene.modules[top.coefs.gata4$Gene]
table(top.coefs.gata4$Module)

mitf.coefs <- subset(coefs.multi, Group == "MITF" & Gene %in% names(gene.modules))$cf
names(mitf.coefs) <- subset(coefs.multi, Group == "MITF" & Gene %in% names(gene.modules))$Gene
mitf.module.effects <- tapply(mitf.coefs, factor(gene.modules[names(mitf.coefs)]), mean)
print(sort(mitf.module.effects, decreasing = T))

gata4.coefs <- subset(coefs.endo, Group == "GATA4" & Gene %in% names(gene.modules))$cf
names(gata4.coefs) <- subset(coefs.endo, Group == "GATA4" & Gene %in% names(gene.modules))$Gene
gata4.module.effects <- tapply(gata4.coefs, factor(gene.modules[names(gata4.coefs)]), mean)
print(sort(gata4.module.effects, decreasing = T))

library(liger)
genesets <- LoadGenesets("hallmark_go_bp_genesets.gmt")

genesets.use <- FilterGenesets(genesets, names(mitf.coefs), min.size = 20, max.size = 400)
mitf.enrich.genesets <- bulk.gsea(mitf.coefs, genesets.use, power = 1, n.rand = 10000, mc.cores = 16)
mitf.enrich.genesets <- subset(mitf.enrich.genesets[order(mitf.enrich.genesets$sscore, decreasing = T),], q.val < 0.05)

genesets.use <- FilterGenesets(genesets, names(gata4.coefs), min.size = 20, max.size = 400)
gata4.enrich.genesets <- bulk.gsea(gata4.coefs, genesets.use, power = 1, n.rand = 10000, mc.cores = 16)
gata4.enrich.genesets <- subset(gata4.enrich.genesets[order(gata4.enrich.genesets$sscore, decreasing = T),], q.val < 0.05)

# write.table(top.coefs.mitf, file = "MITF_top_hits.tsv", sep = "\t")
# write.table(mitf.module.effects, file = "MITF_module_effects.tsv", sep = "\t")
# write.table(mitf.enrich.genesets, file = "MITF_geneset_enrichment.tsv", sep = "\t")

# write.table(top.coefs.gata4, file = "GATA4_top_hits.tsv", sep = "\t")
# write.table(gata4.module.effects, file = "GATA4_module_effects.tsv", sep = "\t")
# write.table(gata4.enrich.genesets, file = "GATA4_geneset_enrichment.tsv", sep = "\t")