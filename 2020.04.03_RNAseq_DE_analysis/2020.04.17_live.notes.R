# Differential expression analysis

###### Setup ######
#Load packages
library(tidyverse)
library(limma)

#Load data
load("data/AM.clean.RData")

#Set seed
set.seed(927)

##### DE genes #####
# Effects of TB infection on AM gene expression

# limma
# 1. Make a model
# 2. Group by donor (if needed)
# 3. Fit model
# 4. Estimate P-val and effect size

# Make model
model <- model.matrix(~TB, data = dat.AM.pc.abund.norm.voom$targets)
  colnames(model) <- c("(Intercept)", "TB")

# Group by donor
# similar to repeated measure
# Calculate how correlated samples from the same donor are
consensus.corr <- duplicateCorrelation(
  object = dat.AM.pc.abund.norm.voom$E,
  block = dat.AM.pc.abund.norm.voom$targets$ptID,
  design = model)

consensus.corr$consensus.correlation

# Fit model with limma
fit <- lmFit(
  object = dat.AM.pc.abund.norm.voom$E,
  design = model,
  block = dat.AM.pc.abund.norm.voom$targets$ptID,
  correlation = consensus.corr$consensus.correlation)

# Empirical Bayes and P-Val estimation
efit <- eBayes(fit)

# Extract P-val
# By-hand
pval <- topTable(fit = efit,
                      coef="TB",
                      number=nrow(dat.AM.pc.abund.norm.voom),
                      adjust.method = "BH")

# Kim's function
source("https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/limma.extract.pval.R")

extract.pval(model = model,
             voom.dat = dat.AM.pc.abund.norm.voom$E,
             eFit = efit,
             name = "gene_pval",
             summary = TRUE)

#results in
gene_pval
gene_pval.summ

##### DE gene plots #####
#volcano plot
gene_pval %>% 
  filter(group == "TB") %>% 
  mutate(signif.group = ifelse(adj.P.Val <= 0.05 & 
                                 FC.group == "up",
                               "up",
                        ifelse(adj.P.Val <= 0.05 & 
                                 FC.group == "down",
                               "down",
                               "NS"))) %>% 
  
ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = signif.group)) +
  theme_classic() +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_hline(yintercept = -log10(0.05)) +
  lims(x = c(-10.5,10.5)) +
  labs(x = "Log2 fold change", y = "-log10( FDR )",
       color = "Significance")

# gene boxplots
# see https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/RNAseq_boxplot_fxn.R

##### Supervised WGCNA #####
# weighted correlation network analysis
# Groups genes based on correlation of expression
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559

# List signif genes
genes.signif <- gene_pval %>% 
  filter(adj.P.Val <= 0.05 & group == "TB") %>% 
  select(geneName) %>% unlist(use.names = FALSE)

# Run clustering
# use blockwiseModules( )
# Actually use Kim's function
source("https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/RNAseq_module_fxn.R")

library(WGCNA)
dir.create(path = "figs", showWarnings = FALSE)
dir.create(path = "results", showWarnings = FALSE)

make.modules(voom.dat = dat.AM.pc.abund.norm.voom,
             genes.signif = genes.signif,
             Rsq.min = 0.8,
             minModuleSize = 50,
             deepSplit = 3,
             nThread = 3,
             basename = "TBgenes")

mods.net <- read_csv("results/module_TBgenes_deepSplit3_minMod50/TBgenes_genes_in_mod.csv")

#MEAN counts of genes in modules, NOT sum
mods.voom <- read_csv("results/module_TBgenes_deepSplit3_minMod50/TBgenes_mod_voom_counts.csv")

##### Summarize modules #####
mods.net %>% 
  count(module.char)

#How many genes in module 00?
#How many modules?

#Exercise: Make PCA of module expression

##### Linear model of modules #####
# Make model
# Use same as for genes

# Group by donor
# similar to repeated measure
# Calculate how correlated samples from the same donor are
consensus.corr.mods <- duplicateCorrelation(
  object = column_to_rownames(mods.voom, "module"),
  block = dat.AM.pc.abund.norm.voom$targets$ptID,
  design = model)

consensus.corr.mods$consensus.correlation

# Fit model with limma
fit <- lmFit(
  object = column_to_rownames(mods.voom, "module"),
  design = model,
  block = dat.AM.pc.abund.norm.voom$targets$ptID,
  correlation = consensus.corr.mods$consensus.correlation)

# Empirical Bayes and P-Val estimation
efit <- eBayes(fit)

# Extract P-val
extract.pval(model = model,
             voom.dat = column_to_rownames(mods.voom, "module"),
             eFit = efit,
             name = "module_pval",
             summary = TRUE)
#results
module_pval
module_pval.summ

##### Plot modules #####
#volcano plot
#Boxplots of expression