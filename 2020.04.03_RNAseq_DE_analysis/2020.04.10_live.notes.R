# Quality filtering genes

##### Setup #####

#Install Bioconductor packages
#install.packages("BiocManager")
#BiocManager::install("edgeR")
#BiocManager::install("limma")

#Install GitHub packages
#install.packages("remotes")
#remotes::install_github("dgrtwo/drlib")

#Load packages
library(tidyverse)
library(edgeR)
library(limma)
#library(drlib)

#Load functions
# Recall %in% meaning "in" like is a in b?
# Create the reverse with Negate as "not in"
`%notin%` <- Negate(`%in%`)

#Set seed for random number generator
set.seed(927)

##### Load data #####
#Cleaning metrics
meta <- read_csv("data/AM.MDM.data.cleaning.metrics.csv")
#Counts of reads in genes(exons only)
counts <- read_csv("data/AM.MDM.counts.paired.clean.csv")
#Sample metadata extracted from sample ID in metrics file
samp <- meta %>% 
  select(sampID) %>% 
  #Separate sampID into multiple columns
  separate(col = sampID, into = c("ptID", "cell", "TB"), sep="_",
           #Keep original sampID column
           remove = FALSE) %>% 
  #Modify ptID to remove leading "AM" text
  mutate(ptID = gsub("AM", "", ptID))

##### Separate cell types #####
# Filter sample metadata
#List samples IDs in groups of interest
AM.vec <- samp %>% 
  filter(cell == "AM") %>% 
  select(sampID) %>% 
  unlist(use.names = FALSE)

MDM.vec <- samp %>% 
  filter(cell == "MDM") %>% 
  select(sampID) %>% 
  unlist(use.names = FALSE)

#c(samp$sampID)

#Select groups from count data
counts.AM <- counts %>% 
  select(geneName, all_of(AM.vec))
  
counts.MDM <- counts %>% 
  select(geneName, all_of(MDM.vec))
  
dim(counts.AM)
dim(counts.MDM)
  
##### QC of genes - AM samples ######
##### Filter to protein coding genes ###### 
# Load gene key
# Filtering to protein coding
key <- read_tsv("data/EnsemblToHGNC_GRCh38.txt",
                na=c(NA, "", ".", "N/A")) %>% 
  filter(gene_biotype == "protein_coding" &
         !is.na(hgnc_symbol)) %>% 
  #remove duplicates if exist
  distinct() %>% 
  #Sort by name
  rename(geneName = ensembl_gene_id) %>% 
  arrange(geneName) %>% 
  #Combine duplicate annotations
  group_by(geneName, gene_biotype) %>% 
  summarize(hgnc_symbol = paste(hgnc_symbol, collapse="_"))

#Filter count table to protein coding
counts.AM.pc <- counts.AM %>% 
  filter(geneName %in% key$geneName)

##### Make DGEList object ##### 
#What we're making 
counts.AM.pc %>% 
  column_to_rownames("geneName") %>% 
  as.matrix()

as.matrix(column_to_rownames(counts.AM.pc, "geneName"))


dat.AM.pc <- DGEList(
  counts = as.matrix(column_to_rownames(counts.AM.pc, "geneName")),
  samples = as.matrix(filter(samp, cell == "AM")),
  genes = as.matrix(filter(key, geneName %in% counts.AM.pc$geneName))
)

##### Filter rare genes #####
# Plot distribution of raw genes
temp <- voomWithQualityWeights(
  counts = dat.AM.pc,
  design = model.matrix(~TB, data = dat.AM.pc$samples),
  plot=FALSE, save.plot=TRUE
)

data.frame(
  x = temp$voom.xy$x,
  y = temp$voom.xy$y,
  linex = temp$voom.line$x,
  liney = temp$voom.line$y) %>% 
  
  ggplot() +
  geom_point(aes(x=x, y=y), size=0.5) +
  geom_path(aes(x=linex, y=liney), color="red") +
  theme_classic() +
  labs(x="log2( count size + 0.5 )", y="Sqrt (stdev)",
       title="voomQW: Mean-variance trend")

# Built in edgeR filter
#filterByExpr()

# Kim's function
source("https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/RNAseq_rare_gene_filter.R")

rare.gene.filter(dat = dat.AM.pc,
                 min.pct = 20,
                 min.CPM = 1,
                 name = "dat.AM.pc.abund")

#Re-assess plot
temp <- voomWithQualityWeights(
  counts = dat.AM.pc.abund,
  design = model.matrix(~TB, data = dat.AM.pc.abund$samples),
  plot=FALSE, save.plot=TRUE
)

data.frame(
  x = temp$voom.xy$x,
  y = temp$voom.xy$y,
  linex = temp$voom.line$x,
  liney = temp$voom.line$y) %>% 
  
  ggplot() +
  geom_point(aes(x=x, y=y), size=0.5) +
  geom_path(aes(x=linex, y=liney), color="red") +
  theme_classic() +
  labs(x="log2( count size + 0.5 )", y="Sqrt (stdev)",
       title="voomQW: Mean-variance trend")

#How genes were removed?
nrow(dat.AM.pc$genes)
nrow(dat.AM.pc.abund$genes)

##### RNA composition normalization #####
# Takes in account relative abundance

dat.AM.pc.abund.norm <- calcNormFactors(dat.AM.pc.abund)

##### voom normalize #####
# log2 counts per million (CPM)
# quality weights option

dat.AM.pc.abund.norm.voomQW <- voomWithQualityWeights(
  counts = dat.AM.pc.abund.norm,
  design = model.matrix(~TB, 
                        data = dat.AM.pc.abund.norm$samples),
  plot=TRUE)

##### PCA #####
PCA.dat <- as.data.frame(dat.AM.pc.abund.norm.voomQW$E) %>% 
  t() %>% 
  prcomp()

PC1.label <- paste("PC1 (", 
                   summary(PCA.dat)$importance[2,1]*100, 
                   "%)", sep="")
PC2.label <-paste("PC2 (", 
                  summary(PCA.dat)$importance[2,2]*100, 
                  "%)", sep="")
#color by TB status
as.data.frame(PCA.dat$x) %>% 
  rownames_to_column("sampID") %>%
  left_join(as.data.frame(dat.AM.pc.abund.norm.voomQW$targets)) %>% 
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=TB), size=3) +
  theme_classic() +
  labs(x=PC1.label, y=PC2.label)

#color by donor
as.data.frame(PCA.dat$x) %>% 
  rownames_to_column("sampID") %>%
  left_join(as.data.frame(dat.AM.pc.abund.norm.voomQW$targets),
            by="sampID") %>% 
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=ptID), size=3) +
  theme_classic() +
  labs(x=PC1.label, y=PC2.label)
