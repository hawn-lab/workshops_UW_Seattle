#### Install packages ####
install.packages("msigdbr")
install.packages("fgsea")

install.packages("BiocManager")
BiocManager::install("clusterProfiler")

#### Load packages ####
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)

#Check MSigDB version
packageVersion("msigdbr")

#### Load data ####
load("data/RSTR_data_clean_subset.RData")

model.results <- read_csv("data/RSTR.Mtb.model.subset.csv")

#### Enrichment ####
#Get gene set database
H <- msigdbr(species = "Homo sapiens", category = "H")

class(H)

#Define significant genes
##Look at FDR spread
ggplot(model.results, aes(x=FDR)) +
  geom_histogram(bins=100) +
  theme_classic()

table(model.results$FDR == 0)

#Set FDR cutoff and get genes 
signif <- model.results %>% 
  filter(FDR <= 1E-16)

#get entrez ID
signif.entrez <- unique(signif$entrezgene_id)
H.entrez <- select(H, gs_name, entrez_gene)

#run enrichment
enrich.H <- enricher(gene = signif.entrez, TERM2GENE = H.entrez)

#or for ensembl
#get ensembl ID
signif.ensembl <- unique(signif$ensembl_gene_id)
H.ensembl <- select(H, gs_name, ensembl_gene)

#run enrichment
enrich.H <- enricher(gene = signif.ensembl, TERM2GENE = H.ensembl)

#### Extract results ####
class(enrich.H)

head(enrich.H@result)

class(enrich.H@result$GeneRatio)

#format results
enrich.H.df <- enrich.H@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)

#### Visualize results ####
enrich.H.df %>% 
  filter(p.adjust <= 0.05) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = gsub("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Mtb significant genes (FDR < 1E-16)\nenriched in Hallmark gene sets (FDR < 0.05)")



########################

#### GSEA ####
#format gene set database
head(H)

H.ensembl.ls <- H %>% 
  select(gs_name, ensembl_gene) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(ensembl_gene))) %>% 
  deframe()

# Calculate fold change 
#Extract expression data
FC <- as.data.frame(dat$E) %>% 
  #Move gene IDs from rownames to a column
  rownames_to_column("ensembl_gene_id") %>% 
  #Make long format with all expression values in 1 column
  pivot_longer(-ensembl_gene_id, 
               names_to = "libID", values_to = "expression") %>% 
  #Extract RSID and TB condition from libID
  #If this info was not in the libID, we could get it by joining
  # with dat$targets
  separate(libID, into = c("RSID","condition"), sep="_") %>% 
  #Make wide with media and tb expression in separate columns
  pivot_wider(names_from = condition, values_from = expression) %>% 
  #Calculate tb minus media fold change (delta for change)
  #Because expression values are log2, subtraction is the same as division
  mutate(delta = TB-MEDIA) %>% 
  #Calculate mean fold change per gene
  group_by(ensembl_gene_id) %>% 
  summarise(mean.delta = mean(delta, na.rm=TRUE)) %>% 
  #Arrange by descending fold change
  arrange(desc(mean.delta))

##format for gsea 
FC.vec <- FC$mean.delta
names(FC.vec) <- FC$ensembl_gene_id

#set score type
min(FC.vec)
max(FC.vec)
scoreType <- "std"

#Run GSEA
gsea.H <- fgseaSimple(pathways = H.ensembl.ls,
                      stats = FC.vec,
                      scoreType = scoreType,
                      nperm=1000)

#### plot gsea results ####
class(gsea.H)

gsea.H %>% 
  filter(padj <= 0.05) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(pathway = gsub("HALLMARK_","", pathway),
         pathway = gsub("_"," ", pathway)) %>% 
  
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y=NES)) +
  geom_col() +
  theme_classic() +
  #Force equal max min
  lims(y=c(-3.2,3.2)) +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Normalized enrichment score (NES)",
       x="Gene set",
       title = "Hallmark GSEA (FDR < 0.05)\nDown in +Mtb <--         --> Up in +Mtb")
