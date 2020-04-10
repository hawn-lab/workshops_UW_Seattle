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
library(drlib)

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
samp <- sampID %>% 
  select(meta) %>% 
  #Separate sampID into multiple columns
  separate(col = sampID, into = c("ptID", "cell", "TB"), sep="_",
           #Keep original sampID column
           remove = FALSE) %>% 
  #Modify ptID to remove leading "AM" text
  mutate(ptID = gsub("AM", "", ptID))

#gsub uses regular expressions aka regex
#for example, this is to replace x## where the # are any single digit
#mutate(ptID = gsub("x[0-9][0-9]", "", ptID))

##### Data cleaning #####
# See the following for some helpful scripts
# https://github.com/kdillmcfarland/R_bioinformatic_scripts

##### Median CV coverage vs. alignment percentage #####
# Median coefficient of variation (CV) coverage is the variance in coverage across genes in a sample
# You want median CV coverage to be low, indicating similar coverage of all genes in a sample

# You want high alignment % to indicate that the reads in a sample were high-quality and successfully aligned to the genome

# Remember that
# ggplot(meta) is the same as
# meta %>% ggplot()

#Plot the two metrics against each other
meta %>% 
  ggplot(aes(x=MEDIAN_CV_COVERAGE, y=PCT_PF_ALIGNED)) +
  geom_point(size=2) +
  theme_classic() +
  #Force axes to go from 0 to 1
  lims(x=c(0,1), y=c(0,1)) +
  #Add horizontal/vertical lines at desired quality cutoffs
  geom_hline(yintercept = 0.9) +
  geom_vline(xintercept = 0.9)

#### Total final sequences ####
# Also look at the total number of sequences per sample

#Unordered option
meta %>% 
  separate(col = sampID, into = c("ptID", "cell", "TB"), sep="_",
           remove = FALSE) %>% 
  
  ggplot(aes(x=sampID, 
             y=both.align.paired_filter.paired)) +
  geom_col() +
  theme_classic()

#Order by min to max sequences in each facte
meta %>% 
  #Create samp variables for coloring plot
  separate(col = sampID, into = c("ptID", "cell", "TB"), sep="_",
           remove = FALSE) %>% 
  
  ggplot(aes(
    #Order sampID (x axis) by the total seqs within each cell type
    x=reorder_within(sampID,
                     by=both.align.paired_filter.paired,
                     within=cell), 
    #Set y variable
    y=both.align.paired_filter.paired,
    #Fill bars by TB status
    fill=TB)) +
  #Column plot
  geom_col() +
  theme_classic() +
  # Facet to separate cell types into 2 adjoining plots
  # format y ~ x, columns ~ rows
  facet_wrap(~cell, scales = "free_x") +
  geom_hline(yintercept = 500000)

#### Exercise #####
# Continue to improve the above plot such as
# roating the x-axis labels so they are readable
# changing the axes names
# anything else you want to change!

##### PCA outliers #####
# An outlier is generally defined as a sample that is > 3 std deviations away from the group mean on any PC axis. However, this is not a hard rule and generally, you can look at the PCA and easily ID potential outliers

PCA <- counts %>% 
  #Put geneName column into rownames b/c prcomp does not allow non-numeric columns
  column_to_rownames("geneName") %>% 
  #Convert to log2 counts per million
  cpm(log=TRUE) %>% 
  #transpose table
  t() %>% 
  #Calculate the PCA
  prcomp()

#Extract PC values to use in plot 
#For PCA$x to be a data frame so it can work with tidyverse functions
PCA.dat <- as.data.frame(PCA$x) %>% 
  #Take rownames and move to a data column named sampID
  rownames_to_column("sampID") %>% 
  #Merge with sample metadata, matching rows based on sampID
  full_join(samp, by="sampID")

#Plot
PCA.dat %>% 
  #Create color variable to show both cell type and TB status
  mutate(color.var = paste(cell, TB, sep="_")) %>% 
  
  ggplot(aes(x=PC1,y=PC2, color=color.var)) +
  geom_point(size=2) +
  theme_classic()

#Forgot to show this in the tutorial
#Extract the % variation explained by each axis
summary(PCA)$importance

#Use above info to label PCA
PCA.dat %>% 
  #Create color variable to show both cell type and TB status
  mutate(color.var = paste(cell, TB, sep="_")) %>% 
  
  ggplot(aes(x=PC1,y=PC2, color=color.var)) +
  geom_point(size=2) +
  theme_classic() +
  #Change axes labels
  labs(x = "PC1 (33.839%)", y = "PC2 (14.957%)",
       #Change legend label
       color="")

##### Data cleaning conclusions ####
# These data are high quality and no samples need to be removed.
# However, if this were not the case, we would use the 3 plots to define reasonable quality cutoffs and remove any samples that fail to meet those cutoffs.

##### END #####