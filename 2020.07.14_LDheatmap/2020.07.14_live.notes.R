library(tidyverse)
library(genetics)
library(LDheatmap)

##### Load data #####
# Data from N100 Shiny app
## Gene PRKAG2 intron, exon, promoter. No expansion of region
## All N100 samples
dat <- read_csv("Hawn_RSTR_SNPlist.PRKAG2.csv")

##### Clean data ##### 
# Convert 0/0 to A/A format
dat.allele <- dat %>% 
  #Remove SNP metadata
  dplyr::select(-c(rsID:POS)) %>% 
  #Convert to long format
  pivot_longer(-c(snpID:allele.1), 
               names_to="FULLIDNO", values_to="genotype") %>% 
  ##Convert genotype to alleles
  mutate(geno.allele = ifelse(genotype == "0/0",
                              paste(allele.0, allele.0, sep="/"),
                       ifelse(genotype == "0/1",
                              paste(allele.0, allele.1, sep="/"),
                       ifelse(genotype == "1/1",
                              paste(allele.1, allele.1, sep="/"),
                              NA)))) %>% 
  #Convert to wide format with SNPs as columns, donors as rows
  dplyr::select(-allele.0, -allele.1, -genotype) %>% 
  pivot_wider(names_from = snpID, values_from = geno.allele) %>% 
  #Convert A/A format to offical "genotype" class
  mutate(across(-FULLIDNO, genetics::as.genotype)) %>% 
  #Move donor IDs to rownames
  column_to_rownames("FULLIDNO")

#To make things run faster here, keep only 10 SNPs
dat.allele.sub <- dat.allele[,1:10]

##### Calculate LD R^2 ##### 
# This can take several minutes if there are a lot of SNPs
# LDmeasure = "D" also available
LD.map <- LDheatmap(dat.allele.sub, LDmeasure = "r")

#Extract LD values in a matrix
LD.vals <- LD.map$LDmatrix

##### LD heatmap ##### 
#Get SNP positions
allele.POS <- dat %>% 
  #Keep only SNPs in LD matrix
  filter(snpID %in% colnames(LD.vals))

#Check that order of positions matches order of LD values
## If false reorder both to match
identical(allele.POS$snpID, colnames(LD.vals))

# Create LD heatmap
LD.plot <- LDheatmap(LD.vals, 
            #Set genetic distances
            genetic.distances = allele.POS$POS,
            #Put annotation on top of plot
            flip=TRUE,
            #Color white to red
            color=heat.colors(50))

##### Add genomic annotation ##### 
# This step can also take several minutes if the chromosome of interest is large
LD.plot.annot <- LDheatmap.addGenes(LD.plot, 
                                    #Auto generate chromosome name from data
                                    chr=paste("chr", unique(dat$CHR), sep=""),
                                    #Set reference genome version
                                    genome="hg19")
