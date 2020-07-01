#### Install packages ####
# Only need to do once on your machine so I've commented out these lines from the script.

#Install tidyverse
# install.packages("tidyverse", Ncpus=2)
### If errors, try with 1 Ncpus

#Help function
# ?install.packages

# Install limma (Bioconductor package)
# install.packages("BiocManager")
# :: specifies a function in a specific package
# BiocManager::install("limma")

#### Load packages ####
#Need to do everytime you open R/RStudio
#Good idea to place at the top of your script for convenience
library(tidyverse)
library(limma)

#### Set seed ####
#We don't have any analyses that use this here but it is best practices to set a seed so everything is reproducible
set.seed(4389)

#### Data ####
# Data from tables like csv, tsv, txt
snp <- read.table(file = "data/Hawn_RSTR_SNPlist.PRKAG2.csv", 
                   sep = ",", header = TRUE)

snp <- read.csv(file = "data/Hawn_RSTR_SNPlist.PRKAG2.csv")

# Note that tab separated is sep="\t"

# Data from RData
load("data/RSTR_RNAseq_dat.voom.RData")

#### Data types ####
# SIMPLE
class(snp)
dim(snp)

#Extract a column
#Each column is itself a vector
snp$snpID
class(snp$snpID)
class(snp$CHR)

# COMPLEX
# S3 (accessed with $) or S4 (accessed with @)
class(dat.norm.voom)
class(dat.norm.voom$genes)
class(dat.norm.voom$E)

#Got further inside an S3 object with more $
dat.norm.voom$genes$symbol
class(dat.norm.voom$genes$symbol)

#### Basic operations ####
# R is a calculator
4+4

# On vectors of dimension 1
var(dat.norm.voom$targets$lib.size)
# also mean, sd, median, mode

# Greater than
dat.norm.voom$targets$lib.size > 1E9

#List unique values
unique(snp$CHR)

# Be careful of types
var(snp$snpID)

# On data frames /matrices /tables: dimension 2
# noted as [rows , columns]
snp$snpID[5]

snp[5,1]
snp[5, "snpID"]

#R is indexed starting at 1

#### Subsetting data ####
logical.vector <- dat.norm.voom$targets$lib.size > 10E6
logical.vector

dat.norm.voom$targets[logical.vector, "sampID"]
#Equal to
dat.norm.voom$targets[dat.norm.voom$targets$lib.size > 10E6,
                      "sampID"]
#Hard to read but we'll see a better way in the tidyverse!

#### tidyverse ####
#Syntax
# verb(data.frame, parameters)

#### Data 2.0 ####
snp <- read_csv(file = "data/Hawn_RSTR_SNPlist.PRKAG2.csv")
# Correctly names columns starting the numbers (unlike read.csv)

# Codes numbers as "doubles" which are more efficient and exact
str(snp)

#Rdata the same as before
load("data/RSTR_RNAseq_dat.voom.RData")

#### Separate out tables #####
#Sample metadata
meta <- dat.norm.voom$targets

#SNP metadata
#options
meta.snp <- snp[ , 1:9]

## TIDYVERSE SELECT COLUMNS
meta.snp <- select(snp, 
                   snpID, rsID, gene_id, symbol, type,
                   CHR, POS, allele.0, allele.1)
meta.snp <- select(snp,
                   snpID:allele.1)
meta.snp <- select(snp,
                   snpID, rsID, gene_id:allele.1)

#Genotypes
geno <- select(snp, 
               snpID, `84165-1-06`:`94295-1-02`)

#### Wide to long format data ####
## TIDYVERSE PIVOT
geno <- pivot_longer(geno,
                     -snpID, names_to = "FULLIDNO",
                     values_to = "genotype")

#### Mutate 0/0 to 0 ####
# Factors
# We use factors as they are more efficient. If you have data like
# cow   cow   llama, this takes 11 characters of space. But if you 
# re-code to a factor, R saves it as 1  1   2 where is knows that 
# 1=cow and 2=llama. Thus, this takes only 3 characters of space

# Moreover, changing 1=cow to 1=horse is faster because it only has to
# change it in the key for factor levels and the actual data of 
# 1   1   2 stays the same

#create a new variable
#Showing the steps

# Reformat to factor
# geno <- mutate(geno, geno.num = factor(genotype))

# Re-code the levels of the factor
# geno <- mutate(geno, 
#         geno.num2 = recode_factor(geno.num,
#                                   "0/0"=0,
#                                   "0/1"=1,
#                                   "1/1"=2))

# Convert to numeric, going through character first to keep the actual values 0,1,2 instead of the factor levels 1,2,3
# geno <- mutate(geno,
#         geno.num3 = as.numeric(as.character(geno.num2)))

#But actually
geno <- mutate(geno, 
  geno.num = as.numeric(as.character(recode_factor(factor(genotype),
                                 "0/0"=0,
                                 "0/1"=1,
                                 "1/1"=2))))
#And even though this is the tidyverse, this is difficult to read
# later on we'll improve this

#### long to wide format ####
#Select only the data you want in the table
geno <- select(geno, -genotype)
#move the snpID to column names and the numeric genotypes to the values in those columns
geno <- pivot_wider(geno, 
                    names_from = "snpID", values_from = "geno.num")

#### tranpose ####
# If you just needed to transpose a table
# But t() works on matrix not data frame so it does not play well with the tidyverse
#remake orig geno
geno <- select(snp, 
               snpID, `84165-1-06`:`94295-1-02`)
geno.transpose <- t(geno)
