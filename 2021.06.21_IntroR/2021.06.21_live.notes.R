#2021.06.21
#Intro R live notes

#### Install packages ####
#You only need to run this once on your computer
#Install package from CRAN
install.packages("tidyverse", Ncpus=2)

#Install package from Bioconductor
install.packages("BiocManager")
BiocManager::install("limma")

#### Ask R for help ####
#Ask for help
#Using the exact function name
?install
?select
#If you're unsure of the function name and want to search all the text documentation
??select
??pubmed

#### Load packages and set a seed for reproducibility ####
#You need to run this EVERY TIME you open R
library(tidyverse)
library(limma)

#Set seed
set.seed(5421)

#### Load data ####
#Load data from a table file
#Reads into the console
##No formatting
read.table(file="data/RSTR_meta_subset.csv")
##Use first row as column names
read.table(file="data/RSTR_meta_subset.csv", header=TRUE)
##Use first row as column names and split data by the comma separator
read.table(file="data/RSTR_meta_subset.csv", header=TRUE, sep=",")

#Saves to your R environment under the name 'meta'
meta <- read.table(file="data/RSTR_meta_subset.csv", header=TRUE, sep=",")
#Shorter option that assumes header=TRUE and sep="," by default
meta <- read.csv(file="data/RSTR_meta_subset.csv")

#Load data from RData
load("data/RSTR_data_clean_subset.RData")

#### Data types ####
#Simple -> meta
##A data frame has 2 dimensions
class(meta)
dim(meta)

##Get a single column (variable)
meta$FULLIDNO

#As what type of data are in a column
##character
class(meta$FULLIDNO)
class(meta$libID)
##numeric
class(meta$lib.size)

#not in our data
#logical = TRUE or FALSE

#Complex -> dat
#A list of 3 data frames
#It has 3 dimensions: length of list plus rows and columns in data frames
# in the list
class(dat)
#Get each data frame in the list
class(dat$targets)
class(dat$E)
class(dat$genes)
#Get a variable (column) within the data frame within the list
class(dat$targets$lib.size)

#### Working with data vectors ####
#Run simple statistics (functions) on data
mean(dat$targets$lib.size)
var(dat$targets$lib.size)
sd(dat$targets$lib.size)

#ask a question
# is lib.size greater than 10 million
dat$targets$lib.size > 10E6
# is lib.size greater than or equal to 10 million
dat$targets$lib.size >= 10E6

#see all possible value
unique(dat$targets$condition)

#using the correct type of data (class)
#Using the incorrect type of character in a function that expects numeric
#It runs but you get a warning
mean(dat$targets$libID)

#error vs warning
#The above warning still runs the mean function 
#wheras the error below does not run mean because there is no object
# in your environment named 'data'
mean(data$targets$lib.size)

#### Subsetting data ####
#data frame
# [rows, columns]
meta[5,]
meta[,5]
meta[5]
meta[5,5]

#vector (variable column)
meta$libID[5]
meta$libID[5,]

#Subset to rows where condition is "MEDIA
meta$condition == "MEDIA"
logical.vector <- meta$condition == "MEDIA"

meta[logical.vector, ]

#Shorthand where you don't have to save the logical.vector
meta[meta$condition == "MEDIA", ]

#Addtl examples
## rows where library size is greater than or equal to 10 million
meta[meta$lib.size >=10E6, ]
## Two ways to get rows with one of two RSID values
meta[meta$RSID %in% c("RS102306","RS102244"), ]
meta[meta$RSID == "RS102306" | meta$RSID == "RS102244", ]
## Using & means BOTH statement need to be TRUE.
## Which is 0 rows in this case
meta[meta$RSID == "RS102306" & meta$RSID == "RS102244", ]

#### Practive exercises ####
# Using help to identify the necessary arguments for the log function, compute the natural logarithm of 4, base 2 logarithm of 4, and base 4 logarithm of 4.
?log
log(4)
log(4, base=2)
log(4, 2) # note that if you push things in the EXACT order as the defaults
          # you don't need to specify the parameter name like base=
log(4,4)

# Using the meta data frame:
# Using an R function, determine what data type the norm.factors variable is
class(meta$norm.factors)

# Using indexing and the square bracket operator []:
# determine what RSID value occurs in the 20th row
meta[20,"RSID"]
#OR
meta$RSID[20]

# return the cell where lib.size equals 14509963
meta[meta$lib.size == 14509963, ]

meta[meta$lib.size == meta$lib.size[1], ]
meta[meta$lib.size == 4575023, ]

# Subset the data to observations where RSID equals “RS102521” or “RS102484”
meta[meta$RSID %in% c("RS102521","RS102484"), ]
#OR
meta[meta$RSID == "RS102521" | meta$RSID == "RS102484", ]
