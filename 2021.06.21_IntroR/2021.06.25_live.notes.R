#### Load data and packages ####
#load packages
library(tidyverse)

#Load in RData object
load("data/RSTR_data_clean_subset.RData")

# Read in a csv
meta <- read_csv(file="data/RSTR_meta_subset.csv")
#Look at data types
class(meta$lib.size)
str(meta)
spec(meta)

#### Select ####
#Select one or more columns

#Base R
##meta[rows, columns]
meta[ , c("libID","FULLIDNO","condition")]

#Tidy way
## get all the columns you want by name
meta.sub <- select(meta, libID, FULLIDNO, condition)
##get all columns between 2 names with :
meta.sub <- select(meta, libID, FULLIDNO:condition)
##remove columns by name
meta.sub <- select(meta, -lib.size, -norm.factors)

##Look at the result
meta.sub

#if you have odd column names, use ` `
meta.sub <- select(meta, -`lib.size`, -`norm.factors`)

#### Filter ####
#Filter one or more rows were X is TRUE
#Base R
meta[meta$condition=="MEDIA",]

#tidy way
filter(meta, condition=="MEDIA")
#Another example
counts <- as.data.frame(dat$E)

#note that dat$E was already a data frame but in some RNAseq data sets,
# it is a matrix and the above as.data.frame() is necessary
class(dat$E)
class(counts)

#### Ronames ####
#move rownames into a data column
counts.rowname <- rownames_to_column(counts)

#Note you can do the opposite with 
column_to_rownames(counts.rowname)

#### Filter again ####
IFNG <- filter(counts.rowname, rowname=="IFNG")
IFNG
##No rows because the rownames are ENSEMBL ID (ENSG##), not symbols like IFNG

#One way around this is to find the correct ENSEMBL ID "by-hand"
#Find the ID for IFNG
filter(dat$genes, hgnc_symbol=="IFNG")
#Copy-paste into filter function
filter(counts.rowname, rowname=="ENSG00000111537")

#But there is a better (more automated) way!!

#### Joining data frame ####
## Basic syntax
## _join(data.frame1, data.frame2, by="column.name.to.match")

#keep rows present in both data frames
genes.counts <- inner_join(counts.rowname, dat$genes)
## error because there are no columns with the same name in BOTH data frames

#instead tell the join what to match, even if they have a different name
genes.counts <- inner_join(counts.rowname, dat$genes, by=c("rowname"="geneName"))

#Or we could have named to "rownames" column so it matched in the beginning
#Rename rowname column
counts.rowname <- rownames_to_column(counts, "geneName")
#join with default matching
genes.counts <- inner_join(counts.rowname, dat$genes)

#### Try filter again ####
#Now we have a column with the symbol name we'e expecting
IFNG <- filter(genes.counts, hgnc_symbol=="IFNG")

#### Pivot data frames ####
#Moving between long vs wide data
#Wide: 1 column per variable, 1 row per sample
#Long: 1 column for values of all variable, many rows for each sample

#These data are wide
head(IFNG)

#from wide to long
#default new columns are "name" and "value"
pivot_longer(IFNG, -c(geneName, hgnc_symbol:locus_group))

#or we can name them what we actually want them to be
IFNG.long <- pivot_longer(IFNG, -c(geneName, hgnc_symbol:locus_group),
             names_to = "libID", values_to = "IFNG")

#and we could undo this with a long to wide function
pivot_wider(IFNG.long, names_from = "libID", values_from = "IFNG")

#### More joining ####
#now we have a column "libID" to join by
IFNG.meta <- inner_join(IFNG.long, meta.sub)
IFNG.meta

#### Select again ####
#let's only keep the variables we really need
IFNG.meta.sub <- select(IFNG.meta, libID, IFNG, FULLIDNO, condition)
IFNG.meta.sub

#### Pipes ####
#All the above created a lot of trash in our R environment
# we could instead use pipes %>% to connect all our tidy functions!

#let's start over
load("data/RSTR_data_clean_subset.RData")
meta <- read_csv(file="data/RSTR_meta_subset.csv")

#CMD+SHIFT+M is the shortcut for a pipe
meta.sub <- meta %>% 
  select(libID, FULLIDNO, condition)
  
IFNG <- as.data.frame(dat$E) %>% 
  rownames_to_column("geneName") %>% 
  inner_join(dat$genes) %>% 
  filter(hgnc_symbol == "IFNG") %>% 
  pivot_longer(-c(geneName, hgnc_symbol:locus_group), 
               names_to = "libID", values_to = "IFNG") %>% 
  inner_join(meta.sub) %>% 
  select(libID, IFNG, FULLIDNO, condition)

#### plotting in ggplot ####
#boxplot with individual points for each sample
ggplot(data=IFNG, aes(x=condition, y=IFNG)) +
  geom_boxplot() +
  geom_point()

#jitter the points so you can see them all and beautify the plot a bit
ggplot(data=IFNG, aes(x=condition, y=IFNG, color=condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, height=0) +
  labs(y="IFNG normalized log2 expression", x="") +
  theme_classic()

#adding color based on a variable
ggplot(data=IFNG, aes(x=condition, y=IFNG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, height=0, aes(color=condition)) +
  labs(y="IFNG normalized log2 expression", x="") +
  theme_classic() +
  theme(legend.position = "none")

#adding a specific color not based on a variable
## since it's not a variable in IFNG, is should NOT be an aesthetic in aes()
ggplot(data=IFNG, aes(x=condition, y=IFNG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, height=0, aes(color="blue")) +
  labs(y="IFNG normalized log2 expression", x="") +
  theme_classic()
ggplot(data=IFNG, aes(x=condition, y=IFNG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, height=0, color="blue") +
  labs(y="IFNG normalized log2 expression", x="") +
  theme_classic()

#### Facets ####
#Make multiple plots all at once!!

as.data.frame(dat$E) %>% 
  #Move rownames to a column. Unlike before, let's name it 
  #to match what we'll be joining with next
  rownames_to_column("geneName") %>% 
  #Join with gene key to get HGNC symbol
  inner_join(dat$genes) %>% 
  #filter IFNG and TNF expression
  filter(hgnc_symbol %in% c("IFNG", "IFNA1")) %>% 
  #Pivot to long format so can combine with metadata
  pivot_longer(-c(geneName, hgnc_symbol:locus_group), 
               names_to="libID", values_to="expression") %>% 
  #join with library metadata
  inner_join(meta.sub) %>% 
  #select variables of interest. Note we need to keep hgnc_symbol to tell
  #data from the two genes apart
  select(libID, hgnc_symbol, expression, FULLIDNO, condition) %>% 
  
  ggplot(aes(x=condition, y=expression)) +
  #Create boxplots
  ## Set the outlier shape to NA so that you don't get duplicate
  ## dots in the next layer.
  geom_boxplot(outlier.shape=NA) +
  #Add points for each sample. Color by MEDIA/TB and "jitter" left 
  # and right (width) to avoid overlap
  geom_jitter(aes(color=hgnc_symbol), width = 0.2, height=0) +
  #Format axis labels
  labs(y="Normalized log2 expression", x="") + 
  #Change theme to classic to remove grey background and other
  #default aspects of ggplot
  theme_classic() +
  #Add facets
  # !!! THIS IS THE ALL POWERFUL SINGLE LINE !!!
  facet_wrap(~hgnc_symbol, scales = "free_y")

