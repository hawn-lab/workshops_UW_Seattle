###### Setup #####
library(tidyverse)

#Load data
dat <- read_csv('data/counts_sub.csv')
#Create matrix format
dat.mat <- dat %>% 
  column_to_rownames("geneName") %>% 
  as.matrix()

#### Apply family functions ####
# apply() for mat and df
# lapply( ) for list
# mapply( ) for multivariate

##### Apply #####
# Some resources
# https://www.guru99.com/r-apply-sapply-tapply.html
# https://nicercode.github.io/guides/repeating-things/

# Syntax
# apply(data, MARGIN, FUNCTION)
# MARGIN=1 (rows), MARGIN=2 (columns), MARGIN=c(1,2) (both)
### How to remember MARGIN. In R, matrices are mat[rows, columns] with rows listed first (1) and columns second (2)

#find the mean expression of genes (rows) in each sample (columns)
apply(dat.mat, MARGIN = 2, FUN = mean)

#Save the result to a data frame
result <- as.data.frame(
  apply(dat.mat, MARGIN = 2, FUN = mean))

#find the mean expression of samples (columns) for each gene (rows) 
apply(dat.mat, MARGIN = 1, FUN = mean)

#apply to all unique combinations of row and column
#Mean does nothing since the mean of 1 number = that number
apply(dat.mat, MARGIN = c(1,2), FUN = mean)
#a better example, add 2 to every value
apply(dat.mat, MARGIN = c(1,2), FUN = function(x) x+2)

##### Summarize #####
# great when you already need to manipulate your data with the tidyverse
# good for data that is not already formatted to work with apply across rows or columns

#find the mean expression of genes (rows) in each sample (columns)
dat %>% 
  #Convert to long format
  pivot_longer(-geneName, names_to = "sampID", values_to = "count") %>% 
  #calculate mean for each sampID
  group_by(sampID) %>% 
  summarize(mean.count = mean(count))

#find the mean expression of genes (rows) in Media vs TB
dat %>% 
  #Convert to long format
  pivot_longer(-geneName, names_to = "sampID", values_to = "count") %>% 
  #separate sampID names into metadata
  separate(sampID, into=c("ptID","cell","TB"), sep="_") %>% 
  #calculate mean for each TB group (e.g. media vs. TB)
  group_by(TB) %>% 
  summarize(mean.count = mean(count))

##### Other summarize functions #####
#find the mean expression of genes (rows) in each sample (columns)
#Summarize only at specific variables (columns)
dat %>% 
  summarize_at(vars(AM10_AM_Media, AM10_AM_TB), mean)

#summarize all columns except geneName
dat %>% 
  summarize_at(vars(-geneName), mean)

#summarize all numeric columns
dat %>% 
  summarize_if(is.numeric, mean)

##### For loops #####
#Syntax

# for(something in something_else){
# DO STUFF to something
# }

## Simplest loop
#Make blank df to hold loop results
result <- data.frame()

#Loop through each column and calculate mean
for(i in 1:4){
  print(i)
  
  dat.sub <- dat.mat[ , i]
  #print(dat.sub)
  
  result.temp <- mean(dat.sub)
  #print(result.temp)
  
  result <- rbind(result, result.temp)
}
#Add names to result
result$sampID <- colnames(dat.mat)
# DANGER! Can cause mislabels!! You have to be 100% sure the order of the result data frame and whatever you're using as names is the same
# See below for a better way!

## Better naming
result <- data.frame()
for(sample in colnames(dat.mat)){
  print(sample)
  
  dat.sub <- dat.mat[ , sample]
  #print(dat.sub)
  #print(colnames(dat.sub))
  
  result.temp <- mean(dat.sub)
  #print(result.temp)
  
  # create a df with the sampID already in there
  result.temp.df <- data.frame(mean = result.temp,
                              sampID = sample)
  #print(result.temp.df)
  
  result <- rbind(result, result.temp.df)
}

## For real (shorter)
result <- data.frame()
for(sample in colnames(dat.mat)){
  result.temp.df <- data.frame(mean = mean(dat.mat[ , sample]),
                               sampID = sample)

  result <- rbind(result, result.temp.df)
}

##### plotting plots #####
#plot a boxplot of each gene media vs TB

for(gene in dat$geneName){
  plot <- dat %>% 
    #Filter to gene of interest
    filter(geneName == gene) %>% 
    #Convert to long format
    pivot_longer(-geneName) %>% 
    #Separate sample names into metadata
    separate(name, into=c("ptID","cell","TB"), sep="_") %>% 
    
    #plot
    ggplot(aes(x = TB, y = value)) +
    geom_boxplot() +
    #add title with gene name
    ggtitle(gene)
  
  #Save plot
  ## using showWarnings=FALSE means this won't fail if the folder already exists
  dir.create("figures/", showWarnings = FALSE)
  ## Create a dynamic file name based on the gene name
  filename <- paste("figures/", gene, ".pdf", sep="")
  ## save
  ggsave(filename, plot)
}

#### Create custom function ####
# Syntax

# function.name <- function(parameter1=..., parameter2=...){
#   STUFF
#  }

#Use
# source("Rscript.with.function.R")
# function.name(parameter1 = x, parameter2 = y)

#Can even source directly from GitHub. Make sure to use the RAW url
# source("https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/SNP.model.fxn.R")

plot.fxn <- function(input){
  #Make sure the data are formatted correctly
  if(is.data.frame(input)){
    print("Data are data frame. Continuing...")
  } else if(is.matrix(input)){
    print("Data are matrix. Reformatting.")
    
    input <- as.data.frame(input) %>% 
      rownames_to_column("geneName")
  } else {
    stop("Data are not data frame or matrix. Please correct.")
  }
  
  for(gene in input$geneName){
    plot <- dat %>% 
      #Filter to gene of interest
      filter(geneName == gene) %>% 
      #Convert to long format
      pivot_longer(-geneName) %>% 
      #Separate sample names into metadata
      separate(name, into=c("ptID","cell","TB"), sep="_") %>% 
      
      #plot
      ggplot(aes(x = TB, y = value)) +
      geom_boxplot() +
      #add title with gene name
      ggtitle(gene)
    
    #Save plot
    ## using showWarnings=FALSE means this won't fail if the folder already exists
    dir.create("figures/", showWarnings = FALSE)
    ## Create a dynamic file name based on the gene name
    filename <- paste("figures/", gene, ".pdf", sep="")
    ## save
    ggsave(filename, plot)
  
}}

#Now it works on both!
plot.fxn(dat)
plot.fxn(dat.mat)

#### Parallel processing loops ####
#for combining parallel results
library(data.table)
#For parallel computing
library(doParallel)
library(foreach)
registerDoParallel(cores=2)
#For getting model results
library(broom)

#Syntax
# foreach(something = something_else) %dopar% {
# STUFF done to something
#  }

#results as list of data frames

#Create blank list to hold results
result <- list()

result <- foreach(gene = dat$geneName) %dopar% {
  result <- dat %>% 
    #format data
    filter(geneName == gene) %>% 
    pivot_longer(-geneName) %>% 
    separate(name, into = c("ptID","cell","TB"), sep="_") %>% 
    #Run ANOVA
    aov(value ~ TB, data=.) %>% 
    #Get results table
    broom::tidy() %>% 
    #Add column with gene name
    mutate(geneName = gene)
    #Note that I forgot to remove this in the workshop but you don't need it once the foreach loop is directed to "result"
    # bind_rows(., result)
}

#Combine list of dfs into 1 df result
result <- data.frame()

result <- data.table::rbindlist(
  foreach(gene = dat$geneName) %dopar% {
    result <- dat %>% 
      #format data
      filter(geneName == gene) %>% 
      pivot_longer(-geneName) %>% 
      separate(name, into = c("ptID","cell","TB"), sep="_") %>% 
      #Run ANOVA
      aov(value ~ TB, data=.) %>% 
      #Get results table
      broom::tidy() %>% 
      #Add column with gene name
      mutate(geneName = gene)
})

