library(tidyverse)

dat <- read_csv("data/counts_sub.csv")

dat.mat <-dat %>% 
  column_to_rownames("geneName") %>% 
  as.matrix()

#### Apply ####
# syntax
## apply(matrix, MARGIN, FUN)
## MARGIN=1 (rows), MARGIN=2 (cols)

apply(dat.mat, MARGIN = 1, FUN = mean)
apply(dat.mat, MARGIN = 2, FUN = mean)
apply(dat.mat, MARGIN = c(1,2), FUN = mean)

apply(dat.mat, MARGIN = 2, FUN = function(x) sum(x)/length(x))

#### Compare summarize ####
dat %>% 
  summarise_at(vars(-geneName), mean)

dat %>% 
  summarise_if(is.numeric, mean)

dat %>% 
  pivot_longer(-geneName, names_to = "sampID", values_to = "count") %>% 
  group_by(sampID) %>% 
  summarise(mean = mean(count))

dat %>% 
  pivot_longer(-geneName, names_to = "sampID", values_to = "count") %>% 
  separate(sampID, into=c("ptID", "cell", "TB"), sep="_") %>% 
  group_by(TB) %>% 
  summarise(mean = mean(count))

#### Compare loop ####

#ncol(dat.mat)

result <- data.frame()
for(i in 1:4){
  print(i)
  
  dat.sub <- dat.mat[,i]
  #print(dat.sub)
  
  result.temp <- mean(dat.mat[,i])
  print(result.temp)
  
  result <- rbind(result, result.temp)
}#

result
#Add names
result$sampID <- colnames(dat.mat)
result

#####
result <- data.frame()
for(sample in colnames(dat.mat)){
  print(sample)
  
  dat.sub <- dat.mat[,sample]
  #print(dat.sub)
  
  result.temp <- mean(dat.sub)
  #print(result.temp)
  
  result.temp.df <- data.frame(mean = result.temp,
                               sampID = sample)
  #print(result.temp.df)
  
  result <- rbind(result, result.temp.df)
}

result <- data.frame()
for(sample in colnames(dat.mat)){
  print(sample)
  
  result.temp.df <- data.frame(mean = mean(dat.mat[,sample]),
                               sampID = sample)
  #print(result.temp.df)
  
  result <- rbind(result, result.temp.df)
}

#### Apply family functions ####
# apply() for mat and df
# lapply( ) for list
# mapply( ) for multivariate
#
for(gene in dat$geneName){
  dat.format <- dat %>% 
  #plot <- dat %>% 
    filter(geneName == gene) %>% 
    pivot_longer(-geneName) %>% 
    separate(name, into=c("ptID", "cell", "TB"), sep="_") %>% 
  #print(dat.format)
    
    ggplot(aes(x=TB, y=value)) +
    geom_boxplot() +
    ggtitle(gene)
  #print(plot)
  
  dir.create("figures/", showWarnings = FALSE)
  filename <- paste("figures/", gene, ".pdf", sep="")
  ggsave(filename, plot)
}

#### Turn into function ####
plot.fxn <- function(input){
  if(is.data.frame(input)){
    print("Date are data frame.")
  } else{
    stop("Date are not data frame.")
  }
  
  for(gene in input$geneName){
    #dat.format <- input %>% 
    plot <- dat %>% 
      filter(geneName == gene) %>% 
      pivot_longer(-geneName) %>% 
      separate(name, into=c("ptID", "cell", "TB"), sep="_") %>% 
      #print(dat.format)
      
      ggplot(aes(x=TB, y=value)) +
      geom_boxplot() +
      ggtitle(gene)
    #print(plot)
    
    dir.create("figures/", showWarnings = FALSE)
    filename <- paste("figures/", gene, ".pdf", sep="")
    ggsave(filename, plot)
  }
}

plot.fxn(input = dat)
plot.fxn(input = dat.mat)


#####
result <- data.frame()
for(gene in dat$geneName){
  require(broom)
  
  result <- dat %>% 
    filter(geneName == gene) %>% 
    pivot_longer(-geneName) %>% 
    separate(name, into=c("ptID", "cell", "TB"), sep="_") %>% 
    aov(value ~ TB, data=.) %>% 
    
    tidy() %>% 
    mutate(geneName = gene) %>% 
    bind_rows(., result)
}

#### Parallel ####
library(data.table)
library(doParallel)
library(foreach)
library(broom)
registerDoParallel(cores=2)

result <- data.frame()
result <- rbindlist(foreach(gene = dat$geneName) %dopar% {
  result <- dat %>% 
    filter(geneName == gene) %>% 
    pivot_longer(-geneName) %>% 
    separate(name, into=c("ptID", "cell", "TB"), sep="_") %>% 
    aov(value ~ TB, data=.) %>% 
    
    tidy() %>% 
    mutate(geneName = gene) %>% 
    bind_rows(., result)
})
