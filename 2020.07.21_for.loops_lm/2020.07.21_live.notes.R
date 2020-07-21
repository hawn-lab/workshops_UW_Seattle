# Running many linear models or plots with a for loop

"
for(something in something.else){

    DO THINGS
}
"

###### Step 1 ######
#write the code for the things you want to do in the loop
#Make sure you can run for 1 iteration before moving to a loop
library(tidyverse)
library(broom)
meta.prkag2 <- read_csv("data/prkag2.clean.csv")

model <- lm(PRKAG2 ~ condition*JHU_7.151349226 + M0_KCVAGE,
              data=meta.prkag2)
results <- tidy(model)
results

##### Step 2 #####
#loop over a couple of the things you want to use
# this allows you to find any errors before running all you data, 
# which would presumable take longer

#For example, i is often used to denote an integer
#Here were loop through integers 1 through 3 and 
# print the value and the value +100 to the console
for(i in c(1:3)){
  print(i)
  
  x = i+100
  print(x)
}

#Can name the temporary variable (i in previous example) anything
#Just be sure to use that name within the loop
for(magic in c(1:3)){
  print(magic)
  
  x = magic+100
  print(x)
}

#Our data
#List out all the SNPS you want to analyze
#Knowing these data, I know the SNPs are in columns 7 through 276
colnames(meta.prkag2)
colnames(meta.prkag2)[7:276]

# Loop through the first two SNPs
## But won't work because colnames() is character not a variable!
for(snpID in colnames(meta.prkag2)[7:8]){
  print(snpID)
  model <- lm(PRKAG2 ~ condition*snpID + M0_KCVAGE,
              data=meta.prkag2)
  results <- tidy(model)
  print(results)
}

#Instead use get() to tell R to "get the snpID variable" not just paste in the word
# A working loop that prints results to the console
for(snpID in colnames(meta.prkag2)[7:8]){
  print(snpID)
  model <- lm(PRKAG2 ~ condition*get(snpID) + M0_KCVAGE,
              data=meta.prkag2)
  results <- tidy(model)
  print(results)
}

#Or you can use integers and define snpID within the loop
#This is helpful when you have multiple things to input for each interation but only want 
# to have one loop

for(i in c(7:8)){
  #Get SNP of interest
  snpID <- colnames(meta.prkag2)[i]
  print(snpID)
  
  #Fit model
  model <- lm(PRKAG2 ~ condition*get(snpID) + M0_KCVAGE,
              data=meta.prkag2)
  #Estimate p-values
  results <- tidy(model)
  print(results)
}

##### step 3 ##### 
# save results and make self containing
# it is best practices to always have the packages you need within a loop with require()
# This facilitates turning the loop into a stand alone function later on

#Make blank df to hold results
results.all <- data.frame()

for(snp in colnames(meta.prkag2)[7:8]){
  #packages
  require(tidyverse)
  require(broom)
  
  #Print snpID so can see progress
  print(snp)
  
  #Fit model
  model <- lm(PRKAG2 ~ condition*get(snp) + M0_KCVAGE,
              data=meta.prkag2)
  #Estimate p-values
  results <- tidy(model)

  #Save results
  results.all <- results %>% 
    #Make new column with snpID for each iteration
    mutate(snpID = snp) %>% 
    #Merge with other iteration results
    bind_rows(results.all)
}

##### step 4 #####
# run on all data
# Once you know it works for a couple and you get your desired output, run on everything!
# Note that below, by using ncol(meta.prkag2) instead of a hard coded 276 for the last
# snp column, we could add more data later on and the exact same code would still work

#You do still have to be careful of the hard coded 7 below. For your date, check
# which column contains the first snp

#Make blank df to hold results
results.all <- data.frame()

for(snp in colnames(meta.prkag2)[7:ncol(meta.prkag2)]){
  #packages
  require(tidyverse)
  require(broom)
  
  #Print snpID so can see progress
  print(snp)
  
  #Fit model
  model <- lm(PRKAG2 ~ condition*get(snp) + M0_KCVAGE,
              data=meta.prkag2)
  #Estimate p-values
  results <- tidy(model)
  
  #Save results
  results.all <- results %>% 
    #Make new column with snpID for each iteration
    mutate(snpID = snp) %>% 
    #Merge with other iteration results
    bind_rows(results.all)
}

#It is best practices to save ALL model results and then apply multiple comparison 
# correction (FDR) and filter to results of interest after. This way you can always
# go back to the full models

#For example, calculating FDR and getting SNPs significant for the interaction of 
# TB and genotype (i.e. snp)

results.fdr <- results.all %>% 
  #Calculate FDR
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
  #Get significant snps
  filter(term == "conditionTB:get(snp)" &
         fdr <= 0.05)

##### Example 2: plots! #####
library(tidyverse)
meta.prkag2 <- read_csv("data/prkag2.clean.csv")

# to plot or not to plot: is interaction term significant?
# Using the results.fdr from above...

for(signif.snp in results.fdr$snpID){
  #Packages
  reguire(ggplot2)
  #Create plot
  plot <- ggplot(meta.prkag2, aes(x = get(signif.snp), y = PRKAG2)) +
    #boxplot
    geom_boxplot(aes(group = get(signif.snp)), outlier.color="white") +
    #add individual pts for each sample, colored by media/tb
    geom_jitter(aes(color=condition), width=0.2, height=0) +
    #add fit lines for media vs tb
    geom_smooth(aes(color=condition), method="lm", se=FALSE) +
    #Beautify
    labs(y="PRKAG2 normalized log2 expression",
         x=signif.snp) +
    scale_x_continuous(breaks=c(0,1,2)) +
    theme_classic()
  
  #Save each plot to PDF in figures directory
  dir.create("figures/", showWarnings = FALSE)
  plot.file <- paste("figures/", signif.snp, ".pdf", sep="")
  ggsave(filename = plot.file, plot, height=5, width=5)
}


##### TB-media plots #####
#Calculate TB-media
delta <- meta.prkag2 %>% 
  #Wide format with Media and TB values in own columns
  select(-name) %>% 
  pivot_wider(names_from = condition, values_from = PRKAG2) %>% 
  #Calculate TB-media for each sample
  group_by(across(c(sampID, FULLIDNO, JHU_7.151253265:rs116605521))) %>% 
  summarise(TB_media = TB-MEDIA) %>% 
  #Reorder columns
  select(sampID, FULLIDNO, TB_media, everything())

#Make a plot
ggplot(delta, aes(x = JHU_7.151253265, y = TB_media)) +
  #boxplot
  geom_boxplot(aes(group = JHU_7.151253265), outlier.color="white") +
  #add individual pts for each sample, colored by media/tb
  geom_jitter(width=0.2, height=0) +
  #add fit lines for media vs tb
  geom_smooth(color="darkgrey", method="lm", se=FALSE) +
  #Beautify
  labs(y="PRKAG2 TB - MEDIA log2 expression") +
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic()
