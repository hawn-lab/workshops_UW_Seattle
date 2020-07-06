##### Packages #####
library(tidyverse)
library(limma)

##### Data #####
# Downloaded from Shiny app.
# SNP genoytpes
snp <- read_csv("data/Hawn_RSTR_SNPlist.PRKAG2.csv")
# RNA-seq data
load("data/RSTR_RNAseq_dat.voom.RData")

##### Pipes #####
# %>% is a pipe
snp <- select(snp, snpID)
snp <- pivot_longer(snp...)

# is equivalent to

snp <- snp %>% select(snpID) %>% pivot_longer(...)

# Mathematical notation
# f(x) %>% g(y) where y = f(x)
# equivalent to
# g(f(x), y)

##### Data cleaning #####
#Genotype data
geno <- read_csv("data/Hawn_RSTR_SNPlist.PRKAG2.csv") %>% 
  #Select genotype data
  select(snpID, `84165-1-06`:`94295-1-02`) %>% 
  #Convert to long format
  pivot_longer(-snpID, names_to="FULLIDNO", values_to="genotype") %>% 
  #Convert genotype 0/0 to 0 format
  mutate(geno.num = recode_factor(factor(genotype),
                                  "0/0" = 0,
                                  "0/1" = 1,
                                  "1/1" = 2)) %>% 
  mutate(geno.num = as.numeric(as.character(geno.num))) %>% 
  #Convert back to wide format
  select(snpID, FULLIDNO, geno.num) %>% 
  pivot_wider(names_from = "snpID", values_from = "geno.num")

#Metadata
meta <- dat.norm.voom$targets %>% 
  #Keep variables of interest
  select(libID, sampID, FULLIDNO, condition, M0_KCVAGE) %>% 
  #Keeps only with both metadata and geno data
  inner_join(geno, by="FULLIDNO")

#Other variable options
colnames(dat.norm.voom$targets)

#Counts (aka expression)
prkag2 <- as.data.frame(dat.norm.voom$E) %>% 
  #move rownames
  rownames_to_column("symbol") %>% 
  #Filter gene of interest
  filter(symbol == "PRKAG2")

#Merge counts and meta
meta.prkag2 <- prkag2 %>% 
  pivot_longer(-symbol) %>% 
  select(-symbol) %>% 
  inner_join(meta, by=c("name"="libID")) %>% 
  rename(PRKAG2 = value)

#Or

meta.prkag2 <- prkag2 %>% 
  pivot_longer(-symbol, names_to = "libID", values_to = "PRKAG2") %>% 
  select(-symbol) %>% 
  inner_join(meta)

##### Save tables #####
write_csv(meta.prkag2, path = "data/prkag2.clean.csv")

##### Linear model #####
# y ~ x1 + x2...
# gene expression ~ Media/TB * SNP.genotype + age

# Fit model
model <- lm(PRKAG2 ~ condition*JHU_7.151349226 + M0_KCVAGE,
            data=meta.prkag2)
#Estimate p-values
summary(model)

# install.packages("broom")
library(broom)
results <- tidy(model)

###### Assess model ######
par(mfrow=c(2,2))
plot(model)

###### Plotting in ggplot2 ######
#Boxplot
ggplot(meta.prkag2, aes(x=JHU_7.151349226, y=PRKAG2)) +
  geom_boxplot(aes(group=JHU_7.151349226))

#Boxplot with dots for samples
ggplot(meta.prkag2, aes(x=JHU_7.151349226, y=PRKAG2)) +
  geom_boxplot(aes(group=JHU_7.151349226),
               outlier.shape=NULL) +
  geom_jitter(width = 0.2, height = 0)

#Add linear fit
ggplot(meta.prkag2, aes(x=JHU_7.151349226, y=PRKAG2)) +
  geom_boxplot(aes(group=JHU_7.151349226),
               outlier.shape=NULL) +
  geom_jitter(width = 0.2, height = 0) +
  geom_smooth(method="lm")

ggplot(meta.prkag2, aes(x=JHU_7.151349226, y=PRKAG2,
                        color=condition)) +
  #Boxplots grouped by genotype 0,1,2
  geom_boxplot(aes(group=JHU_7.151349226),
               outlier.shape=NULL) +
  #Dots for samples
  geom_jitter(width = 0.2, height = 0) +
  #Linear fit
  geom_smooth(method="lm")


#More customization
ggplot(meta.prkag2,
       #Define x and y variables
       aes(x=JHU_7.151349226, y=PRKAG2)) +
  #Create boxplots
  ## Use group to force gglot to treat genotype like character 
  ## instead of numeric
  ## Set the outlier shape to NULL so that you don't get duplicate
  ## dots in the next layer.
  geom_boxplot(aes(group=JHU_7.151349226), outlier.shape=NULL) +
  #Add points for each sample. Color by MEDIA/TB and "jitter" left 
  # and right (width) to avoid overlap
  geom_jitter(aes(color=condition), width = 0.2, height=0) +
  # Add a linear fit for each condition
  geom_smooth(aes(color=condition), method="lm", se=FALSE) +
  #Add a linear fit for genotype overall
  geom_smooth(color="grey", method="lm", se=FALSE) +
  #Format axis labels
  labs(y="PRKAG2 normalized log2 expression") +
  scale_x_continuous(breaks=c(0,1,2)) +
  #Change theme to classic to remove grey background and other
  #default aspects of ggplot
  theme_classic()

##### Facets #####
meta.prkag2 %>% 
  select(name:rs5017429) %>% 
  pivot_longer(JHU_7.151253265:rs5017429,
               names_to = "snpID", 
               values_to = "genotype") %>% 
  
  ggplot(aes(x=genotype, y=PRKAG2)) +
  geom_boxplot(aes(group=genotype), outlier.shape = NULL) +
  geom_jitter(aes(color=condition), width=0.2, height=0) +
  geom_smooth(aes(color=condition), method="lm", se=FALSE) +
  geom_smooth(color="grey", method="lm", se=FALSE) +
  labs(y="PRKAG2 normalized log2 expression") +
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  facet_wrap(~snpID)


plot <- meta.prkag2 %>% 
  pivot_longer(-c(name:M0_KCVAGE),
               names_to = "snpID", 
               values_to = "genotype") %>% 
  
  ggplot(aes(x=genotype, y=PRKAG2)) +
  geom_boxplot(aes(group=genotype), outlier.shape = NULL) +
  geom_jitter(aes(color=condition), width=0.2, height=0) +
  geom_smooth(aes(color=condition), method="lm", se=FALSE) +
  geom_smooth(color="grey", method="lm", se=FALSE) +
  labs(y="PRKAG2 normalized log2 expression") +
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  facet_wrap(~snpID)

##### Save plots #####
ggsave(filename="PRKAG2.pdf", plot, height=22, width=22)

##### Filter for sample size #####
filter.key <- meta.prkag2 %>% 
  pivot_longer(-c(name:M0_KCVAGE),
               names_to = "snpID", 
               values_to = "genotype") %>% 
  #Calculate number of samples in each SNP genotype
  group_by(snpID, genotype) %>% 
  summarize(total = n()/2) %>% 
  ungroup() %>% 
  #Keep SNPs with at least 2 samples in geno 1 AND 2
  pivot_wider(names_from = genotype, values_from = total) %>% 
  filter(`1` > 3 & `2` > 3)

meta.prkag2 %>% 
  pivot_longer(-c(name:M0_KCVAGE),
               names_to = "snpID", 
               values_to = "genotype") %>% 
  filter(snpID %in% filter.key$snpID) %>% 
  
  ggplot(aes(x=genotype, y=PRKAG2)) +
  geom_boxplot(aes(group=genotype), outlier.shape = NULL) +
  geom_jitter(aes(color=condition), width=0.2, height=0) +
  geom_smooth(aes(color=condition), method="lm", se=FALSE) +
  geom_smooth(color="grey", method="lm", se=FALSE) +
  labs(y="PRKAG2 normalized log2 expression") +
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  facet_wrap(~snpID)


##### Filter for sample proportion #####
filter.key2 <- meta.prkag2 %>% 
  pivot_longer(-c(name:M0_KCVAGE),
               names_to = "snpID", 
               values_to = "genotype") %>% 
  #Calculate number of samples in each SNP genotype
  group_by(snpID, genotype) %>% 
  summarize(total = n()/2) %>% 
  ungroup() %>% 
  #Wide format
  pivot_wider(names_from = genotype, values_from = total) %>% 
  #Calculate total nonNA values and proportion
  mutate(total = `0`+`1`+`2`,
         prop.1 = `1`/total,
         prop.2 = `2`/total) %>% 
  #Keep SNPs with at least 2 samples in geno 1 AND 2
  filter(prop.1 > 0.2 & prop.2 > 0.2)

meta.prkag2 %>% 
  pivot_longer(-c(name:M0_KCVAGE),
               names_to = "snpID", 
               values_to = "genotype") %>% 
  filter(snpID %in% filter.key2$snpID) %>% 
  
  ggplot(aes(x=genotype, y=PRKAG2)) +
  geom_boxplot(aes(group=genotype), outlier.shape = NULL) +
  geom_jitter(aes(color=condition), width=0.2, height=0) +
  geom_smooth(aes(color=condition), method="lm", se=FALSE) +
  geom_smooth(color="grey", method="lm", se=FALSE) +
  labs(y="PRKAG2 normalized log2 expression") +
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_classic() +
  facet_wrap(~snpID)
