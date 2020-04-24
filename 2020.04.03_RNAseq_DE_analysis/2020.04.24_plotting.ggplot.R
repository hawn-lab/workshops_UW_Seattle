#Exploring more plotting in ggplot

##### Setup #####
library(tidyverse)
library(edgeR)
set.seed(678)

load("data/AM.clean.RData")

#Extract separate df of counts and meta
counts <- as.data.frame(dat.AM.pc.abund.norm.voom$E)
meta <- as.data.frame(dat.AM.pc.abund.norm.voom$targets)

##### Boxplot of gene expression #####

counts %>% 
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to = "sampID", 
               values_to = "voom.count") %>% 
  left_join(meta, by="sampID") %>% 
  #Example for TBCD
  filter(geneName == "ENSG00000141556") %>%  
  
  ggplot(aes(x=TB, y=voom.count)) +
  geom_boxplot(outlier.shape = NA, 
               outlier.size = 0) +
  geom_jitter(height=0, width=0.2, aes(color=TB)) +
  theme_classic()

##### Violin plot ####
counts %>% 
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to = "sampID", 
               values_to = "voom.count") %>% 
  left_join(meta, by="sampID") %>%  
  #Example for TBCD
  filter(geneName == "ENSG00000141556") %>%  
  
  ggplot(aes(x=TB, y=voom.count)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  geom_jitter(height=0, width=0.2) +
  theme_classic()

##### more on boxplots #####
## Add a slope line
counts %>% 
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to = "sampID", 
               values_to = "voom.count") %>% 
  left_join(meta, by="sampID") %>% 
  #Change levels of TB variable
  mutate(TB.num = as.numeric(as.character(recode_factor(TB, Media=0, TB=1)))) %>% 
  #Example for TBCD
  filter(geneName == "ENSG00000141556") %>%  
  
  ggplot(aes(x=TB.num, y=voom.count)) +
  geom_boxplot(aes(group=TB.num),
               outlier.shape = NA) +
  geom_jitter(height=0, width=0.2, aes(color=TB)) +
  geom_smooth(method="lm", se=FALSE, color="black") +
  theme_classic()

##### encircling dots on scatterplot ######
# https://luisdva.github.io/rstats/Grouping-points/
base_sp <- counts %>% 
  rownames_to_column("geneName") %>% 
  pivot_longer(-geneName, names_to = "sampID", 
               values_to = "voom.count") %>% 
  left_join(meta, by="sampID") %>% 
  #Example for TBCD
  filter(geneName %in% c("ENSG00000141556", 
                         "ENSG00000188976")) %>% 
  pivot_wider(names_from = geneName, 
              values_from = voom.count) %>% 
  
  ggplot(aes(x=ENSG00000141556, y=ENSG00000188976,
             color=TB)) +
  geom_point( )

#95% confidence interval
base_sp + stat_ellipse( ) 

#Smooth outline
library(ggalt)
base_sp + geom_encircle(expand=0)

#Jagged outline
library(ggforce)
base_sp + geom_mark_hull(expand=0, radius = 0)

#### references ####
#Shapes
# http://sape.inf.usi.ch/quick-reference/ggplot2/shape

#Pick colors
# https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3

# Use math symbols
# https://www.calvin.edu/~rpruim/courses/s341/S17/from-class/MathinRmd.html
# http://csrgxtu.github.io/2015/03/20/Writing-Mathematic-Fomulars-in-Markdown/

#ggplot cheatsheet
# https://rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf