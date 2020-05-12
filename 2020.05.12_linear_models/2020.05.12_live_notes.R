###### BEFORE WE BEGIN PLEASE ######

# 1. Create an R project for this workshop and place the data file in the project directory under data/

# 2. Install packages
#install.packages(c("tidyverse",
#                   "broom",
#                   "plyr",
#                   "lme4",
#                   "car",
#                   "lsmeans",
#                   "MASS"))

##### Setup #####
library(tidyverse)

dat <- read_csv("data/AM.MDM.data.cleaning.metrics.csv") %>% 
  separate(sampID, into=c("donorID","cell","treatment"),
           sep="_", remove = FALSE)
  
#Set global theme for global
theme_set(theme_classic())
  
#### ANOVA ####
# Assumptions
# - normal data
# - balanced data
# - independent x-variables

ggplot(dat, aes(x = cell, y = raw)) +
  geom_jitter(width = 0.2, height = 0) +
  stat_summary(fun = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..), color="black") +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.1)

# Equations
# y ~ x1 + x2...
# OR
# x = , y =

# Run 1-way ANOVA
cell_aov <- aov(raw ~ cell, data = dat)

# equivalent to
aov(dat$raw ~ dat$cell)

# Getting p-value results
## Gross way
cell_aov_summ <- summary(cell_aov)
cell_aov_summ[[1]][ ,'Pr(>F)']

## A better way
library(broom)

cell_aov_summ2 <- tidy(cell_aov)
cell_aov_summ2$p.value

# Run 2-way ANOVA
## Additive
add_aov <- aov(raw ~ cell + treatment, data = dat)
tidy(add_aov)

## Interation
interact_aov <- aov(raw ~ cell * treatment, data = dat)
#equivalent to
interact_aov <- aov(raw ~ cell + treatment + cell:treatment, 
                    data = dat)
tidy(interact_aov)

#Interaction
# does the impact of the second term differ between groups in the first term
# does the impact of treatment differ by cell type

ggplot(dat, aes(x = paste(cell, treatment, sep=":"), 
                y = raw)) +
  geom_jitter(width = 0.2, height = 0) +

  stat_summary(fun = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..), color="black") +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.1)

#Variable order is important
ggplot(dat, aes(x = paste(treatment, cell, sep=":"), 
                y = raw)) +
  geom_jitter(width = 0.2, height = 0) +
  
  stat_summary(fun = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..), color="black") +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.1)

##### Linear model #####
tidy(add_aov)

add_lm <- lm(raw ~ cell + treatment, data = dat)
tidy(add_lm)

# Use factor to force variable order
dat <- dat %>% 
  #A way to order ALL levels of a factor
  mutate(cell_fct = factor(cell, levels=c("MDM", "AM"))) %>% 
  #A way to reorder 1 level of a factor
  mutate(treatment_fct = fct_relevel(
                          factor(treatment),
                          "TB", after=1))

tidy(lm(raw ~ cell_fct + treatment, data=dat))

###### Repeated measure ######
# Repeated measures in donors
tidy(aov(raw ~ cell + Error(donorID), data=dat))

library(lme4)
cell_lme <- lmer(raw ~ cell + (1 | donorID), data=dat)
car::Anova(cell_lme)

library(lmerTest)
#TBD

#### No intercept #####
#Good for 3+ groups without reference group like media
tidy(aov(raw ~ cell, data=dat))
#Though this model is not the right one for these data
tidy(aov(raw ~ 0 + cell, data=dat))

###### Assess your model ######
par(mfrow=c(2,2))
plot(cell_aov)
dev.off()

seq_lm <- lm(align_filter.paired ~ raw, data=dat)
par(mfrow=c(2,2))
plot(seq_lm)
