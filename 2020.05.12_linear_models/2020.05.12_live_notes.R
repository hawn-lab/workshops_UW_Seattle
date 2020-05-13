# Original workshop materials
# https://github.com/EDUCE-UBC/workshops_R/blob/master/stats_models/notes/stats_models_notes.pdf

#### BEFORE WE BEGIN PLEASE #### 

# 1. Create an R project for this workshop and place the data file in the project directory under data/

# 2. Install packages
#install.packages(c("tidyverse", "broom", "plyr", "lme4", "lmerTest", "car", lsmeans", "MASS"))
## We will not use all today but this will prepare you for the next several workshops.

#### Setup #### 
library(tidyverse)

dat <- read_csv("data/AM.MDM.data.cleaning.metrics.csv") %>% 
  separate(sampID, into=c("donorID","cell","treatment"),
           sep="_", remove = FALSE)
  
#Set global theme for plots
theme_set(theme_classic())
  
#### ANOVA ####
# Assumptions
# - normal data
# - balanced data
# - independent x-variables

#visualize data with means and error bars
ggplot(dat, aes(x = cell, y = raw)) +
  geom_jitter(width = 0.2, height = 0) +
  stat_summary(fun = mean, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..), color="black") +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.1)

# Equations
## as formula (preferred)
## y ~ x1 + x2...

## As parameters in function (not allowed in some functions)
## x = , y =

#### 1-way ANOVA #### 
# Do total raw sequences differ by cell type (MDM vs. AM)?
cell_aov <- aov(raw ~ cell, data = dat)

# which is equivalent to
cell_aov <- aov(dat$raw ~ dat$cell)

# Getting p-value results
## Gross way
cell_aov_summ <- summary(cell_aov)
cell_aov_summ[[1]][ ,'Pr(>F)']

## A better way
library(broom)

cell_aov_summ2 <- tidy(cell_aov)
cell_aov_summ2$p.value

#### 2-way ANOVA #### 
## Additive
add_aov <- aov(raw ~ cell + treatment, data = dat)
tidy(add_aov)

## Interation
interact_aov <- aov(raw ~ cell * treatment, data = dat)
# equivalent to
interact_aov <- aov(raw ~ cell + treatment + cell:treatment, 
                    data = dat)
tidy(interact_aov)

# Interaction meaning
# does the impact of the second term differ between groups in the first term
# example here: does the impact of treatment differ by cell type?

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
#Compare to previous 2-way ANOVA
tidy(add_aov)

add_lm <- lm(raw ~ cell + treatment, data = dat)
tidy(add_lm)

## SAME p-values because ANOVA is just a special case of a linear model

# Use factors to force variable order. The first level is the reference in the model
dat <- dat %>% 
  #A way to order ALL levels of a factor. All possible levels must be listed
  mutate(cell_fct = factor(cell, levels=c("MDM", "AM"))) %>% 
  #A way to reorder 1 level of a factor
  # Use additional parameter "after" to move to a place other than first
  mutate(treatment_fct = fct_relevel(
                          factor(treatment),
                          "TB"))

tidy(lm(raw ~ cell_fct + treatment, data=dat))

###### Repeated measure linear model ######
# The underlying math is the same, but each function requires different notation for repeated measures (or any other random effects in the model)

# Repeated measures in donors
# In ANOVA
tidy(aov(raw ~ cell + Error(donorID), data=dat))

# In linear model
library(lme4)
library(car)
cell_lme <- lme4::lmer(raw ~ cell + (1 | donorID), data=dat)
car::Anova(cell_lme)

#Or with another package that automatically gives P-values
library(lmerTest)
cell_lme2 <- lmerTest::lmer(raw ~ cell + (1 | donorID), data=dat)
summary(cell_lme2)

#### Remove intercept in a model #####

#Good for 3+ groups without a reference group
# when you want to later do pairwise comparisons like A-B, B-C, A-C
# Though note that this model is not the right one for these data
tidy(aov(raw ~ 0 + cell, data=dat))

###### Assess your model ######
par(mfrow=c(2,2))
plot(cell_aov)
dev.off()

#example of a true linear model with 2 numeric variables
seq_lm <- lm(align_filter.paired ~ raw, data=dat)
par(mfrow=c(2,2))
plot(seq_lm)
