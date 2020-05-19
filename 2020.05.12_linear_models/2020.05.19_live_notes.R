# Original workshop materials
# https://github.com/EDUCE-UBC/workshops_R/blob/master/stats_models/notes/stats_models_notes.pdf

#### BEFORE WE BEGIN PLEASE #### 

# 1. Create (or open last week's) R project for this workshop and create a new working script

# 2. Install packages, if needed
install.packages(c("tidyverse", "broom", "plyr", 
                   "lme4", "lmerTest", "car", "lsmeans", "MASS",
                   "gapminder", "segmented")) # <--- New this week

#### Setup #### 
# Always...
library(tidyverse)
# The data
library(gapminder)

#### Using R package data
# Some packages (like this one) are hidden and need to be assigned a variable name to see it in the environment
gapminder <- gapminder

#### Explore the data #####
# How does life expectancy change over time?
gapminder %>% 
  ggplot(aes(x = year, y = lifeExp)) +
  geom_point()
# How does GDP relate to life expectancy?
gapminder %>% 
  ggplot(aes(x = gdpPercap, y = lifeExp)) +
  geom_point()

#### Simple linear model ####
## Formula syntax y ~ x
lm.simple <- lm(lifeExp ~ year, data = gapminder)

#look at results
summary(lm.simple)
#A cleaner way
library(broom)
tidy(lm.simple)

#Look at fit
gapminder %>% 
  ggplot(aes(x = year, y = lifeExp)) +
  geom_point() +
  geom_smooth(method = "lm")

par(mfrow=c(2,2))
plot(lm.simple)
# Explanation of these plots https://data.library.virginia.edu/diagnostic-plots/

## Other variable of interest
lm.simple2 <- lm(lifeExp ~ gdpPercap, data = gapminder)
tidy(lm.simple2)
#Assess
par(mfrow=c(2,2))
plot(lm.simple2)
# --> Not a good linear fit as we saw in the original plot

##### Multiple linear regression #####
# e.g. Add more variables!
lm.multi <- lm(lifeExp ~ year + gdpPercap, data = gapminder)
#Results
tidy(lm.multi)
par(mfrow=c(2,2))
plot(lm.multi)


# With log transformed predictor
# sometimes called a "log model" but NOT the same as a logistic model
gapminder %>% 
  ggplot(aes(x=log(gdpPercap), y=lifeExp)) +
  geom_point()

lm.multi2 <- lm(lifeExp ~ year + log(gdpPercap), 
               data = gapminder)

#Same as
## can be used with tidyverse to manipulate multiple variables or perform complex manipulations as needed

lm.multi2 <- gapminder %>% 
  mutate(gdpPercap.log = log(gdpPercap)) %>% 
  
  lm(lifeExp ~ year + gdpPercap.log, data=.)

#Results
tidy(lm.multi2)
par(mfrow=c(2,2))
plot(lm.multi2)

#### Other transformations #####
#polynomial
poly()
# Or code with I in formula to code a specific function

# EXAMPLE
lifeExp ~ year + poly(gdpPercap)
lifeExp ~ year + gdpPercap + I(gdpPercap^2)

#DOES NOT equal below. This just squares the GDP values
lifeExp ~ year + gdpPercap + gdpPercap^2

#### Analysis of covariance ANCOVA ####
## e.g. interaction terms
gapminder %>% 
  ggplot(aes(x=year, y=lifeExp, color = continent)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  geom_smooth(method="lm", se=FALSE, color="black")

lm.ancova <- lm(lifeExp ~ year*continent,
                data=gapminder)
tidy(lm.ancova)

#### Segmented models ####
# e.g. multiple linear lines with different slopes
# good ref https://rpubs.com/MarkusLoew/12164
library(segmented)

lm.segment <- segmented(lm.ancova, 
                        seg.Z = ~year,
                        psi=c(1990)) #Break point
#Alternatively, use npsi= to allow model to esimate a given number of breaks

#Results
tidy(lm.segment)

#List break point(s)
lm.segment$psi

#Plot segments
gapminder %>% 
  mutate(year.group = ifelse(year<=1990, "group1", "group2")) %>% 
  
  ggplot(aes(x=year, y=lifeExp, color = continent)) +
  geom_point() +
  geom_smooth(aes(group=paste(year.group, continent)), 
              method="lm", se=FALSE) +
  geom_smooth(method="lm", se=FALSE, color="black")

#### Linear mixed effects models ####
## see last week's
library(lme)
library(lmerTest)

lmer(y ~ x + (1 | patient))

#### General linear model ####
# logistic regression (family = binomial)
# when outcome is binary (yes/no)

ucb <- as.data.frame(UCBAdmissions) 
#If you wanted to change the reference group, use factor levels like
  # mutate(Gender2 = factor(Gender, levels=c("Female", "Male")))

#counts lines for each group and uses that as frequency
lm.binomial <- glm(Admit ~ Gender * Dept, data=ucb,
                   family=binomial)
#No significance because all Gender+Dept groups have the same number of lines in the data
tidy(lm.binomial)

#uses weights variable as frequency. Each group has 1 line in data
glm.model <- glm(Admit ~ Gender * Dept, data=ucb,
    family=binomial, weights=Freq)

tidy(glm.model)

#Summarize p-values for variables
library(car)
car::Anova(glm.model)

#Simplifying models 
#Is my model improved by dropping the Dept variable?
## You want low AIC
drop1(glm.model, test="Chisq", scope="Dept")

# AIC https://en.wikipedia.org/wiki/Akaike_information_criterion
# BIC https://en.wikipedia.org/wiki/Bayesian_information_criterion

# Getting probabilities for groups
library(lsmeans)

#What is the probability of admittance for M vs F in each dept
glm.summ <- lsmeans(glm.model, ~ Gender + Dept, 
                    type="response")
glm.summ

#Get pairwise comparisons of probabilities
contrast(glm.summ, "pairwise", by="Gender")

#### GLM count data (family=poisson) ####
# Count data
# Not used for sequencing data anymore
# See notes linked at top for an example


