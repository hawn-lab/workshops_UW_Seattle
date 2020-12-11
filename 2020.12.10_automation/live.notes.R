### Key shortcuts ###
# Cmd/ctrl+Enter to run code
# Double click to select 1 word
# Triple click to select 1 line
# Under menu Edit>Folding find the shortcut for collapses 

#### Setup ####
# Load packages
library(tidyverse)
library(palmerpenguins)
library(doParallel)
library(foreach)

# Load data
dat <- penguins

#### Base R ####
# Data frames accessed as dat[rows, columns]
# 1. Get 1 species
dat.filter <- dat[dat$species == "Adelie", ]

# 2. Convert mm to cm
bill_length_cm <- dat.filter$bill_length_mm/10
bill_depth_cm <- dat.filter$bill_depth_mm/10

# 3. Plots
plot(x = dat.filter$sex, y = bill_length_cm,
     ylab = "bill length (cm)", xlab = "", 
     main = "Adelie")
plot(x = dat.filter$sex, y = bill_depth_cm,
     ylab = "bill depth (cm)", xlab = "", 
     main = "Adelie")

# Copy and replace for other 2 species. So much repeated code!

#### for loops ####
# for(something in something.else){
#   Do things
# }

#Print each penguin species
for(i in c("Adelie", "Chinstrap", "Gentoo")){
  print(i)
}

for(magic in c("Adelie", "Chinstrap", "Gentoo")){
  print(magic)
}

# Make R figure out the species for you
for(i in unique(dat$species)){
  print(i)
}
#Or if it's a factor variable
for(i in levels(dat$species)){
  print(i)
}

#Set up plot panel in 3 by 2 to fit all 6 plots
#par(mfrow=c(3,2))

# Loop to make plots
## Copy-paste all base R code and replace "Adelie" with i
for(i in unique(dat$species)){
  # 1. Get 1 species
  dat.filter <- dat[dat$species == i, ]
  
  # 2. Convert mm to cm
  bill_length_cm <- dat.filter$bill_length_mm/10
  bill_depth_cm <- dat.filter$bill_depth_mm/10
  
  # 3. Plots
  plot(x = dat.filter$sex, y = bill_length_cm,
       ylab = "bill length (cm)", xlab = "", 
       main = i)
  plot(x = dat.filter$sex, y = bill_depth_cm,
       ylab = "bill depth (cm)", xlab = "", 
       main = i)
}

#### apply functions ####
# apply() something across all rows (MARGIN=1)
# or across all columns (MARGIN=2)
# Remember which is which like dat[1 rows, 2 columns]

# Error b/c works across ALL columns and some of our data are not numeric
dat.filter.cm <- apply(dat.filter, MARGIN = 2,
                       function(x) x/10)

#Use indexes to pick specific columns
dat.filter.cm <- apply(dat.filter[c(3,4)], 
                       MARGIN = 2,
                       function(x) x/10)
#Or columns names
dat.filter.cm <- apply(X=dat.filter[c("bill_length_mm",
                                    "bill_depth_mm")], 
                       MARGIN = 2,
                       FUN=function(x) {x/10}) 
#Can include {} or not for function inside an apply function

#If you have a long function, it's easier to read if you separate it out like
my.function <- function(x){x/10}
#then use my.function in in apply
dat.filter.cm <- apply(X=dat.filter[c("bill_length_mm",
                                      "bill_depth_mm")], 
                       MARGIN = 2,
                       FUN=my.function) 

#### tidyverse mutate ####
#Change each variable individually
dat.mutate <- mutate(dat.filter, 
                     bill_length_cm = bill_length_mm/10,
                     bill_depth_cm = bill_depth_mm/10)

# across() works like apply but within mutate()
# dat.mutate <- mutate(dat.filter, across(vectors.of.variables, function))
dat.mutate <- mutate(dat.filter,
                     across(c(bill_length_mm, bill_depth_mm),
                            .fns = ~.x/10))
#Note does not change variable names even though the data changed
colnames(dat.mutate)
dat.filter$bill_length_mm[1:5]
dat.mutate$bill_length_mm[1:5]

# Automatic rename in mutate
# Can add suffix "{.col}SUFFIX"
# Or prefix "PREFIX{.col}"
dat.mutate <- mutate(dat.filter,
                     across(c(bill_length_mm, bill_depth_mm),
                            .fns = ~.x/10,
                            .names = "{.col}_cm"))
#see new columns
colnames(dat.mutate)

#Old way in tidyverse (before across was added in R4.0)
# mutate_all, mutate_at, mutate_if
dat.mutate <- mutate_at(dat.filter, 
                  vars(c(bill_length_mm, bill_depth_mm)),
                        ~.x/10)

#### More tidyverse ####
#Mutate can be piped ( %>% ) with other tidyverse functions
par(mfrow=c(1,2))
for(i in unique(dat$species)){
  dat.filter <- dat %>% 
    # Step 1: Filter data to the penguin species
    filter(species == i) %>% 
    # Step 2: Convert the mm measurements to cm
    mutate(across(c(bill_length_mm, bill_depth_mm), 
                  .fns = ~.x/10, .names = "{.col}_cm"))
  
  # Step 3: Plot measurements separated by sex
  plot(x = dat.filter$sex, y = dat.filter$bill_length_mm_cm, 
       ylab = "bill length (cm)", xlab = "", main = i)
  plot(x = dat.filter$sex, y = dat.filter$bill_depth_mm_cm, 
       ylab = "bill depth (cm)", xlab = "", main = i)
}

#### Nested loops ####
# Can be very slow but sometime necessary
# for(){
#   for(){
#     Do things
#   }
# }

for(penguin.species in unique(dat$species)){
  # 1. Get 1 species
  dat.filter <- dat[dat$species == penguin.species, ]
  
  for(variable in c("bill_length_mm", "bill_depth_mm")){
    # 2. Convert mm to cm
    size_cm <- dat.filter[[variable]]/10 #Note use of [[]] to get the bill variable instead of "variable", which is not a column name in our data
    
    # 3. Plots
    plot(x = dat.filter$sex, y = size_cm,
         ylab = variable, xlab = "",
         main = penguin.species)
  }
}

#### Base R vs tidyverse ####
# 1. Use a character string (like "name") as a variable 
## Tidy get(name)
## Base in a data frame dat[[name]] 
## Base as a vector !!name

# 2. Calling a function
## Tidy in mutate use .fns and placeholder .x or .
## Base apply use FUN and placeholder whatever you want (probably x)

# 3. variable names (not in loop)
## Tidy doesn't need surrounding " "

#### foreach ####
foreach(i = unique(dat$species)) %dopar% {
  
}
# is the same but on multiple processors as

for(i in unique(dat$species)) {
  
}

# set number of processors
registerDoParallel(cores=2)
#Lookup how many you have
parallel::detectCores()

foreach(i = unique(dat$species)) %dopar% {
  dat.filter <- dat[dat$species == i, ]
  
  bill_depth_cm <- dat.filter$bill_depth_mm/10
  bill_length_cm <- dat.filter$bill_length_mm/10
  
  filename <- paste(i, "pdf", sep=".")
  #print(filename)
  
  #foreach cannot print to the Plots pane in RStudio so you have to save the outputs instead
  pdf(file=filename)
  par(mfrow=c(1,2))
  plot(x = dat.filter$sex, y = bill_length_cm,
       ylab = "bill length (cm)", xlab = "", main = i)
  plot(x = dat.filter$sex, y = bill_depth_cm, 
       ylab = "bill depth (cm)", xlab = "", main = i)
  dev.off()
}


## Other resources
# https://www.r-bloggers.com/2012/12/using-apply-sapply-lapply-in-r/
# purrr() package for alternatives to apply functions
