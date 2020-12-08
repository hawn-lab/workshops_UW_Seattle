# Notes can be added like this
library(tidyverse) # or like this!
library(palmerpenguins)
library(doParallel)
library(foreach)

# Save data to object in environment
dat <- penguins
# View data
dat

# Base R

# Step 1: Filter data to the first penguin species
## Remember that R calls data frames with [rows, columns]
dat.filter <- dat[dat$species == "Adelie", ]

# Step 2: Convert the mm measurements to cm
bill_length_cm <- dat.filter$bill_length_mm/10
bill_depth_cm <- dat.filter$bill_depth_mm/10

# Step 3: Plot measurements separated by sex
plot(x = dat.filter$sex, y = bill_length_cm, 
     ylab = "bill length (cm)", xlab = "", main = "Adelie")
plot(x = dat.filter$sex, y = bill_depth_cm, 
     ylab = "bill depth (cm)", xlab = "", main = "Adelie")

# `for` loops
## Syntax

# list the species by-hand
for(i in c("Adelie", "Chinstrap", "Gentoo")){
  print(i)
}

# Or have R figure out the species for you!
for(i in unique(dat$species)){
  print(i)
}

# Use any iter term
for(magic in unique(dat$species)){
  print(magic)
}

## Using a loop
#Make 3 by 2 grid for 6 plots
par(mfrow=c(3,2))

for(i in unique(dat$species)){
  # Step 1: Filter data to the penguin species
  dat.filter <- dat[dat$species == i, ]

  # Step 2: Convert the mm measurements to cm
  bill_length_cm <- dat.filter$bill_length_mm/10
  bill_depth_cm <- dat.filter$bill_depth_mm/10

  # Step 3: Plot measurements separated by sex
  plot(x = dat.filter$sex, y = bill_length_cm, 
      ylab = "bill length (cm)", xlab = "", main = i)
  plot(x = dat.filter$sex, y = bill_depth_cm, 
      ylab = "bill depth (cm)", xlab = "", main = i)
}

# `apply` functions
## Using `apply`
# Will give an error
dat.filter.cm <- apply(dat.filter, MARGIN = 2, function(x) x/10)

# Using index numbers
## shorter code but more likely to result in errors if you're not careful
dat.filter.cm <- apply(dat.filter[c(3,4)], 
                       MARGIN = 2, function(x) x/10)

# Using column names
## longer code but more exact and no errors if your column order changes
dat.filter.cm <- apply(dat.filter[c("bill_length_mm","bill_depth_mm")], 
                       MARGIN = 2, function(x) x/10)


# The tidyverse
## mutate multiple variables

#Convert bill length to cm
dat.mutate <- mutate(dat.filter, bill_length_cm = bill_length_mm/10)

#list columns mutated data
colnames(dat.mutate)

#Convert bill length AND depth to cm
dat.mutate <- mutate(dat.filter, across(c(bill_length_mm, bill_depth_mm),
                                        .fns = ~.x/10))

#list columns in mutated data
colnames(dat.mutate)

dat.filter$bill_length_mm[1:5]
dat.mutate$bill_length_mm[1:5]

#Convert bill length AND depth to cm. Rename columns
dat.mutate <- mutate(dat.filter, across(c(bill_length_mm, bill_depth_mm), 
                                        .fns = ~.x/10, .names = "{.col}_cm"))

#list columns in mutated data
colnames(dat.mutate)

## More tidyverse functions

#Make 3 by 2 grid for 6 plots
par(mfrow=c(3,2))

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

# Nested loops

#Make 3 by 2 grid for 6 plots
par(mfrow=c(3,2))

for(penguin.species in unique(dat$species)){
  for(variable in c("bill_length_mm", "bill_depth_mm")){
    dat.filter <- dat %>% 
      # Step 1: Filter data to the penguin species
      filter(species == penguin.species) %>% 
      # Step 2: Convert the mm measurements to cm
      mutate(variable.cm = get(variable)/10)
    
    ## Make y-label from variable name
    y.label <- gsub("mm", "(cm)", variable)
    y.label <- gsub("_", " ", y.label)
    # Step 3: Plot measurements separated by sex
    plot(x = dat.filter$sex, y = dat.filter$variable.cm, 
         ylab = y.label, xlab = "", main = penguin.species)
  }
}

# foreach loops

registerDoParallel(cores=2)

## Using foreach

foreach(i = unique(dat$species)) %dopar% {
  # Step 1: Filter data to the penguin species
  dat.filter <- dat[dat$species == i, ]

  # Step 2: Convert the mm measurements to cm
  bill_length_cm <- dat.filter$bill_length_mm/10
  bill_depth_cm <- dat.filter$bill_depth_mm/10

  # Step 3: Plot measurements separated by sex
  #Save to PDF named as penguin species
  pdf(paste(i, ".pdf", sep=""))
  #Make 1 by 2 grid for 2 plots
  par(mfrow=c(1,2))
  plot(x = dat.filter$sex, y = bill_length_cm,
      ylab = "bill length (cm)", xlab = "", main = i)
  plot(x = dat.filter$sex, y = bill_depth_cm, 
       ylab = "bill depth (cm)", xlab = "", main = i)
  dev.off()
}