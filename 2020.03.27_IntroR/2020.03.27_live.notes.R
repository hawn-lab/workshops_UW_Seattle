#### Setup packages ####
#Only need to install once per computer
#Automatically installs from CRAN https://cran.r-project.org/
#install.packages("ggplot2")

#Other biology packages available at Bioconductor.org

#load package
#Must be done everytime you open R/RStudio
library(ggplot2)

#### Load data ####
#Make data folder and put the data file in it
#Basic R function syntax
# functionName(argument1=1, argument2=TRUE)
# can use spaces for improve readability
# functionName ( argument1 = 1 , argument2 = TRUE )

dat <- read.table(file = "data/AM.MDM.data.cleaning.metrics.csv",
           sep=",", header=TRUE)

#Explore data
dim(dat)
class(dat)

#Extract column
dat$sampID

#Extract row
#format [row , column]
dat[1,]
dat[1,1]

#What types are these data?
class(dat$sampID)
class(dat$raw)
class(dat$PCT_CODING_BASES)

#not seen type
#logical =  a TRUE/FALSE variable

#### Doing math ####
#Calculate the variance of a variable
var(dat$raw)
#be careful of type! variance only works on numeric data
var(dat$sampID)

#is x TRUE for each value of vector?
#example: which values of raw are greater than 1 billion?
dat$raw > 1E9

#transformation
dat$PCT_CODING_BASES *100

#R automatically repeats a vector to make it the same length as the other vector
dat$PCT_CODING_BASES * c(100,1000)
dat$PCT_CODING_BASES[1] * c(100,1000,10000)

#unique / levels
unique(dat$sampID)
levels(dat$sampID)

#### help function ####
?var
## if don't function name
??variance

#### Subset data ####
#Subset 1 dimensional data like a vector
dim(dat$sampID)
length(dat$sampID)
#First value
dat$sampID[1]
dat$raw[1]
#2 through 5
dat$raw[2:5]
#2 through 50 - note that R adds NA when there is no data past 24
dat$raw[2:50]

#Subset 2 dimensions
#[row, column]
dat[1:2, 1:2]

#Subset 2 dimensions by some "rule"
#Get all rows where raw is greater than 100 million
logical.vector <- dat$raw > 100E6
logical.vector

dat[logical.vector, 1:3]

#Or use column name in subsetting
dat[logical.vector, "align.paired"]

#### EXERCISES ####
# 2. 
log(4)
log(4,2)
#OR
log2(4)
log(x=4,base=4)
# 3. 
class(dat$paired)
# 4. a. 
dat$trim[20]
# 4. b. 
dat$raw == 95004980
dat$raw[22]
dat[22, "raw"]

dat[dat$raw == 95004980, "raw"]

#5. 
cats <- dat$Assigned >= 20E6
dat[cats, "Assigned"]

dat[dat$Assigned >= 20E6 & 
    dat$raw < 80E6, 
    1:3]

##### Case sensitive #####
J <- "J"
j <- "j"

J == j

#### Tidyverse ####

#load package
library(tidyverse)

#load data
dat <- read_csv(file = "data/AM.MDM.data.cleaning.metrics.csv")
spec(dat)

#force sampID column to be a factor
dat <- read_csv(file = "data/AM.MDM.data.cleaning.metrics.csv",
                col_types = cols(sampID = col_factor()))

#Other options: read_tsv or read_table

#### dplyr ####
#Copy data so have an unaltered raw copy
raw_dat <- dat

## Select columns
colnames(dat)
dat <- select(.data = dat,
              sampID, raw, trim,
              both.align.paired, 
              both.align.paired_filter,
              Assigned_paired)

## Filter rows
dat <- filter(.data = dat,
              grepl(pattern = "Media", sampID))

#base comparison
#dat[ grepl(pattern = "Media", dat$sampID), ]

#### Pipe ####
#### %>%  ####

dat <- raw_dat %>% 
  select(sampID, raw, trim,
         both.align.paired, 
         both.align.paired_filter,
         Assigned_paired) %>% 
  filter(grepl(pattern = "Media", sampID))

dat <- read_csv("data/AM.MDM.data.cleaning.metrics.csv") %>% 
  select(sampID, raw, trim,
         both.align.paired, 
         both.align.paired_filter,
         Assigned_paired) %>% 
  filter(grepl(pattern = "Media", sampID))

#### tidyr ####
# pivot_longer
# pivot_wider
# Transform from wide to long format, or vice versa

dat_long <- dat %>% 
  pivot_longer(cols = -sampID,
               names_to = "name",
               values_to = "value")

dat_wide <- dat_long %>% 
  pivot_wider(names_from = "name",
              values_from = "value")

#### ggplot ####
# Data
# aesthetics
# geoms - dot plot, bar plot...

dat %>% 
  ggplot()

#is the same as

ggplot(dat)

## geom_point = dot plot

dat %>% 
  ggplot(aes(x = sampID, y = raw)) +
  geom_point()

#### Fancy plot ####
plot1 <- 
  #### Data manipulation #### 
raw_dat %>%
  #Select variables of interest
  select(sampID, raw, trim, 
         both.align.paired,
         both.align.paired_filter,
         Assigned_paired) %>% 
  #Filter to Media samples only
  filter(grepl("Media",sampID)) %>% 
  #Convert data to long format
  pivot_longer(-sampID, names_to = "group", values_to = "sequences") %>% 
  #Convert group into a factor and force the level order
  mutate(group = factor(group, levels=c("raw", "trim", "both.align.paired",
                                        "both.align.paired_filter", 
                                        "Assigned_paired"))) %>% 
  
  #### Basic plot  #### 
#Initiate ggplot of sequences in each sample, 
ggplot(aes(x=sampID,y=sequences, 
           #Fill bar color by the variable group
           #Also tell which variable to group bars by
           fill=group)) +
  #Create a bar plot, grouping all groups from each sample together
  geom_bar(stat="identity", position=position_dodge()) +
  
  #### Customization ####
#Add a theme to change the overall look of the plot
  theme_classic() +
  #Rotate x-axis labels so we can read them
  theme(axis.text.x = element_text(angle = 90)) +
  #Re-label the legend
  scale_fill_discrete(name="", labels=c("Raw",
                                        "Min quality, adapter trimmed",
                                        "Aligned, paired",
                                        "High-quality aligned, paired",
                                        "Assigned to gene")) +
  #Re-label the axes
  labs(x="",y="Total sequences") +
  #add a title
  ggtitle("Sequences retained during RNA-seq clean-up")

#### save plot ####
ggsave(filename = "plot1.png",
       plot = plot1, 
       width=5, height=5)

#### END ####