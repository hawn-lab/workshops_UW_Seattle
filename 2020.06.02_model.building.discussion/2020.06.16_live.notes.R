#STATISTICAL MODELS CON'T

##### Interaction term vs contrast matrix vs delta outcome ####

#Interaction term
# Overall difference between any potential contrast within RSTR:infection groups
`expression ~ infection * RSTR`
## infection:RSTR term asks if the slope from media to TB WITHIN RSTR is different from the slope from media to TB WITHIN LTBI
## does gene expression in response to TB infection change differently in RSTR vs LTBI?
## It is more statistically difficult to say two slopes are different (b/c of variance within the sample) than it is to say 1 slope is different from 0 (b/c 0 has not variance)

"Pros
  * Simple to code
  * Simple interpretation with only 1 interaction

Cons
  * Statistically difficult for interaction term to reach significance
  * Lower power within interaction term than main terms
  * Cannot tell WHICH contrast is significant
  * Interpretation becomes difficult with multiple interaction terms
"

#Contrasts
#lm
# New variable that pastes the RSTR status with infection status like RSTR_media, RSTR_TB...
`expression ~ RSTR_infection`

#limma
make.contrasts( )

#other
#Subset data by variable 1 and run a model on each subset

## contrast of RSTR_media vs LTBI_media asks if the slope from RSTR to LTBI WITHIN media samples if different from 0
## contrast of RSTR_media vs RSTR_TB asks if the slope from media to TB WITHIN RSTR is different from 0
## Remember comparing to 0 is statistically easier

"Pros
  * Simple interpretation
  * Only look at pairwise comparisons of interest
  * Can tell WHICH contrast(s) are significant

Cons
  * More difficult to code - so far no easy use in lmekin or eQTL
  * Lower power in subsets compared to main terms in interaction model
"

# Delta model
#Set the outcome (y-variable) as a change in expression in response to one of the variables
`TB minus media expression ~ RSTR`

## is (TB-media in RSTR) minus (TB-media in LTBI) different from 0?
## is the magnitude of TB - media different between LTBI and RSTR? 
## does gene expression in response to TB infection change differently in RSTR vs LTBI?

"Pros
  * focuses on 2 contrasts of interest

Cons
  * but ONLY 2 contrasts of the 4 we're interested in
"

##### Picking a model from above #####
"
Interaction term: you're interested in the overall effect of both terms including all pairwise comparisons

Contrast: you're only interested in specific pairwise comparisons due to biological relevance or experimental setup

Delta: you're only interested in specific pairwise comparisons that can be coded by 1 main model term

Multiple: if you want a more lenient list to go into a downstream analysis (WGCNA, GSEA, etc), maybe to more than one and compare
"

##### Getting contrasts to work outside limma ######
"
1. Subsetting 
    a. Subset to media or TB samples only. Run model of `outcome ~ RSTR` in 2 subsets. 
    b. Subset to RSTR or LTBI samples only. Run model of `outcome ~ infection` in 2 subsets. 
    c. FDR correct across ALL p-values from 4 models. 
2. Other TBD...
"
