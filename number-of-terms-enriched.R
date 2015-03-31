# Number of GO terms overlapping between days 4, 7, and 14
require(dplyr)
require(reshape2)
require(stringr)
require(ggplot2)
require(pander)
source("getAllDavid.R")
source("getAllInnateDB.R")
source('subsetToOverlappingGoTerms.R')

###############################################################################
# How many GO terms are enriched on days 4, 7, and 14 in the DAVID database? 
###############################################################################
# get data and subset to only terms common to days 4, 7, and 14
david <- getAllDavid(50) %>%
   mutate(Concentration = 50)
   
david <- rbind(
   subsetToOverlappingGoTerms(david, direction = "UP", withDay1 = FALSE),
   subsetToOverlappingGoTerms(david, direction = "DOWN", withDay1 = FALSE)
)

temp <- getAllDavid(500) %>%
   mutate(Concentration = 500)

david <- rbind(david,
   subsetToOverlappingGoTerms(temp, direction = "UP", withDay1 = FALSE),
   subsetToOverlappingGoTerms(temp, direction = "DOWN", withDay1 = FALSE)
)

rm(temp)

# add column indicating whether p < 0.05 for all 3 days and replace Term
# column with GO ids where possible
# note: removing 500 down because data are bad
david <- david %>%
   mutate(Term = ifelse(str_detect(Term,"GO:"), substr(Term,1,10), Term)) %>% 
   filter(!(direction == "DOWN" & Concentration == 500)) %>% # remove once data is good
   group_by(Concentration, direction, Term) %>%
   summarize(Significant = all(Benjamini < 0.05)) %>%
   ungroup()

# generate summary table of DAVID results
davidResult <- david %>%
   group_by(Concentration, direction) %>%
   summarize(NumberOfTerms = length(unique(Term)), 
      Significant = sum(Significant),
      Which = "DAVID", 
      onlyGO = "*contains non-GO terms")  

###############################################################################
# How many GO terms are enriched on days 4, 7, and 14 in the InnateDB database? 
###############################################################################
# Copy the terms and p values into new columns for the InnateDB data. Columns
# need these names for subsetToOverlappingGoTerms() to work
innate <- getAllInnateDB(50) %>%
   mutate(Concentration = 50, Term = Pathway.Id, PValue = Pathway.p.value.corrected)

innate <- rbind(
   subsetToOverlappingGoTerms(innate, direction = "UP", withDay1 = FALSE),
   subsetToOverlappingGoTerms(innate, direction = "DOWN", withDay1 = FALSE)
)

temp <- getAllInnateDB(500) %>%
   mutate(Concentration = 500, Term = Pathway.Id, PValue = Pathway.p.value.corrected)

innate <- rbind(innate,
   subsetToOverlappingGoTerms(temp, direction = "UP", withDay1 = FALSE),
   subsetToOverlappingGoTerms(temp, direction = "DOWN", withDay1 = FALSE)
)

rm(temp)

# add column indicating whether p < 0.05 for all 3 days 
innate <- innate %>%
   group_by(Concentration, direction, Term) %>%
   summarize(Significant = all(Pathway.p.value.corrected < 0.05)) %>%
   ungroup()

# generate summary table of InnateDB results
innateResult <- innate %>%
   group_by(Concentration, direction) %>%
   summarize(NumberOfTerms = length(unique(Term)), 
      Significant = sum(Significant),
      Which = "InnateDB",
      onlyGO = "only GO terms")  

###############################################################################
# Combined summary tables 
###############################################################################
result <- rbind(davidResult, innateResult)

pander(result)

###############################################################################
# Which terms are enriched in both DAVID and InnateDB?
###############################################################################
bothDatabases <- rbind(
   mutate(david, which = "DAVID"), 
   mutate(innate, which = "INNATE")
)

bothDatabases %>%
   filter(!(direction == "DOWN" & Concentration == 500)) %>% # remove once DAVID data is good
   group_by(Concentration, direction, Term) %>%
   summarize(OneSignificant = sum(Significant) > 0,
      AllSignificant = sum(Significant) == 2) %>%
   group_by(Concentration, direction) %>%
   summarize(NumberOfTerms = length(unique(Term)), 
      DavidAndInnateDB = sum(AllSignificant),
      DavidOrInnateDB = sum(OneSignificant),
      Which = "Both") 
