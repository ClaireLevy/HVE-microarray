library(dplyr)
library(stringr)
library(ggplot2)
source("DEG-to-long-form.R") # get "longForm" data frame
source("collapse-probes.R") # get function

path <- "C:\\Users\\smhughes\\Desktop\\"

############################################################################
# day 1 50 uM heat map, only genes that were affected
############################################################################
# get target IDs for genes that are up or down on Day 1 at 50 uM
longForm %>%
   filter(Concentration == 50, Day == 1, Direction %in% c("UP", "DOWN")) %>%
   select(TargetID)  %>% 
   subset(select = 1, drop = TRUE) -> genes

# get data for all those target IDs and and collapse multiple probes per gene
longForm %>%
   group_by(Direction, Concentration, Day) %>%
   collapseProbes() %>%
   filter(TargetID %in% genes, Direction %in% c("UP", "DOWN", "FALSE")) -> d1Low

# sort target IDs by direction
d1Low %>%
   filter(Day == 1, Concentration == 50) %>%
   arrange(Direction) -> targets
targets <- targets[["TargetID"]]

d1Low %>%
   mutate(TargetID = factor(TargetID, levels = as.character(targets))) -> toPlot

p <- ggplot(toPlot, aes(x = factor(Day), y = TargetID)) 
p + geom_tile(aes(fill = Direction)) + facet_wrap( ~ Concentration) +
   scale_fill_brewer(type = "div", palette = 5) +
   ggtitle(paste("Genes differentially regulated on day 1 50 μM")) +
   xlab("Day") + ylab("Gene")

ggsave(filename = paste0(path, "day1.pdf"),
   height = 11, width = 8.5, units = "in")

############################################################################
# all genes affected on any day sorted sequentially
############################################################################

## gene level
x <- collapseProbes(longForm)

# omit genes where probes didn't agree
d <- filter(x, Direction %in% c("UP", "DOWN", "FALSE")) %>%
   ungroup()

# put genes in order
genesInOrder <- d %>%
   filter(Concentration == 50) %>%
   # add leading 0 to days < 10 for better ordering (01, 04, 14 not 1, 14, 4)
   mutate(Day = stringr::str_pad(Day, 2, pad = "0")) %>%
   # collapse to one row per target by targetID
   group_by(TargetID) %>%
   summarize(n = sum(Direction != "FALSE"),
      up = sum(Direction == "UP"),
      down = sum(Direction == "DOWN"),
      which = paste(Day[Direction != "FALSE"], collapse = "."), 
      Direction = ifelse(up > 0 & down > 0, "BOTH", 
         ifelse(up > 0, "UP", 
            ifelse(down > 0, "DOWN", "FALSE")))) %>%
   # put in order of the new variables
   arrange(Direction, n, which) %>%
   # extract TargetID column as vector
   .$TargetID

# arrange data according to genesInOrder
toPlot <- d %>%
   mutate(TargetID = factor(TargetID, levels = as.character(genesInOrder)))

# plot
p <- ggplot(toPlot, aes(x = factor(Day), y = TargetID))
p + geom_tile(aes(fill = Direction)) + facet_wrap( ~ Concentration) +
   scale_fill_brewer(type = "div", palette = 5) + 
   scale_y_discrete(labels = "") +
   ggtitle("Genes affected by tenofovir in HVE cells") +
   xlab("Day") + ylab("Gene")