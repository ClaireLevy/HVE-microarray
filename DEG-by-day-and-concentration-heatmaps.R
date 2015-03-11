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
# all genes affected on any day
############################################################################

## heatmaps
d <- collapseProbes(longForm)

# remove TargetIDs where probes don't agree
d <- filter(d, Direction %in% c("DOWN", "UP", "FALSE"))

heatmapBy <- function(toPlot, d, conc) {
   toPlot %>%
      filter_(~Day == d, ~Concentration == conc) %>%
      arrange(Direction) %>%
      select(TargetID) -> z
   z <- z[["TargetID"]]   
   
   toPlot %>%
      mutate(TargetID = factor(TargetID, levels = as.character(z))) -> toPlot
   
   p <- ggplot(toPlot, aes(x = factor(Day), y = TargetID))
   p + geom_tile(aes(fill = Direction)) + facet_wrap( ~ Concentration) +
      scale_fill_brewer(type = "div", palette = 5) +
      scale_y_discrete(labels = "") +
      ggtitle(paste("Sorted by day", d, conc, "μM")) +
      xlab("Day") + ylab("Gene")
         
}

for (conc in c(50, 500)) {
   for (day in c(1, 4, 7, 14)) {
      heatmapBy(d, day, conc)
      ggsave(filename = paste0(path, day, ".", conc, ".pdf"),
         height = 11, width = 8.5, units = "in")
   }
}


heatmapBy(d, 1, 50)
heatmapBy(d, 4, 50)
heatmapBy(d, 7, 50)
heatmapBy(d, 14, 50)
heatmapBy(d, 1, 500)
heatmapBy(d, 4, 500)
heatmapBy(d, 7, 500)
heatmapBy(d, 14, 500)