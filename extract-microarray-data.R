folder <- "J:\\MacLabUsers\\Claire\\Projects\\HVE-microarray\\microarrayData\\"

file <- readLines(paste0(folder, "FinalReport_exvivo_TNF_HVE.txt"))

# File structure
# Line 1:      [Header]
# 2-7:         Header details

# 8:           [Sample Probe Profile]
# 9:           Column names
# 10-47332:    Microarray data

# 47333:       [Control Gene Profile]
# 47334:       Column names
# 47335-432:   Control gene data

# 47433:       [Excluded and Imputed Probes]
# 47434:       Column names
# 47435-45:     Excluded and imputed data

# 47446:       [Samples Table]
# 47447:       Column names
# 47448-47483: Sample table data

# # save just microarray data as separate file
# fileConn <- file(paste0(folder, "MicroarrayDataExtracted.txt"))
# writeLines(file[1:47332], fileConn)
# close(fileConn)

# get control gene profiles 
ctrl <- read.table(paste0(folder, "FinalReport_exvivo_TNF_HVE.txt"), 
   sep = "\t", skip = 47333, nrows = 98, header = TRUE)

# get excluded and imputed probes
# the following skips over several lines between header and the first gene in
#     the table for some reason????
excludedImputed <- read.table(paste0(folder, "FinalReport_exvivo_TNF_HVE.txt"), 
   sep = "\t", skip = 47433, nrows = 7, header = TRUE)

# get sample table
samples <- read.table(paste0(folder, "FinalReport_exvivo_TNF_HVE.txt"), 
   sep = "\t", skip = 47446, nrows = 36, header = TRUE)
