## working with salivary ACC

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

project <- "FFPE-Normals"
outDir <- "FFPE-Normals"

## cases to include
ffpeInc <- c("N2", "N7", "N8", "N11", "N12", "N15")
defuse <- read.table("merged.defuse.results.txt", fill = TRUE, header = TRUE, sep = "\t")
names(defuse)[1:3] <- c("program", "case", "outputFile")
defuse$fusionName <- paste(defuse$gene_name1, defuse$gene_name2, sep = "--")

## extract FFPE from individual filtered fusions
ffpe.star <- star[star$case %in% ffpeInc, ]
ffpe.map <- map[map$case %in% ffpeInc, ]
ffpe.ee <- ee[ee$case %in% ffpeInc, ]
ffpe.integrate <- integrate[integrate$case %in% ffpeInc, ]
ffpe.defuse <- defuse[defuse$case %in% ffpeInc, ]

write.table(ffpe.star, file = "FFPE-Normals/FFPE_star-fusion.fusion_candidates.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ffpe.map, file = "FFPE-Normals/FFPE_mapsplice.fusion_candidates.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ffpe.ee, file = "FFPE-Normals/FFPE_ericscript.results.filtered.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ffpe.integrate, file = "FFPE-Normals/FFPE_integrate.oncofuse.output.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ffpe.defuse, file = "FFPE-Normals/FFPE_defuse.output.txt", quote = FALSE, row.names = FALSE, sep = "\t")

ffpeFilter <- unique(sort(c(ffpe.star$fusionName, ffpe.map$fusionName, ffpe.ee$fusionName, ffpe.integrate$fusionName, ffpe.defuse$fusionName)))
write(ffpeFilter, file = "FFPE-Normals/ffpeFilter.txt")
