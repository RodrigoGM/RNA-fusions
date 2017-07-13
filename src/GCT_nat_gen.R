## working with GCT

library(ggplot2)
library(reshape2)
library(naturalsort)
source("myLib.R")

load("fusions.data.rda")
# load("GCT/fusions.gct.rda")

project <- "Granular Cell Tumor"
outDir <- "GCT"

## cases to include
gctInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project]))

## subset integrate and defuse
gct.integrate <- integrate[integrate$case %in% gctInc, ]
gct.defuse <- defuse[defuse$case %in% gctInc, ]

tcgaInt <- readLines("TCGA-Normals/tcga_brca_normals_integrate.txt")
tcgaDef <- readLines("TCGA-Normals/tcga_brca_normals_defuse.txt")
ffpeInt <- readLines("FFPE-Normals/ffpe_normals_integrate.txt")
ffpeDef <- readLines("FFPE-Normals/ffpe_normals_defuse.txt")

gct.integrate$inTCGA.Normals <- gct.integrate$fusionName %in% c(tcgaInt, tcgaDef)
gct.integrate$inFFPE.Normals <- gct.integrate$fusionName %in% c(ffpeInt, ffpeDef)

gct.defuse$inTCGA.Normals <- gct.defuse$fusionName %in% c(tcgaInt, tcgaDef)
gct.defuse$inFFPE.Normals <- gct.defuse$fusionName %in% c(ffpeInt, ffpeDef)

write.table(gct.integrate, file = "GCT/GCT_NatGen_integrate_fusions.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(gct.defuse, file = "GCT/GCT_NatGen_defuse_fusions.txt", quote = FALSE, sep = "\t", row.names = FALSE)

gct.int.filtered <- gct.integrate[! gct.integrate$inTCGA.Normals, ]
gct.def.filtered <- gct.defuse[! gct.defuse$inTCGA.Normals, ]

oncofuseCoord <- gct.def.filtered[, c("gene_chromosome1", "genomic_break_pos1", "gene_align_strand1", "gene_chromosome2", "genomic_break_pos2", "gene_align_strand2")]

oncofuseCoord$genomic_break_pos1[oncofuseCoord$gene_align_strand1 == "+"] <- oncofuseCoord$genomic_break_pos1[oncofuseCoord$gene_align_strand1 == "+"] + 1
oncofuseCoord$genomic_break_pos1[oncofuseCoord$gene_align_strand1 == "-"] <- oncofuseCoord$genomic_break_pos1[oncofuseCoord$gene_align_strand1 == "-"] - 1
oncofuseCoord$gene_align_strand2[oncofuseCoord$gene_align_strand2 == "+"] <- oncofuseCoord$genomic_break_pos2[oncofuseCoord$gene_align_strand2 == "+"] - 1
oncofuseCoord$gene_align_strand2[oncofuseCoord$gene_align_strand2 == "-"] <- oncofuseCoord$genomic_break_pos2[oncofuseCoord$gene_align_strand2 == "-"] + 1
oncofuseCoord$tissue <- "EPI"

oncofuseCoord$gene_chromosome1 <- paste("chr", oncofuseCoord$gene_chromosome1, sep = "")
oncofuseCoord$gene_chromosome2 <- paste("chr", oncofuseCoord$gene_chromosome2, sep = "")

write.table(oncofuseCoord[, c("gene_chromosome1", "genomic_break_pos1", "gene_chromosome2", "genomic_break_pos2", "tissue") ], file = "GCT/gct_defuse.oncofuse.input.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
