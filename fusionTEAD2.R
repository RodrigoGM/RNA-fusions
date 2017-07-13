## working with salivary ACC

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

project <- "TEAD2fs"
outDir <- "TEAD2"

if (! file.exists(outDir)) {
    dir.create(outDir)
}

## cases to include
tead2Inc <-  unique(c(star$case[star$project == project],
                    map$case[map$project == project],
                    ee$case[ee$project == project],
                    integrate$case[integrate$project == project]))

## extract TEAD2 unfiltered
tead2.star <- star[star$case %in% tead2Inc,]
tead2.map <- map[map$case %in% tead2Inc,]
tead2.ee <- ee[ee$case %in% tead2Inc,]
tead2.integrate <- integrate[integrate$case %in% tead2Inc,]
tead2.defuse <- defuse[defuse$case %in% tead2Inc,]

## extract TEAD2 from individual filtered fusions
tead2.star.filtered <- star.filtered[star.filtered$project == project, ]
tead2.star.filtered <- tead2.star.filtered[!is.na(tead2.star.filtered$case),]
tead2.star.filtered <- tead2.star.filtered[tead2.star.filtered$case %in% tead2Inc, ]

tead2.map.filtered <- map.filtered[map.filtered$project == project, ]
tead2.map.filtered <- tead2.map.filtered[!is.na(tead2.map.filtered$case),]
tead2.map.filtered <- tead2.map.filtered[tead2.map.filtered$case %in% tead2Inc, ]

tead2.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
tead2.ee.filtered <- tead2.ee.filtered[!is.na(tead2.ee.filtered$case),]
tead2.ee.filtered <- tead2.ee.filtered[tead2.ee.filtered$case %in% tead2Inc, ]

tead2.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
tead2.integrate.filtered <- tead2.integrate.filtered[!is.na(tead2.integrate.filtered$case),]
tead2.integrate.filtered <- tead2.integrate.filtered[tead2.integrate.filtered$case %in% tead2Inc, ]

tead2.defuse.filtered <- defuse.filtered[defuse.filtered$project == project, ]
tead2.defuse.filtered <- tead2.defuse.filtered[!is.na(tead2.defuse.filtered$case),]
tead2.defuse.filtered <- tead2.defuse.filtered[tead2.defuse.filtered$case %in% tead2Inc, ]

## extract TEAD2 from recurrent filtered fusions
tead2.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
tead2.star.recurrent <- tead2.star.recurrent[!is.na(tead2.star.recurrent$case),]
tead2.star.recurrent <- tead2.star.recurrent[tead2.star.recurrent$case %in% tead2Inc, ]

tead2.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
tead2.map.recurrent <- tead2.map.recurrent[!is.na(tead2.map.recurrent$case),]
tead2.map.recurrent <- tead2.map.recurrent[tead2.map.recurrent$case %in% tead2Inc, ]

tead2.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
tead2.ee.recurrent <- tead2.ee.recurrent[!is.na(tead2.ee.recurrent$case),]
tead2.ee.recurrent <- tead2.ee.recurrent[tead2.ee.recurrent$case %in% tead2Inc, ]

tead2.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
tead2.ifr.recurrent <- tead2.ifr.recurrent[!is.na(tead2.ifr.recurrent$case),]
tead2.ifr.recurrent <- tead2.ifr.recurrent[tead2.ifr.recurrent$case %in% tead2Inc, ]

tead2.def.recurrent <- def.recurrent[def.recurrent$project == project, ]
tead2.def.recurrent <- tead2.def.recurrent[!is.na(tead2.def.recurrent$case),]
tead2.def.recurrent <- tead2.def.recurrent[tead2.def.recurrent$case %in% tead2Inc, ]


## pair-wise interesections
tead2SM <- merge(tead2.star, tead2.map, by = c("fusionName", "case"))
tead2SI <- merge(tead2.star, tead2.integrate, by = c("fusionName", "case"))
tead2SE <- merge(tead2.star, tead2.ee, by = c("fusionName", "case"))
tead2SD <- merge(tead2.star, tead2.defuse, by = c("fusionName", "case"))
tead2ME <- merge(tead2.map, tead2.ee, by = c("fusionName", "case"))
tead2MI <- merge(tead2.map, tead2.integrate, by = c("fusionName", "case"))
tead2MD <- merge(tead2.map, tead2.defuse, by = c("fusionName", "case"))
tead2EI <- merge(tead2.ee, tead2.integrate, by = c("fusionName", "case"))
tead2DI <- merge(tead2.defuse, tead2.integrate, by = c("fusionName", "case"))
tead2DE <- merge(tead2.defuse, tead2.ee, by = c("fusionName", "case"))

## three-way interesections
tead2MSE <- merge(tead2SM, tead2.ee, by = c("fusionName", "case"))
tead2MSI <- merge(tead2SM, tead2.integrate, by = c("fusionName", "case"))
tead2MSD <- merge(tead2SM, tead2.defuse, by = c("fusionName", "case"))
tead2MEI <- merge(tead2ME, tead2.integrate, by = c("fusionName", "case"))
tead2MED <- merge(tead2ME, tead2.defuse, by = c("fusionName", "case"))
tead2SEI <- merge(tead2SE, tead2.integrate, by = c("fusionName", "case"))
tead2SED <- merge(tead2SE, tead2.defuse, by = c("fusionName", "case"))
tead2MID <- merge(tead2MI, tead2.defuse, by = c("fusionName", "case"))
tead2IED <- merge(tead2EI, tead2.defuse, by = c("fusionName", "case"))

## all programs
#tead24M <- merge(tead2SM, tead2EI, by = c("fusionName", "case"))

## write individual intersections
##write.table(tead2SM, file = "TEAD2/TEAD2_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(tead2SE, file = "TEAD2/TEAD2_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(tead2ME, file = "TEAD2/TEAD2_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(tead2MSE, file = "TEAD2/TEAD2_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(tead2MSI, file = "TEAD2/TEAD2_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(tead2MEI, file = "TEAD2/TEAD2_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(tead2SEI, file = "TEAD2/TEAD2_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(tead24M, file = "TEAD2/TEAD2_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(tead2.star.recurrent, file = "TEAD2/TEAD2_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.map.recurrent, file = "TEAD2/TEAD2_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.ee.recurrent, file = "TEAD2/TEAD2_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.ifr.recurrent, file = "TEAD2/TEAD2_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.def.recurrent, file = "TEAD2/TEAD2_recurrent_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write Unfiltered
write.table(tead2.star, file = "TEAD2/TEAD2_unfiltered_fusions_STAR-Fusion.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.map, file = "TEAD2/TEAD2_unfiltered_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.ee, file = "TEAD2/TEAD2_unfiltered_fusions_Ericscript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.integrate, file = "TEAD2/TEAD2_unfiltered_fusions_INTEGRATE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.defuse, file = "TEAD2/TEAD2_unfiltered_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

unfilt <- merge(tead2.star, tead2.map, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, tead2.ee, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, tead2.integrate, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, tead2.defuse, by = c("fusionName", "case"), all = TRUE)

write.table(unfilt, file = "TEAD2/TEAD2_unfiltered_fusions_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## count individual filtered by programs
tead2SCount <- countFusions(tead2.star.filtered, "fusionName", "case")
tead2MCount <- countFusions(tead2.map.filtered, "fusionName", "case")
tead2ECount <- countFusions(tead2.ee.filtered, "fusionName", "case")
tead2ICount <- countFusions(tead2.integrate.filtered, "fusionName", "case")
tead2DCount <- countFusions(tead2.defuse.filtered, "fusionName", "case")
## count recurrent filtered by programs
tead2SfrCount <- countFusions(tead2.star.recurrent, "fusionName", "case")
tead2MfrCount <- countFusions(tead2.map.recurrent, "fusionName", "case")
tead2EfrCount <- countFusions(tead2.ee.recurrent, "fusionName", "case")
tead2IfrCount <- countFusions(tead2.ifr.recurrent, "fusionName", "case")
tead2DfrCount <- countFusions(tead2.def.recurrent, "fusionName", "case")

## count intersections
tead2SMCounts <- countFusions(tead2SM, "fusionName", "case")
tead2SECounts <- countFusions(tead2SE, "fusionName", "case")
tead2SICounts <- countFusions(tead2SI, "fusionName", "case")
tead2SDCounts <- countFusions(tead2SD, "fusionName", "case")
tead2MECounts <- countFusions(tead2ME, "fusionName", "case")
tead2MICounts <- countFusions(tead2MI, "fusionName", "case")
tead2MDCounts <- countFusions(tead2MD, "fusionName", "case")
tead2EICounts <- countFusions(tead2EI, "fusionName", "case")
tead2DICounts <- countFusions(tead2DI, "fusionName", "case")
tead2DECounts <- countFusions(tead2DE, "fusionName", "case")


##tead2Counts <- countFusions(tead24M, "fusionName", "case")


fusionColors <- c("Not Present" = "#CFCFCF" , "Predicted" = "#F8766D")


pdf("TEAD2/TEAD2_byProgram.pdf", width = 210/25.4, height = 420/25.4)
ggplot(tead2SCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2MCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2ECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2ICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2DCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("DeFuse predicted") +
    scale_fill_manual(values = fusionColors, name = "")
dev.off()



pdf("TEAD2/TEAD2_Recurrent_byProgram.pdf", width = 210/25.4, height = 420/25.4)
ggplot(tead2SfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2MfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2EfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2IfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2DfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("DeFuse predicted") + 
    scale_fill_manual(values = fusionColors, name = "")
dev.off()


pdf("TEAD2/TEAD2_predictedFusions.pdf", width = 210/25.4, height = 420/25.4)
ggplot(tead2SMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2SECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2SICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2SDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & DeFuse Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2MDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & DeFuse Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2MECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2MICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2EICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("EricScript & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2DICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(tead2DECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors, name = "")
dev.off()

##ggplot(tead2Counts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
##    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
##    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
##    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
##                      labels = c("Not Present", "Predicted"))


tead2.MIN2 <- unique(do.call(rbind, lapply(list(tead2SM, tead2SI, tead2SE, tead2SD, tead2ME, tead2MI, tead2MD, tead2EI, tead2DI, tead2DE), function(X) X[,c("fusionName", "case")])))
tead2.MIN2 <- merge(tead2.MIN2, tead2.ee, all.x = TRUE, by = c("fusionName"))
tead2.MIN2 <- merge(tead2.MIN2, tead2.map, all.x = TRUE, by = c("fusionName"))
tead2.MIN2 <- merge(tead2.MIN2, tead2.star, all.x = TRUE, by = c("fusionName"))
tead2.MIN2 <- merge(tead2.MIN2, tead2.integrate, all.x = TRUE, by = c("fusionName"))
tead2.MIN2 <- merge(tead2.MIN2, tead2.defuse, all.x = TRUE, by = c("fusionName"))
##tead2.MIN2$fusionName <- droplevels(tead2.MIN2$fusionName)

tead2.MIN3 <- unique(do.call(rbind, lapply(list(tead2MSE, tead2MSI, tead2MSD, tead2MEI, tead2MED, tead2SEI, tead2SED, tead2MID, tead2IED), function(X) X[,c("fusionName", "case")])))
tead2.MIN3 <- merge(tead2.MIN3, tead2.ee, all.x = TRUE, by = c("fusionName"))
tead2.MIN3 <- merge(tead2.MIN3, tead2.map, all.x = TRUE, by = c("fusionName"))
tead2.MIN3 <- merge(tead2.MIN3, tead2.star, all.x = TRUE, by = c("fusionName"))
tead2.MIN3 <- merge(tead2.MIN3, tead2.integrate, all.x = TRUE, by = c("fusionName"))
tead2.MIN3 <- merge(tead2.MIN3, tead2.defuse, all.x = TRUE, by = c("fusionName"))
## tead2.MIN3$fusionName <- droplevels(tead2.MIN3$fusionName)

min2Count <- countFusions(tead2.MIN2, caseID = "case.x")
min3Count <- countFusions(tead2.MIN3, caseID = "case.x")
#min4Count <- countFusions(tead2.MIN4, caseID = "case.x")


pdf("TEAD2/TEAD2_MIN2_MIN3.pdf", width = 210/25.4, height = 420/25.4)
ggplot(min2Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10), legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least two tools") +
    scale_fill_manual(values = fusionColors, name = "")
ggplot(min3Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least three tools") +
    scale_fill_manual(values = fusionColors, name = "")
dev.off()


write.table(tead2.MIN2, file = "TEAD2/TEAD2_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tead2.MIN3, file = "TEAD2/TEAD2_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##write.table(tead24M, file = "TEAD2/TEAD2_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("TEAD2/fusions.tead2.rda")

load("TEAD2/fusions.tead2.rda")

## java -jar /opt/src/oncofuse-1.1.1/Oncofuse.jar -a hg19 TEAD2/tead2_min2_oncofuse.input.txt 'coord' EPI TEAD2/tead2_min2_oncofuse.output.txt

Onc <- do.call(rbind, lapply(list.files(pattern = "merged.*oncofuse.output.txt"), read.delim, header = TRUE))
Onc$fusionName <- paste(Onc$X5_FPG_GENE_NAME, Onc$X3_FPG_GENE_NAME, "--")

tead2Onc <- rbind(tead2.ee.recurrent[, c("program", "case", "fusionName")], tead2.map.recurrent[, c("program", "case", "fusionName")], tead2.star.recurrent[, c("program", "case", "fusionName")], tead2.ifr.recurrent[, c("program", "case", "fusionName")])

tead2Onc <- merge(tead2Onc, Onc, by = "fusionName")

pdf("TEAD2/TEAD2_Fusions_DrivProbability.pdf", width = 210/25.4, height = 420/25.4)
ggplot(tead2Onc, aes(case.x, fusionName)) + geom_tile(aes(fill = DRIVER_PROB), width = 0.95, height = 0.95) + ylab("Predicted Fusion") + xlab(paste(project, "Case"))
dev.off()


##recFtx <- c("HMGA2--SMUG1", "CCDC169--XRCC2", "CS--", "HELB--", "RP11−96H19.1−−CCDC38", "HMGA2−−RP11−834C11.5")

##recOnc <- Onc[as.character(Onc$fusionName) %in% recFtx, ]
