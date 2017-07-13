## working with BPA
library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

project <- "Pleomorphic Adenoma"
outDir <- "BPA"

## cases to include
bpaInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project], defuse$case[defuse$project == project]))

## extract BPA unfiltered
bpa.star <- star[star$case %in% bpaInc,]
bpa.map <- map[map$case %in% bpaInc,]
bpa.ee <- ee[ee$case %in% bpaInc,]
bpa.integrate <- integrate[integrate$case %in% bpaInc,]
bpa.defuse <- defuse[defuse$case %in% bpaInc,]

## extract BPA from individual filtered fusions
bpa.star.filtered <- star.filtered[star.filtered$project == project, ]
bpa.star.filtered <- bpa.star.filtered[!is.na(bpa.star.filtered$case),]
bpa.star.filtered <- bpa.star.filtered[bpa.star.filtered$case %in% bpaInc, ]

bpa.map.filtered <- map.filtered[map.filtered$project == project, ]
bpa.map.filtered <- bpa.map.filtered[!is.na(bpa.map.filtered$case),]
bpa.map.filtered <- bpa.map.filtered[bpa.map.filtered$case %in% bpaInc, ]

bpa.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
bpa.ee.filtered <- bpa.ee.filtered[!is.na(bpa.ee.filtered$case),]
bpa.ee.filtered <- bpa.ee.filtered[bpa.ee.filtered$case %in% bpaInc, ]

bpa.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
bpa.integrate.filtered <- bpa.integrate.filtered[!is.na(bpa.integrate.filtered$case),]
bpa.integrate.filtered <- bpa.integrate.filtered[bpa.integrate.filtered$case %in% bpaInc, ]

bpa.defuse.filtered <- defuse.filtered[defuse.filtered$project == project, ]
bpa.defuse.filtered <- bpa.defuse.filtered[!is.na(bpa.defuse.filtered$case),]
bpa.defuse.filtered <- bpa.defuse.filtered[bpa.defuse.filtered$case %in% bpaInc, ]

## extract BPA from recurrent filtered fusions
bpa.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
bpa.star.recurrent <- bpa.star.recurrent[!is.na(bpa.star.recurrent$case),]
bpa.star.recurrent <- bpa.star.recurrent[bpa.star.recurrent$case %in% bpaInc, ]

bpa.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
bpa.map.recurrent <- bpa.map.recurrent[!is.na(bpa.map.recurrent$case),]
bpa.map.recurrent <- bpa.map.recurrent[bpa.map.recurrent$case %in% bpaInc, ]

bpa.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
bpa.ee.recurrent <- bpa.ee.recurrent[!is.na(bpa.ee.recurrent$case),]
bpa.ee.recurrent <- bpa.ee.recurrent[bpa.ee.recurrent$case %in% bpaInc, ]

bpa.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
bpa.ifr.recurrent <- bpa.ifr.recurrent[!is.na(bpa.ifr.recurrent$case),]
bpa.ifr.recurrent <- bpa.ifr.recurrent[bpa.ifr.recurrent$case %in% bpaInc, ]

bpa.def.recurrent <- def.recurrent[def.recurrent$project == project, ]
bpa.def.recurrent <- bpa.def.recurrent[!is.na(bpa.def.recurrent$case),]
bpa.def.recurrent <- bpa.def.recurrent[bpa.def.recurrent$case %in% bpaInc, ]


## pair-wise interesections
bpaSM <- merge(bpa.star, bpa.map, by = c("fusionName", "case"))
bpaSI <- merge(bpa.star, bpa.integrate, by = c("fusionName", "case"))
bpaSE <- merge(bpa.star, bpa.ee, by = c("fusionName", "case"))
bpaSD <- merge(bpa.star, bpa.defuse, by = c("fusionName", "case"))
bpaME <- merge(bpa.map, bpa.ee, by = c("fusionName", "case"))
bpaMI <- merge(bpa.map, bpa.integrate, by = c("fusionName", "case"))
bpaMD <- merge(bpa.map, bpa.defuse, by = c("fusionName", "case"))
bpaEI <- merge(bpa.ee, bpa.integrate, by = c("fusionName", "case"))
bpaDI <- merge(bpa.defuse, bpa.integrate, by = c("fusionName", "case"))
bpaDE <- merge(bpa.defuse, bpa.ee, by = c("fusionName", "case"))

## three-way interesections
bpaMSE <- merge(bpaSM, bpa.ee, by = c("fusionName", "case"))
bpaMSI <- merge(bpaSM, bpa.integrate, by = c("fusionName", "case"))
bpaMSD <- merge(bpaSM, bpa.defuse, by = c("fusionName", "case"))
bpaMEI <- merge(bpaME, bpa.integrate, by = c("fusionName", "case"))
bpaMED <- merge(bpaME, bpa.defuse, by = c("fusionName", "case"))
bpaSEI <- merge(bpaSE, bpa.integrate, by = c("fusionName", "case"))
bpaSED <- merge(bpaSE, bpa.defuse, by = c("fusionName", "case"))
bpaMID <- merge(bpaMI, bpa.defuse, by = c("fusionName", "case"))
bpaIED <- merge(bpaEI, bpa.defuse, by = c("fusionName", "case"))

## all programs
#bpa4M <- merge(bpaSM, bpaEI, by = c("fusionName", "case"))

## write individual intersections
##write.table(bpaSM, file = "BPA/BPA_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(bpaSE, file = "BPA/BPA_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(bpaME, file = "BPA/BPA_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bpaMSE, file = "BPA/BPA_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bpaMSI, file = "BPA/BPA_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bpaMEI, file = "BPA/BPA_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bpaSEI, file = "BPA/BPA_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bpa4M, file = "BPA/BPA_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(bpa.star.recurrent, file = "BPA/BPA_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.map.recurrent, file = "BPA/BPA_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.ee.recurrent, file = "BPA/BPA_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.ifr.recurrent, file = "BPA/BPA_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.def.recurrent, file = "BPA/BPA_recurrent_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write Unfiltered
write.table(bpa.star, file = "BPA/BPA_unfiltered_fusions_STAR-Fusion.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.map, file = "BPA/BPA_unfiltered_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.ee, file = "BPA/BPA_unfiltered_fusions_Ericscript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.int, file = "BPA/BPA_unfiltered_fusions_INTEGRATE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.def, file = "BPA/BPA_unfiltered_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

unfilt <- merge(bpa.star, bpa.map, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, bpa.ee, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, bpa.int, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, bpa.def, by = c("fusionName", "case"), all = TRUE)

write.table(unfilt, file = "BPA/BPA_unfiltered_fusions_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## count individual filtered by programs
bpaSCount <- countFusions(bpa.star.filtered, "fusionName", "case")
bpaMCount <- countFusions(bpa.map.filtered, "fusionName", "case")
bpaECount <- countFusions(bpa.ee.filtered, "fusionName", "case")
bpaICount <- countFusions(bpa.integrate.filtered, "fusionName", "case")
bpaDCount <- countFusions(bpa.defuse.filtered, "fusionName", "case")
## count recurrent filtered by programs
bpaSfrCount <- countFusions(bpa.star.recurrent, "fusionName", "case")
bpaMfrCount <- countFusions(bpa.map.recurrent, "fusionName", "case")
bpaEfrCount <- countFusions(bpa.ee.recurrent, "fusionName", "case")
bpaIfrCount <- countFusions(bpa.ifr.recurrent, "fusionName", "case")
bpaDfrCount <- countFusions(bpa.def.recurrent, "fusionName", "case")

## count intersections
bpaSMCounts <- countFusions(bpaSM, "fusionName", "case")
bpaSECounts <- countFusions(bpaSE, "fusionName", "case")
bpaSICounts <- countFusions(bpaSI, "fusionName", "case")
bpaSDCounts <- countFusions(bpaSD, "fusionName", "case")
bpaMECounts <- countFusions(bpaME, "fusionName", "case")
bpaMICounts <- countFusions(bpaMI, "fusionName", "case")
bpaMDCounts <- countFusions(bpaMD, "fusionName", "case")
bpaEICounts <- countFusions(bpaEI, "fusionName", "case")
bpaDICounts <- countFusions(bpaDI, "fusionName", "case")
bpaDECounts <- countFusions(bpaDE, "fusionName", "case")


##bpaCounts <- countFusions(bpa4M, "fusionName", "case")

fusionColors = c("Not Present" = "#CFCFCF", "Predicted" = "#F8766D")

pdf("BPA/BPA_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(bpaSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(bpaMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(bpaECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(bpaICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(bpaDCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("DeFuse predicted") +
    scale_fill_manual(values = fusionColors)
dev.off()



pdf("BPA/BPA_Recurrent_byProgram.pdf", width = 150/25.4, height = 180/25.4)
ggplot(bpaSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion Recurrent") +
    scale_fill_manual(values = fusionColors)
ggplot(bpaMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(bpaEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(bpaIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(bpaDfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Defuse Recurrent") + 
    scale_fill_manual(values = fusionColors)
dev.off()


pdf("BPA/BPA_PairWiseFusions.pdf", width = 105/25.4, height = 105/25.4)
ggplot(bpaSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaSICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaSDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & DeFuse Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaMDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & DeFuse Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaMICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaEICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("EricScript & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaDICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bpaDECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & EricScript Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
dev.off()

##ggplot(bpaCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
##    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
##    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
##    scale_fill_manual(values = fusionColors)
##                      


bpa.MIN2 <- unique(do.call(rbind, lapply(list(bpaSM, bpaSI, bpaSE, bpaSD, bpaME, bpaMI, bpaMD, bpaEI, bpaDI, bpaDE), function(X) X[,c("fusionName", "case")])))
bpa.MIN2 <- merge(bpa.MIN2, bpa.ee, all.x = TRUE, by = c("fusionName"))
bpa.MIN2 <- merge(bpa.MIN2, bpa.map, all.x = TRUE, by = c("fusionName"))
bpa.MIN2 <- merge(bpa.MIN2, bpa.star, all.x = TRUE, by = c("fusionName"))
bpa.MIN2 <- merge(bpa.MIN2, bpa.integrate, all.x = TRUE, by = c("fusionName"))
bpa.MIN2 <- merge(bpa.MIN2, bpa.defuse, all.x = TRUE, by = c("fusionName"))
##bpa.MIN2$fusionName <- droplevels(bpa.MIN2$fusionName)

bpa.MIN3 <- unique(do.call(rbind, lapply(list(bpaMSE, bpaMSI, bpaMSD, bpaMEI, bpaMED, bpaSEI, bpaSED, bpaMID, bpaIED), function(X) X[,c("fusionName", "case")])))
bpa.MIN3 <- merge(bpa.MIN3, bpa.ee, all.x = TRUE, by = c("fusionName"))
bpa.MIN3 <- merge(bpa.MIN3, bpa.map, all.x = TRUE, by = c("fusionName"))
bpa.MIN3 <- merge(bpa.MIN3, bpa.star, all.x = TRUE, by = c("fusionName"))
bpa.MIN3 <- merge(bpa.MIN3, bpa.integrate, all.x = TRUE, by = c("fusionName"))
bpa.MIN3 <- merge(bpa.MIN3, bpa.defuse, all.x = TRUE, by = c("fusionName"))
## bpa.MIN3$fusionName <- droplevels(bpa.MIN3$fusionName)

min2Count <- countFusions(bpa.MIN2, caseID = "case.x")
min3Count <- countFusions(bpa.MIN3, caseID = "case.x")
#min4Count <- countFusions(bpa.MIN4, caseID = "case.x")


pdf("BPA/BPA_MIN2_MIN3.pdf", width = 110/25.4, height = 297/25.4)
ggplot(min2Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10), legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least two tools") +
    scale_fill_manual(values = fusionColors)
ggplot(min3Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least three tools") +
    scale_fill_manual(values = fusionColors) 
dev.off()


write.table(bpa.MIN2, file = "BPA/BPA_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bpa.MIN3, file = "BPA/BPA_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##write.table(bpa4M, file = "BPA/BPA_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("BPA/fusions.bpa.rda")

load("BPA/fusions.bpa.rda")
## java -jar /opt/src/oncofuse-1.1.1/Oncofuse.jar -a hg19 BPA/bpa_min2_oncofuse.input.txt 'coord' EPI BPA/bpa_min2_oncofuse.output.txt

Onc <- do.call(rbind, lapply(list.files(pattern = "merged.*oncofuse.output.txt"), read.delim, header = TRUE))
Onc$fusionName <- paste(Onc$X5_FPG_GENE_NAME, Onc$X3_FPG_GENE_NAME, "--")

bpaOnc <- rbind(bpa.ee.recurrent[, c("program", "case", "fusionName")], bpa.map.recurrent[, c("program", "case", "fusionName")], bpa.star.recurrent[, c("program", "case", "fusionName")], bpa.ifr.recurrent[, c("program", "case", "fusionName")])

bpaOnc <- merge(bpaOnc, Onc, by = "fusionName")

pdf("BPA/BPA_Fusions_DrivProbability.pdf", width = 148/25.4, height = 210/25.4, useDingbats = FALSE)
ggplot(bpaOnc, aes(case.x, fusionName)) + geom_tile(aes(fill = DRIVER_PROB), width = 0.95, height = 0.95) + ylab("Predicted Fusion") + xlab(paste(project, "Case"))
dev.off()


##recFtx <- c("HMGA2--SMUG1", "CCDC169--XRCC2", "CS--", "HELB--", "RP11−96H19.1−−CCDC38", "HMGA2−−RP11−834C11.5")

##recOnc <- Onc[as.character(Onc$fusionName) %in% recFtx, ]
