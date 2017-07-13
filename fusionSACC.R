## working with salivary ACC

## libraries
library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

project <- "Salivary Gland ACC"
outDir <- "SACC"

## cases to include
saccInc <-  c("R2SG-1", "RSG-2", "RSG-3", "RSG-4", "RSG-5", "SG13T", "SG15T", "SG16T")

## extract SACC unfiltered
sacc.star <- star[star$case %in% saccInc,]
sacc.map <- map[map$case %in% saccInc,]
sacc.ee <- ee[ee$case %in% saccInc,]
sacc.integrate <- integrate[integrate$case %in% saccInc,]
sacc.defuse <- defuse[defuse$case %in% saccInc,]

## extract SACC from individual filtered fusions
sacc.star.filtered <- star.filtered[star.filtered$project == project, ]
sacc.star.filtered <- sacc.star.filtered[sacc.star.filtered$case %in% saccInc, ]
sacc.star.filtered <- sacc.star.filtered[!is.na(sacc.star.filtered$case),]

sacc.map.filtered <- map.filtered[map.filtered$project == project, ]
sacc.map.filtered <- sacc.map.filtered[sacc.map.filtered$case %in% saccInc, ]
sacc.map.filtered <- sacc.map.filtered[!is.na(sacc.map.filtered$case),]

sacc.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
sacc.ee.filtered <- sacc.ee.filtered[sacc.ee.filtered$case %in% saccInc, ]
sacc.ee.filtered <- sacc.ee.filtered[!is.na(sacc.ee.filtered$case),]

sacc.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
sacc.integrate.filtered <- sacc.integrate.filtered[sacc.integrate.filtered$case %in% saccInc, ]
sacc.integrate.filtered <- sacc.integrate.filtered[!is.na(sacc.integrate.filtered$case),]

sacc.defuse.filtered <- defuse.filtered[defuse.filtered$project == project, ]
sacc.defuse.filtered <- sacc.defuse.filtered[sacc.defuse.filtered$case %in% saccInc, ]
sacc.defuse.filtered <- sacc.defuse.filtered[!is.na(sacc.defuse.filtered$case),]

## extract SACC from recurrent filtered fusions
sacc.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
sacc.star.recurrent <- sacc.star.recurrent[sacc.star.recurrent$case %in% saccInc, ]
sacc.star.recurrent <- sacc.star.recurrent[!is.na(sacc.star.recurrent$case),]

sacc.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
sacc.map.recurrent <- sacc.map.recurrent[sacc.map.recurrent$case %in% saccInc, ]
sacc.map.recurrent <- sacc.map.recurrent[!is.na(sacc.map.recurrent$case),]

sacc.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
sacc.ee.recurrent <- sacc.ee.recurrent[sacc.ee.recurrent$case %in% saccInc, ]
sacc.ee.recurrent <- sacc.ee.recurrent[!is.na(sacc.ee.recurrent$case),]

sacc.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
sacc.ifr.recurrent <- sacc.ifr.recurrent[sacc.ifr.recurrent$case %in% saccInc, ]
sacc.ifr.recurrent <- sacc.ifr.recurrent[!is.na(sacc.ifr.recurrent$case),]

sacc.def.recurrent <- def.recurrent[def.recurrent$project == project, ]
sacc.def.recurrent <- sacc.def.recurrent[sacc.def.recurrent$case %in% saccInc, ]
sacc.def.recurrent <- sacc.def.recurrent[!is.na(sacc.def.recurrent$case),]


## pair-wise interesections
saccSM <- merge(sacc.star.filtered, sacc.map.filtered, by = c("fusionName", "case"))
saccSI <- merge(sacc.star.filtered, sacc.integrate.filtered, by = c("fusionName", "case"))
saccSE <- merge(sacc.star.filtered, sacc.ee.filtered, by = c("fusionName", "case"))
saccSD <- merge(sacc.star.filtered, sacc.defuse.filtered, by = c("fusionName", "case"))
saccME <- merge(sacc.map.filtered, sacc.ee.filtered, by = c("fusionName", "case"))
saccMI <- merge(sacc.map.filtered, sacc.integrate.filtered, by = c("fusionName", "case"))
saccMD <- merge(sacc.map.filtered, sacc.defuse.filtered, by = c("fusionName", "case"))
saccEI <- merge(sacc.ee.filtered, sacc.integrate.filtered, by = c("fusionName", "case"))
saccDI <- merge(sacc.defuse.filtered, sacc.integrate.filtered, by = c("fusionName", "case"))
saccDE <- merge(sacc.defuse.filtered, sacc.ee.filtered, by = c("fusionName", "case"))

## three-way interesections
saccMSE <- merge(saccSM, sacc.ee.filtered, by = c("fusionName", "case"))
saccMSI <- merge(saccSM, sacc.integrate.filtered, by = c("fusionName", "case"))
saccMSD <- merge(saccSM, sacc.defuse.filtered, by = c("fusionName", "case"))
saccMEI <- merge(saccME, sacc.integrate.filtered, by = c("fusionName", "case"))
saccMED <- merge(saccME, sacc.defuse.filtered, by = c("fusionName", "case"))
saccSEI <- merge(saccSE, sacc.integrate.filtered, by = c("fusionName", "case"))
saccSED <- merge(saccSE, sacc.defuse.filtered, by = c("fusionName", "case"))
saccMID <- merge(saccMI, sacc.defuse.filtered, by = c("fusionName", "case"))
saccIED <- merge(saccEI, sacc.defuse.filtered, by = c("fusionName", "case"))

## all programs
#sacc4M <- merge(saccSM, saccEI, by = c("fusionName", "case"))

## write individual intersections
##write.table(saccSM, file = "SACC/SACC_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(saccSE, file = "SACC/SACC_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(saccME, file = "SACC/SACC_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(saccMSE, file = "SACC/SACC_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(saccMSI, file = "SACC/SACC_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(saccMEI, file = "SACC/SACC_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(saccSEI, file = "SACC/SACC_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(sacc4M, file = "SACC/SACC_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(sacc.star.recurrent, file = "SACC/SACC_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.map.recurrent, file = "SACC/SACC_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.ee.recurrent, file = "SACC/SACC_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.ifr.recurrent, file = "SACC/SACC_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.def.recurrent, file = "SACC/SACC_recurrent_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write Unfiltered
write.table(sacc.star, file = "SACC/SACC_unfiltered_fusions_STAR-Fusion.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.map, file = "SACC/SACC_unfiltered_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.ee, file = "SACC/SACC_unfiltered_fusions_Ericscript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.int, file = "SACC/SACC_unfiltered_fusions_INTEGRATE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.def, file = "SACC/SACC_unfiltered_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

unfilt <- merge(sacc.star, sacc.map, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, sacc.ee, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, sacc.integrate, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, sacc.defuse, by = c("fusionName", "case"), all = TRUE)

write.table(unfilt, file = "SACC/SACC_unfiltered_fusions_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## count individual filtered by programs
saccSCount <- countFusions(sacc.star.filtered, "fusionName", "case")
saccMCount <- countFusions(sacc.map.filtered, "fusionName", "case")
saccECount <- countFusions(sacc.ee.filtered, "fusionName", "case")
saccICount <- countFusions(sacc.integrate.filtered, "fusionName", "case")
## count recurrent filtered by programs
saccSfrCount <- countFusions(sacc.star.recurrent, "fusionName", "case")
saccMfrCount <- countFusions(sacc.map.recurrent, "fusionName", "case")
saccEfrCount <- countFusions(sacc.ee.recurrent, "fusionName", "case")
saccIfrCount <- countFusions(sacc.ifr.recurrent, "fusionName", "case")

## count intersections
saccSMCounts <- countFusions(saccSM, "fusionName", "case")
saccSECounts <- countFusions(saccSE, "fusionName", "case")
saccSICounts <- countFusions(saccSI, "fusionName", "case")
saccSDCounts <- countFusions(saccSD, "fusionName", "case")
saccMECounts <- countFusions(saccME, "fusionName", "case")
saccMICounts <- countFusions(saccMI, "fusionName", "case")
saccMDCounts <- countFusions(saccMD, "fusionName", "case")
saccEICounts <- countFusions(saccEI, "fusionName", "case")
saccDICounts <- countFusions(saccDI, "fusionName", "case")
saccDECounts <- countFusions(saccDE, "fusionName", "case")


##saccCounts <- countFusions(sacc4M, "fusionName", "case")

fusionColors = c("Not Present" = "#CFCFCF", "Predicted" = "#F8766D")

pdf("SACC/SACC_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(saccSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(saccMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(saccECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(saccICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = fusionColors)
dev.off()



pdf("SACC/SACC_Recurrent_byProgram.pdf", width = 150/25.4, height = 180/25.4)
ggplot(saccSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(saccMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = fusionColors)
ggplot(saccEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = fusionColors)
ggplot(saccIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = fusionColors)
dev.off()


pdf("SACC/SACC_predictedFusions.pdf", width = 105/25.4, height = 105/25.4)
ggplot(saccSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccSICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccSDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & DeFuse Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccMDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & DeFuse Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccMICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccEICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("EricScript & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccDICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(saccDECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
dev.off()

##ggplot(saccCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
##    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
##    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
##    scale_fill_manual(values = fusionColors)
##                      


sacc.MIN2 <- unique(do.call(rbind, lapply(list(saccSM, saccSI, saccSE, saccSD, saccME, saccMI, saccMD, saccEI, saccDI, saccDE), function(X) X[,c("fusionName", "case")])))
sacc.MIN2 <- merge(sacc.MIN2, sacc.ee, all.x = TRUE, by = c("fusionName", "case"))
sacc.MIN2 <- merge(sacc.MIN2, sacc.map, all.x = TRUE, by = c("fusionName", "case"))
sacc.MIN2 <- merge(sacc.MIN2, sacc.star, all.x = TRUE, by = c("fusionName", "case"))
sacc.MIN2 <- merge(sacc.MIN2, sacc.integrate, all.x = TRUE, by = c("fusionName", "case"))
sacc.MIN2 <- merge(sacc.MIN2, sacc.defuse, all.x = TRUE, by = c("fusionName", "case"))
##sacc.MIN2$fusionName <- droplevels(sacc.MIN2$fusionName)

sacc.MIN3 <- unique(do.call(rbind, lapply(list(saccMSE, saccMSI, saccMSD, saccMEI, saccMED, saccSEI, saccSED, saccMID, saccIED), function(X) X[,c("fusionName", "case")])))
sacc.MIN3 <- merge(sacc.MIN3, sacc.ee, all.x = TRUE, by = c("fusionName", "case"))
sacc.MIN3 <- merge(sacc.MIN3, sacc.map, all.x = TRUE, by = c("fusionName", "case"))
sacc.MIN3 <- merge(sacc.MIN3, sacc.star, all.x = TRUE, by = c("fusionName", "case"))
sacc.MIN3 <- merge(sacc.MIN3, sacc.integrate, all.x = TRUE, by = c("fusionName", "case"))
sacc.MIN3 <- merge(sacc.MIN3, sacc.defuse, all.x = TRUE, by = c("fusionName", "case"))
## sacc.MIN3$fusionName <- droplevels(sacc.MIN3$fusionName)

min2Count <- countFusions(sacc.MIN2, caseID = "case")
min3Count <- countFusions(sacc.MIN3, caseID = "case")
#min4Count <- countFusions(sacc.MIN4, caseID = "case.x")


pdf("SACC/SACC_MIN2_MIN3.pdf", width = 110/25.4, height = 297/25.4)
ggplot(min2Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10), legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least two tools") +
    scale_fill_manual(values = fusionColors)
ggplot(min3Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least three tools") +
    scale_fill_manual(values = fusionColors) 
dev.off()


write.table(sacc.MIN2, file = "SACC/SACC_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sacc.MIN3, file = "SACC/SACC_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##write.table(sacc4M, file = "SACC/SACC_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("SACC/fusions.sacc.rda")

load("SACC/fusions.sacc.rda")
## java -jar /opt/src/oncofuse-1.1.1/Oncofuse.jar -a hg19 SACC/sacc_min2_oncofuse.input.txt 'coord' EPI SACC/sacc_min2_oncofuse.output.txt

Onc <- do.call(rbind, lapply(list.files(pattern = "merged.*oncofuse.output.txt"), read.delim, header = TRUE))
Onc$fusionName <- paste(Onc$X5_FPG_GENE_NAME, Onc$X3_FPG_GENE_NAME, "--")

saccOnc <- rbind(sacc.ee.recurrent[, c("program", "case", "fusionName")], sacc.map.recurrent[, c("program", "case", "fusionName")], sacc.star.recurrent[, c("program", "case", "fusionName")], sacc.ifr.recurrent[, c("program", "case", "fusionName")])

saccOnc <- merge(saccOnc, Onc, by = "fusionName")

pdf("SACC/SACC_Fusions_DrivProbability.pdf", width = 148/25.4, height = 210/25.4, useDingbats = FALSE)
ggplot(saccOnc, aes(case.x, fusionName)) + geom_tile(aes(fill = DRIVER_PROB), width = 0.95, height = 0.95) + ylab("Predicted Fusion") + xlab(paste(project, "Case"))
dev.off()


##recFtx <- c("HMGA2--SMUG1", "CCDC169--XRCC2", "CS--", "HELB--", "RP11−96H19.1−−CCDC38", "HMGA2−−RP11−834C11.5")

##recOnc <- Onc[as.character(Onc$fusionName) %in% recFtx, ]
