## working with BMEC
## libraries
library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")
## load("BMEC/fusions.bmec.rda")

project <- "Mucoepidermoid Carcinoma"
outDir <- "BMEC"

## cases to include
bmecInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project], defuse$case[defuse$project == project]))

## extract BMEC unfiltered
bmec.star <- star[star$case %in% bmecInc,]
bmec.map <- map[map$case %in% bmecInc,]
bmec.ee <- ee[ee$case %in% bmecInc,]
bmec.integrate <- integrate[integrate$case %in% bmecInc,]
bmec.defuse <- defuse[defuse$case %in% bmecInc,]

## extract BMEC from individual filtered fusions
bmec.star.filtered <- star.filtered[star.filtered$project == project, ]
bmec.star.filtered <- bmec.star.filtered[!is.na(bmec.star.filtered$case),]
bmec.star.filtered <- bmec.star.filtered[bmec.star.filtered$case %in% bmecInc, ]

bmec.map.filtered <- map.filtered[map.filtered$project == project, ]
bmec.map.filtered <- bmec.map.filtered[!is.na(bmec.map.filtered$case),]
bmec.map.filtered <- bmec.map.filtered[bmec.map.filtered$case %in% bmecInc, ]

bmec.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
bmec.ee.filtered <- bmec.ee.filtered[!is.na(bmec.ee.filtered$case),]
bmec.ee.filtered <- bmec.ee.filtered[bmec.ee.filtered$case %in% bmecInc, ]

bmec.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
bmec.integrate.filtered <- bmec.integrate.filtered[!is.na(bmec.integrate.filtered$case),]
bmec.integrate.filtered <- bmec.integrate.filtered[bmec.integrate.filtered$case %in% bmecInc, ]

bmec.defuse.filtered <- defuse.filtered[defuse.filtered$project == project, ]
bmec.defuse.filtered <- bmec.defuse.filtered[!is.na(bmec.defuse.filtered$case),]
bmec.defuse.filtered <- bmec.defuse.filtered[bmec.defuse.filtered$case %in% bmecInc, ]

## extract BMEC from recurrent filtered fusions
bmec.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
bmec.star.recurrent <- bmec.star.recurrent[!is.na(bmec.star.recurrent$case),]
bmec.star.recurrent <- bmec.star.recurrent[bmec.star.recurrent$case %in% bmecInc, ]

bmec.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
bmec.map.recurrent <- bmec.map.recurrent[!is.na(bmec.map.recurrent$case),]
bmec.map.recurrent <- bmec.map.recurrent[bmec.map.recurrent$case %in% bmecInc, ]

bmec.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
bmec.ee.recurrent <- bmec.ee.recurrent[!is.na(bmec.ee.recurrent$case),]
bmec.ee.recurrent <- bmec.ee.recurrent[bmec.ee.recurrent$case %in% bmecInc, ]

bmec.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
bmec.ifr.recurrent <- bmec.ifr.recurrent[!is.na(bmec.ifr.recurrent$case),]
bmec.ifr.recurrent <- bmec.ifr.recurrent[bmec.ifr.recurrent$case %in% bmecInc, ]

bmec.def.recurrent <- def.recurrent[def.recurrent$project == project, ]
bmec.def.recurrent <- bmec.def.recurrent[!is.na(bmec.def.recurrent$case),]
bmec.def.recurrent <- bmec.def.recurrent[bmec.def.recurrent$case %in% bmecInc, ]


## pair-wise interesections
bmecSM <- merge(bmec.star, bmec.map, by = c("fusionName", "case"))
bmecSI <- merge(bmec.star, bmec.integrate, by = c("fusionName", "case"))
bmecSE <- merge(bmec.star, bmec.ee, by = c("fusionName", "case"))
bmecSD <- merge(bmec.star, bmec.defuse, by = c("fusionName", "case"))
bmecME <- merge(bmec.map, bmec.ee, by = c("fusionName", "case"))
bmecMI <- merge(bmec.map, bmec.integrate, by = c("fusionName", "case"))
bmecMD <- merge(bmec.map, bmec.defuse, by = c("fusionName", "case"))
bmecEI <- merge(bmec.ee, bmec.integrate, by = c("fusionName", "case"))
bmecDI <- merge(bmec.defuse, bmec.integrate, by = c("fusionName", "case"))
bmecDE <- merge(bmec.defuse, bmec.ee, by = c("fusionName", "case"))

## three-way interesections
bmecMSE <- merge(bmecSM, bmec.ee, by = c("fusionName", "case"))
bmecMSI <- merge(bmecSM, bmec.integrate, by = c("fusionName", "case"))
bmecMSD <- merge(bmecSM, bmec.defuse, by = c("fusionName", "case"))
bmecMEI <- merge(bmecME, bmec.integrate, by = c("fusionName", "case"))
bmecMED <- merge(bmecME, bmec.defuse, by = c("fusionName", "case"))
bmecSEI <- merge(bmecSE, bmec.integrate, by = c("fusionName", "case"))
bmecSED <- merge(bmecSE, bmec.defuse, by = c("fusionName", "case"))
bmecMID <- merge(bmecMI, bmec.defuse, by = c("fusionName", "case"))
bmecIED <- merge(bmecEI, bmec.defuse, by = c("fusionName", "case"))

## all programs
#bmec4M <- merge(bmecSM, bmecEI, by = c("fusionName", "case"))

## write individual intersections
##write.table(bmecSM, file = "BMEC/BMEC_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(bmecSE, file = "BMEC/BMEC_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(bmecME, file = "BMEC/BMEC_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bmecMSE, file = "BMEC/BMEC_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bmecMSI, file = "BMEC/BMEC_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bmecMEI, file = "BMEC/BMEC_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bmecSEI, file = "BMEC/BMEC_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(bmec4M, file = "BMEC/BMEC_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(bmec.star.recurrent, file = "BMEC/BMEC_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.map.recurrent, file = "BMEC/BMEC_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.ee.recurrent, file = "BMEC/BMEC_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.ifr.recurrent, file = "BMEC/BMEC_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.def.recurrent, file = "BMEC/BMEC_recurrent_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write Unfiltered
write.table(bmec.star, file = "BMEC/BMEC_unfiltered_fusions_STAR-Fusion.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.map, file = "BMEC/BMEC_unfiltered_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.ee, file = "BMEC/BMEC_unfiltered_fusions_Ericscript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.int, file = "BMEC/BMEC_unfiltered_fusions_INTEGRATE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.def, file = "BMEC/BMEC_unfiltered_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

unfilt <- merge(bmec.star, bmec.map, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, bmec.ee, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, bmec.int, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, bmec.def, by = c("fusionName", "case"), all = TRUE)

write.table(unfilt, file = "BMEC/BMEC_unfiltered_fusions_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## count individual filtered by programs
bmecSCount <- countFusions(bmec.star.filtered, "fusionName", "case")
bmecMCount <- countFusions(bmec.map.filtered, "fusionName", "case")
bmecECount <- countFusions(bmec.ee.filtered, "fusionName", "case")
bmecICount <- countFusions(bmec.integrate.filtered, "fusionName", "case")
bmecDCount <- countFusions(bmec.defuse.filtered, "fusionName", "case")
## count recurrent filtered by programs
bmecSfrCount <- countFusions(bmec.star.recurrent, "fusionName", "case")
bmecMfrCount <- countFusions(bmec.map.recurrent, "fusionName", "case")
bmecEfrCount <- countFusions(bmec.ee.recurrent, "fusionName", "case")
bmecIfrCount <- countFusions(bmec.ifr.recurrent, "fusionName", "case")
bmecDfrCount <- countFusions(bmec.def.recurrent, "fusionName", "case")

## count intersections
bmecSMCounts <- countFusions(bmecSM, "fusionName", "case")
bmecSECounts <- countFusions(bmecSE, "fusionName", "case")
bmecSICounts <- countFusions(bmecSI, "fusionName", "case")
bmecSDCounts <- countFusions(bmecSD, "fusionName", "case")
bmecMECounts <- countFusions(bmecME, "fusionName", "case")
bmecMICounts <- countFusions(bmecMI, "fusionName", "case")
bmecMDCounts <- countFusions(bmecMD, "fusionName", "case")
bmecEICounts <- countFusions(bmecEI, "fusionName", "case")
bmecDICounts <- countFusions(bmecDI, "fusionName", "case")
bmecDECounts <- countFusions(bmecDE, "fusionName", "case")


##bmecCounts <- countFusions(bmec4M, "fusionName", "case")

fusionColors = c("Not Present" = "#CFCFCF", "Predicted" = "#F8766D")

pdf("BMEC/BMEC_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(bmecSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(bmecMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(bmecECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(bmecICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(bmecDCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("DeFuse predicted") +
    scale_fill_manual(values = fusionColors)
dev.off()



pdf("BMEC/BMEC_Recurrent_byProgram.pdf", width = 150/25.4, height = 180/25.4)
ggplot(bmecSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion Recurrent") +
    scale_fill_manual(values = fusionColors)
ggplot(bmecMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(bmecEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(bmecIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(bmecDfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Defuse Recurrent") + 
    scale_fill_manual(values = fusionColors)
dev.off()


pdf("BMEC/BMEC_PairWiseFusions.pdf", width = 105/25.4, height = 105/25.4)
ggplot(bmecSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecSICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecSDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & DeFuse Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecMDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & DeFuse Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecMICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecEICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("EricScript & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecDICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(bmecDECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & EricScript Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
dev.off()

##ggplot(bmecCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
##    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
##    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
##    scale_fill_manual(values = fusionColors)
##                      


bmec.MIN2 <- unique(do.call(rbind, lapply(list(bmecSM, bmecSI, bmecSE, bmecSD, bmecME, bmecMI, bmecMD, bmecEI, bmecDI, bmecDE), function(X) X[,c("fusionName", "case")])))
bmec.MIN2 <- merge(bmec.MIN2, bmec.ee, all.x = TRUE, by = c("fusionName"))
bmec.MIN2 <- merge(bmec.MIN2, bmec.map, all.x = TRUE, by = c("fusionName"))
bmec.MIN2 <- merge(bmec.MIN2, bmec.star, all.x = TRUE, by = c("fusionName"))
bmec.MIN2 <- merge(bmec.MIN2, bmec.integrate, all.x = TRUE, by = c("fusionName"))
bmec.MIN2 <- merge(bmec.MIN2, bmec.defuse, all.x = TRUE, by = c("fusionName"))
##bmec.MIN2$fusionName <- droplevels(bmec.MIN2$fusionName)

bmec.MIN3 <- unique(do.call(rbind, lapply(list(bmecMSE, bmecMSI, bmecMSD, bmecMEI, bmecMED, bmecSEI, bmecSED, bmecMID, bmecIED), function(X) X[,c("fusionName", "case")])))
bmec.MIN3 <- merge(bmec.MIN3, bmec.ee, all.x = TRUE, by = c("fusionName"))
bmec.MIN3 <- merge(bmec.MIN3, bmec.map, all.x = TRUE, by = c("fusionName"))
bmec.MIN3 <- merge(bmec.MIN3, bmec.star, all.x = TRUE, by = c("fusionName"))
bmec.MIN3 <- merge(bmec.MIN3, bmec.integrate, all.x = TRUE, by = c("fusionName"))
bmec.MIN3 <- merge(bmec.MIN3, bmec.defuse, all.x = TRUE, by = c("fusionName"))
## bmec.MIN3$fusionName <- droplevels(bmec.MIN3$fusionName)

min2Count <- countFusions(bmec.MIN2, caseID = "case.x")
min3Count <- countFusions(bmec.MIN3, caseID = "case.x")
#min4Count <- countFusions(bmec.MIN4, caseID = "case.x")


pdf("BMEC/BMEC_MIN2_MIN3.pdf", width = 110/25.4, height = 297/25.4)
ggplot(min2Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10), legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least two tools") +
    scale_fill_manual(values = fusionColors)
ggplot(min3Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least three tools") +
    scale_fill_manual(values = fusionColors) 
dev.off()


write.table(bmec.MIN2, file = "BMEC/BMEC_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bmec.MIN3, file = "BMEC/BMEC_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##write.table(bmec4M, file = "BMEC/BMEC_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("BMEC/fusions.bmec.rda")

load("BMEC/fusions.bmec.rda")
## java -jar /opt/src/oncofuse-1.1.1/Oncofuse.jar -a hg19 BMEC/bmec_min2_oncofuse.input.txt 'coord' EPI BMEC/bmec_min2_oncofuse.output.txt

Onc <- do.call(rbind, lapply(list.files(pattern = "merged.*oncofuse.output.txt"), read.delim, header = TRUE))
Onc$fusionName <- paste(Onc$X5_FPG_GENE_NAME, Onc$X3_FPG_GENE_NAME, "--")

bmecOnc <- rbind(bmec.ee.recurrent[, c("program", "case", "fusionName")], bmec.map.recurrent[, c("program", "case", "fusionName")], bmec.star.recurrent[, c("program", "case", "fusionName")], bmec.ifr.recurrent[, c("program", "case", "fusionName")])

bmecOnc <- merge(bmecOnc, Onc, by = "fusionName")

pdf("BMEC/BMEC_Fusions_DrivProbability.pdf", width = 148/25.4, height = 210/25.4, useDingbats = FALSE)
ggplot(bmecOnc, aes(case.x, fusionName)) + geom_tile(aes(fill = DRIVER_PROB), width = 0.95, height = 0.95) + ylab("Predicted Fusion") + xlab(paste(project, "Case"))
dev.off()


##recFtx <- c("HMGA2--SMUG1", "CCDC169--XRCC2", "CS--", "HELB--", "RP11−96H19.1−−CCDC38", "HMGA2−−RP11−834C11.5")

##recOnc <- Onc[as.character(Onc$fusionName) %in% recFtx, ]
