## working with salivary ACC

## libraries
library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")
## load("LPS/fusions.lps.rda")

project <- "Liposarcoma"
outDir <- "LPS"

## cases to include
lpsInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project], defuse$case[defuse$project == project]))

## extract LPS unfiltered
lps.star <- star[star$case %in% lpsInc,]
lps.map <- map[map$case %in% lpsInc,]
lps.ee <- ee[ee$case %in% lpsInc,]
lps.integrate <- integrate[integrate$case %in% lpsInc,]
lps.defuse <- defuse[defuse$case %in% lpsInc,]

## extract LPS from individual filtered fusions
lps.star.filtered <- star.filtered[star.filtered$project == project, ]
lps.star.filtered <- lps.star.filtered[!is.na(lps.star.filtered$case),]
lps.star.filtered <- lps.star.filtered[lps.star.filtered$case %in% lpsInc, ]

lps.map.filtered <- map.filtered[map.filtered$project == project, ]
lps.map.filtered <- lps.map.filtered[!is.na(lps.map.filtered$case),]
lps.map.filtered <- lps.map.filtered[lps.map.filtered$case %in% lpsInc, ]

lps.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
lps.ee.filtered <- lps.ee.filtered[!is.na(lps.ee.filtered$case),]
lps.ee.filtered <- lps.ee.filtered[lps.ee.filtered$case %in% lpsInc, ]

lps.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
lps.integrate.filtered <- lps.integrate.filtered[!is.na(lps.integrate.filtered$case),]
lps.integrate.filtered <- lps.integrate.filtered[lps.integrate.filtered$case %in% lpsInc, ]

lps.defuse.filtered <- defuse.filtered[defuse.filtered$project == project, ]
lps.defuse.filtered <- lps.defuse.filtered[!is.na(lps.defuse.filtered$case),]
lps.defuse.filtered <- lps.defuse.filtered[lps.defuse.filtered$case %in% lpsInc, ]

## extract LPS from recurrent filtered fusions
lps.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
lps.star.recurrent <- lps.star.recurrent[!is.na(lps.star.recurrent$case),]
lps.star.recurrent <- lps.star.recurrent[lps.star.recurrent$case %in% lpsInc, ]

lps.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
lps.map.recurrent <- lps.map.recurrent[!is.na(lps.map.recurrent$case),]
lps.map.recurrent <- lps.map.recurrent[lps.map.recurrent$case %in% lpsInc, ]

lps.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
lps.ee.recurrent <- lps.ee.recurrent[!is.na(lps.ee.recurrent$case),]
lps.ee.recurrent <- lps.ee.recurrent[lps.ee.recurrent$case %in% lpsInc, ]

lps.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
lps.ifr.recurrent <- lps.ifr.recurrent[!is.na(lps.ifr.recurrent$case),]
lps.ifr.recurrent <- lps.ifr.recurrent[lps.ifr.recurrent$case %in% lpsInc, ]

lps.def.recurrent <- def.recurrent[def.recurrent$project == project, ]
lps.def.recurrent <- lps.def.recurrent[!is.na(lps.def.recurrent$case),]
lps.def.recurrent <- lps.def.recurrent[lps.def.recurrent$case %in% lpsInc, ]


## pair-wise interesections
lpsSM <- merge(lps.star, lps.map, by = c("fusionName", "case"))
lpsSI <- merge(lps.star, lps.integrate, by = c("fusionName", "case"))
lpsSE <- merge(lps.star, lps.ee, by = c("fusionName", "case"))
lpsSD <- merge(lps.star, lps.defuse, by = c("fusionName", "case"))
lpsME <- merge(lps.map, lps.ee, by = c("fusionName", "case"))
lpsMI <- merge(lps.map, lps.integrate, by = c("fusionName", "case"))
lpsMD <- merge(lps.map, lps.defuse, by = c("fusionName", "case"))
lpsEI <- merge(lps.ee, lps.integrate, by = c("fusionName", "case"))
lpsDI <- merge(lps.defuse, lps.integrate, by = c("fusionName", "case"))
lpsDE <- merge(lps.defuse, lps.ee, by = c("fusionName", "case"))

## three-way interesections
lpsMSE <- merge(lpsSM, lps.ee, by = c("fusionName", "case"))
lpsMSI <- merge(lpsSM, lps.integrate, by = c("fusionName", "case"))
lpsMSD <- merge(lpsSM, lps.defuse, by = c("fusionName", "case"))
lpsMEI <- merge(lpsME, lps.integrate, by = c("fusionName", "case"))
lpsMED <- merge(lpsME, lps.defuse, by = c("fusionName", "case"))
lpsSEI <- merge(lpsSE, lps.integrate, by = c("fusionName", "case"))
lpsSED <- merge(lpsSE, lps.defuse, by = c("fusionName", "case"))
lpsMID <- merge(lpsMI, lps.defuse, by = c("fusionName", "case"))
lpsIED <- merge(lpsEI, lps.defuse, by = c("fusionName", "case"))

## all programs
#lps4M <- merge(lpsSM, lpsEI, by = c("fusionName", "case"))

## write individual intersections
##write.table(lpsSM, file = "LPS/LPS_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(lpsSE, file = "LPS/LPS_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(lpsME, file = "LPS/LPS_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(lpsMSE, file = "LPS/LPS_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(lpsMSI, file = "LPS/LPS_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(lpsMEI, file = "LPS/LPS_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(lpsSEI, file = "LPS/LPS_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(lps4M, file = "LPS/LPS_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(lps.star.recurrent, file = "LPS/LPS_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.map.recurrent, file = "LPS/LPS_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.ee.recurrent, file = "LPS/LPS_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.ifr.recurrent, file = "LPS/LPS_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.def.recurrent, file = "LPS/LPS_recurrent_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write Unfiltered
write.table(lps.star, file = "LPS/LPS_unfiltered_fusions_STAR-Fusion.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.map, file = "LPS/LPS_unfiltered_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.ee, file = "LPS/LPS_unfiltered_fusions_Ericscript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.int, file = "LPS/LPS_unfiltered_fusions_INTEGRATE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.def, file = "LPS/LPS_unfiltered_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

unfilt <- merge(lps.star, lps.map, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, lps.ee, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, lps.int, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, lps.def, by = c("fusionName", "case"), all = TRUE)

write.table(unfilt, file = "LPS/LPS_unfiltered_fusions_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## count individual filtered by programs
lpsSCount <- countFusions(lps.star.filtered, "fusionName", "case")
lpsMCount <- countFusions(lps.map.filtered, "fusionName", "case")
lpsECount <- countFusions(lps.ee.filtered, "fusionName", "case")
lpsICount <- countFusions(lps.integrate.filtered, "fusionName", "case")
lpsDCount <- countFusions(lps.defuse.filtered, "fusionName", "case")
## count recurrent filtered by programs
lpsSfrCount <- countFusions(lps.star.recurrent, "fusionName", "case")
lpsMfrCount <- countFusions(lps.map.recurrent, "fusionName", "case")
lpsEfrCount <- countFusions(lps.ee.recurrent, "fusionName", "case")
lpsIfrCount <- countFusions(lps.ifr.recurrent, "fusionName", "case")
lpsDfrCount <- countFusions(lps.def.recurrent, "fusionName", "case")

## count intersections
lpsSMCounts <- countFusions(lpsSM, "fusionName", "case")
lpsSECounts <- countFusions(lpsSE, "fusionName", "case")
lpsSICounts <- countFusions(lpsSI, "fusionName", "case")
lpsSDCounts <- countFusions(lpsSD, "fusionName", "case")
lpsMECounts <- countFusions(lpsME, "fusionName", "case")
lpsMICounts <- countFusions(lpsMI, "fusionName", "case")
lpsMDCounts <- countFusions(lpsMD, "fusionName", "case")
lpsEICounts <- countFusions(lpsEI, "fusionName", "case")
lpsDICounts <- countFusions(lpsDI, "fusionName", "case")
lpsDECounts <- countFusions(lpsDE, "fusionName", "case")


##lpsCounts <- countFusions(lps4M, "fusionName", "case")

fusionColors = c("Not Present" = "#CFCFCF", "Predicted" = "#F8766D")

pdf("LPS/LPS_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(lpsSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(lpsMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(lpsECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(lpsICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = fusionColors)
ggplot(lpsDCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = fusionColors)
dev.off()



pdf("LPS/LPS_Recurrent_byProgram.pdf", width = 150/25.4, height = 180/25.4)
ggplot(lpsSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion Recurrent") +
    scale_fill_manual(values = fusionColors)
ggplot(lpsMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(lpsEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(lpsIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate Recurrent") + 
    scale_fill_manual(values = fusionColors)
ggplot(lpsDfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Defuse Recurrent") + 
    scale_fill_manual(values = fusionColors)
dev.off()


pdf("LPS/LPS_predictedFusions.pdf", width = 105/25.4, height = 105/25.4)
ggplot(lpsSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsSICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsSDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & DeFuse Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsMDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & DeFuse Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsMICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsEICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("EricScript & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsDICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & Integrate Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
ggplot(lpsDECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & EricScript Joint fusion Prediction") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = fusionColors)
dev.off()

##ggplot(lpsCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
##    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
##    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
##    scale_fill_manual(values = fusionColors)
##                      


lps.MIN2 <- unique(do.call(rbind, lapply(list(lpsSM, lpsSI, lpsSE, lpsSD, lpsME, lpsMI, lpsMD, lpsEI, lpsDI, lpsDE), function(X) X[,c("fusionName", "case")])))
lps.MIN2 <- merge(lps.MIN2, lps.ee, all.x = TRUE, by = c("fusionName"))
lps.MIN2 <- merge(lps.MIN2, lps.map, all.x = TRUE, by = c("fusionName"))
lps.MIN2 <- merge(lps.MIN2, lps.star, all.x = TRUE, by = c("fusionName"))
lps.MIN2 <- merge(lps.MIN2, lps.integrate, all.x = TRUE, by = c("fusionName"))
lps.MIN2 <- merge(lps.MIN2, lps.defuse, all.x = TRUE, by = c("fusionName"))
##lps.MIN2$fusionName <- droplevels(lps.MIN2$fusionName)

lps.MIN3 <- unique(do.call(rbind, lapply(list(lpsMSE, lpsMSI, lpsMSD, lpsMEI, lpsMED, lpsSEI, lpsSED, lpsMID, lpsIED), function(X) X[,c("fusionName", "case")])))
lps.MIN3 <- merge(lps.MIN3, lps.ee, all.x = TRUE, by = c("fusionName"))
lps.MIN3 <- merge(lps.MIN3, lps.map, all.x = TRUE, by = c("fusionName"))
lps.MIN3 <- merge(lps.MIN3, lps.star, all.x = TRUE, by = c("fusionName"))
lps.MIN3 <- merge(lps.MIN3, lps.integrate, all.x = TRUE, by = c("fusionName"))
lps.MIN3 <- merge(lps.MIN3, lps.defuse, all.x = TRUE, by = c("fusionName"))
## lps.MIN3$fusionName <- droplevels(lps.MIN3$fusionName)

min2Count <- countFusions(lps.MIN2, caseID = "case.x")
min3Count <- countFusions(lps.MIN3, caseID = "case.x")
#min4Count <- countFusions(lps.MIN4, caseID = "case.x")


pdf("LPS/LPS_MIN2_MIN3.pdf", width = 110/25.4, height = 297/25.4)
ggplot(min2Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10), legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least two tools") +
    scale_fill_manual(values = fusionColors)
ggplot(min3Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least three tools") +
    scale_fill_manual(values = fusionColors) 
dev.off()


write.table(lps.MIN2, file = "LPS/LPS_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lps.MIN3, file = "LPS/LPS_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##write.table(lps4M, file = "LPS/LPS_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("LPS/fusions.lps.rda")

load("LPS/fusions.lps.rda")
## java -jar /opt/src/oncofuse-1.1.1/Oncofuse.jar -a hg19 LPS/lps_min2_oncofuse.input.txt 'coord' EPI LPS/lps_min2_oncofuse.output.txt

Onc <- do.call(rbind, lapply(list.files(pattern = "merged.*oncofuse.output.txt"), read.delim, header = TRUE))
Onc$fusionName <- paste(Onc$X5_FPG_GENE_NAME, Onc$X3_FPG_GENE_NAME, "--")

lpsOnc <- rbind(lps.ee.recurrent[, c("program", "case", "fusionName")], lps.map.recurrent[, c("program", "case", "fusionName")], lps.star.recurrent[, c("program", "case", "fusionName")], lps.ifr.recurrent[, c("program", "case", "fusionName")])

lpsOnc <- merge(lpsOnc, Onc, by = "fusionName")

pdf("LPS/LPS_Fusions_DrivProbability.pdf", width = 148/25.4, height = 210/25.4, useDingbats = FALSE)
ggplot(lpsOnc, aes(case.x, fusionName)) + geom_tile(aes(fill = DRIVER_PROB), width = 0.95, height = 0.95) + ylab("Predicted Fusion") + xlab(paste(project, "Case"))
dev.off()


##recFtx <- c("HMGA2--SMUG1", "CCDC169--XRCC2", "CS--", "HELB--", "RP11−96H19.1−−CCDC38", "HMGA2−−RP11−834C11.5")

##recOnc <- Onc[as.character(Onc$fusionName) %in% recFtx, ]
