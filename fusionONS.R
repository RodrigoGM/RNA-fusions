## working with salivary ACC

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

project <- "Olfactory Neuroblastoma"
outDir <- "ONS"

## cases to include
onsInc <-  unique(c(star$case[star$project == project],
                    map$case[map$project == project],
                    ee$case[ee$project == project],
                    integrate$case[integrate$project == project]))


## extract ONS from individual filtered fusions
ons.star.filtered <- star.filtered[star.filtered$project == project, ]
ons.star.filtered <- ons.star.filtered[!is.na(ons.star.filtered$case),]
ons.star.filtered <- ons.star.filtered[ons.star.filtered$case %in% onsInc, ]

ons.map.filtered <- map.filtered[map.filtered$project == project, ]
ons.map.filtered <- ons.map.filtered[!is.na(ons.map.filtered$case),]
ons.map.filtered <- ons.map.filtered[ons.map.filtered$case %in% onsInc, ]

ons.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
ons.ee.filtered <- ons.ee.filtered[!is.na(ons.ee.filtered$case),]
ons.ee.filtered <- ons.ee.filtered[ons.ee.filtered$case %in% onsInc, ]

ons.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
ons.integrate.filtered <- ons.integrate.filtered[!is.na(ons.integrate.filtered$case),]
ons.integrate.filtered <- ons.integrate.filtered[ons.integrate.filtered$case %in% onsInc, ]

## extract ONS from recurrent filtered fusions
ons.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
ons.star.recurrent <- ons.star.recurrent[!is.na(ons.star.recurrent$case),]
ons.star.recurrent <- ons.star.recurrent[ons.star.recurrent$case %in% onsInc, ]

ons.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
ons.map.recurrent <- ons.map.recurrent[!is.na(ons.map.recurrent$case),]
ons.map.recurrent <- ons.map.recurrent[ons.map.recurrent$case %in% onsInc, ]

ons.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
ons.ee.recurrent <- ons.ee.recurrent[!is.na(ons.ee.recurrent$case),]
ons.ee.recurrent <- ons.ee.recurrent[ons.ee.recurrent$case %in% onsInc, ]

ons.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
ons.ifr.recurrent <- ons.ifr.recurrent[!is.na(ons.ifr.recurrent$case),]
ons.ifr.recurrent <- ons.ifr.recurrent[ons.ifr.recurrent$case %in% onsInc, ]

## pair-wise interesections
onsSM <- merge(ons.star.filtered, ons.map.filtered, by = c("fusionName", "case"))
onsSI <- merge(ons.star.filtered, ons.integrate.filtered, by = c("fusionName", "case"))
onsSE <- merge(ons.star.filtered, ons.ee.filtered, by = c("fusionName", "case"))
onsME <- merge(ons.map.filtered, ons.ee.filtered, by = c("fusionName", "case"))
onsMI <- merge(ons.map.filtered, ons.integrate.filtered, by = c("fusionName", "case"))
onsEI <- merge(ons.ee.filtered, ons.integrate.filtered, by = c("fusionName", "case"))

## three-way interesections
onsMSE <- merge(onsSM, ons.ee.filtered, by = c("fusionName", "case"))
onsMSI <- merge(onsSM, ons.integrate.filtered, by = c("fusionName", "case"))
onsMEI <- merge(onsME, ons.integrate.filtered, by = c("fusionName", "case"))
onsSEI <- merge(onsSE, ons.integrate.filtered, by = c("fusionName", "case"))

## all programs
ons4M <- merge(onsSM, onsEI, by = c("fusionName", "case"))

## write individual intersections
write.table(onsSM, file = "ONS/ONS_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(onsSE, file = "ONS/ONS_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(onsME, file = "ONS/ONS_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(onsMSE, file = "ONS/ONS_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(onsMSI, file = "ONS/ONS_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(onsMEI, file = "ONS/ONS_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(onsSEI, file = "ONS/ONS_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ons4M, file = "ONS/ONS_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(ons.star.recurrent, file = "ONS/ONS_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ons.map.recurrent, file = "ONS/ONS_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ons.ee.recurrent, file = "ONS/ONS_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ons.ifr.recurrent, file = "ONS/ONS_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## count individual filtered by programs
onsSCount <- countFusions(ons.star.filtered, "fusionName", "case")
onsMCount <- countFusions(ons.map.filtered, "fusionName", "case")
onsECount <- countFusions(ons.ee.filtered, "fusionName", "case")
onsICount <- countFusions(ons.integrate.filtered, "fusionName", "case")
## count recurrent filtered by programs
onsSfrCount <- countFusions(ons.star.recurrent, "fusionName", "case")
onsMfrCount <- countFusions(ons.map.recurrent, "fusionName", "case")
onsEfrCount <- countFusions(ons.ee.recurrent, "fusionName", "case")
onsIfrCount <- countFusions(ons.ifr.recurrent, "fusionName", "case")

## count intersections
onsSMCounts <- countFusions(onsSM, "fusionName", "case")
onsSECounts <- countFusions(onsSE, "fusionName", "case")
onsSICounts <- countFusions(onsSI, "fusionName", "case")
onsMECounts <- countFusions(onsME, "fusionName", "case")
onsMICounts <- countFusions(onsMI, "fusionName", "case")
onsEICounts <- countFusions(onsEI, "fusionName", "case")
onsCounts <- countFusions(ons4M, "fusionName", "case")



pdf("ONS/ONS_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(onsSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(onsMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(onsECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(onsICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()



pdf("ONS/ONS_Recurrent_byProgram.pdf", width = 298/25.4, height = 420/25.4)
ggplot(onsSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#F8766D", "#CFCFCF"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Predicted", "Not Present"))
ggplot(onsMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = c("#F8766D", "#CFCFCF"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Predicted", "Not Present"))
ggplot(onsEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = c("#F8766D", "#CFCFCF"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Predicted", "Not Present"))
ggplot(onsIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = c("#F8766D", "#CFCFCF"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Predicted", "Not Present"))
dev.off()


pdf("ONS/ONS_predictedFusions.pdf", width = 105/25.4, height = 105/25.4)
ggplot(onsSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(onsSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(onsMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(onsCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


ons.MIN2 <- unique(do.call(rbind, lapply(list(onsEI, onsME, onsMI, onsSE, onsSI, onsSM), function(X) X[,c("fusionName", "case")])))
ons.MIN2 <- merge(ons.MIN2, ons.ee.filtered, all.x = TRUE, by = c("fusionName"))
ons.MIN2 <- merge(ons.MIN2, ons.map.filtered, all.x = TRUE, by = c("fusionName"))
ons.MIN2 <- merge(ons.MIN2, ons.star.filtered, all.x = TRUE, by = c("fusionName"))
ons.MIN2 <- merge(ons.MIN2, ons.integrate.filtered, all.x = TRUE, by = c("fusionName"))
##ons.MIN2$fusionName <- droplevels(ons.MIN2$fusionName)

ons.MIN3 <- unique(do.call(rbind, lapply(list(onsMEI, onsMSE, onsMSI, onsSEI), function(X) X[,c("fusionName", "case")])))
ons.MIN3 <- merge(ons.MIN3, ons.ee.filtered, all.x = TRUE, by = c("fusionName"))
ons.MIN3 <- merge(ons.MIN3, ons.map.filtered, all.x = TRUE, by = c("fusionName"))
ons.MIN3 <- merge(ons.MIN3, ons.star.filtered, all.x = TRUE, by = c("fusionName"))
ons.MIN3 <- merge(ons.MIN3, ons.integrate.filtered, all.x = TRUE, by = c("fusionName"))
## ons.MIN3$fusionName <- droplevels(ons.MIN3$fusionName)

min2Count <- countFusions(ons.MIN2, caseID = "case.x")
min3Count <- countFusions(ons.MIN3, caseID = "case.x")
#min4Count <- countFusions(ons.MIN4, caseID = "case.x")


pdf("ONS/ONS_MIN2_MIN3_MIN4.pdf", width = 110/25.4, height = 297/25.4)
ggplot(min2Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10), legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least two tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(min3Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least three tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(onsCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in all four tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


write.table(ons.MIN2, file = "ONS/ONS_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ons.MIN3, file = "ONS/ONS_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ons4M, file = "ONS/ONS_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("ONS/fusions.ons.rda")

load("ONS/fusions.ons.rda")
## java -jar /opt/src/oncofuse-1.1.1/Oncofuse.jar -a hg19 ONS/ons_min2_oncofuse.input.txt 'coord' EPI ONS/ons_min2_oncofuse.output.txt

Onc <- read.delim("all.oncofuse.annotated.fusions.txt")

onsOnc <- rbind(ons.ee.recurrent[, c("program", "case", "fusionName")], ons.map.recurrent[, c("program", "case", "fusionName")], ons.star.recurrent[, c("program", "case", "fusionName")], ons.ifr.recurrent[, c("program", "case", "fusionName")])

onsOnc <- merge(onsOnc, Onc, by = "fusionName")

pdf("ONS/ONS_Fusions_DrivProbability.pdf", width = 298/25.4, height = 298/25.4, useDingbats = FALSE)
ggplot(onsOnc, aes(case, fusionName)) + geom_tile(aes(fill = DRIVER_PROB), width = 0.95, height = 0.95) + ylab("Predicted Fusion") + xlab("Liposarcoma Case")
dev.off()


recFtx <- c("HMGA2--SMUG1", "CCDC169--XRCC2", "CS--", "HELB--", "RP11−96H19.1−−CCDC38", "HMGA2−−RP11−834C11.5")

recOnc <- Onc[as.character(Onc$fusionName) %in% recFtx, ]
