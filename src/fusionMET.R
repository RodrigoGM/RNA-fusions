## working with woffian tumors

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")
## load("MET/fusions.met.rda")

project <- "Metaplastic Carcinoma"
outDir <- "MET"
if( ! dir.exists(outDir)) { dir.create(file.path(outDir)) }

## cases to include
metInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project]))


## extract MET from individual filtered fusions
met.star.filtered <- star.filtered[star.filtered$project == project, ]
met.star.filtered <- met.star.filtered[!is.na(met.star.filtered$case),]
met.star.filtered <- met.star.filtered[met.star.filtered$case %in% metInc, ]

met.map.filtered <- map.filtered[map.filtered$project == project, ]
met.map.filtered <- met.map.filtered[!is.na(met.map.filtered$case),]
met.map.filtered <- met.map.filtered[met.map.filtered$case %in% metInc, ]

met.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
met.ee.filtered <- met.ee.filtered[!is.na(met.ee.filtered$case),]
met.ee.filtered <- met.ee.filtered[met.ee.filtered$case %in% metInc, ]

met.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
met.integrate.filtered <- met.integrate.filtered[!is.na(met.integrate.filtered$case),]
met.integrate.filtered <- met.integrate.filtered[met.integrate.filtered$case %in% metInc, ]

## extract MET from recurrent filtered fusions
met.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
met.star.recurrent <- met.star.recurrent[!is.na(met.star.recurrent$case),]
met.star.recurrent <- met.star.recurrent[met.star.recurrent$case %in% metInc, ]

met.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
met.map.recurrent <- met.map.recurrent[!is.na(met.map.recurrent$case),]
met.map.recurrent <- met.map.recurrent[met.map.recurrent$case %in% metInc, ]

met.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
met.ee.recurrent <- met.ee.recurrent[!is.na(met.ee.recurrent$case),]
met.ee.recurrent <- met.ee.recurrent[met.ee.recurrent$case %in% metInc, ]

met.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
met.ifr.recurrent <- met.ifr.recurrent[!is.na(met.ifr.recurrent$case),]
met.ifr.recurrent <- met.ifr.recurrent[met.ifr.recurrent$case %in% metInc, ]

## pair-wise interesections
metSM <- merge(met.star.filtered, met.map.filtered, by = c("fusionName", "case"))
metSI <- merge(met.star.filtered, met.integrate.filtered, by = c("fusionName", "case"))
metSE <- merge(met.star.filtered, met.ee.filtered, by = c("fusionName", "case"))
metME <- merge(met.map.filtered, met.ee.filtered, by = c("fusionName", "case"))
metMI <- merge(met.map.filtered, met.integrate.filtered, by = c("fusionName", "case"))
metEI <- merge(met.ee.filtered, met.integrate.filtered, by = c("fusionName", "case"))

## three-way interesections
metMSE <- merge(metSM, met.ee.filtered, by = c("fusionName", "case"))
metMSI <- merge(metSM, met.integrate.filtered, by = c("fusionName", "case"))
metMEI <- merge(metME, met.integrate.filtered, by = c("fusionName", "case"))
metSEI <- merge(metSE, met.integrate.filtered, by = c("fusionName", "case"))

## all programs
met4M <- merge(metSM, metEI, by = c("fusionName", "case"))

## write individual intersections
write.table(metSM, file = "MET/MET_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metSE, file = "MET/MET_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metME, file = "MET/MET_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metMSE, file = "MET/MET_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metMSI, file = "MET/MET_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metMEI, file = "MET/MET_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metSEI, file = "MET/MET_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(met4M, file = "MET/MET_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(met.star.recurrent, file = "MET/MET_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(met.map.recurrent, file = "MET/MET_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(met.ee.recurrent, file = "MET/MET_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(met.ifr.recurrent, file = "MET/MET_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## count individual filtered by programs
metSCount <- countFusions(met.star.filtered, "fusionName", "case")
metMCount <- countFusions(met.map.filtered, "fusionName", "case")
metECount <- countFusions(met.ee.filtered, "fusionName", "case")
metICount <- countFusions(met.integrate.filtered, "fusionName", "case")

## count recurrent filtered by programs
metSfrCount <- countFusions(met.star.recurrent, "fusionName", "case")
metMfrCount <- countFusions(met.map.recurrent, "fusionName", "case")
metEfrCount <- countFusions(met.ee.recurrent, "fusionName", "case")
metIfrCount <- countFusions(met.ifr.recurrent, "fusionName", "case")

## count intersections
metSMCounts <- countFusions(metSM, "fusionName", "case")
metSECounts <- countFusions(metSE, "fusionName", "case")
metSICounts <- countFusions(metSI, "fusionName", "case")
metMECounts <- countFusions(metME, "fusionName", "case")
metMICounts <- countFusions(metMI, "fusionName", "case")
metEICounts <- countFusions(metEI, "fusionName", "case")
metCounts <- countFusions(met4M, "fusionName", "case")



pdf("MET/MET_byProgram.pdf", width = 150/25.4, height = 297/25.4)

ggplot(metSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(metMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(metECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(metICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

dev.off()



pdf("MET/MET_Recurrent_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(metSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(metMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(metEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(metIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


pdf("MET/MET_predictedFusions.pdf", width = 105/25.4, height = 297/25.4)
ggplot(metSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(metSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(metMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(metCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


met.MIN2 <- unique(do.call(rbind, lapply(list(metEI, metME, metMI, metSE, metSI, metSM), function(X) X[,c("fusionName", "case")])))
met.MIN2 <- merge(met.MIN2, met.ee.filtered, all.x = TRUE, by = c("fusionName"))
met.MIN2 <- merge(met.MIN2, met.map.filtered, all.x = TRUE, by = c("fusionName"))
met.MIN2 <- merge(met.MIN2, met.star.filtered, all.x = TRUE, by = c("fusionName"))
met.MIN2 <- merge(met.MIN2, met.integrate.filtered, all.x = TRUE, by = c("fusionName"))
##met.MIN2$fusionName <- droplevels(met.MIN2$fusionName)

met.MIN3 <- unique(do.call(rbind, lapply(list(metMEI, metMSE, metMSI, metSEI), function(X) X[,c("fusionName", "case")])))
met.MIN3 <- merge(met.MIN3, met.ee.filtered, all.x = TRUE, by = c("fusionName"))
met.MIN3 <- merge(met.MIN3, met.map.filtered, all.x = TRUE, by = c("fusionName"))
met.MIN3 <- merge(met.MIN3, met.star.filtered, all.x = TRUE, by = c("fusionName"))
met.MIN3 <- merge(met.MIN3, met.integrate.filtered, all.x = TRUE, by = c("fusionName"))
## met.MIN3$fusionName <- droplevels(met.MIN3$fusionName)

min2Count <- countFusions(met.MIN2, caseID = "case.x")
min3Count <- countFusions(met.MIN3, caseID = "case.x")
#min4Count <- countFusions(met.MIN4, caseID = "case.x")


pdf("MET/MET_MIN2_MIN3_MIN4.pdf", width = 110/25.4, height = 297/25.4)
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
ggplot(metCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in all four tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


write.table(met.MIN2, file = "MET/MET_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(met.MIN3, file = "MET/MET_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(met4M, file = "MET/MET_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("MET/fusions.met.rda")
