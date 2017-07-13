## working with woffian tumors

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")
## load("WOLF/fusions.wolf.rda")

project <- "Wolffian"
outDir <- "WOLF"

## cases to include
wolfInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project]))


## extract WOLF from individual filtered fusions
wolf.star.filtered <- star.filtered[star.filtered$project == project, ]
wolf.star.filtered <- wolf.star.filtered[!is.na(wolf.star.filtered$case),]
wolf.star.filtered <- wolf.star.filtered[wolf.star.filtered$case %in% wolfInc, ]

wolf.map.filtered <- map.filtered[map.filtered$project == project, ]
wolf.map.filtered <- wolf.map.filtered[!is.na(wolf.map.filtered$case),]
wolf.map.filtered <- wolf.map.filtered[wolf.map.filtered$case %in% wolfInc, ]

wolf.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
wolf.ee.filtered <- wolf.ee.filtered[!is.na(wolf.ee.filtered$case),]
wolf.ee.filtered <- wolf.ee.filtered[wolf.ee.filtered$case %in% wolfInc, ]

wolf.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
wolf.integrate.filtered <- wolf.integrate.filtered[!is.na(wolf.integrate.filtered$case),]
wolf.integrate.filtered <- wolf.integrate.filtered[wolf.integrate.filtered$case %in% wolfInc, ]

## extract WOLF from recurrent filtered fusions
wolf.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
wolf.star.recurrent <- wolf.star.recurrent[!is.na(wolf.star.recurrent$case),]
wolf.star.recurrent <- wolf.star.recurrent[wolf.star.recurrent$case %in% wolfInc, ]

wolf.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
wolf.map.recurrent <- wolf.map.recurrent[!is.na(wolf.map.recurrent$case),]
wolf.map.recurrent <- wolf.map.recurrent[wolf.map.recurrent$case %in% wolfInc, ]

wolf.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
wolf.ee.recurrent <- wolf.ee.recurrent[!is.na(wolf.ee.recurrent$case),]
wolf.ee.recurrent <- wolf.ee.recurrent[wolf.ee.recurrent$case %in% wolfInc, ]

wolf.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
wolf.ifr.recurrent <- wolf.ifr.recurrent[!is.na(wolf.ifr.recurrent$case),]
wolf.ifr.recurrent <- wolf.ifr.recurrent[wolf.ifr.recurrent$case %in% wolfInc, ]

## pair-wise interesections
wolfSM <- merge(wolf.star.filtered, wolf.map.filtered, by = c("fusionName", "case"))
wolfSI <- merge(wolf.star.filtered, wolf.integrate.filtered, by = c("fusionName", "case"))
wolfSE <- merge(wolf.star.filtered, wolf.ee.filtered, by = c("fusionName", "case"))
wolfME <- merge(wolf.map.filtered, wolf.ee.filtered, by = c("fusionName", "case"))
wolfMI <- merge(wolf.map.filtered, wolf.integrate.filtered, by = c("fusionName", "case"))
wolfEI <- merge(wolf.ee.filtered, wolf.integrate.filtered, by = c("fusionName", "case"))

## three-way interesections
wolfMSE <- merge(wolfSM, wolf.ee.filtered, by = c("fusionName", "case"))
wolfMSI <- merge(wolfSM, wolf.integrate.filtered, by = c("fusionName", "case"))
wolfMEI <- merge(wolfME, wolf.integrate.filtered, by = c("fusionName", "case"))
wolfSEI <- merge(wolfSE, wolf.integrate.filtered, by = c("fusionName", "case"))

## all programs
wolf4M <- merge(wolfSM, wolfEI, by = c("fusionName", "case"))

## write individual intersections
write.table(wolfSM, file = "WOLF/WOLF_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolfSE, file = "WOLF/WOLF_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolfME, file = "WOLF/WOLF_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolfMSE, file = "WOLF/WOLF_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolfMSI, file = "WOLF/WOLF_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolfMEI, file = "WOLF/WOLF_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolfSEI, file = "WOLF/WOLF_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolf4M, file = "WOLF/WOLF_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(wolf.star.recurrent, file = "WOLF/WOLF_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolf.map.recurrent, file = "WOLF/WOLF_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolf.ee.recurrent, file = "WOLF/WOLF_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolf.ifr.recurrent, file = "WOLF/WOLF_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## count individual filtered by programs
wolfSCount <- countFusions(wolf.star.filtered, "fusionName", "case")
wolfMCount <- countFusions(wolf.map.filtered, "fusionName", "case")
wolfECount <- countFusions(wolf.ee.filtered, "fusionName", "case")
wolfICount <- countFusions(wolf.integrate.filtered, "fusionName", "case")

## count recurrent filtered by programs
wolfSfrCount <- countFusions(wolf.star.recurrent, "fusionName", "case")
wolfMfrCount <- countFusions(wolf.map.recurrent, "fusionName", "case")
wolfEfrCount <- countFusions(wolf.ee.recurrent, "fusionName", "case")
wolfIfrCount <- countFusions(wolf.ifr.recurrent, "fusionName", "case")

## count intersections
wolfSMCounts <- countFusions(wolfSM, "fusionName", "case")
wolfSECounts <- countFusions(wolfSE, "fusionName", "case")
wolfSICounts <- countFusions(wolfSI, "fusionName", "case")
wolfMECounts <- countFusions(wolfME, "fusionName", "case")
wolfMICounts <- countFusions(wolfMI, "fusionName", "case")
wolfEICounts <- countFusions(wolfEI, "fusionName", "case")
wolfCounts <- countFusions(wolf4M, "fusionName", "case")



pdf("WOLF/WOLF_byProgram.pdf", width = 150/25.4, height = 297/25.4)

ggplot(wolfSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(wolfMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(wolfECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(wolfICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

dev.off()



pdf("WOLF/WOLF_Recurrent_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(wolfSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#F8766D", "#CFCFCF"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Predicted", "Not Present"))
ggplot(wolfMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(wolfEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(wolfIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


pdf("WOLF/WOLF_predictedFusions.pdf", width = 105/25.4, height = 297/25.4)
ggplot(wolfSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(wolfSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(wolfMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(wolfCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


wolf.MIN2 <- unique(do.call(rbind, lapply(list(wolfEI, wolfME, wolfMI, wolfSE, wolfSI, wolfSM), function(X) X[,c("fusionName", "case")])))
wolf.MIN2 <- merge(wolf.MIN2, wolf.ee.filtered, all.x = TRUE, by = c("fusionName"))
wolf.MIN2 <- merge(wolf.MIN2, wolf.map.filtered, all.x = TRUE, by = c("fusionName"))
wolf.MIN2 <- merge(wolf.MIN2, wolf.star.filtered, all.x = TRUE, by = c("fusionName"))
wolf.MIN2 <- merge(wolf.MIN2, wolf.integrate.filtered, all.x = TRUE, by = c("fusionName"))
##wolf.MIN2$fusionName <- droplevels(wolf.MIN2$fusionName)

wolf.MIN3 <- unique(do.call(rbind, lapply(list(wolfMEI, wolfMSE, wolfMSI, wolfSEI), function(X) X[,c("fusionName", "case")])))
wolf.MIN3 <- merge(wolf.MIN3, wolf.ee.filtered, all.x = TRUE, by = c("fusionName"))
wolf.MIN3 <- merge(wolf.MIN3, wolf.map.filtered, all.x = TRUE, by = c("fusionName"))
wolf.MIN3 <- merge(wolf.MIN3, wolf.star.filtered, all.x = TRUE, by = c("fusionName"))
wolf.MIN3 <- merge(wolf.MIN3, wolf.integrate.filtered, all.x = TRUE, by = c("fusionName"))
## wolf.MIN3$fusionName <- droplevels(wolf.MIN3$fusionName)

min2Count <- countFusions(wolf.MIN2, caseID = "case.x")
min3Count <- countFusions(wolf.MIN3, caseID = "case.x")
#min4Count <- countFusions(wolf.MIN4, caseID = "case.x")


pdf("WOLF/WOLF_MIN2_MIN3_MIN4.pdf", width = 110/25.4, height = 297/25.4)
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
ggplot(wolfCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in all four tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


write.table(wolf.MIN2, file = "WOLF/WOLF_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolf.MIN3, file = "WOLF/WOLF_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wolf4M, file = "WOLF/WOLF_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("WOLF/fusions.wolf.rda")
