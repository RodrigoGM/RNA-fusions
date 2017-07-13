## working with papillary adenomas

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")
## load("PAP/fusions.pap.rda")

project <- "Papillary"
outDir <- "PAP"

## cases to include
papInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project]))


## extract PAP from individual filtered fusions
pap.star.filtered <- star.filtered[star.filtered$project == project, ]
pap.star.filtered <- pap.star.filtered[!is.na(pap.star.filtered$case),]
pap.star.filtered <- pap.star.filtered[pap.star.filtered$case %in% papInc, ]

pap.map.filtered <- map.filtered[map.filtered$project == project, ]
pap.map.filtered <- pap.map.filtered[!is.na(pap.map.filtered$case),]
pap.map.filtered <- pap.map.filtered[pap.map.filtered$case %in% papInc, ]

pap.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
pap.ee.filtered <- pap.ee.filtered[!is.na(pap.ee.filtered$case),]
pap.ee.filtered <- pap.ee.filtered[pap.ee.filtered$case %in% papInc, ]

pap.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
pap.integrate.filtered <- pap.integrate.filtered[!is.na(pap.integrate.filtered$case),]
pap.integrate.filtered <- pap.integrate.filtered[pap.integrate.filtered$case %in% papInc, ]

## extract PAP from recurrent filtered fusions
pap.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
pap.star.recurrent <- pap.star.recurrent[!is.na(pap.star.recurrent$case),]
pap.star.recurrent <- pap.star.recurrent[pap.star.recurrent$case %in% papInc, ]

pap.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
pap.map.recurrent <- pap.map.recurrent[!is.na(pap.map.recurrent$case),]
pap.map.recurrent <- pap.map.recurrent[pap.map.recurrent$case %in% papInc, ]

pap.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
pap.ee.recurrent <- pap.ee.recurrent[!is.na(pap.ee.recurrent$case),]
pap.ee.recurrent <- pap.ee.recurrent[pap.ee.recurrent$case %in% papInc, ]

pap.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
pap.ifr.recurrent <- pap.ifr.recurrent[!is.na(pap.ifr.recurrent$case),]
pap.ifr.recurrent <- pap.ifr.recurrent[pap.ifr.recurrent$case %in% papInc, ]

## pair-wise interesections
papSM <- merge(pap.star.filtered, pap.map.filtered, by = c("fusionName", "case"))
papSI <- merge(pap.star.filtered, pap.integrate.filtered, by = c("fusionName", "case"))
papSE <- merge(pap.star.filtered, pap.ee.filtered, by = c("fusionName", "case"))
papME <- merge(pap.map.filtered, pap.ee.filtered, by = c("fusionName", "case"))
papMI <- merge(pap.map.filtered, pap.integrate.filtered, by = c("fusionName", "case"))
papEI <- merge(pap.ee.filtered, pap.integrate.filtered, by = c("fusionName", "case"))

## three-way interesections
papMSE <- merge(papSM, pap.ee.filtered, by = c("fusionName", "case"))
papMSI <- merge(papSM, pap.integrate.filtered, by = c("fusionName", "case"))
papMEI <- merge(papME, pap.integrate.filtered, by = c("fusionName", "case"))
papSEI <- merge(papSE, pap.integrate.filtered, by = c("fusionName", "case"))

## all programs
pap4M <- merge(papSM, papEI, by = c("fusionName", "case"))

## write individual intersections
write.table(papSM, file = "PAP/PAP_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(papSE, file = "PAP/PAP_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(papME, file = "PAP/PAP_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(papMSE, file = "PAP/PAP_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(papMSI, file = "PAP/PAP_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(papMEI, file = "PAP/PAP_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(papSEI, file = "PAP/PAP_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pap4M, file = "PAP/PAP_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(pap.star.recurrent, file = "PAP/PAP_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pap.map.recurrent, file = "PAP/PAP_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pap.ee.recurrent, file = "PAP/PAP_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pap.ifr.recurrent, file = "PAP/PAP_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## count individual filtered by programs
papSCount <- countFusions(pap.star.filtered, "fusionName", "case")
papMCount <- countFusions(pap.map.filtered, "fusionName", "case")
papECount <- countFusions(pap.ee.filtered, "fusionName", "case")
papICount <- countFusions(pap.integrate.filtered, "fusionName", "case")

## count recurrent filtered by programs
papSfrCount <- countFusions(pap.star.recurrent, "fusionName", "case")
papMfrCount <- countFusions(pap.map.recurrent, "fusionName", "case")
papEfrCount <- countFusions(pap.ee.recurrent, "fusionName", "case")
papIfrCount <- countFusions(pap.ifr.recurrent, "fusionName", "case")

## count intersections
papSMCounts <- countFusions(papSM, "fusionName", "case")
papSECounts <- countFusions(papSE, "fusionName", "case")
papSICounts <- countFusions(papSI, "fusionName", "case")
papMECounts <- countFusions(papME, "fusionName", "case")
papMICounts <- countFusions(papMI, "fusionName", "case")
papEICounts <- countFusions(papEI, "fusionName", "case")
papCounts <- countFusions(pap4M, "fusionName", "case")



pdf("PAP/PAP_byProgram.pdf", width = 150/25.4, height = 297/25.4)

ggplot(papSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(papMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(papECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(papICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

dev.off()



pdf("PAP/PAP_Recurrent_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(papSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(papMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(papEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(papIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


pdf("PAP/PAP_predictedFusions.pdf", width = 105/25.4, height = 297/25.4)
ggplot(papSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(papSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(papMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(papCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


pap.MIN2 <- unique(do.call(rbind, lapply(list(papEI, papME, papMI, papSE, papSI, papSM), function(X) X[,c("fusionName", "case")])))
pap.MIN2 <- merge(pap.MIN2, pap.ee.filtered, all.x = TRUE, by = c("fusionName"))
pap.MIN2 <- merge(pap.MIN2, pap.map.filtered, all.x = TRUE, by = c("fusionName"))
pap.MIN2 <- merge(pap.MIN2, pap.star.filtered, all.x = TRUE, by = c("fusionName"))
pap.MIN2 <- merge(pap.MIN2, pap.integrate.filtered, all.x = TRUE, by = c("fusionName"))
##pap.MIN2$fusionName <- droplevels(pap.MIN2$fusionName)

pap.MIN3 <- unique(do.call(rbind, lapply(list(papMEI, papMSE, papMSI, papSEI), function(X) X[,c("fusionName", "case")])))
pap.MIN3 <- merge(pap.MIN3, pap.ee.filtered, all.x = TRUE, by = c("fusionName"))
pap.MIN3 <- merge(pap.MIN3, pap.map.filtered, all.x = TRUE, by = c("fusionName"))
pap.MIN3 <- merge(pap.MIN3, pap.star.filtered, all.x = TRUE, by = c("fusionName"))
pap.MIN3 <- merge(pap.MIN3, pap.integrate.filtered, all.x = TRUE, by = c("fusionName"))
## pap.MIN3$fusionName <- droplevels(pap.MIN3$fusionName)

min2Count <- countFusions(pap.MIN2, caseID = "case.x")
min3Count <- countFusions(pap.MIN3, caseID = "case.x")
#min4Count <- countFusions(pap.MIN4, caseID = "case.x")


pdf("PAP/PAP_MIN2_MIN3_MIN4.pdf", width = 110/25.4, height = 297/25.4)
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
ggplot(papCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in all four tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


write.table(pap.MIN2, file = "PAP/PAP_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pap.MIN3, file = "PAP/PAP_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pap4M, file = "PAP/PAP_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("PAP/fusions.pap.rda")
