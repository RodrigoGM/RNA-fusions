## working with thyroid adenomas

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

project <- "Trabecular Adenoma Thyroid"
outDir <- "THT"

## cases to include
thtInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project]))

## extract THT from individual filtered fusions
tht.star.filtered <- star.filtered[star.filtered$project == project, ]
tht.star.filtered <- tht.star.filtered[!is.na(tht.star.filtered$case),]
tht.star.filtered <- tht.star.filtered[tht.star.filtered$case %in% thtInc, ]

tht.map.filtered <- map.filtered[map.filtered$project == project, ]
tht.map.filtered <- tht.map.filtered[!is.na(tht.map.filtered$case),]
tht.map.filtered <- tht.map.filtered[tht.map.filtered$case %in% thtInc, ]

tht.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
tht.ee.filtered <- tht.ee.filtered[!is.na(tht.ee.filtered$case),]
tht.ee.filtered <- tht.ee.filtered[tht.ee.filtered$case %in% thtInc, ]

tht.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
tht.integrate.filtered <- tht.integrate.filtered[!is.na(tht.integrate.filtered$case),]
tht.integrate.filtered <- tht.integrate.filtered[tht.integrate.filtered$case %in% thtInc, ]

## extract THT from recurrent filtered fusions
tht.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
tht.star.recurrent <- tht.star.recurrent[!is.na(tht.star.recurrent$case),]
tht.star.recurrent <- tht.star.recurrent[tht.star.recurrent$case %in% thtInc, ]

tht.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
tht.map.recurrent <- tht.map.recurrent[!is.na(tht.map.recurrent$case),]
tht.map.recurrent <- tht.map.recurrent[tht.map.recurrent$case %in% thtInc, ]

tht.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
tht.ee.recurrent <- tht.ee.recurrent[!is.na(tht.ee.recurrent$case),]
tht.ee.recurrent <- tht.ee.recurrent[tht.ee.recurrent$case %in% thtInc, ]

tht.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
tht.ifr.recurrent <- tht.ifr.recurrent[!is.na(tht.ifr.recurrent$case),]
tht.ifr.recurrent <- tht.ifr.recurrent[tht.ifr.recurrent$case %in% thtInc, ]

## pair-wise interesections
thtSM <- merge(tht.star.filtered, tht.map.filtered, by = c("fusionName", "case"))
thtSI <- merge(tht.star.filtered, tht.integrate.filtered, by = c("fusionName", "case"))
thtSE <- merge(tht.star.filtered, tht.ee.filtered, by = c("fusionName", "case"))
thtME <- merge(tht.map.filtered, tht.ee.filtered, by = c("fusionName", "case"))
thtMI <- merge(tht.map.filtered, tht.integrate.filtered, by = c("fusionName", "case"))
thtEI <- merge(tht.ee.filtered, tht.integrate.filtered, by = c("fusionName", "case"))

## three-way interesections
thtMSE <- merge(thtSM, tht.ee.filtered, by = c("fusionName", "case"))
thtMSI <- merge(thtSM, tht.integrate.filtered, by = c("fusionName", "case"))
thtMEI <- merge(thtME, tht.integrate.filtered, by = c("fusionName", "case"))
thtSEI <- merge(thtSE, tht.integrate.filtered, by = c("fusionName", "case"))

## all programs
tht4M <- merge(thtSM, thtEI, by = c("fusionName", "case"))

## write individual intersections
write.table(thtSM, file = "THT/THT_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(thtSE, file = "THT/THT_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(thtME, file = "THT/THT_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(thtMSE, file = "THT/THT_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(thtMSI, file = "THT/THT_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(thtMEI, file = "THT/THT_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(thtSEI, file = "THT/THT_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tht4M, file = "THT/THT_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(tht.star.recurrent, file = "THT/THT_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tht.map.recurrent, file = "THT/THT_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tht.ee.recurrent, file = "THT/THT_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tht.ifr.recurrent, file = "THT/THT_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## count individual filtered by programs
thtSCount <- countFusions(tht.star.filtered, "fusionName", "case")
thtMCount <- countFusions(tht.map.filtered, "fusionName", "case")
thtECount <- countFusions(tht.ee.filtered, "fusionName", "case")
thtICount <- countFusions(tht.integrate.filtered, "fusionName", "case")
## count recurrent filtered by programs
thtSfrCount <- countFusions(tht.star.recurrent, "fusionName", "case")
thtMfrCount <- countFusions(tht.map.recurrent, "fusionName", "case")
thtEfrCount <- countFusions(tht.ee.recurrent, "fusionName", "case")
thtIfrCount <- countFusions(tht.ifr.recurrent, "fusionName", "case")

## count intersections
thtSMCounts <- countFusions(thtSM, "fusionName", "case")
thtSECounts <- countFusions(thtSE, "fusionName", "case")
thtSICounts <- countFusions(thtSI, "fusionName", "case")
thtMECounts <- countFusions(thtME, "fusionName", "case")
thtMICounts <- countFusions(thtMI, "fusionName", "case")
thtEICounts <- countFusions(thtEI, "fusionName", "case")
thtCounts <- countFusions(tht4M, "fusionName", "case")



pdf("THT/THT_byProgram.pdf", width = 150/25.4, height = 297/25.4)

ggplot(thtSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(thtMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(thtECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(thtICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

dev.off()



pdf("THT/THT_Recurrent_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(thtSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(thtMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(thtEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(thtIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


pdf("THT/THT_predictedFusions.pdf", width = 105/25.4, height = 297/25.4)
ggplot(thtSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(thtSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(thtMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(thtCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


tht.MIN2 <- unique(do.call(rbind, lapply(list(thtEI, thtME, thtMI, thtSE, thtSI, thtSM), function(X) X[,c("fusionName", "case")])))
tht.MIN2 <- merge(tht.MIN2, tht.ee.filtered, all.x = TRUE, by = c("fusionName"))
tht.MIN2 <- merge(tht.MIN2, tht.map.filtered, all.x = TRUE, by = c("fusionName"))
tht.MIN2 <- merge(tht.MIN2, tht.star.filtered, all.x = TRUE, by = c("fusionName"))
tht.MIN2 <- merge(tht.MIN2, tht.integrate.filtered, all.x = TRUE, by = c("fusionName"))
##tht.MIN2$fusionName <- droplevels(tht.MIN2$fusionName)

tht.MIN3 <- unique(do.call(rbind, lapply(list(thtMEI, thtMSE, thtMSI, thtSEI), function(X) X[,c("fusionName", "case")])))
tht.MIN3 <- merge(tht.MIN3, tht.ee.filtered, all.x = TRUE, by = c("fusionName"))
tht.MIN3 <- merge(tht.MIN3, tht.map.filtered, all.x = TRUE, by = c("fusionName"))
tht.MIN3 <- merge(tht.MIN3, tht.star.filtered, all.x = TRUE, by = c("fusionName"))
tht.MIN3 <- merge(tht.MIN3, tht.integrate.filtered, all.x = TRUE, by = c("fusionName"))
## tht.MIN3$fusionName <- droplevels(tht.MIN3$fusionName)

min2Count <- countFusions(tht.MIN2, caseID = "case.x")
min3Count <- countFusions(tht.MIN3, caseID = "case.x")
#min4Count <- countFusions(tht.MIN4, caseID = "case.x")


pdf("THT/THT_MIN2_MIN3_MIN4.pdf", width = 110/25.4, height = 297/25.4)
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
ggplot(thtCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in all four tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


write.table(tht.MIN2, file = "THT/THT_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tht.MIN3, file = "THT/THT_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tht4M, file = "THT/THT_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("THT/fusions.tht.rda")
