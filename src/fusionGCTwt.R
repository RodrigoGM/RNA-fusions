## working with GCT

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

project <- "Granular Cell Tumor"
outDir <- "GCTwt"

## cases to include
## gctInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project]))
gctInc <- c("GCT8", "GCT11", "GCT12", "GCT14", "GCT20")

## extract GCTwt from individual filtered fusions
gctwt.star.filtered <- star.filtered[star.filtered$project == project, ]
gctwt.star.filtered <- gctwt.star.filtered[!is.na(gctwt.star.filtered$case),]
gctwt.star.filtered <- gctwt.star.filtered[gctwt.star.filtered$case %in% gctInc, ]

gctwt.map.filtered <- map.filtered[map.filtered$project == project, ]
gctwt.map.filtered <- gctwt.map.filtered[!is.na(gctwt.map.filtered$case),]
gctwt.map.filtered <- gctwt.map.filtered[gctwt.map.filtered$case %in% gctInc, ]

gctwt.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
gctwt.ee.filtered <- gctwt.ee.filtered[!is.na(gctwt.ee.filtered$case),]
gctwt.ee.filtered <- gctwt.ee.filtered[gctwt.ee.filtered$case %in% gctInc, ]

gctwt.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
gctwt.integrate.filtered <- gctwt.integrate.filtered[!is.na(gctwt.integrate.filtered$case),]
gctwt.integrate.filtered <- gctwt.integrate.filtered[gctwt.integrate.filtered$case %in% gctInc, ]

## extract GCTwt from recurrent filtered fusions
gctwt.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
gctwt.star.recurrent <- gctwt.star.recurrent[!is.na(gctwt.star.recurrent$case),]
gctwt.star.recurrent <- gctwt.star.recurrent[gctwt.star.recurrent$case %in% gctInc, ]

gctwt.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
gctwt.map.recurrent <- gctwt.map.recurrent[!is.na(gctwt.map.recurrent$case),]
gctwt.map.recurrent <- gctwt.map.recurrent[gctwt.map.recurrent$case %in% gctInc, ]

gctwt.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
gctwt.ee.recurrent <- gctwt.ee.recurrent[!is.na(gctwt.ee.recurrent$case),]
gctwt.ee.recurrent <- gctwt.ee.recurrent[gctwt.ee.recurrent$case %in% gctInc, ]

gctwt.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
gctwt.ifr.recurrent <- gctwt.ifr.recurrent[!is.na(gctwt.ifr.recurrent$case),]
gctwt.ifr.recurrent <- gctwt.ifr.recurrent[gctwt.ifr.recurrent$case %in% gctInc, ]

## pair-wise interesections
gctwtSM <- merge(gctwt.star.filtered, gctwt.map.filtered, by = c("fusionName", "case"))
gctwtSI <- merge(gctwt.star.filtered, gctwt.integrate.filtered, by = c("fusionName", "case"))
gctwtSE <- merge(gctwt.star.filtered, gctwt.ee.filtered, by = c("fusionName", "case"))
gctwtME <- merge(gctwt.map.filtered, gctwt.ee.filtered, by = c("fusionName", "case"))
gctwtMI <- merge(gctwt.map.filtered, gctwt.integrate.filtered, by = c("fusionName", "case"))
gctwtEI <- merge(gctwt.ee.filtered, gctwt.integrate.filtered, by = c("fusionName", "case"))

## three-way interesections
gctwtMSE <- merge(gctwtSM, gctwt.ee.filtered, by = c("fusionName", "case"))
gctwtMSI <- merge(gctwtSM, gctwt.integrate.filtered, by = c("fusionName", "case"))
gctwtMEI <- merge(gctwtME, gctwt.integrate.filtered, by = c("fusionName", "case"))
gctwtSEI <- merge(gctwtSE, gctwt.integrate.filtered, by = c("fusionName", "case"))

## all programs
gctwt4M <- merge(gctwtSM, gctwtEI, by = c("fusionName", "case"))

## write individual intersections
write.table(gctwtSM, file = "GCTwt/GCTwt_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwtSE, file = "GCTwt/GCTwt_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwtME, file = "GCTwt/GCTwt_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwtMSE, file = "GCTwt/GCTwt_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwtMSI, file = "GCTwt/GCTwt_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwtMEI, file = "GCTwt/GCTwt_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwtSEI, file = "GCTwt/GCTwt_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwt4M, file = "GCTwt/GCTwt_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(gctwt.star.recurrent, file = "GCTwt/GCTwt_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwt.map.recurrent, file = "GCTwt/GCTwt_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwt.ee.recurrent, file = "GCTwt/GCTwt_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwt.ifr.recurrent, file = "GCTwt/GCTwt_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## count individual filtered by programs
gctwtSCount <- countFusions(gctwt.star.filtered, "fusionName", "case")
gctwtMCount <- countFusions(gctwt.map.filtered, "fusionName", "case")
gctwtECount <- countFusions(gctwt.ee.filtered[gctwt.ee.filtered$], "fusionName", "case")
gctwtICount <- countFusions(gctwt.integrate.filtered, "fusionName", "case")

## count recurrent filtered by programs
gctwtSfrCount <- countFusions(gctwt.star.recurrent, "fusionName", "case")
gctwtMfrCount <- countFusions(gctwt.map.recurrent, "fusionName", "case")
gctwtEfrCount <- countFusions(gctwt.ee.recurrent, "fusionName", "case")
gctwtIfrCount <- countFusions(gctwt.ifr.recurrent, "fusionName", "case")

## count intersections
gctwtSMCounts <- countFusions(gctwtSM, "fusionName", "case")
gctwtSECounts <- countFusions(gctwtSE, "fusionName", "case")
gctwtSICounts <- countFusions(gctwtSI, "fusionName", "case")
gctwtMECounts <- countFusions(gctwtME, "fusionName", "case")
gctwtMICounts <- countFusions(gctwtMI, "fusionName", "case")
gctwtEICounts <- countFusions(gctwtEI, "fusionName", "case")
gctwtCounts <- countFusions(gctwt4M, "fusionName", "case")



pdf("GCTwt/GCTwt_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(gctwtSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()



pdf("GCTwt/GCTwt_Recurrent_byProgram.pdf", width = 150/25.4, height = 200/25.4)
ggplot(gctwtSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


pdf("GCTwt/GCTwt_predictedFusions.pdf", width = 105/25.4, height = 297/25.4)
ggplot(gctwtSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctwtCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


gctwt.MIN2 <- unique(do.call(rbind, lapply(list(gctwtEI, gctwtME, gctwtMI, gctwtSE, gctwtSI, gctwtSM), function(X) X[,c("fusionName", "case")])))
gctwt.MIN2 <- merge(gctwt.MIN2, gctwt.ee.filtered, all.x = TRUE, by = c("fusionName"))
gctwt.MIN2 <- merge(gctwt.MIN2, gctwt.map.filtered, all.x = TRUE, by = c("fusionName"))
gctwt.MIN2 <- merge(gctwt.MIN2, gctwt.star.filtered, all.x = TRUE, by = c("fusionName"))
gctwt.MIN2 <- merge(gctwt.MIN2, gctwt.integrate.filtered, all.x = TRUE, by = c("fusionName"))
##gctwt.MIN2$fusionName <- droplevels(gctwt.MIN2$fusionName)

gctwt.MIN3 <- unique(do.call(rbind, lapply(list(gctwtMEI, gctwtMSE, gctwtMSI, gctwtSEI), function(X) X[,c("fusionName", "case")])))
gctwt.MIN3 <- merge(gctwt.MIN3, gctwt.ee.filtered, all.x = TRUE, by = c("fusionName"))
gctwt.MIN3 <- merge(gctwt.MIN3, gctwt.map.filtered, all.x = TRUE, by = c("fusionName"))
gctwt.MIN3 <- merge(gctwt.MIN3, gctwt.star.filtered, all.x = TRUE, by = c("fusionName"))
gctwt.MIN3 <- merge(gctwt.MIN3, gctwt.integrate.filtered, all.x = TRUE, by = c("fusionName"))
## gctwt.MIN3$fusionName <- droplevels(gctwt.MIN3$fusionName)

min2Count <- countFusions(gctwt.MIN2, caseID = "case.x")
min3Count <- countFusions(gctwt.MIN3, caseID = "case.x")
#min4Count <- countFusions(gctwt.MIN4, caseID = "case.x")


pdf("GCTwt/GCTwt_MIN2_MIN3_MIN4.pdf", width = 110/25.4, height = 110/25.4)
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
ggplot(gctwtCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in all four tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


write.table(gctwt.MIN2, file = "GCTwt/GCTwt_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwt.MIN3, file = "GCTwt/GCTwt_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctwt4M, file = "GCTwt/GCTwt_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("GCTwt/fusions.gctwt.rda")
