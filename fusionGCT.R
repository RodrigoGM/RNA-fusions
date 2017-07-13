## working with GCT

library(ggplot2)
library(reshape2)
library(naturalsort)
source("myLib.R")

# load("fusions.data.rda")
# load("GCT/fusions.gct.rda")

project <- "Granular Cell Tumor"
outDir <- "GCT"

## cases to include
gctInc <-  unique(c(star$case[star$project == project], map$case[map$project == project], ee$case[ee$project == project], integrate$case[integrate$project == project]))


## extract GCT from individual filtered fusions
gct.star.filtered <- star.filtered[star.filtered$project == project, ]
gct.star.filtered <- gct.star.filtered[!is.na(gct.star.filtered$case),]
gct.star.filtered <- gct.star.filtered[gct.star.filtered$case %in% gctInc, ]

gct.map.filtered <- map.filtered[map.filtered$project == project, ]
gct.map.filtered <- gct.map.filtered[!is.na(gct.map.filtered$case),]
gct.map.filtered <- gct.map.filtered[gct.map.filtered$case %in% gctInc, ]

gct.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
gct.ee.filtered <- gct.ee.filtered[!is.na(gct.ee.filtered$case),]
gct.ee.filtered <- gct.ee.filtered[gct.ee.filtered$case %in% gctInc, ]

gct.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
gct.integrate.filtered <- gct.integrate.filtered[!is.na(gct.integrate.filtered$case),]
gct.integrate.filtered <- gct.integrate.filtered[gct.integrate.filtered$case %in% gctInc, ]

## extract GCT from recurrent filtered fusions
gct.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
gct.star.recurrent <- gct.star.recurrent[!is.na(gct.star.recurrent$case),]
gct.star.recurrent <- gct.star.recurrent[gct.star.recurrent$case %in% gctInc, ]

gct.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
gct.map.recurrent <- gct.map.recurrent[!is.na(gct.map.recurrent$case),]
gct.map.recurrent <- gct.map.recurrent[gct.map.recurrent$case %in% gctInc, ]

gct.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
gct.ee.recurrent <- gct.ee.recurrent[!is.na(gct.ee.recurrent$case),]
gct.ee.recurrent <- gct.ee.recurrent[gct.ee.recurrent$case %in% gctInc, ]

gct.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
gct.ifr.recurrent <- gct.ifr.recurrent[!is.na(gct.ifr.recurrent$case),]
gct.ifr.recurrent <- gct.ifr.recurrent[gct.ifr.recurrent$case %in% gctInc, ]

## pair-wise interesections
gctSM <- merge(gct.star.filtered, gct.map.filtered, by = c("fusionName", "case"))
gctSI <- merge(gct.star.filtered, gct.integrate.filtered, by = c("fusionName", "case"))
gctSE <- merge(gct.star.filtered, gct.ee.filtered, by = c("fusionName", "case"))
gctME <- merge(gct.map.filtered, gct.ee.filtered, by = c("fusionName", "case"))
gctMI <- merge(gct.map.filtered, gct.integrate.filtered, by = c("fusionName", "case"))
gctEI <- merge(gct.ee.filtered, gct.integrate.filtered, by = c("fusionName", "case"))

## three-way interesections
gctMSE <- merge(gctSM, gct.ee.filtered, by = c("fusionName", "case"))
gctMSI <- merge(gctSM, gct.integrate.filtered, by = c("fusionName", "case"))
gctMEI <- merge(gctME, gct.integrate.filtered, by = c("fusionName", "case"))
gctSEI <- merge(gctSE, gct.integrate.filtered, by = c("fusionName", "case"))

## all programs
gct4M <- merge(gctSM, gctEI, by = c("fusionName", "case"))

## write individual intersections
write.table(gctSM, file = "GCT/GCT_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctSE, file = "GCT/GCT_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctME, file = "GCT/GCT_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctMSE, file = "GCT/GCT_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctMSI, file = "GCT/GCT_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctMEI, file = "GCT/GCT_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gctSEI, file = "GCT/GCT_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gct4M, file = "GCT/GCT_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(gct.star.recurrent, file = "GCT/GCT_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gct.map.recurrent, file = "GCT/GCT_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gct.ee.recurrent, file = "GCT/GCT_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gct.ifr.recurrent, file = "GCT/GCT_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## count individual filtered by programs
gctSCount <- countFusions(gct.star.filtered, "fusionName", "case")
gctMCount <- countFusions(gct.map.filtered, "fusionName", "case")
gctECount <- countFusions(gct.ee.filtered, "fusionName", "case")
gctICount <- countFusions(gct.integrate.filtered, "fusionName", "case")

## count recurrent filtered by programs
gctSfrCount <- countFusions(gct.star.recurrent, "fusionName", "case")
gctMfrCount <- countFusions(gct.map.recurrent, "fusionName", "case")
gctEfrCount <- countFusions(gct.ee.recurrent, "fusionName", "case")
gctIfrCount <- countFusions(gct.ifr.recurrent, "fusionName", "case")

## count intersections
gctSMCounts <- countFusions(gctSM, "fusionName", "case")
gctSECounts <- countFusions(gctSE, "fusionName", "case")
gctSICounts <- countFusions(gctSI, "fusionName", "case")
gctMECounts <- countFusions(gctME, "fusionName", "case")
gctMICounts <- countFusions(gctMI, "fusionName", "case")
gctEICounts <- countFusions(gctEI, "fusionName", "case")
gctCounts <- countFusions(gct4M, "fusionName", "case")



pdf("GCT/GCT_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(gctSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()



pdf("GCT/GCT_Recurrent_byProgram.pdf", width = 150/25.4, height = 200/25.4)
ggplot(gctSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


pdf("GCT/GCT_predictedFusions.pdf", width = 105/25.4, height = 297/25.4)
ggplot(gctSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


gct.MIN2 <- unique(do.call(rbind, lapply(list(gctEI, gctME, gctMI, gctSE, gctSI, gctSM), function(X) X[,c("fusionName", "case")])))
gct.MIN2 <- merge(gct.MIN2, gct.ee.filtered, all.x = TRUE, by = c("fusionName"))
gct.MIN2 <- merge(gct.MIN2, gct.map.filtered, all.x = TRUE, by = c("fusionName"))
gct.MIN2 <- merge(gct.MIN2, gct.star.filtered, all.x = TRUE, by = c("fusionName"))
gct.MIN2 <- merge(gct.MIN2, gct.integrate.filtered, all.x = TRUE, by = c("fusionName"))
##gct.MIN2$fusionName <- droplevels(gct.MIN2$fusionName)

gct.MIN3 <- unique(do.call(rbind, lapply(list(gctMEI, gctMSE, gctMSI, gctSEI), function(X) X[,c("fusionName", "case")])))
gct.MIN3 <- merge(gct.MIN3, gct.ee.filtered, all.x = TRUE, by = c("fusionName"))
gct.MIN3 <- merge(gct.MIN3, gct.map.filtered, all.x = TRUE, by = c("fusionName"))
gct.MIN3 <- merge(gct.MIN3, gct.star.filtered, all.x = TRUE, by = c("fusionName"))
gct.MIN3 <- merge(gct.MIN3, gct.integrate.filtered, all.x = TRUE, by = c("fusionName"))
## gct.MIN3$fusionName <- droplevels(gct.MIN3$fusionName)

min2Count <- countFusions(gct.MIN2, caseID = "case.x")
min3Count <- countFusions(gct.MIN3, caseID = "case.x")
#min4Count <- countFusions(gct.MIN4, caseID = "case.x")


pdf("GCT/GCT_MIN2_MIN3_MIN4.pdf", width = 110/25.4, height = 110/25.4)
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
ggplot(gctCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 8)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in all four tools") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


write.table(gct.MIN2, file = "GCT/GCT_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gct.MIN3, file = "GCT/GCT_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gct4M, file = "GCT/GCT_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("GCT/fusions.gct.rda")


####

## FILTERED PLOTS

ee.blast <- read.table("GCT/EricSCript.nhom.FILTER.txt")
ee.fusionName <- unique(paste(ee.blast$V1, ee.blast$V2, sep = "--"))

int.blast <- read.table("GCT/Integrate.FILTER.txt")
int.fusionName <- unique(paste(int.blast$V1, int.blast$V2, sep = "--"))

star.fusionName <- readLines("GCT/STAR.maunualFilter.txt")


## extras
gctECount.nh <- countFusions(gct.ee.filtered[gct.ee.filtered$homology == "",], "fusionName", "case")
gctECount.nh$case <- naturalfactor(gctECount.nh$case)
gctECount.nh <- gctECount.nh[! gctECount.nh$fusionName %in% ee.fusionName , ]

gctSCount$case <- naturalfactor(gctSCount$case)

gctICount$case <- naturalfactor(gctICount$case)
gctICount <- gctICount[! gctICount$fusionName %in% int.fusionName , ]

write.table(gctSCount$fusionName, file = "STAR-Fusion.fusions.txt", quote = FALSE, row.names = FALSE)
write.table(gctECount.nh$fusionName, file = "EricSCript.nhom.fusions.txt", quote = FALSE, row.names = FALSE)
write.table(gctICount$fusionName, file = "Integrate.fusions.txt", quote = FALSE, row.names = FALSE)



pdf("GCT/GCT_fusion_suppl.pdf", width = 148/25.4, height = 298/25.4)
ggplot(gctSCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6.9)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion Predictions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctECount.nh, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript Predictions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(gctICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6.9)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript Predictions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()



pdf("GCT/STAR-fusion_ManualFilter.pdf", width = 148/25.4, height = 210/25.4)
ggplot(gctSCount[! as.character(gctSCount$fusionName) %in% as.character(star.fusionName), ], aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6.9)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion Predictions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()
