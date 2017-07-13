## working with salivary ACC

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

project <- "Pancreatic Sem."
outDir <- "PNC"

## cases to include
pncInc <-  unique(c(star$case[star$project == project],
                    map$case[map$project == project],
                    ee$case[ee$project == project],
                    integrate$case[integrate$project == project]))

## extract PNC unfiltered
pnc.star <- star[star$case %in% pncInc,]
pnc.map <- map[map$case %in% pncInc,]
pnc.ee <- ee[ee$case %in% pncInc,]
pnc.integrate <- integrate[integrate$case %in% pncInc,]
pnc.defuse <- defuse[defuse$case %in% pncInc,]

## extract PNC from individual filtered fusions
pnc.star.filtered <- star.filtered[star.filtered$project == project, ]
pnc.star.filtered <- pnc.star.filtered[!is.na(pnc.star.filtered$case),]
pnc.star.filtered <- pnc.star.filtered[pnc.star.filtered$case %in% pncInc, ]

pnc.map.filtered <- map.filtered[map.filtered$project == project, ]
pnc.map.filtered <- pnc.map.filtered[!is.na(pnc.map.filtered$case),]
pnc.map.filtered <- pnc.map.filtered[pnc.map.filtered$case %in% pncInc, ]

pnc.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
pnc.ee.filtered <- pnc.ee.filtered[!is.na(pnc.ee.filtered$case),]
pnc.ee.filtered <- pnc.ee.filtered[pnc.ee.filtered$case %in% pncInc, ]

pnc.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
pnc.integrate.filtered <- pnc.integrate.filtered[!is.na(pnc.integrate.filtered$case),]
pnc.integrate.filtered <- pnc.integrate.filtered[pnc.integrate.filtered$case %in% pncInc, ]

pnc.defuse.filtered <- defuse.filtered[defuse.filtered$project == project, ]
pnc.defuse.filtered <- pnc.defuse.filtered[!is.na(pnc.defuse.filtered$case),]
pnc.defuse.filtered <- pnc.defuse.filtered[pnc.defuse.filtered$case %in% pncInc, ]

## extract PNC from recurrent filtered fusions
pnc.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
pnc.star.recurrent <- pnc.star.recurrent[!is.na(pnc.star.recurrent$case),]
pnc.star.recurrent <- pnc.star.recurrent[pnc.star.recurrent$case %in% pncInc, ]

pnc.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
pnc.map.recurrent <- pnc.map.recurrent[!is.na(pnc.map.recurrent$case),]
pnc.map.recurrent <- pnc.map.recurrent[pnc.map.recurrent$case %in% pncInc, ]

pnc.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
pnc.ee.recurrent <- pnc.ee.recurrent[!is.na(pnc.ee.recurrent$case),]
pnc.ee.recurrent <- pnc.ee.recurrent[pnc.ee.recurrent$case %in% pncInc, ]

pnc.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
pnc.ifr.recurrent <- pnc.ifr.recurrent[!is.na(pnc.ifr.recurrent$case),]
pnc.ifr.recurrent <- pnc.ifr.recurrent[pnc.ifr.recurrent$case %in% pncInc, ]

pnc.def.recurrent <- def.recurrent[def.recurrent$project == project, ]
pnc.def.recurrent <- pnc.def.recurrent[!is.na(pnc.def.recurrent$case),]
pnc.def.recurrent <- pnc.def.recurrent[pnc.def.recurrent$case %in% pncInc, ]


## pair-wise interesections
pncSM <- merge(pnc.star, pnc.map, by = c("fusionName", "case"))
pncSI <- merge(pnc.star, pnc.integrate, by = c("fusionName", "case"))
pncSE <- merge(pnc.star, pnc.ee, by = c("fusionName", "case"))
pncSD <- merge(pnc.star, pnc.defuse, by = c("fusionName", "case"))
pncME <- merge(pnc.map, pnc.ee, by = c("fusionName", "case"))
pncMI <- merge(pnc.map, pnc.integrate, by = c("fusionName", "case"))
pncMD <- merge(pnc.map, pnc.defuse, by = c("fusionName", "case"))
pncEI <- merge(pnc.ee, pnc.integrate, by = c("fusionName", "case"))
pncDI <- merge(pnc.defuse, pnc.integrate, by = c("fusionName", "case"))
pncDE <- merge(pnc.defuse, pnc.ee, by = c("fusionName", "case"))

## three-way interesections
pncMSE <- merge(pncSM, pnc.ee, by = c("fusionName", "case"))
pncMSI <- merge(pncSM, pnc.integrate, by = c("fusionName", "case"))
pncMSD <- merge(pncSM, pnc.defuse, by = c("fusionName", "case"))
pncMEI <- merge(pncME, pnc.integrate, by = c("fusionName", "case"))
pncMED <- merge(pncME, pnc.defuse, by = c("fusionName", "case"))
pncSEI <- merge(pncSE, pnc.integrate, by = c("fusionName", "case"))
pncSED <- merge(pncSE, pnc.defuse, by = c("fusionName", "case"))
pncMID <- merge(pncMI, pnc.defuse, by = c("fusionName", "case"))
pncIED <- merge(pncEI, pnc.defuse, by = c("fusionName", "case"))

## all programs
#pnc4M <- merge(pncSM, pncEI, by = c("fusionName", "case"))

## write individual intersections
##write.table(pncSM, file = "PNC/PNC_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(pncSE, file = "PNC/PNC_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
##write.table(pncME, file = "PNC/PNC_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(pncMSE, file = "PNC/PNC_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(pncMSI, file = "PNC/PNC_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(pncMEI, file = "PNC/PNC_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(pncSEI, file = "PNC/PNC_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(pnc4M, file = "PNC/PNC_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(pnc.star.recurrent, file = "PNC/PNC_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.map.recurrent, file = "PNC/PNC_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.ee.recurrent, file = "PNC/PNC_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.ifr.recurrent, file = "PNC/PNC_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.def.recurrent, file = "PNC/PNC_recurrent_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write Unfiltered
write.table(pnc.star, file = "PNC/PNC_unfiltered_fusions_STAR-Fusion.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.map, file = "PNC/PNC_unfiltered_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.ee, file = "PNC/PNC_unfiltered_fusions_Ericscript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.int, file = "PNC/PNC_unfiltered_fusions_INTEGRATE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.def, file = "PNC/PNC_unfiltered_fusions_DeFuse.txt", sep = "\t", quote = FALSE, row.names = FALSE)

unfilt <- merge(pnc.star, pnc.map, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, pnc.ee, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, pnc.int, by = c("fusionName", "case"), all = TRUE)
unfilt <- merge(unfilt, pnc.def, by = c("fusionName", "case"), all = TRUE)

write.table(unfilt, file = "PNC/PNC_unfiltered_fusions_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## count individual filtered by programs
pncSCount <- countFusions(pnc.star.filtered, "fusionName", "case")
pncMCount <- countFusions(pnc.map.filtered, "fusionName", "case")
pncECount <- countFusions(pnc.ee.filtered, "fusionName", "case")
pncICount <- countFusions(pnc.integrate.filtered, "fusionName", "case")
## count recurrent filtered by programs
pncSfrCount <- countFusions(pnc.star.recurrent, "fusionName", "case")
pncMfrCount <- countFusions(pnc.map.recurrent, "fusionName", "case")
pncEfrCount <- countFusions(pnc.ee.recurrent, "fusionName", "case")
pncIfrCount <- countFusions(pnc.ifr.recurrent, "fusionName", "case")

## count intersections
pncSMCounts <- countFusions(pncSM, "fusionName", "case")
pncSECounts <- countFusions(pncSE, "fusionName", "case")
pncSICounts <- countFusions(pncSI, "fusionName", "case")
pncSDCounts <- countFusions(pncSD, "fusionName", "case")
pncMECounts <- countFusions(pncME, "fusionName", "case")
pncMICounts <- countFusions(pncMI, "fusionName", "case")
pncMDCounts <- countFusions(pncMD, "fusionName", "case")
pncEICounts <- countFusions(pncEI, "fusionName", "case")
pncDICounts <- countFusions(pncDI, "fusionName", "case")
pncDECounts <- countFusions(pncDE, "fusionName", "case")


##pncCounts <- countFusions(pnc4M, "fusionName", "case")



pdf("PNC/PNC_byProgram.pdf", width = 150/25.4, height = 297/25.4)
ggplot(pncSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()



pdf("PNC/PNC_Recurrent_byProgram.pdf", width = 150/25.4, height = 180/25.4)
ggplot(pncSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + 
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


pdf("PNC/PNC_predictedFusions.pdf", width = 105/25.4, height = 105/25.4)
ggplot(pncSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncSICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncSDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & DeFuse Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncMDCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & DeFuse Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncMICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncEICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("EricScript & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncDICounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & Integrate Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(pncDECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("DeFuse & EricScript Joint fusion Prediction") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()

##ggplot(pncCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
##    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") +
##    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
##    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
##                      labels = c("Not Present", "Predicted"))


pnc.MIN2 <- unique(do.call(rbind, lapply(list(pncSM, pncSI, pncSE, pncSD, pncME, pncMI, pncMD, pncEI, pncDI, pncDE), function(X) X[,c("fusionName", "case")])))
pnc.MIN2 <- merge(pnc.MIN2, pnc.ee, all.x = TRUE, by = c("fusionName"))
pnc.MIN2 <- merge(pnc.MIN2, pnc.map, all.x = TRUE, by = c("fusionName"))
pnc.MIN2 <- merge(pnc.MIN2, pnc.star, all.x = TRUE, by = c("fusionName"))
pnc.MIN2 <- merge(pnc.MIN2, pnc.integrate, all.x = TRUE, by = c("fusionName"))
pnc.MIN2 <- merge(pnc.MIN2, pnc.defuse, all.x = TRUE, by = c("fusionName"))
##pnc.MIN2$fusionName <- droplevels(pnc.MIN2$fusionName)

pnc.MIN3 <- unique(do.call(rbind, lapply(list(pncMSE, pncMSI, pncMSD, pncMEI, pncMED, pncSEI, pncSED, pncMID, pncIED), function(X) X[,c("fusionName", "case")])))
pnc.MIN3 <- merge(pnc.MIN3, pnc.ee, all.x = TRUE, by = c("fusionName"))
pnc.MIN3 <- merge(pnc.MIN3, pnc.map, all.x = TRUE, by = c("fusionName"))
pnc.MIN3 <- merge(pnc.MIN3, pnc.star, all.x = TRUE, by = c("fusionName"))
pnc.MIN3 <- merge(pnc.MIN3, pnc.integrate, all.x = TRUE, by = c("fusionName"))
pnc.MIN3 <- merge(pnc.MIN3, pnc.defuse, all.x = TRUE, by = c("fusionName"))
## pnc.MIN3$fusionName <- droplevels(pnc.MIN3$fusionName)

min2Count <- countFusions(pnc.MIN2, caseID = "case.x")
min3Count <- countFusions(pnc.MIN3, caseID = "case.x")
#min4Count <- countFusions(pnc.MIN4, caseID = "case.x")


pdf("PNC/PNC_MIN2_MIN3.pdf", width = 110/25.4, height = 297/25.4)
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
dev.off()


write.table(pnc.MIN2, file = "PNC/PNC_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pnc.MIN3, file = "PNC/PNC_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##write.table(pnc4M, file = "PNC/PNC_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("PNC/fusions.pnc.rda")

load("PNC/fusions.pnc.rda")
## java -jar /opt/src/oncofuse-1.1.1/Oncofuse.jar -a hg19 PNC/pnc_min2_oncofuse.input.txt 'coord' EPI PNC/pnc_min2_oncofuse.output.txt

Onc <- do.call(rbind, lapply(list.files(pattern = "merged.*oncofuse.output.txt"), read.delim, header = TRUE))
Onc$fusionName <- paste(Onc$X5_FPG_GENE_NAME, Onc$X3_FPG_GENE_NAME, "--")

pncOnc <- rbind(pnc.ee.recurrent[, c("program", "case", "fusionName")], pnc.map.recurrent[, c("program", "case", "fusionName")], pnc.star.recurrent[, c("program", "case", "fusionName")], pnc.ifr.recurrent[, c("program", "case", "fusionName")])

pncOnc <- merge(pncOnc, Onc, by = "fusionName")

pdf("PNC/PNC_Fusions_DrivProbability.pdf", width = 148/25.4, height = 210/25.4, useDingbats = FALSE)
ggplot(pncOnc, aes(case.x, fusionName)) + geom_tile(aes(fill = DRIVER_PROB), width = 0.95, height = 0.95) + ylab("Predicted Fusion") + xlab(paste(project, "Case"))
dev.off()


##recFtx <- c("HMGA2--SMUG1", "CCDC169--XRCC2", "CS--", "HELB--", "RP11−96H19.1−−CCDC38", "HMGA2−−RP11−834C11.5")

##recOnc <- Onc[as.character(Onc$fusionName) %in% recFtx, ]
