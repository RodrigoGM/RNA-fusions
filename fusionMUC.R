## working with salivary MUCINOUS

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

load("MUC/fusions.muc.rda")
project <- "Mucinus"
outDir <- "MUC"

## cases to include
mucInc <-  unique(c(star$case[star$project == project],
                    map$case[map$project == project],
                    ee$case[ee$project == project],
                    integrate$case[integrate$project == project]))

mucType <- read.delim("MUC/mucinus_types.txt", colClasses = c("character", "character", "character"))
mucType$tumor.type <- factor(mucType$tumor.type, levels = c("Mucinous A", "DCIS", "Mucinous B", "Mucinous M", "IDC"))


## extract MUC from individual filtered fusions
muc.star.filtered <- star.filtered[star.filtered$project == project, ]
muc.star.filtered <- muc.star.filtered[!is.na(muc.star.filtered$case),]
muc.star.filtered <- muc.star.filtered[muc.star.filtered$case %in% mucInc, ]
muc.star.filtered <- merge(muc.star.filtered, mucType, by = "case")
muc.star.filtered$case <- factor(muc.star.filtered$case, levels = c(mucType[mucType$tumor.type == "Mucinous A", "case"], mucType[mucType$tumor.type == "Mucinous B", "case"], mucType[mucType$tumor.type == "Mucinous M", "case"], mucType[mucType$tumor.type == "IDC", "case"]))

muc.map.filtered <- map.filtered[map.filtered$project == project, ]
muc.map.filtered <- muc.map.filtered[!is.na(muc.map.filtered$case),]
muc.map.filtered <- muc.map.filtered[muc.map.filtered$case %in% mucInc, ]
muc.map.filtered <- merge(muc.map.filtered, mucType, by = "case")
muc.map.filtered$case <- factor(muc.map.filtered$case, levels = c(mucType[mucType$tumor.type == "Mucinous A", "case"], mucType[mucType$tumor.type == "Mucinous B", "case"], mucType[mucType$tumor.type == "Mucinous M", "case"], mucType[mucType$tumor.type == "IDC", "case"]))

muc.ee.filtered <- ee.filtered[ee.filtered$project == project, ]
muc.ee.filtered <- muc.ee.filtered[!is.na(muc.ee.filtered$case),]
muc.ee.filtered <- muc.ee.filtered[muc.ee.filtered$case %in% mucInc, ]
muc.ee.filtered <- merge(muc.ee.filtered, mucType, by = "case")
muc.ee.filtered$case <- factor(muc.ee.filtered$case, levels = c(mucType[mucType$tumor.type == "Mucinous A", "case"], mucType[mucType$tumor.type == "Mucinous B", "case"], mucType[mucType$tumor.type == "Mucinous M", "case"], mucType[mucType$tumor.type == "IDC", "case"]))
muc.ee.filtered <- merge(muc.ee.filtered, mucType, by = "case")

muc.integrate.filtered <- integrate.filtered[integrate.filtered$project == project, ]
muc.integrate.filtered <- muc.integrate.filtered[!is.na(muc.integrate.filtered$case),]
muc.integrate.filtered <- muc.integrate.filtered[muc.integrate.filtered$case %in% mucInc, ]
muc.integrate.filtered <- merge(muc.integrate.filtered, mucType, by = "case")
muc.integrate.filtered$case <- factor(muc.integrate.filtered$case, levels = c(mucType[mucType$tumor.type == "Mucinous A", "case"], mucType[mucType$tumor.type == "Mucinous B", "case"]))

## extract MUC from recurrent filtered fusions
muc.star.recurrent <- star.recurrent[star.recurrent$project == project, ]
muc.star.recurrent <- muc.star.recurrent[!is.na(muc.star.recurrent$case),]
muc.star.recurrent <- muc.star.recurrent[muc.star.recurrent$case %in% mucInc, ]
muc.star.recurrent <- merge(muc.star.recurrent, mucType, by = "case")
muc.star.recurrent$case <- factor(muc.star.recurrent$case, levels = c(mucType[mucType$tumor.type == "Mucinous A", "case"], mucType[mucType$tumor.type == "Mucinous B", "case"]))

muc.map.recurrent <- map.recurrent[map.recurrent$project == project, ]
muc.map.recurrent <- muc.map.recurrent[!is.na(muc.map.recurrent$case),]
muc.map.recurrent <- muc.map.recurrent[muc.map.recurrent$case %in% mucInc, ]
muc.map.recurrent <- merge(muc.map.recurrent, mucType, by = "case")
muc.map.recurrent$case <- factor(muc.map.recurrent$case, levels = c(mucType[mucType$tumor.type == "Mucinous A", "case"], mucType[mucType$tumor.type == "Mucinous B", "case"], mucType[mucType$tumor.type == "Mucinous M", "case"], mucType[mucType$tumor.type == "IDC", "case"]))


muc.ee.recurrent <- ee.recurrent[ee.recurrent$project == project, ]
muc.ee.recurrent <- muc.ee.recurrent[!is.na(muc.ee.recurrent$case),]
muc.ee.recurrent <- muc.ee.recurrent[muc.ee.recurrent$case %in% mucInc, ]
muc.ee.recurrent <- merge(muc.ee.recurrent, mucType, by = "case")
muc.ee.recurrent$case <- factor(muc.ee.recurrent$case, levels = c(mucType[mucType$tumor.type == "Mucinous A", "case"], mucType[mucType$tumor.type == "Mucinous B", "case"], mucType[mucType$tumor.type == "Mucinous M", "case"], mucType[mucType$tumor.type == "IDC", "case"]))


muc.ifr.recurrent <- ifr.recurrent[ifr.recurrent$project == project, ]
muc.ifr.recurrent <- muc.ifr.recurrent[!is.na(muc.ifr.recurrent$case),]
muc.ifr.recurrent <- muc.ifr.recurrent[muc.ifr.recurrent$case %in% mucInc, ]
muc.ifr.recurrent <- merge(muc.ee.recurrent, mucType, by = "case")
muc.ifr.recurrent$case <- factor(muc.ifr.recurrent$case, levels = c(mucType[mucType$tumor.type == "Mucinous A", "case"], mucType[mucType$tumor.type == "Mucinous B", "case"], mucType[mucType$tumor.type == "Mucinous M", "case"], mucType[mucType$tumor.type == "IDC", "case"]))


## pair-wise interesections
mucSM <- merge(muc.star.filtered, muc.map.filtered, by = c("fusionName", "case"))
mucSI <- merge(muc.star.filtered, muc.integrate.filtered, by = c("fusionName", "case"))
mucSE <- merge(muc.star.filtered, muc.ee.filtered, by = c("fusionName", "case"))
mucME <- merge(muc.map.filtered, muc.ee.filtered, by = c("fusionName", "case"))
mucMI <- merge(muc.map.filtered, muc.integrate.filtered, by = c("fusionName", "case"))
mucEI <- merge(muc.ee.filtered, muc.integrate.filtered, by = c("fusionName", "case"))

## Three-way interesections
mucMSE <- merge(mucSM, muc.ee.filtered, by = c("fusionName", "case"))
mucMSI <- merge(mucSM, muc.integrate.filtered, by = c("fusionName", "case"))
mucMEI <- merge(mucME, muc.integrate.filtered, by = c("fusionName", "case"))
mucSEI <- merge(mucSE, muc.integrate.filtered, by = c("fusionName", "case"))

## all programs
muc4M <- merge(mucSM, mucEI, by = c("fusionName", "case"))

## write individual intersections
write.table(mucSM, file = "MUC/MUC_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mucSE, file = "MUC/MUC_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mucME, file = "MUC/MUC_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mucMSE, file = "MUC/MUC_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mucMSI, file = "MUC/MUC_fusions_MSI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mucMEI, file = "MUC/MUC_fusions_MEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mucSEI, file = "MUC/MUC_fusions_SEI.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(muc4M, file = "MUC/MUC_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## write recurrent per program
write.table(muc.star.recurrent, file = "MUC/MUC_recurrent_fusions_STAR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(muc.map.recurrent, file = "MUC/MUC_recurrent_fusions_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(muc.ee.recurrent, file = "MUC/MUC_recurrent_fusions_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(muc.ifr.recurrent, file = "MUC/MUC_recurrent_fusions_Integrate.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## count individual filtered by programs
mucSCount <- countFusions(muc.star.filtered, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucMCount <- countFusions(muc.map.filtered, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucECount <- countFusions(muc.ee.filtered, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucICount <- countFusions(muc.integrate.filtered, "fusionName", "case", merge.type = TRUE, type.df = mucType)
## count recurrent filtered by programs
mucSfrCount <- countFusions(muc.star.recurrent, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucMfrCount <- countFusions(muc.map.recurrent, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucEfrCount <- countFusions(muc.ee.recurrent, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucIfrCount <- countFusions(muc.ifr.recurrent, "fusionName", "case", merge.type = TRUE, type.df = mucType)

## count intersections
mucSMCounts <- countFusions(mucSM, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucSECounts <- countFusions(mucSE, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucSICounts <- countFusions(mucSI, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucMECounts <- countFusions(mucME, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucMICounts <- countFusions(mucMI, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucEICounts <- countFusions(mucEI, "fusionName", "case", merge.type = TRUE, type.df = mucType)
mucCounts <- countFusions(muc4M, "fusionName", "case", merge.type = TRUE, type.df = mucType)



pdf("MUC/MUC_byProgram.pdf", width = 210/25.4, height = 298/25.4)
ggplot(mucSCount, aes(case, fusionName)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 3)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucMCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 3)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucECount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 3)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucICount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 3)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()



pdf("MUC/MUC_Recurrent_byProgram.pdf", width = 210/25.4, height = 298/25.4)
ggplot(mucSfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion predicted") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucMfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("MAPSplice predicted") +  facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucEfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("EricScript predicted") +  facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucIfrCount, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 4)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Integrate predicted") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
##    scale_fill_manual(values = c("#F8766D", "#CFCFCF"), name = "Fusion", breaks = c("Not Present", "Predicted"),
##                      labels = c("Predicted", "Not Present"))
dev.off()


pdf("MUC/MUC_predictedFusions.pdf", width = 210/25.4, height = 298/25.4)
ggplot(mucSMCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") + facet_grid( ~ tumor.type, scales = "free_x") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucSECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("STAR-Fusion & EricScript Joint fusion Prediction") + facet_grid( ~ tumor.type, scales = "free_x") +
    theme(legend.position = "none", axis.text.y = element_text(size = 4), axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucMECounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("MAPSplice & EricScript Joint fusion Prediction") + facet_grid( ~ tumor.type, scales = "free_x") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mucCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    xlab(paste(project, "Case")) + ylab("All 3 Joint fusion Prediction") + facet_grid( ~ tumor.type, scales = "free_x") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6)) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()


muc.MIN2 <- unique(do.call(rbind, lapply(list(mucEI, mucME, mucMI, mucSE, mucSI, mucSM), function(X) X[,c("fusionName", "case")])))
muc.MIN2 <- merge(muc.MIN2, muc.ee.filtered, all.x = TRUE, by = c("fusionName"))
muc.MIN2 <- merge(muc.MIN2, muc.map.filtered, all.x = TRUE, by = c("fusionName"))
muc.MIN2 <- merge(muc.MIN2, muc.star.filtered, all.x = TRUE, by = c("fusionName"))
muc.MIN2 <- merge(muc.MIN2, muc.integrate.filtered, all.x = TRUE, by = c("fusionName"))
##muc.MIN2$fusionName <- droplevels(muc.MIN2$fusionName)

muc.MIN3 <- unique(do.call(rbind, lapply(list(mucMEI, mucMSE, mucMSI, mucSEI), function(X) X[,c("fusionName", "case")])))
muc.MIN3 <- merge(muc.MIN3, muc.ee.filtered, all.x = TRUE, by = c("fusionName"))
muc.MIN3 <- merge(muc.MIN3, muc.map.filtered, all.x = TRUE, by = c("fusionName"))
muc.MIN3 <- merge(muc.MIN3, muc.star.filtered, all.x = TRUE, by = c("fusionName"))
muc.MIN3 <- merge(muc.MIN3, muc.integrate.filtered, all.x = TRUE, by = c("fusionName"))
## muc.MIN3$fusionName <- droplevels(muc.MIN3$fusionName)

min2Count <- countFusions(muc.MIN2, caseID = "case.x")
min2Count <- merge(min2Count, mucType, by.x = "case.x", by.y = "case", all.x = TRUE)
min3Count <- countFusions(muc.MIN3, caseID = "case.x")
min3Count <- merge(min3Count, mucType, by.x = "case.x", by.y = "case", all.x = TRUE)

#min4Count <- countFusions(muc.MIN4, caseID = "case.x")


pdf("MUC/MUC_MIN2_MIN3_MIN4.pdf", width = 210/25.4, height = 298/25.4)
ggplot(min2Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6), legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least two tools") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(min3Count, aes(case.x, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in at least three tools") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

ggplot(mucCounts, aes(case, fusionName)) + geom_tile(aes(fill = factor(qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6)) +   ## , legend.position = "none") +
    xlab(paste(project, "Case")) + ylab("Fusion predictions comon in all four tools") + facet_grid( ~ tumor.type, scales = "free_x") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))

dev.off()


write.table(muc.MIN2, file = "MUC/MUC_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(muc.MIN3, file = "MUC/MUC_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(muc4M, file = "MUC/MUC_Unique_Fusions_in_4_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("MUC/fusions.muc.rda")

load("MUC/fusions.muc.rda")
## java -jar /opt/src/oncofuse-1.1.1/Oncofuse.jar -a hg19 MUC/muc_min2_oncofuse.input.txt 'coord' EPI MUC/muc_min2_oncofuse.output.txt

Onc <- read.delim("all.oncofuse.annotated.fusions.txt")

mucOnc <- rbind(muc.ee.recurrent[, c("program", "case", "fusionName")], muc.map.recurrent[, c("program", "case", "fusionName")], muc.star.recurrent[, c("program", "case", "fusionName")], muc.ifr.recurrent[, c("program", "case", "fusionName")])

mucOnc <- merge(mucOnc, Onc, by = "fusionName")
mucOnc <- merge(mucOnc, mucType, by = "case")

pdf("MUC/MUC_Fusions_DrivProbability.pdf", width = 210/25.4, height = 298/25.4, useDingbats = FALSE)
ggplot(mucOnc, aes(case, fusionName)) + geom_tile(aes(fill = DRIVER_PROB), width = 0.95, height = 0.95) + ylab("Predicted Fusion") + xlab("Liposarcoma Case") + facet_grid( ~ tumor.type) + theme(axis.text.x  = element_text(angle = 90, hjust = 1))
dev.off()

## recFtx <- c("HMGA2--SMUG1", "CCDC169--XRCC2", "CS--", "HELB--", "RP11−96H19.1−−CCDC38", "HMGA2−−RP11−834C11.5")
## recOnc <- Onc[as.character(Onc$fusionName) %in% recFtx, ]
