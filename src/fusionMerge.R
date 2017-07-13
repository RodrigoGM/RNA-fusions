## working with salivary ACC

library(ggplot2)
library(reshape2)
source("myLib.R")

load("fusions.data.rda")

mmSM <- merge(star.filtered, map.filtered, by = c("fusionName", "case"))
mmSI <- merge(star.filtered, integrate, by = c("fusionName", "case"))
mmSE <- merge(star.filtered, ee.filtered, by = c("fusionName", "case"))
mmME <- merge(map.filtered, ee.filtered, by = c("fusionName", "case"))
mmMI <- merge(map.filtered, integrate, by = c("fusionName", "case"))
mmEI <- merge(ee.filtered, integrate, by = c("fusionName", "case"))

mmMSE <- merge(mmSM, ee.filtered, by = c("fusionName", "case"))
mmMSI <- merge(mmSM, integrate, by = c("fusionName", "case"))
mmMEI <- merge(mmME, integrate, by = c("fusionName", "case"))
mmSEI <- merge(mmSE, integrate, by = c("fusionName", "case"))

mm4M <- merge(mmSM, mmEI, by = c("fusionName", "case"))


write.table(mmSM, file = "MM_fusions_STAR_MAPSplice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mmSE, file = "MM_fusions_STAR_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mmME, file = "MM_fusions_MAPSplice_EricScript.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mmMSE, file = "MM_fusions_MSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mm4M, file = "MM_fusions_4M.txt", sep = "\t", quote = FALSE, row.names = FALSE)


mmSMCounts <- melt(table(mmSM$fusionName, mmSM$case))
mmSMCounts$value[mmSMCounts$value == 0] <- NA
mmSMCounts$qual <- mmSMCounts$value
mmSMCounts$qual[mmSMCounts$qual == 0] <- NA
mmSMCounts$qual[!is.na(mmSMCounts$qual)] <- 1

mmSECounts <- melt(table(mmSE$fusionName, mmSE$case))
mmSECounts$value[mmSECounts$value == 0] <- NA
mmSECounts$qual <- mmSECounts$value
mmSECounts$qual[mmSECounts$qual == 0] <- NA
mmSECounts$qual[!is.na(mmSECounts$qual)] <- 1

mmMECounts <- melt(table(mmME$fusionName, mmME$case))
mmMECounts$value[mmMECounts$value == 0] <- NA
mmMECounts$qual <- mmMECounts$value
mmMECounts$qual[mmMECounts$qual == 0] <- NA
mmMECounts$qual[!is.na(mmMECounts$qual)] <- 1

mmCounts <- melt(table(mmMSE$fusionName, mmMSE$case))
mmCounts$value[mmCounts$value == 0] <- NA
mmCounts$qual <- mmCounts$value
mmCounts$qual[mmCounts$qual == 0] <- NA
mmCounts$qual[!is.na(mmCounts$qual)] <- 1


MIN2 <- unique(do.call(rbind, lapply(list(mmSM, mmSE, mmSI, mmME, mmMI, mmEI), function(X) X[,c("fusionName", "case")])))
MIN2 <- merge(MIN2, ee.filtered, all.x = TRUE, by = c("fusionName"))
MIN2 <- merge(MIN2, map.filtered, all.x = TRUE, by = c("fusionName"))
MIN2 <- merge(MIN2, star.filtered, all.x = TRUE, by = c("fusionName"))
MIN2 <- merge(MIN2, integrate, all.x = TRUE, by = c("fusionName"))

min2Count <- melt(table(MIN2$fusionName, MIN2$case.x))
min2Count$value[min2Count$value == 0] <- NA
min2Count$qual <- min2Count$value
min2Count$qual[!is.na(min2Count$qual)] <- "Predicted"
min2Count$qual <- setKnown(min2Count, caseID = "Var.2", fusionName = "Var.1")
min2Count$project <- fuseProject(min2Count, caseID = "Var.2")

MIN3 <- unique(do.call(rbind, lapply(list(mmMEI, mmMSE, mmMSI, mmSEI), function(X) X[,c("fusionName", "case")])))
MIN3 <- merge(MIN3, ee.filtered, all.x = TRUE, by = c("fusionName"))
MIN3 <- merge(MIN3, map.filtered, all.x = TRUE, by = c("fusionName"))
MIN3 <- merge(MIN3, star.filtered, all.x = TRUE, by = c("fusionName"))
MIN3 <- merge(MIN3, integrate, all.x = TRUE, by = c("fusionName"))

min3Count <- melt(table(MIN3$fusionName, MIN3$case.x))
min3Count$value[min3Count$value == 0] <- NA
min3Count$qual <- min3Count$value
min3Count$qual[!is.na(min3Count$qual)] <- "Predicted"
min3Count$qual <- setKnown(min3Count, caseID = "Var.2", fusionName = "Var.1")
min3Count$project <- fuseProject(min3Count, caseID = "Var.2")


min4Count <- melt(table(mm4M$fusionName, mm4M$case))
min4Count$value[min4Count$value == 0] <- NA
min4Count$qual <- min4Count$value
min4Count$qual[!is.na(min4Count$qual)] <- "Predicted"
min4Count$qual <- setKnown(min4Count, caseID = "Var.2", fusionName = "Var.1")
min4Count$project <- fuseProject(min4Count, caseID = "Var.2")

pdf("MM_MIN2_MIN3_MIN4.pdf", width = 297/25.4, height = 210/25.4)
ggplot(min2Count, aes(Var.2, Var.1)) + geom_tile(aes(fill = (qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90), legend.position = "none", axis.text.y = element_text(size = 4)) +
    xlab("Case") + ylab("Fusion predictions comon in at least two tools") +
    facet_grid(. ~ project, scales = "free")

ggplot(min3Count, aes(Var.2, Var.1)) + geom_tile(aes(fill = (qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90), legend.position = "none", axis.text.y = element_text(size = 10)) +
    xlab("Case") + ylab("Fusion predictions comon in at least three tools") +
    facet_grid(. ~ project, scales = "free")

ggplot(min4Count, aes(Var.2, Var.1)) + geom_tile(aes(fill = (qual)), width = 0.95, height = 0.95) +
    theme(axis.text.x  = element_text(angle = 90), legend.position = "none", axis.text.y = element_text(size = 10)) +
    xlab("Salivary Gland ACC Case") + ylab("Fusion predictions comon in all four tools") +
    facet_grid(. ~ project, scales = "free")
dev.off()

write.table(MIN2, file = "MM_Unique_Fusions_in_2_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(MIN3, file = "MM_Unique_Fusions_in_3_tools.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## fusionPlots.R is for plotting, here just be aware it's not what you may want
## ggplot(mmCounts, aes(Var.2, Var.1)) + geom_tile(aes(fill = factor(value)))
## pdf("PredictedFusions.pdf", width = 297/25.4, height = 210/25.4)
## ggplot(mmSMCounts, aes(Var.2, Var.1)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) + xlab("Case") + ylab("STAR-Fusion & MAPSplice Joint Fusion Prediction") + theme(legend.position = "none")
## ggplot(mmSECounts[!is.na(mmSECounts$value), ], aes(Var.2, Var.1)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) + xlab("Case") + ylab("STAR-Fusion & EricScript Joint fusion Prediction") + theme(legend.position = "none", axis.text.y = element_text(size = 4))
## ggplot(mmMECounts[!is.na(mmMECounts$value), ], aes(Var.2, Var.1)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) + xlab("Case") + ylab("MAPSplice & EricScript Joint fusion Prediction") + theme(legend.position = "none")
## ggplot(mmCounts, aes(Var.2, Var.1)) + geom_tile(aes(fill = qual), width = 0.95, height = 0.95) + xlab("Case") + xlab("Case") + ylab("All 3 Joint fusion Prediction") + theme(legend.position = "none")
## dev.off()


