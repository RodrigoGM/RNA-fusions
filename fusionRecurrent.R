## takes data to identify recurrent fusions per program

## As if any are recurrent prior to intersecting
## STAR-Fusion -- All fusions -- Are there any recurrent ?
## EricScript  -- All fusions -- Are there any recurrent ?
## MAPSplice   -- All fusions -- Are there any recurrent ?
## Integrate   -- All fusions -- Are there any recurrent ?
## deFuse      -- All fusions -- Are there any recurrent ?

library(reshape2)
library(ggplot2)
source("myLib.R")

## load data  -- make sure that the myLib.R is the updtodate otheriwse functions get re-written
load("fusions.data.rda")

## count fusions by case
sfr <- getRecurrent(star.filtered, caseID = "case", fusionName = "fusionName")
mfr <- getRecurrent(map.filtered, caseID = "case", fusionName = "fusionName")
efr <- getRecurrent(ee.filtered, caseID = "case", fusionName = "fusionName")
ifr <- getRecurrent(integrate.filtered, caseID = "case", fusionName = "fusionName")
dfr <- getRecurrent(defuse.filtered, caseID = "case", fusionName = "fusionName")

## keep only those that appear in more than one case
star.recurrent <- star.filtered[star.filtered$fusionName %in% rownames(sfr), ]
map.recurrent <- map.filtered[map.filtered$fusionName %in% rownames(mfr), ]
ee.recurrent <- ee.filtered[ee.filtered$fusionName %in% rownames(efr), ]
ifr.recurrent <- integrate.filtered[integrate.filtered$fusionName %in% rownames(ifr), ]
def.recurrent <- defuse.filtered[defuse.filtered$fusionName %in% rownames(ifr), ]

## remove cross projecct fusions from the recurrent fusions
sfpr <- getCrossProject(star.recurrent)
star.recurrent <- star.recurrent[! star.recurrent$fusionName %in% rownames(sfpr), ]

mfpr <- getCrossProject(map.recurrent)
map.recurrent <- map.recurrent[! map.recurrent$fusionName %in% rownames(mfpr), ]

efpr <- getCrossProject(ee.recurrent)
ee.recurrent <- ee.recurrent[! ee.recurrent$fusionName %in% rownames(efpr), ]

ifpr <- getCrossProject(ifr.recurrent)
ifr.recurrent <- ifr.recurrent[! ifr.recurrent$fusionName %in% rownames(ifpr), ]

dfpr <- getCrossProject(def.recurrent)
def.recurrent <- def.recurrent[! def.recurrent$fusionName %in% rownames(dfpr), ]


## count recurrent fusions
sfrCount <- countFusions(star.recurrent)
sfrCount$project <- fuseProject(sfrCount)

mfrCount <- countFusions(map.recurrent)
mfrCount$project <- fuseProject(mfrCount)

efrCount <- countFusions(ee.recurrent)
efrCount$project <- fuseProject(efrCount)

ifrCount <- countFusions(ifr.recurrent)
ifrCount$project <- fuseProject(ifrCount)

dfrCount <- countFusions(def.recurrent)
dfrCount$project <- fuseProject(dfrCount)

## write the output
write.table(star.recurrent, file = "recurrent.star.fusion_candidates.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(map.recurrent, file = "recurrent.mapsplice.fusion_candidates.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ee.recurrent, file = "recurrent.ericscript.results.filtered.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ifr.recurrent, file = "recurrent.integrate.oncofuse.output.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(def.recurrent, file = "recurrent.defuse.results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## save R image
save.image("fusions.data.rda")

## plot for QC
pdf("RecurrentFusions_byProgram.pdf", width = 594/25.4, height = 420/25.4)
ggplot(sfrCount, aes(case, fusionName)) + geom_tile(aes(fill = qual)) +
    facet_grid(. ~ project, scales = "free_x") +
    xlab("Case by Tumor Type") + ylab("STAR-Fusion Recurrent Fusions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mfrCount, aes(case, fusionName)) + geom_tile(aes(fill = qual)) +
    facet_grid(. ~ project, scales = "free_x") +
    xlab("Case by Tumor Type") + ylab("MAPSplice Recurrent Fusions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(efrCount, aes(case, fusionName)) + geom_tile(aes(fill = qual)) +
    facet_grid(. ~ project, scales = "free_x") +
    xlab("Case by Tumor Type") + ylab("EricScript Recurrent Fusions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(ifrCount, aes(case, fusionName)) + geom_tile(aes(fill = qual)) +
    facet_grid(. ~ project, scales = "free_x") +
    xlab("Case by Tumor Type") + ylab("Integrate Recurrent Fusions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(dfrCount, aes(case, fusionName)) + geom_tile(aes(fill = qual)) +
    facet_grid(. ~ project, scales = "free_x") +
    xlab("Case by Tumor Type") + ylab("Defuse Recurrent Fusions") +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
dev.off()
