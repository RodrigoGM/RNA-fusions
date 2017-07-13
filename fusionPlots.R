## code for plotting fusion data

library(ggplot2)
library(reshape2)
library(RColorBrewer)

source("myLib.R")

load("fusions.data.rda")

pdf("FusionPrediction_per_program.pdf", width = 594/25.4, height = 420/25.4)
ggplot(sfCount[!is.na(sfCount$qual),], aes(as.character(case), as.character(fusionName))) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("STAR-Fusion Predictions") + theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mapCount[!is.na(mapCount$qual),], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("MAPSplice Predictions") + theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(eeCount[!is.na(eeCount$qual),], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("EricScript Predictions") + theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
lapply(unique(sfCount$project), function(PR) {
    ggplot(sfCount[sfCount$project == PR & !is.na(sfCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free_x") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("STAR-Fusion Predictions") + theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
})
lapply(unique(mapCount$project), function(PR) {
    ggplot(mapCount[mapCount$project == PR & !is.na(mapCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free_x") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("MAPSplice Predictions") + theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
})
lapply(unique(eeCount$project), function(PR) {
    ggplot(eeCount[eeCount$project == PR & !is.na(eeCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free_x") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("Ericscript Predictions")  + theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
})
lapply(unique(intCount$project), function(PR) {
    ggplot(intCount[intCount$project == PR & !is.na(intCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free_x") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("Ericscript Predictions")  + theme(axis.text.y=element_blank()) +
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
})
dev.off()


pdf("RecurrentFusions_per_program.pdf", width = 594/25.4, height = 420/25.4)
ggplot(sfrCount[!is.na(sfrCount$qual),], aes(as.character(case), as.character(fusionName))) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("STAR-Fusion Predictions") +  ## + theme(axis.text.y=element_blank())
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(mfrCount[!is.na(mfrCount$qual),], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("MAPSplice Predictions") +  ## + theme(axis.text.y=element_blank())
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(efrCount[!is.na(efrCount$qual),], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("EricScript Predictions") +  ## + theme(axis.text.y=element_blank())
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
ggplot(ifrCount[!is.na(ifrCount$qual),], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("Integrate Predictions") +  ## + theme(axis.text.y=element_blank())
    scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                      labels = c("Not Present", "Predicted"))
lapply(unique(sfrCount$project), function(PR) {
    ggplot(sfrCount[sfrCount$project == PR & !is.na(sfrCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free_x") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("STAR-Fusion Predictions") +  ## + theme(axis.text.y=element_blank())
        scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                          labels = c("Not Present", "Predicted"))
})
lapply(unique(mfrCount$project), function(PR) {
    ggplot(mfrCount[mfrCount$project == PR & !is.na(mfrCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free_x") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("MAPSplice Predictions") +  ## + theme(axis.text.y=element_blank())
        scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                          labels = c("Not Present", "Predicted"))
})
lapply(unique(efrCount$project), function(PR) {
    ggplot(efrCount[efrCount$project == PR & !is.na(efrCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free_x") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("Ericscript Predictions") +  ## + theme(axis.text.y=element_blank())
        scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                          labels = c("Not Present", "Predicted"))
})
lapply(unique(ifrCount$project), function(PR) {
    ggplot(ifrCount[ifrCount$project == PR & !is.na(ifrCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free_x") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("Integrate Predictions") +  ## + theme(axis.text.y=element_blank())
        scale_fill_manual(values = c("#CFCFCF", "#F8766D"), name = "Fusion", breaks = c("Not Present", "Predicted"),
                          labels = c("Not Present", "Predicted"))
})
dev.off()




pdf("FusionPrediction_wOncofuse_per_program.pdf", width = 594/25.4, height = 420/25.4)
ggplot(sfOncCount[!is.na(sfOncCount$qual),], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("Oncofuse Annotated STAR-Fusion Predictions")  + theme(axis.text.y=element_blank())
ggplot(mapOncCount[!is.na(mapOncCount$qual),], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("Oncofuse Annotated MAPSplice Predictions")  + theme(axis.text.y=element_blank())
ggplot(eeOncCount[!is.na(eeOncCount$qual),], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual)), width = .95, height = .95) +
    theme(axis.text.x  = element_text(angle = 90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
    xlab("Case") + ylab("Oncofuse Annotated EricScript Predictions") + theme(axis.text.y=element_blank())
dev.off()


pdf("FusionPrediction_by_project_wOncofuse_per_program.pdf", width = 594/25.4, height = 420/25.4)
lapply(unique(sfOncCount$project), function(PR) {
    ggplot(sfOncCount[sfOncCount$project == PR & !is.na(sfOncCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("STAR-Fusion Predictions") + theme(axis.text.y=element_blank())
})
lapply(unique(mapOncCount$project), function(PR) {
    ggplot(mapOncCount[mapOncCount$project == PR & !is.na(mapOncCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("MAPSplice Predictions") + theme(axis.text.y=element_blank())
})
lapply(unique(eeOncCount$project), function(PR) {
    ggplot(eeOncCount[eeOncCount$project == PR & !is.na(eeOncCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("EricScript Predictions") + theme(axis.text.y=element_blank())
})
dev.off()

## plot 3-way merged output
pdf("Merged_Unfiltered_Map_Er_StF.pdf", width = 420/25.4, height = 297/25.4)
ggplot(mmCount, aes(case, fusionName, fill = factor(value))) + geom_tile(width = .95, height = .95) +
    theme( axis.text.x  = element_text(angle=90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse=TRUE, title = "Fusion\nIterations")) +
    xlab("Case") + ylab("Shared Fusion Predictions") 
#    scale_fill_discrete(breaks=c())
ggplot(mmCount, aes(case, fusionName, fill = factor(qual))) + geom_tile(width = .95, height = .95) +
    theme( axis.text.x  = element_text(angle=90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse=TRUE, title = "Fusion")) +
    xlab("Case") + ylab("Shared Fusion Predictions") 
lapply(unique(mmCount$project)[-c(5,9)], function(PR) {
    ggplot(mmCount[mmCount$project == PR & !is.na(mmCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("Shared Fusion Predictions")
})
dev.off()

## merged filtered
pdf("Merged_Filtered_Map_Er_StF.pdf", width = 420/25.4, height = 297/25.4)
ggplot(mmfCount, aes(case, fusionName, fill = factor(value))) + geom_tile(width = .95, height = .95) +
    theme( axis.text.x  = element_text(angle=90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse=TRUE, title = "Fusion\nIterations")) +
    xlab("Case") + ylab("Shared Fusion Predictions") 
ggplot(mmfCount, aes(case, fusionName, fill = factor(qual))) + geom_tile(width = .95, height = .95) +
    theme( axis.text.x  = element_text(angle=90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse=TRUE, title = "Fusion")) +
    xlab("Case") + ylab("Shared Fusion Predictions") 
lapply(unique(mmfCount$project)[-c(5,9)], function(PR) {
    ggplot(mmfCount[mmfCount$project == PR & !is.na(mmfCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("Shared Fusion Predictions")
})
dev.off()

## plot 3-way merged oncofuse annotated
pdf("Merged_Unfiltered_Oncofuse_Map_Er_StF.pdf", width = 420/25.4, height = 297/25.4)
ggplot(mmOncCount, aes(case, fusionName, fill = factor(value))) + geom_tile(width = .95, height = .95) +
    theme( axis.text.x  = element_text(angle=90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse=TRUE, title = "Fusion\nIterations")) +
    xlab("Case") + ylab("Oncofuse Annotated Shared Fusion Predictions") 
#    scale_fill_discrete(breaks=c())
ggplot(mmOncCount, aes(case, fusionName, fill = factor(qual))) + geom_tile(width = .95, height = .95) +
    theme( axis.text.x  = element_text(angle=90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse=TRUE, title = "Fusion")) +
    xlab("Case") + ylab("Oncofuse Annotated Shared Fusion Predictions") 
lapply(unique(mmOncCount$project)[-c(5, 11)], function(PR) {
    ggplot(mmOncCount[mmOncCount$project == PR & !is.na(mmOncCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
        facet_grid(. ~ project, scales = "free") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("Shared Fusion Predictions")
})
dev.off()

pdf("Merged_Filtered_Oncofuse_Map_Er_StF.pdf", width = 420/25.4, height = 297/25.4)
ggplot(mmfOncCount, aes(case, fusionName, fill = factor(value))) + geom_tile(width = .95, height = .95) +
    theme( axis.text.x  = element_text(angle=90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse=TRUE, title = "Fusion\nIterations")) +
    xlab("Case") + ylab("Oncofuse Annotated Shared Fusion Predictions") 
ggplot(mmfOncCount, aes(case, fusionName, fill = factor(qual))) + geom_tile(width = .95, height = .95) +
    theme( axis.text.x  = element_text(angle=90)) + facet_grid(. ~ project, scales = "free_x") +
    guides(fill = guide_legend(reverse=TRUE, title = "Fusion")) +
    xlab("Case") + ylab("Oncofuse Annotated Shared Fusion Predictions") 
lapply(unique(mmfOncCount$project)[-c(5, 11)], function(PR) {
    ggplot(mmfOncCount[mmfOncCount$project == PR & !is.na(mmfOncCount$qual), ], aes(case, fusionName)) + geom_tile(aes(fill = as.factor(qual))) +
#        facet_grid(. ~ project, scales = "free") +
        guides(fill = guide_legend(reverse = TRUE, title = "Fusion")) +
        xlab("Case") + ylab("Shared Fusion Predictions")
})
dev.off()


pdf("Known_Fusions_summary.pdf", width = 297/25.4, height = 148/25.4)
ggplot(known, aes(variable, program)) + geom_tile(aes(fill = (value)), width = 0.95, height = 0.95) + xlab("Case and Known Fusion") + ylab("Program") + guides(fill = guide_legend(reverse=TRUE, title = ""))
ggplot(oncoKnown, aes(variable, program)) + geom_tile(aes(fill = (value)), width = 0.95, height = 0.95) + xlab("Case and Known Fusion") + ylab("Program w/Oncofuse annotation") + guides(fill = guide_legend(reverse=TRUE, title = ""))
ggplot(bk, aes(variable, program)) + geom_tile(aes(fill = (value)), width = 0.95, height = 0.95) + facet_grid(. ~ annotation) + xlab("Case and Known Fusion") + ylab("Program w/Oncofuse annotation") + guides(fill = guide_legend(reverse=TRUE, title = ""))
dev.off()


## Venn-Euler diagrams
## vve <- venneuler(c(STARFusion =  nrow(star), MapSplice = nrow(map), EricScript = nrow(ee),
##                    "STARFusion&MapSplice" = sum(as.character(star$fusionName) %in% as.character(map$fusionName)),
##                    "STARFusion&EricScript" = sum(as.character(star$fusionName) %in% as.character(ee$fusionName)),
##                    "MapSplice&EricScript" = sum(as.character(map$fusionName) %in% as.character(ee$fusionName)),
##                    "STARFusion&MapSplice&EricScript" = nrow(mm))
##                  )

## vvef <- venneuler(c(STARFusion =  nrow(star.filtered), MapSplice = nrow(map.filtered), EricScript = nrow(ee.filtered),
##                     "STARFusion&MapSplice" = sum(as.character(star.filtered$fusionName) %in% as.character(map.filtered$fusionName)),
##                     "STARFusion&EricScript" = sum(as.character(star.filtered$fusionName) %in% as.character(ee.filtered$fusionName)),
##                     "MapSplice&EricScript" = sum(as.character(map.filtered$fusionName) %in% as.character(ee.filtered$fusionName)),
##                     "STARFusion&MapSplice&EricScript" = nrow(mmOnc))
##                   )

## adjust plot in inkscape
## pdf("SharedFusions_filtered.pdf")
## plot(vve)
## text(vve$centers, labels = c(nrow(star.filtered), nrow(map.filtered), nrow(ee.filtered)), pos = 1)
## plot(vvef)
## text(vvef$centers, labels = c(nrow(star.filtered), nrow(map.filtered), nrow(ee.filtered)), pos = 1)
## dev.off()

## vvo <- venneuler(c(STARFusion = 1483, MapSplice = 757, EricScript = 14493, "STARFusion&MapSplice" = 103, "STARFusion&EricScript" = 162, "MapSplice&EricScript" = 112, "STARFusion&MapSplice&EricScript" = 61))
## pdf("SharedFusions_oncofuse.pdf")
## plot(vvo)
## dev.off()
##options( java.parameters = "-Xmx16g" )
##vvm <- data.frame(program = c(star.filtered$program, ee.filtered$program, map.filtered$program), fusionName = c(star.filtered$fusionName, ee.filtered$fusionName, map.filtered$fusionName))
##vv <- venneuler(vvm)
