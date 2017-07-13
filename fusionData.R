## merging and analyzing the fusion analysis
## libraries used
library(reshape2)
library(ggplot2)
source("myLib.R")

## tcgaFilter <- unique(readLines("TCGA-Normals/tcga_normals.filter.txt"))
tcgaFilter <- unique(readLines("TCGA-Normals/filterList.txt"))
ffpeFilter <- unique(readLines("FFPE-Normals/ffpeFilter.txt"))
blastFilter <- unique(readLines("FusionFilter/pairwise.tx.blast.txt"))
chimerDB <- unique(readLines("ChimerDB/ChimerDB3_ChimerKB.txt"))
                    
## loading STAR-fusion data
star <- prepStar("merged.star-fusion.fusion_candidates.txt", tcga = tcgaFilter, ffpe = ffpeFilter, knownFusions = chimerDB,  
                 stringsAsFactors = FALSE)
## attach project name
star$project <- fuseProject(star, "case")
## & keep fusions not in TCGA.Normals , i.e. keep all FALSE
star.filtered <- star[!star$inTCGA.Normals, ]
star.filtered <- star.filtered[!star.filtered$inFFPE.Normals, ]
star.filtered <- star.filtered[!star.filtered$inBLAST, ]

## melt of star.filtered fusions
sfCount <- countFusions(star)
sfCount$qual <- setKnown(sfCount)
sfCount$project <- fuseProject(sfCount, "case")

## loading EricScript data
ee <- prepEricScript("merged.ericscript.results.filtered.txt", tcga = tcgaFilter, ffpe = ffpeFilter, knownFusions = chimerDB,  
                     header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
ee$project <- fuseProject(ee, "case")
ee.filtered <- ee[!ee$inTCGA.Normals, ]
ee.filtered <- ee.filtered[!ee.filtered$inFFPE.Normals, ]
ee.filtered <- ee.filtered[!ee.filtered$inBLAST, ]

eeCount <- countFusions(ee)
eeCount$qual <- setKnown(eeCount)
eeCount$project <- fuseProject(eeCount, "case")

## loading MAPSplice data
## run on bash or use system()
## awk '{OFS="\t"; print $4, $6, $6+1, ".", ".", $10}' merged.mapsplice.fusions_candidates.txt | sed -e 's/++/+/' -e 's/+-/+/' -e 's/-+/-/' -e 's/--/-/' | sed '/chrom/d' > mapsplice.out.bed
## awk '{OFS="\t"; print $5, $7, $7+1, ".", ".", $10}' merged.mapsplice.fusions_candidates.txt | sed -e 's/++/+/' -e 's/+-/-/' -e 's/-+/+/' -e 's/--/-/' | sed '/chrom/d' >> mapsplice.out.bed
## bedSort mapsplice.out.bed mapsplice.out.bed
## awk '$6 == "+"' mapsplice.out.bed > mapsplice.pos.bed
## awk '$6 == "-"' mapsplice.out.bed > mapsplice.neg.bed
## not run: sed -e 's/^chr//' -e 's/E.*gene_name \"/ /' -e 's/\".*//' ~/genomes/homo_sapiens/Ensembl/Grch37.p13/Annotation/Genes/gencode.v19.annotation.bed  > gencode.v19.names.tmp
## not run: awk '{OFS="\t"; print "987", $6}' ~/genomes/homo_sapiens/Ensembl/Grch37.p13/Annotation/Genes/gencode.v19.annotation.bed | paste gencode.v19.names.tmp - > gencode.v19.names.bed
## not run: awk '$6 == "+"' gencode.v19.names.bed > gencode.v19.pos.names.bed
## not run: awk '$6 == "-"' gencode.v19.names.bed > gencode.v19.neg.names.bed
## bedmap --delim "\t" --echo --echo-map-id-uniq mapsplice.pos.bed ~/genomes/homo_sapiens/Ensembl/Grch37.p13/Annotation/Genes/gencode.v19.pos.names.bed | awk '{OFS="\t"; print $1, $2, $3, $7, $5, $6}' > mapsplice.annotation.out.bed
## bedmap --delim "\t" --echo --echo-map-id-uniq mapsplice.neg.bed ~/genomes/homo_sapiens/Ensembl/Grch37.p13/Annotation/Genes/gencode.v19.neg.names.bed | awk '{OFS="\t"; print $1, $2, $3, $7, $5, $6}' >> mapsplice.annotation.out.bed

## merging annotation bed files to the mapsplice fusions output
preAnnot <- read.table("merged.mapsplice.fusions_candidates.txt", header = FALSE, sep = "\t", na.strings = "-,", stringsAsFactors = FALSE)[,1:66]
colnames(preAnnot) <- maphh
##preAnnot$doner_strand <- substr(as.character(preAnnot$strand), 1, 1)
##preAnnot$acceptor_strand <- substr(as.character(preAnnot$strand), 2, 2)
annot <- unique(read.table("mapsplice.annotation.out.bed", header = FALSE, fill = TRUE, stringsAsFactors = FALSE))

invisible(
sapply(1:nrow(preAnnot), function(i) preAnnot$annotated_gene_donor[i] <<- as.character(annot$V4[ preAnnot$doner_chrom[i] == annot$V1 & preAnnot$doner_end[i] == annot$V2]))
)
invisible(
sapply(1:nrow(preAnnot), function(i) preAnnot$annotated_gene_acceptor[i] <<- as.character(annot$V4[ preAnnot$acceptor_chrom[i] == annot$V1 & preAnnot$acceptor_start[i] == annot$V2]))
)
write.table(preAnnot, file = "merged.mapsplice.fusions_candidates.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

map <- prepMapsplice("merged.mapsplice.fusions_candidates.txt", tcga = tcgaFilter, ffpe = ffpeFilter, knownFusions = chimerDB,  
                     header = FALSE, fill = TRUE, sep = "\t", stringsAsFactors = FALSE)
map$project <- fuseProject(map, "case")
map.filtered <- map[!map$inTCGA.Normals, ]
map.filtered <- map.filtered[!map.filtered$inFFPE.Normals, ]
map.filtered <- map.filtered[!map.filtered$inBLAST, ]


mapCount <- countFusions(map)
mapCount$qual <- setKnown(mapCount)
mapCount$project <- fuseProject(mapCount, "case")


## loading DEFUSE-fusion data
defuse <- prepDefuse("merged.defuse.results.txt", tcga = tcgaFilter, ffpe = ffpeFilter, knownFusions = chimerDB,  
                     stringsAsFactors = FALSE, na.strings = c(NA, "", "."))
## attach project name
defuse$project <- fuseProject(defuse, "case")
## & keep fusions not in TCGA.Normals , i.e. keep all FALSE
defuse.filtered <- defuse[!defuse$inTCGA.Normals, ]
defuse.filtered <- defuse.filtered[!defuse.filtered$inFFPE.Normals, ]
defuse.filtered <- defuse.filtered[!defuse.filtered$inBLAST, ]

## melt of defuse.filtered fusions
deCount <- countFusions(defuse)
deCount$qual <- setKnown(deCount)
deCount$project <- fuseProject(deCount, "case")


integrate <- prepIntegrate("merged.integrate.oncofuse.txt", tcga = tcgaFilter, ffpe = ffpeFilter, knownFusions = chimerDB,  
                           header = FALSE, skip = 1, fill = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
integrate$project <- fuseProject(integrate, "case")
integrate.filtered <- integrate[!integrate$inTCGA.Normals, ]
integrate.filtered <- integrate.filtered[!integrate.filtered$inFFPE, ]
integrate.filtered <- integrate.filtered[!integrate.filtered$inBLAST, ]

intCount <- countFusions(integrate)
intCount$qual <- setKnown(intCount)
intCount$project <- fuseProject(intCount, "case")


## known fusions
known <- melt(read.table("raw.known.fusions.detected.txt", sep = "\t", header = TRUE), id.vars = "program")
known$value <- factor(known$value, levels = c("No", "Yes"))
oncoKnown <- melt(read.table("oncofuse.annotated.known.fusions.detected.txt", sep = "\t", header = TRUE), id.vars = "program")
oncoKnown$value <- factor(oncoKnown$value, levels = c("No", "Yes"))
known$annotation <- "Standalone"
oncoKnown$annotation <- "Oncofuse"

bk <- rbind(known, oncoKnown)
bk$annotation <- factor(bk$annotation, levels = c("Standalone", "Oncofuse"))
bk$program <- factor(bk$program, levels = rev(c("Integrate", "MAPSplice", "EricScript", "STAR-Fusion")))
bk$value <- factor(bk$value, levels = c("No", "Yes"))

save.image("fusions.data.rda")


##
#### loading STAR-fusion with oncofuse annotation data
##starOnc <- prepOnc("merged.star-fusion.oncofuse.output.txt", tcga = tcgaFilter, ffpe = ffpeFilter, knownFusions = chimerDB,  
##                   header = TRUE, row.names = NULL, fill = TRUE, sep = "\t", stringsAsFactors = FALSE)
##starOnc$project <- fuseProject(starOnc, "case")
##starOnc.filtered <- starOnc[!starOnc$inTCGA.Normals, ]
##starOnc.filtered <- starOnc.filtered[!starOnc.filtered$inFFPE.Normals, ]
##starOnc.filtered <- starOnc.filtered[!starOnc.filtered$inBLAST, ]
##
##sfOncCount <- countFusions(starOnc)
##sfOncCount$qual <- setKnown(sfOncCount)
##sfOncCount$project <- fuseProject(sfOncCount, "case")
##
#### loading EricScript with oncofuse annotation data
##eeOnc <- prepOnc("merged.ericscript.oncofuse.output.txt", tcga = tcgaFilter, ffpe = ffpeFilter, knownFusions = chimerDB,  
##                 header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
##eeOnc$project <- fuseProject(eeOnc, "case")
##eeOnc.filtered <- eeOnc[!eeOnc$inTCGA.Normals, ]
##eeOnc.filtered <- eeOnc.filtered[!eeOnc.filtered$inFFPE.Normals, ]
##eeOnc.filtered <- eeOnc.filtered[!eeOnc.filtered$inBLAST, ]
##
##eeOncCount <- countFusions(eeOnc)
##eeOncCount$qual <- setKnown(eeOncCount)
##eeOncCount$project <- fuseProject(eeOncCount)
##
#### loading MAPSplice with oncofuse annotation
##mapOnc <- prepOnc("merged.mapsplice.oncofuse.output.txt", tcga = tcgaFilter, ffpe = ffpeFilter, knownFusions = chimerDB, 
##                  header = TRUE, sep = "\t", fill = TRUE)
##mapOnc$project <- fuseProject(mapOnc, "case")
##mapOnc.filtered <- mapOnc[!mapOnc$inTCGA.Normals, ]
##mapOnc.filtered <- mapOnc.filtered[!mapOnc.filtered$inFFPE.Normals, ]
##mapOnc.filtered <- mapOnc.filtered[!mapOnc.filtered$inBLAST, ]
##mapOncCount <- countFusions(mapOnc)
##mapOncCount$qual <- setKnown(mapOncCount)
##mapOncCount$project <- fuseProject(mapOncCount, "case")
##
#### read in integrate.oncofuse results
##intOnc <- prepOnc("merged.integrate.oncofuse.output.txt", header = TRUE, fill = TRUE, sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
##intOnc$fusionName <- paste(intOnc$X5_FPG_GENE_NAME, intOnc$X3_FPG_GENE_NAME, sep = "--")
##intOnc$project <- fuseProject(intOnc, "case")
##intOnc.filtered <- intOnc[!intOnc$inTCGA.Normals, ]
##intOnc.filtered <- intOnc.filtered[!intOnc.filtered$inFFPE, ]
##intOnc.filtered <- intOnc.filtered[!intOnc.filtered$inBLAST, ]
##
##intOncCount <- countFusions(intOnc)
##intOncCount$qual <- setKnown(intOncCount)
##intOncCount$project <- fuseProject(intOncCount, "case")
##
##
#### merging oncofuse outputs
## mmOnc <- merge(mapOnc, eeOnc, by = "fusionName")
## mmOnc <- merge(mmOnc, starOnc, by = "fusionName")
## mmOnc <- merge(mmOnc, starOnc, by = "fusionName")
## write.table(table(mmOnc$fusionName, mmOnc$case), file = "intersect.oncofuse.star.ericScript.mapsplice.counts.txt", quote = FALSE, sep = "\t")
##
## mmOncCount <- countFusions(mmOnc)
## mmOncCount$qual <- setKnown(mmOncCount)
## mmOncCount$project <- fuseProject(mmOncCount, "case")

#### merging raw outputs
## mm <- merge(ee, map, by = "fusionName")
## mm <- merge(mm, star, by = "fusionName")
## mm <- merge(mm, defuse, by = "fusionName")
## write.table(table(mm$fusionName, mm$case), file = "intersect.star.ericScript.mapsplice.counts.txt", quote = FALSE, sep = "\t")
## 
## mmCount <- countFusions(mm)
## mmCount$qual <- setKnown(mmCount)
## mmCount$project <- fuseProject(mmCount, "case")
##
## 
### merging filtered raw outputs
## mmf <- merge(star.filtered, ee.filtered, by = "fusionName")
## mmf <- merge(mmf, map.filtered, by = "fusionName")
## write.table(table(mmf$fusionName, mmf$case), file = "intersect.tcgaFilter.ffpeFilter.star.ericScript.mapsplice.counts.txt", quote = FALSE, sep = "\t")
## 
## mmfCount <- countFusions(mmf)
## mmfCount$qual <- setKnown(mmfCount)
## mmfCount$project <- fuseProject(mmfCount, "case")
## 
## merging oncofuse filtered output
## mmfOnc <- merge(starOnc.filtered, eeOnc.filtered, by = "fusionName")
## mmfOnc <- merge(mmfOnc, mapOnc.filtered, by = "fusionName")
## mmfOnc <- merge(mmfOnc, defuse.filtered, by = "fusionName")
## mmfOnc <- merge(mmfOnc, integrate.filtered, by = "fusionName")
## write.table(table(mmfOnc$fusionName, mmfOnc$case), file = "intersect.tcga.oncofuse.star.ericScript.mapsplice.integrate.defuse.counts.txt", quote = FALSE, sep = "\t")
## 
## mmfOncCount <- countFusions(mmfOnc)
## mmfOncCount$qual <- setKnown(mmfOncCount)
## mmfOncCount$project <- fuseProject(mmfOncCount, "case")

