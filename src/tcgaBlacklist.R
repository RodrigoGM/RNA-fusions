load("fusions.data.rda")
source("myLib.R")

getNormals()

sapply(list(star.normals$case, starOnc.normals$case, ee.normals$case, eeOnc.normals$case, mapOnc.normals$case, map.normals$case), unique)

filterList <- c(star.normals$fusionName, starOnc.normals$fusionName, ee.normals$fusionName, eeOnc.normals$fusionName, mapOnc.normals$fusionName, map.normals$fusionName)
uf <- unique(filterList)

write(filterList, file = "TCGA-Normals/filterList.txt")
