library(lattice)
library(RColorBrewer)
library(gplots)
library(hexbin)

kDat <- read.table("GSE4051_MINI.txt")
str(kDat)
table(kDat$devStage)

table(kDat$gType)
with(kDat,table(devStage,gType))

prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)

prDes <- read.csv("GSE4051_design.csv")
str(prDes)

set.seed(540)
(yo <- sample(1:nrow(prDat), size = 20,replace=TRUE))
hDat <- prDat[yo, ]
str(hDat)
hDat <- as.matrix(t(hDat))
rownames(hDat) <- with(prDes,
                       paste(devStage, gType, sidChar, sep="_"))
str(hDat)
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8))
heatmap(hDat, Rowv = NA, Colv = NA, col = cm.colors(256),
        scale="none", margins = c(5, 8))
jGraysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8),
        col = jGraysFun(256))

heatmap.2(hDat, col = jGraysFun, trace = "none")

heatmap.2(cor(t(hDat)), col = jGraysFun, trace = "none")

set.seed(540)
(yo <- sample(1:ncol(prDat), size = 2))
pairDat <- subset(prDat, select = yo)
str(pairDat)
pairs(pairDat)

splom(pairDat, panel = panel.smoothScatter, raster = TRUE)
hexplom(pairDat)
                 