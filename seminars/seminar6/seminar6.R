library(limma)
library(lattice)
library(plyr)

prDat <- read.table("GSE4051_data.tsv",header=T)
prDes <- readRDS("GSE4051_design.rds")
str(prDat, max.level = 0)
str(prDes)

prepareData <- function(myGenes) {
  miniDat <- t(wtDat[myGenes, ])
  miniDat <- data.frame(gExp = as.vector(miniDat),
                        gene = factor(rep(colnames(miniDat), each =
                                            nrow(miniDat)), levels = colnames(miniDat)))
  miniDat <- suppressWarnings(data.frame(wtDes, miniDat))
  miniDat
}
stripplotIt <- function(myData, ...) {
  stripplot(gExp ~ devStage | gene, myData,
            jitter.data = TRUE,
            auto.key = TRUE, type = c('p', 'a'), grid = TRUE, ...)
}

## Reproduce all seminar 6 by following instructions
m <- 1000
n <- 3
x <- matrix(rnorm(m * n,0,3), nrow = m)
obsVars <- apply(x, 1, var)
summary(obsVars)
mean(obsVars < 1/3)
densityplot(~obsVars, n = 200)

wtDes <- subset(prDes, gType == "wt")
str(wtDes)

wtDat <- subset(prDat, select = prDes$gType == "wt")
str(wtDat, max.level = 0)

wtDesMat <- model.matrix(~devStage, wtDes)
str(wtDesMat)

wtFit <- lmFit(wtDat, wtDesMat)
wtEbFit <- eBayes(wtFit)
topTable(wtEbFit)
topTable(wtEbFit, coef = 2:5)  # cryptic! error-prone!
colnames(coef(wtEbFit)) 
(dsHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit)))))
stripplotIt(prepareData(rownames(dsHits)[c(3, 6, 9)]))

#dat <- prepareData(rownames(dsHits)[3])
#Fit <- lm(gExp~devStage,dat)
#summary(Fit)

cutoff <- 1e-05
dsHits <- topTable(wtEbFit,
                   coef = grep("devStage", colnames(coef(wtEbFit))),
                   p.value = cutoff, n = Inf)
(numBHhits <- nrow(dsHits))

dsHits[63, c("F", "adj.P.Val", "devStageP6")]

P2Hits <- topTable(wtEbFit, coef = "devStageP2", n = Inf, sort = "none")
P10Hits <- topTable(wtEbFit, coef = "devStageP10", n = Inf, sort = "none")
xyplot(P10Hits$t ~ P2Hits$t, aspect = 1,
       xlab = "t-statistic for P2 effect",
       ylab = "t-statistic for P10 effect",
       xlim = c(-20, 16), ylim = c(-20, 16),
       panel = function(x, y, ...) {
         panel.smoothScatter(x, y, nbin = 100, ...)
         panel.abline(a = 0, b = 1, col = "orange")
       })


densityplot(~ P10Hits$adj.P.Val + P2Hits$adj.P.Val, auto.key = TRUE,
            plot.points = FALSE, n = 300)

cutoff <- 1e-03
foo <- data.frame(P2 = P2Hits$adj.P.Val < cutoff,
                  P10 = P10Hits$adj.P.Val < cutoff)
addmargins(with(foo, table(P2, P10)))

P10pVals <- data.frame(raw = P10Hits$P.Value,
                       BH = P10Hits$adj.P.Val,
                       BY = p.adjust(P10Hits$P.Value, method = "BY"))
splom(P10pVals,
      panel = function(x, y, ... ) {
        panel.xyplot(x, y, pch = ".", ...)
        panel.abline(a = 0, b = 1, col = "orange")
      })


colnames(wtDesMat)

(cont.matrix <- makeContrasts(P10VsP6 = devStageP10 - devStageP6, fourweeksVsP10 = devStage4_weeks - 
                                devStageP10, levels = wtDesMat))
wtFitCont <- contrasts.fit(wtFit, cont.matrix)
wtEbFitCont <- eBayes(wtFitCont)
topTable(wtEbFitCont)

foo <- topTable(wtEbFitCont)
stripplotIt(prepareData(rownames(foo)[1:4]))

cutoff <- 1e-04
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)

(hits1 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] < 0)])
stripplotIt(prepareData(hits1))
(hits2 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] < 0)])
stripplotIt(prepareData(hits2[1:4]))
intersect(hits1, hits2)
(hits3 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] > 0)])
stripplotIt(prepareData(hits3[1:4]))
intersect(hits1, hits3)
intersect(hits2, hits3)

cutoff <- 0.01
nHits <- 8
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)

hits1 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] < 0)]
stripplotIt(prepareData(hits1[1:nHits]))
hits2 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] < 0)]
stripplotIt(prepareData(hits2[1:nHits]))
hits3 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] > 0)]
stripplotIt(prepareData(hits3[1:nHits]))
hits4 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] > 0)]
stripplotIt(prepareData(hits4[1:nHits]))

vennDiagram(wtResCont)
hits5 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] != 0 & wtResCont[, "fourweeksVsP10"] != 
                                 0)]
stripplotIt(prepareData(hits5))

hits6 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] > 0 & wtResCont[, "fourweeksVsP10"] < 
                                 0)]
stripplotIt(prepareData(hits6))

## Exercise 1
## Make the above simulation more realistic with two (or more) groups, different data-generating 
## means and group differences, different data-generating gene-wise variances, etc.
m <- 1000
n <- 3
x <- matrix(rnorm(m/2 * n,0,3), nrow = m/2)
y <- matrix(rnorm(m/2 * n,3,1), nrow = m/2)
z <- rbind(x,y)

obsVars <- apply(z, 1, var)
summary(obsVars)
mean(obsVars < 1/3)
densityplot(~obsVars, n = 200)

## Take Home
#colnames(wtDesMat)
(cont.matrixM <- makeContrasts(P2VsE16 = devStageP2 - Intercept,
                               P6VsP2 = devStageP6 - devStageP2, 
                               P10VsP6 = devStageP10 - devStageP6,
                               fourweeksVsP10 = devStage4_weeks - devStageP10,
                               levels = wtDesMat))
wtFitCont1 <- contrasts.fit(wtFit, cont.matrixM)
wtEbFitCont1 <- eBayes(wtFitCont1)
topTable(wtEbFitCont1)
cutoff <- 1e-04
wtResCont1 <- decideTests(wtEbFitCont1, p.value = cutoff, method = "global")
summary(wtResCont1)

hits7 <- rownames(prDat)[which(wtResCont1[,"P2VsE16"] < 0 
                               & wtResCont1[,"P6VsP2"] > 0 &
                                 wtResCont1[, "P10VsP6"] == 0 & 
                                 wtResCont1[, "fourweeksVsP10"] == 0)]

stripplotIt(prepareData(hits7))