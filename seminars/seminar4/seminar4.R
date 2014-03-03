library(lattice)
library(plyr)
library(xtable)
prDat <- read.table("GSE4051_data.tsv",header=T)
prDes <- readRDS("GSE4051_design.rds")
str(prDat, max.level = 0)


set.seed(540)
KeepGenes <- sample(1:nrow(prDat), 100)
KeepGenes <- rownames(prDat)[KeepGenes]
miniDat <- subset(prDat, rownames(prDat) %in% KeepGenes)
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                      gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                    levels = rownames(miniDat)))
miniDat <- suppressWarnings(data.frame(prDes, miniDat))
str(miniDat)

pvalue <- ddply(miniDat, ~ gene, function(z) {
  t <- suppressWarnings(t.test(gExp ~ gType, z))
  wc <- suppressWarnings(wilcox.test(gExp ~ gType, z))
  ks <- suppressWarnings(ks.test(z$gExp[z$gType == "NrlKO"],
        z$gExp[z$gType == "wt"]))
  c(pvalue_t=t$p.value,pvalue_ws=wc$p.value,pvalue_ks=ks$p.value)
})
#head(pvalue)
pvalue

## get all columns except names
p_value <- pvalue[,-1]

##change to binary
binary <- p_value <= 0.05
binary <- as.data.frame(binary)
rownames(binary) <- pvalue[,1]

(sig_by_test <- apply(binary, 2, sum)) # 2 means column sum

## Choose interesting Genes (sig in all tests)
interest <- rownames(subset(binary, apply(binary,1,sum)==3))

miniDat <- subset(prDat, rownames(prDat) %in% interest)
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                      gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                    levels = interest))
miniDat <- suppressWarnings(data.frame(prDes, miniDat))
str(miniDat)

stripplot(gType ~ gExp | gene, miniDat,
          scales = list(x = list(relation = "free")),
          group = gType, auto.key = TRUE)
