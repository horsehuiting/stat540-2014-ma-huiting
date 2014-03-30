library(edgeR)
dat <- read.table("bottomly_count_table.tsv", header = TRUE, 
                  row.names = 1)
des <- read.table("bottomly_phenodata.tsv", header = TRUE, 
                  row.names = 1)
str(dat)
show(des)
all(rownames(des) == colnames(dat))
with(des, table(strain))

group <- factor(c(rep("1", 10), rep("2", 11)))
group

# this produces an object of type DGEList with can be manipulated in a
# similar way to any other list object in R
dge.glm <- DGEList(counts = dat, group = group)
names(dge.glm)
dge.glm[["samples"]]
nrow(dge.glm[[1]])
ncol(dge.glm[[1]])
str(dge.glm)
design <- model.matrix(~group)
design
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp, design)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)
# plot the tagwise dispersion against log2-CPM (counts per million)
plotBCV(dge.glm.tag.disp)
fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)
nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])
interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]
summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))
# plotting the tagwise log fold changes against log-cpm
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")






library(DESeq)
# reading in the same count table data and grouping information
deSeqDat <- newCountDataSet(dat, group)
head(counts(deSeqDat))
deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)
deSeqDat <- estimateDispersions(deSeqDat)
# plotting the estimated dispersions against the mean normalized counts
plotDispEsts(deSeqDat)
## this takes a minute or so for JB
results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)
plotMA(results)



library(limma)
norm.factor <- calcNormFactors(dat)
dat.voomed <- voom(dat, design, plot = TRUE, 
                   lib.size = colSums(dat) * norm.factor)
dat.voomed
fit <- lmFit(dat.voomed, design)
fit <- eBayes(fit)
topTable(fit)


## Take Home Problem

## Choose a specific threshold for the adjusted p value, 
## find the genes identified as differentially expressed using each of edgeR, 
## DESeq and voom+limma. Compare the number of genes in these 3 lists, 
## and draw a venn digram demonstrating the overlap (if any!).
cut <- 1e-4
edgeR <- subset(tt.glm$table, FDR < cut)
edgerR <- na.omit(edgeR)
DeSeq <- subset(results, padj < cut)
DeSeq <- na.omit(DeSeq)
hits.voomed <- topTable(fit, coef ="group2",  n = Inf)
voomed <- subset(hits.voomed, adj.P.Val < cut)
voomed<- na.omit(voomed)

library(VennDiagram)

edger.genes <- rownames(edgeR)
deseq.genes <- DeSeq$id
voom.limma.genes <- rownames(voomed)

# Put the things you want to plot in a list. The names in the list will be
# put on the plot.
de.genes <- list(edger = edger.genes, deseq = deseq.genes, voom.limma = voom.limma.genes)

# Start a new plot
plot.new()

# Draw the Venn diagram. Note the argument `filename=NULL` tells it to
# create a plot object instead of outputting to file.
venn.plot <- venn.diagram(de.genes, filename = NULL, fill = c("red", "blue", "green"))

# Draw the plot on the screen.
grid.draw(venn.plot)




## Mini Exercise

##Redo the above analysis but first filter the data and remove any gene that has: 
## 1. count equal tot zero across all samples 
## 2. count equal to zero in at least one sample in each genotype group
dat <- read.table("bottomly_count_table.tsv", header = TRUE, 
                  row.names = 1)
des <- read.table("bottomly_phenodata.tsv", header = TRUE, 
                  row.names = 1)

dat <- subset(dat, rowSums(dat) == 0 & (all(dat[, 1:10] != 0)) & (all(dat[, 11:21]!=0)))
all(rownames(des) == colnames(dat))
with(des, table(strain))

group <- factor(c(rep("1", 10), rep("2", 11)))
group

# this produces an object of type DGEList with can be manipulated in a
# similar way to any other list object in R
dge.glm <- DGEList(counts = dat, group = group)
names(dge.glm)
dge.glm[["samples"]]
nrow(dge.glm[[1]])
ncol(dge.glm[[1]])
str(dge.glm)
design <- model.matrix(~group)
design
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp, design)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)
# plot the tagwise dispersion against log2-CPM (counts per million)
plotBCV(dge.glm.tag.disp)
fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)
nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])
interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]
summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))
# plotting the tagwise log fold changes against log-cpm
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")
