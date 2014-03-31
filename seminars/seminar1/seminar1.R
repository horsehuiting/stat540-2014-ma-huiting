prDat <- read.table("GSE4051_MINI.txt", header = TRUE, row.names = 1)
str(prDat)
head(prDat)
rownames(prDat)

# Exercise: Why do we use the header = and row.names = arguments above upon import? 
# Submit the command without these arguments and note any difference in the result. 
# Form the habit of reading error messages carefully and working the problem. 
# Mastering the arguments of read.table() and friends is time well spent.

prDat1 <- read.table("GSE4051_MINI.txt")
str(prDat1)
head(prDat1)
rownames(prDat1)
## I did not notice any difference, it might due to read.table default. Let's try make header = FALSE

prDat2 <- read.table("GSE4051_MINI.txt", header = FALSE)
str(prDat2)
head(prDat2)

## can't read it now