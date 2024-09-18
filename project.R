library(edgeR)
library(baySeq)
library(RMySQL)
cleaned <- read.csv("cleaned.csv")
View(cleaned)
head(cleaned)
mobDataGroups <- c("Mock","Mock","Mock","dCC_MeCP2","dCC_MeCP2","dCC_MeCP2","dCC_scramble","dCC_scramble","dCC_scramble")
data(mobAnnotation)
head(mobAnnotation)
d <- DGEList(counts=cleaned,group=factor(mobDataGroups))
d
dim(d)
d.full <- d # keep the old one in case we mess up
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum) # total gene counts per sample
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d)
d
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

'''
#not run
require(DESeq)
cds <- newCountDataSet( data.frame(d$counts), d$samples$group )
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cds <- estimateDispersions( cds , method="blind")
plotDispEsts(cds)
#not run
'''

et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
et13 <- exactTest(d1, pair=c(1,3)) # compare groups 1 and 3
et23 <- exactTest(d1, pair=c(2,3)) # compare groups 2 and 3
#Exporting the data sets
write.csv(et12, file = "et12.csv")
write.csv(et13, file = "et13.csv")
write.csv(et23, file = "et23.csv")

#1-dCC_MeCP2
#2-dCC_scramble
#3-Mock

top13<-topTags(et13, n=2745)
write.csv(top13,file="top13.csv")
top12<-topTags(et12, n=2745)
write.csv(top12,file="top12.csv")
top23<-topTags(et23, n=2745)
write.csv(top23,file="top23.csv")
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
de1 <- decideTestsDGE(et13, adjust.method="BH", p.value=0.05)
summary(de1)
de1tags12 <- rownames(d1)[as.logical(de1)]
de1tags13 <- rownames(d1)[as.logical(de1)]

plotSmear(et13, de.tags=de1tags13)
abline(h = c(-2, 2), col = "blue")
design.mat
fit <- glmFit(d2, design.mat)
# compare (group 1 - group 2) to 0:
# this is equivalent to comparing group 1 to group 2
lrt12 <- glmLRT(fit, contrast=c(1,-1,0))
lrt13 <- glmLRT(fit, contrast=c(1,0,-1))
lrt23 <- glmLRT(fit, contrast=c(0,1,-1))
topTags(lrt12, n=10)
#topTags(lrt13, n=10)
#topTags(lrt23, n=10)
de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")


