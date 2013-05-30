# 
# Author: Jonathan Rosenblatt<john.ros.work@gmail.com>
###############################################################################

install.packages("FPF", repos="http://R-Forge.R-project.org")
require(FPF)
# Note: manuscript plots were done with Version 0.54 of FPF package. 

require(lattice)

#### Generate Simulation data and prevalence estimates ####
# Warning: long run to generate data and estimates.
source(file='GenerateToyExampleData.R')
levelplot(.result2["p3.1",,], xlab="", ylab="", main="Estimated Signal Prevalence (overlap)" )



wilcoxon.toy.pvals2<- apply(noisy.prevalence2, c(1,2), function(x) wilcox.test(x, mu=0, alternative="two.sided", exact=FALSE)$p.value 	)
# FDR threshing:
wilcoxon.toy.adjusted.pvals2 <- array(p.adjust(wilcoxon.toy.pvals2, method = 'BH'), dim=dim(wilcoxon.toy.pvals2))
plot(wilcoxon.toy.adjusted.pvals2 ~ wilcoxon.toy.pvals2)
FDR <- 0.05
wilcoxon.toy.mask2 <- wilcoxon.toy.adjusted.pvals2 <= FDR
image(wilcoxon.toy.mask2)


wilcoxon.toy.pvals3<- apply(noisy.prevalence3, c(1,2), function(x) wilcox.test(x, mu=0, alternative="two.sided", exact=FALSE)$p.value 	)
wilcoxon.toy.adjusted.pvals3 <- array(p.adjust(wilcoxon.toy.pvals3, method = 'BH'), dim=dim(wilcoxon.toy.pvals3))
plot(wilcoxon.toy.adjusted.pvals3 ~ wilcoxon.toy.pvals3)
wilcoxon.toy.mask3 <- wilcoxon.toy.adjusted.pvals3 <= FDR




#### Create a single image of *six* plots####
plot1 <- levelplot(raw.signal4, xlab="", ylab="", main="A- Single Brain Activation Example", scales=list(draw=FALSE), useRaster=TRUE )
plot2 <- levelplot(signal.prevalence, xlab="", ylab="", main="B- True Activation Prevalence", scales=list(draw=FALSE), useRaster=TRUE)
plot4<- levelplot(.result2["p3.1",,], xlab="", ylab="", main="D- Estimates (misspecified) ", scales=list(draw=FALSE), useRaster=TRUE)
plot5<- levelplot(.result2["p3.1",,] * wilcoxon.toy.mask2, xlab="", ylab="", main="F-Masked Estimates (misspecified)", scales=list(draw=FALSE), useRaster=TRUE)
plot6<- levelplot(.result3["p3.1",,], xlab="", ylab="", main="C-  Estimates (assumed) ", scales=list(draw=FALSE), useRaster=TRUE)
plot7<- levelplot(.result3["p3.1",,] * wilcoxon.toy.mask3, xlab="", ylab="", main="E-Masked Estimates (assumed)", scales=list(draw=FALSE), useRaster=TRUE)
print(plot1, position=c(0, 0.66, 0.5, 1), more=TRUE)
print(plot2, position=c(0.5, 0.66, 1, 1), more=TRUE)
print(plot6, position=c(0, 0.33, 0.5, 0.66), more=TRUE)
print(plot4, position=c(0.5, 0.33, 1, 0.66), more=TRUE)
print(plot5, position=c(0.5, 0, 1, 0.33), more=TRUE)
print(plot7, position=c(0, 0, 0.5, 0.33), more=FALSE)









































