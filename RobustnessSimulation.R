# 
# Author: johnros
###############################################################################


#### load required functions ####
rm(list=ls())
source(file='RobustnessSimulationFunctions.R')

## Generate simulation configurations:
prevalences <- 	makePrevalences(prevalences = c(0.05, 0.2, 0.5, 0.8, 0.95), prop.ratio=c(0, 0.05, 0.1)) 
parameters <- list(
    		sample.size=c(15,20,30,50,75,100,150,200),
		sd1s = c(sqrt(0.15)),
		sd2s = c(1),
		sd3s = c(sqrt(0.6)),
		sd4s = c(sqrt(1.8)),# robustness is achieved with mixing proportions
		mu= c(0.8,1,2) 
		)
library(reshape)
configurations <- expand.grid.df(do.call(expand.grid, parameters), prevalences)



## Compute MSE for each configuration.
install.packages("FPF", repos="http://R-Forge.R-project.org")
library(FPF)


# Main function. Will require some time to run....
MSEs <- computeMSE(configurations, replications = 200 )


# Now might be a good time to save() MSEs object.

# In case you want an email sent when simulation finishes (linux only. requires ssmpt)
system("ssmtp myemail@something.com < someTextFile.txt")




#### Analyze results ####
require(ggplot2)
MSEs$kurtosis<-  with(MSEs, ScaleMixtureKurtosis(p3/(p3+p4), sd3s, sd4s, standardized=TRUE))
MSEs$p <- with(MSEs, p3+p4)
MSEs$xi <- with(MSEs, round(p4/(p3+p4),2))
MSEs$RMSE<- with(MSEs, sqrt(MSE))
MSEs$RejectedRMSE<- with(MSEs, sqrt(RejectedMSE))
MSEs$RejectedRMSE2<- with(MSEs, sqrt(RejectedMSE2))
MSEs$root.sample.size<- with(MSEs, sqrt(sample.size))
MSEs$oracle<- with(MSEs, sqrt(prevalence*(1-prevalence)/sample.size))




myplot<- function(y, data){
  base.plot<- ggplot(aes_string(y=as.character(match.call()[2]) , x="sample.size"), data=data) + xlab("Sample Size")
  base.plot.2<- base.plot + geom_point() + facet_grid(mu+p ~  xi, labeller = label_both)
  base.plot.2 + ylim(c(0, 0.6)) + scale_x_sqrt() + geom_vline(xintercept = 32, lty=2, col='darkgrey') #+ stat_smooth(se = FALSE, method=loess, span = 0.9)  
}

# Choose the required mu value:
table(MSEs$mu)
.subset<- subset(MSEs, mu==2 )

# Plots used in manuscript.
myplot(RMSE, subset(MSEs, mu==1 ))
myplot(RMSE, subset(MSEs, mu==2 ))

# Plots not used. Demostrating selective accuracy.
myplot(RejectedRMSE, subset(MSEs, mu==1 ))
myplot(RejectedRMSE2, subset(MSEs, mu==1 ))




