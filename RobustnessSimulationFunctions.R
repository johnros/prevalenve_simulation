# Make prevalence values to be tested 
makePrevalences <-function (prevalences, prop.ratio){
  #prevalences <- seq(0,1,0.2)
	#prop.ratio <- seq(0,1,length=10)
	base <- c(p1=0.5, p2=0.5) 
	combinations <- expand.grid(prevalence=prevalences, prop.ratio=prop.ratio)
	ps <- matrix(NA, ncol=4, nrow=nrow(combinations),dimnames=list(NULL,c('p1','p2','p3','p4')))
	for(i in 1:nrow(combinations)){
		ps[i,] <- with(combinations[i,], c((1-prevalence)*base, prevalence*c((1-prop.ratio), prop.ratio)))
	}
	
	return(ps)    
}

# Computes the square error given a sample and true parameter value.
getPrevalenceSquaredError <- function(data, true.prevalence){
  #data <- .data[1,] # .data is created in the following
	estimates <- FPF:::pointWiseMixtureFitFast(data, generateMixtureControl())
	return(	(estimates[['p3.1']]-true.prevalence)^2 )
}

# Get MSE for rejected locations only
getRejectedPrevalenceSquaredError <- function(data, true.prevalence, alpha=0.05){
	#data <- .data[1,] # .data is created in the following
	rejection <- wilcox.test(data, alternative='two.sided')$p.value < alpha
	if(!rejection)	return(	( 0-true.prevalence)^2 )
	else{
		estimates <- FPF:::pointWiseMixtureFitFast(data, generateMixtureControl())
		return(	(estimates[['p3.1']]-true.prevalence)^2 )
	}
}


getRejectedPrevalenceSquaredError2 <- function(data, true.prevalence, alpha=0.05){
	#data <- .data[1,] # .data is created in the following
	rejection <- wilcox.test(data, alternative='two.sided')$p.value < alpha
	if(!rejection)	return(NA)
	else{
		estimates <- FPF:::pointWiseMixtureFitFast(data, generateMixtureControl())
		return(	(estimates[['p3.1']]-true.prevalence)^2 )
	}
}




# Compute Kurtosis of Active Component
ScaleMixtureKurtosis<- function(p, sd1, sd2, standardized){ 
  # Fourth Momeht:
  fourth.moment<- 3*p*sd1^4 + 3*(1 - p)*sd2^4
  
  #Squared Variance:
  squared.variance<- ( p * sd1^2 + (1-p) * sd2^2 )^2
  
  # Fourth Standardized Moment:
  fourth.standardized.moment <- fourth.moment / squared.variance 
  
  if(standardized) return(fourth.standardized.moment)
  else return(fourth.moment)
  
}
##Testing:
#curve(ScaleMixtureKurtosis(0.5, 1, x, standardized=FALSE), 1, 10)



# Main workhorse. Note that for reasons of speed and accuracy, the same data is used for all configurations.
computeMSE <- function(configurations, replications){
	#replications <- 10
	result <- data.frame(configurations, replications, MSE=NA, MSE.SE=NA, RejectedMSE=NA, RejectedMSE.SE=NA, RejectedMSE2=NA, RejectedMSE2.SE=NA  )	
	result$kurtosis <- with(result, ScaleMixtureKurtosis(p3/(p3+p4), sd3s, sd4s, standardized=TRUE))
	  
	max.n <- max(configurations$sample.size)
	base.data <- matrix(rnorm(max.n*replications), nrow=replications, ncol=max.n, dimnames=list(replications=NULL, sample=NULL))
	
	for(i in seq(length=nrow(result))){
	  
		#i <- 101
		.parameters <- configurations[i,]
		.n <- .parameters[["sample.size"]]
		true.prevalence <- sum(.parameters[c("p3","p4")])
		
		# Thin n
		.base.data <- base.data[,sample(1:max.n, size=.n, replace=FALSE)]
    
		
		# Apply parameters to base.data:
		group.array <- array(
				sample(x=c(1,2,3,4), size=.n*replications, replace=TRUE, prob=.parameters[c("p1","p2","p3","p4")]),
				dim=dim(.base.data)
		)
		mean.array <- (group.array == 3 | group.array == 4) * .parameters$mu
		sd.array <-  as.numeric(group.array == 1) * .parameters$sd1s +  as.numeric(group.array == 2) * .parameters$sd2s + as.numeric(group.array == 3) * .parameters$sd3s + 	as.numeric(group.array == 4) * .parameters$sd4s
		.data <- .base.data * sd.array + mean.array
		
				
		# Compute squared eror per row
		.MSEs <- apply(.data, 1, getPrevalenceSquaredError, true.prevalence=true.prevalence)		
		.RejectedMSEs <- apply(.data, 1, getRejectedPrevalenceSquaredError, true.prevalence=true.prevalence)
		.RejectedMSEs2 <- apply(.data, 1, getRejectedPrevalenceSquaredError2, true.prevalence=true.prevalence)	
		# Write Result:
		result[i,'MSE'] <- mean(.MSEs, na.rm=TRUE)
		result[i,'MSE.SE'] <- sd(.MSEs, na.rm=TRUE)
    result[i,'RejectedMSE'] <- mean(.RejectedMSEs, na.rm=TRUE)
		result[i,'RejectedMSE.SE'] <- sd(.RejectedMSEs, na.rm=TRUE)
		result[i,'RejectedMSE2'] <- mean(.RejectedMSEs2, na.rm=TRUE)
		result[i,'RejectedMSE2.SE'] <- sd(.RejectedMSEs2, na.rm=TRUE)

	}
  return(result)
}
