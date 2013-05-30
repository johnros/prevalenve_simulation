#
# Note: Typically, should not be run as the data has already been created and loaded in PreparePaperArtUsingFPF2.R
#
###############################################################################




# Generate signal image
cols<- 256
rows<- 256
raw.signal<- matrix(0, ncol=256, nrow=256)
signal.constuctor<- function(x, y, x.center, y.center, radius.x, radius.y ) {
  as.numeric( ((x-x.center)/radius.x)^2 + ((y-y.center)/radius.y)^2 < 1 )
}
  
	
	raw.signal<- 
			outer(1:cols, 1:rows, FUN = signal.constuctor, x.center=cols*1/3, y.center=rows/2, radius.x=cols/8, radius.y=cols/8  ) + 
			outer(1:cols, 1:rows, FUN = signal.constuctor, x.center=cols*2/3, y.center=rows/2, radius.x=cols/8, radius.y=cols/8  )
	image(raw.signal)
	
	
	
	
	
signal.constructor4<- function(x, y, x.center, y.center, radius.x, radius.y, theta) {
		as.numeric(
				(((x-x.center)*cos(theta) + (y-y.center)*sin(theta)) / radius.x)^2  + (((x-x.center)*sin(theta) + (y-y.center)*cos(theta)) / radius.y)^2  < 1	
		)
		#(x-x.center)^2 + (y-y.center)^2 <= (radius.x*radius.y)^2 / ( (radius.y*cos(theta))^2 + (radius.x*sin(theta))^2 ) 
	}	
	
	raw.signal4<- 
			outer(1:cols, 1:rows, FUN = signal.constructor4, x.center=cols*1/3, y.center=rows/2, radius.x=cols/18, radius.y=cols/14, theta=1 ) + 
			outer(1:cols, 1:rows, FUN = signal.constructor4, x.center=cols*2/3, y.center=rows/2, radius.x=cols/10, radius.y=cols/18, theta=pi  )
	image(raw.signal4)
	
	
# Perturbate image
perturbations<- 100
image.perturbations<- array(dim = c(cols, rows, perturbations))
contour.perturbations<- list()
# Perturbation parameters:
sigma.center<- 5
sigma.radius.x<- 5
sigma.radius.y<- 5

for (i in 1:perturbations){
		#i<- 1
		x.center1 <- cols*1/3 + rnorm(1, sd = sigma.center)
		x.center2 <- cols*2/3 + rnorm(1, sd = sigma.center)
		y.center1<- rows/2 + rnorm(1, sd = sigma.center) 
		y.center2<- rows/2 + rnorm(1, sd = sigma.center)
		x.radius1<-cols/10+ rnorm(1, sd = sigma.radius.x)
		x.radius2<- cols/10+ rnorm(1, sd = sigma.radius.x*2)
		y.radius1<- cols/10 + rnorm(1, sd = sigma.radius.y)
		y.radius2<- cols/10 + rnorm(1, sd = sigma.radius.y*2)
		theta1<- runif(1, 0,pi/2)
		theta2<- 0
		
		image.perturbations[,,i]<-
				outer(1:cols, 1:rows, FUN = signal.constructor4, x.center= x.center1, y.center=y.center1, radius.x=x.radius1, radius.y=y.radius1, theta=theta1 ) + 
				outer(1:cols, 1:rows, FUN = signal.constructor4, x.center= x.center2, y.center=y.center2, radius.x=x.radius2, radius.y=y.radius2, theta=theta2 )
	}
	
	# Show percent overlap
	signal.prevalence<- apply(image.perturbations, c(1,2), mean)
require(lattice)
levelplot(signal.prevalence)
	
	
	
	# Estimate signal prevalence (new version):
	pointwise2<- function(y){
		clean.y<- na.omit(y)
		EMfit<- unlist(FPF:::pointWiseMixtureFitFast(clean.y, generateMixtureControl()))			
	}
	
	





#### Start generating estimate maps ####
sigma.noise<- 0.5



#### Add heavytailed noise to image ####
rmixednorm2<-function(draws,p1, p2, p3, mu, A, B, C){
  group.count<- rmultinom(1, draws, prob=c(p1,p2,p3))
  mus<- c(0,0,mu)
  sds<- sqrt(c(A,B,C))
  result<- c(rnorm(group.count[1], mus[1], sds[1]),rnorm(group.count[2], mus[2], sds[2]) ,rnorm(group.count[3], mus[3], sds[3]) )  
  return(result)
} 
## Testing:
#.p<- 0.4
#stem(rmixednorm2(draws = 1000, p1=(1-.p)/2, p2=(1-.p)/2 , p3=0.2 , mu= 5, A=0.5, B=0.7,  C=0.4))




noisy.prevalence2<- image.perturbations + array(
		rmixednorm2(draws = length(image.perturbations),p1=0.88, p2=0.12, p3=0, A = 0.15, B=1, C=0.25 , mu=0),
		dim=dim(image.perturbations))
# Warning: long run (several hours)
.result2 <- apply(noisy.prevalence2, c(1,2), pointwise2)


#### Add true model assumptions to image #### 
sigma.noise<- 0.5
noisy.prevalence3<- image.perturbations +  
		image.perturbations * array(rnorm(length(image.perturbations), sd=sigma.noise), dim=dim(image.perturbations))+
		(1-image.perturbations) * array(rmixednorm2(draws = length(image.perturbations),p1=0.88, p2=0.12, p3=0, A = 0.15, B=1, C=0.25 , mu=0),dim=dim(image.perturbations))
# Warning: long run (several hours)
.result3 <- apply(noisy.prevalence3, c(1,2), pointwise2)



