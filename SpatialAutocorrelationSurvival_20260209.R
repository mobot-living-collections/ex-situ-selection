
######################################################################################################
######################################################################################################
##
## TITLE
##
## Provenance selection and opportunity for selection in an ex situ collection of Quercus arkansana
## at MBG.
##
## INTRODUCTION
##
## This script is part of a project on provenance selection and opportunity for selection in ex
## situ collections for plant conservation. The project uses as a study system a living collection
## of Quercus arkansana (Fagaceae) maintained in the Plant Nursery of MBG.
##
## DATA FILES REQUIRED
##
## 1. "QaSurDat_2024Jun21_122839.csv" [Survival data]
##
## CONTENTS
##
## 1. Read and examine survival data.
## 2. Examine spatial autocorrelation in individual plant survival.
## 3. Examine spatial autocorrelation in deviance residuals from the logistic regression model
##    of individual plant survival as a function of maternal line.
##
######################################################################################################
######################################################################################################


#############################################################################################
# 1. Read and examine survival data.
#############################################################################################

#read and examine survival data
Qa.surv.data <- read.table("QaSurDat_2024Jun21_122839.csv", header = T, sep =",")
dim(Qa.surv.data)
colnames(Qa.surv.data)
str(Qa.surv.data)
head(Qa.surv.data)


#############################################################################################
# 2. Examine spatial autocorrelation in individual plant survival.
#############################################################################################

#calculate spatial autocorrelation function for survival data
autocor.model.1 <- acf(as.numeric(Qa.surv.data$Dead[Qa.surv.data$Position]), lag.max=30)
str(autocor.model.1)
autocor.model.1$acf[,1,1]

#create null distribution of the spatial autocorrelation function for survival data
max.lag <- 30
k.sim <- 5000
null.acf <- matrix(NA, nrow=k.sim, ncol=max.lag+1)
for(i in 1:k.sim){
	null.Dead <- sample(as.numeric(Qa.surv.data$Dead[Qa.surv.data$Position]))
	null.acf[i,] <- acf(null.Dead, plot=F, lag.max=30)$acf[1:(max.lag+1),1,1]
}

#obtain and examine p-values
p.values <- rep(NA, times = (max.lag+1))
for(j in 1:(max.lag+1)){
	p.values[j] <- min(sum(autocor.model.1$acf[j,1,1] >= null.acf[,j])/k.sim,
		sum(autocor.model.1$acf[j,1,1] <= null.acf[,j])/k.sim)
}
p.values
which(p.values <= 0.025)

#obtain 95% confidence intervals for the null distribution of the spatial autocorrelation function
#of the spatial autocorrelation function for survival data 
LL.95CI.model.1 <- apply(null.acf, MARGIN=2, FUN=function(x) quantile(x, probs=0.025))
UL.95CI.model.1 <- apply(null.acf, MARGIN=2, FUN=function(x) quantile(x, probs=0.975))

#plot spatial autocorrelation function for survival data 
#and 95% confidence interval for the null model
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 5, 4, 2) + 0.1)
plot(autocor.model.1, ci=0,  ylim=c(-0.17, 1),
	type="o", pch=19,
	main="",
	xlab="Lag in plant position",
	ylab="Autocorrelation coefficient",
	cex.axis=1.5, cex.lab=1.5)
points(1:30, LL.95CI.model.1[2:31], type="l", lty=2, col="gray")
points(1:30, UL.95CI.model.1[2:31], type="l", lty=2, col="gray")
#mtext(side=2, "a)", cex=1.5, las=1, at=1.15, line=2.5)

#############################################################################################
# 3. Examine spatial autocorrelation in deviance residuals from the logistic regression model
#    of individual plant survival as a function of maternal line.
#############################################################################################

#Fit logistic regression model of individual plant survival as a function of  maternal line
model.2 <- glm(Qa.surv.data$Dead ~ Qa.surv.data$AccessionNumber, family = binomial(link="logit"))
summary(model.2)

#examine the distribution of deviance residuals
hist(residuals(model.2, type="deviance"), breaks=seq(-2.2, 2.2, 0.1))

#examine relationship between deviance residuals and sequential plant position
plot(Qa.surv.data$Position, residuals(model.2, type="deviance")[Qa.surv.data$Position])

#calculate the spatial autocorrelation function for deviance residuals
autocor.model.2 <- acf(residuals(model.2, type="deviance")[Qa.surv.data$Position], lag.max=30)

#obtain 95% confidence intervals for the null distribution of the spatial autocorrelation function
LL.95CI.model.2 <- qnorm(0.025)/sqrt(autocor.model.2$n.used) 
UL.95CI.model.2 <- qnorm(0.975)/sqrt(autocor.model.2$n.used) 

#plot the spatial autocorrelation function for deviance residuals
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 5, 4, 2) + 0.1)
plot(autocor.model.2, ci=0, ylim=c(-0.17, 1),
	type="o", pch=19, ci.col="gray",
	main="",
	xlab="Lag in plant position",
	ylab="Autocorrelation coefficient",
	cex.axis=1.5, cex.lab=1.5)
points(1:30, rep(LL.95CI.model.2, 30), type="l", lty=2, col="gray")
points(1:30, rep(UL.95CI.model.2, 30), type="l", lty=2, col="gray")
#mtext(side=2, "b)", cex=1.5, las=1, at=1.15, line=2.5)


