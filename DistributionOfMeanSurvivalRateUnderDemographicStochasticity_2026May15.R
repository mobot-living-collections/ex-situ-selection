
######################################################################################################
######################################################################################################
##
## TITLE
##
## Distribution of mean survival rate under the null model of demographic stochastocity.
##
## INTRODUCTION
##
## This script is part of the supplementary material for:
##
## Sucher et al. Measuring maternal line selection driven by differential survival in ex situ
## collections for plant conservation. Accepted in Conservation Biology.
##
## This code implements two simulations described in Appendix S1 of Sucher et al., designed to explore
## properties of the distribution of mean survival rate of individual plants across maternal lines
## sourced from a given region "r", a quantity described in equation 1 of the main text and
## hereafter referred to as "Srt.hat". The first simulation explores how, under the null model for
## provenance selection, the mean, variance and 95% confidence interval (CI) of "Srt.hat" are afected
## by survival probability, initial number of plants per maternal line and number of maternal lines.
## The second simulation focuses on how the mean, variance and 95% CI of "Srt.hat" are affected by
## variation in the initial number of plants per maternal line, also under the null model for
## provenance selection.
##
## DATA FILES REQUIRED
##
## None
##
## CONTENTS
##
## 1. First simulation: effects of survival probability, initial number of plants per maternal line
##    and number of maternal lines on the mean, variance and 95% confidence interval of "Srt.hat".
##
## 2. Second simulation: effect of variation in the initial number of plants per maternal line on
##    the mean, variance and 95% confidence interval of "Srt.hat".
##
######################################################################################################
######################################################################################################


######################################################################################################
# 1. First simulation: effects of survival probability, initial number of plants per maternal line
#    and number of maternal lines on the mean, variance and 95% confidence interval of "Srt.hat".
######################################################################################################

#create dataframe with values of survival probability ("overall.survival"), number of maternal lines
#from region "r" ("nr") and the initial number of plants per maternal line from region "r" ("Nri.0").  
PS.parA <- expand.grid(overall.survival = seq(0.1, 0.9, 0.1), nr = c(10, 20), Nri.0=c(2, 10, 20))

#create vectors to store simulated values of mean, variance and 95% confidence interval of "Srt.hat"
mean.Srt.hat <- rep(NA, times=nrow(PS.parA))
var.Srt.hat <- rep(NA, times=nrow(PS.parA))
length.95CI.Srt.hat <- rep(NA, times=nrow(PS.parA))

#define the initial number of plants in the whole ex situ collection (N) and the total number of
#simulation iterations (k)
N <- 800 #initial number of plants in the whole ex situ collection
k <- 10000 #simulation iterations

#begin simulation
for(i in 1:nrow(PS.parA)){
	Srt.hat <- rep(NA, times=k) #mean survival across maternal lines
	for(j in 1:k){
		Nr <- PS.parA$nr[i] * PS.parA$Nri[i] 
		nr <- PS.parA$nr[i] #10, 20		
		overall.survival <- PS.parA$overall.survival[i] #overall survival probability: range from 0.1 to 0.9 every 0.1
		K <-  N*overall.survival
		#K
		
		Nri.0 <- rep(PS.parA$Nri[i], times=nr)
		#Nri.0
		#sum(Nri.0)
		#Nr

		survival <- sample(rep(c(0,1), times=c(N-K,K)), size=Nr)
		#sum(survival)
		#K
		#length(survival)
		#N

		#maternal line index, lower bound
		mli.L <- c(1,(1+cumsum(Nri.0)[-nr]))
		#maternal line index, upper bound 
		mli.U <- cumsum(Nri.0)

		Nri.t <- rep(NA, times=nr)
		for(w in 1:nr){
			Nri.t[w] <- sum(survival[mli.L[w]:mli.U[w]])
		}
		#Nri.t
		#sum(Nri.t)
		#length(Nri.t)

		Srt.hat[j] <- mean(Nri.t/Nri.0) #mean survival across maternal lines
	}
	mean.Srt.hat[i] <- mean(Srt.hat, na.rm=T)
	var.Srt.hat[i] <- var(Srt.hat, na.rm=T)
	length.95CI.Srt.hat[i] <- quantile(Srt.hat, 0.975, na.rm=T) - quantile(Srt.hat, 0.025, na.rm=T) 
}
#end simulation

#create vector with colors to plot results
#hcl.pals()
#nr.col <- hcl.colors(length(unique(PS.parA$nr)), palette = "Viridis")
#nr.col <- nr.col[length(nr.col):1]
#Oct 2, 2025 fixing the yellow color
#nr.col <- viridis::viridis(10, option = "D")[c(3,6)]
nr.col <- c("#D81B60", "#1E88E5")

#create vector with symbol size to plot results
nr.point.cex <- c(1.5, 1)
nr.pch <- c(21, 19)

#select a value for the initial number of plants per maternal line ("Nri.0"),
#the plots below will show the results for that value of "Nri.0"
Nri.0.to.plot <- 2 #select a value from unique(PS.parA$Nri.0)
Nri.0.to.plot <- 10 #select a value from unique(PS.parA$Nri.0)
Nri.0.to.plot <- 20 #select a value from unique(PS.parA$Nri.0)

#plot results for the mean of "Srt.hat" (i.e., mean.Srt.hat)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(PS.parA$overall.survival[PS.parA$Nri.0==Nri.0.to.plot], mean.Srt.hat[PS.parA$Nri.0==Nri.0.to.plot], xlim=c(0.1,0.9), ylim=c(0.1,0.9),
	bty="n", type="n", asp=1,
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("E[ ", italic(bar(S[r])(t)), " ]", sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(PS.parA$nr))){
	points(PS.parA$overall.survival[PS.parA$Nri.0==Nri.0.to.plot & PS.parA$nr==unique(PS.parA$nr)[i]],
	mean.Srt.hat[PS.parA$Nri.0==Nri.0.to.plot & PS.parA$nr==unique(PS.parA$nr)[i]], col=nr.col[i], type="o", pch=nr.pch[i],
	cex=nr.point.cex[i])
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
axis(2, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[ri](0)) == .(Nri.0.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K)/italic(N), sep="")), side=1,
	line=3.5, cex=1.5)
#mtext(side=2, "a)", at=1, line=3, cex=1.5, las=2)
#mtext(side=2, "c)", at=1, line=3, cex=1.5, las=2)
#mtext(side=2, "e)", at=1, line=3, cex=1.5, las=2)
legend(0.1, 0.85, col=nr.col, lty=1, pch=nr.pch,
	cex=1.5, pt.cex=nr.point.cex,
	c(expression(paste(italic(n[r]), "=10")),
	expression(paste(italic(n[r]), "=20"))))

#plot results for the variance of "Srt.hat" (i.e., var.Srt.hat)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(PS.parA$overall.survival[PS.parA$Nri.0==Nri.0.to.plot], var.Srt.hat[PS.parA$Nri.0==Nri.0.to.plot], xlim=c(0.1,0.9), ylim=c(0,0.013),
	bty="n", type="n",
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("Var [ ", italic(bar(S[r])(t)), " ]", sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(PS.parA$nr))){
	points(PS.parA$overall.survival[PS.parA$Nri.0==Nri.0.to.plot & PS.parA$nr==unique(PS.parA$nr)[i]],
	var.Srt.hat[PS.parA$Nri.0==Nri.0.to.plot & PS.parA$nr==unique(PS.parA$nr)[i]], col=nr.col[i], type="o", pch=nr.pch[i],
	cex=nr.point.cex[i])
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[ri](0)) == .(Nri.0.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K)/italic(N), sep="")), side=1,
	line=3.5, cex=1.5)
#mtext(side=2, "b)", at=0.0145, line=3, cex=1.5, las=2)
#mtext(side=2, "d)", at=0.0145, line=3, cex=1.5, las=2)
#mtext(side=2, "f)", at=0.0145, line=3, cex=1.5, las=2)

#plot results for the lenght of the 95% confidence interval for "Srt.hat" (i.e., length.95CI.Srt.hat)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(PS.parA$overall.survival[PS.parA$Nri==Nri.0.to.plot], length.95CI.Srt.hat[PS.parA$Nri==Nri.0.to.plot], xlim=c(0.1,0.9), ylim=c(0,1),
	bty="n", type="n",
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("Length of 95% CI for ", italic(bar(I[r])), sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(PS.parA$nr))){
	points(PS.parA$overall.survival[PS.parA$Nri==Nri.0.to.plot & PS.parA$nr==unique(PS.parA$nr)[i]],
	length.95CI.Srt.hat[PS.parA$Nri==Nri.0.to.plot & PS.parA$nr==unique(PS.parA$nr)[i]], col=nr.col[i], type="o", pch=nr.pch[i],
	cex=nr.point.cex[i])
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[ri](0)) == .(Nri.0.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K[r])/italic(N[r]), sep="")), side=1,
	line=3.5, cex=1.5)
legend(0.1, 0.99, col=nr.col, lty=1, pch=nr.pch,
	cex=1.5, pt.cex=nr.point.cex,
	c(expression(paste(italic(n[r]), "=10")),
	expression(paste(italic(n[r]), "=20"))))


######################################################################################################
# 2. Second simulation: effect of variation in the initial number of plants per maternal line on the
#    mean, variance and 95% confidence interval of "Srt.hat".
######################################################################################################

#create dataframe with values of survival probability ("overall.survival"), number of maternal
#lines from region "r" ("nr"), the initial number of plants from region "r" represented in the
#ex situ collection ("Nr") and whether the distribution of the initial number of plants from
#region "r" is evenly distributed across maternal lines or not ("even.Nri.0"). 
PS.parB <- expand.grid(overall.survival = seq(0.1, 0.9, 0.1), nr = c(10, 50), even.Nri.0=c(T,F), Nr = c(100, 200, 400))

#create vectors to store simulated values of the mean, variance and 95% confidence interval of "Srt.hat"
mean.Srt.hat <- rep(NA, times=nrow(PS.parB))
var.Srt.hat <- rep(NA, times=nrow(PS.parB))
length.95CI.Srt.hat <- rep(NA, times=nrow(PS.parB))

#define the initial number of plants in the whole ex situ collection (N) and the total number of
#simulation iterations (k)
N <- 800 #initial number of plants in the whole ex situ collection
k <- 10000 #simulation iterations

#begin simulation
for(i in 1:nrow(PS.parB)){
	Srt.hat <- rep(NA, times=k) #mean survival across maternal lines
	for(j in 1:k){
		Nr <- PS.parB$Nr[i] #200, 400, 100
		nr <- PS.parB$nr[i] #10, 50		
		overall.survival <- PS.parB$overall.survival[i] #overall survival probability: range from 0.1 to 0.9 every 0.1
		K <-  N*overall.survival
		#K
		
		#distribution of plants among maternal lines
		if(PS.parB$even.Nri.0[i] == T) Nri.0 <- rep(Nr/nr, times=nr)
		if(PS.parB$even.Nri.0[i] == F) Nri.0 <- rep(c(0.75*Nr, 0.25*Nr)/(nr*0.5), each=c(nr*0.5))
		#Nri.0
		#sum(Nri.0)
		#Nr

		survival <- sample(rep(c(0,1), times=c(N-K,K)), size=Nr)
		#sum(survival)
		#K
		#length(survival)
		#N

		#maternal line index, lower bound
		mli.L <- c(1,(1+cumsum(Nri.0)[-nr]))
		#maternal line index, upper bound 
		mli.U <- cumsum(Nri.0)

		Nri.t <- rep(NA, times=nr)
		for(w in 1:nr){
			Nri.t[w] <- sum(survival[mli.L[w]:mli.U[w]])
		}
		#Nri.t
		#sum(Nri.t)
		#length(Nri.t)

		Srt.hat[j] <- mean(Nri.t/Nri.0) #mean survival across maternal lines
	}
	mean.Srt.hat[i] <- mean(Srt.hat, na.rm=T)
	var.Srt.hat[i] <- var(Srt.hat, na.rm=T)
	length.95CI.Srt.hat[i] <- quantile(Srt.hat, 0.975, na.rm=T) - quantile(Srt.hat, 0.025, na.rm=T)
}
#end simulation

#create vector with colors to plot results
#hcl.pals()
#nr.col <- hcl.colors(length(unique(PS.parB$nr)), palette = "Viridis")
#nr.col <- nr.col[length(nr.col):1]
#Oct 2, 2025 fixing the yellow color
#nr.col <- viridis::viridis(10, option = "D")[c(3,6)]
nr.col <- c("#D81B60", "#1E88E5")

#create vector with symbol type and size to plot results
nr.point.cex <- c(1.5, 1)
nr.pch.even <- c(21, 19)
nr.pch.uneven <- c(2, 17)

#select a value for the initial number of plants from region r ("Nr") in the ex situ collection,
#the plots below will show the results for that value of "Nr"
Nr.to.plot <- 100 #select a value from unique(PS.parB$Nr)
Nr.to.plot <- 200 #select a value from unique(PS.parB$Nr)
Nr.to.plot <- 400 #select a value from unique(PS.parB$Nr)

#plot results for the mean of "Srt.hat" (i.e., mean.Srt.hat)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot], mean.Srt.hat[PS.parB$Nr==Nr.to.plot], xlim=c(0.1,0.9), ylim=c(0.1,0.9),
	bty="n", type="n", asp=1,
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("E[ ", italic(bar(S[r])(t)), " ]", sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(PS.parB$nr))){
	points(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==T & PS.parB$nr==unique(PS.parB$nr)[i]],
	mean.Srt.hat[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==T & PS.parB$nr==unique(PS.parB$nr)[i]], col=nr.col[i], type="o", pch=nr.pch.even[i],
	cex=nr.point.cex[i])
}
for(i in 1:length(unique(PS.parB$nr))){
	points(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==F & PS.parB$nr==unique(PS.parB$nr)[i]],
	mean.Srt.hat[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==F & PS.parB$nr==unique(PS.parB$nr)[i]], col=nr.col[i], type="o", pch=c(2,17),
	cex=nr.point.cex[i])
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
axis(2, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[r]) == .(Nr.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K)/italic(N), sep="")), side=1,
	line=3.5, cex=1.5)
#mtext(side=2, "a)", at=1, line=3, cex=1.5, las=2)
#mtext(side=2, "c)", at=1, line=3, cex=1.5, las=2)
#mtext(side=2, "e)", at=1, line=3, cex=1.5, las=2)
legend(0.1, 0.88, col=nr.col, lty=1, pch=c(21, 19, 2, 17),
	cex=1.5, pt.cex=rep(c(1.5, 1), times=2),
	c(expression(paste(italic(n[r]), "=10", ", even")),
	expression(paste(italic(n[r]), "=50", ", even")),
	expression(paste(italic(n[r]), "=10", ", uneven")),
	expression(paste(italic(n[r]), "=50", ", uneven"))))

##plot results for the variance of "Srt.hat" (i.e., var.Srt.hat)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot], var.Srt.hat[PS.parB$Nr==Nr.to.plot], xlim=c(0.1,0.9), ylim=c(0,0.0032),
	bty="n", type="n",
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("Var [ ", italic(bar(S[r])(t)), " ]", sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(PS.parB$nr))){
	points(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==T & PS.parB$nr==unique(PS.parB$nr)[i]],
	var.Srt.hat[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==T & PS.parB$nr==unique(PS.parB$nr)[i]], col=nr.col[i], type="o", pch=nr.pch.even[i],
	cex=nr.point.cex[i])
}
for(i in 1:length(unique(PS.parB$nr))){
	points(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==F & PS.parB$nr==unique(PS.parB$nr)[i]],
	var.Srt.hat[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==F & PS.parB$nr==unique(PS.parB$nr)[i]], col=nr.col[i], type="o", pch=nr.pch.uneven[i],
	cex=nr.point.cex[i])
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[r]) == .(Nr.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K)/italic(N), sep="")), side=1,
	line=3.5, cex=1.5)
#mtext(side=2, "b)", at=0.0035, line=3, cex=1.5, las=2)
#mtext(side=2, "d)", at=0.0035, line=3, cex=1.5, las=2)
#mtext(side=2, "f)", at=0.0035, line=3, cex=1.5, las=2)

#plot results for the lenght of the 95% CI for "Irt" (i.e., length.95CI.Srt.hat)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot], length.95CI.Srt.hat[PS.parB$Nr==Nr.to.plot], xlim=c(0.1,0.9), ylim=c(0,1),
	bty="n", type="n",
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("Length of 95% CI for ", italic(bar(I[r])), sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(PS.parB$nr))){
	points(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==T & PS.parB$nr==unique(PS.parB$nr)[i]],
	length.95CI.Srt.hat[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==T & PS.parB$nr==unique(PS.parB$nr)[i]], col=nr.col[i], type="o", pch=nr.pch.even[i],
	cex=nr.point.cex[i])
}
for(i in 1:length(unique(PS.parB$nr))){
	points(PS.parB$overall.survival[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==F & PS.parB$nr==unique(PS.parB$nr)[i]],
	length.95CI.Srt.hat[PS.parB$Nr==Nr.to.plot & PS.parB$even.Nri.0==F & PS.parB$nr==unique(PS.parB$nr)[i]], col=nr.col[i], type="o", pch=nr.pch.uneven[i],
	cex=nr.point.cex[i])
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[r]) == .(Nr.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K[r])/italic(N[r]), sep="")), side=1,
	line=3.5, cex=1.5)
legend(0.1, 0.99, col=nr.col, lty=1, pch=c(21, 19, 2, 17),
	cex=1.5, pt.cex=rep(c(1.5, 1), times=2),
	c(expression(paste(italic(n[r]), "=10", ", even")),
	expression(paste(italic(n[r]), "=50", ", even")),
	expression(paste(italic(n[r]), "=10", ", uneven")),
	expression(paste(italic(n[r]), "=50", ", uneven"))))


