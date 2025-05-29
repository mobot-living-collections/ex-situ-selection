
######################################################################################################
######################################################################################################
##
## TITLE
##
## Simulations described in Appendix S1 to explore properties of the variance in relative fitness
## among maternal lines sourced from region "r" (equation 3 in main text), under the null model
## for opportunity for selection, based on demographic stochasticity.  
##
## INTRODUCTION
##
## This script is part of a project on provenance selection and opportunity for selection in ex
## situ collections for plant conservation.
##
## The code below implements a simulation (described in detail in Appendix S1) that explores
## properties of the distribution of the variance in relative fitness among maternal lines sourced
## from region "r", a quantity described in equation 3 of the main text and hereafter referred to
## as "Irt". The simulation explores how, under the null model for opportunity for selection, the
## mean and variance of "Irt" are affected by survival probability, initial number of plants from
## region "r", the number of maternal lines and variation in the initial number of plants per
## maternal line.
##
## DATA FILES REQUIRED
##
## None
##
######################################################################################################
######################################################################################################

#create dataframe with values of survival probability ("overall.survival"), number of maternal lines
#from region "r" ("nr"), the initial number of plants from region "r" ("Nr") and whether the
#distribution of the initial number of plants from region "r" is evenly distributed across maternal
#lines or not ("even.Nri.0"). 
OS.par <- expand.grid(overall.survival = seq(0.1, 0.9, 0.1), nr = c(10, 50), even.Nri.0=c(T,F), Nr = c(100, 200, 400))

#create vectors to store simulated values of the mean survival rate of individual plants across
#maternal lines sourced from a given region "r", "Srt.hat", the mean and variance of "Irt" and th 95%
#confidence interval for "Irt".
mean.Srt.hat <- rep(NA, times=nrow(OS.par))
mean.Irt <- rep(NA, times=nrow(OS.par))
var.Irt <- rep(NA, times=nrow(OS.par))
length.95CI.Irt <- rep(NA, times=nrow(OS.par))

#define number of simulation iterations
k <- 10000 #iterations

#begin simulation
for(i in 1:nrow(OS.par)){
	Srt.hat <- rep(NA, times=k) #mean survival across maternal lines
	Irt <- rep(NA, times=k) #opportunity for selection
	for(j in 1:k){
		Nr <- OS.par$Nr[i]
		nr <- OS.par$nr[i]		
		overall.survival <- OS.par$overall.survival[i] #overall survival probability: range from 0.1 to 0.9 every 0.1
		Kr <-  Nr*overall.survival
		#Kr
		
		#distribution of plants among maternal lines
		if(OS.par$even.Nri.0[i] == T) Nri.0 <- rep(Nr/nr, times=nr)
		if(OS.par$even.Nri.0[i] == F) Nri.0 <- rep(c(0.75*Nr, 0.25*Nr)/(nr*0.5), each=c(nr*0.5))
		#Nri.0
		#sum(Nri.0)
		#Nr

		survival <- sample(rep(c(0,1), times=c(Nr-Kr,Kr)))
		#sum(survival)
		#Kr
		#length(survival)
		#Nr

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
		Irt[j] <- var((Nri.t/Nri.0)/Srt.hat[j]) #opportunity for selection
	}
	mean.Srt.hat[i] <- mean(Srt.hat, na.rm=T)
	mean.Irt[i] <- mean(Irt, na.rm=T)
	var.Irt[i] <- var(Irt, na.rm=T)
	length.95CI.Irt[i] <- quantile(Irt, 0.975, na.rm=T) - quantile(Irt, 0.025, na.rm=T) 
}
#end simulation

#create vector with colors to plot results
nr.col <- hcl.colors(length(unique(OS.par$nr)), palette = "Viridis")
#nr.col <- nr.col[length(nr.col):1]
#hcl.pals()

#select a value for the initial number of plants from region "r" ("Nr"),
#the plots below will show the results for that value of "Nr"
Nr.to.plot <- 100 #select a value from unique(OS.par$Nr)

#plot results for the mean of "Srt.hat" (i.e., mean.Srt.hat)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(OS.par$overall.survival[OS.par$Nr==Nr.to.plot], mean.Srt.hat[OS.par$Nr==Nr.to.plot], xlim=c(0.1,0.9), ylim=c(0.1,0.9),
	bty="n", type="n", asp=1,
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("E[ ", italic(bar(S[r])(t)), " ]", sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(OS.par$nr))){
	points(OS.par$overall.survival[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==T & OS.par$nr==unique(OS.par$nr)[i]],
	mean.Srt.hat[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==T & OS.par$nr==unique(OS.par$nr)[i]], col=nr.col[i], type="o", pch=19)
}
for(i in 1:length(unique(OS.par$nr))){
	points(OS.par$overall.survival[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==F & OS.par$nr==unique(OS.par$nr)[i]],
	mean.Srt.hat[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==F & OS.par$nr==unique(OS.par$nr)[i]], col=nr.col[i], type="o", pch=17)
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
axis(2, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[r]) == .(Nr.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K[r])/italic(N[r]), sep="")), side=1,
	line=3.5, cex=1.5)
legend(0.1, 0.85, col=nr.col, lty=1, pch=rep(c(19, 17), each=3),
	cex=1.5, pt.cex=1,
	c(expression(paste(italic(n[r]), "=10", ", even")),
	expression(paste(italic(n[r]), "=50", ", even")),
	expression(paste(italic(n[r]), "=10", ", uneven")),
	expression(paste(italic(n[r]), "=50", ", uneven"))))

#plot results for the mean of "Irt" (i.e., mean.Irt)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(OS.par$overall.survival[OS.par$Nr==Nr.to.plot], mean.Irt[OS.par$Nr==Nr.to.plot], xlim=c(0.1,0.9), ylim=c(0,6.1),
	bty="n", type="n",
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("E [ ", italic(I[r](t)), " ]", sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(OS.par$nr))){
	points(OS.par$overall.survival[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==T & OS.par$nr==unique(OS.par$nr)[i]],
	mean.Irt[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==T & OS.par$nr==unique(OS.par$nr)[i]], col=nr.col[i], type="o", pch=19)
}
for(i in 1:length(unique(OS.par$nr))){
	points(OS.par$overall.survival[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==F & OS.par$nr==unique(OS.par$nr)[i]],
	mean.Irt[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==F & OS.par$nr==unique(OS.par$nr)[i]], col=nr.col[i], type="o", pch=17)
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[r]) == .(Nr.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K[r])/italic(N[r]), sep="")), side=1,
	line=3.5, cex=1.5)
#mtext(side=2, "a)", at=6.8, line=3, cex=1.5, las=2)
#mtext(side=2, "c)", at=6.8, line=3, cex=1.5, las=2)
mtext(side=2, "e)", at=6.8, line=3, cex=1.5, las=2)
legend(0.5, 5.5, col=nr.col, lty=1, pch=rep(c(19, 17), each=3),
	cex=1.5, pt.cex=1,
	c(expression(paste(italic(n[r]), "=10", ", even")),
	expression(paste(italic(n[r]), "=50", ", even")),
	expression(paste(italic(n[r]), "=10", ", uneven")),
	expression(paste(italic(n[r]), "=50", ", uneven"))))

#plot results for the variance of "Irt" (i.e., var.Irt)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(OS.par$overall.survival[OS.par$Nr==Nr.to.plot], var.Irt[OS.par$Nr==Nr.to.plot], xlim=c(0.1,0.9), ylim=c(0,0.4),
	bty="n", type="n",
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("Var [ ", italic(I[r](t)), " ]", sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(OS.par$nr))){
	points(OS.par$overall.survival[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==T & OS.par$nr==unique(OS.par$nr)[i]],
	var.Irt[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==T & OS.par$nr==unique(OS.par$nr)[i]], col=nr.col[i], type="o", pch=19)
}
for(i in 1:length(unique(OS.par$nr))){
	points(OS.par$overall.survival[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==F & OS.par$nr==unique(OS.par$nr)[i]],
	var.Irt[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==F & OS.par$nr==unique(OS.par$nr)[i]], col=nr.col[i], type="o", pch=17)
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[r]) == .(Nr.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K[r])/italic(N[r]), sep="")), side=1,
	line=3.5, cex=1.5)
mtext(side=2, "b)", at=0.45, line=3, cex=1.5, las=2)
#mtext(side=2, "d)", at=0.45, line=3, cex=1.5, las=2)
#mtext(side=2, "f)", at=0.45, line=3, cex=1.5, las=2)
legend(0.5, 5.5, col=nr.col, lty=1, pch=rep(c(19, 17), each=3),
	cex=1.5, pt.cex=1,
	c(expression(paste(italic(n[r]), "=10", ", even")),
	expression(paste(italic(n[r]), "=50", ", even")),
	expression(paste(italic(n[r]), "=20", ", uneven")),
	expression(paste(italic(n[r]), "=50", ", uneven"))))

#plot results for the lenght of the 95% CI for "Irt" (i.e., length.95CI.Irt)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5.5, 5, 4, 2) + 0.1)
plot(OS.par$overall.survival[OS.par$Nr==Nr.to.plot], length.95CI.Irt[OS.par$Nr==Nr.to.plot], xlim=c(0.1,0.9), ylim=c(0,2.2),
	bty="n", type="n",
	#xlab=expression(paste("Survival probability for plants from region ", italic(r), " (", italic(S[or]), ")", sep="")),
	xlab="",
	ylab=expression(paste("Length of 95% CI for ", italic(bar(I[r])), sep="")),
	cex.lab=1.5, cex.axis=1.5)
for(i in 1:length(unique(OS.par$nr))){
	points(OS.par$overall.survival[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==T & OS.par$nr==unique(OS.par$nr)[i]],
	length.95CI.Irt[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==T & OS.par$nr==unique(OS.par$nr)[i]], col=nr.col[i], type="o", pch=19)
}
for(i in 1:length(unique(OS.par$nr))){
	points(OS.par$overall.survival[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==F & OS.par$nr==unique(OS.par$nr)[i]],
	length.95CI.Irt[OS.par$Nr==Nr.to.plot & OS.par$even.Nri.0==F & OS.par$nr==unique(OS.par$nr)[i]], col=nr.col[i], type="o", pch=17)
}
axis(1, at=seq(0.1, 0.9, 0.1), cex.axis=1.5)
mtext(side=3, bquote(italic(N[r](0)) == .(Nr.to.plot)), cex=1.5)
mtext(expression(paste("Survival probability, ", italic(K[r])/italic(N[r]), sep="")), side=1,
	line=3.5, cex=1.5)
legend(0.5, 5.5, col=nr.col, lty=1, pch=rep(c(19, 17), each=3),
	cex=1.5, pt.cex=1,
	c(expression(paste(italic(n[r]), "=10", ", even")),
	expression(paste(italic(n[r]), "=20", ", even")),
	expression(paste(italic(n[r]), "=10", ", uneven")),
	expression(paste(italic(n[r]), "=50", ", uneven"))))





