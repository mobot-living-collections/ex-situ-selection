
######################################################################################################
######################################################################################################
##
## TITLE
##
## Provenance selection and opportunity for selection in an ex situ collection of Quercus arkansana
## at the Missouri Botanical Garden.
##
## INTRODUCTION
##
## This script is part of a project on provenance selection and opportunity for selection in ex
## situ collections for plant conservation. The project uses as a study system a living collection
## of Quercus arkansana (Fagaceae) maintained in the Oertli Hardy Plant Nursery of the
## Missouri Botanical Garden.
##
## DATA FILES REQUIRED
##
## 1. "QaSurDat_2024Jun21_122839.csv" [Survival data]
## 2. "QaAccCoor_2024Jul24_104558.csv" [Provenance data]
##
## CONTENTS
##
## 1. Read and examine data files.
## 2. Calculate selection coefficients for provenances (i.e., geographic regions).
## 3. Relationship between the number of maternal lines and the selection coefficients
#     across provenances (i.e., geographic regions).
## 4. Calculate opportunity for selection within each provenance (i.e., geographic region).
##
######################################################################################################
######################################################################################################


#############################################################################################
# 1. Read and examine data files.
#############################################################################################

#set working directory

#read and examine survival data
Qa.surv.data <- read.table("QaSurDat_2024Jun21_122839.csv", header = T, sep =",")
dim(Qa.surv.data) 
colnames(Qa.surv.data)
str(Qa.surv.data)
head(Qa.surv.data)

#read and examine provenance data
Qa.prov.data <- read.table("QaAccCoor_2024Jul24_104558.csv", header = T, sep =",")
dim(Qa.prov.data) 
colnames(Qa.prov.data)
str(Qa.prov.data)
head(Qa.prov.data)
 
#make sure each accession number in the survival data is also in the provenance data
match(Qa.surv.data$AccessionNumber, Qa.prov.data$AccessionNumber)
sum(is.na(match(Qa.surv.data$AccessionNumber, Qa.prov.data$AccessionNumber)))

#make sure each accession number in the provenance data is also in the survival data
match(Qa.prov.data$AccessionNumber, Qa.surv.data$AccessionNumber)
sum(is.na(match(Qa.prov.data$AccessionNumber, Qa.surv.data$AccessionNumber)))


#############################################################################################
# 2. Calculate selection coefficients for provenances (i.e., geographic regions).
#############################################################################################

#############################################################################################
# 2.1. Observed selection coefficients for provenances (i.e., geographic regions).

#count plants per maternal line (i.e., accession)
plants.per.accession <- tapply(Qa.surv.data$Dead, Qa.surv.data$AccessionNumber, FUN=length)
plants.per.accession
length(plants.per.accession)
class(plants.per.accession)
str(plants.per.accession)
barplot(table(plants.per.accession), ylab="Maternal lines", xlab="Plants")

#count deaths per maternal line (i.e., accession)
deaths.per.accession <- tapply(Qa.surv.data$Dead, Qa.surv.data$AccessionNumber, FUN=sum)
deaths.per.accession
length(deaths.per.accession)
class(deaths.per.accession)
str(deaths.per.accession)
barplot(table(deaths.per.accession), ylab="Maternal lines", xlab="Deaths")

#calculate survival rate per maternal line (i.e., accession)
survival.rate.per.accession <- deaths.per.accession / plants.per.accession
survival.rate.per.accession
length(survival.rate.per.accession)
class(survival.rate.per.accession)
str(survival.rate.per.accession)
hist(survival.rate.per.accession, breaks=seq(0, 1, 0.01))

#count maternal lines per provenance (i.e., geographic region)
identical(names(survival.rate.per.accession), Qa.prov.data$AccessionNumber) #make sure accessions are ordered in the same way
maternal.lines.per.provenance <- tapply(survival.rate.per.accession, Qa.prov.data$GeographicRegion, FUN=length) 
maternal.lines.per.provenance
length(maternal.lines.per.provenance)
class(maternal.lines.per.provenance)
str(maternal.lines.per.provenance)

#count plants per provenance (i.e., geographic region)
identical(names(plants.per.accession), Qa.prov.data$AccessionNumber) #make sure accessions are ordered in the same way
plants.per.provenance <- tapply(plants.per.accession, Qa.prov.data$GeographicRegion, FUN=sum) 
plants.per.provenance
length(plants.per.provenance)
class(plants.per.provenance)
str(plants.per.provenance)

#examine the number of plants per accession in each provenance (i.e., geographic region)
plants.per.accession[Qa.prov.data$GeographicRegion==2]
plants.per.accession[Qa.prov.data$GeographicRegion==3]
plants.per.accession[Qa.prov.data$GeographicRegion==4]
plants.per.accession[Qa.prov.data$GeographicRegion==5]
plants.per.accession[Qa.prov.data$GeographicRegion==8]

#calculate average survival rate per provenance (i.e., geographic region)
identical(names(survival.rate.per.accession), Qa.prov.data$AccessionNumber) #make sure accessions are ordered in the same way
mean.survival.rate.per.provenance <- tapply(survival.rate.per.accession, Qa.prov.data$GeographicRegion, FUN=mean) 
mean.survival.rate.per.provenance
length(mean.survival.rate.per.provenance)
class(mean.survival.rate.per.provenance)
str(mean.survival.rate.per.provenance)

#Calculate selection coefficients for provenances (i.e., geographic regions)
#region 2
(mean.survival.rate.per.provenance["2"]/mean.survival.rate.per.provenance["2"]) - 1
#region 3
(mean.survival.rate.per.provenance["3"]/mean.survival.rate.per.provenance["2"]) - 1
#region 4
(mean.survival.rate.per.provenance["4"]/mean.survival.rate.per.provenance["2"]) - 1
#region 5
(mean.survival.rate.per.provenance["5"]/mean.survival.rate.per.provenance["2"]) - 1
#region 8
(mean.survival.rate.per.provenance["8"]/mean.survival.rate.per.provenance["2"]) - 1

observed.selection.coefficients <- (mean.survival.rate.per.provenance/mean.survival.rate.per.provenance["2"]) - 1

#############################################################################################
# 2.2. Null distribution of selection coefficients for provenances (i.e., geographic regions).

#set up null model
k <- 10000 #null model iterations
null.selection.coefficients <- matrix(NA, nrow=k, ncol=5) #matrix to store null selection coefficients
colnames(null.selection.coefficients) <- sort(unique(Qa.prov.data$GeographicRegion))
dim(null.selection.coefficients) #matrix dimensions should equal the result of running the next line of code 
c(k,length(unique(Qa.prov.data$GeographicRegion)))
head(null.selection.coefficients)

#run null model
for(i in 1:k){
	null.dead <- sample(x=Qa.surv.data$Dead, size=length(Qa.surv.data$Dead))
	null.deaths.per.accession <- tapply(null.dead, Qa.surv.data$AccessionNumber, FUN=sum)
	null.survival.rate.per.accession <- null.deaths.per.accession / plants.per.accession
	null.mean.survival.rate.per.provenance <- tapply(null.survival.rate.per.accession, Qa.prov.data$GeographicRegion, FUN=mean)
	null.reference.region <- names(which.max(plants.per.provenance[null.mean.survival.rate.per.provenance>0]))
	null.selection.coefficients[i,1] <- (null.mean.survival.rate.per.provenance["2"]/null.mean.survival.rate.per.provenance[null.reference.region]) - 1
	null.selection.coefficients[i,2] <- (null.mean.survival.rate.per.provenance["3"]/null.mean.survival.rate.per.provenance[null.reference.region]) - 1
	null.selection.coefficients[i,3] <- (null.mean.survival.rate.per.provenance["4"]/null.mean.survival.rate.per.provenance[null.reference.region]) - 1
	null.selection.coefficients[i,4] <- (null.mean.survival.rate.per.provenance["5"]/null.mean.survival.rate.per.provenance[null.reference.region]) - 1
	null.selection.coefficients[i,5] <- (null.mean.survival.rate.per.provenance["8"]/null.mean.survival.rate.per.provenance[null.reference.region]) - 1
}

#examine null model results
summary(null.selection.coefficients)


#############################################################################################
# 3. Relationship between the number of maternal lines and the selection coefficients
#    across provenances (i.e., geographic regions).
#############################################################################################

#calculate median of null selection coefficients
null.median.selection.coefficients <- apply(X=null.selection.coefficients, MARGIN=2, FUN=quantile, probs=0.5)
#calculate mean of null selection coefficients
null.mean.selection.coefficients <- apply(X=null.selection.coefficients, MARGIN=2, FUN=mean)
#calculate lower 0.025 quantile of null selection coefficients
null.LL.selection.coefficients <- apply(X=null.selection.coefficients, MARGIN=2, FUN=quantile, probs=0.025)
#calculate upper 0.975 quantile of null selection coefficients
null.UL.selection.coefficients <- apply(X=null.selection.coefficients, MARGIN=2, FUN=quantile, probs=0.975)

#plot number of maternal lines for each provenance ((i.e., geographic region)in the horizontal axis
#and the respective observed and null selection coefficients in the vertical axis 
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 5, 4, 2) + 0.1)
#examine range of values in the horizontal axis
range(maternal.lines.per.provenance)
#examine range of values in the vertical axis
range(observed.selection.coefficients, null.LL.selection.coefficients, null.UL.selection.coefficients)
#set up plot
plot(maternal.lines.per.provenance, null.mean.selection.coefficients, xlim=c(0, 16), ylim=c(-1, 1.2),
	bty="n", cex.axis=1.5, cex.lab=1.5, xlab="Maternal lines", ylab=expression(paste("Selection coefficient (", italic(S[r]), " )", sep="")), type="n")
arrows(0, 0, 15, 0, code=0, lty=2, col="gray80", lwd=2) #add horizontal line marking selection coefficients equal to zero
#graph mean of null selection coefficients
points(maternal.lines.per.provenance, null.mean.selection.coefficients, pch=19, col="gray", cex=1)
#graph median of null selection coefficients
points(maternal.lines.per.provenance, null.median.selection.coefficients, pch=21, col="gray", cex=2.3)
#graph 95% confidence interval of null selection coefficients
arrows(maternal.lines.per.provenance, null.LL.selection.coefficients, maternal.lines.per.provenance, null.UL.selection.coefficients,
	length=0.1, angle=90, code=3, col="gray")
#graph observed selection coefficients
points(maternal.lines.per.provenance, observed.selection.coefficients, pch=19, col="black", cex=1.5)

#add simple version of legend
#legend("topright", c("Observed", "Null model mean", "Null model median", "Null model 95% CI"), pch=c(19,19,21, NA),
#	pt.cex=c(1.5, 1, 2.3, NA), col=c("black", "gray", "gray", "gray"), lty=c(NA, NA, NA, 1))

#add fidgety version of legend
legend("topright", c("Observed", "Null model mean", "Null model median", "Null model 95% CI"), pch=c(19,19,21,NA),
	pt.cex=c(1.5, 1, 2.3, 1), col=c("black", "gray", "gray", "gray"))
arrows(11.8, 0.88, 11.8, 0.945, length=0.1, angle=90, code=3, col="gray")

#label regions
text(15.6, 0.01, "2", cex=1.5, col="#D81B60")
text(3.5, 0.6, "3", cex=1.5, col="#1E88E5")
text(8.5, 0.1, "4", cex=1.5, col="#FFC107")
text(12.5, -0.1, "5", cex=1.5, col="#42FD7F")
text(0.5, -0.62, "8", cex=1.5, col="#E7A0FD")


#############################################################################################
# 4. Calculate opportunity for selection within each provenance (i.e., geographic region).
#############################################################################################

#############################################################################################
# 4.1. Observed opportunity for selection in each provenance (i.e., geographic region).

#calculate the standard deviation of survival rate per provenance (i.e., geographic region)
identical(names(survival.rate.per.accession), Qa.prov.data$AccessionNumber) #make sure accessions are ordered in the same way
standard.deviation.survival.rate.per.provenance <- tapply(survival.rate.per.accession, Qa.prov.data$GeographicRegion, FUN=sd)
standard.deviation.survival.rate.per.provenance

#Calculate opportunity for selection coefficients for provenances (i.e., geographic regions)
#region 2
(standard.deviation.survival.rate.per.provenance["2"]/mean.survival.rate.per.provenance["2"])^2
#region 3
(standard.deviation.survival.rate.per.provenance["3"]/mean.survival.rate.per.provenance["3"])^2
#region 4
(standard.deviation.survival.rate.per.provenance["4"]/mean.survival.rate.per.provenance["4"])^2
#region 5
(standard.deviation.survival.rate.per.provenance["5"]/mean.survival.rate.per.provenance["5"])^2
#region 8
(standard.deviation.survival.rate.per.provenance["8"]/mean.survival.rate.per.provenance["8"])^2

#############################################################################################
# 4.2. Compare null and observed opportunity for selection for region 2.

#set up null model for region 2
k <- 10000 #null model iterations
null.opportunity.for.selection.region2 <- rep(NA, times=k) #vector to store null opportunity for selection values
#extract survival data for region 2
index.region2 <- !is.na(match(Qa.surv.data$AccessionNumber, Qa.prov.data$AccessionNumber[Qa.prov.data$GeographicRegion==2]))
Qa.surv.data.region2 <- Qa.surv.data[index.region2,]
#examine result
dim(Qa.surv.data.region2)
head(Qa.surv.data.region2)
#the next two lines of code should produce the information (although in objects of different class)
plants.per.accession[Qa.prov.data$GeographicRegion==2]
table(Qa.surv.data.region2$AccessionNumber)
#check the match
identical(names(plants.per.accession[Qa.prov.data$GeographicRegion==2]), names(table(Qa.surv.data.region2$AccessionNumber)))
identical(as.vector(plants.per.accession[Qa.prov.data$GeographicRegion==2]), as.vector(table(Qa.surv.data.region2$AccessionNumber)))

#run null model
for(i in 1:k){
	null.dead.region2 <- sample(x=Qa.surv.data.region2$Dead, size=length(Qa.surv.data.region2$Dead))
	null.deaths.per.accession.region2 <- tapply(null.dead.region2, Qa.surv.data.region2$AccessionNumber, FUN=sum)
	null.survival.rate.per.accession.region2 <- null.deaths.per.accession.region2 / plants.per.accession[Qa.prov.data$GeographicRegion==2]
	null.mean.survival.rate.region2 <- mean(null.survival.rate.per.accession.region2)
	null.standard.deviation.survival.rate.region2 <- sd(null.survival.rate.per.accession.region2)
	null.opportunity.for.selection.region2[i] <- (null.standard.deviation.survival.rate.region2/null.mean.survival.rate.region2)^2
}

#examine null model results
summary(null.opportunity.for.selection.region2)

#define observed value of opportunity for selection
observed.opportunity.for.selection.region2 <- (standard.deviation.survival.rate.per.provenance["2"]/mean.survival.rate.per.provenance["2"])^2

#graphically contrast observed and null model results
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 5, 4, 2) + 0.1)
#examine range of values in the horizontal axis
range(observed.opportunity.for.selection.region2, null.opportunity.for.selection.region2) # 
#histogram with null values
hist(null.opportunity.for.selection.region2, breaks=seq(0, 0.6, 0.01),
	main="Region 2", col.main="black", xlab=expression(paste("Opportunity for selection (", italic(I[2]), ")", sep="")),
	ylab="Null model iterations", font.main=1, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,
	border="gray", col="gray90")

###########################################################################################################
# COLOR OPTIONS FOR REGIONS #
# IVAN'S COLOR SCHEME "#D81B60", "#1E88E5", "#FFC107", "#42FD7F", "#E7A0FD"
# VIRIDIS PACKAGE COLOR SCHEME OPTION "D" OR "#440154FF" "#3B528BFF" "#21908CFF" "#5DC863FF" "#FDE725FF"
###########################################################################################################

#add observed value
abline(v=observed.opportunity.for.selection.region2, col="#440154FF", lty=1, lwd=2)

#add legend
legend("topright", c("Observed", "Null model"), fill=c(NA, "gray90"), border=c(NA, "gray"))
#add symbol for observed value, a bit fidgety
arrows(0.504, 630, 0.504, 670, length=0, angle=90, col="#440154FF", lwd=2)

#calculate p-values
#lower tail
sum(null.opportunity.for.selection.region2 <= observed.opportunity.for.selection.region2)/length(null.opportunity.for.selection.region2)
#upper tail
sum(null.opportunity.for.selection.region2 >= observed.opportunity.for.selection.region2)/length(null.opportunity.for.selection.region2)

#############################################################################################
# 4.3. Compare null and observed opportunity for selection for region 3.

#set up null model for region 3
k <- 10000 #null model iterations
null.opportunity.for.selection.region3 <- rep(NA, times=k) #vector to store null opportunity for selection values
#extract survival data for region 3
index.region3 <- !is.na(match(Qa.surv.data$AccessionNumber, Qa.prov.data$AccessionNumber[Qa.prov.data$GeographicRegion==3]))
Qa.surv.data.region3 <- Qa.surv.data[index.region3,]
#examine result
dim(Qa.surv.data.region3)
head(Qa.surv.data.region3)
#the next two lines of code should produce the information (although in objects of different class)
plants.per.accession[Qa.prov.data$GeographicRegion==3]
table(Qa.surv.data.region3$AccessionNumber)
#check the match
identical(names(plants.per.accession[Qa.prov.data$GeographicRegion==3]), names(table(Qa.surv.data.region3$AccessionNumber)))
identical(as.vector(plants.per.accession[Qa.prov.data$GeographicRegion==3]), as.vector(table(Qa.surv.data.region3$AccessionNumber)))

#run null model
for(i in 1:k){
	null.dead.region3 <- sample(x=Qa.surv.data.region3$Dead, size=length(Qa.surv.data.region3$Dead))
	null.deaths.per.accession.region3 <- tapply(null.dead.region3, Qa.surv.data.region3$AccessionNumber, FUN=sum)
	null.survival.rate.per.accession.region3 <- null.deaths.per.accession.region3 / plants.per.accession[Qa.prov.data$GeographicRegion==3]
	null.mean.survival.rate.region3 <- mean(null.survival.rate.per.accession.region3)
	null.standard.deviation.survival.rate.region3 <- sd(null.survival.rate.per.accession.region3)
	null.opportunity.for.selection.region3[i] <- (null.standard.deviation.survival.rate.region3/null.mean.survival.rate.region3)^2
}

#examine null model results
summary(null.opportunity.for.selection.region3)

#define observed value of opportunity for selection
observed.opportunity.for.selection.region3 <- (standard.deviation.survival.rate.per.provenance["3"]/mean.survival.rate.per.provenance["3"])^2

#graphically contrast observed and null model results
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 5, 4, 2) + 0.1)
#examine range of values in the horizontal axis
range(observed.opportunity.for.selection.region3, null.opportunity.for.selection.region3) # 
#histogram with null values
hist(null.opportunity.for.selection.region3, breaks=seq(0, 4, 0.1),
	main="Region 3", col.main="black", xlab=expression(paste("Opportunity for selection (", italic(I[3]), ")", sep="")),
	ylab="Null model iterations", font.main=1, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,
	border="gray", col="gray90")
#add observed value
abline(v=observed.opportunity.for.selection.region3, col="#3B528BFF", lty=1, lwd=2)

#add legend
legend("topright", c("Observed", "Null model"), fill=c(NA, "gray90"), border=c(NA, "gray"))
#add symbol for observed value, a bit fidgety
arrows(3.37, 5000, 3.37, 5350, length=0, angle=90, col="#3B528BFF", lwd=2)

#calculate p-values
#lower tail
sum(null.opportunity.for.selection.region3 <= observed.opportunity.for.selection.region3)/length(null.opportunity.for.selection.region3)
#upper tail
sum(null.opportunity.for.selection.region3 >= observed.opportunity.for.selection.region3)/length(null.opportunity.for.selection.region3)

#############################################################################################
# 4.4. Compare null and observed opportunity for selection for region 4.

#set up null model for region 4
k <- 10000 #null model iterations
null.opportunity.for.selection.region4 <- rep(NA, times=k) #vector to store null opportunity for selection values
#extract survival data for region 4
index.region4 <- !is.na(match(Qa.surv.data$AccessionNumber, Qa.prov.data$AccessionNumber[Qa.prov.data$GeographicRegion==4]))
Qa.surv.data.region4 <- Qa.surv.data[index.region4,]
#examine result
dim(Qa.surv.data.region4)
head(Qa.surv.data.region4)
#the next two lines of code should produce the information (although in objects of different class)
plants.per.accession[Qa.prov.data$GeographicRegion==4]
table(Qa.surv.data.region4$AccessionNumber)
#check the match
identical(names(plants.per.accession[Qa.prov.data$GeographicRegion==4]), names(table(Qa.surv.data.region4$AccessionNumber)))
identical(as.vector(plants.per.accession[Qa.prov.data$GeographicRegion==4]), as.vector(table(Qa.surv.data.region4$AccessionNumber)))

#run null model
for(i in 1:k){
	null.dead.region4 <- sample(x=Qa.surv.data.region4$Dead, size=length(Qa.surv.data.region4$Dead))
	null.deaths.per.accession.region4 <- tapply(null.dead.region4, Qa.surv.data.region4$AccessionNumber, FUN=sum)
	null.survival.rate.per.accession.region4 <- null.deaths.per.accession.region4 / plants.per.accession[Qa.prov.data$GeographicRegion==4]
	null.mean.survival.rate.region4 <- mean(null.survival.rate.per.accession.region4)
	null.standard.deviation.survival.rate.region4 <- sd(null.survival.rate.per.accession.region4)
	null.opportunity.for.selection.region4[i] <- (null.standard.deviation.survival.rate.region4/null.mean.survival.rate.region4)^2
}

#examine null model results
summary(null.opportunity.for.selection.region4)

#define observed value of opportunity for selection
observed.opportunity.for.selection.region4 <- (standard.deviation.survival.rate.per.provenance["4"]/mean.survival.rate.per.provenance["4"])^2

#graphically contrast observed and null model results
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 5, 4, 2) + 0.1)
#examine range of values in the horizontal axis
range(observed.opportunity.for.selection.region4, null.opportunity.for.selection.region4) # 
#histogram with null values
hist(null.opportunity.for.selection.region4, breaks=seq(0, 1.2, 0.01),
	main="Region 4", col.main="black", xlab=expression(paste("Opportunity for selection (", italic(I[4]), ")", sep="")),
	ylab="Null model iterations", font.main=1, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,
	border="gray", col="gray90")
#add observed value
abline(v=observed.opportunity.for.selection.region4, col="#21908CFF", lty=1, lwd=2)

#add legend
legend("topright", c("Observed", "Null model"), fill=c(NA, "gray90"), border=c(NA, "gray"))
#add symbol for observed value, a bit fidgety
arrows(1.01, 475, 1.01, 505, length=0, angle=90, col="#21908CFF", lwd=2)

#calculate p-values
#lower tail
sum(null.opportunity.for.selection.region4 <= observed.opportunity.for.selection.region4)/length(null.opportunity.for.selection.region4)
#upper tail
sum(null.opportunity.for.selection.region4 >= observed.opportunity.for.selection.region4)/length(null.opportunity.for.selection.region4)

#############################################################################################
# 4.5. Compare null and observed opportunity for selection for region 5.

#set up null model for region 5
k <- 10000 #null model iterations
null.opportunity.for.selection.region5 <- rep(NA, times=k) #vector to store null opportunity for selection values
#extract survival data for region 5
index.region5 <- !is.na(match(Qa.surv.data$AccessionNumber, Qa.prov.data$AccessionNumber[Qa.prov.data$GeographicRegion==5]))
Qa.surv.data.region5 <- Qa.surv.data[index.region5,]
#examine result
dim(Qa.surv.data.region5)
head(Qa.surv.data.region5)
#the next two lines of code should produce the information (although in objects of different class)
plants.per.accession[Qa.prov.data$GeographicRegion==5]
table(Qa.surv.data.region5$AccessionNumber)
#check the match
identical(names(plants.per.accession[Qa.prov.data$GeographicRegion==5]), names(table(Qa.surv.data.region5$AccessionNumber)))
identical(as.vector(plants.per.accession[Qa.prov.data$GeographicRegion==5]), as.vector(table(Qa.surv.data.region5$AccessionNumber)))

#run null model
for(i in 1:k){
	null.dead.region5 <- sample(x=Qa.surv.data.region5$Dead, size=length(Qa.surv.data.region5$Dead))
	null.deaths.per.accession.region5 <- tapply(null.dead.region5, Qa.surv.data.region5$AccessionNumber, FUN=sum)
	null.survival.rate.per.accession.region5 <- null.deaths.per.accession.region5 / plants.per.accession[Qa.prov.data$GeographicRegion==5]
	null.mean.survival.rate.region5 <- mean(null.survival.rate.per.accession.region5)
	null.standard.deviation.survival.rate.region5 <- sd(null.survival.rate.per.accession.region5)
	null.opportunity.for.selection.region5[i] <- (null.standard.deviation.survival.rate.region5/null.mean.survival.rate.region5)^2
}

#examine null model results
summary(null.opportunity.for.selection.region5)

#define observed value of opportunity for selection
observed.opportunity.for.selection.region5 <- (standard.deviation.survival.rate.per.provenance["5"]/mean.survival.rate.per.provenance["5"])^2

#graphically contrast observed and null model results
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 5, 4, 2) + 0.1)
#examine range of values in the horizontal axis
range(observed.opportunity.for.selection.region5, null.opportunity.for.selection.region5) # 
#histogram with null values
hist(null.opportunity.for.selection.region5, breaks=seq(0, 0.99, 0.01),
	main="Region 5", col.main="black", xlab=expression(paste("Opportunity for selection (", italic(I[5]), ")", sep="")),
	ylab="Null model iterations", font.main=1, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,
	border="gray", col="gray90")
#add observed value
abline(v=observed.opportunity.for.selection.region5, col="#5DC863FF", lty=1, lwd=2)

#add legend
legend("topright", c("Observed", "Null model"), fill=c(NA, "gray90"), border=c(NA, "gray"))
#add symbol for observed value, a bit fidgety
arrows(0.835, 450, 0.835, 485, length=0, angle=90, col="#5DC863FF", lwd=2)

#calculate p-values
#lower tail
sum(null.opportunity.for.selection.region5 <= observed.opportunity.for.selection.region5)/length(null.opportunity.for.selection.region5)
#upper tail
sum(null.opportunity.for.selection.region5 >= observed.opportunity.for.selection.region5)/length(null.opportunity.for.selection.region5)


