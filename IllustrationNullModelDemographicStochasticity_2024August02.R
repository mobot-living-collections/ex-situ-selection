
######################################################################################################
######################################################################################################
##
## TITLE
##
## Illustration of the null model, representing demographic stochasticity, used to gauge provenance
## selection and opportunity for selection in an ex situ collection of Quercus arkansana
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
## 2. Distribution of initial plants per maternal line per provenance (i.e., geographic region). 
## 3. Distribution of surviving plants per maternal line per provenance (i.e., geographic region).
## 4. Demographic stochasticity null model.
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
# 2. Distribution of initial plants per maternal line per provenance (i.e., geographic region). 
#############################################################################################

#############################################################################################
# 2.1. Count the number of initial plants per maternal line (i.e., accession).

plants.per.accession <- tapply(Qa.surv.data$Dead, Qa.surv.data$AccessionNumber, FUN=length)
plants.per.accession
length(plants.per.accession)
class(plants.per.accession)
str(plants.per.accession)
barplot(table(plants.per.accession), ylab="Maternal lines", xlab="Plants")

#############################################################################################
# 2.2. Graph the number of initial plants per accession, in increasing order for each provenance (i.e., geographic region).

identical(names(plants.per.accession), Qa.prov.data$AccessionNumber) #make sure accessions are ordered in the same way
region.plants.order <- order(Qa.prov.data$GeographicRegion, plants.per.accession)
#define colors for bar plot
bar.col <- rep(c("#D81B60", "#1E88E5", "#FFC107", "#42FD7F", "#E7A0FD"),
	times=c(sum(Qa.prov.data$GeographicRegion==2), sum(Qa.prov.data$GeographicRegion==3),
	sum(Qa.prov.data$GeographicRegion==4), sum(Qa.prov.data$GeographicRegion==5),
	sum(Qa.prov.data$GeographicRegion==8)))
#par(mar=c(5, 4, 4, 2) + 0.1) #default margins
par(mar=c(5, 4.5, 4, 2) + 0.1)
barplot(plants.per.accession[region.plants.order], ylim=c(0,11), names.arg="", ylab="Initial plants", cex.lab=1.5, cex.axis=1.5, col=bar.col)
#edit vertical axis
axis(2, at=seq(, 10, 1), labels=T, cex.axis=1.5)
#add axis 1 for region 2
ml.r2 <- sum(Qa.prov.data$GeographicRegion==2)
axis(1, at=c(0.2, ml.r2*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("2", side=1, at=0.2 + ml.r2*(1+0.2)/2, line=1.2, cex=1.5) 
#add axis 1 for region 3
ml.r3 <- sum(Qa.prov.data$GeographicRegion==3)
axis(1, at=c(0.2 + ml.r2*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("3", side=1, at=ml.r2*(1+0.2) + 0.2 + (ml.r3*(1+0.2))/2, line=1.2, cex=1.5) 
#add axis 1 for region 4
ml.r4 <- sum(Qa.prov.data$GeographicRegion==4)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("4", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + 0.2 + (ml.r4*(1+0.2))/2, line=1.2, cex=1.5) 
#add axis 1 for region 5
ml.r5 <- sum(Qa.prov.data$GeographicRegion==5)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("5", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + 0.2 + (ml.r5*(1+0.2))/2, line=1.2, cex=1.5) 
#add axis 1 for region 8
ml.r8 <- sum(Qa.prov.data$GeographicRegion==8)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2) + ml.r8*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("8", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2) + 0.2 + (ml.r8*(1+0.2))/2, line=1.2, cex=1.5) 
#add axis 1 title
mtext("Region", side=1, line=3, cex=1.5) 
#add legend
legend(16.5*(1+0.2),10.5, paste("Maternal lines region", c(2,3,4,5,8)), fill=c("#D81B60", "#1E88E5", "#FFC107", "#42FD7F", "#E7A0FD"), cex=1.2)


#############################################################################################
# 3. Distribution of surviving plants per maternal line per provenance (i.e., geographic region).
#############################################################################################

#############################################################################################
# 3.1. Count the number of surviving plants per maternal line (i.e., accession).

surviving.plants.per.accession <- tapply(Qa.surv.data$Dead == F, Qa.surv.data$AccessionNumber, FUN=sum)
surviving.plants.per.accession
length(surviving.plants.per.accession)
class(surviving.plants.per.accession)
str(surviving.plants.per.accession)
barplot(table(surviving.plants.per.accession), ylab="Maternal lines", xlab="Plants")

#############################################################################################
# 3.2. Graph the number of surviving plants per accession, in the order used in section 2.

identical(names(surviving.plants.per.accession), Qa.prov.data$AccessionNumber) #make sure accessions are ordered in the same way
region.plants.order <- order(Qa.prov.data$GeographicRegion, plants.per.accession)
#define colors for bar plot, see: https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40-%23664ef8-%2342fd7f
bar.col <- rep(c("#D81B60", "#1E88E5", "#FFC107", "#42FD7F", "#E7A0FD"),
	times=c(sum(Qa.prov.data$GeographicRegion==2), sum(Qa.prov.data$GeographicRegion==3),
	sum(Qa.prov.data$GeographicRegion==4), sum(Qa.prov.data$GeographicRegion==5),
	sum(Qa.prov.data$GeographicRegion==8)))
#par(mar=c(5, 4, 4, 2) + 0.1) #default margins
par(mar=c(5, 4.5, 4, 2) + 0.1)
barplot(surviving.plants.per.accession[region.plants.order], ylim=c(0,11), names.arg="", ylab="Surviving plants", cex.lab=1.5, cex.axis=1.5, col=bar.col)
#edit vertical axis
axis(2, at=seq(, 10, 1), labels=T, cex.axis=1.5)
#add horizontal axis for region 2
ml.r2 <- sum(Qa.prov.data$GeographicRegion==2)
axis(1, at=c(0.2, ml.r2*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("2", side=1, at=0.2 + ml.r2*(1+0.2)/2, line=1.2, cex=1.5) 
#add horizontal axis for region 3
ml.r3 <- sum(Qa.prov.data$GeographicRegion==3)
axis(1, at=c(0.2 + ml.r2*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("3", side=1, at=ml.r2*(1+0.2) + 0.2 + (ml.r3*(1+0.2))/2, line=1.2, cex=1.5) 
#add horizontal axis for region 4
ml.r4 <- sum(Qa.prov.data$GeographicRegion==4)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("4", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + 0.2 + (ml.r4*(1+0.2))/2, line=1.2, cex=1.5) 
#add horizontal axis for region 5
ml.r5 <- sum(Qa.prov.data$GeographicRegion==5)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("5", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + 0.2 + (ml.r5*(1+0.2))/2, line=1.2, cex=1.5) 
#add horizontal axis for region 8
ml.r8 <- sum(Qa.prov.data$GeographicRegion==8)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2) + ml.r8*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("8", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2) + 0.2 + (ml.r8*(1+0.2))/2, line=1.2, cex=1.5) 
#add horizontal axis title
mtext("Region", side=1, line=3, cex=1.5) 
#add legend
legend(16.5*(1+0.2),10.5, paste("Maternal lines region", c(2,3,4,5,8)), fill=c("#D81B60", "#1E88E5", "#FFC107", "#42FD7F", "#E7A0FD"), cex=1.2)

#############################################################################################
# 4. Demographic stochasticity null model.
#############################################################################################

#############################################################################################
# 4.1. Generate null distribution of surviving plants.

#set up null model
k <- 10000 #null model iterations
null.surviving.plants.per.accession <- matrix(NA, nrow=k, ncol=length(unique(Qa.surv.data$AccessionNumber)))
identical(unique(Qa.surv.data$AccessionNumber), names(surviving.plants.per.accession)) #make sure accessions are ordered in the same way
colnames(null.surviving.plants.per.accession) <- unique(Qa.surv.data$AccessionNumber)
dim(null.surviving.plants.per.accession) #matrix dimensions should equal the result of running the next line of code 
c(k,length(unique(Qa.surv.data$AccessionNumber)))
head(null.surviving.plants.per.accession)

#run null model
for(i in 1:k){
	null.dead <- sample(x=Qa.surv.data$Dead, size=length(Qa.surv.data$Dead))
	null.surviving.plants.per.accession[i,] <- tapply(null.dead == F, Qa.surv.data$AccessionNumber, FUN=sum)
}

#examine null model results
summary(null.surviving.plants.per.accession)

#############################################################################################
# 4.2. Graph the null and observed distributions of surviving plants per accession, in the order used in section 2.

region.plants.order <- order(Qa.prov.data$GeographicRegion, plants.per.accession)
#define colors for bar plot, see: https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40-%23664ef8-%2342fd7f
bar.col <- rep(c("#D81B60", "#1E88E5", "#FFC107", "#42FD7F", "#E7A0FD"),
	times=c(sum(Qa.prov.data$GeographicRegion==2), sum(Qa.prov.data$GeographicRegion==3),
	sum(Qa.prov.data$GeographicRegion==4), sum(Qa.prov.data$GeographicRegion==5),
	sum(Qa.prov.data$GeographicRegion==8)))
#par(mar=c(5, 4, 4, 2) + 0.1) #default margins
par(mar=c(5, 4.5, 4, 2) + 0.1)
barplot(surviving.plants.per.accession[region.plants.order], ylim=c(0,11), names.arg="", ylab="Surviving plants", cex.lab=1.5, cex.axis=1.5, col=bar.col)
#edit vertical axis
axis(2, at=seq(, 10, 1), labels=T, cex.axis=1.5)
#add horizontal axis for region 2
ml.r2 <- sum(Qa.prov.data$GeographicRegion==2)
axis(1, at=c(0.2, ml.r2*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("2", side=1, at=0.2 + ml.r2*(1+0.2)/2, line=1.2, cex=1.5) 
#add horizontal axis for region 3
ml.r3 <- sum(Qa.prov.data$GeographicRegion==3)
axis(1, at=c(0.2 + ml.r2*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("3", side=1, at=ml.r2*(1+0.2) + 0.2 + (ml.r3*(1+0.2))/2, line=1.2, cex=1.5) 
#add horizontal axis for region 4
ml.r4 <- sum(Qa.prov.data$GeographicRegion==4)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("4", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + 0.2 + (ml.r4*(1+0.2))/2, line=1.2, cex=1.5) 
#add horizontal axis for region 5
ml.r5 <- sum(Qa.prov.data$GeographicRegion==5)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("5", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + 0.2 + (ml.r5*(1+0.2))/2, line=1.2, cex=1.5) 
#add horizontal axis for region 8
ml.r8 <- sum(Qa.prov.data$GeographicRegion==8)
axis(1, at=c(0.2 + ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2), ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2) + ml.r8*(1+0.2)),
	labels=F, line=0.5, lwd=1.5)
mtext("8", side=1, at=ml.r2*(1+0.2) + ml.r3*(1+0.2) + ml.r4*(1+0.2) + ml.r5*(1+0.2) + 0.2 + (ml.r8*(1+0.2))/2, line=1.2, cex=1.5) 
#add horizontal axis title
mtext("Region", side=1, line=3, cex=1.5) 

#add demographic stochasticity null model
LL.null.surviving.plants.per.accession <- apply(X=null.surviving.plants.per.accession, MARGIN=2, FUN=quantile, probs=0.05)[region.plants.order]
UL.null.surviving.plants.per.accession <- apply(X=null.surviving.plants.per.accession, MARGIN=2, FUN=quantile, probs=0.975)[region.plants.order]
x.coor.bars <- seq(0.2 + 0.5, 42*(1+0.2)-0.5 , 1.2)
arrows(x.coor.bars, LL.null.surviving.plants.per.accession, x.coor.bars, UL.null.surviving.plants.per.accession,
	length=0.02, angle=90, code=3, col="black")

#add legend
legend(5*(1+0.2), 10.5, "Demographic stochasticity null model", pch=NA, cex=1.2)
#add symbol (a bit fidgety)
arrows(7.8, 9.7, 7.8, 10.3, length=0.02, angle=90, code=3, col="black")

