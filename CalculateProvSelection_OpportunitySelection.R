
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
## 1. [Survival data]
## 2. [Provenance data]
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
# 2. Calculate mean survival rates for provenances (i.e., geographic regions).
#############################################################################################

#############################################################################################
# 2.1. Observed mean survival rates for provenances (i.e., geographic regions).

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

#count surviving plants per maternal line (i.e., accession)
survive.per.accession <- tapply(Qa.surv.data$Dead, Qa.surv.data$AccessionNumber, FUN=function(x) sum(x==FALSE))
survive.per.accession
length(survive.per.accession)
class(survive.per.accession)
str(survive.per.accession)
barplot(table(survive.per.accession), ylab="Maternal lines", xlab="Surviving plants")

#calculate survival rate per maternal line (i.e., accession)
survival.rate.per.accession <- survive.per.accession / plants.per.accession
survival.rate.per.accession
length(survival.rate.per.accession)
class(survival.rate.per.accession)
str(survival.rate.per.accession)
hist(survival.rate.per.accession, breaks=seq(0, 1, 0.05))

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

#calculate mean survival rate per provenance (i.e., geographic region)
identical(names(survival.rate.per.accession), Qa.prov.data$AccessionNumber) #make sure accessions are ordered in the same way
mean.survival.rate.per.provenance <- tapply(survival.rate.per.accession, Qa.prov.data$GeographicRegion, FUN=mean)
mean.survival.rate.per.provenance
length(mean.survival.rate.per.provenance)
class(mean.survival.rate.per.provenance)
str(mean.survival.rate.per.provenance)

#############################################################################################
# 2.2. Null distribution of mean survival rates for provenances (i.e., geographic regions).

#set up null model
k <- 10000 #null model iterations
null.mean.survival.rate.per.provenance <- matrix(NA, nrow=k, ncol=5) #matrix to store null mean survival rates for provenances
colnames(null.mean.survival.rate.per.provenance) <- sort(unique(Qa.prov.data$GeographicRegion))
dim(null.mean.survival.rate.per.provenance) #matrix dimensions should equal the result of running the next line of code
c(k,length(unique(Qa.prov.data$GeographicRegion)))
head(null.mean.survival.rate.per.provenance)

#run null model
for(i in 1:k){
  null.dead <- sample(x=Qa.surv.data$Dead, size=length(Qa.surv.data$Dead))
  null.survival.per.accession <- tapply(null.dead, Qa.surv.data$AccessionNumber, FUN=function(x) sum(x==FALSE))
  null.survival.rate.per.accession <- null.survival.per.accession / plants.per.accession
  null.mean.survival.rate.per.provenance[i,] <- tapply(null.survival.rate.per.accession, Qa.prov.data$GeographicRegion, FUN=mean)
}

#examine null model results
summary(null.mean.survival.rate.per.provenance)

#############################################################################################
# 2.3. Calculate p-values for mean survival rates for provenances (i.e., geographic regions).

#region 2 lower tail
sum(null.mean.survival.rate.per.provenance[,"2"] <= mean.survival.rate.per.provenance["2"])/length(null.mean.survival.rate.per.provenance[,"2"])
#region 2 upper tail
sum(null.mean.survival.rate.per.provenance[,"2"] >= mean.survival.rate.per.provenance["2"])/length(null.mean.survival.rate.per.provenance[,"2"])
hist(null.mean.survival.rate.per.provenance[,"2"])
abline(v=mean.survival.rate.per.provenance["2"])

#region 3 lower tail
sum(null.mean.survival.rate.per.provenance[,"3"] <= mean.survival.rate.per.provenance["3"])/length(null.mean.survival.rate.per.provenance[,"3"])
#region 3 upper tail
sum(null.mean.survival.rate.per.provenance[,"3"] >= mean.survival.rate.per.provenance["3"])/length(null.mean.survival.rate.per.provenance[,"3"])
hist(null.mean.survival.rate.per.provenance[,"3"])
abline(v=mean.survival.rate.per.provenance["3"])

#region 4 lower tail
sum(null.mean.survival.rate.per.provenance[,"4"] <= mean.survival.rate.per.provenance["4"])/length(null.mean.survival.rate.per.provenance[,"4"])
#region 4 upper tail
sum(null.mean.survival.rate.per.provenance[,"4"] >= mean.survival.rate.per.provenance["4"])/length(null.mean.survival.rate.per.provenance[,"4"])
hist(null.mean.survival.rate.per.provenance[,"4"])
abline(v=mean.survival.rate.per.provenance["4"])

#region 5 lower tail
sum(null.mean.survival.rate.per.provenance[,"5"] <= mean.survival.rate.per.provenance["5"])/length(null.mean.survival.rate.per.provenance[,"5"])
#region 5 upper tail
sum(null.mean.survival.rate.per.provenance[,"5"] >= mean.survival.rate.per.provenance["5"])/length(null.mean.survival.rate.per.provenance[,"5"])
hist(null.mean.survival.rate.per.provenance[,"5"])
abline(v=mean.survival.rate.per.provenance["5"])

#region 8 lower tail
sum(null.mean.survival.rate.per.provenance[,"8"] <= mean.survival.rate.per.provenance["8"])/length(null.mean.survival.rate.per.provenance[,"8"])
#region 8 upper tail
sum(null.mean.survival.rate.per.provenance[,"8"] >= mean.survival.rate.per.provenance["8"])/length(null.mean.survival.rate.per.provenance[,"8"])
hist(null.mean.survival.rate.per.provenance[,"8"])
abline(v=mean.survival.rate.per.provenance["8"])


#############################################################################################
# 3. Relationship between the number of maternal lines and mean survival rates
#    across provenances (i.e., geographic regions).
#############################################################################################

#calculate median of the null mean survival rates for provenances (i.e., geographic regions)
median.null.mean.survival.rate.per.provenance <- apply(X=null.mean.survival.rate.per.provenance, MARGIN=2, FUN=quantile, probs=0.5)
#calculate mean of the mean survival rates for provenances (i.e., geographic regions)
mean.null.mean.survival.rate.per.provenance <- apply(X=null.mean.survival.rate.per.provenance, MARGIN=2, FUN=mean)
#calculate lower 0.025 quantile of the null mean survival rates for provenances (i.e., geographic regions)
LL.null.mean.survival.rate.per.provenance <- apply(X=null.mean.survival.rate.per.provenance, MARGIN=2, FUN=quantile, probs=0.025)
#calculate upper 0.975 quantile of null selection coefficients
UL.null.mean.survival.rate.per.provenance <- apply(X=null.mean.survival.rate.per.provenance, MARGIN=2, FUN=quantile, probs=0.975)

#plot number of maternal lines for each provenance (i.e., geographic region) in the horizontal axis
#and the respective observed and null mean survival rate in the vertical axis
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 5, 4, 2) + 0.1)
#examine range of values in the horizontal axis
range(maternal.lines.per.provenance)
#examine range of values in the vertical axis
range(mean.survival.rate.per.provenance, LL.null.mean.survival.rate.per.provenance, UL.null.mean.survival.rate.per.provenance)
#set up plot
plot(maternal.lines.per.provenance, mean.null.mean.survival.rate.per.provenance, xlim=c(0, 16), ylim=c(0, 1),
     bty="n", cex.axis=1.5, cex.lab=1.5, xlab=expression("Number of maternal lines, "* italic(n[i])), ylab=expression(paste("Mean survival rate, ", italic(bar(S[r])(t)), sep="")), type="n")
#arrows(0, 0, 15, 0, code=0, lty=2, col="gray80", lwd=2) #add horizontal line marking selection coefficients equal to zero
#graph mean of null selection coefficients
points(maternal.lines.per.provenance, mean.null.mean.survival.rate.per.provenance, pch=19, col="gray70", cex=1)
#graph median of null selection coefficients
points(maternal.lines.per.provenance, median.null.mean.survival.rate.per.provenance, pch=21, col="gray70", cex=3)
#graph 95% confidence interval of null selection coefficients
arrows(maternal.lines.per.provenance, LL.null.mean.survival.rate.per.provenance, maternal.lines.per.provenance, UL.null.mean.survival.rate.per.provenance,
       length=0.1, angle=90, code=3, col="gray70")
#graph observed selection coefficients
#define shapes and colors
shapes <- c(21,22,23,24,25)
colors <- viridis::viridis(5, option = "D")
points(maternal.lines.per.provenance, mean.survival.rate.per.provenance, pch=shapes, col="black", bg = colors, cex=2.0)

#a version of the legend with p-values, separating observed and null values
#observed values
legend(x=6, y=1.2,
       legend = c(expression(paste("Observed region 2, ", italic(p), " = 0.456")),
                  expression(paste("Observed region 3, ", italic(p), " = 0.104")),
                  expression(paste("Observed region 4, ", italic(p), " = 0.374")),
                  expression(paste("Observed region 5, ", italic(p), " = 0.293")),
                  expression(paste("Observed region 8, ", italic(p), " = 0.141"))),
       pch=c(21,22,23,24,25),
       pt.cex=2,
       col=c("black", "black", "black", "black", "black"),
       pt.bg = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"),
       ncol = 1,
       cex = 1.22,
       xpd = TRUE)
#null values
legend(x=8, y=0.35,
       legend = c("Null model mean", "Null model median", "Null model 95% CI"),
       pch=c(19,21,NA),
       pt.cex=c(1,3,2.3),
       col=c("gray70", "gray70", "gray70"),
       pt.bg = NA,
       ncol = 1,
       cex = 1.22,
       xpd = TRUE)
arrows(8.4, 0.1, 8.4, 0.15, length=0.1, angle=90, code=3, col="gray70")

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
#region 3
(standard.deviation.survival.rate.per.provenance["4"]/mean.survival.rate.per.provenance["4"])^2
#region 3
(standard.deviation.survival.rate.per.provenance["5"]/mean.survival.rate.per.provenance["5"])^2
#region 3
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
  null.survive.per.accession.region2 <- tapply(null.dead.region2, Qa.surv.data.region2$AccessionNumber, FUN=function(x) sum(x==FALSE))
  null.survival.rate.per.accession.region2 <- null.survive.per.accession.region2 / plants.per.accession[Qa.prov.data$GeographicRegion==2]
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
hist(null.opportunity.for.selection.region2, breaks=seq(0, 0.7, 0.01),
     main="", col.main="black", xlab=expression(paste("Opportunity for selection, ", italic(I[2]), sep="")),
     ylab="Null model iterations (bars)", font.main=1, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,
     border="gray", col="gray90")
#add observed value
abline(v=observed.opportunity.for.selection.region2, col="#440154FF", lty=1, lwd=2)
text(observed.opportunity.for.selection.region2, 620, "Observed region 2", cex=1.5, xpd=TRUE, col="#440154FF")
mtext(side=2, "a)", cex=1.5, las=2,  at=650, line=4, adj=0)
text(0.5, 400, expression(paste(italic(p), " = 0.116")), cex=1.5)

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
  null.survive.per.accession.region3 <- tapply(null.dead.region3, Qa.surv.data.region3$AccessionNumber, FUN=function(x) sum(x==FALSE))
  null.survival.rate.per.accession.region3 <- null.survive.per.accession.region3 / plants.per.accession[Qa.prov.data$GeographicRegion==3]
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
     main="", col.main="black", xlab=expression(paste("Opportunity for selection, ", italic(I[3]), sep="")),
     ylab="Null model iterations (bars)", font.main=1, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,
     border="gray", col="gray90")
#add observed value
abline(v=observed.opportunity.for.selection.region3, col="#3B528BFF", lty=1, lwd=2)
text(observed.opportunity.for.selection.region3-0.5, 5750, "Observed region 3", cex=1.5, xpd=TRUE, col="#3B528BFF")
mtext(side=2, "b)", cex=1.5, las=2,  at=6000, line=4, adj=0)
text(2.5, 3700, expression(paste(italic(p), " = 0.018")), cex=1.5)

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
  null.survive.per.accession.region4 <- tapply(null.dead.region4, Qa.surv.data.region4$AccessionNumber, FUN=function(x) sum(x==FALSE))
  null.survival.rate.per.accession.region4 <- null.survive.per.accession.region4 / plants.per.accession[Qa.prov.data$GeographicRegion==4]
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
     main="", col.main="black", xlab=expression(paste("Opportunity for selection, ", italic(I[4]), sep="")),
     ylab="Null model iterations (bars)", font.main=1, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,
     border="gray", col="gray90")
#add observed value
abline(v=observed.opportunity.for.selection.region4, col="#21908CFF", lty=1, lwd=2)
text(observed.opportunity.for.selection.region4, 750, "Observed region 4", cex=1.5, xpd=TRUE, col="#21908CFF")
mtext(side=2, "c)", cex=1.5, las=2,  at=800, line=4, adj=0)
text(0.8, 500, expression(paste(italic(p), " = 0.321")), cex=1.5)

#add legend
legend("topright", c("Observed", "Null model"), fill=c(NA, "gray90"), border=c(NA, "gray"))
#add symbol for observed value, a bit fidgety
arrows(1.01, 465, 1.01, 500, length=0, angle=90, col="#21908CFF", lwd=2)

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
  null.survive.per.accession.region5 <- tapply(null.dead.region5, Qa.surv.data.region5$AccessionNumber, FUN=function(x) sum(x==FALSE))
  null.survival.rate.per.accession.region5 <- null.survive.per.accession.region5 / plants.per.accession[Qa.prov.data$GeographicRegion==5]
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
hist(null.opportunity.for.selection.region5, breaks=seq(0, 0.85, 0.01),
     main="", col.main="black", xlab=expression(paste("Opportunity for selection, ", italic(I[5]), sep="")),
     ylab="Null model iterations (bars)", font.main=1, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,
     border="gray", col="gray90")
#add observed value
abline(v=observed.opportunity.for.selection.region5, col="#5DC863FF", lty=1, lwd=2)
text(observed.opportunity.for.selection.region5, 600, "Observed region 5", cex=1.5, xpd=TRUE, col="#5DC863FF")
mtext(side=2, "d)", cex=1.5, las=2,  at=630, line=4, adj=0)
text(0.6, 400, expression(paste(italic(p), " = 0.409")), cex=1.5)

#add legend
legend("topright", c("Observed", "Null model"), fill=c(NA, "gray90"), border=c(NA, "gray"))
#add symbol for observed value, a bit fidgety
arrows(0.715, 468, 0.715, 500, length=0, angle=90, col="#5DC863FF", lwd=2)

#calculate p-values
#lower tail
sum(null.opportunity.for.selection.region5 <= observed.opportunity.for.selection.region5)/length(null.opportunity.for.selection.region5)
#upper tail
sum(null.opportunity.for.selection.region5 >= observed.opportunity.for.selection.region5)/length(null.opportunity.for.selection.region5)


