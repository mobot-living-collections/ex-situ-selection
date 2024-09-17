
######################################################################################################
######################################################################################################
##
## TITLE
##
## Quercus arkansana survival data preparation and validation.
##
## INTRODUCTION
##
## This script is part of a project on provenance selection and opportunity for selection in ex
## situ collections for plant conservation. The project uses as a study system a living collection
## of Quercus arkansana (Fagaceae) maintained in the Oertli Hardy Plant Nursery of the
## Missouri Botanical Garden.
##
## DATA FILE REQUIRED
##
## 1. "QaSurvivalData20230922.csv".
##
## CONTENTS
##
## 1. Read and examine data files.
## 2. Edit the data frame "Qa.surv.data" and create a new data frame named "Qa.sd".
## 3. Examine survival survey dates.
## 4. Compare survival data between the two surveys.
## 5. Examine the distribution of plantings among accessions.
## 6. Examine the distribution of planting deaths among accessions.
## 7. Save the data frame named "Qa.sd" as a text file.
##
######################################################################################################
######################################################################################################


######################################################################################################
## 1. Read and examine data files.
######################################################################################################

#set working directory

#read and examine survival data
Qa.surv.data <- read.table(file="QaSurvivalData20230922.csv", header=TRUE, sep=",")
class(Qa.surv.data)
dim(Qa.surv.data)
head(Qa.surv.data)
summary(Qa.surv.data)

######################################################################################################
## 2. Edit the data frame "Qa.surv.data" and create a new data frame named "Qa.sd".
######################################################################################################

#remove (non-relevant) data on the latitude and longitude of the Oertli Nursery
Qa.sd <- Qa.surv.data[,-c(4,5)]

#format survival survey dates
Qa.sd$SurveyDatePossiblyDead <- as.Date(Qa.surv.data$SurveyDatePossiblyDead, format="%m/%d/%Y %H:%M")
Qa.sd$SurveyDateDead <- as.Date(Qa.surv.data$SurveyDateDead, format="%m/%d/%Y")

#order the rows by accession number and planting number
o.AccessionPlanting <- order(Qa.sd$AccessionNumber, Qa.sd$PlantingNumber)
Qa.sd <- Qa.sd[o.AccessionPlanting,]

#examine the new data frame
head(Qa.sd, n=20)
tail(Qa.sd, n=20)
summary(Qa.sd)

#################
#annotations below, describing each of the variables (i.e., columns) in the
#data frame Qa.sd.   
#################

#Accession number represents a single mother plant.
#Plantings are siblings of the same mother.
#Position is the number used to randomize the position of a plant in the garden.
#SurveyDataPossiblyDead is the date that the plant was observed as likely being dead.
#PossiblyDead is whether or not the plant was observed as possibly dead.
#SurveyDateDead is the data that plant death was confirmed.
#Dead is whether or not the plant was actually dead.

######################################################################################################
## 3. Examine survival survey dates.
######################################################################################################

#there are three dates for the first survey:
unique(Qa.sd$SurveyDatePossiblyDead)

#graph the distribution of survey data among the three dates
PossiblyDead.per.date <- tapply(Qa.sd$PossiblyDead, Qa.sd$SurveyDatePossiblyDead, FUN=sum)
#par(mar=c(5, 4, 4, 2) + 0.1) #default
par(mar=c(5, 4.5, 4, 2) + 0.1)
barplot(table(Qa.sd$SurveyDatePossiblyDead), col="gray70", xlab="Date", ylab="Plantings", cex.lab=1.5, cex.axis=1.5, cex.names=1.5)
barplot(PossiblyDead.per.date, add=T, density=15, col="black", axes=F, axisnames=F)
#################
#add a legend to the barplot above
legend("topright", c("Alive", "Possibly Dead"), col=c("gray70", "gray70"), fill=c("gray70", "black"), cex=1.1, density=c(NA,15))
#################

#there is only one date for the second survey:
unique(Qa.sd$SurveyDateDead)

######################################################################################################
## 4. Compare survival data between the two surveys.
######################################################################################################

#these rows have different survival data in the two surveys
Qa.sd[Qa.sd$PossiblyDead != Qa.sd$Dead,]

######################################################################################################
## 5. Examine the distribution of plantings among accessions.
######################################################################################################

#examine unique accessions
unique(Qa.sd$AccessionNumber)
length(unique(Qa.sd$AccessionNumber)) #the number of unique accessions

#examine the distribution of plantings per accession,
table(Qa.sd$AccessionNumber)
#ordering accessions according to decreasing number of plantings
sort(table(Qa.sd$AccessionNumber), decreasing=T)

#graph the distribution of plantings per accession

barplot(sort(table(Qa.sd$AccessionNumber), decreasing=T), las=3) 

barplot(sort(table(Qa.sd$AccessionNumber), decreasing=T), las=3, ylim=c(0, 11), ylab="Plantings", xlab="Accession Number")
axis(side = 2, at = seq(0,11))

######################################################################################################
## 6. Examine the distribution of planting deaths among accessions.
######################################################################################################

#overall mortality
summary(Qa.sd$Dead)
sum(Qa.sd$Dead)/length(Qa.sd$Dead) #proportion dead plantings

#obtain and examine the distribution of deaths among accessions
deaths.per.accession <- tapply(Qa.sd$Dead, Qa.sd$AccessionNumber, FUN=sum)
deaths.per.accession

#graph the distribution of deaths among accessions
barplot(deaths.per.accession, las=3)

#graph the distribution of deaths among accessions,
#this time ordering accessions according to decreasing number of plantings;
#first use the function "order" to order accessions: 
o.FrequencyPlantings <- order(table(Qa.sd$AccessionNumber), decreasing=T)
#make sure the order of accessions is correct
table(Qa.sd$AccessionNumber)[o.FrequencyPlantings]
deaths.per.accession[o.FrequencyPlantings]
identical(names(table(Qa.sd$AccessionNumber)[o.FrequencyPlantings]), names(deaths.per.accession[o.FrequencyPlantings]))
#now create a barplot
barplot(deaths.per.accession[o.FrequencyPlantings], las=3)

#graph the distribution of deaths among accessions on top of the number of plantings per accession

barplot(table(Qa.sd$AccessionNumber)[o.FrequencyPlantings], las=3, col="transparent")
barplot(deaths.per.accession[o.FrequencyPlantings], add=T, col="gray60", axisnames=F)

######################################################################################################
## 7. Save the data frame named "Qa.sd" as a text file.
######################################################################################################

write.csv(Qa.sd, file=paste("QaSurDat_", format(Sys.time(), "%Y%b%d_%H%M%S"), ".csv", sep=""), row.names = F)


