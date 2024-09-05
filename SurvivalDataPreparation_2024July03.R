
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
#setwd("C:/_transfer/Papers/OpportunityForSelection/Data") #Ivan's directory

#################
#Homework: set your working directory below and please don't erase mine above. We will be sharing versions
#of this script and it is convenient not having to retype the working directory. Also, having the working
#directories in the code, as oppossed to using the "file.choose" function, bypasses the need to manually
#navigate to a directory. That will save time, given we will run each script many times and open and
#save several files each time.    
#################
#setwd("") #Lauren's directory

"C:/Users/laurd/OneDrive/Documents/MBG REU Internship/REU R Environment/REU Project"

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
#Homework: write annotations below, describing each of the variables (i.e., columns) in the
#data frame Qa.sd.   
#################

#Accession number represents a single mother plant.
#Plantings are siblings of the same mother.
#Position is the number used to randomize the position of a plant in the garden.
#SurverDataPossiblyDead is the date that the plant was observed as likely being dead.
#PossiblyDead is whether or not the plant was observed as possibly dead.
#SurveryDateDead is the data that plant death was confirmed.
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
#Homework: add a legend to the barplot above. Read the help page for function "legend".

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
#################
#Homework: modify the code below to edit the barplot and produce a high quality figure for
#publication in Conservation Biology. Read the journal author guidelines regarding figures and see
#examples in published papers. Focus on features of the barplot (including axis labels, axis annotation,
#font size, etc.) as it appears in the graphical window. Do not worry for now about saving the figure
#in a high resolution file. Read the help pages for functions "barplot" and "par".   
#################
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
#################
#Homework: modify the code below to edit the barplot and produce a high quality figure for
#publication in Conservation Biology. See indications for previous homework. Again, read the help pages
#for functions "barplot" and "par".
#################
barplot(table(Qa.sd$AccessionNumber)[o.FrequencyPlantings], las=3, col="transparent")
barplot(deaths.per.accession[o.FrequencyPlantings], add=T, col="gray60", axisnames=F)


######################################################################################################
## 7. Save the data frame named "Qa.sd" as a text file.
######################################################################################################

#set working directory
#setwd("C:/_transfer/Papers/OpportunityForSelection/Data") #Ivan's directory
#################
#Homework: set your working directory below and don't erase mine above. We will be sharing future versions
#of this script and it is convenient not having to retype the working directories.   
#################
#setwd("") #Lauren
write.csv(Qa.sd, file=paste("QaSurDat_", format(Sys.time(), "%Y%b%d_%H%M%S"), ".csv", sep=""), row.names = F)


