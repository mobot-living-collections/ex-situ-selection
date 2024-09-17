
#############################################################################################
#############################################################################################
##
## TITLE
##
## Preparation and validation of occurrence and provenance data for Quercus arkansana.
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
## Three data files are required:
##
## 1. "ExpeditionDataALL_Static2023-06-21.txt" [Latitude and longitude (i.e., provenance) data for collected accessions]
## 2. "OccurrenceDatabaseSubset2023-06-21.txt" [Latitude and longitude data for all occurences of Quercus arkansana]
## 3. "QaSurDat_2024Jun21_122839.csv" [Survival data]
## 4. "gadm40_USA_1.shp" [shape file with USA state boundaries]
##
## CONTENTS
##
## 1. Read and examine data files.
## 2. Obtain occurrence and provenance data for Quercus arkansana.
## 3. Map known occurrences of Q. arkansana and collection localities of the unique
##    accessions in the common garden experiment.
## 4. Assign accessions to geographic regions, following Thomas et al. (2023, Biological
#     Conservation 283: 110052).
## 5. Save provenance data as a text file named "QaAccCoor.csv".
##
#############################################################################################
#############################################################################################


#############################################################################################
# 1. Read and examine data files.
#############################################################################################

#set working directory

accession.provenance <- read.table("ExpeditionDataALL_Static2023-06-21.txt", header = T, sep ="\t")
dim(accession.provenance) 
colnames(accession.provenance)
str(accession.provenance)
head(accession.provenance)

occurrence <- read.table("OccurrenceDatabaseSubset2023-06-21.txt", header = T, sep = "\t")
dim(occurrence)
str(occurrence)
head(occurrence)

Qa.sd <- read.table("QaSurDat_2024Jun21_122839.csv", header = T, sep = ",")
dim(Qa.sd) 
str(Qa.sd) 
head(Qa.sd)


#############################################################################################
# 2. Obtain occurrence and provenance data for Quercus arkansana.
#############################################################################################

#extract occurrence data for Q. arkansana
unique(occurrence$Taxon)
Qa.occurrence <- occurrence[occurrence$Taxon == "Quercus arkansana",]
dim(Qa.occurrence) #should have 792 rows and 18 columns
str(Qa.occurrence)
head(Qa.occurrence)

#extract accession provenace data for MBG accessions of Q. arkansana
sort(unique(accession.provenance$Taxon))
Qa.provenance <- data.frame(accession.provenance[accession.provenance$Taxon == "Quercus arkansana" | 
                            accession.provenance$Taxon == "Quercus arkansana " |
                            accession.provenance$Taxon == "Quercus arkansana?",])
dim(Qa.provenance) #should have 317 rows and 44 columns
str(Qa.provenance)
head(Qa.provenance)

#extract coordinates for the collection locality of each unique accession number in the common garden experiment
Qa.accessions <- sort(unique(Qa.sd$AccessionNumber)) 
Qa.accessions.coor <- Qa.provenance[match(Qa.accessions, Qa.provenance$AccessionNumber),c("DDLatitude", "DDLongitude")]


#############################################################################################
# 3. Map known occurrences of Q. arkansana and collection localities of the unique accessions
#    in the common garden experiment.
#############################################################################################

#load package "terra"
library(terra)

#load state boundaries
States <- vect("gadm40_USA_1.shp")

#examine range of geographic coordinates
range(Qa.occurrence$DDLongitud)
range(Qa.occurrence$DDLatitude)

#map 1: shows the whole range of the species
windows(14,14) #for windows machines only
plot(States, xlim=c(-95, -81.5), ylim=c(30, 34.2), axes=F, box=T, mar=c(3.1,6,2.1,1))
axis(1, line=-12.2, seq(-94, -82, 1), labels=F, cex.axis=2)
axis(1, line=-12.2, seq(-94, -82, 2), labels=T, cex.axis=2)
axis(1, line=-12.2, cex.axis=2)
axis(2, line=0, seq(30, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-9)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud, Qa.occurrence$DDLatitude, type = "p", pch = 19, col = "blue")
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 24, col = "red", cex=2)

#map 2: focused on the western side of the range of the species, showing all collection localities
#of the accessions in the commmon garden.
windows(14,14)
plot(States, xlim=c(-95, -89), ylim=c(30, 34.2), axes=F, box=T, mar=c(6,6,2.1,1))
axis(1, line=0, seq(-95, -89, 1), labels=T, cex.axis=2)
axis(2, line=-0.65, seq(30, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=3.7)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud, Qa.occurrence$DDLatitude, type = "p", pch = 19, col = "blue")
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 24, col = "red", cex=2)

#map 3: focused on the northwestern side of the range of the species, showing all collection localities
#of the accessions in the commmon garden.
windows(14,14)
plot(States, xlim=c(-94.5, -89.2), ylim=c(31.5, 34.2), axes=F, box=T, mar=c(6,6,2.1,1))
axis(1, line=-4.9, seq(-94, -90), labels=T, cex.axis=2)
axis(2, line=0, seq(32, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-1)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud[Qa.occurrence$DDLatitude>31.4], Qa.occurrence$DDLatitude[Qa.occurrence$DDLatitude>31.4], type = "p", pch = 19, col = "blue")
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 24, col = "red", cex=2)


#############################################################################################
# 4. Assign accessions to geographic regions, following Thomas et al. (2023, Biological
#    Conservation 283: 110052).
#############################################################################################

geographic.region <- rep(NA, length=nrow(Qa.accessions.coor))

#geographic region 2
#windows(14,14)
plot(States, xlim=c(-94.5, -89.2), ylim=c(31.5, 34.2), axes=F, box=T, mar=c(6,6,2.1,1))
axis(1, line=-4.9, seq(-94, -90), labels=T, cex.axis=2)
axis(2, line=0, seq(32, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-1)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud[Qa.occurrence$DDLatitude>31.4], Qa.occurrence$DDLatitude[Qa.occurrence$DDLatitude>31.4], type = "p", pch = 19, col = "blue")
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 24, col = "red", cex=2)
abline(v = c(-92.5, -92), lty=3, lwd=2)
geographic.region[Qa.accessions.coor$DDLongitude > -92.5 & Qa.accessions.coor$DDLongitude < -92] <- 2
points(Qa.accessions.coor$DDLongitude[geographic.region==2], Qa.accessions.coor$DDLatitude[geographic.region==2], type = "p", pch = 21, col = "orange", cex=3, lwd=2)

#geographic region 3
#windows(14,14)
plot(States, xlim=c(-94.5, -89.2), ylim=c(31.5, 34.2), axes=F, box=T, mar=c(6,6,2.1,1))
axis(1, line=-4.9, seq(-94, -90), labels=T, cex.axis=2)
axis(2, line=0, seq(32, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-1)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud[Qa.occurrence$DDLatitude>31.4], Qa.occurrence$DDLatitude[Qa.occurrence$DDLatitude>31.4], type = "p", pch = 19, col = "blue", lwd=2)
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 24, col = "red", cex=2)
abline(v = - 93.5, lty=3, lwd=2)
abline(h = 33.5, lty=3, lwd=2)
geographic.region[Qa.accessions.coor$DDLongitude < -93.5 & Qa.accessions.coor$DDLatitude < 33.5] <- 3
points(Qa.accessions.coor$DDLongitude[geographic.region==3], Qa.accessions.coor$DDLatitude[geographic.region==3], type = "p", pch = 21, col = "yellow", cex=3, lwd=2)

#geographic region 4
#windows(14,14)
plot(States, xlim=c(-94.5, -89.2), ylim=c(31.5, 34.2), axes=F, box=T, mar=c(6,6,2.1,1))
axis(1, line=-4.9, seq(-94, -90), labels=T, cex.axis=2)
axis(2, line=0, seq(32, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-1)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud[Qa.occurrence$DDLatitude>31.4], Qa.occurrence$DDLatitude[Qa.occurrence$DDLatitude>31.4], type = "p", pch = 19, col = "blue", lwd=2)
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 24, col = "red", cex=2)
abline(v = - 93.5, lty=3, lwd=2)
abline(h = 33.5, lty=3, lwd=2)
geographic.region[Qa.accessions.coor$DDLongitude > -93.5 & Qa.accessions.coor$DDLatitude > 33.5] <- 4
points(Qa.accessions.coor$DDLongitude[geographic.region==4], Qa.accessions.coor$DDLatitude[geographic.region==4], type = "p", pch = 21, col = "lawngreen", cex=3, lwd=2)

#geographic region 5
#windows(14,14)
plot(States, xlim=c(-94.5, -89.2), ylim=c(31.5, 34.2), axes=F, box=T, mar=c(6,6,2.1,1))
axis(1, line=-4.9, seq(-94, -90), labels=T, cex.axis=2)
axis(2, line=0, seq(32, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-1)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud[Qa.occurrence$DDLatitude>31.4], Qa.occurrence$DDLatitude[Qa.occurrence$DDLatitude>31.4], type = "p", pch = 19, col = "blue", lwd=2)
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 24, col = "red", cex=2)
abline(v = - 93.5, lty=3, lwd=2)
abline(h = 33.5, lty=3, lwd=2)
geographic.region[Qa.accessions.coor$DDLongitude < -93.5 & Qa.accessions.coor$DDLatitude > 33.5] <- 5
points(Qa.accessions.coor$DDLongitude[geographic.region==5], Qa.accessions.coor$DDLatitude[geographic.region==5], type = "p", pch = 21, col = "green", cex=3, lwd=2)

#geographic region 8
windows(14,14)
plot(States, xlim=c(-94.5, -89.2), ylim=c(31.5, 34.2), axes=F, box=T, mar=c(6,6,2.1,1))
axis(1, line=-4.9, seq(-94, -90), labels=T, cex.axis=2)
axis(2, line=0, seq(32, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-1)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud[Qa.occurrence$DDLatitude>31.4], Qa.occurrence$DDLatitude[Qa.occurrence$DDLatitude>31.4], type = "p", pch = 19, col = "blue")
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 24, col = "red", cex=2)
abline(v = -90, lty=3, lwd=2)
geographic.region[Qa.accessions.coor$DDLongitude > -90] <- 8
points(Qa.accessions.coor$DDLongitude[geographic.region==8], Qa.accessions.coor$DDLatitude[geographic.region==8], type = "p", pch = 21, col = "blue", cex=3)

#examine result
geographic.region
length(geographic.region)

#visualize the result using a map of the northwestern side of the range of the species,
#showing all collection localities of the accessions in the commmon garden.
#define color of each accession accordng to the respective geographic region
accession.by.geographic.region.col <- rep(NA, times=length(geographic.region))
accession.by.geographic.region.col[geographic.region==2] <- "#D81B60"
accession.by.geographic.region.col[geographic.region==3] <- "#1E88E5"
accession.by.geographic.region.col[geographic.region==4] <- "#FFC107"
accession.by.geographic.region.col[geographic.region==5] <- "#42FD7F"
accession.by.geographic.region.col[geographic.region==8] <- "#E7A0FD"
#create map
plot(States, xlim=c(-94.5, -89.2), ylim=c(31.5, 34.2), axes=F, box=T, mar=c(6,6,2.1,1))
axis(1, line=-4.9, seq(-94, -90), labels=T, cex.axis=2)
axis(2, line=0, seq(32, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-1)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud[Qa.occurrence$DDLatitude>31.4], Qa.occurrence$DDLatitude[Qa.occurrence$DDLatitude>31.4], type = "p", pch = 19, col = "gray")
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 17, col = accession.by.geographic.region.col, cex=2)
text(-92.2, 32.9, "2", cex=1.5, col="#D81B60")
text(-94.2, 33, "3", cex=1.5, col="#1E88E5")
text(-92.8, 33.7, "4", cex=1.5, col="#FFC107")
text(-93.5, 33.8, "5", cex=1.5, col="#42FD7F")
text(-89.5, 31.8, "8", cex=1.5, col="#E7A0FD")

#visualize the result using a map of the whole range of the species
windows(14,14)
plot(States, xlim=c(-95, -81.5), ylim=c(30, 34.2), axes=F, box=T, mar=c(3.1,6,2.1,1))
axis(1, line=-12.2, seq(-94, -82, 1), labels=F, cex.axis=2)
axis(1, line=-12.2, seq(-94, -82, 2), labels=T, cex.axis=2)
axis(1, line=-12.2, cex.axis=2)
axis(2, line=0, seq(30, 34, 1), cex.axis=2)
mtext("Longitude (decimal degrees)", side=1, cex=2, line=-9)
mtext("Latitude (decimal degrees)", side=2, cex=2, line=3.7)
points(Qa.occurrence$DDLongitud[Qa.occurrence$DDLatitude>31.4], Qa.occurrence$DDLatitude[Qa.occurrence$DDLatitude>31.4], type = "p", pch = 19, col = "gray")
points(Qa.accessions.coor$DDLongitude, Qa.accessions.coor$DDLatitude, type = "p", pch = 17, col = accession.by.geographic.region.col, cex=2)
text(-91.9, 32.9, "2", cex=1.5, col="#D81B60")
text(-94.5, 33, "3", cex=1.5, col="#1E88E5")
text(-92.5, 33.7, "4", cex=1.5, col="#FFC107")
text(-93.5, 34.1, "5", cex=1.5, col="#42FD7F")
text(-89.5, 32, "8", cex=1.5, col="#E7A0FD")

#############################################################################################
# 5. Save provenance data as a text file named "QaAccCoor.csv".
#############################################################################################

#create data frame with provenance data for each accession
Qa.prov.data <- data.frame(Qa.accessions, Qa.accessions.coor, geographic.region)
colnames(Qa.prov.data)[1] <- "AccessionNumber"
colnames(Qa.prov.data)[4] <- "GeographicRegion"
#examine the result
dim(Qa.prov.data)
str(Qa.prov.data)
head(Qa.prov.data)

write.csv(Qa.prov.data, file=paste("QaAccCoor_", format(Sys.time(), "%Y%b%d_%H%M%S"), ".csv", sep=""), row.names = F)


