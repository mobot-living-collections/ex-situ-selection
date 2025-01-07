######################################################################################################
######################################################################################################
##
## TITLE
##
## Preparation and mapping of occurrence and provenance data for Quercus arkansana.
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
## 1. [Latitude and longitude (i.e., provenance) data for collected accessions]
## 2. [Latitude and longitude data for all occurrences of Quercus arkansana]
##
## CONTENTS
##
## 1. Read and examine data files
## 2. Extract occurrence data for Q. arkansas
## 3. Download rnaturalearth vector data for base map
## 4. Build maps with ggplot by layering provenance data over occurrence data
##
##
######################################################################################################
######################################################################################################

library(rnaturalearth)
library(raster)
library(ggrepel)
library(showtext)
library(ggplot2)
library(ggspatial)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(viridis)
library(tmap)
library(ggpattern)

#read in data files
occurrence <- read.table("OccurrenceDatabaseSubset2023-06-21.txt", header = T, sep = "\t")
Qa.prov.data <- read.table("QaAccCoor_2024Jul24_104558.csv", header = T, sep =",")

#extract occurrence data for Q. arkansana
unique(occurrence$Taxon)
Qa.occurrence <- occurrence[occurrence$Taxon == "Quercus arkansana",]

#read vector data using rnaturalearth
land <- ne_download(
  scale = 50,
  type = "land",
  category = "physical",
  returnclass = "sf"
)
states <- ne_download(type = 'states',
                      scale = 50,
                      returnclass = 'sf')

#start building the map with ggplot2
p1 <- ggplot() +
  #first layer is the land mass outline
  geom_sf(data = land,
          fill = "#dfdfdf", # land fill
          color = "grey50", # border color
          lwd = 1) +
  #crop the whole thing to size
  coord_sf(xlim = c(-100,-75),
           ylim = c(25,45)) +
  # make the background ocean blue and remove grid lines
  theme(panel.background = element_rect(fill = "#a6cee3"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(p1)

#add states with slightly thinner line width
p2 <- p1 +
  geom_sf(data = states,
          color = "white",
          fill = "#dfdfdf",
          size = 0.2) +
  #crop the whole thing to size
  coord_sf(xlim = c(-100,-75),
           ylim = c(28,45))

print(p2)

#add Quercus arkansas occurrence points with no fill and transposition
p3 <- p2 +
  geom_point(data = Qa.occurrence,
             aes(DDLongitude,DDLatitude),
             colour = "grey70",
             fill = "grey80",
             shape = 21,
             alpha = 0.5,
             size = 2
  ) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf(xlim = c(-96,-81),
           ylim = c(29,35))
print(p3)

#add Quercus arkansas collection points with provenance "GeographicRegion"
#color and shape by Region

p4 <- p3 +
  geom_point(data = Qa.prov.data,
             aes(DDLongitude, DDLatitude,
                 fill = factor(GeographicRegion),
                 shape = factor(GeographicRegion)),
             color = "grey20",
             size = 3) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_viridis(name = expression("Region, " * italic(r)),
                     discrete = TRUE,
                     option = "D",
                     labels = c("2 (15,106)",
                                "3 (4,8)",
                                "4 (9,50)",
                                "5 (13,73)",
                                "8 (1,6)")) +
  scale_shape_manual(name = expression("Region, " * italic(r)),
                     values = c(21, 22, 23, 24, 25),
                     labels = c("2 (15,106)",
                                "3 (4,8)",
                                "4 (9,50)",
                                "5 (13,73)",
                                "8 (1,6)")) +

  scale_color_manual(values = "grey20", guide = "none") +

  theme(legend.key = element_rect(fill = NA),
        legend.background = element_rect(color = "gray60", fill = NA),
        legend.title = element_text(hjust = 0.5)) +
  #scale bar
  annotation_scale(location = "br",
                   width_hint = 0.2,
                   text_col = "gray40",
                   style = "ticks",
                   pad_x = unit(3, "cm"))+

  #compass
  annotation_north_arrow(location = "br",
                         which_north = "true",
                         style = north_arrow_fancy_orienteering(
                           fill = c("white", "gray70"),
                           line_col = "gray40",
                           text_size = 7,
                           text_col = "gray40"),
                         pad_y = unit(1.7, "cm"),
                         pad_x = unit(0.03, "cm"),
                         height = unit(0.9, "cm"),
                         width = unit(0.9, "cm")
                          )

print(p4)

#save the plot
ggsave("v8_viridis_QAProvOverOccurrences.png", dpi = 300)
