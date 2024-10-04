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
## 1. "QaAccCoor_2024Jul24_104558.csv" [Latitude and longitude (i.e., provenance) data for collected accessions]
## 2. "OccurrenceDatabaseSubset2023-06-21.txt" [Latitude and longitude data for all occurrences of Quercus arkansana]
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
#use ggplot2 to build layers of a map - base layer (p1)
# + states (p2) + points (p3)...etc.
library(tidyverse)
library(sf)
library(RColorBrewer)
library(viridis)
library(tmap)
library(ggpattern)

#optional for aesthetics
#ggrepel and showtext provide fancy formatting for text
#font_add_google("Lato",
               #regular.wt = 300,
               #bold.wt = 700)

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
  #first layer is the land mass outline as a grey border (fancy)
  geom_sf(data = land,
          fill = NA,
          color = "grey50",
          fill = "#dfdfdf",
          lwd = 1) +
  #crop the whole thing to size
  coord_sf(xlim = c(-100,-75),
           ylim = c(25,45))
print(p1)
#add states and slightly thinner line width
p2 <- p1 +
  geom_sf(data = states,
          color = "white",
          fill = "#dfdfdf",
          size = 0.2) +
  #crop the whole thing to size
    #we can zoom further in later
  coord_sf(xlim = c(-100,-75),
           ylim = c(28,45))
print(p2)

#add Quercus arkansas occurrence points
p3 <- p2 +
  geom_point(data = Qa.occurrence,
             aes(DDLongitude,DDLatitude),
             colour = "grey20",
             fill = "grey50",
             shape = 21,
             alpha = 0.6,
             size = 2
                ) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf(xlim = c(-97,-80),
           ylim = c(29,36))
#add Quercus arkansas collection points with provenance "GeographicRegion"  
p4 <- p3 +
  geom_point(data = Qa.prov.data,
             aes(DDLongitude, DDLatitude,
                 fill = factor(GeographicRegion)),
                 color = "grey20",
             shape = 24,
             size = 3) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_discrete(name = "Region") +
    scale_color_manual(values = "grey20", guide = "none")
print(p4)  

#save the plot
ggsave("QAProvOverOccurrences.png", dpi = 300)

#zoom in on provenances
p5 <- p4 + coord_sf(xlim = c(-96,-88.75),
              ylim = c(30.5,34.5))
print(p5)

#save the plot
ggsave("QAProvOverOccurrences_B.png", dpi = 300)

########################################################
## Option with buffers ##
########################################################

#add 50km buffers to occurrence points
#shaded buffers around provenance points
#lighter buffers around occurrence data


# Convert occurrence and prov data to sf objects, if they aren't already
# The remove = FALSE argument in st_as_sf() ensures that columns like 
# GeographicRegion are kept in the new sf object.

Qa.occurrence_sf <- st_as_sf(Qa.occurrence, coords = c("DDLongitude", "DDLatitude"), crs = 4326, remove = FALSE)
Qa.prov.data_sf <- st_as_sf(Qa.prov.data, coords = c("DDLongitude", "DDLatitude"), crs = 4326, remove = FALSE)

# Buffer the points (50 km buffer, transform to an appropriate projection for distance calculation)
Qa.occurrence_buffer <- st_transform(Qa.occurrence_sf, crs = 3395) %>%  # Use a suitable projection for accurate distance
  st_buffer(dist = 50000) %>%    # 50 km buffer
  st_transform(crs = 4326)       # Back to lat/long

Qa.prov_buffer <- st_transform(Qa.prov.data_sf, crs = 3395) %>% 
  st_buffer(dist = 50000) %>% 
  st_transform(crs = 4326)

#start building the map
p4_with_pattern_buffers <- p4 +
  # First, add the buffer for Qa.occurrence with lighter transparency (below the points)
  geom_sf(data = Qa.occurrence_buffer, 
          fill = "grey80",    # Light, subdued fill color
          color = NA,         # No border for buffer
          alpha = 0.3) +      # Transparent
  
  # Then, add the buffer for Qa.prov.data with a diagonal stroke pattern
  geom_sf_pattern(data = Qa.prov_buffer,
                  pattern = "stripe",   # Set pattern to diagonal stripes
                  pattern_angle = 45,   # Diagonal stroke at 45 degrees
                  pattern_density = 0.025, # Density of stripes
                  pattern_spacing = 0.02, # Adjust spacing of the stripes
                  pattern_fill = NA,    # No fill color for pattern
                  pattern_color = "grey65", # Color of the pattern
                  fill = NA,            # No solid fill
                  color = "grey65",       # Border outline color
                  linetype = "solid",  # Keep the outline
                  alpha = 0.5,          # Semi-transparent
                  size = 0.5) +         # Adjust size if necessary
  
  # Add the occurrence points (Qa.occurrence) on top of the buffers
  geom_point(data = Qa.occurrence,
             aes(DDLongitude, DDLatitude),
             colour = "grey20",
             fill = "grey50",
             shape = 21,
             alpha = 0.6,
             size = 2) +
  
  # Add the Qa.prov.data points on top of the buffers
  geom_point(data = Qa.prov.data,
             aes(DDLongitude, DDLatitude,
                 fill = factor(GeographicRegion)),
             color = "grey20",
             shape = 24,
             size = 4) +
  
  # Add axis labels and the fill legend for regions
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_discrete(name = "Region") + # add legend for scaled fill of the GeographicRegion points
  scale_color_manual(values = "grey20", guide = "none") +  # Suppress outline color legend
  coord_sf(xlim = c(-97,-80),
           ylim = c(29,36))
  
# Display the plot
print(p4_with_pattern_buffers)

#save the plot
ggsave("QAProvOverOccurrences_wBuffer.png", dpi = 300)

#zoom in on provenances
p5_with_buffers <- p4_with_pattern_buffers + coord_sf(xlim = c(-96,-88.75),
                    ylim = c(30.5,34.5))
print(p5_with_buffers)

#save the plot
ggsave("QAProvOverOccurrences_B_wBuffer.png", dpi = 300)


##############################################################
# Map edits #
# merge buffer circles
# make each region a different shape
# make occurrence points not filled
# try to make transparency not stacked
# remove gray grid background and replace with blank or ocean blue
# make lat long axes more balanced
# can there be a zoomed in map of the # of maternal lines to show concentrated of trees where they were collected?
##############################################################
#Oct 3, 2024

#start over with a different background on the map

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
             fill = NA,
             shape = 21,
             
             size = 2
  ) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf(xlim = c(-96,-81),
           ylim = c(29,35)) 
  # + scale_x_continuous(breaks = seq(-97, -80, by = 1))
print(p3)

#add Quercus arkansas collection points with provenance "GeographicRegion" 
#have colors and different shapes corresponding to Regions
p4 <- p3 +
  geom_point(data = Qa.prov.data,
             aes(DDLongitude, DDLatitude,
                 fill = factor(GeographicRegion),
                 shape = factor(GeographicRegion)),
             color = "grey20",
             size = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_discrete(name = "Region") +
  scale_shape_manual(name = "Region",
                     values = c(21, 22, 23, 24, 25)) + #shapes
  scale_color_manual(values = "grey20", guide = "none") +
  theme(legend.key = element_rect(fill = NA))
print(p4)  

#save the plot
ggsave("v2QAProvOverOccurrences.png", dpi = 300)

#zoom in on provenances
p5 <- p4 + coord_sf(xlim = c(-95,-89),
                    ylim = c(30,35))
print(p5)

#save the plot
ggsave("v2QAProvOverOccurrences_zoomed.png", dpi = 300)

###########################################################################

# add buffers but make them so the circle outlines don't overlap, they connect, by using st_union

###########################################################################

# Buffer the points (50 km buffer, transform to an appropriate projection for distance calculation)
Qa.occurrence_buffer <- st_transform(Qa.occurrence_sf, crs = 3395) %>%  # Use a suitable projection for accurate distance
  st_buffer(dist = 50000) %>%    # 50 km buffer
  st_transform(crs = 4326)       # Back to lat/long

Qa.prov_buffer <- st_transform(Qa.prov.data_sf, crs = 3395) %>% 
  st_buffer(dist = 50000) %>% 
  st_transform(crs = 4326)

# Combine overlapping buffers for Qa.occurrence
Qa.occurrence_buffer_union <- st_union(Qa.occurrence_buffer)

# Combine overlapping buffers for Qa.prov.data
Qa.prov_buffer_union <- st_union(Qa.prov_buffer)

#start building the map
p4_with_pattern_buffers <- p4 +
  # First, add the unioned buffer for Qa.occurrence with lighter transparency
  geom_sf(data = Qa.occurrence_buffer_union, 
          fill = "grey80",    # Light, subdued fill color
          color = NA,         # No border for buffer
          alpha = 0.3) +      # Transparent
  # Then, add the unioned buffer for Qa.prov.data with a diagonal stroke pattern
  geom_sf_pattern(data = Qa.prov_buffer_union,
                  pattern = "stripe",   # Set pattern to diagonal stripes
                  pattern_angle = 45,   # Diagonal stroke at 45 degrees
                  pattern_density = 0.025, # Density of stripes
                  pattern_spacing = 0.02, # Adjust spacing of the stripes
                  pattern_fill = NA,    # No fill color for pattern
                  pattern_color = "grey65", # Color of the pattern
                  fill = NA,            # No solid fill
                  color = "grey65",       # Border outline color
                  linetype = "solid",  # Keep the outline
                  alpha = 0.5,          # Semi-transparent
                  size = 0.5)   +        # Adjust size if necessary

# Add the occurrence points (Qa.occurrence) on top of the buffers
  geom_point(data = Qa.occurrence,
             aes(DDLongitude,DDLatitude),
             colour = "grey70",
             fill = NA,
             shape = 21,
             
             size = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf(xlim = c(-96,-81),
           ylim = c(29,35)) +
  #add Quercus arkansas collection points with provenance "GeographicRegion" 
  #have colors and different shapes corresponding to Regions
  geom_point(data = Qa.prov.data,
             aes(DDLongitude, DDLatitude,
                 fill = factor(GeographicRegion), #colors based on region
                 shape = factor(GeographicRegion)), #shapes based on region
             color = "grey20",
             size = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_discrete(name = "Region") +
  scale_shape_manual(name = "Region",
                     values = c(21, 22, 23, 24, 25)) + #shapes
  scale_color_manual(values = "grey20", guide = "none") +
  theme(legend.key = element_rect(fill = NA))
print(p4_with_pattern_buffers) 
#save the plot
ggsave("v2QAProvOverOccurrences_buffer_join.png", dpi = 300)
#zoom in on provenances
p5 <- p4_with_pattern_buffers + coord_sf(xlim = c(-95,-89),
                    ylim = c(30,35))
#save the plot
print(p5)
ggsave("v2QAProvOverOccurrences_buffer_join_zoomed.png", dpi = 300)


################ Oct 4, 2024 ##########################
# note: Ivan uses these colors in other plot
# "#D81B60", "#1E88E5", "#FFC107", "#42FD7F", "#E7A0FD" 
#######################################################

