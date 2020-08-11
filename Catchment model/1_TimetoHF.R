#Script to produce travel time rasters for MDG using friction surface  
setwd(paste0(Sys.getenv("HOME"), "/Prevalence-Madagascar/Catchment model"))
# clear workspace
rm(list = ls())

## Required Packages
require(gdistance)
library(raster)
library(rgdal)

#pull in polygon to get extent
poly <- readOGR(dsn = "/home/suzanne/Shapefiles/Malareo_District.shp", 
                layer = "Malareo_District")
e <- extent(poly)

# Input Files
friction.surface.filename <- "friction_surface_2015_v1.tif"


#  Define the spatial information from the friction surface
friction <- raster(friction.surface.filename)
fs1 <- crop(friction, e)
plot(fs1)
## Read in the points table.
load("../Incidence/sample_incidence_data.RData")
n_hf <- dim(hf_coords)[1]
points <- cbind(1:n_hf, hf_coords)    ## only contains those validated (2801)
names(points)[1] <- "FID"
names(points)[2] <- "X"
names(points)[3] <- "Y"

head(points)

#Loop through points
HF_list <- unique(points$FID)
filepath <- ("MDG_time_travel_raster/") 

for (i in seq_along(HF_list)) { 

  output.filename <- paste(filepath, HF_list[i], "HF.access.tif", sep='')
  
  T <- transition(fs1, function(x) 1/mean(x), 8) # RAM intensive, can be very slow for large areas
  T.GC <- geoCorrection(T)                    
  
  HF.coords <- c(points$X[i],points$Y[i])
  HF.raster <- accCost(T.GC, HF.coords)
  
  writeRaster(HF.raster, output.filename)
  print(output.filename)
}
