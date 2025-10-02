# Script to generate the 10 sqkm resolution grid to avoid the "gird is empty changing translation" problem in feems

library(fs)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(tmaptools)
library(tmap)
library(dggridR)

world <- ne_countries(
  scale = "medium",   # "small", "medium", or "large"
  returnclass = "sf"  # return as an sf object
)

# Save to shapefile
st_write(world, "world.shp", append = FALSE)

grid_resolution <- 10

dggs <- dgconstruct(
  res = grid_resolution,
  precision = 30,
  projection = "ISEA",
  aperture = 4,
  topology = "TRIANGLE"
)

hisp_grid_file_name_as_written <- dgshptogrid(
  dggs, "world.shp",
  cellsize = 0.01, savegrid = "world_10resolution.shp"
)
