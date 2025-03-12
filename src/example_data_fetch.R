library(tidyverse)
library(tigris)
library(sf)
library(here)

# each census tract in mult co
ct <- tigris::tracts(state = "OR", county = "Multnomah", cb=TRUE)
place <- tigris::places(state = "OR", cb=TRUE)

pdx <- place %>% filter(NAME == "Portland")

# filter to just in portland, only centroids of tracts
ct <- ct %>% st_centroid %>% mutate(
  in_portland = as.integer(st_intersects(ct, pdx))
)
pdx_tracts <- ct %>% filter(in_portland == 1)

# write to /data

st_write(pdx_tracts, paste(here(),"data","pdx_tracts","pdx_tracts.shp", sep = "/"))