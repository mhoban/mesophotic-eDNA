library(tidyverse)
library(sf)
library(maps)
library(mapdata)
library(here)
library(cowplot)

sample_data <- read_csv(here::here("data","sample_data.csv"))
spots <- sample_data %>% 
  filter(project == "mesophotic") %>%
  distinct(station_grouping,lat,lon,.keep_all = T) %>%
  st_as_sf(coords=c("lon","lat"),crs=4326,remove=F)

box <- st_bbox(spots)


pbm <- cbind(c(box['xmin'],box['xmax'],box['xmax'],box['xmin'],box['xmin']),
                 c(box['ymin'],box['ymin'],box['ymax'],box['ymax'],box['ymin']))
boxpoly <- st_polygon(list(pbm))
boxsfc <- st_sfc(boxpoly,crs=4326)
plotbox <- st_sf(name="study area",geometry=boxsfc)

hbox <- c('xmin'=-156.178, 'xmax'=-154.687, 'ymin'=18.846, 'ymax'=20.334)
world <- st_as_sf(map("worldHires",plot=F,fill=T),crs=4326)
hawaii <- world %>% filter(ID == "Hawaii")
x <- 0.01
ggplot() + 
  geom_sf(data=hawaii, fill='gray20',color="lightgrey",size=0.07) +
  geom_sf(data=spots, shape=25,size=4,color="black",fill="firebrick") +
  coord_sf(xlim = c(box['xmin']-x,box['xmax']+x),
           ylim = c(box['ymin']-x,box['ymax']+x))

ggplot() + 
  geom_sf(data=hawaii, fill='gray20',color="lightgrey",size=0.07) +
  geom_sf(data=plotbox, color="red", fill="transparent") + 
  coord_sf(xlim = c(hbox['xmin']-x,hbox['xmax']+x),
           ylim = c(hbox['ymin']-x,hbox['ymax']+x)) +
  theme_void()
  
