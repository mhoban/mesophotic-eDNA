library(tidyverse)
library(sf)
library(maps)
library(mapdata)
library(here)
library(cowplot)
library(beyonce)
library(ggrepel)

sample_data <- read_csv(here::here("data","sample_data.csv"))
spots <- sample_data %>% 
  filter(
    project == "mesophotic",
    !is.na(lon) & !is.na(lat)
  ) %>%
  distinct(station_grouping,lat,lon,.keep_all = T) %>% 
  # let's cheat the point over a bit to make it easier to see
  mutate(
    lat = case_when(
      str_detect(station,"Outhouse") ~  19.636009,
      TRUE ~ lat
    ),
    lon = case_when(
      str_detect(station,"Outhouse") ~  -156.005123,
      TRUE ~ lon
    )
  ) %>%
  st_as_sf(coords=c("lon","lat"),crs=4326,remove=F)

box <- st_bbox(spots)
box['xmin'] <- box['xmin']-0.2
box['xmax'] <- box['xmax']+0.2
box['ymin'] <- box['ymin'] - 0.05
box['ymax'] <- box['ymax'] + 0.05


box <- c(ymax=20.337941429661853, xmin=-156.62139587608522, ymin=18.84548284767367, xmax=-155.4117057257134)


colors <- beyonce_palette(54)
colors <- beyonce_palette(2)

pbm <- cbind(c(box['xmin'],box['xmax'],box['xmax'],box['xmin'],box['xmin']),
                 c(box['ymin'],box['ymin'],box['ymax'],box['ymax'],box['ymin']))
boxpoly <- st_polygon(list(pbm))
boxsfc <- st_sfc(boxpoly,crs=4326)
plotbox <- st_sf(name="study area",geometry=boxsfc)


hbox <- c("ymin"=18.63314470051002, "xmax"=-154.6469259428272, "ymax"=22.434801224940056, "xmin"=-160.76622913125863)
world <- st_as_sf(map("worldHires",plot=F,fill=T),crs=4326)
hawaii <- world %>% filter(ID == "Hawaii")
x <- 0.01


if (TRUE) {
  mainmap <- ggplot() + 
  geom_sf(data=hawaii, fill="grey57",color="black",size=0.5) +
  # geom_sf(data=hawaii, fill=colors[5],color="black",size=0.5) +
  geom_sf(data=spots, shape=25,size=6,color="black",fill="firebrick") +
  coord_sf(xlim = c(box['xmin']-x,box['xmax']+x),
           ylim = c(box['ymin']-x,box['ymax']+x)) +
  # geom_sf_label(data=spots,aes(label=station_grouping),position="dodge") + 
  geom_label_repel(
    data=spots,
    aes(label=station_grouping,geometry=geometry),
    point.padding = 12,
    stat = "sf_coordinates",
    min.segment.length = 0
  ) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  theme(
    axis.text = element_text(color="black", size=14),
    # panel.background = element_rect(fill=colors[2]) ,
    panel.background = element_rect(fill="white") ,
    plot.background = element_rect(fill="white",color=NA),
    # text=element_text(color="grey45")      
    text=element_text(color="black")      
  ) #+
  # xlab("Longitude (decimal degrees)") + 
  # ylab("Latitude (decimal degrees)")
mainmap

inset <- ggplot() + 
  geom_sf(data=hawaii, fill="grey57",color="black",size=0.1) +
  geom_sf(data=plotbox, color="red", fill="transparent") + 
  coord_sf(xlim = c(hbox['xmin']-x,hbox['xmax']+x),
           ylim = c(hbox['ymin']-x,hbox['ymax']+x)) +
  theme_void() +
  theme(
    panel.background = element_rect(fill="white", color="black") 
  )

# box2 <- c(xmin = -156.156018, ymax =   19.734519, xmax = -155.945894, ymin =   19.608708)
# inset2 <- ggplot() + 
#   geom_sf(data=hawaii,fill="grey57",color="black",size=0.1) + 
#   geom_sf(data=spots, shape=25,size=6,color="black",fill="firebrick") +
#   coord_sf(xlim = c(box2['xmin']-x,box2['xmax']+x),
#            ylim = c(box2['ymin']-x,box2['ymax']+x)) +
#   theme_void() +
#   theme(
#     panel.background = element_rect(fill="white", color="black") 
#   )
# inset2
}
# inset
# inset
(
  p <- ggdraw(mainmap) + 
    draw_plot(
      inset,
      x = -0.130,
      y = 0.35,
      scale=0.4
    ) 
)

ggsave(p,file=here("output/map.svg"),device="svg",width = 8.5, height = 11, units = "in")
ggsave(p,file=here("output/map.pdf"),device=cairo_pdf(),width = 8.5, height = 11, units = "in")
  