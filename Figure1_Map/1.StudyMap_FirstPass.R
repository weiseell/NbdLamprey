#Map with all sampling sites colorized by collection time
#Created by: Ellie Weise

#setwd
setwd("~/OneDrive - Michigan State University/Documents/Sea_Lamprey_MS_project/exp_bioinformatics/Figures_Tables/")
#libraries
library(tidyverse)
library(sp)
library(sf)
library(rnaturalearth)
library(rgdal)
library(rgeos)
library(ggrepel)
library(cowplot)
#load data
df <- read.csv(file = "site_coordinates.csv")

#getting lake, river and state shapes
points_sf <- st_as_sf(df, coords = c("Long", "Lat"), crs = 4326)

lakes110 <- ne_download(scale = 110, type = 'lakes', category = 'physical', returnclass = 'sf')
statesus <- ne_states(country = 'united states of america', returnclass = 'sf')
statescan <- ne_states(country = 'canada', returnclass = 'sf')

#putting data into proper sf format
points_sf <- st_as_sf(df, coords = c("Long", "Lat"), crs = 4326)
points_sf <- cbind(points_sf,chYear_collect = as.character(df$Year_collect))
class(points_sf)
rnaturalearth::ne_states(geounit = "australia", returnclass = "sf")


#selecting specific lakes, states and provedences for map
lakes <- lakes110 %>% 
  filter(name_alt == "Great Lakes") %>% 
  #select(name) %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84")

USstates <- statesus %>% 
  #select(name) %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84")
  #filter(name %in% c('Michigan','Wisconsin','Ohio','New York','Pennsylvania'))

CANprov <- statescan %>% 
  #filter(name %in% c('Ontario')) %>% 
  select(name) %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84")

#map plot
ggplot()+
  geom_sf(data = lakes, fill = "white")+
  geom_sf(data = USstates)+
  geom_sf(data = CANprov)+
  geom_sf(data = points_sf,
          show.legend = "point",
          size = 2)+
  scale_fill_discrete(c("blue","light blue"))+
  coord_sf(xlim = c(-94, -74), ylim = c(40, 49), expand = FALSE)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(color = "Collection \nType",shape = "Year \nCollected")
#making a new plot with just the expansion location
#selecting specific lakes, states and provedences for map
lakes2 <- lakes110 %>% 
  filter(name_alt == "Great Lakes") %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84") %>% 
  filter(name %in% c('Lake Michigan','Lake Huron'))

USstates2 <- statesus %>% 
  #select(name) %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84") %>% 
  filter(name %in% c('Michigan'))

points_sf2 <- points_sf %>% 
  filter(Type == "Barrier")

inset <- ggplot()+
  geom_sf(data = lakes, fill = "white")+
  geom_sf(data = USstates)+
  geom_sf(data = CANprov)+
  geom_sf(data = points_sf2,
          aes(color = Type, shape=chYear_collect),
          size = 2, color = "black")+
  coord_sf(xlim = c(-91, -79), ylim = c(41.5, 47.6), expand = FALSE)+
  geom_rect(aes(xmin=-84.7,xmax=-83.5,ymin=45,ymax=46),
            fill="transparent",color = "black")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))

#creating a smaller map with just the study locations (zoomed in)
#adding river shape files
river_shapes <- st_read(dsn = "USA_Rivers_and_Streams-shp/")
lake_shapes <- st_read(dsn = "USA_Detailed_Water_Bodies/")
#need the following object IDs

#BMR: 56604
bmr <- river_shapes %>% filter(OBJECTID==56604)
bmr$Name <- as.factor("Black Mallard River")
#BMR Lake: 294985
bmr_lake <- lake_shapes %>% filter(OBJECTID==294985)
#CHE:
che <- river_shapes %>% filter(OBJECTID==58322)
#OCQ: 58137, 58508,58561
ocq <- river_shapes %>% filter(OBJECTID==58137|OBJECTID==58508|OBJECTID==58561)
river_shapes %>% filter(Name == "Ocqueoc River")
#OCQ Lake: 22370
ocq_lake <- lake_shapes %>% filter(OBJECTID==22370)

rivers <- rbind(bmr,che,ocq)

#MI Rivers
MI_rivers <- river_shapes %>% filter(State=="MI")
MI_lakes <- lake_shapes %>% filter(State=="MI")

p2 <- ggplot()+
  #geom_sf(data = lakes, fill = "white")+
  geom_sf(data = USstates)+
  geom_sf(data = CANprov)+
  geom_sf(data = MI_rivers, col = "lightblue")+
  geom_sf(data = rivers,col = "black")+
  geom_sf(data = bmr_lake, col = "black")+
  geom_sf(data = ocq_lake, col = "white")+
  coord_sf(xlim = c(-84.7, -83.5), ylim = c(45, 46), expand = FALSE)+
  theme_classic()

tiff(filename = "StudyMap.tiff",height = 6,width = 5,units = "in",res = 200)
ggdraw()+
  draw_plot(p2)+
  draw_plot(inset,x = 0.50, y = 0.55, width = 0.5, height = 0.5)+
  draw_plot_label(label = "Black \nMallard",x=0.37,y=0.57,size = 8)+
  draw_plot_label(label = "Ocqueoc",x=0.55,y=0.45,size = 8)+
  draw_plot_label(label = "Pigeon",x=0.2,y=0.45,size = 8)
dev.off()


