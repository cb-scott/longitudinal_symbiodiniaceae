
########################
# from https://seethedatablog.wordpress.com/2016/12/31/r-shift-centralprime-meridian-of-world-map/

library("data.table")
library("ggplot2")
library("maps")
library("maptools")
library("ggsn") #ggplot scalebars

# Load world map as map object
worldmap <- map("world", fill=TRUE, plot=FALSE)

# Convert map object to SpatialPolygons object (check ?map2SpatialPolygons)
WGS84 <- CRS("+proj=longlat +datum=WGS84")
worldmapPolys <- map2SpatialPolygons(worldmap, 
IDs=sapply(strsplit(worldmap$names, ":"), "[", 1L), 
proj4string=WGS84)

# shift central/prime meridian towards west – positive values only
shift <- 180 - 20

# transform map in a data table that ggplot can use
XY <- data.table(map_data(as(worldmapPolys, "SpatialPolygonsDataFrame")))

# Shift coordinates (data.table way)
XY[, long.new := long + shift]
XY[, long.new := ifelse(long.new > 180, long.new - 360, long.new)]
# split up coordinates of polygons that differ too much/across globe
XY[, to.split := sum(diff(long.new) > 300, na.rm=TRUE) > 0, by=group]
XY[, gr.split := ifelse(to.split & long.new < 0, paste0(group, ".", 1), group)]

# ADD FP's

library("rgdal")

dat<-readOGR("data/FaunalProvs_shp/FaunalProvinces_Global.shp")
coords <- spTransform(dat, WGS84)
FP <- data.table(map_data(as(coords, "SpatialPolygonsDataFrame")))
head(FP)
FP[, long.new := long + shift]
FP[, long.new := ifelse(long.new > 180, long.new - 360, long.new)]
FP[, to.split := sum(diff(long.new) > 300, na.rm=TRUE) > 0, by=group]
FP[, gr.split := ifelse(to.split & long.new < 0, paste0(group, ".", 1), group)]


head(FP)
FP <- FP[!FP$group %in% c(6,10,16,15),]
info <- read.csv("data/info.csv")
mapcols<-info$col 
names(mapcols)<-info$shp.n

# plot shifted map
map<-ggplot() + 
geom_polygon(data=FP, aes(x=long.new, y=lat, group=group,fill=as.factor(group)))+   
    geom_polygon(data=XY,  aes(x=long.new, y=lat, group=gr.split), 
colour="black", fill="grey65", size = 0.5)+
    geom_polygon(data=XY,  aes(x=long.new, y=lat, group=gr.split), 
colour="grey65", fill="grey65", size = 0.1)+
scale_fill_manual(values=mapcols)+
coord_map( ylim=c(-27, 35)) +
 scalebar(x.min = -110, x.max = -100, y.min = -28, y.max = -24, dist = 2000, dist_unit = "km",st.size=2,st.dist=0.5, box.fill=c("black", "white"),st.bottom = FALSE, transform = TRUE, model = "WGS84", border.size=0.1, height=0.2) +
guides(fill="none")+
theme_void()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1), plot.margin = unit(c(1,1,1,1), "mm"))
map


