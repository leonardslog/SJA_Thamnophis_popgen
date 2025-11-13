library(terra)
library(scales)

# SJA plot area
x <- c(-123.5, -122.5)
y <- c(48.4, 48.8)
sja_ext <- ext(-122.9, -122.5, 48.4, 48.8)

# conversion to m for albers projection
sja_ext_m <- ext(-2e+06, -1.8e+06, 1.2e+06, 1.5e+06)
inset_crs <- "EPSG:4326"

north_america <- vect("~/Documents/sji_thamnophis/NorthAmerica_shape/North_America.shp")

crop_ext <- c(-123, -59, 23, 54)
sja_region <- crop(north_america, sja_ext)

north_america_inset <- north_america
test <- terra::project(north_america_inset, "EPSG:102005")


plot(test)
plot(test, xlim=c(-130, -50),ylim=c(15, 60),xlab="Latitude", ylab="Longitude")


# albers equal area projection for inset
na_albers <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"


# Reproject 
projected_na <- project(north_america, na_albers)
projected_icesheet <- project(ice_sheet, na_albers)
# projected_sja <- project(sja_region, na_albers)

# Plot the projected data for inset
plot(projected_na, xlim=c(-1.88e+06, -1.83e+06), ylim=c(1.25e+06, 1.3e+06))

sja_ext_m <- ext(-1.88e+06, -1.83e+06, 1.25e+06, 1.3e+06)

north_america_inset <- ext(-2.5e+06, 3e+06, -2e+06, 2e+06)
crop_test <- crop(x=projected_na, y=north_america_inset, ext=TRUE)

sja_box_m <- vect(ext(-1.875e+06, -1.83e+06, 1.25e+06, 1.3e+06))

plot(projected_na, xlim=c(-2.5e+06, -1e+06), ylim=c(5e+05, 2e+06))
plot(sja_box_m, border="red", lwd=2, add=TRUE)
plot(projected_icesheet, col="blue", lwd=2)



main_ext <- ext(-2e+06, -1.8e+06, 1.1e+06, 1.4e+06)

plot(sja_box_m, border="red", lwd=2, add=TRUE)
plot(projected_icesheet, col=alpha("lightblue",alpha = 0.4),border=NA, add=TRUE)
plot(grat, col="lightgray", lab.loc=c(1,2), add=TRUE)


# lat lon graticule
grat <- graticule(lon=seq(-125, -122, by=.2), lat=seq(48, 49, by=.2), crs=na_albers)
cropped_grat <- crop(grat, sja_extent)

terra::plot(projected_na, ext=sja_extent, buffer=FALSE, axes=FALSE, box=TRUE)
terra::plot(cropped_grat, col="lightgray", lab.loc=c(1,2), lty=3, 
            add=TRUE, lab.lat=1:4, lab.lon=2:10, retro=FALSE, tickmarks=TRUE)
# north(xy="bottomleft", type=1, label="N", angle=0, head=0.1, xpd=TRUE)
# sbar(d=50, xy="bottomleft", type="bar", divs=2, below="km", labels = c(0,25,50),
#      lonlat=NULL, lwd=1.5, xpd=TRUE, ticks=FALSE)


inset(crop_test, offset = 0, loc="bottomleft", box=sja_extent,
      pbox=list(col="red", lwd=.5),
      scale=0.45, add=TRUE) # offset required to get rid of inset's inner margin 


# olympic peninsula, vancouver island, and puget sound extent
plot(projected_na, xlim=c(-2e+06, -1.75e+06), ylim=c(1.1e+06, 1.4e+06))

########################
# San Juan Archipelago #
########################

adj_extent = ext(-123.25, -122.73, 48.4, 48.8)

pdf(file = "sja_sampling_site.pdf", 
    width = 2.63,       
    height = 2.95,
    onefile="onefile",
    bg="white",
    pointsize = 10)

terra::plot(north_america, ext=adj_extent, box=TRUE, oma=c(0,0,0,0),
            lwd=0.5, mar=c(4, 3, 1, 0)+0.15, xaxt="n", yaxt="n", axes=FALSE)
axis(side=2, outer=FALSE, lwd=0.5, line=0, cex.axis=0.7, padj=1)
axis(side=1, outer=FALSE, at=c(-123.2,-123,-122.8), cex.axis=0.7, padj=-1.5, 
     lwd=0.5, line=0.95)
mtext("Longitude", side = 1, line = 2.5)
mtext("Latitude", side = 2, line = 2)
dev.off()

##################################
# Broader extent with glaciation #
##################################

# from https://github.com/awickert/North-American-Ice-Sheets
ice_sheet <- vect("North-American-Ice-Sheets-master/WGS84/ice015000/ice015000.shp")
sja_ext <- ext(-123.3, -122.5, 48.4, 48.8)

pdf(file = "sja_location_ice_sheet.pdf", 
    width = 2.63,       
    height = 2.95,
    onefile="onefile",
    bg="white",
    pointsize = 10)

plot(north_america, xlim=c(-125, -121), ylim=c(47, 50), oma=c(0,0,0,0),
     col='darkgrey',border="black",lwd=0.5,
     mar=c(2,1,0,3),
     pax=list(side=c(1,4)))
mtext("Longitude", side = 1, line = 3) 
mtext("Latitude", side = 4, line = 1)
plot(ice_sheet, col=alpha("white",alpha = 0.25), lty="dotted", 
     border="black", add=TRUE)
plot(sja_ext, lwd=1.5, border="black", add=TRUE)
north(xy="bottomleft", type=1, label="N", angle=0, head=0.1, xpd=TRUE)
sbar(d=100, xy=c(-122.5,49.7), type="bar", divs=2, below="km", 
     labels = c(0,50,100), cex=0.8, lonlat=NULL, lwd=1, xpd=TRUE, ticks=FALSE)
dev.off()


################################
# SJA version 2 - sampling map #
################################

## read sampling data

specimens <- read.csv("sample_coords.csv")

ele <- subset(specimens,species=="ele")
ord <- subset(specimens,species=="ord")
sir <- subset(specimens,species=="sir")

inset_ext <- ext(-125, -121, 47, 50)
test_inset <- crop(north_america, inset_ext)
ice_sheet_inset <- crop(ice_sheet, inset_ext)

pdf(file = "test_samplingmap_2.pdf", 
    width = 2.92,
    height = 3.55,
    onefile="onefile",
    bg="white",
    pointsize = 10)
terra::plot(north_america, ext=adj_extent, box=TRUE, ylab="Latitude", 
            oma=c(0,0,0,0), mar=c(2,2,0,0),
            xlab="Longitude", lwd=0.5, axes=TRUE)
# points(ele[,5],ele[,4], pch=21,bg=alpha("#ffff99",0.6), cex=1.2)
# points(sir[,5],sir[,4], pch=21,bg=alpha("#1f78b4",0.6), cex=1.2)
# points(ord[,5],ord[,4], pch=21,bg=alpha("#d73027",0.6), cex=1.2)
terra::inset(test_inset,
             box = adj_extent,
             background="white", scale = 0.25, offset = 0,
             col = "darkgrey", border = NA, loc = "bottomleft", add = TRUE)
# terra::inset(ice_sheet_inset, col=alpha("white",alpha = 0.25), perimeter=FALSE, 
#              background=NA, lty="dotted", offset = 0, border="black", 
#              loc = "bottomleft", add=TRUE)
dev.off()








terra::inset(ice_sheet, col=alpha("white",alpha = 0.25), lty="dotted", 
                  border="black", add=TRUE)
plot(sja_ext, lwd=1.5, border="black", add=TRUE)
north(xy="bottomleft", type=1, label="N", angle=0, head=0.1, xpd=TRUE)
sbar(d=100, xy=c(-122.5,49.7), type="bar", divs=2, below="km", 
     labels = c(0,50,100), cex=0.8, lonlat=NULL, lwd=1, xpd=TRUE, ticks=FALSE)
