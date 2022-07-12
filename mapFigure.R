################################  FILE LICENSE  ################################
#
#	This file is copyright (C) 2022 Matthew Kosnik and Yvette Bauder
#
#	This program is free software; you can redistribute it and/or modify it 
#	under the terms of version 3 the GNU General Public License as published 
#	by the Free Software Foundation.
#
#	This program is distributed in the hope that it will be useful, but WITHOUT
#	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
#	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
#	more details.
#
#	To view a copy of the license go to:
#	http://www.gnu.org/licenses/licenses.html#GPL
#	
################################################################################

## USE THE SAME PLOT PAR ACROSS ALL PLOTS
source('plotPreferences.R')
fontFamily <- 'Helvetica'	

############################  Packages needed as well as base R  #################
{my_packages <- c("jpeg")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)}

## NEEDED TO INCLUDE A JPEG IN A PLOT
library(jpeg) 

## LOAD THE GOOGLE MAPS SATELITE IMAGE INTO R
mapImage <- readJPEG("./data/MapRawBits/reefMap2.jpg") 
mapSize <- dim(mapImage)

## SET PAGE HEIGHT SO THAT THE MAP RATIO IS CORRECT
pageHeight <-  (pageWidthTwo / mapSize[2])* mapSize[1]

## COORDINATES TO DOUBLE CHECK THAT IT IS WORKING!
## SYKES - IMAGE SPLICE CORNER:-23.439425° x 152.004504°
##	Wistari - image edge @ ripple: -23.467471° x 151.841413°
##	Wistari - image dark corner @ front edge: -23.485113° x 151.888326°
##	ORTRS accom bloc: -23.506982° x 152.091938°
##	Heron Jetty: -23.442115° x 151.910359°
##	Core: -23.49677 ° x 152.06587 °
testpoints <- data.frame(x=c(152.004504, 151.841413, 151.888326, 152.091938, 151.910359, 152.06587), y=c(-23.439425, -23.467471, -23.485113, -23.506982, -23.442115, -23.49677))

## GET SITE COORINATES FROM FILE
gp <- read.csv('./data/Mamo_Water_Depth.csv')
gp$y <- gp$lat
gp$x <- gp$lng

## OPEN A PDF FILE...
pdf("./figs/Fig2-map.pdf", width=pageWidthTwo, height= pageHeight, family=fontFamily, paper=pagePaper)

## SET PLOTTING PARAMETERS
par(oma=c(0,0,0,0), mar=c(0,0,0,0), bty='u', yaxs='i', xaxs='i', mgp=c(2,0.6,0), cex=0.75)

## SET PLOTTING AREA IN LAT / LONG
xLim <- c(151.8365, 152.1033)
yLim <- c(-23.5237, -23.4081)

## MAKE OFFSETS FOR LABELS
x10 <- xLim[1] + (xLim[2] - xLim[1])*0.042
y10 <- yLim[1] + (yLim[2] - yLim[1])*0.05
x20 <- xLim[1] + (xLim[2] - xLim[1])*0.08
y20 <- yLim[1] + (yLim[2] - yLim[1])*0.08

## PLOT MAP AREA
plot(1:1, type='n', ann=FALSE, axes=FALSE, ylim=yLim, yaxs='i', xlim=xLim, xaxs='i', las=1, lwd=0.5)

##	map image
rasterImage(mapImage, xleft=xLim[1], ybottom=yLim[1], xright= xLim[2], ytop=yLim[2])

cexReef <- 1.2
cexLabs <- 0.9

##	ADD: lat / lng lines
segments(x20,-23.5, xLim[2],-23.5, col='white', lwd=1.2)
segments(xLim[1],-23.441666666666667, xLim[2],-23.441666666666667, col='white', lwd=1.2, lty=2)
segments(152, y20,152, yLim[2], col='white', lwd=1.2)

##	ADD: lat / lng labels
text(151.9, y10, '151.9°', col='white', cex= cexLabs)
text(152.0, y10, '152.0°', col='white', cex= cexLabs)
text(x10, -23.45, '23.45°', col='white', cex= cexLabs)
text(x10, -23.50, '23.50°', col='white', cex= cexLabs)
text(152.05, -23.441666666666667, 'Tropic of Capricorn', col='white', cex= cexLabs, pos=3, offset=0.2)

# ADD: SCALE LABEL
text(151.92, -23.511, '5km', col='white', cex= cexLabs)

## ADD: lat ticks
par(mgp=c(3,1,0))
lat <- seq(-24, -22 ,by=0.05)
axis(2, at=lat, labels=NA, tck=0.01, col='white', lwd=1.2)

## ADD: lng ticks
lng <- seq(150, 153, by=0.05)
axis(1, at=lng, tck=0.01, col='white', col.lab='white', col.sub='white', lwd=1.2)

## ADD: REEF LABELS
#text(151.86, -23.457, 'Wistari', col='white', cex= cexReef)
text(151.9, -23.49, 'Wistari', col='white', cex= cexReef)
text(151.97, -23.429, 'Heron', col='white', cex= cexReef)
text(152.06, -23.480, 'One Tree', col='white', cex= cexReef)

gp$Site2 <- substring(gp$Site,2,nchar(gp$Site))
gp$Site2[27] <- 'C'

# ADD: site points & labels
#points(gp$Lng, gp$Lat, pch=21, col='black', bg='orange', lwd=0.6, cex=1.2)
#text(gp$Lng, gp$Lat, gp$Site2, col='black', bg='orange', lwd=0.5, cex=0.45, offset=0.0)

pGrab <- gp[1:26,]
points(pGrab$Lng, pGrab$Lat, pch=22, col='black', bg='orange', lwd=0.6, cex=1.4)
text(pGrab$Lng, pGrab$Lat, pGrab$Site2, col='black', bg='orange', lwd=0.5, cex=0.45, offset=0.0)

pCore <- gp[27,]
points(pCore $Lng, pCore $Lat, pch=21, col='black', bg='yellow', lwd=0.6, cex=1.3)
text(pCore $Lng, pCore $Lat, pCore $Site2, col='black', bg='yellow', lwd=0.5, cex=0.5, offset=0.0)


dev.off()
