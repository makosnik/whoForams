################################  FILE LICENSE  ################################
#
#	This file is copyright (C) 2022 Yvette Bauder and Matthew Kosnik
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

################################################################################
## Yvette Bauder, Briony Mamo, Glenn A. Brock, Matthew A. Kosnik. 2022 (In revision). 
## One Tree Reef Foraminifera: a relic of the pre-colonial Great Barrier Reef
## Conservation Palaeobiology of Marine Ecosystems: Concepts and Applications
## Geological Society of London Special Publication.
##
## Contact: Matthew Kosnik, Matthew.Kosnik@mq.edu.au or mkosnik@alumni.uchicago.edu
################################################################################

############################  About this script file  ##########################
## This script will set the preferences for the plots

############################  Packages needed as well as base R  #################
{my_packages <- c("dichromat","viridis")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)}

library(dichromat)
library(viridis)

##If there is no figure directory create one

PATH <- './'
PATH.OUT <- paste0(PATH,'output','/')
PATH.FIG <- paste0(PATH,'figs','/')

if(!paste(PATH.OUT) %in% list.dirs(PATH)) dir.create(PATH.OUT)#
if(!paste(PATH.FIG) %in% list.dirs(PATH)) dir.create(PATH.FIG)#

##Set page width, paper type and font
#pageWidthOne <- 3.386	
#pageWidthTwo <- 6.929		
#pageHeight <- 9.291		

## SET GEOSOC LONDON PAGE SIZES
pageWidthOne <- 65 / 25.4	
pageWidthTwo <- 135 / 25.4		
pageHeight <- 203 / 25.4		


pagePaper <- 'A4'		
fontFamily <- 'Times'	

##Set colours for time intervals 
colMod <- 'hotpink'
colCol <- 'green1'
colPre <- 'dodgerblue'

colModTr <- col2rgb(colMod)/255
colColTr <- col2rgb(colCol)/255
colPreTr <- col2rgb(colPre)/255

colModTr <- rgb(colModTr['red',1],colModTr['green',1],colModTr['blue',1],0.2)
colColTr <- rgb(colColTr['red',1],colColTr['green',1],colColTr['blue',1],0.2)
colPreTr <- rgb(colPreTr['red',1],colPreTr['green',1],colPreTr['blue',1],0.2)

## FROM AGE AGREEMENT PLOT
par(mar=c(3,3,1,1), las=0, mgp=c(2,0.75,0), cex=1.0)

intervalMeta <- data.frame(yrStart=c(1945, 1788, 1440),yrStop=c(2012,1945,1788),col=c(colMod,colCol,colPre),tCol=c(colModTr,colColTr,colPreTr))
intervalMeta$col <- as.character(intervalMeta$col)
intervalMeta$tCol <- as.character(intervalMeta$tCol)

timePolygons <- function(meta, yLim) {
  yLim[1] <- -10
  yLim[2] <- yLim[2]*1.04
  for (i in 1:nrow(meta)) {
    polygon(c(meta[i,1],meta[i,1],meta[i,2],meta[i,2]),c(yLim[1],yLim[2],yLim[2],yLim[1]),col=meta[i,'tCol'],lty=0)
  }
}


## PLOT PREFERENCES / STYLE
##Use where black is required in line graphs
col063 <- 'black'
col125 <- 'black'
col250 <- 'black'
fractionMeta <- data.frame(fraction=c('Total','063','125','250'), pch=c(5,0,6,2), col=c('black',col063,col125,col250), lty=1:4)
fractionMeta$name <- paste(fractionMeta$fraction,'\u03BCm')
fractionMeta$name <- paste(fractionMeta$fraction,'um')
fractionMeta[1,'name'] <- 'Total'

                        
                        

pCol <- c('skyblue','dodgerblue','navy','hotpink','forestgreen','red')
pCol <- c('skyblue','dodgerblue','navy','gold','forestgreen','red')
GroupNames <- c('pre-Colonial','Colonial','Modern','OneTree','Heron','Wistari')
#pch = 21: filled circle, 
#pch = 22: filled square, 
#pch = 23: filled diamond, 
#pch = 24: filled triangle point-up, 
#pch = 25: filled triangle point down. 

pCh <- c(21,21,21,23,24,25)

panelLabelSize <- 0.8
panelLabelSize <- 0.7
panelLabelLine <- -1.5
pValSide <- 3
pValAdj <- 1
pValLines <- c(-5,-4,-3)
pValLines <- c(0,-1,-2,-3)

pMeta <- data.frame(GroupNames= GroupNames, pCh=pCh, pCol=pCol)

                    
##########################################################################################
plotPointsWithErrors <- function (colName, panelLabelText, sumM,sumE,pCh,pCol,panelLabelSize, panelLabelLine, statSumText, pValSide, pValLines, pValAdj, yRange) {
	
	if (missing(yRange))
		yRange <- range(sumM[,colName]+sumE[,colName], sumM[,colName]-sumE[,colName])

	plot(sumM[,colName], ylim=yRange, axes=FALSE, ann=FALSE, pch=pCh, col=pCol, bg=pCol)
	box(bty='l')
	arrows(1:nrow(sumM), sumM[,colName]+sumE[,colName], 1:nrow(sumM), sumM[,colName]-sumE[,colName], code=3, angle=90, length=0.05, lwd=2, col= pCol)

	mtext(panelLabelText, side=1, line=panelLabelLine, adj=0.1, cex=panelLabelSize)
	axis(2)
	axis(1, labels=FALSE)

	if (is.data.frame(statSumText)) {
		mtext(statSumText[colName,'allGroups'], side=pValSide, line=pValLines[1], adj=pValAdj, cex=panelLabelSize)
		mtext(statSumText[colName,'sCollect'], side=pValSide, line=pValLines[2], adj=pValAdj, cex=panelLabelSize)
		mtext(statSumText[colName,'sTime'], side=pValSide, line=pValLines[3], adj=pValAdj, cex=panelLabelSize)
		mtext(statSumText[colName,'sReef'], side=pValSide, line=pValLines[4], adj=pValAdj, cex=panelLabelSize)
	}
}
      
      
plotBarWhisker <- function (dataFrame, panelLabelText, panelLabelSize = 1, panelLabelSide=3, panelLabelLine=-1, statSum, pValSide, pValLines, pValAdj, pCh, pCol, yRange) {
	
	dataFrame[(dataFrame$Groups == 'Pre-Colonial'),'Groups'] <- 'pre-Colonial'
	dataFrame$Groups <- factor(dataFrame$Groups, levels=GroupNames)
	
	if (missing(yRange))
		yRange <-range(dataFrame$ProSM)

	plot(dataFrame$Groups, dataFrame$ProSM, ylim=yRange, axes=FALSE, ann=FALSE, pch=pCh, col=pCol, bg=dataFrame$pCol)
	box(bty='l')
	points(dataFrame$Groups, dataFrame$ProSM, pch=dataFrame$pCh, bg=dataFrame$pCol, cex=0.8)
	
	mtext(panelLabelText, side=panelLabelSide, line= panelLabelLine, adj=0.1, cex=panelLabelSize)
	axis(2, las=1)
	axis(1, labels=FALSE)

#	statSum$p.adj <- paste('p =', round(statSum$p.adj,4))

	if (pValSide == 1)
		pValLines <- pValLines[c(4,3,2,1)] - 2

	comps <- c('aGrp','Type','Time','Reef')

	if (is.data.frame(statSum)) {
		for (i in 1:length(comps)) {
			if (statSum[(statSum$comp == comps[i]),'p.adj'] < 0.05) { 
				col <- 'black'
			} else {
				col <- 'grey'				
			}
			mtext(statSum[(statSum$comp == comps[i]),'p.txt'], side=pValSide, line=pValLines[i], adj=pValAdj, cex=panelLabelSize, col=col)
		}
	}
}



plotBarWhiskerPVB <- function (dataFrame, panelLabelText, panelLabelSize = 1, panelLabelSide=3, panelLabelLine=-1, statSum, pValSide, pValLines, pValAdj, pCh, pCol, yRange, addNtext=FALSE) {
	
	dataFrame[(dataFrame$Groups == 'Pre-Colonial'),'Groups'] <- 'pre-Colonial'
	dataFrame$Groups <- factor(dataFrame$Groups, levels=GroupNames)
	if (missing(addNtext))
		addNtext <- FALSE
	
#	if (missing(yRange))
	yRange <- range(dataFrame$ProSM)
	yRng <- (yRange[2]-yRange[1])
	yRange[1] <- yRange[1] - yRng*0.17
	yRange[2] <- yRange[2] + yRng*0.02
	if (addNtext == TRUE)
		yRange[2] <- yRange[2] + yRng*0.13
		
	plot(dataFrame$Groups, dataFrame$ProSM, ylim=yRange, axes=FALSE, ann=FALSE, pch=pCh, col=pCol, bg=dataFrame$pCol)
	box(bty='l')
	points(dataFrame$Groups, dataFrame$ProSM, pch=dataFrame$pCh, bg=dataFrame$pCol, cex=0.8)
	
	axis(2, las=1)
	axis(1, labels=FALSE)

#	statSum$p.adj <- paste('p =', round(statSum$p.adj,4))

	y1 <- c(yRange[1]+yRng*0.015, yRange[1]+yRng*0.055,yRange[1]+yRng*0.095,yRange[1]+yRng*0.095)

	comps <- c('aGrp','Type','Time','Reef')
	cX1 <- c(1,3,1,4) - 0.5
	cX2 <- c(6,4,3,6) + 0.5

	if (is.data.frame(statSum)) {
		for (i in 1:length(comps)) {
			if (statSum[(statSum$comp == comps[i]),'p.adj'] < 0.05) { 
				col <- 'black'
			} else {
				col <- 'grey'				
			}
			segments(cX1[i],y1[i],cX2[i],y1[i],lwd=2,col=col)
		}
	}
	if (panelLabelSide == 4) {
		mtext(panelLabelText, side=1, line=panelLabelLine, adj=0.9, cex=panelLabelSize)
	} else {
		mtext(panelLabelText, side=panelLabelSide, line=panelLabelLine, adj=0.1, cex=panelLabelSize)
	}
	
}
               
                        
