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

################################################################################
## Yvette Bauder, Briony Mamo, Glenn A. Brock, Matthew A. Kosnik. 2022 (In revision). 
## One Tree Reef Foraminifera: a relic of the pre-colonial Great Barrier Reef
## Conservation Palaeobiology of Marine Ecosystems: Concepts and Applications
## Geological Society of London Special Publication.
##
## Contact: Matthew Kosnik, Matthew.Kosnik@mq.edu.au or mkosnik@alumni.uchicago.edu
################################################################################

# SET PATH TO LOCATION OF FILES
source('./plotPreferences.R')

# FILES CONTAINING DATA
FILE.Pb <- paste(PATH,'data','oti-lead.csv',sep='/')
FILE.Sed <- paste(PATH,'data','All_sed_results.csv',sep='/')

## LOAD DATA
dataPb <- read.table(FILE.Pb, header=TRUE, sep='\t', na.strings=c("\\N","NA"))

datedPb210 <- dataPb$top_cm

foramLayers <- read.csv('./data/OTR_Core_pb210_CIC_Ages.csv')
a <- lm(fit ~ (top_mm), data= foramLayers)
k <- coefficients(a)
keyAges <- c(1945,1788)
keyDepths <- (keyAges - k[1])/k[2]
axisAges <- c(2000,1900,1800,1700,1600)
axisDepths <- (axisAges - k[1])/k[2]
axisAges <- seq(1500,2050,by=50)
axisDepths <- (axisAges - k[1])/k[2]


dataSed <- read.csv(FILE.Sed, header=TRUE, na.strings=c("\\N","NA"))

dataSed <- dataSed[(dataSed[,'sample_type']!='Average'),]
dataSed <- dataSed[!is.na(dataSed[,'sample_top']),]
dataSed <- dataSed[(dataSed[,'Sample_Name']!='Kos62'),]

cols <- colnames(dataSed)
grains <- cols[33:length(cols)]
grainsNum <- as.numeric(substr(grains,2,20))

sedData <- dataSed[,grains]
cumData <- dataSed[,grains]

for (i in 2:length(grains)) {	
	cumData[,(i)] <- cumData[,(i-1)]+ cumData[,(i)]
}

cumData <- cbind(dataSed[,'sample_top'],cumData)
sedData <- cbind(dataSed[,'sample_top'],sedData)

startCol <-35
endCol <- 89

grainWindow <- grainsNum[startCol:endCol]

plot(grainWindow, cumData[i,(startCol:endCol)], type='n', log='x', ylim=c(160,-100))
for (i in 1:nrow(cumData))
	lines(grainWindow,(cumData[i,1]-(cumData[i,(startCol:endCol)])))

volmean <- cbind(dataSed[,c('Sample_Name','D..4..3....Volume.weighted.mean')],dataSed[,'sample_top'])
volord <- sort(dataSed[,'sample_top'])

plot(grainWindow, sedData[i,(startCol:endCol)], type='n', log='x', ylim=c(170,-60), ylab='core depth (cm)', xlab='grain size (um)')
for (i in 1:nrow(sedData))
	lines(grainWindow,(sedData[i,1]-(10*sedData[i,(startCol:endCol)])), col='grey')
points(dataSed[,'D..4..3....Volume.weighted.mean'],dataSed[,'sample_top'], lty=1, lwd=2)

tops <- sort(unique(dataSed[,'sample_top']))

a <- matrix(nrow=length(tops),ncol=11)
a[,1] <- tops
rownames(a) <- tops
for (d in tops) {
	dataTop <- dataSed[(dataSed['sample_top']==d),]
	a[as.character(d),2] <- nrow(dataTop)
	a[as.character(d),3] <- min(dataTop[,'d..0.1.']) 
	a[as.character(d),4] <- median(dataTop[,'d..0.1.']) 
	a[as.character(d),5] <- max(dataTop[,'d..0.1.']) 
	a[as.character(d),6] <- min(dataTop[,'d..0.5.']) 
	a[as.character(d),7] <- median(dataTop[,'d..0.5.']) 
	a[as.character(d),8] <- max(dataTop[,'d..0.5.']) 
	a[as.character(d),9] <- min(dataTop[,'d..0.9.']) 
	a[as.character(d),10] <- median(dataTop[,'d..0.9.']) 
	a[as.character(d),11] <- max(dataTop[,'d..0.9.']) 
}
mid <- a[,1]-1
a <- cbind(a,mid)

xmin <- min(dataSed[,'d..0.1.'])
xmax <- max(dataSed[,'d..0.9.'])
mm <- quantile(a[,7], c(0.025,0.5,0.975))
ml <- quantile(a[,4], c(0.025,0.5,0.975))
mh <- quantile(a[,10], c(0.025,0.5,0.975))
dataSed[,'sample_top'] <- dataSed[,'sample_top']-2

#######################################################################################################
pdf(paste(PATH.FIG,'Fig02-CoreSed&Chron.pdf',sep='/'), width=pageWidthTwo, height=pageHeight/2)

colLine <- 'white'
colText <- 'white'
colPoly <- 'grey50'

colLine <- 'black'
colText <- 'black'
colOutline <- 'white'
colPoly <- 'grey75'
cexDefault <- 0.8


par (mfcol=c(1,2), cex=cexDefault, fg=colLine, mar=c(3,3.25,1,0), mgp=c(2.2,0.7,0), oma=c(1,0.5,1,1))

ytop <- -2
ybot <- 180

plot(dataSed[,'d..0.5.'],dataSed[,'sample_top'], log='x', ylim=c(ybot,ytop), xlim=c(xmin,xmax), ylab='Hand core depth (m)', xlab=expression(paste("Grain size (diameter (", mu, "m))")), type='n', main='', yaxs='i', yaxt='n', col=colLine, col.axis=colLine,col.sub=colLine,col.main=colText,col.lab=colText)
axis(2,at=seq(0,160, by=20),labels=formatC(seq(0,1.6, by=0.20),format='f',digits=1),col=colLine, col.axis=colText, las=1)

mtext('A.',adj=0, line=0.5)
mtext('Grain size', line=0.5)

#atop <- 0
#abot <- 22
#polygon(c(1,1,2000,2000), c(atop,abot, abot, atop), col= colModTr, lty=0)

x0 <- rep(1,2)
x1 <- rep(2000,2)

segments(x0,keyDepths/10,x1, keyDepths/10, lty=2, col='blue')

#atop <- 27
#abot <- 80
#polygon(c(1,1,2000,2000), c(atop,abot, abot, atop), col= colColTr, lty=0)

#atop <- 90
#abot <- 180
#polygon(c(1,1,2000,2000), c(atop,abot, abot, atop), col= colPreTr, lty=0)



polygon(c(mm[1],mm[1],mm[3],mm[3]), c(ytop, ybot, ybot, ytop), col=colPoly, lty=0)
segments(mm[2], ytop, mm[2], ybot-10, lwd=2, lty=2)
text(mm[2], 175, round(mm[2],0), font=2)

polygon(c(ml[1], ml[1], ml[3], ml[3]),c(ytop, ybot, ybot, ytop), col=colPoly, lty=0)
segments(ml[2], ytop, ml[2], ybot-10, lwd=1, lty=2)
text(ml[2], 175, round(ml[2],0), font=2)

polygon(c(mh[1], mh[1], mh[3], mh[3]), c(ytop, ybot, ybot, ytop), col=colPoly, lty=0)
segments(mh[2], ytop, mh[2], ybot-10, lwd=1, lty=2)
text(mh[2], 175, round(mh[2],0), font=2)

segments(a[,4],a[,'mid'],a[,10],a[,'mid'], col='grey')

points(dataSed[,'d..0.1.'], dataSed[,'sample_top']+1, cex=0.4, pch=19, col=colLine)
points(dataSed[,'d..0.5.'], dataSed[,'sample_top']+1, cex=0.4, pch=19, col=colLine)
points(dataSed[,'d..0.9.'], dataSed[,'sample_top']+1, cex=0.4, pch=19, col=colLine)

points(a[,7],a[,'mid'])

#for (i in 1:length(datedLayers))
	#polygon(c(10,17,17,10),c(datedLayers[i],datedLayers[i],datedLayers[i]+5,datedLayers[i]+5), col=colText)

#for (i in 1:length(datedPb210))
#	polygon(c(10,16,16,10),c(datedPb210[i], datedPb210[i], datedPb210[i]+2, datedPb210[i]+2), col=colText)

fl <- foramLayers$top_mm/10
for (i in 1:length(fl)) {
	polygon(c(10,17,17,10),c(fl[i], fl[i], fl[i]+1, fl[i]+1), col=rgb(0,0,0,0.1), lwd=1)
}

pbDepths <- dataPb[,'top_cm']+0.5
pbAges <- dataPb[,'CICage']

axis(4, at=(axisDepths/10)+0.5, labels=rep("",length(axisDepths)), las=1)

#axis(4, at=pbDepths, labels=pbAges, las=1)
## LEAD

#axis(4, at=(foramLayers$top_mm/10)+0.5, labels=foramLayers$fit, las=1)

depthAdjust <- -1.5

plot(dataPb[,'Pb_total'],dataPb[,'top_cm']+depthAdjust, ylim=c(ybot,ytop), xlim=c(0.1,11), type='n', log='x', xlab=expression(""^{210} ~ "Pb (Bq/Kg)"), ylab='', yaxs='i', yaxt='n')
#axis(2,at=seq(0,60, by=10), labels=c('0.0','0.1','0.2','0.3','0.4','0.5','0.6'), las=1, col.axis='black')

mtext('B.',adj=0, line=0.5)
mtext('Lead activity', line=0.5)

#axis(2, at=(foramLayers$top_mm/10)+0.5, labels=rep('',9), las=1)

#atop <- 0
#abot <- 22
#yRng <- c(atop, abot, abot, atop)
xRng <- c(0.01,0.01,20,20)
#polygon(xRng, yRng, col= colModTr, lty=0)

#atop <- 27
#abot <- 80
#yRng <- c(atop, abot, abot, atop)
#polygon(xRng, yRng, col= colColTr, lty=0)

#atop <- 90
#abot <- 180
#yRng <- c(atop, abot, abot, atop)
#polygon(xRng, yRng, col= colPreTr, lty=0)

lines(dataPb[,'Pb_total'], dataPb[,'top_cm']+depthAdjust, lty=2,col='grey')
points(dataPb[,'Pb_total'], dataPb[,'top_cm']+depthAdjust, pch=21, col='grey', bg='grey')
segments(dataPb[,'Pb_total']+ dataPb[,'Pb_total_err'], dataPb[,'top_cm']+depthAdjust, dataPb[,'Pb_total']-dataPb[,'Pb_total_err'], dataPb[,'top_cm']+depthAdjust,col='grey')

points(dataPb[4:9,'Pb_unsupported'], dataPb[4:9,'top_cm']+depthAdjust, pch=21, col='black', bg='white', cex=1.2, lwd=0.4)

lines(dataPb[,'Pb_unsupported'], dataPb[,'top_cm']+depthAdjust, lty=2,col='black')
points(dataPb[,'Pb_unsupported'], dataPb[,'top_cm']+depthAdjust, pch=21, col='black', bg='black', cex=0.6)
segments(dataPb[,'Pb_unsupported']+ dataPb[,'Pb_unsupported_err'], dataPb[,'top_cm']+depthAdjust, dataPb[,'Pb_unsupported']-dataPb[,'Pb_unsupported_err'], dataPb[,'top_cm']+depthAdjust,col='black')
lines(c(0.01,0.3),c(52+depthAdjust,52+depthAdjust),col='black')



#legend(0.1,0,legend=c(expression(" "^{210} ~ "Pb"["total"]), expression(" "^{210} ~ "Pb"["unsupported"])),pch=21,col=c('grey','black'), lty=2, bty='n', cex=0.8)
legend(0.1,90,legend=c("Total","Unsupported"),pch=21,col=c('grey','black'),pt.bg=c('grey','black'), lty=2, bty='n')

#pbAge <- lm(CICage ~ top_cm, data=dataPb)
#abline(pbAge)


axis(2, at=(axisDepths/10)+0.5, labels= axisAges, las=1)
#axis(2, at=(keyDepths/10)+0.5, labels= keyAges, las=1)
#axis(2, at=(foramLayers$top_mm/10)+0.5, labels=foramLayers$fit, las=1)
#axis(2, at=(foramLayers$top_mm/10)+0.5, labels=rep("",9), las=1)


#mtext(paste0('age = ',round(k[2],2),' x depth(mm) + ',round(k[1])),side=1, line = -2, cex=0.8)

dev.off()
