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

############################  About this script file  ##########################
#	This script consolidates all the stats mentioned in the manuscript into
#	one place.


TAXONLEVEL
DROP63
COMBINEFRACTIONS
DROPUNIDENTIFIED
DROPJUVENILE
DROPRARETAXA

# SAMPLE SUMMARY
( bob <- summary(as.factor(substring(rownames(dataToUseT),0,1))) )

## Number of taxa
nrow(dataCore)		# taxa in core
nrow(dataMamo)		# taxa in grabs
dataMamoOTR <- (dataMamo[,which((substring(colnames(dataMamo),0,1) == 'O'))])
dataMamoOTR <- dataMamoOTR[(rowSums(dataMamoOTR) > 0),]
nrow(dataMamoOTR)	# taxa in grabs - at OTR
nrow(dataToUse)		# taxa total dataset

## Number of specimens
sum(data63[7:ncol(data63)], na.rm=TRUE)  # Core w/o 63um, w/ juv. and unident
sum(dataCore[,-1])	# core (only ID > 125um)
sum(dataMamo[,-1])	# grabs
sum(dataToUse[,-1])	# whole dataset

## Number of taxa with 1 or 2 specimens
length(which(rowSums(dataCore[,-1]) < 3)) / length(rowSums(dataCore[,-1]))
length(which(rowSums(dataMamo[,-1]) < 3)) / length(rowSums(dataMamo[,-1]))
taxonSums <- rowSums(dataToUse[,-1])
length(which(taxonSums == 1))
length(which(taxonSums == 2))
length(which(taxonSums < 3)) / length(taxonSums)

## Stats table...
#Sum.Table

## STATS by category table...
# sColl, sTime, sReef
statSumText
statSum

## RARE
mamoSum <- data.frame(taxon=dataMamo[,1],sum=rowSums(dataMamo[,-1]))
mamoTot <- sum(mamoSum$sum)
( mamo1pc <- mamoTot * 0.01 )
mamoRare <- (mamoSum[(mamoSum$sum < mamo1pc),])
( mamoProRare <- nrow(mamoRare)/nrow(mamoSum) )

coreSum <- data.frame(taxon=dataCore[,1],sum=rowSums(dataCore[,-1]))
coreTot <- sum(coreSum $sum)
( core1pc <- coreTot * 0.01 )
coreRare <- (coreSum[(coreSum$sum < core1pc),])
( coreProRare <- nrow(coreRare)/nrow(coreSum) )

makeCategoryProportions <- function(fsiMeta, fsiColumn, dataToUse, TAXONLEVEL) {

	sumAll <- colSums(dataToUse[,-1])
	
	taxSM.Y <- fsiMeta[(fsiMeta[,fsiColumn] == TRUE),TAXONLEVEL]
	matSM.Y <- dataToUse[(dataToUse$Group.1 %in% taxSM.Y),]
	sumSM.Y <- colSums(matSM.Y[,-1])
	proSM.Y <- sumSM.Y / sumAll
	
	dfSM.Y <- data.frame(Groups = meta$grpNames, ProSM = proSM.Y*100)
	return(dfSM.Y)
	
}

makeCategoryAggregate <- function(proSM.Y, TAXONLEVEL) {
	proSM.m <- aggregate(proSM.Y[,2], by=list(proSM.Y[,1]), mean)
	colnames(proSM.m) <- c(TAXONLEVEL,'mean')
	proSM.e <- aggregate(proSM.Y[,2], by=list(proSM.Y[,1]), sd)
	colnames(proSM.e) <- c(TAXONLEVEL,'sd')
		
	catSum <- merge(proSM.m, proSM.e)
	catSum <- catSum[c(5,1,3,4,2,6),]
	catSum$mean <- round(catSum$mean,3)
	catSum$sd <- round(catSum$sd,3)
	
	return(catSum)
}





##  FSI
# FSI

plot(FSI ~ as.factor(meta$grpNames))

kwFSI <- metaCompareKW(FSI, meta, 'FSI')
(kwFSI <- kwFSI[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )


##  MILIOLIDS
# fsiMeta$smallMiliolid
MiliolidRaw <- makeCategoryProportions(fsiMeta, 'smallMiliolid', dataToUse, TAXONLEVEL)
( MiliolidAgg <- makeCategoryAggregate(MiliolidRaw, TAXONLEVEL) )

plot(ProSM ~ as.factor(Groups), data = MiliolidRaw)

kwMil <- metaCompareKW(MiliolidRaw$ProSM, meta, 'Miliolid')
(kwMil <- kwMil[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )

##  INFAUNAL
# fsiMeta$infaunal
infaunalRaw <- makeCategoryProportions(fsiMeta, 'infaunal', dataToUse, TAXONLEVEL)
infaunalAgg <- makeCategoryAggregate(infaunalRaw, TAXONLEVEL)

plot(ProSM ~ as.factor(Groups), data = infaunalRaw)

kwInf <- metaCompareKW(infaunalRaw$ProSM, meta, 'infaunal')
(kwInf <- kwInf[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )


##  LBF
# fsiMeta$LBF
lbfRaw <- makeCategoryProportions(fsiMeta, 'LBF', dataToUse, TAXONLEVEL)
lbfAgg <- makeCategoryAggregate(lbfRaw, TAXONLEVEL)

plot(ProSM ~ as.factor(Groups), data = lbfRaw)

kwLBF <- metaCompareKW(lbfRaw$ProSM, meta, 'LBF')
(kwLBF <- kwLBF[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )


##  Bolivinids
# fsiMeta$Bolivinids
bolivinidRaw <- makeCategoryProportions(fsiMeta, 'Bolivinids', dataToUse, TAXONLEVEL)
bolivinidAgg <- makeCategoryAggregate(bolivinidRaw, TAXONLEVEL)

plot(ProSM ~ as.factor(Groups), data = bolivinidRaw)

kwBol <- metaCompareKW(bolivinidRaw$ProSM, meta, 'Bolivinids')
(kwBol <- kwBol[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )

##  Agglutinated
# fsiMeta$Agglutinated
agglutRaw <- makeCategoryProportions(fsiMeta, 'Agglutinated', dataToUse, TAXONLEVEL)
agglutAgg <- makeCategoryAggregate(agglutRaw, TAXONLEVEL)

plot(ProSM ~ as.factor(Groups), data = agglutRaw)

kwAgg <- metaCompareKW(agglutRaw$ProSM, meta, 'Agglutinated')
(kwAgg <- kwAgg[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )

kwSum2 <- as.data.frame(rbind(kwFSI, kwAgg, kwBol, kwInf, kwLBF, kwMil))

write.csv(kwSum2, file=paste0(PATH.OUT, 'Table2-Fig3-', TAXONLEVEL, DROPRARETAXA,'-KW-summary.csv'))

## p.value = NA means that they are found at sites in all groups.
isaGrps$sign[is.na(isaGrps$sign$p.value),]

##	A == positive predictive value == probability of site in group if species found.
##	B == fidelity or sensitivity == probability of finding taxon in sites belonging to the group

summary(isaGrps, indvalcomp=TRUE, alpha=1)
summary(isaType, indvalcomp=TRUE)
summary(isaTime, indvalcomp=TRUE)
summary(isaReef, indvalcomp=TRUE)


pdf(file=paste0(PATH.FIG,'Fig04-', DROPRARETAXA,'-traits.pdf'), width=pageWidthTwo, height=pageHeight)
par(mfrow=c(5,2), oma=c(6,3,1,1), mar=c(1,1,1,2), las=1)

FSI2 <- data.frame(Groups=meta$grpNames,ProSM=FSI,pCh=meta$pCh, pCol=meta$pCol)
letterline <- -0.5
par(mar=c(1,1,1,2))
plotBarWhiskerPVB(FSI2,"Foram Stress Index", panelLabelSize, panelLabelSide=4, panelLabelLine=-2.5, kwFSI, pValSide=1, pValLines, pValAdj, pCh, pCol, yRange=c(8.2,10), addNtext=TRUE)
mtext('A',side=3, line=letterline, cex=1.2, adj=0.05)

nv <- aggregate(FSI2$Groups,by=list(FSI2$Groups),FUN=length)
mv <- aggregate(FSI2$ProSM,by=list(FSI2$Groups),FUN=max)
sq <- c(5,1,3,4,2,6)
text(x=1:6,y=mv$x[sq],paste0('n=',nv$x[sq]), pos=3)

par(mar=c(1,2,1,1))
plotBarWhiskerPVB(agglutRaw,"% Agglutinated", panelLabelSize, panelLabelSide=3, panelLabelLine=-2, kwAgg, pValSide, pValLines, pValAdj, pCh, pCol, yRange=c(0,5))
mtext('B',side=3, line= letterline, cex=1.2, adj=0.05)

par(mar=c(1,1,1,2))
plotBarWhiskerPVB(bolivinidRaw,"% Bolivinid", panelLabelSize, panelLabelSide=1, panelLabelLine=-2.5, kwBol, pValSide, pValLines, pValAdj, pCh, pCol, yRange=c(0,20))
mtext('C',side=3, line= letterline, cex=1.2, adj=0.05)

par(mar=c(1,2,1,1))
plotBarWhiskerPVB(infaunalRaw,"% Infaunal", panelLabelSize, panelLabelSide=1, panelLabelLine=-2.5, kwInf, pValSide=1, pValLines, pValAdj, pCh, pCol, yRange=c(45,95))
mtext('D',side=3, line= letterline, cex=1.2, adj=0.05)

par(mar=c(1,1,1,2))
plotBarWhiskerPVB(lbfRaw,"% Larger benthic forams", panelLabelSize, panelLabelSide=3, panelLabelLine=-2, kwLBF, pValSide, pValLines, pValAdj, pCh, pCol, yRange=c(0,40))
mtext('E',side=3, line= letterline, cex=1.2, adj=0.05)

axis(1, at=1:6, labels=GroupNames, las=2)
mtext('OTR Core', side=1, line=6.5, adj=0.14, cex=0.8)
mtext('Surface', side=1, line=5.3, adj=0.82, cex=0.8)
mtext('Grabs', side=1, line=6.5, adj=0.79, cex=0.8)



par(mar=c(1,2,1,1))
plotBarWhiskerPVB(MiliolidRaw,"% Miliolid", panelLabelSize, panelLabelSide=3, panelLabelLine=-2, kwMil, pValSide=1, pValLines, pValAdj, pCh, pCol, yRange=c(45,85))
mtext('F',side=3, line= letterline, cex=1.2, adj=0.05)

axis(1, at=1:6, labels=GroupNames, las=2)
mtext('OTR Core', side=1, line=6.5, adj=0.14, cex=0.8)
mtext('Surface', side=1, line=5.3, adj=0.82, cex=0.8)
mtext('Grabs', side=1, line=6.5, adj=0.79, cex=0.8)

box('outer')
dev.off()
