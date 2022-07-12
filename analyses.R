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

##	DATA IMPORT / IMPLEMENT CLEANING / SUBSETTING / GROUPING
##########################################################################################

{my_packages <- c("indicspecies","vegan")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)}

library(indicspecies)
library(vegan)


##	TAXONOMIC LEVEL TO USE FOR ANALYSES
TAXONLEVEL <- 'Family'
TAXONLEVEL <- 'Genus'

TAXONLEVELPlural <- 'Genera'
if (TAXONLEVEL == 'Family')
	TAXONLEVELPlural <- 'Families'


##	HOW TO SUBSET / COMBINE DATA
DROP63 <- TRUE
COMBINEFRACTIONS <- TRUE
DROPUNIDENTIFIED <- TRUE
DROPJUVENILE <- TRUE

##	DROP RARE TAXA FROM ANALYSES
##	0 means nothing is dropped
##	otherwise a numerical value to be multiplied by 1%
DROPRARETAXA <- 0.2
DROPRARETAXA <- 2
DROPRARETAXA <- 0.5
DROPRARETAXA <- 1
DROPRARETAXA <- 0.1
DROPRARETAXA <- 0

source('plotPreferences.R')

##	SCRIPT TO IMPORT / IMPLEMENT CLEANING / SUBSETTING / GROUPING
##	dataCore == data from sediment core
##	dataMamo == data from surface samples
##	dataToUse == the merged dataset 
source('./dataImport.R')

#tail(dataCore)
#tail(dataMamo)
#tail(dataToUse)

#tail(dataCoreT)
#tail(dataMamoT)
#tail(dataToUseT)


## STATS META FUNCTION
#################
statsMe <- function(Sum.Table,Categories) {
	
	keep <- which(!is.na(Categories))
	
	Sum.Table <- Sum.Table[,keep]
	Categories <- Categories[keep]
		
	p <- rep(NA,length=nrow(Sum.Table))
	c <- rep(NA,length=nrow(Sum.Table))
	for (r in 1:nrow(Sum.Table)) {
		kt <- kruskal.test(Sum.Table[r,] ~ Categories)
		p[r] <- kt$p.value
		c[r] <- kt$statistic
	}
	statsSum <- data.frame(var= rownames(Sum.Table), p.value = p, p.adj = NA, chisq = round(c,2))
	return(statsSum)
}


adjust.p.value <- function(statSum, row) {
	
	pCols <- which(colnames(statSum) %in% "p.value")
	pVec <- as.vector(statSum[row,pCols], mode='double')

	aCols <- which(colnames(statSum) %in% "p.adj")
	statSum[row,aCols] <- p.adjust(pVec,'holm')
	
	return(statSum)
}

##	RUN 4 KRUSKAL-WALLIS / WILCOXON'S RANK SUM TESTS AND ADJUST P VALUES
metaCompareKW <- function(metricData, meta, metaName,sigDigits=2) {

	if (is.data.frame(metricData)) {
		metricData <- metricData[,2]
	}

	kwG <- kruskal.test(metricData, meta$grpNames)
	kwC <- kruskal.test(metricData, meta$compareCollect)
	kwT <- kruskal.test(metricData, meta$compareTime)
	kwR <- kruskal.test(metricData, meta$compareReef)

	kwMil <- list(kwG,kwC,kwT,kwR)
	kwMil <- as.data.frame(do.call(rbind,kwMil))
	kwMil$p.adj <- round(p.adjust(kwMil$p.value), sigDigits)
	kwMil$p.value <- round(unlist(kwMil$p.value), sigDigits)
	kwMil$comp <- c('aGrp','Type','Time','Reef')
	kwMil$metaName <- rep(metaName,4)
	kwMil$statistic <- unlist(kwMil$statistic)
	kwMil$parameter <- unlist(kwMil$parameter)
	kwMil$p.txt <- paste('p =', kwMil$p.adj)
	kwMil[(kwMil$p.txt == 'p = 0'),'p.txt'] <- paste('p <',1/(10^sigDigits))
	kwMil[(kwMil$p.txt == 'p = 1'),'p.txt'] <- 'p = 1.00'
	
	kwMil[(nchar(kwMil$p.txt) < (6 + sigDigits)),'p.txt'] <- paste0(kwMil[(nchar(kwMil$p.txt) < (6 + sigDigits)),'p.txt'],0)
	
	return(kwMil)	
}



##	CALCULATE DIVERSITY METRICS / INDICES
##########################################################################################

Chao.Table <- apply(dataToUseT, 1, estimateR)
Chao.Table <- Chao.Table[1:2,]

##	Richness
SObs <- data.frame(Groups=meta$grpNames, ProSM=Chao.Table['S.obs',], pCh=meta$pCh, pCol= meta$pCol)
kwObs <- metaCompareKW(SObs, meta, 'Observed Richness')
(kwObs <- kwObs[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )

## Chao 
Chao <- data.frame(Groups=meta$grpNames, ProSM=Chao.Table['S.chao1',], pCh= meta$pCh, pCol= meta$pCol)
kwCha <- metaCompareKW(Chao, meta, 'Chao1 Richness')
(kwCha <- kwCha[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )

mean(Chao$ProSM - SObs$ProSM)
sd(Chao$ProSM - SObs$ProSM)

## Shannon 
Shannon <- diversity(dataToUseT, index='shannon')
kwShn <- metaCompareKW(Shannon, meta, 'Shannon')
(kwShn <- kwShn[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )

# Shannon equitability index - 0 is perfectly even,  1 is uneven.
SEq <- data.frame(Groups=meta$grpNames, ProSM=(Shannon/log(Chao.Table[1,])), pCh= meta$pCh, pCol= meta$pCol)
kwSEq <- metaCompareKW(SEq, meta, 'Shannon equitability')
(kwSEq <- kwSEq[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )

## FISHER'S ALPHA
Fishers <- data.frame(Groups=meta$grpNames, ProSM=fisher.alpha(dataToUseT), pCh=meta$pCh, pCol= meta$pCol)
kwFis <- metaCompareKW(Fishers, meta, "Fisher's alpha")
(kwFis <- kwFis[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )

## SIMPSON'S D
Simpson <- diversity(dataToUseT, index='simpson')
kwSim <- metaCompareKW(Simpson, meta, "Simpson's D")
(kwSim <- kwSim[,c('metaName','comp','p.adj','p.txt','p.value','statistic','parameter')] )



pdf(file=paste0(PATH.FIG,'Fig03-', TAXONLEVEL, DROPRARETAXA,'-3X-diversity.pdf'), width=pageWidthOne, height=pageHeight)
par(mfcol=c(5,1), oma=c(8,3,1,1), mar=c(1,1,1,1), las=1)

panelLabelSide <- 1
panelLabelLine <- -2
plotBarWhiskerPVB(SObs,"Observed richness", panelLabelSize, panelLabelSide, panelLabelLine, kwObs, pValSide=3, pValLines, pValAdj, pCh, pCol, addNtext=TRUE)
mtext('A',side=3, line=-0.5, cex=1.2, adj=0.05)
mtext(TAXONLEVELPlural, side=2, line=2.2, las=0, cex=0.8)

nv <- aggregate(SObs$Groups,by=list(SObs$Groups),FUN=length)
mv <- aggregate(SObs$ProSM,by=list(SObs$Groups),FUN=max)
sq <- c(5,1,3,4,2,6)
text(x=1:6,y=mv$x[sq],paste0('n=',nv$x[sq]), pos=3)

plotBarWhiskerPVB(Chao,"Chao1", panelLabelSize, panelLabelSide, panelLabelLine, kwCha, pValSide=3, pValLines, pValAdj, pCh, pCol)
mtext('B',side=3, line=-0.5, cex=1.2, adj=0.05)
mtext(TAXONLEVELPlural, side=2, line=2.2, las=0, cex=0.8)

plotBarWhiskerPVB(SEq,"Shannon equitability", panelLabelSize, panelLabelSide, panelLabelLine, kwSEq, pValSide=3, pValLines, pValAdj, pCh, pCol, yRange=c(0.5,1))
mtext('C',side=3, line=-0.5, cex=1.2, adj=0.05)

plotBarWhiskerPVB(Fishers,"Fisher's alpha", panelLabelSize, panelLabelSide, panelLabelLine, kwFis, pValSide=3, pValLines, pValAdj, pCh, pCol)
mtext('D',side=3, line=-0.5, cex=1.2, adj=0.05)

axis(1, at=1:6, labels=GroupNames, las=2)
box('outer')

mtext('OTR Core', side=1, line=6.5, adj=0.14, cex=0.8)
mtext('Surface', side=1, line=5.3, adj=0.82, cex=0.8)
mtext('Grabs', side=1, line=6.5, adj=0.79, cex=0.8)

dev.off()









## Foram Stress Index
fsiMeta[duplicated(fsiMeta$Genus),]

ProportionSI <- dataToUseT[,-1]/rowSums(dataToUseT[,-1])

fsi.Sensitive.taxa <- fsiMeta[fsiMeta$FSI == 'Sensitive', TAXONLEVEL]
fsi.Tolerant.taxa <- fsiMeta[fsiMeta$FSI == 'Tolerant', TAXONLEVEL]

fsi.Sensitive.Cols <- which(colnames(ProportionSI) %in% fsi.Sensitive.taxa)
colnames(ProportionSI)[fsi.Sensitive.Cols]

fsi.Tolerant.Cols <- which(colnames(ProportionSI) %in% fsi.Tolerant.taxa)
colnames(ProportionSI)[fsi.Tolerant.Cols]

FSI <- rowSums(ProportionSI[, fsi.Sensitive.Cols])*10 + rowSums(ProportionSI[, fsi.Tolerant.Cols])

## FORAM index
foram.Hetero.taxa <- fsiMeta[fsiMeta$Functional_Type == 'Heterotrophic', TAXONLEVEL]
foram.Symbio.taxa <- fsiMeta[fsiMeta$Functional_Type == 'Symbiont', TAXONLEVEL]
foram.Opport.taxa <- fsiMeta[fsiMeta$Functional_Type == 'Opportunistic', TAXONLEVEL]

foram.Hetero.Cols <- which(colnames(ProportionSI) %in% foram.Hetero.taxa)
foram.Symbio.Cols <- which(colnames(ProportionSI) %in% foram.Symbio.taxa)
foram.Opport.Cols <- which(colnames(ProportionSI) %in% foram.Opport.taxa)

FORAM <- rowSums(ProportionSI[, foram.Symbio.Cols])*10 + rowSums(ProportionSI[, foram.Opport.Cols]) + rowSums(ProportionSI[, foram.Hetero.Cols])*2


##	SUMMMARY TABLE
#Sum.Table <- rbind(Chao.Table, Shannon, SEq, Fishers, Simpson, FSI, FORAM)
Sum.Table <- rbind(Chao.Table, SEq=SEq$ProSM, Fishers=Fishers$ProSM)
write.csv(Sum.Table, file=paste0(PATH.OUT, TAXONLEVEL, DROPRARETAXA,'-Table_DiversityStats-sample.csv'))

##	COMPARE DIVERSITY METRICS
##########################################################################################

##	STATS FOR SUMMARY TABLE
sAGrp <- statsMe(Sum.Table, meta$grpNames)
sColl <- statsMe(Sum.Table, meta$compareCollect)
sTime <- statsMe(Sum.Table, meta$compareTime)
sReef <- statsMe(Sum.Table, meta$compareReef)

statSum <- cbind(sAGrp,sColl, sTime, sReef)

for (r in 1:nrow(statSum))
	statSum <- adjust.p.value(statSum,r)
	
statSumText <- as.data.frame(statSum[, which(colnames(statSum) %in% "p.adj")])
rownames(statSumText) <- statSum$var
colnames(statSumText) <- c('allGroups','sCollect','sTime','sReef')
statSumText <- round(statSumText,3)
statSumText$allGroups <- paste('p =', statSumText$allGroups)
statSumText$sCollect <- paste('p =', statSumText$sCollect)
statSumText$sTime <- paste('p =', statSumText$sTime)
statSumText$sReef <- paste('p =', statSumText$sReef)
statSumText[(statSumText == 'p = 0')] <- 'p < 0.001'

write.csv(statSumText, file=paste0(PATH.OUT, 'Table1-Fig2-', TAXONLEVEL, DROPRARETAXA,'-Diversity-kw-summary.csv'))

##	AGGREGATE STATS BY REEF / TIME PERIOD
SumT <- t(Sum.Table[])
Grps <- substring(colnames(Chao.Table),0,1)
Grps[1:3] <- 'M'
Grps[4:6] <- 'C'
Grps[7:9] <- 'P'
sumM <- aggregate(SumT,by=list(meta$grpNames),mean)
sumE <- aggregate(SumT,by=list(meta$grpNames),sd)

##	put in preferred order
sumM <- sumM[c(5,1,3,4,2,6),]
sumE <- sumE[c(5,1,3,4,2,6),]

write.csv(rbind(sumM, sumE), file=paste0(PATH.OUT, TAXONLEVEL, DROPRARETAXA,'-Table_DiversityStats-summary.csv'))


##	INDICATOR SPECIES ANALYSIS
##########################################################################################

d2uType <- dataToUseT[!is.na(meta$compareCollect),]
d2uTime <- dataToUseT[!is.na(meta$compareTime),]
d2uReef <- dataToUseT[!is.na(meta$compareReef),]

d2uTypeCat <- meta$type[!is.na(meta$compareCollect)]
d2uTimeCat <- meta$age[!is.na(meta$compareTime)]
d2uReefCat <- meta$reef[!is.na(meta$compareReef)]

nperm <- 99999
#nperm <- 9999

for (i in 1:5) {

isaGrps <- multipatt(dataToUseT, meta$grpNames, control = how(nperm= nperm))
isaType <- multipatt(d2uType, d2uTypeCat, control = how(nperm= nperm))
isaReef <- multipatt(d2uReef, d2uReefCat, control = how(nperm= nperm))
isaTime <- multipatt(d2uTime, d2uTimeCat, control = how(nperm= nperm))

#d2uPA <- ifelse(dataToUseT>0,1,0)
#phi <- multipatt(d2uPA, meta$grpNames, func = "r.g", control = how(nperm=999))
#summary(phi)
#round(head(phi$str),3)

#summary(isaGrps, indvalcomp=TRUE)
#summary(isaType, indvalcomp=TRUE)
#summary(isaTime, indvalcomp=TRUE)
#summary(isaReef, indvalcomp=TRUE)

out <- capture.output(summary(isaGrps, indvalcomp=TRUE, alpha=0.05))
write(out, file=paste0(PATH.OUT, 'Table3', TAXONLEVEL, DROPRARETAXA,'-isaGrps','-',i,'.txt'))

out <- capture.output(summary(isaType, indvalcomp=TRUE, alpha=0.05))
write(out, file=paste0(PATH.OUT, 'Table4', TAXONLEVEL, DROPRARETAXA,'-isaType','-',i,'.txt'))

out <- capture.output(summary(isaTime, indvalcomp=TRUE, alpha=0.05))
write(out, file=paste0(PATH.OUT, 'Table5', TAXONLEVEL, DROPRARETAXA,'-isaTime','-',i,'.txt'))

out <- capture.output(summary(isaReef, indvalcomp=TRUE, alpha=0.05))
write(out, file=paste0(PATH.OUT, 'Table6', TAXONLEVEL, DROPRARETAXA,'-isaReef','-',i,'.txt'))

}


##	PLOTS
##########################################################################################


##	NONMETRIC MULTIDIMENSIONAL SCALING
##########################################################################################
nmds.dataToUseT <- metaMDS(dataToUseT)
stressplot(nmds.dataToUseT)
nmds.dataToUseT

##	pull siteScores so enable more nuanced plotting...
siteScores <- as.data.frame(nmds.dataToUseT$points[,c('MDS1','MDS2')])

##	make columns for point style (point type & colour)
siteScores$pCol <- NA
siteScores$pCh <- NA
siteScores$grp <- NA

##	set point style and colour by reef / time period
siteScores[rownames(siteScores,0,1) %in% c('L110','L150','L140'),'pCol'] <- pCol[1]
siteScores[rownames(siteScores,0,1) %in% c('L060','L032','L037'),'pCol'] <- pCol[2]
siteScores[rownames(siteScores,0,1) %in% c('L002','L007','L012'),'pCol'] <- pCol[3]
siteScores[substring(rownames(siteScores),0,1)=='O','pCol'] <- pCol[4]
siteScores[substring(rownames(siteScores),0,1)=='H','pCol'] <- pCol[5]
siteScores[substring(rownames(siteScores),0,1)=='W','pCol'] <- pCol[6]

siteScores[rownames(siteScores,0,1) %in% c('L002','L007','L012'),'pCh'] <- pCh[1]
siteScores[rownames(siteScores,0,1) %in% c('L032','L037','L060'),'pCh'] <- pCh[2]
siteScores[rownames(siteScores,0,1) %in% c('L110','L140','L150'),'pCh'] <- pCh[3]
siteScores[substring(rownames(siteScores),0,1)=='O','pCh'] <- pCh[4]
siteScores[substring(rownames(siteScores),0,1)=='H','pCh'] <- pCh[5]
siteScores[substring(rownames(siteScores),0,1)=='W','pCh'] <- pCh[6]

GroupNames2 <- c('Core','OneTree','Heron','Wistari')

siteScores[rownames(siteScores,0,1) %in% c('L110','L150','L140'),'grp'] <- GroupNames2[1]
siteScores[rownames(siteScores,0,1) %in% c('L060','L032','L037'),'grp'] <- GroupNames2[1]
siteScores[rownames(siteScores,0,1) %in% c('L002','L007','L012'),'grp'] <- GroupNames2[1]
siteScores[substring(rownames(siteScores),0,1)=='O','grp'] <- GroupNames[4]
siteScores[substring(rownames(siteScores),0,1)=='H','grp'] <- GroupNames[5]
siteScores[substring(rownames(siteScores),0,1)=='W','grp'] <- GroupNames[6]


##	NMDS PLOT
####################
pdf(file=paste0(PATH.FIG,'Fig05-', TAXONLEVEL, DROPRARETAXA,'-nmdsPlot.pdf'),width = pageWidthTwo, height= pageWidthTwo)
par(mar=c(3,3,1,1), las=0, mgp=c(2,0.75,0), cex=1.0, las=1)

xRange <- range(siteScores$MDS1)
yRange <- range(siteScores$MDS2)

ordiplot(nmds.dataToUseT, type="n", xlim=xRange, ylim=yRange)

ordiellipse(nmds.dataToUseT, siteScores$grp, col= names(table(siteScores$pCol)), alpha=0.2, lwd=1, draw='polygon', conf=0.95)

points(siteScores$MDS1, siteScores$MDS2, pch=as.integer(siteScores$pCh), bg=siteScores$pCol, col='black', cex=2)


#text(siteScores[,'MDS1'], siteScores[,'MDS2'], labels= row.names(siteScores), pos=3)

#orditorp(nmds.dataToUseT, display="species", col="black", air=0.8, font=3, cex=0.56)

#mtext(TAXONLEVEL,side=3, line=-1.25, adj=0.95)
#mtext(paste0("> ", DROPRARETAXA,"%"),side=3, line=-2, adj=0.95)

#legend('bottomleft',legend= GroupNames, pch=pCh, col=pCol, pt.bg=pCol)
legend('bottomright',legend= GroupNames[c(1,4,2,5,3,6)], pch=pCh[c(1,4,2,5,3,6)], col=pCol[c(1,4,2,5,3,6)], pt.bg=pCol[c(1,4,2,5,3,6)], ncol =3, bty='n', cex=0.9)

mtext('OTR Core:', adj=0.02, side=1, line=-2.35, cex=0.9)
mtext('Surface Grabs:',adj=0.02, side=1, line=-1.5, cex=0.9)


box('outer')
dev.off()


##	NMDS AXIS SCORES
####################
write.csv(siteScores, file=paste0(PATH.OUT, TAXONLEVEL, DROPRARETAXA,'-Table_nMDS_siteScores.csv'))
write.csv(nmds.dataToUseT$species, file=paste0(PATH.OUT, TAXONLEVEL, DROPRARETAXA,'-Table_nMDS_taxonScores.csv'))


nmds.dataToUseT$species[order(nmds.dataToUseT$species[,'MDS1']),]
plot(nmds.dataToUseT$species[order(nmds.dataToUseT$species[,'MDS1']),'MDS1'])

nmds.dataToUseT$points[order(nmds.dataToUseT$points[,'MDS1']),]
plot(nmds.dataToUseT$points[order(nmds.dataToUseT$points[,'MDS1']),'MDS1'])


source('./textStats.R')
source('./plot-Sed&Lead.R')


