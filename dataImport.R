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
## This script is designed to be called from analyses.R
##
## This script will import the original data files:
##	Mamo_Lagoon_Forams.csv containing published foram count data from Mamo 2016.
##	OTR_Core_Forams.csv containing the foram count data from Bauder 2019 (MRes Thesis).
##	OTR_Core_pb210_CIC_Ages.csv containing published pb-210 age data from Kosnik 2015.
##
## This script will subset the data and remove the data not used in these analyses.
##
## This script enables the options to:
##	1) Drop 63 fraction data, 
##	2) Combine fractions into single samples, 
##	3) Drop unidentfied taxa 
##	4) Drop juveniles.
##	5) Recreate the dataset at arbritary taxonomic levels (only tested Genus and Family)
##	6) Drop taxa with fewer specimens than a specified cutoff (as proportion of dataset realtive to 1%) 

############################  Packages needed as well as base R  #################
{my_packages <- c("stringr")
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed)}

library(stringr)

## USE THE SAME PLOT PAR ACROSS ALL PLOTS
#source('plotPreferences.R')

## READ MAMO SURFACE SAMPLES - FORAM DATA
dataMamo <- read.csv("./data/Mamo_Lagoon_Forams.csv", as.is=TRUE, blank.lines.skip=TRUE)
tail(dataMamo)

## READ OTR CORE - AGE DATA
dataAge <- read.csv("./data/OTR_Core_pb210_CIC_Ages.csv", row.names=1)
head(dataAge)

## READ OTR CORE - FORAM DATA
dataCore <- read.csv("./data/OTR_Core_Forams.csv", as.is=TRUE, blank.lines.skip=TRUE)
tail(dataCore)

### DATA META
###############################
#fsiMeta <- rbind(dataMamo[,c(TAXONLEVEL,'Functional_Type','Morphotype')], dataCore[,c(TAXONLEVEL,'Functional_Type','Morphotype')])
fsiMeta <- rbind(dataMamo[,c('Family','Genus','Functional_Type','Morphotype')], dataCore[,c('Family','Genus','Functional_Type','Morphotype')])
for (colNam in colnames(fsiMeta)) {
	fsiMeta[, colNam] <- trimws(fsiMeta[, colNam], which='both', whitespace = "[\\h\\v]")	
}
fsiMeta <- fsiMeta[!duplicated(fsiMeta),]
fsiMeta <- fsiMeta[!is.na(fsiMeta[,TAXONLEVEL]),]

summary(as.factor(fsiMeta$Morphotype))

fsiMeta$FSI <- NA
fsiMeta[(fsiMeta$Functional_Type == 'Heterotrophic'),'FSI']	<- 'Sensitive'
fsiMeta[(fsiMeta$Functional_Type == 'Symbiont'),'FSI']		<- 'Sensitive'
fsiMeta[(fsiMeta$Functional_Type == 'Opportunistic'),'FSI'] <- 'Tolerant'

summary(as.factor(fsiMeta$Functional_Type))

fsiMeta$infaunal <- NA
fsiMeta[(fsiMeta$Morphotype %in% c('Rounded.planispiral','Flattened.ovoid','Tapered.cylindrical','Flattened.tapered')),'infaunal'] <- 'TRUE'
fsiMeta[(fsiMeta$Morphotype %in% c('Biconvex.trochospiral','Rounded.trochospiral','Milioline','Plano.convex')),'infaunal'] <- 'FALSE'

summary(as.factor(fsiMeta$infaunal))

fsiMeta$LBF <- FALSE
fsiMeta[(fsiMeta$Family %in% c('Alveolinidae', 'Amphisteginidae', 'Calcarinidae', 'Nummulitidae', 'Soritidae', 'Peneroplidae')),'LBF'] <- TRUE

fsiMeta$smallMiliolid <- FALSE
fsiMeta[(fsiMeta$Family %in% c('Spiroloculinidae', 'Nubeculariidae', 'Riveroinidae', 'Ophthalmidiidae', 'Cornuspiridae', 'Hauerinidae', 'Cribrolinoididae')),'smallMiliolid'] <- TRUE

fsiMeta$Bolivinids <- FALSE
fsiMeta[(fsiMeta$Family %in% c('Tortoplectellidae', 'Bolivinitidae', 'Bolivinellidae', 'Buliminoididae','Uvigerinidae','Reussellidae')),'Bolivinids'] <- TRUE

fsiMeta$Agglutinated <- FALSE
fsiMeta[(fsiMeta$Family %in% c('Textulariidae', 'Pseudogaudryinidae', 'Eggerellidae', 'Olgiidae','Reophacidae')),'Agglutinated'] <- TRUE

write.csv(fsiMeta, file=paste0('./output/Supplement2-FSIcategories-',TAXONLEVEL,'.csv'))

### DATA CLEANING & SUBSETTING -- OTR CORE
###############################

## DROP ALL THE COLUMNS WITH 63 MICRON DATA
if (DROP63 == TRUE) {
  subSet <- (substr(colnames(dataCore),8,9) != 63)
  data63 <- dataCore[,subSet]
}

#	CLEAN UP THE RELEVANT TAXON / GROUPING COLUMN
dataCore <- dataCore[!is.na(dataCore[, TAXONLEVEL]),]
dataCore <- dataCore[(dataCore[,TAXONLEVEL] != ''),]
dataCore[,TAXONLEVEL] <- trimws(dataCore[,TAXONLEVEL], which='both', whitespace = "[\\h\\v]")

# SET NA <- 0
dataCore[is.na(dataCore)] <- 0
dataCore <- dataCore[,which(substr(colnames(dataCore),0,1)!='X')]

colNames <- names(dataCore)
layNames <- colNames[(substr(colnames(dataCore),0,1) == 'L')]

## KEEP ONLY COLUMNS NEEDED FOR GROUPING & DATA
dataSub <- dataCore[,c(TAXONLEVEL,layNames)]
head(dataSub)

## DROP ALL THE COLUMNS WITH 63 MICRON DATA
if (DROP63 == TRUE) {
  subSet <- (substr(colnames(dataSub),8,9) != 63)
  dataSub <- dataSub[,subSet]
}

## COMBINE FRACTIONS INTO A SINGLE SAMPLE 
if (COMBINEFRACTIONS == TRUE) {
  head(dataSub)

  # UNIQUE layer names (MINUS THE GROUPING COLUMN)
  dslayers <- substr(colnames(dataSub),0,4)
  usdlayers <- unique(dslayers)[-1]
  # NEW DATA FRAME
  datalayer <- data.frame(Genus=dataSub[,1], x=matrix(nrow=nrow(dataSub), ncol=length(usdlayers)))
  colnames(datalayer) <- c(TAXONLEVEL,usdlayers)
  usdlayers
  for(u in usdlayers){
    temp <- dataSub[,(dslayers==u)]
    datalayer[,u] <- rowSums(temp)
  }
  head(dataSub)
  dataSub <- datalayer
}


##	drop rows with unidentified specimens
if (DROPUNIDENTIFIED == TRUE) {
	groups <- dataSub[,1]
	dropIfFound <- "Unidentified"
	groupsTrunc <- substring(groups,0, nchar(dropIfFound))
	selectMe <- (groupsTrunc != dropIfFound)
	dataSub <- dataSub[selectMe,]
}

##	drop rows with juvenile specimens
if (DROPJUVENILE == TRUE) {
  groups <- dataSub[,1]
  dropIfFound <- "Juvenile"
  groupsTrunc <- substring(groups,0, nchar(dropIfFound))
  selectMe <- (groupsTrunc != dropIfFound)
  dataSub <- dataSub[selectMe,]
}

##	DROP TAXA WITH NO SPECIMENS...
dataSub <- dataSub[(rowSums(dataSub[,2:ncol(dataSub)])>0),]

## Group by first / GROUPING column
dataGrp <- aggregate(dataSub[,2:ncol(dataSub)], by=list(dataSub[,TAXONLEVEL]), FUN=sum)

## number of forams per sample and per taxon
sampleTotals <- colSums(dataGrp[,2:ncol(dataGrp)])
groupTotals <- rowSums(dataGrp[,2:ncol(dataGrp)])

taxonSummary <- data.frame(Group.1=dataGrp[,1], groupTotals)
taxonSummary[order(taxonSummary[,'groupTotals'], decreasing=TRUE),]

###Calculate 1%
sampleTotals <- colSums(dataGrp[,2:ncol(dataGrp)])
onePercent <- sum(sampleTotals)/100
onePercent

# DROP TAXA WITH TOO FEW SPECIMENS
if (DROPRARETAXA != 0) {
	(THRESHOLD <- floor(onePercent * DROPRARETAXA) )
	dataGrp <- dataGrp[(rowSums(dataGrp[,2:ncol(dataGrp)]) > THRESHOLD),]
}

### DATA ORDERING

## MAKE A COLUMN IN THE AGE DATA THAT WILL MATCH COLUMN NAMES IN FAUNA DATA
ageMeta <- data.frame(dataAge)
ageMeta$name <- paste0('L', format(ageMeta$top,digits=3))
ageMeta$name <- str_replace(ageMeta$name,' ','0')
ageMeta$name <- str_replace(ageMeta$name,' ','0')
ageMeta$name <- substring(ageMeta$name,0,4)

## order oldest to youngest, bigest to smallest BY MERGING WITH AGE DATA
cn <- colnames(dataGrp)
sMeta <- data.frame(cols=cn[2:length(cn)], seq=1:(length(cn)-1))
sMeta$name <- substring(sMeta$cols,0,4)
sMeta$frac <- substring(sMeta$cols,7,9)
sMeta <- merge(sMeta, ageMeta)
cn3 <- as.character(sMeta$cols)
dataGrp <- dataGrp[,c('Group.1', cn3)]
#dataGrp<-dataGrp[-1,]

## TRANSPOSE AND ORDER DATA FOR VEGAN
dataGrpT <- t(dataGrp[,-1])
colnames(dataGrpT) <- dataGrp[,1]

dataCore <- dataGrp
dataCoreT <- dataGrpT


### DATA CLEANING & SUBSETTING -- MAMO SURFACE SAMPLES
###############################
head(dataMamo)

dataMamo <- dataMamo[!is.na(dataMamo[, TAXONLEVEL]),]
dataMamo <- dataMamo[(dataMamo[,TAXONLEVEL] != ''),]
dataMamo[,TAXONLEVEL] <- trimws(dataMamo[,TAXONLEVEL], which='both', whitespace = "[\\h\\v]")

# SET NA <- 0
dataMamo[is.na(dataMamo)] <- 0
dataMamo <- dataMamo[,which(substr(colnames(dataMamo),0,1)!='X')]

## ONLY KEEP NEEDED COLUMNS
tCol <- which(names(dataMamo) == TAXONLEVEL)
dataMamo <- dataMamo[,c(tCol,7:ncol(dataMamo))]
head(dataMamo)

## Group by first / GROUPING column
dataMamo <- aggregate(dataMamo[,2:ncol(dataMamo)], by=list(dataMamo[,TAXONLEVEL]), FUN=sum)

taxonSummary <- data.frame(Group.1= dataMamo[,1], groupTotals= rowSums(dataMamo[,2:ncol(dataMamo)]))
taxonSummary[order(taxonSummary[,'groupTotals'], decreasing=TRUE),]

###Calculate 1%
sampleTotals <- colSums(dataMamo[,2:ncol(dataMamo)])
onePercent <- sum(sampleTotals)/100
onePercent

# DROP TAXA WITH TOO FEW SPECIMENS
if (DROPRARETAXA != 0) {
	
	(THRESHOLD <- floor(onePercent * DROPRARETAXA) )
	dataMamo <- dataMamo[(rowSums(dataMamo[,2:ncol(dataMamo)]) > THRESHOLD),]
}

# abbreviate all reefs with one letter
newNames <- str_replace(names(dataMamo),'OTI','O')
names(dataMamo) <- newNames

## TRANSPOSE AND ORDER DATA FOR VEGAN
dataMamoT <- t(dataMamo[,-1])
colnames(dataMamoT) <- dataMamo[,1]
#dataGrpT <- dataGrpT[,-1]


## COMBINE
###############################
dataToUse <- merge(dataCore, dataMamo, all=TRUE)
dataToUse[is.na(dataToUse)] <- 0

dataToUseT <- t(dataToUse[,-1])
colnames(dataToUseT) <- dataToUse[,1]

#metaDM <- dataMamo
#metaDC <- dataMamo

## META 
###############################
meta <- data.frame(sample= rownames(dataToUseT),reef=NA,type=NA,age=NA, grpNames=NA)
meta$seq <- 1:nrow(meta)

meta[(substring(meta$sample,0,1) %in% c('L','O')),'reef'] <- 'OneTree'
meta[(substring(meta$sample,0,1) == 'W'),'reef'] <- 'Wistari'
meta[(substring(meta$sample,0,1) == 'H'),'reef'] <- 'Heron'

meta[(substring(meta$sample,0,1) == 'L'),'type'] <- 'Core'
meta[(substring(meta$sample,0,1) != 'L'),'type'] <- 'Grab'

meta[(meta$reef == 'OneTree') & (meta$type == 'Grab'),'age'] <- 'Modern'
meta[(meta$reef == 'OneTree') & (substring(meta$sample,0,1) == 'L'),'age'] <- 'Colonial'
meta[(meta$reef == 'OneTree') & as.numeric(substring(meta$sample,2,4)) < 15,'age'] <- 'Modern'
meta[(meta$reef == 'OneTree') & as.numeric(substring(meta$sample,2,4)) > 100,'age'] <- 'pre-Colonial'

meta$grpNames <- meta$reef
meta[1:9,'grpNames'] <- meta[1:9,'age']

meta$compareCollect <- NA
meta[(substring(meta$sample,0,1) == 'L') & (meta$age == 'Modern'),'compareCollect'] <- 'Core'
meta[(substring(meta$sample,0,1) == 'O'),'compareCollect'] <- 'Grab'

meta$compareTime <- meta$grpNames
meta[(meta$grpNames == 'Heron'),'compareTime'] <- NA
meta[(meta$grpNames == 'Wistari'),'compareTime'] <- NA
meta[(meta$grpNames == 'OneTree'),'compareTime'] <- NA

meta$compareReef <- meta$reef
meta[(meta$type == 'Core'),'compareReef'] <- NA

meta <- merge(meta,pMeta, by.x='grpNames',by.y='GroupNames')
meta <- meta[order(meta$seq),]
