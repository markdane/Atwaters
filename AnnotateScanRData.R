#12/2015 preprocess scanr data
library("MEMA")
library("limma")
library("data.table")
library("parallel")
library(XLConnect)
library(ggplot2)

medianDT <- function(x) x/median(x,na.rm=TRUE)

#TODO Move these functions into the MEMA package


kmeansClusterEdUPositiveWellPlate <- function(x,value){
  EdUPos <- kmeansCluster(log2(x[[value]]))
  edUPositiveThresh <- min(x[[value]][EdUPos==2])
  clusters <- rep.int(0,nrow(x))
  clusters[x[[value]]>edUPositiveThresh] <- 1
  return(clusters)
}


fixedThreshEdUPositive <- function(x,value,thresh){
  clusters <- rep.int(1,nrow(x))
  clusters[x[[value]]>thresh] <- 2
  return(clusters)
}


readWellMetadata <- function (xlsFile) 
{
  sheetList <- sapply(gdata::sheetNames(path.expand(xlsFile)), 
                      gdata::read.xls, xls = path.expand(xlsFile), simplify = FALSE, 
                      stringsAsFactors = TRUE, check.names = FALSE, row.names = "Row/Column")
  nrRows <- dim(sheetList[[1]])[1]
  nrCols <-  dim(sheetList[[1]])[2]
  nrWells = nrRows * nrCols
  sheetDF <- data.frame(lapply(sheetList, function(df, nrCols) {
    dfM <- matrix(t(df[, 1:nrCols]), byrow = TRUE)
  }, nrCols = nrCols), WellIndex = 1:nrWells, Well = wellAN(nrRows, 
                                                            nrCols), check.names = TRUE, stringsAsFactors = FALSE)
  return(sheetDF)
}

# Function from Roland on Stack Overflow
#http://stackoverflow.com/questions/23018056/convert-cartesian-angles-to-polar-compass-cardinal-angles-in-r
calcTheta <- function(x,y) {
  z <- x + 1i * y
  res <- 90 - Arg(z) / pi * 180
  res %% 360
}

#\code{findPerimeterCell} Determine the perimeter cell in wedge
#
# @param x A datatable or dataframe with a RadialPosition column
# @return A logical vector the length of x with a TRUE value for the Perimeter cell
#
#
findPerimeterCell <- function(x){
  if(!nrow(x)==0){
    perimeterLogicals <- vector(length=nrow(x))
    perimeterLogicals[which.max(x$RadialPosition)] <- TRUE
  }
  return(perimeterLogicals)
}

#Classify Outer cells
#
#\code{labelOuterCells} Determine if a cell is an Outer cell in the spot
#
# @param x A numeric vector of RadialPositions
# @param thresh A quantile value between 0 and 1 used to threshold x. The returned logical will be TRUE for cells with x values in quantiles greater than thresh.
# @return A logical vector the length of x with a TRUE value for the Outer cells
#
#
labelOuterCells <- function(x, thresh=.75){
  outerLogicals <- NULL
  if(!length(x)==0){
    outerLogicals <- x>quantile(x,probs = thresh, na.rm=TRUE)
  }
  return(outerLogicals)
}


wellPositionParms <- function (DT, densityRadius = 160, outerThresh = 0.2, wedges = 18, sparseThresh = 0.8) 
{
  lDT <- copy(DT)
  lDT <- lDT[, `:=`(XLocal, X - median(X)), by = "Barcode,Well"]
  lDT <- lDT[, `:=`(YLocal, Y - median(Y)), by = "Barcode,Well"]
  lDT <- lDT[, `:=`(RadialPosition, sqrt(XLocal^2 + YLocal^2))]
  lDT <- lDT[, `:=`(Theta, calcTheta(XLocal, YLocal))]
  lDT <- lDT[, `:=`(Density, spotCellDensities(.SD, radius = densityRadius) * 
                      10000), by = "Barcode,Well"]
  lDT <- lDT[, `:=`(Sparse, Density < sparseThresh)]
  wedgeAngs <- 360/wedges
  lDT <- lDT[, `:=`(Wedge, ceiling(Theta/wedgeAngs))]
  lDT <- lDT[, `:=`(OuterCell, labelOuterCells(RadialPosition, 
                                               thresh = outerThresh)), by = "Barcode,Well"]
  denseOuterDT <- lDT[!lDT$Sparse & lDT$OuterCell]
  denseOuterDT <- denseOuterDT[, `:=`(Perimeter, findPerimeterCell(.SD)), 
                               by = "Barcode,Well,Wedge"]
  setkey(lDT, Barcode, Well, ObjectID)
  setkey(denseOuterDT, Barcode, Well, ObjectID)
  lDT <- denseOuterDT[, list(Barcode, Well, ObjectID, 
                             Perimeter)][lDT]
  lDT$Perimeter[is.na(lDT$Perimeter)] <- FALSE
  return(lDT[, list(Barcode, Well, ObjectID, XLocal, 
                    YLocal, RadialPosition, Theta, Wedge, Density, Sparse, 
                    OuterCell, Perimeter)])
}

#Calculate Robust Z Scores of a numeric vector
robustZScores <- function(x){
  xMedian <- median(x,na.rm=TRUE)
  xMad <- mad(x,na.rm=TRUE)
  robustZScores <- (x-xMedian)/xMad
}

#Get all detail on the data files to be stitched together
getWellSubsetFileNames <- function(path){
  rdFiles <- data.frame(FilePaths = dir(path, full.names = TRUE), stringsAsFactors = FALSE)
  splits <- strsplit2(rdFiles$FilePaths,"_")
  rdFiles$Barcode <- gsub(".*/","",splits[,1])
  #Todo: read these values from the google doc instead of the file names
  rdFiles$Wells <- gsub("W","",splits[,4])
  rdFiles$Location <- gsub(".txt","",splits[,5])
  rdFiles$CellLine <- splits[,2]
  rdFiles$PlateType <- as.integer(strsplit2(rdFiles$Barcode,"")[,7])
  return(rdFiles)
}

stitchWellData <- function(fs){
  dtL <-apply(fs[grepl("_Main",fs$FilePaths),], 1, function(fa){
    mdt <- fread(fa["FilePaths"])
    cdt <- fread(gsub("Main","Cyto",fa["FilePaths"]))
    mdt$Barcode <- fa["Barcode"]
    mdt$CellLine <- fa["CellLine"]
    mdt <- convertColumnNames(mdt)
    cdt <- convertColumnNames(cdt)
    setnames(cdt,grep("Intensity",colnames(cdt),value=TRUE), paste0(grep("Intensity",colnames(cdt),value=TRUE),"Cyto"))
    setkey(mdt,ObjectID)
    setkey(cdt,ParentObjectIDMO)
    mcdt <- mdt[cdt]
    mcdt <- mcdt[,grep("^i[.]|ObjectID|ParentTraceID",colnames(mcdt), value=TRUE,invert=TRUE), with=FALSE]
    return(mcdt)
  })
  dt <- rbindlist(dtL)
  #TODO: Add check for redundant well data
}

#########

#Get the raw data file names
rdf <- getWellSubsetFileNames("RawData")

#Combine raw data from each plate
cDT <- rbindlist(lapply(unique(rdf$Barcode), function(barcode, df){
  bdf <- df[df$Barcode==barcode,]
  bDT <- stitchWellData(fs=bdf)
}, df=rdf))

densityThresh <- 0.4
outerThresh <- 0.5
ss <- "lineageEdU"

#Count the cells at each well
cDT<-cDT[,WellCellCount := .N, by="Barcode,Well"]

#Filter out debris based on nuclear area
nuclearAreaThresh <- 50
cDT <- cDT[cDT$Area >nuclearAreaThresh,]

cDT$TotalIntensityDAPI <- cDT$Area*cDT$MeanIntensityDAPI

#Add row and column indices
cDT$Row <- ceiling(cDT$Well/24)
cDT$Column <- (cDT$Well-1) %% 24 + 1

# #Read in anotations for the gene class
# geneAnn <- data.table(read.xls("2015June_pilot-sirna-screen-v4-1.xlsx"),key="Gene")
# setkey(cDT,GeneSymbol)
# cDT <- cDT[geneAnn]

normToGeneSymbol <- function(DT, value, baseGeneSymbol) {
  if(!c("GeneSymbol") %in% colnames(DT)) stop(paste("DT must contain a GeneSymbol column."))
  valueMedian <- median(unlist(DT[grepl(baseGeneSymbol, DT$GeneSymbol),value, with=FALSE]), na.rm = TRUE)
  normedValues <- DT[,value,with=FALSE]/valueMedian
  return(normedValues)
}

# 
#   cDT <- cDT[,MeanIntensityAlexa488WellNorm := normToGeneSymbol(.SD,value = "MeanIntensityAlexa488Cyto",baseGeneSymbol = "NegCtrl"), by="CellLine"]
#   cDT <- cDT[,CytoMeanIntensityAlexa555WellNorm := normToGeneSymbol(.SD,value = "MeanIntensityAlexa555Cyto",baseGeneSymbol = "NegCtrl"), by="CellLine"]
#   cDT <- cDT[,MeanIntensityAlexa647WellNorm := normToGeneSymbol(.SD,value = "MeanIntensityAlexa647",baseGeneSymbol = "NegCtrl"), by="CellLine"]
#   cDT$LineageRatio <- cDT$CytoMeanIntensityAlexa488/cDT$CytoMeanIntensityAlexa555
#   #Classify each cell based on the plate EdU levels
#   #TODO: add control gene ID to crate gate value
#   cDT <- cDT[,EdUPositive := kmeansClusterEdUPositiveWellPlate(.SD, value="MeanIntensityAlexa647"), by="Barcode"]
#   #Calculate the EdU Positive Proportion at each well
#   cDT <- cDT[,EdUPositiveProportion := sum(EdUPositive)/length(EdUPositive),by="Barcode,Well"]

#Summarize cell data to well level by taking the medians of these parameters
wNames<-grep(pattern="(Total|Mean|Area|EdUPositiveProportion|Cyto|Density|Z|X53BP1|WellCellCount|Barcode|^Spot$|^Well$)",x=names(cDT),value=TRUE)
#Remove the well normalized values as the median is 1 by definition
wNames <- grep("WellNorm",x=wNames, value=TRUE, invert=TRUE)
wKeep<-cDT[,wNames,with=FALSE]
wDT<-wKeep[,lapply(.SD,numericMedian),keyby="Barcode,Well"]

#Merge back in the well and plate metadata
mDT <- cDT[,grep("Barcode|^Well$|CellLine|WellIndex|Row|Column|GeneSymbol|EndpointDAPI|Endpoint488|Endpoint555|Endpoint647|Class|Drugs..FDA.Approved.|Drugs..Kinase.Inhibitors.|Drugs..Perturbagen.Signatures.|Pathways..Hallmark.Gene.Sets.|Pathways..KEGG.Pathways.|Pathways..Nature.PID.Pathways.",colnames(cDT),value=TRUE), with=FALSE]
setkey(mDT,Barcode,Well)
wDT <- mDT[wDT, mult="first"]

ggplot(wDT, aes(x=Column, y=Row, colour=WellCellCount))+
  geom_point()+
  facet_wrap(~Barcode)

# #Calculate CVs for each set of replicates in the ScanR data
# cvNames<-grep(pattern="(Intensity|Area|SpotCellCount|Population|EdUPositivePercent|Barcode|^Name$|^Well$)",x=names(slDT),value=TRUE)
# cvKeep<-slDT[,cvNames,with=FALSE]
# repSpots<-c('Name','Well','Barcode')
# cv<-cvKeep[,lapply(.SD,CV),by=repSpots]
# data.table::setnames(cv,colnames(cv), paste0("CV.",colnames(cv)))
# data.table::setkey(cv,CV.Well,CV.Name, CV.Barcode)
# data.table::setkey(slDT,Well,Name,Barcode)
# slDT <- slDT[cv]

#Summarize well data to replicate level by taking the medians of these parameters
repNames<-grep(pattern="(Total|Mean|Area|EdUPositiveProportion|Cyto|Density|Z|X53BP1|WellCellCount|Barcode|^Spot$|GeneSymbol)",x=names(slDT),value=TRUE)
repKeep<-slDT[,repNames,with=FALSE]
repDT<-repKeep[,lapply(.SD,numericMedian),keyby="Barcode,GeneSymbol"]

#Merge back in the well and plate metadata
mDT <- slDT[,list(CellLine,WellIndex,Row,Column,GeneSymbol,EndpointDAPI,Endpoint488,Endpoint555,Endpoint647,Class,Drugs..FDA.Approved.,Drugs..Kinase.Inhibitors.,Drugs..Perturbagen.Signatures.,Pathways..Hallmark.Gene.Sets.,Pathways..KEGG.Pathways.,Pathways..Nature.PID.Pathways.),keyby="Barcode,GeneSymbol"]
repDT <- mDT[repDT, mult="first"]

repDT <- repDT[,WellCellCountRobustZScore := robustZScores(WellCellCount), by=Barcode]
repDT <- repDT[,CytoMeanIntensityAlexa488RobustZScore := robustZScores(CytoMeanIntensityAlexa488), by=Barcode]
repDT <- repDT[,CytoMeanIntensityAlexa555RobustZScore := robustZScores(CytoMeanIntensityAlexa555), by=Barcode]

# #Calculate CVs for each set of replicates in the ScanR data
# cvNames<-grep(pattern="(Intensity|Area|SpotCellCount|Population|EdUPositivePercent|Barcode|^Name$|^Well$)",x=names(slDT),value=TRUE)
# cvKeep<-slDT[,cvNames,with=FALSE]
# repSpots<-c('Name','Well','Barcode')
# cv<-cvKeep[,lapply(.SD,CV),by=repSpots]
# data.table::setnames(cv,colnames(cv), paste0("CV.",colnames(cv)))
# data.table::setkey(cv,CV.Well,CV.Name, CV.Barcode)
# data.table::setkey(slDT,Well,Name,Barcode)
# slDT <- slDT[cv]
# 
# Cleanup and write
# write.table(cDT, paste0(ss,"/Cell/Annotated Data/",ss,"_","CellAnn.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
# write.table(wDT, paste0(ss,"/Cell/Annotated Data/",ss,"_","CellWellAnn.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
# write.table(repDT, paste0(ss,"/Cell/Annotated Data/",ss,"_","CellRepAnn.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
