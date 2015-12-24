#12/2015 preprocess scanr data
library("MEMA")
library("limma")
library("data.table")
library("parallel")
library(XLConnect)
library(ggplot2)

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
  rdFiles$Plate <- as.integer(strsplit2(rdFiles$Barcode,"")[,7])
  return(rdFiles)
}

stitchWellData <- function(fs){
  dtL <-apply(fs[grepl("_Main",fs$FilePaths),], 1, function(fa){
    mdt <- fread(fa["FilePaths"])
    cdt <- fread(gsub("Main","Cyto",fa["FilePaths"]))
    mdt$Barcode <- fa["Barcode"]
    mdt$CellLine <- fa["CellLine"]
    mdt$Plate <- fa["Plate"]
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

annotateCellSeedWells <- function(dt){
  # First and last columns with cells 1:2 serial dilution as 
  #A  CSCtrlA
  #B(1/2 #A) CSCtrlB
  #C(1/4 #A) CSCtrlC
  #D(1/8 #A) CSCtrlD
  #the rest of wells with same as A. 
  CSCtrlAWells <- c(paste0("A0",1:9),paste0("A",10:24),
                    paste0("P0",2:9),paste0("P",10:23),
                    "E01","E24","I01","I24","M01","M24")
  
  CSCtrlBWells <- c("B01","B24","F01","F24","J01","J24","N01","N24")
  CSCtrlCWells <- c("C01","C24","G01","G24","K01","K24","O01","O24")
  CSCtrlDWells <- c("D01","D24","H01","H24","L01","L24","P01","P24")
  
  dt$GeneSymbol[dt$Well %in% CSCtrlAWells] <- "CSCtrlA"
  dt$GeneSymbol[dt$Well %in% CSCtrlBWells] <- "CSCtrlB"
  dt$GeneSymbol[dt$Well %in% CSCtrlCWells] <- "CSCtrlC"
  dt$GeneSymbol[dt$Well %in% CSCtrlDWells] <- "CSCtrlD"
  return(dt$GeneSymbol)
}

#########
rawDataVersion <- "v1.0"

#Get the raw data file names
rdf <- getWellSubsetFileNames("RawData")

#Combine raw data from each plate
cDT <- rbindlist(lapply(unique(rdf$Barcode), function(barcode, df){
  bdf <- df[df$Barcode==barcode,]
  bDT <- stitchWellData(fs=bdf)
}, df=rdf))

#Add row and column indices
cDT$Row <- ceiling(cDT$Well/24)
cDT$Column <- (cDT$Well-1) %% 24 + 1

#Convert well index to an alphanumeric label
cDT$Well <- wellAN(16,24)[cDT$Well]

densityThresh <- 0.4
outerThresh <- 0.5
ss <- "lineageEdU"

#Count the cells at each well
cDT<-cDT[,WellCellCount := .N, by="Barcode,Well"]

#Filter out debris based on nuclear area
nuclearAreaThresh <- 50
cDT <- cDT[cDT$Area >nuclearAreaThresh,]

cDT$TotalIntensityDAPI <- cDT$Area*cDT$MeanIntensityDAPI

#Read in siRNA location and annotations
siRNAs <- readWorksheetFromFile("./Metadata/LH_2015Oct384 G-CUSTOM-185752.xls", sheet="Sequences_by_Ord_w384ShippingRa",startRow=3, header=TRUE)
siRNAs <- convertColumnNames(data.table(siRNAs))
siRNAs$Plate <- gsub("Plate ","",siRNAs$Plate)
siRNAs <- unique(siRNAs[,list(Plate,Well,GeneSymbol,GENEID,GeneAccession,GINumber,PoolCatalogNumber)])
siRNAs$GeneSymbol[grepl("ON-TARGET",siRNAs$GeneSymbol)] <- "NegCtrl"

#Add siRNA annotations to the cell level data
setkey(siRNAs,Plate,Well)
setkey(cDT,Plate,Well)
cDT <- siRNAs[cDT]

#Add labels to cell seeding wells
cDT$GeneSymbol <- annotateCellSeedWells(cDT[,list(GeneSymbol,Well)])

#Create a lineageRatio signal for each cell KRT19/KRT5 luminal/basal
cDT$LineageRatio <- log2((cDT$MeanIntensityAlexa555Cyto+1)/(cDT$MeanIntensityAlexa488Cyto+1))

#Summarize cell data to well level by taking the medians of these parameters
wNames<-grep(pattern="(Intensity|Area|ElongationFactor|Perimeter|Lineage|EdUPositiveProportion|Density|WellCellCount|Barcode|^Well$)",x=names(cDT),value=TRUE)
#Remove the well normalized values as the median is 1 by definition
wNames <- grep("WellNorm",x=wNames, value=TRUE, invert=TRUE)
wKeep<-cDT[,wNames,with=FALSE]
wDT<-wKeep[,lapply(.SD,numericMedian),keyby="Barcode,Well"]

#Merge back in the well and plate metadata
mDT <- cDT[,c("Barcode","Well",setdiff(colnames(cDT),wNames)), with=FALSE]
setkey(mDT,Barcode,Well)
wDT <- mDT[wDT, mult="first"]

#Summarize well data to replicate level by taking the medians of these parameters
repNames<-grep(pattern="(Intensity|Area|ElongationFactor|Perimeter|Lineage|EdUPositiveProportion|Density|WellCellCount|Barcode|GeneSymbol)",x=names(wDT),value=TRUE)
repKeep<-wDT[,repNames,with=FALSE]
repDT<-repKeep[,lapply(.SD,numericMedian),keyby="Barcode,GeneSymbol"]

#Merge back in the well and plate metadata
mDT <- wDT[,c("Barcode","GeneSymbol",setdiff(colnames(wDT),repNames)), with=FALSE]
setkey(mDT,Barcode,GeneSymbol)
repDT <- mDT[repDT, mult="first"]

#Write data to disk
write.table(cDT, paste0("AtwatersCell_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
write.table(wDT, paste0("AtwatersWell_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
write.table(repDT, paste0("AtwatersRep_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
