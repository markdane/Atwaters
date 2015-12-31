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
  CSCtrlAWells <- unique(dt$Well[is.na(dt$GeneSymbol)])
  CSCtrlBWells <- c("B01","B24","F01","F24","J01","J24","N01","N24")
  CSCtrlCWells <- c("C01","C24","G01","G24","K01","K24","O01","O24")
  CSCtrlDWells <- c("D01","D24","H01","H24","L01","L24","P01","P24")
  
  dt$GeneSymbol[dt$Well %in% CSCtrlAWells] <- "CSCtrlA"
  dt$GeneSymbol[dt$Well %in% CSCtrlBWells] <- "CSCtrlB"
  dt$GeneSymbol[dt$Well %in% CSCtrlCWells] <- "CSCtrlC"
  dt$GeneSymbol[dt$Well %in% CSCtrlDWells] <- "CSCtrlD"
  return(dt$GeneSymbol)
}

# normToNegCtrl <- function(dt){
#   WCCNegCtrl <- median(dt$WellCellCount[dt$GeneSymbol=="NegCtrl"],na.rm=TRUE)+1
#   return(dt$WellCellCount/WCCNegCtrl)
# }

normToNegCtrl <- function(dt){
  negCtrlMedian <- median(unlist(dt[,1,with=FALSE][dt$GeneSymbol=="NegCtrl"]),na.rm=TRUE)
  return(dt[,1,with=FALSE]/negCtrlMedian)
}

normToLogNegCtrl <- function(dt){
  negCtrlMedian <- median(unlist(dt[,1,with=FALSE][dt$GeneSymbol=="NegCtrl"]),na.rm=TRUE)+.01
  return(dt[,1,with=FALSE]-negCtrlMedian)
}

#Todo Move kmeansDNACluster fix into MEMA package 
kmeansDNACluster <- function (x, centers = 2) 
{
  if(length(x) == 1) return(1L)
  x <- data.frame(x)
  xkmeans <- kmeans(x, centers = centers)
  if (centers == 2) {
    if (xkmeans$centers[1] > xkmeans$centers[2]) {
      tmp <- xkmeans$cluster == 1
      xkmeans$cluster[xkmeans$cluster == 2] <- 1L
      xkmeans$cluster[tmp] <- 2L
    }
  }
  return(xkmeans$cluster)
}

#Count the number of neighboring cells
countNeighbors <- function(dt, radius){
  distMatrix <- as.matrix(dist(dt))
  count <- apply(distMatrix, 2, function(x) {
    sum(x <= radius) - 1
  })
  return(count)
}

labelPerimeterCells <- function(x){
  #browser()
  if(!length(x)==0){
    perimeterLogicals <- vector(length=length(x))
    perimeterLogicals[which.max(x)] <- TRUE
  }
  return(perimeterLogicals)
}

#########
rawDataVersion <- "v1.0"
calcNeighbors <- FALSE
neighborsThresh <- 5
wedgeAngs <- 18

densityThresh <- 0.4
outerThresh <- 0.5
ss <- "lineageEdU"

#Get the raw data file names
rdf <- getWellSubsetFileNames("RawData")

if(calcNeighbors){
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
  
  #Label cells DNA 2N, 4N and EdU state
  #Gate each well DAPI signal independently
  #Set 2N and 4N DNA status
  cDT <- cDT[,DNA2N := kmeansDNACluster(TotalIntensityDAPI), by="Barcode,Well"]
  
  cDT <- cDT[,DNA2NProportion := calc2NProportion(DNA2N),by="Barcode,Well"]
  cDT$DNA4NProportion <- 1-cDT$DNA2NProportion
  
  #Logit transform DNA Proportions
  #logit(p) = log[p/(1-p)]
  if(any(grepl("DNA2NProportion",colnames(cDT)))){
    DNA2NImpute <- cDT$DNA2NProportion
    DNA2NImpute[DNA2NImpute==0] <- .01
    DNA2NImpute[DNA2NImpute==1] <- .99
    cDT$DNA2NProportionLogit <- log2(DNA2NImpute/(1-DNA2NImpute))
  }
  
  if(any(grepl("DNA4NProportion",colnames(cDT)))){
    DNA4NImpute <- cDT$DNA4NProportion
    DNA4NImpute[DNA4NImpute==0] <- .01
    DNA4NImpute[DNA4NImpute==1] <- .99
    cDT$DNA4NProportionLogit <- log2(DNA4NImpute/(1-DNA4NImpute))
  }
  
  cDT <- cDT[,EduPositive := kmeansDNACluster(MeanIntensityAlexa647)-1L, by="Barcode,Well"]
  #Calculate the EdU Positive Percent at each spot
  cDT <- cDT[,EduPositiveProportion := sum(EduPositive)/length(EduPositive),by="Barcode,Well"]
  #Logit transform EduPositiveProportion
  #logit(p) = log[p/(1-p)]
  EdUppImpute <- cDT$EduPositiveProportion
  EdUppImpute[EdUppImpute==0] <- .01
  EdUppImpute[EdUppImpute==1] <- .99
  cDT$EduPositiveProportionLogit <- log2(EdUppImpute/(1-EdUppImpute))
  
  #Add in imageIDs
  
  #Add in adjacency parameters
  cDT <- cDT[,XLocal := (X - median(X, na.rm=TRUE)), by="Barcode,Well"]
  cDT <- cDT[,YLocal := (Y - median(Y, na.rm=TRUE)), by="Barcode,Well"]
  cDT <- cDT[,RadialPosition := sqrt(XLocal^2 + YLocal^2), by="Barcode,Well"]
  cDT <- cDT[,Theta := calcTheta(XLocal, YLocal), by="Barcode,Well"]
  
  
  cDTL <- mclapply(unique(cDT$Barcode), function(barcode, dt, nrRadii=5){
    setkey(dt,Barcode)
    bdt <- dt[barcode]
    neighborhoodNucleiRadii <-sqrt(median(bdt$Area/pi, na.rm=TRUE))
    bdt <- bdt[,Neighbors := countNeighbors(.SD, radius = nrRadii*neighborhoodNucleiRadii), by = "Well", .SDcols=c("XLocal","YLocal")]
  }, dt=cDT, mc.cores=detectCores())
  cDT <- rbindlist(cDTL)
} else {
  cat("loading cDT from Disk...")
  load("cDT.RData")
}

#Rules for classifying perimeter cells
cDT <- cDT[,Sparse := Neighbors < neighborsThresh]

#Add a local wedge ID to each cell based on conversations with Michel Nederlof
cDT <- cDT[,Wedge:=ceiling(Theta/wedgeAngs)]

#Define the perimeter cell if it exists in each wedge
#Classify cells as outer if they have a radial position greater than a thresh
cDT <- cDT[,OuterCell := labelOuterCells(RadialPosition, thresh=outerThresh),by="Barcode,Well"]

#Classify the perimeter cells
cDT <- cDT[,PerimeterCell:=labelPerimeterCells(RadialPosition),by="Barcode,Well,Wedge"]
#Require a perimeter cell not be in a sparse region
cDT$PerimeterCell[cDT$Sparse] <- FALSE

#Add a cell cycle state and summarize proportions
cDT$CellCycleState <- 2^(cDT$DNA2N)
cDT$CellCycleState[cDT$EduPositive==1] <- 3
cDT$CellCycleState <- cDT$CellCycleState-1
cDT$CellCycleState <- factor(cDT$CellCycleState,levels=c(1,2,3), labels=c("G0G1","S","G2"))

#Summarize cell data to well level by taking the medians of these parameters
wNames<-grep(pattern="(Intensity|Area|ElongationFactor|^Perimeter$|Lineage|Proportion|Neighbors|Density|WellCellCount|Barcode|^Well$)",x=names(cDT),value=TRUE)
wKeep<-cDT[,wNames,with=FALSE]
wDT<-wKeep[,lapply(.SD,numericMedian),keyby="Barcode,Well"]

#Add proportions for cell cycle state
wDT$G0G1Proportion <- cDT[,sum(CellCycleState=="G0G1")/.N,by="Barcode,Well"][,V1]
wDT$SProportion <- cDT[,sum(CellCycleState=="S")/.N,by="Barcode,Well"][,V1]
wDT$G2Proportion <- cDT[,sum(CellCycleState=="G2")/.N,by="Barcode,Well"][,V1]

#Merge back in the well and plate metadata
mDT <- unique(cDT[,c("Barcode","Well",setdiff(colnames(cDT),c(wNames,"X","Y","Position","DNA2N","EduPositive","CellCycleState","XLocal","YLocal","RadialPosition","Theta","Sparse","Wedge","OuterCell","PerimeterCell"))), with=FALSE], by=NULL)
setkey(mDT,Barcode,Well)
wDT <- mDT[wDT]

#Normalize WellCellCount to the NegCtrls
wDT <- wDT[,WellCellCountNorm := normToNegCtrl(.SD), by="Barcode", .SDcols=c("WellCellCount","GeneSymbol")]

#Normalize LineageRatio to the NegCtrls
wDT <- wDT[,LineageRatioNorm := normToLogNegCtrl(.SD), by="Barcode", .SDcols=c("LineageRatio","GeneSymbol")]

#Normalize EduPositiveProportionLogit to the NegCtrls
wDT <- wDT[,EduPositiveProportionLogitNorm := normToLogNegCtrl(.SD), by="Barcode", .SDcols=c("EduPositiveProportionLogit","GeneSymbol")]

#Summarize well data to replicate level by taking the medians of these parameters
repNames<-grep(pattern="(Intensity|Area|ElongationFactor|Perimeter|Lineage|Proportion|Neighbors|Density|WellCellCount|Barcode|GeneSymbol)",x=names(wDT),value=TRUE)
repKeep<-wDT[,repNames,with=FALSE]
repDT<-repKeep[,lapply(.SD,numericMedian),keyby="Barcode,GeneSymbol"]

#Merge back in the well and plate metadata
mDT <- unique(wDT[,c("Barcode","GeneSymbol",setdiff(colnames(wDT),c(repNames,"Well","Row","Column"))), with=FALSE],by=NULL)
#Combine negative controls from different pools
setkey(mDT,Barcode,GeneSymbol)
repDT <- mDT[repDT, mult="first"]

#Write data to disk
write.table(cDT, paste0("AtwatersCell_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
write.table(wDT, paste0("AtwatersWell_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
write.table(repDT, paste0("AtwatersRep_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
