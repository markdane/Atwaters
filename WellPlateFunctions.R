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
  #Identify NA wells on a plate basis
  dt <- rbindlist(lapply(unique(dt$Plate), function(plate){
    setkey(dt,Plate)
    pdt <- dt[plate]
    CSCtrlAWells <- unique(pdt$Well[is.na(pdt$GeneSymbol)])
    CSCtrlBWells <- c("B01","B24","F01","F24","J01","J24","N01","N24")
    CSCtrlCWells <- c("C01","C24","G01","G24","K01","K24","O01","O24")
    CSCtrlDWells <- c("D01","D24","H01","H24","L01","L24","P01","P24")
    
    pdt$GeneSymbol[pdt$Well %in% CSCtrlAWells] <- "CSCtrlA"
    pdt$GeneSymbol[pdt$Well %in% CSCtrlBWells] <- "CSCtrlB"
    pdt$GeneSymbol[pdt$Well %in% CSCtrlCWells] <- "CSCtrlC"
    pdt$GeneSymbol[pdt$Well %in% CSCtrlDWells] <- "CSCtrlD"
    return(pdt)
  }))
  return(dt)
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

numericMedian <- function(x) as.numeric(median(x))

evalMedians <- function(values, reps){
  tmp <- median(values[reps], na.rm=TRUE)
}

create96WellHistograms <- function(DT, pr, prDisplay, binwidth = diff(quantile(DT[[pr]],probs = c(0,.98),na.rm=TRUE))/50, upperProb = .99, ncol = 4) {
  p <- ggplot(DT, aes_string(x=pr))+
    geom_histogram(binwidth = binwidth)+
    scale_x_continuous(limits = quantile(DT[[pr]],probs = c(0,upperProb),na.rm=TRUE))+
    ggtitle(paste(prDisplay,"in",unique(DT$CellLine), "cells in plate",unique(DT$Barcode)))+
    xlab(prDisplay)+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(.5)), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=rel(1)), plot.title = element_text(size = rel(.5)),legend.text=element_text(size = rel(.3)),legend.title=element_text(size = rel(.3)))
}


normWithinPlate <- function (DT, value, base) 
{
  if (!c("GeneSymbol") %in% colnames(DT)) 
    stop(paste("DT must contain a GeneSymbol column."))
  valueMedian <- median(unlist(DT[(grepl(base, DT$GeneSymbol)), value, with = FALSE]),na.rm = TRUE)
  normedValues <- DT[, value, with = FALSE]/valueMedian
  return(normedValues)
}

#Compute Z' Factors
calcZPrime <- function(dt, method="original", posCtrl="PLK1"){
  imageHTS::zprime(dt[GeneSymbol=="NegCtrl",WellCellCount],dt[GeneSymbol==posCtrl,WellCellCount], method)
}

#Loess normalization of WellCellCount
loessNormalize <- function(dt, signalNames, span=1){
  #Normalize a data.table of signals 
  #There must be barcode and well columns
  if(!(any(grepl("Barcode",colnames(dt)))&any(grepl("Row",colnames(dt)))&any(grepl("Column",colnames(dt))))) stop("There must be Barcode, Column and Row and columns")
  
  dtbList <-lapply(unique(dt$Barcode),function(barcode){
    setkey(dt,Barcode)
    dtb <- dt[barcode]
    #ToDo: Either make the code general or specific
    sigNormedList <- lapply(signalNames, function(signalName, dtb){
      #get mtarget values for each gene symbol
      setnames(dtb,signalName,"Value")
      dtb <- dtb[,mt := as.numeric(median(Value, na.rm=TRUE)), by=GeneSymbol]
      dtb <- dtb[,WCCRes := as.numeric(Value-mt)]
      dtb$WCCResFitted <-loess(WCCRes~Row+Column, dtb, span=span)$fitted
      dtb <- dtb[,WellCellCountLoessNorm := as.numeric(Value-WCCResFitted)]
      setnames(dtb,"Value",signalName)
      return(dtb)
    },dtb=dtb)
    dtb<- rbindlist(sigNormedList)
  })
  dt <-rbindlist(dtbList)
  return(dt)
}

#Calculate SSMD
SSMDPlate <- function(df){
  if(length(unique(df)))
    mup <- mean(df$WellCellCount[df$GeneSymbol=="PLK1"], na.rm = TRUE)
  sdp <- sd(df$WellCellCount[df$GeneSymbol=="PLK1"], na.rm = TRUE)
  mun <- mean(df$WellCellCount[df$GeneSymbol=="NegCtrl"], na.rm = TRUE)
  sdn <- sd(df$WellCellCount[df$GeneSymbol=="NegCtrl"], na.rm = TRUE)
  SSMD <- (mup-mun)/sqrt(sdp^2+sdn^2)
}

#Calculate SSMD on the loess normed well cell count
SSMDPlateLoess <- function(df){
  if(length(unique(df)))
    mup <- mean(df$WellCellCountLoessNorm[df$GeneSymbol=="PLK1"], na.rm = TRUE)
  sdp <- sd(df$WellCellCountLoessNorm[df$GeneSymbol=="PLK1"], na.rm = TRUE)
  mun <- mean(df$WellCellCountLoessNorm[df$GeneSymbol=="NegCtrl"], na.rm = TRUE)
  sdn <- sd(df$WellCellCountLoessNorm[df$GeneSymbol=="NegCtrl"], na.rm = TRUE)
  SSMD <- (mup-mun)/sqrt(sdp^2+sdn^2)
}

calcZPrimeLoess <- function(dt, method="original", posCtrl="PLK1"){
  imageHTS::zprime(dt[GeneSymbol=="NegCtrl",WellCellCountLoessNorm],dt[GeneSymbol==posCtrl,WellCellCountLoessNorm], method)
}


localMinima <- function(x, probs=c(.2,.8)){
  #browser()
  #Finds the local minima between the probs quantiles
  #x numeric vector
  #probs interval limits on where to search for the minima
  h <- hist(x,breaks=300, plot=FALSE)
  if(length(h$mids)<2) return(max(x))
  f <- approxfun(h$mids, h$counts)
  o <- optimise(f, interval=quantile(x, probs))
  if(length(o)>2) stop()
  return(o$minimum)
}


gateOnlocalMinima <- function(x, ...){
  thresh <- localMinima(x, ...)
  cluster <- rep.int(1,times=length(x))
  cluster[x>thresh] <- 2
  return(cluster)
}

addQAMetrics <- function(dt){
  
  #Remove rows from multiple fields/images of each well
  colNames <- grep("Omero|Image|Field",colnames(dt),value=TRUE,invert=TRUE)
  DT <- unique(dt[,colNames, with=FALSE], by=c("Barcode","Well"))
  
  DT <- DT[,ZPrimeRobust := calcZPrime(.SD, method="robust"), by="Barcode", .SDcols=c("GeneSymbol","WellCellCount")]
  DT <- DT[,SSMD := SSMDPlate(.SD), by="Barcode", .SDcols=c("GeneSymbol","WellCellCount")]
  DT <- DT[,ZPrimeRobustLoess := calcZPrimeLoess(.SD, method="robust"), by="Barcode", .SDcols=c("GeneSymbol","WellCellCountLoessNorm")]
  DT <- DT[,SSMDLoess := SSMDPlateLoess(.SD), by="Barcode", .SDcols=c("GeneSymbol","WellCellCountLoessNorm")]
  
  DT <- unique(DT[,list(Barcode,ZPrimeRobust,SSMD,ZPrimeRobustLoess,SSMDLoess)])
  dt <- merge(dt,DT, by="Barcode")
  
  return(dt)
}
