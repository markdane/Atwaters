#12/2015 preprocess scanr data
library("MEMA")
library("limma")
library("data.table")
library("parallel")
library(XLConnect)
library(ggplot2)

source("WellPlateFunctions.R")

#########
rawDataVersion <- "v1.1"
loadcDT <- FALSE
calcNeighbors <- FALSE
mergeOmeroIDs <- TRUE
neighborsThresh <- 5
wedgeAngs <- 36

densityThresh <- 0.4
outerThresh <- 0.5
ss <- "lineageEdU"

#Get the raw data file names
rdf <- getWellSubsetFileNames("RawData")
#rdf <- rdf[grep("MDAMB134VI",rdf$FilePaths),]
imageURLFiles <- grep("imageIDs",dir(paste0("./Metadata"),full.names = TRUE), value=TRUE)

if(!loadcDT){
  #Combine raw data from each plate
  cL <- mclapply(unique(rdf$Barcode), function(barcode, df){
    bdf <- df[df$Barcode==barcode,]
    bDT <- stitchWellData(fs=bdf)
    if(mergeOmeroIDs){
      #Read in and merge the Omero URLs
      omeroIndex <- fread(grep(barcode, imageURLFiles, value=TRUE))[,list(PlateID,Row,Column,Field,ImageID)]
      omeroIndex$Well <- omeroIndex$Column + (omeroIndex$Row-1)*24
      bDT$Field <- bDT$Position-1
      bDT <- merge(bDT,omeroIndex,by=c("Well","Field"))
    }
    return(bDT)
  }, df=rdf, mc.cores = detectCores())
  
  cDT <- rbindlist(cL)
  #Convert well index to an alphanumeric label
  cDT$Well <- wellAN(16,24)[cDT$Well]
  
  #Filter out debris based on nuclear area
  nuclearAreaThresh <- 200
  cDT <- cDT[cDT$Area >nuclearAreaThresh,]
  #Count the cells at each well
  cDT<-cDT[,WellCellCount := .N, by="Barcode,Well"]
  
  cDT$TotalIntensityDAPI <- cDT$Area*cDT$MeanIntensityDAPI
  
  #Read in siRNA location and annotations
  siRNAs <- readWorksheetFromFile("./Metadata/LH_2015Oct384 G-CUSTOM-185752.xls", sheet="Sequences_by_Ord_w384ShippingRa",startRow=3, header=TRUE)
  siRNAs <- siRNAs[-which(is.na(siRNAs)[,1]),]
  siRNAs <- convertColumnNames(data.table(siRNAs))
  siRNAs$Plate <- gsub("Plate ","",siRNAs$Plate)
  siRNAs <- unique(siRNAs[,list(Plate,Well,GeneSymbol,GENEID,GeneAccession,GINumber,PoolCatalogNumber)])
  siRNAs$GeneSymbol[grepl("ON-TARGET",siRNAs$GeneSymbol)] <- "NegCtrl"
  
  #Add siRNA annotations to the cell level data
  setkey(siRNAs,Plate,Well)
  setkey(cDT,Plate,Well)
  cDT <- siRNAs[cDT]
  
  #Add labels to cell seeding wells
  GSdt <- unique(annotateCellSeedWells(cDT[,list(GeneSymbol,Well,Plate)]))
  setkey(GSdt,Plate,Well)
  setkey(cDT,Plate,Well)
  cDT <- cDT[,GeneSymbol:=NULL]
  cDT <- merge(cDT,GSdt)
  
  #Create a lineageRatio signal for each cell KRT19/KRT5 luminal/basal
  cDT$LineageRatio <- log2((cDT$MeanIntensityAlexa555Cyto+1)/(cDT$MeanIntensityAlexa488Cyto+1))
  
  #Label cells DNA 2N, 4N and EdU state
  #Gate each well DAPI signal independently
  #Set 2N and 4N DNA status
  #   cDT <- cDT[,DNA2N := kmeansDNACluster(TotalIntensityDAPI), by="Barcode"]
  #   cDT <- cDT[,DNA2NProportion := calc2NProportion(DNA2N),by="Barcode"]
  #   cDT$DNA4NProportion <- 1-cDT$DNA2NProportion
  cDT <- cDT[,DNA2N := gateOnlocalMinima(TotalIntensityDAPI,probs=c(.05,.8)), by="Barcode"]
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
  
  #cDT <- cDT[,EduPositive := kmeansDNACluster(MeanIntensityAlexa647)-1L, by="Barcode"]
  cDT <- cDT[,EduPositive := gateOnlocalMinima(MeanIntensityAlexa647, probs=c(.05,.95))-1, by="Barcode"]
  #Calculate the EdU Positive Percent at each spot
  cDT <- cDT[,EduPositiveProportion := sum(EduPositive)/length(EduPositive),by="Barcode,Well"]
  #Logit transform EduPositiveProportion
  #logit(p) = log[p/(1-p)]
  EdUppImpute <- cDT$EduPositiveProportion
  EdUppImpute[EdUppImpute==0] <- .01
  EdUppImpute[EdUppImpute==1] <- .99
  cDT$EduPositiveProportionLogit <- log2(EdUppImpute/(1-EdUppImpute))
  
  
  #Add in adjacency parameters
  cDT <- cDT[,XLocal := (X - median(X, na.rm=TRUE)), by="Barcode,Well"]
  cDT <- cDT[,YLocal := (Y - median(Y, na.rm=TRUE)), by="Barcode,Well"]
  cDT <- cDT[,RadialPosition := sqrt(XLocal^2 + YLocal^2), by="Barcode,Well"]
  cDT <- cDT[,Theta := calcTheta(XLocal, YLocal), by="Barcode,Well"]
  
  if(calcNeighbors){
    cDTL <- lapply(unique(cDT$Barcode), function(barcode, dt, nrRadii=5){
      setkey(dt,Barcode)
      bdt <- dt[barcode]
      neighborhoodNucleiRadii <-sqrt(median(bdt$Area/pi, na.rm=TRUE))
      bdt <- bdt[,Neighbors := countNeighbors(.SD, radius = nrRadii*neighborhoodNucleiRadii), by = "Well", .SDcols=c("XLocal","YLocal")]
    }, dt=cDT)
    cDT <- rbindlist(cDTL)
  }
  
} else {
  cat("loading cDT from Disk...")
  load("cDT.RData")
}


#Rules for classifying perimeter cells
if("Neighbors" %in% colnames(cDT)) cDT <- cDT[,Sparse := Neighbors < neighborsThresh]

#Add a local wedge ID to each cell based on conversations with Michel Nederlof
cDT <- cDT[,Wedge:=ceiling(Theta/wedgeAngs)]

#Define the perimeter cell if it exists in each wedge
#Classify cells as outer if they have a radial position greater than a thresh
cDT <- cDT[,OuterCell := labelOuterCells(RadialPosition, thresh=outerThresh),by="Barcode,Well"]

#Classify the perimeter cells
cDT <- cDT[,PerimeterCell:=labelPerimeterCells(RadialPosition),by="Barcode,Well,Wedge"]
#Require a perimeter cell not be in a sparse region
if("Sparse" %in% colnames(cDT)) cDT$PerimeterCell[cDT$Sparse] <- FALSE

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
setkey(cDT,Barcode,Well)
setkey(wDT,Barcode,Well)
wDT$G0G1Proportion <- cDT[,sum(CellCycleState=="G0G1")/.N,by="Barcode,Well"][,V1]
wDT$SProportion <- cDT[,sum(CellCycleState=="S")/.N,by="Barcode,Well"][,V1]
wDT$G2Proportion <- cDT[,sum(CellCycleState=="G2")/.N,by="Barcode,Well"][,V1]

#Merge back in the well and plate metadata
mDT <- unique(cDT[,c("Barcode","Well",setdiff(colnames(cDT),c(wNames,"X","Y","Position","DNA2N","EduPositive","CellCycleState","XLocal","YLocal","RadialPosition","Theta","Sparse","Wedge","OuterCell","PerimeterCell"))), with=FALSE], by=NULL)
setkey(mDT,Barcode,Well)
wDT <- mDT[wDT]

#Lormalize the signals within each plate by subtracting the fitted loess value
wDT <- loessNormalize(dt=wDT,signalNames=c("WellCellCount"), span=.05)

#Normalize WellCellCount to the NegCtrls
wDT <- wDT[,WellCellCountNorm := normToNegCtrl(.SD), by="Barcode", .SDcols=c("WellCellCountLoessNorm","GeneSymbol")]

#Normalize LineageRatio to the NegCtrls
wDT <- wDT[,LineageRatioNorm := normToLogNegCtrl(.SD), by="Barcode", .SDcols=c("LineageRatio","GeneSymbol")]

#Normalize EduPositiveProportionLogit to the NegCtrls
wDT <- wDT[,EduPositiveProportionLogitNorm := normToLogNegCtrl(.SD), by="Barcode", .SDcols=c("EduPositiveProportionLogit","GeneSymbol")]

#Summarize well data to replicate level by taking the medians of these parameters
repNames<-grep(pattern="(Intensity|Area|ElongationFactor|Perimeter|Lineage|Proportion|Neighbors|Density|WellCellCount|Barcode|GeneSymbol)",x=names(wDT),value=TRUE)
repKeep<-wDT[,repNames,with=FALSE]
repDT<-repKeep[,lapply(.SD,numericMedian),keyby="Barcode,GeneSymbol"]

#Merge back in the well and plate metadata
mDT <- unique(wDT[,c("Barcode","GeneSymbol",setdiff(colnames(wDT),c(repNames,"Well","Row","Column","ImageID"))), with=FALSE],by=NULL)
#Combine negative controls from different pools
setkey(mDT,Barcode,GeneSymbol)
repDT <- mDT[repDT, mult="first"]

#Write data to disk
write.table(cDT, paste0("AtwatersCell_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
write.table(wDT, paste0("AtwatersWell_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
write.table(repDT, paste0("AtwatersRep_",rawDataVersion,".txt"), sep = "\t",row.names = FALSE, quote=FALSE)
