#Read cell level data from disk
rawDataVersion <- "v1.0"
load("cDTL.RData")


#Add in adjacency parameters
cDT <- cDT[,XLocal := (X - median(X, na.rm=TRUE)), by="Barcode,Well"]
cDT <- cDT[,YLocal := (Y - median(Y, na.rm=TRUE)), by="Barcode,Well"]
cDT <- cDT[,RadialPosition := sqrt(XLocal^2 + YLocal^2), by="Barcode,Well"]
cDT <- cDT[,Theta := calcTheta(XLocal, YLocal), by="Barcode,Well"]

if(calcNeighbors){
  cDTL <- mclapply(unique(cDT$Barcode), function(barcode, dt, nrRadii=5){
    setkey(dt,Barcode)
    bdt <- dt[barcode]
    neighborhoodNucleiRadii <-sqrt(median(bdt$Area/pi, na.rm=TRUE))
    bdt <- bdt[,Neighbors := countNeighbors(.SD, radius = nrRadii*neighborhoodNucleiRadii), by = "Well", .SDcols=c("XLocal","YLocal")]
  }, dt=cDT, mc.cores=detectCores())
  
  cDT <- rbindlist(cDTL)
}
