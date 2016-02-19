#Create analysis Output File
#Mark Dane 1/2015
library(data.table)
wDT <- fread("AtwatersWell_v1.1.txt")
sDT <- wDT[wDT$CellLine %in% c("CAMA1","HCC1143","HCC1419","CAL51","SKBR3","JIMT1","Hs578T","SUM149PT","SUM159PT","MDAMB157","MDAMB231","ZR751","ZR75B","AU565","T47D","HCC38","HCC1806"),]
sDT <- unique(sDT[,list(Barcode,Well,Plate,GeneSymbol,GeneAccession,GINumber,PoolCatalogNumber,CellLine,PlateID,Row,Column,MeanIntensityDAPI,MeanIntensityAlexa647,ElongationFactor,Perimeter,Area,WellCellCount,TotalIntensityDAPI,DNA2NProportion,DNA4NProportion,DNA2NProportionLogit,DNA4NProportionLogit,EduPositiveProportion,EduPositiveProportionLogit,G0G1Proportion,SProportion,G2Proportion,WellCellCountNorm,EduPositiveProportionLogitNorm)])

write.table(sDT,file="AtwatersSelectedWell.txt",quote = FALSE,sep = "\t",row.names = FALSE)
