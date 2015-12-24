
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