# Create Docplots:
setwd("/Users/yk1/Documents/Docplots_for_Liz/")

chrs<-unique(tumors.pcf[,1])

#new.cnvs: cnv table

for(y in unique(new.cnvs$Total.Samples)){
  par(mfrow=c(4,2))
  pdf(paste0(y, ".pdf"), height = 10, width = 30)
  par(mfrow=c(4,2))
  y.pcf<-tumors.pcf[,c(1,2, which(gsub("X", "", colnames(tumors.pcf)) %in% y))]
  
  y.cnv.index<-new.cnvs[which(new.cnvs$Total.Samples %in% y),]
  plot.new()
  
  for(z in chrs){
    plot(y.pcf$pos[which(y.pcf$chrom == z)], 
         y.pcf[which(y.pcf$chrom == z),3], 
         ylab = "LogR", xlab = paste0(y, " Chrom", z),
         ylim = c(-2,2),
         pch = 20, 
         cex = .3)
    
    if(length(which(y.cnv.index$Chrom==z))>0){
      abline(v = y.cnv.index$Start[which(y.cnv.index$Chrom==z)],
             col = "red"
      )
      abline(v = y.cnv.index$End[which(y.cnv.index$Chrom==z)],
             col = "blue"
      )
      #par(opar)
      text(x = y.cnv.index$Start[which(y.cnv.index$Chrom==z)], cex = .7, 
           y = par("usr")[3] - .3, srt = 45, adj = 1,
           labels = y.cnv.index$Numeric[which(y.cnv.index$Chrom==z)], xpd = TRUE)
    }
    
  }
  dev.off()
}
