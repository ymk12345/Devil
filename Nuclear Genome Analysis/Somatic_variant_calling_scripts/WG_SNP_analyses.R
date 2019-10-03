

data<-read.table("test.genotyped.vcf", stringsAsFactors = FALSE, quote="\"")
colnames(data)<-read.table("test.genotyped.vcf", stringsAsFactors = FALSE, skip = 47, header = F, nrow = 1)


# FILTER TUMORS BASED ON HOST
## Take NV field for each variant (row) and sample (col)
variants.supp.reads = as.data.frame(apply(datafile[,10:ncol(datafile)], 2, function(var){as.numeric(str_split_fixed(var, ":", 6)[,6])}))
colnames(variants.supp.reads)<-paste0(colnames(datafile)[10:ncol(datafile)])


# Remove any variants observed in the host:
merged.datafile<-cbind(datafile[,1:9], variants.supp.reads)


# Exclude any variants that is also observed in the host
host.idx<-grep("H", colnames(merged.datafile)[10:ncol(merged.datafile)])+9
merged.datafile[,host.idx]<-as.numeric(as.character(merged.datafile[,host.idx]))

merged.datafile<-merged.datafile[-which(merged.datafile[,host.idx]>3),]













