# Argument
matching<-commandArgs(trailingOnly=TRUE)[1]

# LIBRARIES
require(stringr)

datafile <- read.table("~/Documents/Tasmanian_Devil/WG_Analyses_variants/100_genotyped_contigs.revised.vcf", stringsAsFactors = FALSE, quote="\"");
colnames(datafile) <- read.delim("~/Documents/Tasmanian_Devil/WG_Analyses_variants/100_genotyped_contigs.revised.vcf", stringsAsFactors = FALSE, skip = 47, header = F, nrow = 1);
colnames(datafile)[1]<-"CHROM"

# FILTER TUMORS BASED ON HOST
## Take NV field for each variant (row) and sample (col)
variants.supp.reads = as.data.frame(apply(datafile[,10:ncol(datafile)], 2, function(var){as.numeric(str_split_fixed(var, ":", 6)[,6])}))
colnames(variants.supp.reads)<-paste0(colnames(datafile)[10:ncol(datafile)])


# Remove any variants observed in the host:
merged.datafile<-cbind(datafile[,1:9], variants.supp.reads)


# Exclude any variants that is also observed in the host
host.idx<-grep("H", colnames(merged.datafile)[10:ncol(merged.datafile)])+9
merged.datafile[,host.idx]<-as.numeric(as.character(merged.datafile[,host.idx]))

merged.datafile<-merged.datafile[-which(merged.datafile[,host.idx]>0),]
merged.datafile<-merged.datafile[-which(as.numeric(as.character(merged.datafile$QUAL))<10),]
merged.datafile<-merged.datafile[-grep("SC|badReads|MQ", as.character(merged.datafile$FILTER)),]

# Write out vcf file
directory<-commandArgs(trailingOnly=TRUE)[2]
setwd(directory)
write.table(datafile, file = paste0("filtered.tumor.host", matching), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

