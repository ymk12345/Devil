# Variant Matching ID identification

vcf_filtered[,2]<-gsub("_", "", 
                   gsub("platypusVariants_", "", 
                   gsub("-Devil.*", "", as.character(vcf_filtered[,1]))))

vcf_filtered[,3]<-"Tumor"
vcf_filtered[grep("H", vcf_filtered[,2]),3]<-"Host"
vcf_filtered$ID<-gsub("[[:alpha:]].*", "", vcf_filtered[,2])

hosts<-vcf_filtered[grep("H", vcf_filtered[,2]),]
tumours<-vcf_filtered[grep("T", vcf_filtered[,2]),]


tumours$matching.host.exists <- match(gsub("[[:alpha:]].*", "", tumours[,2]), gsub("[[:alpha:]].*", "", hosts[,2]))

unmatched<-unique(tumours[which(is.na(tumours$matching.host.exists)),4])
matched<-unique(tumours[-which(is.na(tumours$matching.host.exists)),4])

write.table(unmatched, file = "unmatched.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


# already matched:
already_genotyped[,2]<-gsub(".gz", "", gsub("final_regenotyped_", "", gsub("_genotyped*.*", "", already_genotyped[,1])))
tumours.still.need.to.genotype<-matched[-which(matched %in% already_genotyped[,2])]
write.table(tumours.still.need.to.genotype, file = "need.to.genotype.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

unmatched.tumours.still.need.to.genotype<-unmatched[-which(unmatched %in% already_genotyped[,2])]
write.table(unmatched.tumours.still.need.to.genotype, file = "unmatched.need.to.genotype.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
