potentially.somatic<-DFT1.DFT2.potentially.somatic.SNVs[,1:4]

potentially.somatic$STR<-"."

potentially.somatic<-potentially.somatic[,c(1,2,5,3,4)]

write.table(potentially.somatic, file = "potentially.somatic.txt", row.names = FALSE,
            col.names=FALSE, sep = "\t", quote = FALSE)