
# 1. Import and Load Relevant Data 

  # a) Metadata
    load("~/Desktop/Devil/drive-download-20190919T030557Z-001/manifest_total_2_10_18.RData")
  
  # b) Allele Count Data
    load("~/Desktop/Devil/drive-download-20190919T030557Z-001/allele.counter.data.updated.final.3.12.18.RData")
    load("~/Desktop/Devil/allele.counter.data.comb_2_10_18.RData")
  
  # c) Probe Metadata
    load("~/Desktop/Devil/mips_rebalanced_vep.RData")
    
  # d) Haplotype Database
    load("haplos.RData")

  # e) Libraries
    library(stringr)
    library(data.table)
    library(reshape2)

    
# 2. Transform Allele Count Data    
    allele.counter.data.comb$TAG.INDEX<-str_split_fixed(str_split_fixed(allele.counter.data.comb$Sample, "#", 2)[,2],  "[.]", 2)[,1]
    allele.counter.data.comb$seq_run<-str_split_fixed(str_split_fixed(allele.counter.data.comb$Sample, "#", 2)[,1], "_", 2)[,1]


    
# 3 Merge Allele Count and Manifest Metadata   
    #a. Merge Old Allele Count Data with Sample and MIPs Metadata
      # i) Sample Manifests
      allele.counter.data.comb1<-merge(allele.counter.data.comb, manifest_total_2_10_18)
      colnames(allele.counter.data.comb1)[7]<-"ALT"

      # ii). Merge Allele Count with MIPs
      allele.counter.data.comb1<-merge(allele.counter.data.comb1, mips_rebalanced_vep)
      allele.counter.data.comb1$Sample.TCG_ID<-gsub("[-]Devil", "", allele.counter.data.comb1$TCG_ID)
      allele.counter.data.comb1$TCG_ID<-gsub("[-]Devil", "", allele.counter.data.comb1$TCG_ID)


    #b. Merge New Allele Count with Manifests; reformat so the tables are the same.
      allele.counter.data.updated1<-allele.counter.data.updated1[,-c(31:32)]
      mip.ac.rebalanced1<-merge(allele.counter.data.updated1, mips_rebalanced_vep)
      mip.ac.rebalanced1$Sample.TCG_ID<-mip.ac.rebalanced1$Supplier.name
      mip.ac.rebalanced1$TCG_ID<-gsub("[-]Devil", "", mip.ac.rebalanced1$Sample.TCG_ID)
      mip.ac.rebalanced1$Sample.TCG_ID<-gsub("[-]Devil", "", mip.ac.rebalanced1$Sample.TCG_ID)
      mip.ac.rebalanced1$INTERMEDIATE.PLATE.ID<-""

    
      
# 5. Merge Old and New Allele Count Data; and reformat table so it's compatible with the haplotype data table
    allele.counter.data.comb1<-rbind(allele.counter.data.comb1[,which(colnames(allele.counter.data.comb1) %in% colnames(mip.ac.rebalanced1))], 
                                 mip.ac.rebalanced1[,which(colnames(mip.ac.rebalanced1) %in% colnames(allele.counter.data.comb1))])

    allele.counter.data.comb1$VAF<-as.numeric(as.character(allele.counter.data.comb1$value))/as.numeric(as.character(allele.counter.data.comb1$Depth))
    allele.counter.data.comb1$Group<-unlist(stringr::str_extract_all(allele.counter.data.comb1$mip_name, "\\([^()]+\\)"))
    allele.counter.data.comb1$Fraction<-paste0("'", allele.counter.data.comb1$value, "/", allele.counter.data.comb1$Depth)
    allele.counter.data.comb1$unique.ids<-allele.counter.data.comb1$TCG_ID
    allele.counter.data.comb1$id<-allele.counter.data.comb1$TCG_ID
    
    
    
# 6.  Merge Allele Count Data and Haplotype
    allele.counter.data.comb1<-cbind(allele.counter.data.comb1, 
                                     haplos[match(allele.counter.data.comb1$Sample.TCG_ID, haplos$Sample.TCG_ID),])
    
    
    
# 7. Check Samples for Quality Check    
    for(x in unique(allele.counter.data.comb1$seq_run)){
  pdf(paste0(x, "ac_check_hist.pdf"), height = 8, width = 20)
  #par(mfrow=c(2,1))
  temp<-allele.counter.data.comb1[which(allele.counter.data.comb1$seq_run ==x),]
  temp$TAG.INDEX<-as.numeric(as.character(temp$TAG.INDEX))
  temp<-temp[order(temp$TAG.INDEX),]
  temp<-temp[grep("T", temp$TCG_ID),]
  
  for(y in unique(temp$TCG_ID)){
    
    temp.samp<-temp[which(temp$TCG_ID %in% y),]
    med.vaf<- as.numeric(supp.table1[which(supp.table1[,1] %in% y),2])/2
    
    
    hist(temp.samp$VAF[which(temp.samp$VAF>0)], breaks = seq(0, 1, by = .01), xlab = "proportion of reads in plate", main = paste0(y, "  Median VAF: ", med.vaf))
    abline(v=med.vaf, col = "Red")
    
    
    
    
    #print(ggplot(temp.samp, aes(Group, VAF, group = Group))+
    #       labs(x = "", title = paste0(y, "  ", temp.samp$TAG.INDEX[1], "  ", haplos$Groups[match(temp.samp$TCG_ID[1], haplos$Sample.TCG_ID)], ":  ", 
    #                                    haplos$Complete.Haplotype[match(temp.samp$TCG_ID[1], haplos$Sample.TCG_ID)]))+
    #       geom_violin(fill = NA)+
    #      geom_jitter(width = .2)+
    #     theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylim(c(0,1))+
    #    geom_hline(aes(yintercept = med.vaf), col = "red"))
    
    
  }
  
  dev.off()
}



    
# 8. Combine Read Data for Samples that where libraries were recreated
    #a. Isolate Duplicate Samples from Sequencing Runs
      samples.to.combine<-data.frame(samp=unique(paste0(allele.counter.data.comb1$Sample.TCG_ID, "_", allele.counter.data.comb1$seq_run)))
      samples.to.combine$Sample.TCG_ID<-gsub("-Devil", "", str_split_fixed(samples.to.combine$samp, "_", 2)[,1])
      samples.to.combine<-samples.to.combine$Sample.TCG_ID[which(duplicated(samples.to.combine$Sample.TCG_ID))]
      samples.to.combine<-samples.to.combine[grep("T|H", samples.to.combine)]


    #b. Find Samples in Multiple Sequencing Run; Merge Data for Duplicated Samples
      allele.counter.data.comb1$Sample.TCG_ID<-gsub("-Devil", "", allele.counter.data.comb1$Sample.TCG_ID)
      duplicated.samples1<-allele.counter.data.comb1[which(allele.counter.data.comb1$Sample.TCG_ID %in% samples.to.combine),]
      duplicated.samples1.combined<-NULL
      for(y in samples.to.combine){
        temp.dup<-allele.counter.data.comb1[which(allele.counter.data.comb1$Sample.TCG_ID %in% y),]
        temp<-temp.dup
        t.val<-aggregate(temp[,c("value", "Depth")], by = list(temp$mip_name), FUN=sum)
        temp.dup$value<-t.val$value[match(temp.dup$mip_name, t.val$Group.1)]
        temp.dup$Depth<-t.val$Depth[match(temp.dup$mip_name, t.val$Group.1)]
        
        temp.dup$VAF<-temp.dup$value/temp.dup$Depth
        temp.dup$Fraction<-paste0("'", temp.dup$value, "/", temp.dup$Depth)
        #temp.dup$PLATE_ID<-paste(unique(temp.dup$PLATE_ID), collapse = "_")
        temp.dup$seq_run<-paste(unique(temp.dup$seq_run), collapse = "_")
        temp.dup$INTERMEDIATE.PLATE.ID<-paste(unique(temp.dup$INTERMEDIATE.PLATE.ID), collapse = "_")
        
        temp.dup<-unique(temp.dup)
        duplicated.samples1.combined<-rbind(duplicated.samples1.combined, temp.dup)
        
      }
      duplicated.samples1.combined<-unique(duplicated.samples1.combined)

      
    #c. Merge Samples Genotyped in Single and Multiple Runs
      allele.counter.data.comb1<-allele.counter.data.comb1[-which(allele.counter.data.comb1$Sample.TCG_ID %in% samples.to.combine),]
      mip.ac.rebalanced.d1.comb<-rbind(allele.counter.data.comb1, duplicated.samples1.combined)

      
    #d. Add Clade and Group Info to Allele Count
      mip.ac.rebalanced.d1.comb$Clade<-haplos$Clade[match(mip.ac.rebalanced.d1.comb$Sample.TCG_ID, haplos$Sample.TCG_ID)]
      mip.ac.rebalanced.d1.comb$Groups<-haplos$Groups[match(mip.ac.rebalanced.d1.comb$Sample.TCG_ID, haplos$Sample.TCG_ID)]
      mip.ac.rebalanced.d1.comb<-mip.ac.rebalanced.d1.comb[grep("T", mip.ac.rebalanced.d1.comb$Sample.TCG_ID),]
  
    
    #e. Reformat Allele Count Data and Calculate VAF
      mip.ac.rebalanced.d1.only<-mip.ac.rebalanced.d1.comb
      mip.ac.rebalanced.d1.only$mip_name_group<-unlist(stringr::str_extract_all(mip.ac.rebalanced.d1.comb$mip_name, "\\([^()]+\\)"))
      mip.ac.rebalanced.d1.only$Unique_IDs<-paste0(mip.ac.rebalanced.d1.only$TCG_ID, ":",
                                                   mip.ac.rebalanced.d1.only$PLATE_ID)
      mip.ac.rebalanced.d1.only<-mip.ac.rebalanced.d1.only[order(mip.ac.rebalanced.d1.only$TCG_ID),]
      mip.ac.rebalanced.d1.only$MIP_fraction<-with(mip.ac.rebalanced.d1.only, 
                                                   paste0("'", value, "/", Depth))
      
      mip.ac.rebalanced.d1.only$VAF<-as.numeric(mip.ac.rebalanced.d1.only$value)/as.numeric(mip.ac.rebalanced.d1.only$Depth)

    
    
# 9. Calculate Average VAF
    setDT(mip.ac.rebalanced.d1.only)[, Mode_VAF:= Mode(VAF, 0.001), by=Sample.TCG_ID][]
    setDT(mip.ac.rebalanced.d1.only)[, Median_VAF:= median(VAF[which(VAF>0.001)]), by=Sample.TCG_ID][]

  
  
# 10. Remove Empty Samples and Only Select Tumors
    mip.ac.rebalanced.d1.only<-mip.ac.rebalanced.d1.only[which(mip.ac.rebalanced.d1.only$Sample.TCG_ID!="BLANK"),]
    mip.ac.rebalanced.d1.only<-mip.ac.rebalanced.d1.only[grep("T", mip.ac.rebalanced.d1.only$Sample.TCG_ID),]
    mip.ac.rebalanced.d1.only$Clade<-haplos$Clade[match(mip.ac.rebalanced.d1.only$Sample.TCG_ID, haplos$Sample.TCG_ID)]

  
    
# 11. Generate MIPS Table
    mip.ac.rebalanced.d1.only.reshaped<-dcast(unique(mip.ac.rebalanced.d1.only[,c("mip_name", "Sample.TCG_ID", "Clade",
                                                                                "Mode_VAF", "Median_VAF", "MIP_fraction")]),
                                            Sample.TCG_ID  + Mode_VAF + Median_VAF+Clade  ~ mip_name, value.var = "MIP_fraction" )
    mip.ac.rebalanced.d1.only.reshaped<-t(mip.ac.rebalanced.d1.only.reshaped)


    
# 12. Identify Allele Count of MIPS that do not support the Sample Group
    tumors.unique<-unique(haplos$Sample.TCG_ID)
    only.unsupporting.vars<-NULL
    for(y in tumors.unique){
      
      temp<-mip.ac.rebalanced.d1.comb[which(mip.ac.rebalanced.d1.comb$Sample.TCG_ID %in% y),]
      temp$value[grep(paste0(temp$Clade[1], "|All_groups|VEP"), temp$mip_name)]<-NA
      
      temp$value[grep(paste0(temp$Clade[1], "|All_groups|VEP"), temp$mip_name)]<-NA
      temp$value[which(temp$value==0)]<-NA
      
      
      only.unsupporting.vars<-rbind(only.unsupporting.vars, temp)
      
    }
    
    only.unsupporting.vars$sample.ids<-paste0(only.unsupporting.vars$Sample.TCG_ID)
    only.unsupporting.vars.reshape<-dcast(unique(only.unsupporting.vars[,c("mip_name", "Sample.TCG_ID","Groups", "Clade", "value")]), 
                                          Sample.TCG_ID + Groups + Clade~mip_name, value.var = "value")
    only.unsupporting.vars.reshape<-as.data.frame(t(only.unsupporting.vars.reshape))
  
    write.csv(only.unsupporting.vars.reshape, file = "only.unsupporting.vars.reshape.tumors.duplicates.csv")
  



# 13. Reformat Data for Power Analysis
    mip.ac.rebalanced.d1.comb$id<-with(mip.ac.rebalanced.d1.comb, paste0(Sample.TCG_ID))
    mip.ac.rebalanced.tumors.only<-mip.ac.rebalanced.d1.only
    mip.ac.rebalanced.tumors.only$TCG_ID<-mip.ac.rebalanced.tumors.only$Sample.TCG_ID
    mip.ac.rebalanced.tumors.only$pseudoid<-paste0(mip.ac.rebalanced.tumors.only$Sample.TCG_ID)
    tumors.psuedo.id<-unique(mip.ac.rebalanced.tumors.only$pseudoid[grep("T", mip.ac.rebalanced.tumors.only$Sample.TCG_ID)])
    
    only.unsupporting.vars.reshape<-as.data.frame(only.unsupporting.vars.reshape)
    pseudo.id.for.temp.max<-paste0(as.character(unlist(only.unsupporting.vars.reshape[1,])), "_", as.character(unlist(only.unsupporting.vars.reshape[2,])))
    mip.ac.rebalanced.d1.comb$VAF<-mip.ac.rebalanced.d1.comb$value/mip.ac.rebalanced.d1.comb$Depth


# 14. Power Analysis and Contamination Calculation for Tumors
  pwr.tumors<-NULL
  for(y in unique(mip.ac.rebalanced.tumors.only$pseudoid)){
    
    temp<-mip.ac.rebalanced.tumors.only[which(mip.ac.rebalanced.tumors.only$Sample.TCG_ID %in% gsub("__1", "", y)),]
    temp.val<-temp$value[match(mips_rebalanced_vep$mip_name, temp$mip_name)]
    temp.Depth<-temp$Depth[match(mips_rebalanced_vep$mip_name, temp$mip_name)]
    if(length(which(as.character(unlist(only.unsupporting.vars.reshape[1,]))==y))>0){
      temp.max<-max(as.numeric(as.character(only.unsupporting.vars.reshape[c(4:(which(rownames(only.unsupporting.vars.reshape) == "Chrx_supercontig_000000040:210607_(DFT1_0_Buckland_exclusive)_0597")-1),
                                                                             (which(rownames(only.unsupporting.vars.reshape) == "Chrx_supercontig_000000040:210607_(DFT1_0_Buckland_exclusive)_0597")+1):nrow(only.unsupporting.vars.reshape)),
                                                                           which(as.character(unlist(only.unsupporting.vars.reshape[1,]))==y)])),na.rm=TRUE)}else{temp.max<-3}
    
    temp.cont.vars<-temp$VAF[which(temp$value>temp.max & temp$Depth>=10 & temp$VAF>0)]
    temp.cont.orig<-as.numeric(as.character(median(temp$VAF[which(temp$value>temp.max & temp$Depth>=10 & temp$VAF>0)], na.rm = TRUE)))*2
    
    if(length(temp.cont.vars)>=10){
      temp.cont<-temp.cont.orig
    }else{temp.cont<-NA
    temp.cont.orig<-as.numeric(as.character(median(temp$VAF[which(temp$Depth>=10)], na.rm = TRUE)))*2
    }
    
    read.depth<-median(temp.Depth, na.rm = TRUE)
    
    if(length(which(temp$VAF>0))==0 & length(temp.cont.vars)>=10){temp.cont<-0}
    
    #temp<-merge(temp, mips_rebalanced_vep[,c("mip_name")], all.y = TRUE)
    read.thres<-ifelse(temp.max %in% c(0,1,2), 5, ifelse(temp.max==3, 6, ifelse(temp.max==4, 7, ifelse(temp.max==5, 8, ifelse(temp.max ==6, 9, ifelse(temp.max >=7 & temp.max<=47, temp.max+3, 50))))))
    #temp.cont<-as.numeric(mip.ac.rebalanced.tumors.only$Median_VAF[match(y, mip.ac.rebalanced.tumors.only$pseudoid)])*2
    
    temp.dat<-c(temp$Sample.TCG_ID[1], temp.cont, temp.max, read.thres, read.depth, (temp.val - as.numeric(temp.max))/(temp.Depth * temp.cont.orig *0.5))
    pwr.tumors<-cbind(pwr.tumors, temp.dat)
    
  }
  
  colnames(pwr.tumors)<-pwr.tumors[1,]
  rownames(pwr.tumors)<-c("TCG_ID", "Median_VAF", "Max_contam_reads", "Read_threshold", "Median Read Depth", mips_rebalanced_vep$mip_name)
  
  supp.table1<-t(pwr.tumors[1:5,])
  pwr.tumors<-as.data.frame(pwr.tumors)

  tumors_metadata<-as.data.frame(supp.table1)


  # Filter out contaminating variants (per sample):
  mip.ac.rebalanced.d1.comb$unique.ids<-with(mip.ac.rebalanced.d1.comb, paste0(Sample.TCG_ID))
  mip.ac.rebalanced.d1.only.reshaped<-as.data.frame(t(mip.ac.rebalanced.d1.only.reshaped))

  mip.depth.pres.abs<-data.frame(mip_name = as.character(colnames(mip.ac.rebalanced.d1.only.reshaped))[4:ncol(mip.ac.rebalanced.d1.only.reshaped)])
  mip.depth.pres.abs$REF<-mips_rebalanced_vep$REF[match(mip.depth.pres.abs$mip_name, mips_rebalanced_vep$mip_name)]
  mip.depth.pres.abs$ALT<-mips_rebalanced_vep$ALT[match(mip.depth.pres.abs$mip_name, mips_rebalanced_vep$mip_name)]
  
  tumors_metadata<-tumors_metadata[grep("T", tumors_metadata$TCG_ID),]
  tumors_metadata<-unique(tumors_metadata)
  tumor_list<-tumors_metadata$TCG_ID
  tumors_metadata$Read_threshold<-as.numeric(as.character(tumors_metadata$Read_threshold))


  for(u in 1:nrow(tumors_metadata)){
    x_md<-tumors_metadata[u,]
    x = as.character(x_md$TCG_ID[1])
    
    x.mips<-rep("N", nrow(mip.depth.pres.abs))
    
    # First step: identify those with expected reads > 10
    reads<-NA
    reads<-mip.ac.rebalanced.d1.only.reshaped[which(mip.ac.rebalanced.d1.only.reshaped$Sample.TCG_ID==x),4:ncol(mip.ac.rebalanced.d1.only.reshaped)]
    
    x_rd_depth<-as.numeric(stringr::str_split_fixed(unlist(reads), "[/]", 2)[,2])
    x_rc<-as.numeric(gsub("'", "", stringr::str_split_fixed(unlist(reads), "[/]", 2)[,1]))
    exp_reads<-x_rd_depth*as.numeric(as.character(x_md$Median_VAF[1]))*0.5
    
    
    pab.index<-NA
    pab.index<-which(exp_reads>=10 & x_rc>=as.numeric(as.character(x_md$Read_threshold[1])))
    
    x.mips[pab.index]<-mip.depth.pres.abs$ALT[pab.index]
    
    
    
    # Second step: identify number of reads >2 but less than the read_threshold and call them ambiguous
    amb.index<-NA
    amb.index<-which(x_rc>2 & x_rc<x_md$Read_threshold[1])
    
    x.mips[amb.index]<-"N"
    
    
    # Third step: identify those that were called present, but with VAF<.01 # Use varying threshold
    x_vaf<-x_rc/x_rd_depth
    vaf.index<-which(x_vaf<.01 & exp_reads>=10 & x_rc>=as.numeric(as.character(x_md$Read_threshold[1])))
    x.mips[vaf.index]<-"N"
    
    
    # Third step: identify those absent if they have 0, 1, 2 reads
    abs.index<-NA
    abs.index<-which(x_rc <= 2 & exp_reads>=10)
    
    x.mips[abs.index]<-mip.depth.pres.abs$REF[abs.index]
    
    # Fourth step: expected reads ambiguous
    abs.e.index<-NA
    abs.e.index<-which(exp_reads<5)
    
    x.mips[abs.e.index]<-"N"
    
    # Fifth step: ambiguous
    amb.e.index<-NA
    amb.e.index<-which(exp_reads>=5 & exp_reads<10 & x_rc<x_md$Read_threshold[1])
    
    x.mips[amb.e.index]<-"N"
    
    # Sixth step: Present
    p.e.index<-NA
    p.e.index<-which(exp_reads>=5 & exp_reads<10 & x_rc>=x_md$Read_threshold[1])
    
    x.mips[p.e.index]<-mip.depth.pres.abs$ALT[p.e.index]  
    
    
    # Place in order
    mip.depth.pres.abs<-cbind(mip.depth.pres.abs, x.mips)
    colnames(mip.depth.pres.abs)[ncol(mip.depth.pres.abs)]<-as.character(x)
  }

  
# 15. Calculate the number of variants that cannot be assigned (per sample)
  mip.counts<-as.data.frame(apply(mip.depth.pres.abs[,4:ncol(mip.depth.pres.abs)], 2, 
                                function(x){c(length(which(x=="A")), length(which(x=="C")), 
                                              length(which(x=="T")), length(which(x=="G")),
                                              length(which(x=="N")), (length(which(x=="N"))/552))}))

# 16. Remove all samples with greater than 50% missing data
  samp.to.remove<-which(mip.counts[6,] > .5)
  samp.to.remove.table<-tumors_metadata[which(mip.counts[6,] > .6),]
  samp.to.remove.table$TCG_ID<-rownames(supp.table1)[which(mip.counts[6,] > .5)]

  supp.table1<-cbind(supp.table1[grep("T", rownames(supp.table1)),], t(mip.counts[6,]))
  colnames(supp.table1)[6]<-"N_prop"
  
  
# 17. Calculate the representation (missing fraction) per probe
  probes.to.remove<-apply(mip.depth.pres.abs[,-c(1:3, samp.to.remove+3)], 1, function(x){length(which(x=="N"))/(length(tumor_list))})
  probes.N.count<-data.frame(mip_name = mip.depth.pres.abs$mip_name,
                           N_prop = probes.to.remove)
  
# 18. Identify probes that have greater than 70% missing
  probes.lowN<-probes.N.count$mip_name[which(probes.N.count$N_prop>.7)]
  probes.lowN<-c(as.character(probes.lowN), host.probes)

  
# 19. Identify Manually  Germline mutations through the QC
  host.probes<-c( "Chr5_supercontig_000000478:61461_(All_groups_ancestral)_0535",
                  "Chr1_supercontig_000000015:1677191_(All_groups_excluding_377T1)_0004")

  
# 20. Adding Mitochondria info to to the presence/absence table
  MT_variants$mip_name<-with(MT_variants, paste0(POS, REF, ">", ALT))
  MT_pres_abs<-matrix(nrow = nrow(MT_variants), ncol = ncol(mip.depth.pres.abs)-3)
  MT_pres_abs<-cbind(MT_variants[,c("mip_name", "REF", "ALT")], MT_pres_abs)
  colnames(MT_pres_abs)<-colnames(mip.depth.pres.abs)

  for(z in 4:ncol(MT_pres_abs)){
  x = colnames(MT_pres_abs)[z]
  mt_samp<-NULL
  if(length(which(haplos$Sample.TCG_ID %in% x))>0){
    for(y in 1:nrow(MT_variants)){
      mt.var.samp<-ifelse(grepl(MT_pres_abs$mip_name[y], haplos$Defining.variants[match(x, haplos$Sample.TCG_ID)]),
                          as.character(MT_pres_abs$ALT[y]), as.character(MT_pres_abs$REF[y]))
      mt_samp<-c(mt_samp, mt.var.samp)
    }}else{
      mt_samp<-rep("N", nrow(MT_variants))
    }
  MT_pres_abs[,z]<-mt_samp
  
}


# 21. Merge the MIPS and Mitochondrial data tables:
  p.a.tot<-rbind(mip.depth.pres.abs,MT_pres_abs)
  p.a.tot[,ncol(p.a.tot)]<-p.a.tot$REF
  colnames(p.a.tot)[ncol(p.a.tot)]<-"reference"
  p.a.tot1<-p.a.tot[,-c(samp.to.remove+3)]
  samp.to.remove<-colnames(p.a.tot)[c(samp.to.remove+3)]
  
  
# 22. Add 0 (sample carry ref only) and 1 (sample carries at least one alt allele) and NA (undetermined) to create a presence-absence table
  p.a.table<-t(apply(p.a.tot, 1, function(x){
    idx<-which(as.character(x[4:length(x)])==x[2])+3
    a.idx<-which(as.character(x[4:length(x)])==x[3])+3
    x[idx]<-0
    x[a.idx]<-1
    
    return(x)
    
  }))
  
  
# 23. Reformat table to include haplogroup information as a Supplementary Table
  p.a.tot_pres.abs<-rbind(p.a.table, haplos$Groups[match(colnames(p.a.table), haplos$Sample.TCG_ID)])
  p.a.tot_pres.abs<-rbind(p.a.tot_pres.abs, haplos$Complete.Haplotype[match(colnames(p.a.tot_pres.abs), haplos$Sample.TCG_ID)])
  p.a.tot_pres.abs<-rbind(p.a.tot_pres.abs, haplos$Clade[match(colnames(p.a.tot_pres.abs), haplos$Sample.TCG_ID)])
  p.a.tot_pres.abs<-p.a.tot_pres.abs[-1,]
  rownames(p.a.tot_pres.abs)[626:628]<-c("Group", "MT Haplogroup", "Clade")

  p.a.tot_pres.abs<-p.a.tot_pres.abs[,c(1:3,order(p.a.tot_pres.abs[626,4:ncol(p.a.tot_pres.abs)])+3)]
  p.a.tot_pres.abs<-p.a.tot_pres.abs[,c(1:3,order(p.a.tot_pres.abs[627,4:ncol(p.a.tot_pres.abs)])+3)]
  p.a.tot_pres.abs<-p.a.tot_pres.abs[,c(1:3,order(p.a.tot_pres.abs[628,4:ncol(p.a.tot_pres.abs)])+3)]
  
  p.a.tot_pres.abs<-p.a.tot_pres.abs[c(626:628, 1:625),]

  p.a.tot_pres.abs.meta<-p.a.tot_pres.abs[1:3,]
  p.a.tot_pres.abs.data<-p.a.tot_pres.abs[4:nrow(p.a.tot_pres.abs),]
  MT<-p.a.tot_pres.abs[-grep("Chr", p.a.tot_pres.abs[,1]),]

  p.a.tot_pres.abs.data<-p.a.tot_pres.abs.data[order(unlist(str_match_all(p.a.tot_pres.abs.data[,1], "(?<=\\().+?(?=\\))"))),]
  p.a.tot_pres.abs.data<-rbind(p.a.tot_pres.abs.meta, p.a.tot_pres.abs.data, MT[-c(1:3),])
  p.a.tot_pres.abs.data<-p.a.tot_pres.abs.data[,c(1,2,3, order(p.a.tot_pres.abs.data[2,4:ncol(p.a.tot_pres.abs.data)])+3)]
  p.a.tot_pres.abs.data<-p.a.tot_pres.abs.data[,c(1,2,3, order(p.a.tot_pres.abs.data[1,4:ncol(p.a.tot_pres.abs.data)])+3)]
  p.a.tot_pres.abs.data<-p.a.tot_pres.abs.data[,c(1,2,3, order(p.a.tot_pres.abs.data[3,4:ncol(p.a.tot_pres.abs.data)])+3)]

  
# 24. Filter table removing samples and probes
  filtered.table<-p.a.tot_pres.abs.data[-which(p.a.tot_pres.abs.data[,1] %in% probes.lowN), -which(colnames(p.a.tot_pres.abs.data) %in% samp.to.remove)]
  removed.table<-p.a.tot_pres.abs.data[c(1:3, which(p.a.tot_pres.abs.data[,1] %in% probes.lowN)), c(1:3, which(colnames(p.a.tot_pres.abs.data) %in% samp.to.remove))]

  
# 25. Recalculate missing fraction in the probes and also in the samples, and remove any probes that is represented in less than 50% 
  probes.lowN.added<-apply(filtered.table[-c(1:3),4:ncol(filtered.table)], 1, function(x){
    
    y<-length(which(x=="N"))/length(x)
    return(y)
  })
  names(probes.lowN.added)<-filtered.table[-c(1:3),1]
  probes.lowN.added<-names(probes.lowN.added)[which(probes.lowN.added>.5)]

  
# 26. Recalculate missing fraction in the probes and also in the samples, and remove any probes that are called as reference across all filtered samples
  probes.nopres<-apply(filtered.table[-c(1:3),4:ncol(filtered.table)], 1, function(x){
    
    y<-length(which(x=="1"))/length(x)
    return(y)
  })
  names(probes.nopres)<-filtered.table[-c(1:3),1]
  probes.nopres<-names(probes.nopres)[which(probes.nopres==0)]
  
  filtered.table.postfilter<-filtered.table[-which(filtered.table[,1] %in% c(probes.lowN.added, probes.nopres)),]



# 27. Remove Probes with >50% N within the target group (group that sample is classified in) using group-representing probes
  probe.names<-filtered.table.postfilter[-c(1:3),1]
  probe.names<-unlist(str_match_all(probe.names, "(?<=\\().+?(?=\\))"))
  probe.names<-probe.names[-which(probe.names %in% c("All_groups_excluding_377T1", "All_groups_ancestral",  "VEP_ATR", "VEP_MET", "VEP_MSH6"))]

  groups<-t(filtered.table.postfilter)
  colnames(groups)[4:ncol(groups)]<-groups[1,4:ncol(groups)]
  groups<-groups[-1,]
  groups<-as.data.frame(groups)
  groups<-groups[-c(1:2),]
  groups$Group<-haplos$Groups[match(rownames(groups), haplos$Sample.TCG_ID)]
  groups$`MT Haplogroup`<-haplos$Complete.Haplotype[match(rownames(groups), haplos$Sample.TCG_ID)]
  

  # Updates to sample classification:
  groups$`MT Haplogroup`[grep("358T3.1", rownames(groups))]<-"DFT1_0"
  groups$Group[which(rownames(groups)=="426T1")]<-"DFT1_3_Channel"
  groups$Group[which(rownames(groups)=="390T1")]<-"DFT1_0_Freycinet"
  groups$Group[which(rownames(groups)=="219T")]<-"DFT1_0_Freycinet"
  
  haplos$Groups[which(haplos$Sample.TCG_ID=="426T1")]<-"DFT1_3_Channel"
  haplos$Groups[which(haplos$Sample.TCG_ID=="390T1")]<-"DFT1_0_Freycinet"
  haplos$Groups[which(haplos$Sample.TCG_ID=="219T")]<-"DFT1_0_Freycinet"
  haplos$Groups[which(haplos$Sample.TCG_ID=="1049T1")]<-"DFT1_2_Kempton"



  groups<-groups[-which(rownames(groups)=="reference"),]
  group.calls<-NULL
  for(y in 1:3){
    for(x in 4:ncol(groups)){
      temp<-table(groups[,c(y, x)]) 
      temp<-melt(temp)
      temp$mip_name<-names(temp)[2]
      names(temp)<-c("Grouping", "Call", "Value", "mip_name")
      group.calls<-rbind(group.calls, temp)
    }
  }
  
  
# 28. Import MIPs grouping and Haplotype info
  
  mip_rebalanced_group_match <- read.csv("~/Downloads//mip_rebalanced_group_match1.csv", stringsAsFactors=FALSE)
  
  mip_rebalanced_group_match$Groups[which(mip_rebalanced_group_match$Groups=="")]<-mip_rebalanced_group_match$Clades[which(mip_rebalanced_group_match$Groups=="")]
  
  mip_rebalanced_group_match.filtered<-mip_rebalanced_group_match[which(mip_rebalanced_group_match$mip_name %in% filtered.table.postfilter[,1]),]
  
  group.mips.match<-NULL
  for(y in 1:nrow(mip_rebalanced_group_match.filtered)){
    
    list.groups<-gsub(" ", "", unlist(str_split(mip_rebalanced_group_match.filtered$Targeted.Groups[y], ";")))
    filtered.temp<-group.calls[which(gsub(" ", "", group.calls$Grouping) %in% list.groups & group.calls$mip_name == mip_rebalanced_group_match.filtered$mip_name[y]),]
    filtered.temp<-filtered.temp[which(filtered.temp[,2] %in% c("0", "1", "N")),]
    N.prop<-sum(filtered.temp$Value[which(filtered.temp$Call=="N")])/sum(filtered.temp$Value)
    tot.counts<-sum(filtered.temp$Value)
    calls<-c(filtered.temp$mip_name[1],N.prop, tot.counts)
    group.mips.match<-rbind(group.mips.match, calls)
  }
  
  group.mips.match<-group.mips.match[-which(is.na(group.mips.match[,1])),]
  group.mips.match<-as.data.frame(group.mips.match)
  group.mips.match$Probe_set<-unlist(stringr::str_extract_all(group.mips.match[,1], "\\([^()]+\\)"))
  group.mips.match.table<-unique(group.mips.match[,c("Probe_set", "V3")])


# 29. Identify group-representing probes which have greater than 50% missingness across samples that are in the group
  probes.lowN.added2<-group.mips.match[which(as.numeric(as.character(group.mips.match[,2]))>.5),1]
  filtered.table.postfilter2<-filtered.table.postfilter[-which(filtered.table.postfilter[,1] %in% probes.lowN.added2),]
  
  
# 30. Remove Probes with all 0s/Ns, 0s, or Ns
  f.z<-apply(filtered.table.postfilter2[-c(1:3),4:ncol(filtered.table.postfilter2)], 1, function(x){
    y<-length(which(x %in% c("0", "N")))/length(x)
    return(y)
    
  })
  f.z<-filtered.table.postfilter[-c(1:3),1][which(f.z==1)]
  
  f.N<-apply(filtered.table.postfilter2[-c(1:3),4:ncol(filtered.table.postfilter2)], 1, function(x){
    y<-length(which(x=="N"))/length(x)
    return(y)
    
  })
  f.N<-filtered.table.postfilter[-c(1:3),1][which(f.N==1)]

  write.csv(filtered.table.postfilter2, file = "filtered.table.postfilter2.duplicates.50.csv")



# 31. Identify probes and samples that remain from the above filtering
  tot.probes.removed<-mips_rebalanced_vep$mip_name[-which(mips_rebalanced_vep$mip_name %in% filtered.table.postfilter2[,1])]
  samp.removed<-supp.table1[-which(supp.table1[,1] %in% colnames(filtered.table.postfilter2)),]


  filtered.table.postfilter3<-filtered.table.postfilter2

  mip.ac.rebalanced.d1.only.reshaped<-dcast(unique(mip.ac.rebalanced.d1.only[,c("mip_name", "Sample.TCG_ID", 
                                                                                "Mode_VAF", "Median_VAF", "MIP_fraction")]),
                                            Sample.TCG_ID  + Mode_VAF + Median_VAF  ~ mip_name, value.var = "MIP_fraction" )

  mip.ac.rebalanced.d1.only.reshaped$Median_VAF<-as.character(mip.ac.rebalanced.d1.only.reshaped$Median_VAF)
  mip.ac.rebalanced.d1.only.reshaped$Median_VAF<-supp.table1[match(unlist(mip.ac.rebalanced.d1.only.reshaped$Sample.TCG_ID), supp.table1[,1]), 2]


# 32. Calculate the number of expected reads to call presence for variants representing groups.
      # If number of reads is lower than number of expected reads, or there is less than 30 total reads across the variant, then assign "N" for ambiguity

  exp.reads.df<-NULL
  for(x in colnames(filtered.table.postfilter3)[-c(1:3, ncol(filtered.table.postfilter3))]){
    temp<- unlist(mip.ac.rebalanced.d1.only.reshaped[which(mip.ac.rebalanced.d1.only.reshaped$Sample.TCG_ID %in% x),])
    temp<-as.character(temp)
    temp.reads<- as.numeric(str_split_fixed(temp[4:length(temp)], "/", 2)[,2])*as.numeric(as.character(temp[3]))*0.5
    names(temp.reads)<-colnames(mip.ac.rebalanced.d1.only.reshaped)[4:ncol(mip.ac.rebalanced.d1.only.reshaped)]
    temp.probes<-names(temp.reads)[which(temp.reads<30)]
    temp.Ns<-filtered.table.postfilter3[which(filtered.table.postfilter3[,1] %in% temp.probes),x]
    filtered.table.postfilter3[which(filtered.table.postfilter3[,1] %in% temp.probes),x][which(temp.Ns=="0")]<-"N"
    
    temp.complete<-c(as.character(temp[1:3]), temp.reads)
  
    exp.reads.df<-rbind(exp.reads.df, temp.complete)
  }


# 33. Remove any probes that are either all 0/Ns or Ns
  f.z<-apply(filtered.table.postfilter3[-c(1:3),4:ncol(filtered.table.postfilter3)], 1, function(x){
    y<-length(which(x %in% c("0", "N")))/length(x)
    return(y)
    
  })
  f.z<-filtered.table.postfilter3[-c(1:3),1][which(f.z==1)]
  
  f.N<-apply(filtered.table.postfilter3[-c(1:3),4:ncol(filtered.table.postfilter3)], 1, function(x){
    y<-length(which(x=="N"))/length(x)
    return(y)
    
  })
  f.N<-filtered.table.postfilter3[-c(1:3),1][which(f.N==1)]


# 34. Remove probes and samples with Greater than 50% Ns
  f.N<-apply(filtered.table.postfilter3[-c(1:3),4:ncol(filtered.table.postfilter3)], 1, function(x){
    y<-length(which(x=="N"))/length(x)
    return(y)
    
  })
  f.N<-filtered.table.postfilter3[-c(1:3),1][which(f.N>.5)]
  filtered.table.postfilter3<-filtered.table.postfilter3[-which(filtered.table.postfilter3[,1] %in% f.N),]


  f.s<-apply(filtered.table.postfilter3[-c(1:3),4:ncol(filtered.table.postfilter3)], 2, function(x){
    y<-length(which(x %in% c("N")))/length(x)
    return(y)
    
  })
  f.s<-colnames(filtered.table.postfilter3[-c(1:3),4:ncol(filtered.table.postfilter3)])[which(f.s>.5)]
  filtered.table.postfilter3<-filtered.table.postfilter3[,-which(colnames(filtered.table.postfilter3) %in% f.s)]



# 35. Reformat table to create the fasta file
  dt.mat.samples.orig<-filtered.table.postfilter3[4:nrow(filtered.table.postfilter3),-c(1:3)]
  rownames(dt.mat.samples.orig)<-rownames(filtered.table.postfilter3)[4:nrow(filtered.table.postfilter3)]



  d.knn.data.tree.seq.orig<-apply(dt.mat.samples.orig, 2, 
                                  function(x){
                                    y<-x
                                    y[which(x==0)]<-filtered.table.postfilter2.dist[(which(x==0)+3),2]
                                    y[which(x>0)]<-filtered.table.postfilter2.dist[(which(x>0)+3),3]
                                    y[which(is.na(x))]<-ifelse(grepl("GT|TG", paste0(filtered.table.postfilter2.dist[which(is.na(x)),2],filtered.table.postfilter2.dist[which(is.na(x)),3])), "K", ifelse(
                                      grepl("AC|CA", paste0(filtered.table.postfilter2.dist[which(is.na(x)),2],filtered.table.postfilter2.dist[which(is.na(x)),3])), "M", ifelse(
                                        grepl("CG|GC", paste0(filtered.table.postfilter2.dist[which(is.na(x)),2],filtered.table.postfilter2.dist[which(is.na(x)),3])), "S", ifelse( 
                                          grepl("AT|TA", paste0(filtered.table.postfilter2.dist[which(is.na(x)),2],filtered.table.postfilter2.dist[which(is.na(x)),3])), "W", ifelse( 
                                            grepl("AG|GA", paste0(filtered.table.postfilter2.dist[which(is.na(x)),2],filtered.table.postfilter2.dist[which(is.na(x)),3])), "R", ifelse( 
                                              grepl("CT|TC", paste0(filtered.table.postfilter2.dist[which(is.na(x)),2],filtered.table.postfilter2.dist[which(is.na(x)),3])), "Y", "N"))))))
                                    return(y)
                                    
                                  })
  
  
  fas.seq<-apply(d.knn.data.tree.seq.orig, 2, function(x){
    paste(x, collapse = "")
  })
  names(fas.seq)<-paste0(names(fas.seq), "_", seq(1, length(fas.seq), 1))
  fas.seq<-DNAStringSet(fas.seq)

