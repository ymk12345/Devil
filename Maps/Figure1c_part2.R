Combined<-haplos
tumors<-Combined[which(Combined$Duplicate.or.Time.Course.biopsy.==""),]

DFT1_tumors<-tumors[grep("DFT1", tumors$Haplogroup),]
DFT1_tumors<-DFT1_tumors[-which(duplicated(DFT1_tumors[,c("TCG_ID", "Complete.Haplotype")])),]
DFT1_tumors$Numeric_haplogroup<-gsub("DFT1_", "", DFT1_tumors$Haplogroup)
DFT1_tumors$Numeric_Complete.Haplotype<-gsub("DFT1_", "", DFT1_tumors$Complete.Haplotype)

DFT1_tumors$Haplogroup<-factor(DFT1_tumors$Haplogroup, levels = unique(DFT1_tumors$Haplogroup[order(as.numeric(DFT1_tumors$Numeric_haplogroup))]))
DFT1_tumors$Complete.Haplotype<-factor(DFT1_tumors$Complete.Haplotype, levels = unique(DFT1_tumors$Complete.Haplotype[order(as.numeric(DFT1_tumors$Numeric_haplogroup))]))
DFT1_tumors$Location[grep("West Pencil", DFT1_tumors$Location)]<-"West Pencil Pine"

# West Pencil Pine:

p1<-ggplot(DFT1_tumors[which(DFT1_tumors$Location=="West Pencil Pine"),], 
           aes(as.numeric(Year), fill = Haplogroup))+
  geom_bar(width = 0.9, color = "black")+
  theme_bw() + 
  scale_fill_manual(values = c("grey", "#315195", "#F1C33B", "#FDBE85", "#A5153E", "#4D1313"))+
  xlab("")+
  ylab("No. of Devils")+
  labs(fill = "Haplogroup")+
  ggtitle("West Pencil Pine")+
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(2002, 2018), breaks = seq(2002, 2018, 1))+
  theme(axis.line = element_line(colour = "black"),
        text=element_text(size =15, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1.2, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(vjust = 1.1, size = 20),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


# Forestier   
p2<-ggplot(DFT1_tumors[which(DFT1_tumors$Location=="Forestier"),], 
           aes(as.numeric(Year), fill = Complete.Haplotype))+
  geom_bar(width = 0.9, color = "black"
           #color = "black"
  )+
  theme_bw() + 
  scale_x_continuous(limits = c(2002, 2018), breaks = seq(2002, 2018, 1))+
  scale_fill_brewer(palette="PRGn")+
  #scale_fill_manual(values = rev(c("grey", "#315195", "#F1C33B", "#8A3CD5", "#A5153E")))+
  xlab("")+
  ylab("No. of Devils")+
  labs(fill = "Haplogroup")+
  ggtitle("Forestier")+
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.line = element_line(colour = "black"),
        text=element_text(size =15, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1.2, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(vjust = 1.1, size = 20),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# Bronte
p3<-ggplot(DFT1_tumors[which(DFT1_tumors$Location=="Bronte Park"),], 
           aes(as.numeric(Year), fill = Complete.Haplotype))+
  geom_bar(width = 0.9, color = "black")+
  theme_bw() + 
  #ggthemes::scale_fill_ptol()+
  scale_fill_manual(values = c("grey", "#315195", rev(brewer.pal(4, "Blues"))))+
  xlab("")+
  ylab("No. of Devils")+
  labs(fill = "Haplogroup")+
  ggtitle("Bronte Park")+
  scale_x_continuous(limits = c(2002, 2018), breaks = seq(2002, 2018, 1))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 10, 2)) +
  theme(axis.line = element_line(colour = "black"),
        text=element_text(size =15, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1.2, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(vjust = 1.1, size = 20),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


# Fentonbury

p4<-ggplot(DFT1_tumors[which(DFT1_tumors$Location=="Fentonbury"),], 
           aes(as.numeric(Year), fill = Haplogroup))+
  geom_bar(width = 0.9, color = "black")+
  theme_bw() + 
  #ggthemes::scale_fill_ptol()+
  scale_fill_manual(values = c("grey", "#315195", rev(brewer.pal(4, "PuOr")[1:2])))+
  xlab("")+
  ylab("No. of Devils")+
  labs(fill = "Haplogroup")+
  ggtitle("Fentonbury")+
  scale_x_continuous(limits = c(2002, 2018), breaks = seq(2002, 2018, 1))+
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.line = element_line(colour = "black"),
        text=element_text(size =15, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1.2, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(vjust = 1.1, size = 20),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

pdf("lineage_dynamics_part2_devils.pdf", height = 18, width = 20)
grid.arrange(p1, p2, p3, p4, ncol = 1)
dev.off()
