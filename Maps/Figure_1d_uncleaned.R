
library(plyr)
library(dplyr)
library(ggplot2)

load("~/6_4_18_haplos.RData")

# DotPlot: 
Combined<-haplos
Combined<-Combined[which(Combined$Duplicate.or.Time.Course.biopsy.==""),]
Combined<-Combined[duplicated(Combined[,c("TCG_ID")]) | duplicated(Combined[,c("TCG_ID")], fromLast=TRUE),]

Combined$IDs<-with(Combined, paste0(TCG_ID, ":", Complete.Haplotype))
Combined$col<-""
cols<-brewer.pal(3,"Set3")
for(y in unique(Combined$TCG_ID)){
  test<-Combined[which(Combined$TCG_ID == y),]
  test.count<-plyr::count(test$Complete.Haplotype)
  if(nrow(test.count)==1){
    Combined$col[which(Combined$TCG_ID == y)]<-cols[1]
  }
  if(nrow(test.count)>1){
    test.count$col<-cols[1:nrow(test.count)]
    test.count$TCG_ID<-unique(test$TCG_ID)
    test.count$Complete.Haplotype<-test.count$x
    test.count<-test.count[order(test.count$freq),]
    test.count<-test.count[,-c(1:2)]
    test.count$IDs<-with(test.count, paste0(TCG_ID, ":", Complete.Haplotype))
    for(i in 1:nrow(test.count)){
      Combined$col[which(Combined$IDs == test.count$IDs[i] & 
                           Combined$TCG_ID==y)]<-test.count$col[i]
    }
  }
}

# determine dot positions

x<-data %>% group_by(TCG_ID, col) %>% summarize(n=n()) %>% arrange(desc(n))
xt<-data %>% group_by(TCG_ID) %>% filter(n()>1) %>% summarize(total=n()) %>% arrange(desc(total))
x<-merge(x, xt)

xc.plot<-matrix(ncol = 5, nrow = 0, dimnames = list(NULL, c(colnames(x), "y")))
for(z in 1:nrow(x)){
  temp<-data.frame(x[z,] %>% slice(rep(row_number(), x$n[z])))
  if(temp$col[1]=="#8DD3C7"){
    temp$y<-seq(1, nrow(temp), 1)
  }
  if(temp$col[1]=="#FFFFB3"){
    temp$y<-seq(temp$total[1]-temp$n[1]+1, temp$total[1], 1)
  }
  
  xc.plot<-rbind(xc.plot, temp)
}

xc.plot<-xc.plot[(order(xc.plot$n, decreasing = TRUE)),]
xc.plot<-xc.plot[(order(xc.plot$col)),]
xc.plot<-xc.plot[order(xc.plot$total),]

xc.id<-data.frame(unique(xc.plot$TCG_ID))
colnames(xc.id)<-"TCG_ID"
xc.id$id<-seq(1, nrow(xc.id), 1)
xc.plot<-merge(xc.plot,xc.id)
xc.plot$col<-as.numeric(as.factor(xc.plot$col))
xc.plot$col[which(xc.plot$col==1)]<-"Predominant DFTD Strain"
xc.plot$col[which(xc.plot$col==2)]<-"Secondary DFTD Strain"


p<-ggplot(xc.plot, aes(id,y, color = col))+
  geom_line(aes(id, y, group = id), lwd = 1,  color = "#48403D", alpha = .3)+
  geom_point(size = 5)+
  scale_color_manual(values = c("black", "red"))+
  labs(x = "Individual Devils", y = "Number of Tumors", main = "Different Tumor Types in Devils")+
  scale_y_continuous(breaks = seq(0,8,1))+
  scale_x_continuous(expand = c(0.01, 0))+
  theme(axis.line = element_line(colour = "black"),
        text=element_text(size =15, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        axis.title.y = element_text(vjust = 1.1, size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

bar<-matrix(ncol=2,nrow = 2, dimnames = list(NULL, c("Infection", "Devils")))
bar[,1]<-c("1 DFT1 Infection", ">1 DFT1 Infection")
bar[,2]<-c(length(unique(x$TCG_ID))-length(unique(x$TCG_ID[which(x$col=="#FFFFB3")])), 
           length(unique(x$TCG_ID[which(x$col=="#FFFFB3")])))
bar<-as.data.frame(bar)
bar$Infection<-factor(bar$Infection, levels = bar$Infection[order(bar$Devils, decreasing = TRUE)])

p1<-ggplot(bar, aes(Infection, as.numeric(as.character(Devils)), fill = Infection))+geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#685A5A","#A73232"))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 100, 20), limits = c(0,90))+
  labs(x = "", y = "Number of Devils")+
  theme(legend.position = c(0.75, 0.85))+
  theme(axis.line = element_line(colour = "black"),
        text=element_text(size =15, colour = "black"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(vjust = 1, hjust = 1, size = 16),
        axis.title.y = element_text(vjust = 1.1, size = 18),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf("Figure1c_proposed.pdf", height = 3, width = 24)
print(p)
dev.off()

pdf("Figure1c.1_proposed.pdf", height = 8, width = 10)
print(p1)
dev.off()
