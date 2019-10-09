library(xlsx);
library(png);
library(ggplot2);
library(grid);
library(gridExtra);

give.n <- function(x){
  return(c(y = median(x)+.08, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


load("~/6_4_18_haplos.RData")
Combined<-haplos
Combined<-Combined[grep("Yes", Combined$Include.for.Spread.Analyses.),]
Combined<-Combined[-which(Combined$Clade %in% c("Unassigned", "", "Non-DFTD", "377T1", "DFT2")),]
Combined$Clade[grep("B", Combined$Clade)]<-"Clade_B"
Combined$Clade[grep("A1", Combined$Clade)]<-"Clade_A1"
Combined$Clade[grep("A2", Combined$Clade)]<-"Clade_A2"
Combined$Clade[grep("Clade C", Combined$Clade)]<-"Clade_C"



# Include the Mapping Data
img<-readPNG("~/Desktop/Tasmania.png", native = TRUE)
g<-rasterGrob(img, interpolate=TRUE)

unique.list<-unique(as.numeric(as.character(Combined$Latitude)))
unique.list<-unique.list[!is.na(unique.list)]
Mapping.table<-matrix(ncol = 8, nrow = 0, dimnames = list(NULL, c("Sample.TCG_ID","Haplotype.New","Complete.Haplotype", "Clade", "Location","Latitude", "Longitude", "Year")))
for (x in 1:length(unique.list)){
  no.sample.10<-trunc(nrow(Combined[which(Combined$Latitude == unique.list[x]),])/10);
  data1<-Combined[which(Combined$Latitude == unique.list[x]),]
  data1$Number<-NA
  data1$Map<-NA
  data<-data1[order(data1$Clade),]
  if(nrow(data)<=10){
    data$Number[1:nrow(data)]<- -1
    if(nrow(data)<=4){
      data$Map[1:nrow(data)]<-'dot'
    }
    if(nrow(data)>4){
      data$Map[1:nrow(data)]<-'box'
      
    }
  }
  if(nrow(data)>10){
    data$Map[1:nrow(data)]<-'box'
    for(i in 1:trunc(nrow(data)/10)){
      if(i == 1){
        data$Number[1:((i*10)+1)]<- -1
        if(i == max(trunc(nrow(data)/10)) & nrow(data)>(max(trunc(nrow(data)/10))*10)){
          data$Number[((max(trunc(nrow(data)/10))*10)+1):((max(trunc(nrow(data)/10))*10)+(nrow(data)-(max(trunc(nrow(data)/10))*10)))]<- -(max(trunc(nrow(data)/10))+1)
        }
      }
      
      else{
        if(i>1 & i< max(trunc(nrow(data)/10))){
          data$Number[(((i-1)*10)+1):(i*10)]<- -i}
        if(i == max(trunc(nrow(data)/10)) & nrow(data)==(max(trunc(nrow(data)/10))*10)){
          data$Number[(((i-1)*10)+1):nrow(data)]<- -i
        }
        if(i == max(trunc(nrow(data)/10)) & nrow(data)>(max(trunc(nrow(data)/10))*10)){
          data$Number[(((i-1)*10)+1):(i*10)]<- -i
          data$Number[((max(trunc(nrow(data)/10))*10)+1):((max(trunc(nrow(data)/10))*10)+(nrow(data)-(max(trunc(nrow(data)/10))*10)))]<- -(max(trunc(nrow(data)/10))+1)
        }
      }}
  }
  Mapping.table<-rbind(Mapping.table, data)
}


# Colors
#colorHaplotype<-c('#01A9DB', '#01BE82', '#DF7401', '#86B404', '#013ADF', '#7401DF', '#B43104', '#31B404', '#0174DF', '#FFBF00', '#299571', '#2AD5C1', '#B833AB', '#55BFF3', '#F355A7', '#1ADCC8', '#DF0101', '#04B45F', '#FD026F', '#F7AF13', '#A102FD', '#FD860E', '#16A726', "#2E86C1", "#A569BD", "#16A085", "#E67E22", "#6C3483", "#CB4335", "#CB4335", "#CB4335", "#CB4335")
colorHaplotype<-c("#F80B0B", "#F477A0","#68A9CC","#8AB948") # Clade A1, Clade A2, Clade B, Clade C
Haplogrp.list<-unique(Combined$Clade)
shapes<-c("a", "b", "c", "d")

colorHaplotype<-cbind(colorHaplotype, Haplogrp.list, shapes)
#--------------------------------------------------------------------------------------------#

pdf("broad_prelim_clades_per_year_maps_not_filled_dots_v1.pdf", height = 10, width = 10.5)
ggplot(data = Mapping.table, aes(x = as.numeric(as.character(Longitude)), y = as.numeric(as.character(Latitude)), color = Clade), size = 3)+
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_bw()+
  ylim(-43.676477, -40.5479)+xlim(144.45642, 148.48056)+
  scale_shape_discrete(solid=F) +
  facet_wrap(~Year, ncol = 4)+
  geom_jitter(aes(x = as.numeric(as.character(Longitude)), y = as.numeric(as.character(Latitude)), shape=Clade), size = 1.8, alpha = 1, width = .08, height = .08)+
  #geom_point(data = Dots.Final, aes(x=as.numeric(as.character(PseudoLong)), 
  #                                  y = as.numeric(as.character(PseudoLat))),
  #           colour = Dots$Color, size = 5.7, alpha = .5, show.legend = TRUE)+
  #scale_shape_manual(values = c(19,1))+
  scale_color_manual(values = colorHaplotype)+
  theme(
    #legend.position = "none",
    strip.background = element_rect(colour = "#B0B0B6"),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  labs(x = "", y = "")
dev.off()



# Per Clade:

for(y in Haplogrp.list){

#temp.shape<-colorHaplotype[which(colorHaplotype[,2]==y),3]
temp.color<-colorHaplotype[which(colorHaplotype[,2]==y),1]

  
pdf(paste0("total_", y, "_clade.pdf"), height = 10, width = 10)
print(ggplot(data = Mapping.table[which(Mapping.table$Clade==y),], aes(x = as.numeric(as.character(Longitude)), y = as.numeric(as.character(Latitude))), size = 3)+
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_bw()+
  ylim(-43.676477, -40.5479)+xlim(144.45642, 148.48056)+
  scale_shape_discrete(solid=F) +
  geom_jitter(aes(x = as.numeric(as.character(Longitude)), y = as.numeric(as.character(Latitude))), shape = 4, color = temp.color,  size = 6, alpha = 1, width = .08, height = .08)+
  scale_color_manual(values = colorHaplotype)+
  theme(
    legend.position = "none",
    strip.background = element_rect(colour = "#B0B0B6"),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  labs(x = "", y = "", title = y))
dev.off()
}





