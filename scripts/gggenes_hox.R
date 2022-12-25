#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gggenes))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(eoffice))


species="Res"
species="Cja"
HGT<-read.delim(paste0(species,"_hox.txt"), header = TRUE)

#check strand start stop
for (i in 1:nrow(HGT)) {
  if (HGT$strand[i] == "-") {
    if (HGT$start[i]<HGT$stop[i]) {
      a<-HGT$stop[i]
      HGT$stop[i]<-HGT$start[i]
      HGT$start[i]<-a
    } else {
      
    }
  } else {
    if (HGT$start[i]>HGT$stop[i]) {
      a<-HGT$stop[i]
      HGT$stop[i]<-HGT$start[i]
      HGT$start[i]<-a
    } else {
      
    }
  }
}

HGT<-HGT[order(HGT$chr,HGT$start),]
HGT$genelength<-abs(HGT$stop-HGT$start)
row.names(HGT)<-1:nrow(HGT)

##intergenic distance
intergenic<-as.data.frame(matrix(NA, nrow = nrow(HGT), ncol = 8))
names(intergenic)<-c("gene.u", "gene.d",
                     "chr", "chr2",
                     "start", "stop",
                     "genecount","distance")
for (i in 1:(nrow(HGT)-1)) {
  intergenic$gene.u[i]<-HGT$Gene[i]
  intergenic$gene.d[i]<-HGT$Gene[i+1]
  intergenic$chr[i]<-HGT$chr[i]
  intergenic$chr2[i]<-HGT$chr[i+1]
  intergenic$start[i]<-max(HGT$stop[i],HGT$start[i])
  intergenic$stop[i]<-min(HGT$start[i+1], HGT$stop[i+1])
}
intergenic<-intergenic[intergenic$chr==intergenic$chr2 & is.na(intergenic$chr)==F,]
intergenic$distance<-intergenic$stop-intergenic$start
intergenic<-intergenic[,c("gene.u", "gene.d",
                          "chr",
                          "start", "stop",
                          "genecount","distance")]
row.names(intergenic)<-1:nrow(intergenic)

gff<-read.delim(paste0(species,"_gene.gff3"), header = FALSE)
gff<-gff[!duplicated(gff$V9),]
gff<-gff[,c(1,4,5,7,9)]
gff$l<-gff$V5-gff$V4
names(gff)<-c("chr", "start", "stop", "strand", "gene", "genelength")

for (i in 1:nrow(intergenic)) {
  gff.sub<-subset(gff, chr == intergenic$chr[i] & start > intergenic$start[i] & stop < intergenic$stop[i])
  intergenic$genecount[i]<-nrow(gff.sub)
}

#intergenic.min<-intergenic[intergenic$genecount>0,] %>% group_by(chr) %>% top_n(-1, distance)
#intergenic.min<-intergenic.min[,c("chr", "distance")]
#intergenic.min$distance[intergenic.min$distance>10^7]=10^7
#intergenic.min$distance[intergenic.min$distance>10^6]=10^6
#names(intergenic.min)[names(intergenic.min)=="distance"]<-"unit"

intergenic$distance1<-NA
intergenic$distance1[intergenic$genecount==0]=0
intergenic$distance1[intergenic$genecount>0]=180
intergenic$distance1[intergenic$genecount>10]=180*1.5
intergenic$distance1[intergenic$genecount>50]=180*2
intergenic$distance1[intergenic$genecount>100]=180*3
intergenic$distance1[intergenic$genecount>1000]=180*5

intergenic.new<-intergenic[,c("chr", "start", "genecount", "distance", "distance1")]
intergenic.new$Genegroup="intergenic"
intergenic.new$distance<-as.numeric(intergenic.new$distance)

for (i in 1:nrow(intergenic.new)) {
  if (intergenic.new$distance[i] < 10^5) {
    intergenic.new$Gene[i]<-paste0(intergenic.new$genecount[i],"/", 
                                   round(intergenic.new$distance[i]/10^3,2), "Kb")
  } else {
    intergenic.new$Gene[i]<-paste0(intergenic.new$genecount[i],"/", 
                                   round(intergenic.new$distance[i]/10^6,2), "Mb")
  }
}

intergenic.new$GeneID<-NA
intergenic.new$strand<-"."
names(intergenic.new)[names(intergenic.new)=="distance1"]<-"genelength"
HGT.new<-HGT[,c("Genegroup","Gene","GeneID",
                "strand","chr","start","genelength")]
HGT.new$genelength<-180
intergenic.new<-intergenic.new[intergenic.new$genelength>0,c("Genegroup","Gene","GeneID",
                                  "strand","chr","start","genelength")]
HGT.new<-rbind(HGT.new,intergenic.new)
HGT.new<-HGT.new[order(HGT.new$chr,HGT.new$start),]
row.names(HGT.new)<-1:nrow(HGT.new)
HGT.new$stop<-NA

HGT.new.final<-HGT.new[0,]
chrlist<-unique(HGT.new$chr)

for (j in 1:length(chrlist)) {
  HGT.sub<-HGT.new[HGT.new$chr==chrlist[j],]

  if (HGT.sub$strand[1]=="-") {
    HGT.sub$start[1]<-360
    HGT.sub$genelength[1]=0-HGT.sub$genelength[1]
  } else {
    HGT.sub$start[1]<-180
  }
  HGT.sub$stop[1]<-HGT.sub$start[1]+HGT.sub$genelength[1]
  if (nrow(HGT.sub)>1) {
    for (i in 2:nrow(HGT.sub)) {
      
      if (HGT.sub$strand[i]=="-") {
        HGT.sub$start[i]<-max(HGT.sub$stop[i-1], HGT.sub$start[i-1])+180
        HGT.sub$genelength[i]=0-HGT.sub$genelength[i]
      } else {
        HGT.sub$start[i]<-max(HGT.sub$stop[i-1], HGT.sub$start[i-1])
        HGT.sub$genelength[i]=HGT.sub$genelength[i]
      }
      HGT.sub$stop[i]<-HGT.sub$start[i]+HGT.sub$genelength[i]
    }
  } else {
    
  }
  
  HGT.new.final<-rbind(HGT.new.final,HGT.sub)
}

HGT.new.final$start<-as.numeric(HGT.new.final$start)
HGT.new.final$stop<-as.numeric(HGT.new.final$stop)

yscale<-1
#xlab = "(MB)"

Genegrouporder<-c("Hox","Nk","NKL","MOX","Gsx",
                  "Dlx","Evx",
                  "BarxL","Msx","Noto",
                  "XloxL","EmxL",
                  "Hhex",
                  "unique","intergenic")

colorset3<-brewer.pal(12, "Set3")
genecolor=c("#B3DE69","#80B1D3","#BEBADA", "#8DD3C7","#FFFFB3",
            "#FB8072", "#FDB462",
            "#FCCDE5", "gray60", "#BC80BD",
            "#CCEBC5","#A65628","red","gray90","white")


a<-ggplot(HGT.new.final, 
          aes(xmin = start/yscale, xmax = stop/yscale, y = chr, 
              fill = factor(Genegroup, levels = Genegrouporder)))+
  geom_gene_arrow(arrowhead_height = grid::unit(12, "mm"),
                  arrow_body_height = grid::unit(10, "mm"))+
  #geom_gene_arrow(data = subset(HGT.new.final, HGT.new.final$Genegroup != "intergenic"))+
  #geom_gene_arrow(data = subset(HGT.new.final, HGT.new.final$Genegroup == "intergenic"),
  #                arrowhead_width = grid::unit(0, "mm"),
  #                arrowhead_height = grid::unit(3, "mm"))+
  geom_gene_label(aes(xmin = start/yscale, xmax = stop/yscale, y = chr, label = Gene, angle = 30),
                  HGT.new.final,
                  align = "centre",
                  grow = F, reflow = F,
                  min.size = 0,
                  padding.x = grid::unit(0, "mm"),
                  padding.y = grid::unit(0, "lines"),
                  height = grid::unit(20, "mm"))+
  facet_wrap(~ chr, scales = "free_y", ncol = 1)+
  labs(x= "", y = "", fill = "Gene groups")+
  scale_fill_manual(values = genecolor)+
  theme_genes()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5, face = "plain"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "right"
  )
a

ggsave(paste0(species,"_unscale_full.png"),
       width = (max(HGT.new.final$stop,HGT.new.final$start)/360+3), height = 1*length(unique(HGT$chr)), units = "in", dpi = 300, limitsize = F)
f=paste0(species,"_unscale_full.pptx")
topptx(a,f, width = (max(HGT.new.final$stop,HGT.new.final$start)/360+3), height = 1*length(unique(HGT$chr)), units = "in")


ggsave(paste0(species,"_unscale.png"),
       width = (max(HGT.new.final$stop,HGT.new.final$start)/180+3), height = 1*length(unique(HGT$chr)), units = "in", dpi = 300, limitsize = F)
f=paste0(species,"_unscale.pptx")
topptx(a,f, width = (max(HGT.new.final$stop,HGT.new.final$start)/180+3), height = 1*length(unique(HGT$chr)), units = "in")

