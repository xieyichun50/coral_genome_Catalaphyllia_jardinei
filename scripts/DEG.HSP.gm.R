setwd("D:/coral/DEG_HSP_GM/")

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(eoffice)
library(pheatmap)
library(pathview)
library(clusterProfiler)
DEG<-read.delim("Cjar_gene_count_matrix.csv.R23_vs_R29.edgeR.DE_results", header = TRUE)
##Filter DEGs →
DEG$logFC<- (0-DEG$logFC)
DEG$Groups<-"23ºC>29ºC"
DEG$change<-"ns"
DEG$change[DEG$logFC > 1 & DEG$FDR < 0.05]<-"Up-regulated"
DEG$change[DEG$logFC < 1 & DEG$FDR < 0.05]<-"Down-regulated"

write.table(DEG, file = "DEG.allpairs.txt",
            sep = "\t", quote = F, row.names = T)

DEG.BK<-DEG

#pie chart
DEG.summary<-as.data.frame(xtabs(~change+sampleA+sampleB,DEG.BK))
piep<-ggplot(DEG.summary, aes(x="", y=Freq, fill=change))+
  geom_bar(stat = "identity", position = position_fill())+
  scale_fill_manual(values=c("blue", "gray", "red"))+
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), color = "white")+
  labs(fill = "25°C>32°C")+
  coord_polar(theta = "y", start=0)+
  theme_void()

piep

f = "DEG.pie.pptx"
topptx(piep,f, width = 5, height = 5, units = "in")

write.table(DEG.summary, "DEG.summary.txt", row.names = F, quote = F, sep = "\t")

##Volcano plot
a<-ggplot(data = DEG, aes(x = logFC, y = -log10(FDR), 
                          colour = change))+
  geom_point(shape = 20, size = 1)+
  scale_colour_manual(values = c("blue", "grey", "red"))+
  scale_x_continuous(limits = c(-12,12), breaks = seq(-12,12,4))+
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2))+
  labs(x = "log2(Fold change)", y = "-log10(FDR)", legend = "")+
  geom_hline(yintercept = -log10(0.05), color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = -1, color = "grey30", linetype = "dashed")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 12))
a
ggsave("Volcano.png", width = 5, height = 5, units = "in", dpi = 300)
ggsave("Volcano.tiff", width = 5, height = 5, units = "in", dpi = 300)

##Functional annotation
addline_format <- function(x,...){ 
  gsub('_','\n',x) 
} 

GenesKOGpair.1v1<-read.delim("D:/coral/eggnog/KOG.1v1.txt", header = TRUE)
GenesKEGGpair.1v1<-read.delim("D:/coral/eggnog/KEGG.1v1.txt", header = TRUE)
GenesGOpair.1v1<-read.delim("D:/coral/eggnog/GO.1v1.txt", header = TRUE)
Geneskopair.1v1<-read.delim("D:/coral/eggnog/KO.1v1.txt", header = TRUE)

Geneskopair.1v1$ko<-gsub("ko:","",Geneskopair.1v1$ko)

kog2name<-read.delim("D:/3enrichment/kog2name.txt", header = TRUE)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

ko2name<-read.delim("D:/3enrichment/ko2name.txt")
kegg2name<-read.delim("D:/3enrichment/kegg2name.txt")

kegg2ont<- read.delim("D:/3enrichment/kegglevel.AC.txt", 
                      sep = "\t", colClasses = "character")
names(kegg2ont)[2]="ID"

kegg2koID<-read.delim("D:/3enrichment/kegg2koID.txt",
                      sep = "\t", colClasses = "character")

go2name<-read.delim("D:/3enrichment/go2name.txt", 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

DEG<-DEG[DEG$change != "ns",]
DEG$Genes<-row.names(DEG)
ref="gm."

####KOG enrich
{
  KOG.all.1<-compareCluster(Genes ~ Groups+change, 
                            data = DEG, 
                            fun = 'enricher',
                            TERM2GENE = GenesKOGpair.1v1,
                            TERM2NAME = kog2name,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 1,
                            minGSSize = 1,
                            maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KOG.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"
  
  plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"KOGenrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  list(plotdata$Group)
  plotdata$ratio2[plotdata$change=="Down-regulated"]<- (0-plotdata$ratio2[plotdata$change=="Down-regulated"])
  a<-plotdata %>%
    #mutate(Groups = factor(Groups, levels = group.order.sub)) %>% 
    ggplot(aes(x = ratio2, y = Description, 
               #size = abs(ratio2), 
               size = Count, 
               colour = p.adjust))+
    labs(title = "", 
         size = "Gene counts", colour = "Adjust P value",
         x = "Down-regulated  <<     Gene ratio     >>  Up-regulated", 
         y = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_x_continuous(limits = c(-max(abs(plotdata$ratio2)), max(abs(plotdata$ratio2))))+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    geom_vline(xintercept = 0, color = "black", linetype = "solid")+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "black", size = 12), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 12), 
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave(paste0(ref,"KOG.tiff"), width = 10, height = 6, units = "in", dpi = 300)
  ggsave(paste0(ref,"KOG.png"), width = 10, height = 6, units = "in", dpi = 300)
  
  f = paste0(ref,"KOG.pptx")
  topptx(a,f, width = 10, height = 6, units = "in")
}


####KEGG enrich
{
  KEGG.all.1<-compareCluster(Genes ~ Groups+change, 
                             data = DEG, 
                             fun = 'enricher',
                             TERM2GENE = GenesKEGGpair.1v1,
                             TERM2NAME = kegg2name,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             minGSSize = 1,
                             maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KEGG.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)
  
  plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"KEGG.enrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  plotdata<-read.delim(paste0(ref,"KEGG.enrich.txt"),header = T)
  plotdata$ONTOLOGY<-gsub(" ", "\n", plotdata$ONTOLOGY)
  list(plotdata$Group)
  
  plotdata$ratio2[plotdata$change=="Down-regulated"]<- (0-plotdata$ratio2[plotdata$change=="Down-regulated"])
  plotdata1<-subset(plotdata, abs(ratio2)>0 & p.adjust< 1 & ONTOLOGY != "Human\nDiseases")
  a<-plotdata1 %>%
    ggplot(aes(x = ratio2, y = Description, size = Count, colour = p.adjust))+
    labs(title = "", 
         size = "Gene count", colour = "Adjust P value",
         x = "Down-regulated  <<     Gene ratio     >>  Up-regulated", 
         y = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_x_continuous(limits = c(-max(abs(plotdata1$ratio2)), max(abs(plotdata1$ratio2))))+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    geom_vline(xintercept = 0, color = "black", linetype = "solid")+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "black", size = 12), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 12), 
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 8, angle = 0),
          strip.text.y = element_text(size = 8, angle = 0))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")
  
  a
  ggsave(paste0(ref,"KEGG.tiff"), width = 10, height = 10, units = "in", dpi = 300)
  ggsave(paste0(ref,"KEGG.png"), width = 10, height = 10, units = "in", dpi = 300)
  f = paste0(ref,"KEGG.pptx")
  topptx(a,f, width = 10, height = 10, units = "in")
  
  a<-plotdata1 %>%
    ggplot(aes(x = ratio2, y = Description, size = Count, colour = Count))+
    labs(title = "", 
         size = "Gene count", colour = "Gene count",
         x = "Down-regulated  <<     Gene ratio     >>  Up-regulated", 
         y = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_x_continuous(limits = c(-max(abs(plotdata1$ratio2)), max(abs(plotdata1$ratio2))))+
    scale_color_gradient(high = "#FF0000", low = "#0000FF")+
    geom_vline(xintercept = 0, color = "black", linetype = "solid")+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "black", size = 12), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 12), 
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 8, angle = 0),
          strip.text.y = element_text(size = 8, angle = 0))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")
  
  a
  ggsave(paste0(ref,"KEGG.color.tiff"), width = 10, height = 10, units = "in", dpi = 300)
  ggsave(paste0(ref,"KEGG.color.png"), width = 10, height = 10, units = "in", dpi = 300)
  f = paste0(ref,"KEGG.color.pptx")
  topptx(a,f, width = 10, height = 10, units = "in")
}

####ko enrich
{
  ko.all.1<-compareCluster(Genes ~ Groups+change, 
                           data = DEG, 
                           fun = 'enricher',
                           TERM2GENE = Geneskopair.1v1,
                           TERM2NAME = ko2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(ko.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"ko.enrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  list(plotdata$Group)
  
  plotdata$ratio2[plotdata$change=="Down-regulated"]<- (0-plotdata$ratio2[plotdata$change=="Down-regulated"])
  plotdata1<-subset(plotdata, abs(ratio2)>0.01 & p.adjust< 1 )
  a<-plotdata1 %>%
    ggplot(aes(x = ratio2, y = Description, size = Count, colour = p.adjust))+
    labs(title = "", 
         size = "Gene count", colour = "Adjust P value",
         x = "Down-regulated  <<     Gene ratio     >>  Up-regulated", 
         y = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_x_continuous(limits = c(-max(abs(plotdata$ratio2)), max(abs(plotdata$ratio2))))+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    geom_vline(xintercept = 0, color = "black", linetype = "solid")+
    guides(size = guide_legend(order = 1))+
    facet_grid(~Groups, scales = "free", space = "free")+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "black", size = 11), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 8, angle = 0),
          strip.text.y = element_text(size = 8, angle = 0))
  
  a
  ggsave(paste0(ref,"ko.tiff"), width = 12, height = 8, units = "in", dpi = 300)
  ggsave(paste0(ref,"ko.png"), width = 12, height = 8, units = "in", dpi = 300)
  
}

####GO enrich
{
  GO.all.1<-compareCluster(Genes ~ Groups+change, 
                           data = DEG, 
                           fun = 'enricher',
                           TERM2GENE = GenesGOpair.1v1,
                           TERM2NAME = go2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 2000000)
  
  #Plot
  plotin<-as.data.frame(GO.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[4]<-"goClass"
  go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)
  
  plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"GOenrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  plotdata1<-plotdata
  GO.L4<-read.delim("D:/3enrichment/GOlevel4.txt")
  plotdata1<-plotdata[-which(plotdata$goClass %in% GO.L4$goClass),]
  
  GO.L5<-read.delim("D:/3enrichment/GOlevel5.txt")
  plotdata1<-plotdata1[-which(plotdata1$goClass %in% GO.L5$goClass),]
  
  GO.L6<-read.delim("D:/3enrichment/GOlevel6.txt")
  plotdata1<-plotdata[-which(plotdata1$goClass %in% GO.L6$goClass),]
  
  GO.L7<-read.delim("D:/3enrichment/GOlevel7.txt")
  plotdata1<-plotdata[which(plotdata$goClass %in% GO.L7$goClass),]
  
  GO.Ln<-rbind(GO.L5, GO.L6, GO.L7)
 # plotdata1<-plotdata[which(plotdata$goClass %in% GO.Ln$goClass),]
  
  #top 5 function
  plotdata.top5<-plotdata1 %>% group_by(Cluster, ONTOLOGY) %>% top_n(-5, p.adjust)
  #plotdata.top5<-plotdata1 %>% group_by(Cluster) %>% arrange(-ratio2/p.adjust) %>% slice_head(n=5)
  write.table(plotdata.top5, 
              paste0(ref,"GO.top5bypadj.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  ##all GO
  plotdata1<-plotdata.top5
  plotdata1<-subset(plotdata1, abs(ratio2)>0 & p.adjust< 1 )
  plotdata1$ratio2[plotdata1$change=="Down-regulated"]<- (0-plotdata1$ratio2[plotdata1$change=="Down-regulated"])
  #plotdata1<-plotdata1[-which(plotdata1$goClass == "GO:0048856") ,]
  
  a<-plotdata1 %>%
    ggplot(aes(x = ratio2, y = Description, size = abs(ratio2), colour = p.adjust))+
    labs(title = "", 
         size = "Gene ratio", colour = "Adjust P value",
         x = "Down-regulated  <<     Gene ratio     >>  Up-regulated", 
         y = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_x_continuous(limits = c(-max(abs(plotdata1$ratio2)), max(abs(plotdata1$ratio2))))+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    geom_vline(xintercept = 0, color = "black", linetype = "solid")+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "black", size = 11), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  
  ggsave(paste0(ref,"GO.tiff"), width = 16, height = 16, units = "in", dpi = 300)
  ggsave(paste0(ref,"GO.png"), width = 16, height = 16, units = "in", dpi = 300)
  
  
}


####Read in TMM expression matrix
library(edgeR)
rawcount<-read.delim("gene_count_matrix.csv.R23_vs_R29.edgeR.count_matrix", header = TRUE)
TMM<-rawcount
data=as.matrix(TMM) 

####Read in sample group
trait<-read.delim("samples_n_reads_decribed.txt", header = FALSE)
group <- factor(trait[,2])
y<-DGEList(counts=data,group=group)
y <- calcNormFactors(y, method = "TMM")
TMM.norm<-cpm(y)

write.table(TMM.norm, file = "gene_count_TMMmatrix.txt", quote = F, sep = "\t")

TMM.norm.log2<-log2(TMM.norm+1)
write.table(TMM.norm.log2, file = "gene_count_log2TMMmatrix.txt", quote = F, sep = "\t")

TMM.norm.log2<-as.data.frame(TMM.norm.log2)
TMM.norm.log2<-read.delim("gene_count_log2TMMmatrix.txt", header = TRUE)

heatmap.input<-TMM.norm.log2[row.names(TMM.norm.log2) %in% unique(DEG$Genes[DEG$change!="ns"]),]

names(heatmap.input)=c("23°C.rep1","23°C.rep2","23°C.rep3","29°C.rep1","29°C.rep2","29°C.rep3")
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
names(heatmap.input)<-gsub("X","",names(heatmap.input))

#heatmap.input<-scale(heatmap.input, center = T, scale = T)
bk.limit<-ceiling(max(abs(heatmap.input)))
bk.limit<-0.75*bk.limit
bk.limit<-floor(max(abs(heatmap.input)))
bk.limit<-0.8*bk.limit
bk.limit<-1
bk <- c(seq(-bk.limit,-0.1,by=0.1),seq(0,bk.limit,by=0.1))
library(pheatmap)
heatp<-pheatmap(heatmap.input,
                #color = c(colorRampPalette(colors = c("blue","yellow"))(length(bk)/2),
                #          colorRampPalette(colors = c("yellow","red"))(length(bk)/2)),
                #legend_breaks=seq(-bk.limit,bk.limit,1), breaks=bk,
                legend = TRUE,
                scale = "none",
                cluster_rows = TRUE, cluster_cols = FALSE,
                treeheight_row = 100,
                border_color = NA,
                cellwidth = 50, cellheight = 5, 
                angle_col = c("0"),
                width = 8, height = 10,
                show_colnames = TRUE,
                show_rownames = FALSE,
                filename = paste0(ref,"heatmap.hsp.png"))
#filename = "heatmap.DEG.all.png")
#filename = "heatmap.DEG.mean.png")

##pathview mapping
GenesKOpair.1v1<-read.delim("D:/coral/eggnog/KO.1v1.txt", header = TRUE)
DEG<-merge(DEG, GenesKOpair.1v1, all.x = TRUE, by = "Genes")
DEG<-DEG[is.na(DEG$ko)==F,]
KEGG2koID<-read.delim("D:/3enrichment/kegg2koID.txt",header = TRUE)
DEG<-merge(DEG, KEGG2koID, by = "ko", all.x = TRUE)
pathways<-as.data.frame(unique(DEG$KEGG))
names(pathways)[1]="KEGG"
pathways<-merge(pathways, kegg2name, by = "KEGG", all.x = TRUE)
names(pathways)[1]="pathways"
write.table(pathways, "pathways.txt", row.names = F, quote = F, sep = "\t")

DEG.matrix<-as.data.frame(DEG[,c("Genes","logFC")])

DEG.matrix<-merge(DEG.matrix, Geneskopair.1v1, all.x = TRUE, by = "Genes")
DEG.ko.matrix<-as.data.frame(unique(DEG.matrix$ko))
names(DEG.ko.matrix)[1]="ko"

DEG.matrix.sub<-DEG.matrix[,c(2,3)]
names(DEG.matrix.sub)[1]="FC"
DEG.ko.matrix.sub<-DEG.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(FC))
names(DEG.ko.matrix.sub)[2]=names(DEG.matrix)[2]
DEG.ko.matrix<-merge(DEG.ko.matrix, DEG.ko.matrix.sub, all.x = TRUE, by = "ko")  

rm(DEG.matrix.sub, DEG.ko.matrix.sub)

DEG.ko.matrix<-DEG.ko.matrix[is.na(DEG.ko.matrix$ko)==F,]
row.names(DEG.ko.matrix)<-DEG.ko.matrix$ko
DEG.ko.matrix<-subset(DEG.ko.matrix, select = c("logFC"))

write.table(DEG.ko.matrix, paste0(ref,"DEG.ko.matrix.txt"),
            sep = "\t", quote = F, row.names = T)
setwd("pv/")
library(pathview)
for (i in 1:nrow(pathways)) {
  pathwayid<-pathways$pathways[i]
  #pathwayid<-"ko00480"
  pathview(gene.data = DEG.ko.matrix,
           pathway.id = pathwayid, 
           species = "ko", 
           gene.idtype = "KEGG", 
           limit = list(gene = 5), 
           bins = list(gene=10), 
           multi.state = TRUE, 
           na.col="transparent", 
           out.suffix = ref)
}
