setwd("D:/coral/DEG/")
save.image("DEG.RData")
load("DEG.RData")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(eoffice)
library(pheatmap)
library(pathview)
library(clusterProfiler)
DEG<-read.delim("regeneration_conclude.txt")
##Filter DEGs →
DEG.BK<-DEG

##Functional annotation
addline_format <- function(x,...){ 
  gsub('_','\n',x) 
} 

GenesKOGpair.1v1<-read.delim("D:/coral/eggnog/KOG.1v1.txt", header = TRUE)
GenesKEGGpair.1v1<-read.delim("D:/coral/eggnog/KEGG.1v1.txt", header = TRUE)
GenesGOpair.1v1<-read.delim("D:/coral/eggnog/GO.1v1.txt", header = TRUE)
Geneskopair.1v1<-read.delim("D:/coral/eggnog/KO.1v1.txt", header = TRUE)

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

go2name<-read.delim("D:/3enrichment/go2name.txt", 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

DEG.bk<-DEG

group.order=c("0vs6", "6vs18", "0vs18")
ref="all"

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
  a<-plotdata %>%
    mutate(Groups = factor(Groups, levels = group.order)) %>% 
    ggplot(aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave(paste0(ref,"KOG.tiff"), width = 12, height = 8, units = "in", dpi = 300)
  ggsave(paste0(ref,"KOG.png"), width = 12, height = 8, units = "in", dpi = 300)
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
  plotdata$ONTOLOGY<-gsub(" ", "\n", plotdata$ONTOLOGY)
  list(plotdata$Group)
  
  plotdata1<-subset(plotdata, ratio2>0 & p.adjust< 1 & ONTOLOGY != "Human\nDiseases")
  a<-plotdata1 %>%
    mutate(Groups = factor(Groups, levels = group.order)) %>% 
    ggplot(aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 11, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))+
    facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")
  
  a
  ggsave(paste0(ref,"KEGG.tiff"), width = 12, height = 12, units = "in", dpi = 300)
  ggsave(paste0(ref,"KEGG.png"), width = 12, height = 12, units = "in", dpi = 300)
  
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
  plotdata$Description<-paste0(plotdata$Description," (", plotdata$ID, ")")
  plotdata1<-subset(plotdata, ratio2>0.01 & p.adjust< 1 )
  a<-plotdata1 %>%
    mutate(Groups = factor(Groups, levels = group.order)) %>% 
    ggplot(aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    facet_grid(~Groups, scales = "free", space = "free")+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 11, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))
  
  a
  ggsave(paste0(ref,"ko.tiff"), width = 12, height = 12, units = "in", dpi = 300)
  ggsave(paste0(ref,"ko.png"), width = 12, height = 12, units = "in", dpi = 300)
  
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
  
  GO.L4<-read.delim("D:/3enrichment/GOlevel4.txt")
  plotdata1<-plotdata[-which(plotdata$goClass %in% GO.L4$goClass),]
  
  GO.L5<-read.delim("D:/3enrichment/GOlevel5.txt")
  plotdata1<-plotdata1[-which(plotdata1$goClass %in% GO.L5$goClass),]
  
  GO.L6<-read.delim("D:/3enrichment/GOlevel6.txt")
  plotdata1<-plotdata[which(plotdata$goClass %in% GO.L6$goClass),]
  
  GO.L7<-read.delim("D:/3enrichment/GOlevel7.txt")
  plotdata1<-plotdata[which(plotdata$goClass %in% GO.L7$goClass),]
  
  ##all GO
  plotdata1<-subset(plotdata, ratio2>0.05 & p.adjust< 0.2 & ONTOLOGY = "Cellular Component")
  plotdata1$Description<-gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced ascorbate as one donor, and incorporation of one atom of oxygen",
                              "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen,\nreduced ascorbate as one donor, and incorporation of one atom of oxygen",
                              plotdata1$Description)
  a<-plotdata1[plotdata1$ONTOLOGY == "Biological Process" & plotdata1$ratio2>0.1,] %>%
  #a<-plotdata1[plotdata1$ONTOLOGY == "Molecular Function",] %>%
    mutate(Groups = factor(Groups, levels = group.order)) %>% 
    ggplot(aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  
  ggsave(paste0(ref,"GO.BP.tiff"), width = 16, height = 40, units = "in", dpi = 300)
  ggsave(paste0(ref,"GO.BP.png"), width = 16, height = 40, units = "in", dpi = 300)
  
  #top 5 function
  plotdata.top5<-plotdata1 %>% group_by(Cluster, ONTOLOGY) %>% top_n(-5, p.adjust)
  #plotdata.top5<-plotdata1 %>% group_by(Cluster) %>% arrange(-ratio2/p.adjust) %>% slice_head(n=5)
  write.table(plotdata.top5, 
              paste0(ref,"GO.top5bypadj.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
}

##Gene name
gene2name<-read.delim("D:/coral/eggnog/gene2name.txt", header = TRUE)
##DEG freq
DEG.freq<-as.data.frame(xtabs(~Genes, DEG))
DEG.freq<-merge(DEG.freq, gene2name, by = "Genes", all.x = T)

DEG.freq.ud<-as.data.frame(xtabs(~Genes+change, DEG))
DEG.freq.ud<-merge(DEG.freq.ud, gene2name, by = "Genes", all.x = T)

#continuous changes
DEG<-DEG.bk
TMM.norm.log2<-read.delim("RSEM.gene.counts.matrix.TMM_log2.xls", header = TRUE)

heatmap.input<-TMM.norm.log2[row.names(TMM.norm.log2) %in% unique(DEG$Genes),
                             c("LC_0h.1","RC_0h.1","NL_6h.1","NR_6h.1","NL_18h.1","NR_18h.1")]
names(heatmap.input)<-c("0h.rep1", "0h.rep2",
                        "6h.rep1", "6h.rep2",
                        "18h.rep1", "18h.rep2")
#heatmap.input<-TMM.mean.log2[row.names(TMM.mean.log2) %in% unique(DEG$Genes),]



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
heatp<-pheatmap(heatmap.input,
                #color = c(colorRampPalette(colors = c("blue","yellow"))(length(bk)/2),
                #          colorRampPalette(colors = c("yellow","red"))(length(bk)/2)),
                #legend_breaks=seq(-bk.limit,bk.limit,1), breaks=bk,
                legend = TRUE,
                scale = "none",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = 6, treeheight_row = 100,
                border_color = NA,
                cellwidth = 50, cellheight = 6, fontsize = 6,
                angle_col = c("0"),
                width = 10, height = 40,
                show_colnames = TRUE,
                show_rownames = TRUE,
                filename = "test6L.png")
#filename = "heatmap.DEG.all.png")
#filename = "heatmap.DEG.mean.png")

heatp.clust<-cbind(heatmap.input, 
                   order = as.list(heatp$tree_row["order"]),
                   cluster = cutree(heatp$tree_row, k=6))
write.table(heatp.clust, "test6.clust.txt", sep = "\t", row.names = TRUE, quote = FALSE)

heatp.clust<-read.delim("test6.clust.txt", header = TRUE)
heatp.clust$Genes<-row.names(heatp.clust)
DEG<-heatp.clust[,c("Genes", "cluster")]
names(DEG)[2]="Groups"
DEG$Groups<-paste0("cluster", DEG$Groups)
ref= "heatmap.test6"

####KOG enrich
{
  KOG.all.1<-compareCluster(Genes ~ Groups, 
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
  a<-ggplot(plotdata, aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave(paste0(ref,"KOG.tiff"), width = 12, height = 8, units = "in", dpi = 300)
  ggsave(paste0(ref,"KOG.png"), width = 12, height = 8, units = "in", dpi = 300)
}


####KEGG enrich
{
  KEGG.all.1<-compareCluster(Genes ~ Groups, 
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
  plotdata$ONTOLOGY<-gsub(" ", "\n", plotdata$ONTOLOGY)
  list(plotdata$Group)
  
  plotdata1<-subset(plotdata, ratio2>0 & p.adjust< 0.2 & ONTOLOGY != "Human\nDiseases")
  plotdata1$ONTOLOGY<-gsub("Environmental\nInformation\nProcessing", 
                           "Environmental\nInformation Processing",
                           plotdata1$ONTOLOGY)
  plotdata1$ONTOLOGY<-gsub("Genetic\nInformation\nProcessing", 
                           "Genetic\nInformation Processing",
                           plotdata1$ONTOLOGY)
  a<-ggplot(plotdata1, aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 11, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")
  
  a
  ggsave(paste0(ref,"KEGG.tiff"), width = 12, height = 12, units = "in", dpi = 300)
  ggsave(paste0(ref,"KEGG.png"), width = 12, height = 12, units = "in", dpi = 300)
  
}

####ko enrich
{
  ko.all.1<-compareCluster(Genes ~ Groups, 
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
  
  plotdata1<-subset(plotdata, ratio2>0.01 & p.adjust< 0.2 )
  a<-ggplot(plotdata1, aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 11, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))
  
  a
  ggsave(paste0(ref,"ko.tiff"), width = 12, height = 12, units = "in", dpi = 300)
  ggsave(paste0(ref,"ko.png"), width = 12, height = 12, units = "in", dpi = 300)
  
}

####GO enrich
{
  GO.all.1<-compareCluster(Genes ~ Groups, 
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
  
  colnames(plotinsep)[names(plotinsep)=="ID"]<-"goClass"
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
  
  GO.L4<-read.delim("D:/3enrichment/GOlevel4.txt")
  plotdata1<-plotdata[-which(plotdata$goClass %in% GO.L4$goClass),]
  
  GO.L5<-read.delim("D:/3enrichment/GOlevel5.txt")
  plotdata1<-plotdata1[-which(plotdata$goClass %in% GO.L5$goClass),]
  
  GO.L6<-read.delim("D:/3enrichment/GOlevel6.txt")
  plotdata1<-plotdata1[-which(plotdata$goClass %in% GO.L6$goClass),]
  
  GO.L7<-read.delim("D:/3enrichment/GOlevel7.txt")
  plotdata1<-plotdata1[-which(plotdata$goClass %in% GO.L7$goClass),]
  
  ##all GO
  plotdata1<-subset(plotdata, ratio2>0.1 & p.adjust< 0.2)
  plotdata1$Description<-gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced ascorbate as one donor, and incorporation of one atom of oxygen",
                              "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen,\nreduced ascorbate as one donor, and incorporation of one atom of oxygen",
                              plotdata1$Description)
  a<-ggplot(subset(plotdata1, ONTOLOGY=="Molecular Function"),
  #a<-ggplot(subset(plotdata1, ONTOLOGY=="Biological Process"),
            aes(x = Groups, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  
  ggsave(paste0(ref,"GO.MF.tiff"), width = 20, height = 16, units = "in", dpi = 300)
  ggsave(paste0(ref,"GO.MF.png"), width = 20, height = 16, units = "in", dpi = 300)
  
  #top 5 function
  plotdata.top5<-plotdata1 %>% group_by(Groups, ONTOLOGY) %>% top_n(-10, p.adjust)
  #plotdata.top5<-plotdata1 %>% group_by(Cluster) %>% arrange(-ratio2/p.adjust) %>% slice_head(n=5)
  write.table(plotdata.top5, 
              paste0(ref,"GO.topbypadj.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
}



##Gene groups
Up.0to6<-DEG[DEG$Groups=="0h > 6h" & DEG$change=="Up",]
Down.after6<-DEG[DEG$sampleA=="48" & DEG$sampleB !="0" & DEG$change=="Down",]
Genelist<-Up.0to6[Up.0to6$Genes %in% Down.after6$Genes, "Genes"]
TMM.mean.log2$Genes<-rownames(TMM.mean.log2)
TMM.mean.log2<-merge(TMM.mean.log2, gene2name, by = "Genes", all.x = TRUE)
TMM.mean.log2.sub<-TMM.mean.log2[TMM.mean.log2$Genes %in% Genelist,]



continue.up<-subset(TMM.mean.log2, `0h` < `6h` & `6h` < `12h` & 
                      `12h` < `24h` & `24h` < `48h`)
continue.up<-continue.up[row.names(continue.up) %in% DEG$Genes,]


##pathview mapping
GenesKOpair.1v1<-read.delim("D:/coral/eggnog/KO.1v1.txt", header = TRUE)
DEG.sex<-merge(DEG.sex, GenesKOpair.1v1, all.x = TRUE, by = "Genes")
DEG.sex<-DEG.sex[is.na(DEG.sex$ko)==F,]
KEGG2koID<-read.delim("D:/3enrichment/kegg2koID.txt",header = TRUE)
DEG.sex<-merge(DEG.sex, KEGG2koID, by = "ko", all.x = TRUE)
pathways<-as.data.frame(unique(DEG.sex$KEGG))
names(pathways)[1]="KEGG"
pathways<-merge(pathways, kegg2name, by = "KEGG", all.x = TRUE)
names(pathways)[1]="pathways"
write.table(pathways, "pathways.txt", row.names = F, quote = F, sep = "\t")

DEG.matrix<-as.data.frame(unique(DEG$Genes))
names(DEG.matrix)[1]="Genes"

for (i in 1:nrow(DEG.sum.sub)) { 
  DEG.sub<-subset(DEG, sampleA == DEG.sum.sub$sampleA[i] & sampleB == DEG.sum.sub$sampleB[i],
                  select = c("Genes", "logFC"))
  names(DEG.sub)[2]=paste0(DEG.sum.sub$sampleA[i],".", DEG.sum.sub$sampleB[i])
  DEG.sub<-unique(DEG.sub)
  DEG.matrix<-merge(DEG.matrix, DEG.sub, all.x = TRUE, by = "Genes")
}

DEG.matrix<-merge(DEG.matrix, GenesKOpair.1v1, all.x = TRUE, by = "Genes")
DEG.matrix<-DEG.matrix[,-1]
DEG.ko.matrix<-as.data.frame(unique(DEG.matrix$ko))
names(DEG.ko.matrix)[1]="ko"

for (i in 1:(ncol(DEG.matrix)-2)){
  DEG.matrix.sub<-DEG.matrix[,c(i,9)]
  names(DEG.matrix.sub)[1]="FC"
  DEG.ko.matrix.sub<-DEG.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(FC))
  names(DEG.ko.matrix.sub)[2]=names(DEG.matrix)[i]
  DEG.ko.matrix<-merge(DEG.ko.matrix, DEG.ko.matrix.sub, all.x = TRUE, by = "ko")
}
rm(DEG.matrix.sub, DEG.ko.matrix.sub)

DEG.ko.matrix<-DEG.ko.matrix[is.na(DEG.ko.matrix$ko)==F,]
row.names(DEG.ko.matrix)<-DEG.ko.matrix$ko
DEG.ko.matrix<-DEG.ko.matrix[,-1]
DEG.ko.matrix<-DEG.ko.matrix[,c(1,2,3,6,4,7,5,8)]
write.table(DEG.ko.matrix, "DEG.ko.matrix.sex.txt",
            sep = "\t", quote = F, row.names = T)
setwd("pv/")
for (i in 1:nrow(pathways)) {
  pathwayid<-pathways$pathways[i]
  pathwayid<-"ko04150"
  pathview(gene.data = DEG.ko.matrix,
           pathway.id = pathwayid, 
           species = "ko", 
           gene.idtype = "KEGG", 
           limit = list(gene = 1), 
           bins = list(gene=10), 
           multi.state = TRUE, 
           na.col="transparent", 
           out.suffix = ref)
}
