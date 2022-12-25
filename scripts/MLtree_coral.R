library(ape)
library(ggplot2)
library(tidytree)
library(phylotools)
library(phangorn)
library(ggtree)
library(flextable)
library(tidyr)
library(dplyr)
library(stringr)
library(svglite)
library(ggplotify)
library(eoffice)

gene<-"coral_hox_ClustalW_alignment_10_species"
gene<-"coral_hox_alignment"
gene<-"coral_Nk_alignment"
gene<-"coral_Nk_alignment_v1"
tree<-read.tree(paste(gene,".align.fa.treefile.txt", sep = ""))

node.order<-data.frame(node=1:Nnode(tree) + Ntip(tree), node.label = tree$node.label)
node.order$bootstrap=node.order$node.label

node.order<-separate(node.order, node.label,
                     c("SH-aLRT", "bootstrap"), 
                     sep = "/", convert = FALSE, 
                     remove = FALSE)

node.order$`SH-aLRT`<-as.numeric(node.order$`SH-aLRT`)
node.order$bootstrap<-as.numeric(node.order$bootstrap)

##unroot
{
  mytree<-full_join(tree, node.order, by = "node")
  
  ggtree(mytree, layout = "circular", size = 1)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)
  
  ggsave(paste("ML.",gene,".circular.png", sep = ""), 
         width = 50, height = 50, units = "in", dpi = 300, limitsize = FALSE)
  
  
  ggtree(mytree, layout = "circular", size = 1)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)+
    geom_label2(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                    color = bootstrap), 
                nudge_x = 0, nudge_y = 0, label.size=0.5, 
                size=4, color = "black")
  
  ggsave(paste("ML.",gene,".circular.label.png", sep = ""), 
         width = 50, height = 50, units = "in", dpi = 300, limitsize = FALSE)
  
  
  ggtree(mytree, layout = "circular", size = 1)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)+
    geom_nodelab(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                     hjust= 0, vjust = 0),
                 nudge_x = -0.03, nudge_y = 0,
                 size=6, color = "red")
  
  ggsave(paste("ML.",gene,".circular.nodelab.png", sep = ""), 
         width = 50, height = 50, units = "in", dpi = 300, limitsize = FALSE)
  
  
  ggtree(mytree, layout = "circular", size = 1, branch.length = "none")+
    geom_tiplab(size=6)+
    geom_label2(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                    color = bootstrap), 
                nudge_x = 0, nudge_y = 0, label.size=0.5, 
                size=4, color = "black")
  
  ggsave(paste("ML.",gene,".circular-nonbranch.png", sep = ""),
         width = 50, height = 50, units = "in", dpi = 300, limitsize = FALSE)
  
  a<-ggtree(mytree, size = 1, layout = "fan", open.angle = 10)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)+
    geom_nodelab(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                     hjust= 0, vjust = 0),
                 nudge_x = -0.03, nudge_y = 0,
                 size=6, color = "black")
  a
  ggsave(paste("ML.",gene,".fan.png", sep = ""), width = 30, height = 30, units = "in", dpi = 300, limitsize = FALSE)
  
  
  a<-ggtree(mytree, size = 1)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)+
    #geom_text2(aes(subset=!isTip,label=bootstrap, hjust=1.5, vjust = -1),size=5, color = "black")+
    geom_nodelab(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                     hjust= 0, vjust = 0),
                 nudge_x = -0.3, nudge_y = 0.1,
                 size=6, color = "black")+
    xlim(0,20)
  a
  ggsave(paste("ML.",gene,".rectangle.png", sep = ""), width = 20, height = 15, units = "in", dpi = 300)
  f = paste("ML.",gene,".rectangle.pptx", sep = "")
  topptx(a,f, width = 20, height = 15, units = "in")
}

##Midpoint
{
  test<-reroot(tree, )
  
  test<-midpoint(tree)

  mytree<-full_join(test, node.order, by = "node")
  
  ggtree(mytree, layout = "circular", size = 1)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)
  
  ggsave(paste("ML.",gene,".circular.mid.png", sep = ""), width = 40, height = 40, units = "in", dpi = 300, limitsize = FALSE)
  
  
  ggtree(mytree, layout = "circular", size = 1)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)+
    geom_label2(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                    color = bootstrap), 
                nudge_x = 0, nudge_y = 0, label.size=0.5, 
                size=4, color = "black")
  
  ggsave(paste("ML.",gene,".circular.labelmid.png", sep = ""), 
         width = 50, height = 50, units = "in", dpi = 300, limitsize = FALSE)
  
  
  a<-ggtree(mytree, layout = "circular", size = 1)+
    geom_tiplab(size=7, align=TRUE, linesize=.5)+
    geom_nodelab(aes(subset=!isTip & bootstrap >=0, label=bootstrap,
                     hjust= 0, vjust = 0),
                 nudge_x = -0.03, nudge_y = 0,
                 size=6, color = "red")
  a
  ggsave(paste("ML.",gene,".circular.nodelabmid.png", sep = ""), 
         width = 50, height = 50, units = "in", dpi = 300, limitsize = FALSE)
  f = paste("ML.",gene,".circular.nodelabmid.pptx", sep = "")
  topptx(a,f, width = 50, height = 50, units = "in")
  
  ggtree(mytree, layout = "circular", size = 1, branch.length = "none")+
    geom_tiplab(size=6)+
    geom_label2(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                    color = bootstrap), 
                nudge_x = 0, nudge_y = 0, label.size=0.5, 
                size=4, color = "black")
  
  ggsave(paste("ML.",gene,".circular-nonbranchmid.png", sep = ""), 
         width = 50, height = 50, units = "in", dpi = 300, limitsize = FALSE)
  
  a<-ggtree(mytree, size = 1, layout = "fan", open.angle = 10)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)+
    geom_nodelab(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                     hjust= 0, vjust = 0),
                 nudge_x = -0.03, nudge_y = 0,
                 size=6, color = "black")
  a
  ggsave(paste("ML.",gene,".fanmid.png", sep = ""), width = 30, height = 30, units = "in", dpi = 300, limitsize = FALSE)
  
  a<-ggtree(mytree, size = 1)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)+
    #geom_text2(aes(subset=!isTip,label=bootstrap, hjust=1.5, vjust = -1),size=5, color = "black")+
    geom_nodelab(aes(subset=!isTip & bootstrap >=0, label=bootstrap,
                     hjust= 0, vjust = 0),
                 nudge_x = -0.3, nudge_y = 0.1,
                 size=6, color = "black")+
    geom_treescale()+
    xlim(0,20)
  a
  ggsave(paste("ML.",gene,".rectanglemid.png", sep = ""), width = 20, height = 20, units = "in", dpi = 300)
  f = paste("ML.",gene,".rectanglemid.pptx", sep = "")
  topptx(a,f, width = 20, height = 20, units = "in")
  
  a<-ggtree(mytree, size = 1)+
    geom_tiplab(size=8, align=TRUE, linesize=.5)+
    #geom_text2(aes(subset=!isTip,label=bootstrap, hjust=1.5, vjust = -1),size=5, color = "black")+
    geom_nodelab(aes(subset=!isTip & bootstrap >=80, label=bootstrap,
                     hjust= 0, vjust = 0),
                 nudge_x = -0.03, nudge_y = 0,
                 size=6, color = "black")
  a
  ggsave(paste("ML.",gene,".rectanglemid.png", sep = ""), width = 33, height = 44, units = "in", dpi = 300)
}

