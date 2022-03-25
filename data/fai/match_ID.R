setwd("F://Coral/fai/")
fai<-read.delim("Montipora_capitata.fa.fai.bk", header = FALSE)
names(fai)[1]="NCBI"
assrp<-read.delim("Montipora_capitata.match", header = FALSE)
names(assrp)[1]="db"
names(assrp)[2]="NCBI"
fai.new<-merge(fai, assrp, by = "NCBI", all.x = TRUE)
fai.new<-subset(fai.new, select = c("db", "V2","V3","V4","V5"))
write.table(fai.new, file = "Montipora_capitata.fa.fai", 
            row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE)