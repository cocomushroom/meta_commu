## This script extract representatives of OTUs based on their ranking of percentages summarized across replicates
## Use OTU table as input

setwd("/Users/kc178/Documents/HiSeq1269/newtaxa_combine_0119/taxa_sum/")

OTU_raw<-read.delim("culture_otutable.txt", sep="\t", header=T)
OTU<-OTU_raw[,2:10]
#row.names(OTU)=OTU_raw$OTUId
row.names(OTU)=make.names(OTU_raw$Genus,unique=TRUE)
name <- row.names(OTU)

OTU_2 <- apply(OTU, 2, prop.table)
  
colIndex<- vector(mode="list")
colIndex$b<- grep("b.", names(OTU))
colIndex$m<- grep("m.", names(OTU))
colIndex$t<- grep("t.", names(OTU))


OTU3<-t( apply(OTU_2, 1, function(x){
  c(mean(x[colIndex$b]/3),
    mean(x[colIndex$m]/3),
    mean(x[colIndex$t]/3) )
})  )
colnames(OTU3)<- c("B","M","T")

OTU_sum<- apply(OTU3, 1, sum)

otu_order<- order(OTU_sum, decreasing=T)

col_set<-c(b="chocolate4", m="#556B2F", t="#ADFF2F")
otu_order_subset<- otu_order[1:30]
## for culture-table it's 14, the rest 15
OTU_raw[otu_order_subset,14]
write.table(OTU_raw[otu_order_subset,14], row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", "top30_seq_culture_taxa.txt")
write.table(OTU_raw[otu_order_subset,1], row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", "top30_seq_culture_OTUId.txt")



maxY <- 0.45
png("culture.png", width=1000, height=900)
par(mar=c(20,8,6,3), font.axis=3)
barplot(t(OTU3[otu_order_subset,]), main="strains isolated on 2% MEA", 
        cex.main= 2.8, cex.names= 2.0, cex.axis=2.0, 
        names.arg=name[otu_order_subset], las=2, col=col_set, ylim=c(0, maxY), yaxt="n")
axis(2, at=seq(0,maxY,by=0.05), labels=seq(0,maxY,by=0.05), cex.axis=2.0)
dev.off()


maxY <- 0.25
png("culture_all.png", width=1500, height=900)
par(mar=c(19,8,6,3), font.axis=3)
barplot(t(OTU3[otu_order,]), main="culture_all", 
        cex.main= 2.8, cex.names= 0.9, cex.axis=0.7, 
        names.arg=name[otu_order], las=2, col=col_set, ylim=c(0, maxY), yaxt="n")
axis(2, at=seq(0,maxY,by=0.2), labels=seq(0,maxY,by=0.2), cex.axis=2.0)
dev.off()


######
maxY <- 0.25
png("culture_OTUname.png", width=1000, height=900)
par(mar=c(36,5,3,2), font.axis=3)
barplot(t(OTU3[otu_order_subset,]), main="culture_OTUname", 
        cex.main= 2.8, cex.names= 1.5, cex.axis=1.5, 
        names.arg=name[otu_order_subset], las=2, col=col_set, ylim=c(0, maxY), yaxt="n")
axis(2, at=seq(0,maxY,by=0.2), labels=seq(0,maxY,by=0.2), cex.axis=2.0)
dev.off()

