setwd("/Users/kc178/Downloads/all_GB+customized_silva_wLRORLR3_1208/cutfirst")
library("DESeq2", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/resources/library")
test <-read.csv("allcount.csv")
countC <- test[,2:10]
row.names(countC) <- test[, 1]
colnames(countC) <- c("t1", "m1", "b1", "t2", "m2", "b2", "t3", "m3", "b3")
#load ("ddsC.RData")
#load ("C_resTM.RData")
#load ("C_resTB.RData")
#load ("C_resMB.RData")

log_con_mean <- matrix(nrow=NROW(countC), ncol=3)
colnames(log_con_mean) <- c("top", "middle", "bottom")
row.names(log_con_mean) <- row.names(countC)
con_index<- list()
con_index[[1]]<- c(1,4,7)
con_index[[2]]<- con_index[[1]]+1
con_index[[3]]<- con_index[[1]]+2

for (i in 1:NROW(countC)) {
  rr <- as.numeric(countC[i,])
  for (k in 1:3) {
    tt <- rr[con_index[[k]]]
    log_con_mean[i,k] <- log2(mean(tt)+1) 
  }
}

heatmap.2(log_con_mean, col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", dendrogram="row", trace="none", margin=c(10, 6))


type <- read.delim("condition_type", header=TRUE)
a <- type[,2:3]
row.names(a) <- as.character(type[,1])
as.character(type[,1])
row.names(a) <- as.character(type[,1])



print(type)
dim(type)
ddsC <- DESeqDataSetFromMatrix(countData = countC,
                               colData = a,
                               design = ~ condition)
ddsC <- ddsC[, c(1,4,7,2,5,8,3,6,9)]
ddsC <- DESeq(ddsC)

#the relvel, droplevel commands are important when changing the comparison groups
#ddsC$condition <- relevel(ddsC$condition, "middle")
#ddsC$condition <- relevel(ddsC$condition, "bottom")
#ddsC$condition <- droplevels(ddsC$condition, "top")

C_resTB<- results(ddsC)
C_resTB
summary(C_resTB)
#change which condition to compare with
ddsC$condition <- relevel(ddsC$condition, "middle")
# C_resTM gives the C_results of top vs. middle group
C_resTM <- results(ddsC)
C_resTM
summary(C_resTM)
# C_resMB gives the C_results of middle vs. bottom group
# sometimes run "top" first and change it to bottom will work 
ddsC$condition <- relevel(ddsC$condition, "top")
ddsC$condition <- relevel(ddsC$condition, "bottom")
C_resMB <- results(ddsC)
C_resMB
summary(C_resMB)
ddsC$condition <- relevel(ddsC$condition, "top")

#C_res <- C_results(ddsC)
#C_res
C_resTBOrdered <- C_resTB[order(C_resTB$padj),]
head(C_resTBOrdered)
summary(C_resTB)
plotMA(C_resTB, main="DESeq2")
C_resMLE <- C_results(ddsC, addMLE=TRUE)
head(C_resMLE, 4)

#Before all the visualization, let's transform and normalize the read countC
rld <- rlog(ddsC)
vsd <- varianceStabilizingTransformation(ddsC)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

#setting up the color & plot library
library("RColorBrewer")
library("gplots")
select <- plist
select <- order(rowMeans(countCs(ddsC,normalized=TRUE)),decreasing=TRUE)
hmcol <- colorRampPalette(brewer.pal(9, "YlOrBr"))(100)

#a <- log2(C_res$baseMean[plist])

#heatmap.2(a, col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", dendrogram="row", trace="none", margin=c(10, 6))

#Normal heatmap
heatmap.2(assay(rld)[C_resTBP_T_index,], col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", dendrogram="row", trace="none", margin=c(10, 6))
heatmap.2(countCs(ddsC,normalized=TRUE)[select,], col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", dendrogram="row", trace="none", margin=c(10,6))
heatmap.2(assay(vsd)[taxonomy_index,], col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", dendrogram="row", trace="none", margin=c(10, 6))

#clustering heatmap
distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(ddsC), paste(condition, type, sep=" : "))
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(blue_yellow), margin=c(10, 13))

#PCA plot

dataPCA <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
data<-dataPCA
percentVar <- round(100 * attr(data, "percentVar"))
data$condition<- factor(data$condition, levels=c("top",    "middle", "bottom") )

ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  #scale_fill_manual(values=c("greenyellow","chocolate4","darkolivegreen")) +
  #scale_fill_manual(values=plotColor) +
  #scale_color_hue(l=40, c=35) +
  geom_point(colour="black",size=4.5)+
  scale_color_manual(values=col2hex(c("greenyellow","darkolivegreen","chocolate4"))) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  ggtitle("LR0R-MP variance")

###
### MDS ###
##MDS plot
sampleDists <- dist( t( assay(vsd) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$condition, vsd$type, sep="-" )
#colnames(sampleDistMatrix) <- c("MDS1", "MDS2","MDS3", "MDS4","MDS5", "MDS6", "MDS7")
colnames(sampleDistMatrix) <- NULL
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
mds$condition<- factor(mds$condition, levels=c("top",    "middle", "bottom") )


#ggplot(mds, aes(X1,X2,color=condition,shape=type)) + geom_point(size=3)
png("LR0R_MP_MDS.png", width=900, height=800)


ggplot(mds, aes(X1, X2, color=condition, shape=type)) +
  geom_point(colour="black",size=8)+
  scale_color_manual(values=col2hex(c("greenyellow","darkolivegreen","chocolate4"))) +
  geom_point(size=7.5) +
  xlab(paste0("MDS1")) +
  ylab(paste0("MDS2"))+
  ggtitle("LR0R-MP MDS")+
  theme(text= element_text(size=33.5),
        axis.text = element_text(size=30), panel.background=element_rect(fill="white",colour="black")) 

dev.off()




##ggplot (PCA)
png("LR0R_MP_PCA.png", width=900, height=800)
ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(colour="black",size=8)+
  scale_color_manual(values=col2hex(c("greenyellow","darkolivegreen","chocolate4"))) +
  geom_point(size=7.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  ggtitle("LR0R-MP PCA")+
  theme(text= element_text(size=33.5),
        axis.text = element_text(size=30), panel.background=element_rect(fill="white",colour="black")) 

dev.off()









##list of genes with p< certain thC_reshold  in Top vs. Bottom comparison
C_resTBP <- C_resTB$padj
C_resTBP_T <- C_resTBP < 0.05
#c: "which" can give the index of whatever are true
#list of genes with p< certain thC_reshold  in Top vs. Bottom comparison
C_resTBP_T_index  <- which(C_resTBP_T)
length(C_resTBP_T_index)

C_resTMP <- C_resTM$padj
C_resTMP_T <- C_resTMP < 0.05
#c: "which" can give the index of whatever are true
#list of genes with p< certain thC_reshold  in Top vs. Bottom comparison
C_resTMP_T_index  <- which(C_resTMP_T)
length(C_resTMP_T_index)

C_resMBP <- C_resMB$padj
C_resMBP_T <- C_resMBP < 0.05
#c: "which" can give the index of whatever are true
#list of genes with p< certain thC_reshold  in Top vs. Bottom comparison
C_resMBP_T_index  <- which(C_resMBP_T)
length(C_resMBP_T_index)

#
## add taxonomy to the heatmap
taxonomy <- read.delim("LR0R_MP_taxa_list.txt", header=TRUE, stringsAsFactors=F)
contig <- taxonomy$ID
genus <- taxonomy$genus
phylum <- taxonomy$phylum
#require("gdata")
genus <- trim(genus)
#make the contig vs. genus relationship dictionary-like 
genus <- as.list(genus)
names(genus) <- contig
genus <- trim(genus)
phylum <- trim(phylum)
#make the contig vs. genus relationship dictionary-like 
genus <- as.list(genus)
phylum<- as.list(phylum)
names(genus) <- contig
names(phylum) <- contig

# C_res & ddsC are both class 4 datatype and therefore will need to use "slots" to check the contents (?slots)
C_resTBP_rownames <- C_resTB@rownames[C_resTBP_T_index]
C_resTMP_rownames <- C_resTM@rownames[C_resTMP_T_index]
C_resMBP_rownames <- C_resMB@rownames[C_resMBP_T_index]
all_plist <- c(C_resTBP_rownames, C_resTMP_rownames, C_resMBP_rownames)
all_plist <- unique(all_plist)

#write.csv(all_plist, row.names = FALSE, "all_plist.csv")

#common listc
common <- abs(C_resTB$log2FoldChange) < log2(2)
commonnames <- C_resTB@rownames[common]

commonP <- C_resTB$padj > 0.05
commonP[is.na(commonP)] <- F
commonPnames <- C_resTB@rownames[commonP]

common2P <- common & commonP
common2Pnames <- C_resTB@rownames[common2P]

##check up- down- regulated ones
sigTB <- C_resTB$padj < 0.05
TBup <- C_resTB$log2FoldChange > log2(4) ## fold change >1=> logfoldchange >0
sigTBup <- sigTB & TBup
sigTBup[is.na(sigTBup)] <- F ## this truns NA to false so only the true ones will be printed out
TBnamesup <- C_resTB@rownames[sigTBup]
gup <- genus[TBnamesup]
write.csv(TBnamesup, row.names = FALSE, "TBnamesup.csv")
length(TBnamesup)

TBdown <- C_resTB$log2FoldChange < -log2(4)
sigTBdown <- sigTB & TBdown
sigTBdown[is.na(sigTBdown)] <- F 
check_fungi <- C_resTB@rownames[sigTBdown] 
test <- check_fungi[check_fungi %in% contig]
want <- C_resTB@rownames %in% test
TBnamesdown <- C_resTB@rownames[want]
length(TBnamesdown)

###

sig <- all_plist[all_plist %in% contig] #remember some are non-fungi (Node282)
nonsig <- contig[! contig %in% all_plist]
taxonomy_index <- C_resTB@rownames %in% sig
nonsig_tax_index <- C_resTB@rownames %in% nonsig
nonsig_tax <- C_resTB@rownames[C_resTB@rownames %in% nonsig]
contig_want <- C_resTB@rownames[C_resTB@rownames %in% sig]
all_genus <-  C_resTB@rownames %in% contig
all_genus_index <- C_resTB@rownames[C_resTB@rownames %in% contig]

####




#Italic
par(font.axis=3)
#change main title size
par(cex.main=0.8)
#genus_index <- genus[contig %in% sig]
heatmap.2(assay(rld)[taxonomy_index,c(1,2,3,7,8,9)], labRow=genus[contig_want], main = "taxa FDR<0.05 rlog", col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(5, 5), cexCol = 1, cexRow=0.8, srtCol=70, keysize=5)
heatmap.2(assay(vsd)[contig_want,c(1,2,3,7,8,9)], labRow=genus[contig_want], main = c("taxa differentially presented:", "LROR-MP"), col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 8), cexCol = 1, cexRow=0.9, srtCol=60, keysize=1.5, density.info="none")

heatmap.2(assay(vsd)[all_genus,], labRow=genus[all_genus_index], main = "taxa all vsd", col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 8), cexCol = 1, cexRow=0.7, srtCol=70, keysize=1)
heatmap.2(assay(vsd)[all_genus,c(1,2,3,7,8,9)], labRow=genus[all_genus_index], main = "taxa all vsd", col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 8), cexCol = 1, cexRow=0.7, srtCol=70, keysize=1)
heatmap.2(assay(vsd)[TBnamesup,], labRow=genus[TBnamesup], main = "top up vsd", col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 13), cexCol = 1, cexRow=1, srtCol=70, keysize=2)
heatmap.2(assay(vsd)[TBnamesup,c(1,2,3,7,8,9)], labRow=genus[TBnamesup], main = "top up vsd", col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 13), cexCol = 1, cexRow=1, srtCol=70, keysize=1)
heatmap.2(assay(vsd)[TBnamesdown,], labRow=genus[TBnamesdown], main = "down up vsd", col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 13), cexCol = 1, cexRow=1, srtCol=70, keysize=1.5)
heatmap.2(assay(vsd)[common,], labRow=genus[commonnames], main = "log2fold < 2", col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(3, 13), cexCol = 1, cexRow=0.8, srtCol=70, keysize=1)
heatmap.2(assay(vsd)[commonP,], labRow=genus[commonPnames], main = "non-significant (FDR <0.05) vsd", col = hmcol, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(2.5, 13), cexCol = 1, cexRow=0.6, srtCol=70, keysize=1)
heatmap.2(assay(vsd)[common2P,], labRow=genus[common2Pnames], main = c("taxa equally presented:", "LROR-MP"), col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram= "row", trace="none", margin=c(3, 8), cexCol = 1, cexRow=0.9, srtCol=70, keysize=1.5, density.info="none")

#test

#test

png("LR0R_MP_2P.png", width=1200, height=1000)
par(font.axis=3)#make label italic#
par(cex.main=3) #change main title size
heatmap.2(assay(vsd)[common2P,], labRow=genus[common2Pnames], 
          main = c("taxa equally presented: LR0R-MP"), 
          col = blue_yellow, Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(6, 24.2), 
          cexCol = 3.2, cexRow=2.9, srtCol=70, keysize=0.75, 
          density.info="none", key.par=list(cex.lab=2.2, cex.axis=2.2, cex.main=2.5),
          RowSideColors=unlist(pcol3[common2Pnames],
                               
          )  )
dev.off()

png("LR0R_MP_FDR.png", width=1200, height=1000)
par(font.axis=3)#make label italic#
par(cex.main=3) #change main title size
heatmap.2(assay(vsd)[contig_want,c(1,2,3,7,8,9)], labRow=genus[contig_want], 
          main = c("taxa differentially presented: LR0R-MP"), col = blue_yellow, 
          Rowv = TRUE, Colv = FALSE, scale="none", 
          dendrogram="row", trace="none", margin=c(6, 24.2), 
          cexCol = 3.2, cexRow=2.9, srtCol=70, keysize=0.75, density.info="none",
          key.par=list(cex.lab=2.2, cex.axis=2.2, cex.main=2.5),
          RowSideColors=unlist(pcol3[contig_want])  )
dev.off()






#find row number & row names that has padj < 0.1 to "plist"
plist <- {}
pnames<- {}
for (i in 1:NROW(C_res)) {
  p <- C_res$padj[i]
  if (!is.na(p)) {
    if (p<0.1) {
      plist <- c(plist, i)
      pnames <- c(pnames, rownames(C_res)[i])
    } 
  }  
}
length(plist)

#use the padj C_results in "C_res" to select rownames to look at heatmap
df <- countC
for  (i in 1:NROW(countC)) {
  name <- rownames(countC)[i]
  if (name %in% pnames) {
    df <- rbind(df, countC[i,])
  } 
}  


df <- df[3828:3956,]
head(df)

plot(C_res$baseMean+1, -log10(C_res$pvalue),
     log="x", xlab="mean of normalized countCs",
     ylab=expC_ression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
save(genus, file="genus.Rdata")
save(ddsC, file="ddsC.RData")
save(C_resTM, file="C_resTM.RData")
save(C_resTB, file="C_resTB.RData")
save(C_resMB, file="C_resMB.RData")


pcol<- unlist(phylum)

trNameToCol<- matrix(c(
  "Basidiomycota", "blue",
  "Ascomycota", "red",
  "Chytridiomycota", "orange",
  "Glomeromycota", "purple",
  "Entomophthoromycota", "yellow", 
  "Mortierellomycotina", "green", 
  "Mucoromycotina", "brown"), byrow=T, ncol=2)

pcol2<-pcol
for(i in 1:NROW(trNameToCol)){  
  pcol2<- gsub(trNameToCol[i,1], trNameToCol[i,2], pcol2)
}
pcol3<-as.list(pcol2)

legend("topright",      
       legend = trNameToCol[,1],
       col = trNameToCol[,2], 
       lty= 1,             
       lwd = 5,           
       cex=.7
)


#check genus list
df_g<- genus[contig_want]
common2p_g <- genus[common2Pnames]