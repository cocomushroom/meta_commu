#Class level barplot
setwd("/Users/kc178/Documents/HiSeq1269/newtaxa_combine_0119/taxa_sum/class")
LR3_95_MA <- read.delim("LR3_MA_taxa_L3.txt", sep="\t", header=TRUE)
culture_Ori<- read.delim("moss_culture_list_otutable_tax_L3.txt", sep="\t", header=TRUE)
LR0R_amp_Ori <-read.delim("LR0R_amp_taxa_L3.txt", sep="\t", header=TRUE)
LR0R_amp_Ori <- LR0R_amp_Ori[,c(1,4,3,2,7,6,5)]
LR3_amp_Ori <-read.delim("LR3_amp_taxa_L3.txt", sep="\t", header=TRUE)
LR3_amp_Ori <- LR3_amp_Ori[,c(1,4,3,2,7,6,5)]
LR0R_MA_Ori <-read.delim("LR0R_MA_taxa_L3.txt", sep="\t", header=TRUE)
LR0R_MP_Ori <-read.delim("LR0R_MP_taxa_L3.txt", sep="\t", header=TRUE)
LR3_MA_Ori <-read.delim("LR3_MA_taxa_L3.txt", sep="\t", header=TRUE)
LR3_MP_Ori <-read.delim("LR3_MP_taxa_L3.txt", sep="\t", header=TRUE)


class_names <- c(levels(LR0R_amp_Ori$ID), levels(LR3_amp_Ori$ID), levels(LR0R_MA_Ori$ID), levels(LR0R_MP_Ori$ID), 
                 levels(LR3_MA_Ori$ID), levels(LR3_MP_Ori$ID),levels(culture_Ori$ID))
u_class_names <- unique(class_names)

row.names(LR3_95_MA) <- LR3_95_MA$ID
LR3_95_MA<-LR3_95_MA[,2:10]

row.names(culture_Ori) <- culture_Ori$ID
culture<-culture_Ori[,2:10]


row.names(LR0R_MA_Ori) <- LR0R_MA_Ori$ID
LR0R_MA<-LR0R_MA_Ori[,2:10]

row.names(LR0R_MP_Ori) <- LR0R_MP_Ori$ID
LR0R_MP<-LR0R_MP_Ori[,2:10]

row.names(LR3_MP_Ori) <- LR3_MP_Ori$ID
LR3_MP <-LR3_MP_Ori[,2:10]

row.names(LR3_MA_Ori) <- LR3_MA_Ori$ID
LR3_MA <-LR3_MA_Ori[,2:10]

row.names(LR0R_amp_Ori) <- LR0R_amp_Ori$ID
LR0R_amp<-LR0R_amp_Ori[,2:7]

row.names(LR3_amp_Ori) <- LR3_amp_Ori$ID
LR3_amp<-LR3_amp_Ori[,2:7]





#cols <- colorRampPalette(brewer.pal(9, "Set1"))
# to get the real color codes: cols(9)
#ngroups <- length(u_phyla_names) # not working

col_class <- c(Agaricomycetes="#FFFF33", Agaricostilbomycetes="#FFFFB3", Archaeorhizomycetes="#BEBADA", 
               Arthoniomycetes="#FB8072", Ascomycota="#80B1D3", Chytridiomycetes="#FDB462", 
               Cystobasidiomycetes="#B3DE69", Leotiomyceta="blue", Dothideomycetes="deepskyblue", 
               Entomophthoromycetes="#D95F02", Eurotiomycetes="#1B9E77", Exobasidiomycetes="#7570B3", 
               Glomeromycetes="pink", Lecanoromycetes="#8DD3C7", Leotiomycetes="#E7298A", 
               Microbotryomycetes="#E6AB02", Pezizomycotina="#A6761D", 
               Mortierellomycotina="#E41A1C", Mucoromycotina="#377EB8", Basidiomycota="black", 
               Orbiliomycetes="#666666", Saccharomycetes="#984EA3", Sordariomycetes="#FF7F00", 
               Tremellomycetes="#66A61E", Ustilaginomycetes="#A65628", Sordariomyceta="#F781BF", 
               Schizosaccharomycetes="#999999", Malasseziomycetes="gray49", Pucciniomycetes="springgreen1", 
               Fungi="darkgreen", Pezizomycetes="purple")


#par(mfrow=c(3,2))
png("LR3_amp.png", width=1200, height=800)
barplot(as.matrix(LR3_95_MA), col=col_class[row.names(LR3_95_MA)], border=NA, xlim=c(0,14), main="LR3 MA 95", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR3_95_MA), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend
legend(11, 0.98,
       legend =row.names(LR3_95_MA),
       col=col_class[row.names(LR3_95_MA)],
       lty= 2,             
       lwd = 12,           
       cex=1.5
)

dev.off()

png("culture.png", width=1200, height=800)
barplot(as.matrix(culture), col=col_class[row.names(culture)], border=NA, xlim=c(0,14), main="Strains isolated on 2% MEA", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_MA), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend
dev.off()

png("LR0R_MA.png", width=1200, height=800)
barplot(as.matrix(LR0R_MA), col=col_class[row.names(LR0R_MA)], border=NA, xlim=c(0,14), main="LR0R MA", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_MA), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend
dev.off()

png("LR0R_MP.png", width=1200, height=800)
barplot(as.matrix(LR0R_MP), col=col_class[row.names(LR0R_MP)], border=NA, xlim=c(0,14), main="LR0R MP", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_MP), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend
dev.off()

png("LR3_MA.png", width=1200, height=800)
barplot(as.matrix(LR3_MA), col=col_class[row.names(LR3_MA)], border=NA, xlim=c(0,14), main="LR3 MA", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR3_MA), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend
dev.off()

png("LR3_MP.png", width=1200, height=800)
barplot(as.matrix(LR3_MP), col=col_class[row.names(LR3_MP)], border=NA, xlim=c(0,14), main="LR3 MP", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR3_MP), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend
dev.off()

png("LR0R_amp.png", width=1200, height=800)
barplot(as.matrix(LR0R_amp), col=col_class[row.names(LR0R_amp)], border=NA, xlim=c(0,14), main="LR0R amplicon sequencing", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_amp), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend
dev.off()

png("LR3_amp.png", width=1200, height=800)
barplot(as.matrix(LR3_amp), col=col_class[row.names(LR3_amp)], border=NA, xlim=c(0,14), main="LR3 amplicon sequencing", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR3_amp), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend
dev.off()

png("LR0R_with_legend.png", width=1200, height=1500)
barplot(as.matrix(LR0R_amp), col=col_class[row.names(LR0R_amp)], border=NA, xlim=c(0,17), main="LR0R amplicon sequencing", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_amp), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend



legend(11, 0.98,
       legend =row.names(LR0R_amp),
       col=col_class[row.names(LR0R_amp)],
       lty= 2,             
       lwd = 12,           
       cex=2.0
)
dev.off()

png("culture_with_legend.png", width=1200, height=1500)
barplot(as.matrix(culture), col=col_class[row.names(culture)], border=NA, xlim=c(0,17), main="LR0R amplicon sequencing", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_amp), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend




## to produce legend
setwd("/Users/kc178/Documents/HiSeq1269/newtaxa_combine_0119/taxa_sum/class")
mock <-read.delim("mock_forlegend.txt", sep="\t", header=TRUE)
row.names(mock) <- mock$ID
mock<-mock[,2:7]
png("mock.png", width=1200, height=1500)
barplot(as.matrix(mock), col=col_class[row.names(mock)], border=NA, xlim=c(0,14), 
        main="mock", cex.main=2.2, cex.axis=2.4, 
        names.arg = colnames(mock), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend

legend(9, 22,
       legend =row.names(mock),
       col=col_class[row.names(mock)],
       lty= 2,             
       lwd = 12,           
       cex=2.0
)
dev.off()