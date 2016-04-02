#Phyla level barplot
setwd("/Users/kc178/Documents/HiSeq1269/newtaxa_combine_0119/taxa_sum/phyla")

culture_Ori<- read.delim("moss_culture_list_otutable_tax_L2.txt", sep="\t", header=TRUE)
LR0R_amp <-read.delim("LR0R_amp_taxa_L2.txt", sep="\t", header=TRUE)

LR3_amp <-read.delim("LR3_amp_taxa_L2.txt", sep="\t", header=TRUE)

LR0R_MA <-read.delim("LR0R_MA_taxa_L2.txt", sep="\t", header=TRUE)
LR0R_MP <-read.delim("LR0R_MP_taxa_L2.txt", sep="\t", header=TRUE)
LR3_MA <-read.delim("LR3_MA_taxa_L2.txt", sep="\t", header=TRUE)
LR3_MP <-read.delim("LR3_MP_taxa_L2.txt", sep="\t", header=TRUE)
LR3_95_MA <- read.delim("LR3_MA_taxa_L2.txt", sep="\t", header=TRUE)



phyla_names <- c(levels(LR0R_amp$ID), levels(LR3_amp$ID), levels(LR0R_MA$ID), levels(LR0R_MP$ID), 
                 levels(LR3_MA$ID), levels(LR3_MP$ID), levels(culture_Ori$ID))
u_phyla_names<-unique(phyla_names)

#color_set <- c(Ascomycota="red", Basidiomycota="blue", Chytridiomycota="orange", 
#               Mortierellomycotina="purple", Mucoromycotina="white", Neocallimastigomycota="cyan",
#               Entomophthoromycota="green", Glomeromycota="yellow", Microsporidiaincertaesedis="black"
#               )
row.names(LR3_95_MA) <- LR3_95_MA$ID
LR3_95_MA<-LR3_95_MA[,2:10]



row.names(culture_Ori) <- culture_Ori$ID
culture<-culture_Ori[,2:10]


row.names(LR0R_MA) <- LR0R_MA$ID
LR0R_MA<-LR0R_MA[,2:10]

row.names(LR0R_MP) <- LR0R_MP$ID
LR0R_MP<-LR0R_MP[,2:10]

row.names(LR3_MP) <- LR3_MP$ID
LR3_MP <-LR3_MP [,2:10]

row.names(LR3_MA) <- LR3_MA$ID
LR3_MA <-LR3_MA [,2:10]

row.names(LR0R_amp) <- LR0R_amp$ID
LR0R_amp<-LR0R_amp[,2:7]
LR0R_amp <- LR0R_amp[,c(3,2,1,6,5,4)]
row.names(LR3_amp) <- LR3_amp$ID
LR3_amp<-LR3_amp[,2:7]
LR3_amp <- LR3_amp[,c(3,2,1,6,5,4)]


#cols <- colorRampPalette(brewer.pal(9, "Set1"))
# to get the real color codes: cols(9)
#ngroups <- length(u_phyla_names) # not working

color_set_1 <- c(Ascomycota="#E41A1C", Basidiomycota="#377EB8", Chytridiomycota="#4DAF4A", 
               Mortierellomycotina="#FF7F00", Mucoromycotina="#FFFF33", Neocallimastigomycota="#A65628",
               Entomophthoromycota="#F781BF", Glomeromycota="#984EA3", Microsporidiaincertaesedis="#999999", Fungi= "darkgreen")

#par(mfrow=c(3,2))
png("LR3_95_MA.png", width=1200, height=800)
barplot(as.matrix(LR3_95_MA), col=color_set_1[row.names(LR3_95_MA)], border=NA, xlim=c(0,14), main="LR3 MA 95", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR3_95_MA), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend

legend(11, 0.98,
       legend =row.names(LR3_95_MA),
       col=color_set_1[row.names(LR3_95_MA)],
       lty= 2,             
       lwd = 12,           
       cex=1.5
)
dev.off()

png("culture.png", width=1200, height=800)
barplot(as.matrix(culture), col=color_set_1[row.names(culture)], border=NA, xlim=c(0,14), main="Strains isolated on 2% MEA", cex.main=2.2, cex.axis=2.4, names.arg = colnames(culture), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend

dev.off()

png("LR0R_MA.png", width=1200, height=800)
barplot(as.matrix(LR0R_MA), col=color_set_1[row.names(LR0R_MA)], border=NA, xlim=c(0,14), main="LR0R MA", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_MA), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend
dev.off()

png("LR0R_MP.png", width=1200, height=800)
barplot(as.matrix(LR0R_MP), col=color_set_1[row.names(LR0R_MP)], border=NA, xlim=c(0,14), main="LR0R MP", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_MP), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend
dev.off()

png("LR3_MA.png", width=1200, height=800)
barplot(as.matrix(LR3_MA), col=color_set_1[row.names(LR3_MA)], border=NA, xlim=c(0,14), main="LR3 MA", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR3_MA), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend
dev.off()

png("LR3_MP.png", width=1200, height=800)
barplot(as.matrix(LR3_MP), col=color_set_1[row.names(LR3_MP)], border=NA, xlim=c(0,14), main="LR3 MP", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR3_MP), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend
dev.off()

png("LR0R_amp.png", width=1200, height=800)
barplot(as.matrix(LR0R_amp), col=color_set_1[row.names(LR0R_amp)], border=NA, xlim=c(0,14), main="LR0R amplicon sequencing", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR0R_amp), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend
dev.off()

png("LR3_amp.png", width=1200, height=800)
barplot(as.matrix(LR3_amp), col=color_set_1[row.names(LR3_amp)], border=NA, xlim=c(0,14), main="LR3 amplicon sequencing", cex.main=2.2, cex.axis=2.4, names.arg = colnames(LR3_amp), cex.names=2.5) # have to match color_set back to the row.names; border=NA, xlim set can spare space on the right to put the legend
dev.off()









#legend not right




mock <-read.delim("mock_forlegend.txt", sep="\t", header=TRUE)
row.names(mock) <- mock$ID
mock<-mock[,2:7]
png("mock.png", width=1200, height=1500)
barplot(as.matrix(mock), col=color_set_1[row.names(mock)], border=NA, xlim=c(0,14), 
        main="mock", cex.main=2.2, cex.axis=2.4, 
        names.arg = colnames(mock), cex.names=2.5) # have to match color_set back to the row.names; xlim set can spare space on the right to put the legend

legend(5, 4,
       legend =row.names(mock),
       col=color_set_1[row.names(mock)],
       lty= 2,             
       lwd = 12,           
       cex=2.0
)
dev.off()