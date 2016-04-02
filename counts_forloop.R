setwd("/Users/kc178/Downloads/all_GB+customized_silva_wLRORLR3_1208/cutfirst")
list_tab <- dir(pattern = "*IdxStats.tabular") 
#make sure the files are read in the correct order
list_tab <- list_tab[c(3,2,1,6,5,4,9,8,7)]
ldf <- list()
present <- list()
for (k in 1:length(list_tab)) {
  ldf[[k]] <- read.delim(list_tab[k], header=FALSE)
  names(ldf[[k]]) <- c("seqID", "length", "mapped", "unmapped")
}

for (a in 1:length(ldf)) {
  data0 <- ldf[[a]]
  present[[a]] <- data0[data0$mapped>0, c(1,3)]
  present[[a]]$seqID <- as.character(present[[a]]$seqID)
}

#ll<-vector(length=9, mode="list")

#combine count from all condition/replicate
allcount <- present[[1]][,c(1,2)]
#allcount <- present1[,c(1,2)]
# this conmmand can make sure exactly how many "0" are assigned; otherwise in some case only 1 "0" is assigned
#note that "$new column name" is enough to add a new column
for (i in 3:10) {
  allcount[,i] <- rep(0, NROW(allcount))  
}
names(allcount) <- c("seqID", "t1", "m1", "b1", "t2", "m2", "b2", "t3", "m3", "b3")


numCol <- ncol(allcount)
numColZero <- rep(0, numCol)
for (q in 2:length(present)) {
  for (i in 1:NROW(present[[q]])) {
    index<- present[[q]]$seqID[i]==allcount$seqID
    if(  any(index)){
      allcount[index, 1+q] <- present[[q]]$mapped[i]
    }else {
      allcount <- rbind(allcount, numColZero)
      rownumber <- NROW(allcount)
      allcount[rownumber, 1] <- present[[q]]$seqID[i]
      allcount[rownumber, q+1] <- present[[q]]$mapped[i]
      #    temp <- c(present2$seqID[i], 0, present2$mapped[i])
    }  #this syntax will only return "TRUE"
  }
}

write.csv(allcount, row.names = FALSE, "allcount.csv")

#remove contaminated sequences checked by BLAST & MEGAN
not_fungi <- read.delim("/Users/kc178/Documents/HiSeq1269/community_RDP_NC_J_CF/velvet91/vel91_c95_meganNonfungi", header=FALSE)
not_fungi <- as.character(not_fungi$V1)
remove <- allcount$seqID %in% not_fungi

write.csv(allcount[!remove,], row.names = FALSE, "allcount_removenonfungi.csv")