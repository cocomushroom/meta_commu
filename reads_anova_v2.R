###this is to test any read/percentage question, 1) reads mapped 2) ascomycetes 3) eurotiomycetes
setwd("/Users/kc178/Documents/HiSeq1269/newtaxa_combine_0119/taxa_sum/")

raw_read<-read.delim("LR3_MP_otutable.txt", sep="\t", header=TRUE)
raw_read_data<- raw_read[,2:10]

# read_N<- vector(mode="list", length=9)
read_N<- apply(raw_read_data, 2, sum)

tmbName <- c("t","m","b")
tmbList <- sapply(tmbName, function(x){
    grep(x, names(read_N))
}, simplify=F)



table_as<- list()
table_sum<- list()

table_as[["total"]] <- read_N # table_as$total <- read_N #t

group_names<- trimws( levels(raw_read$Phylum) )
for(i in 1:length(group_names) ){
    index<- grep(group_names[i], raw_read$Phylum)
    sub_data<- raw_read_data[index, ]
    as_sum <- apply(sub_data, 2, sum)

    table_as[[group_names[i]]] <- as_sum
    table_as[[paste0("non_",group_names[i])]] <- table_as$total - as_sum
    table_sum[[ group_names[i] ]] <- sapply(tmbList, function(x){ 
        sum(table_as[[group_names[i] ]][x])
    })
}


table_asTMB <- list()
for(j in 1:length(names(table_as))){
  table_asTMB[[names(table_as)[j]]] <- sapply(tmbList, function(x){
    sum(table_as[[names(table_as)[j]]][x])
  })  
}

bi<-as.data.frame(table_asTMB)


A_LR3_MA<-chisq.test(bi[,c(2,3)], simulate.p.value = TRUE)
B_LR3_MA<-chisq.test(bi[,c(4,5)], simulate.p.value = TRUE)
C_LR3_MA<-chisq.test(bi[,c(6,7)], simulate.p.value = TRUE)
G_LR3_MA<-chisq.test(bi[,c(8,9)], simulate.p.value = TRUE)
Mo_LR3_MA <-chisq.test(bi[,c(10,11)], simulate.p.value = TRUE)

A_LR3_MP<-chisq.test(bi[,c(2,3)],simulate.p.value = TRUE)
B_LR3_MP<-chisq.test(bi[,c(4,5)],simulate.p.value = TRUE)
C_LR3_MP<-chisq.test(bi[,c(6,7)],simulate.p.value = TRUE)
G_LR3_MP<-chisq.test(bi[,c(8,9)],simulate.p.value = TRUE)
Mo_LR3_MP <-chisq.test(bi[,c(10,11)],simulate.p.value = TRUE)
Mu_LR3_MP <-chisq.test(bi[,c(12,13)],simulate.p.value = TRUE)


A_LR0R_MP<-chisq.test(bi[,c(2,3)], simulate.p.value = TRUE)
B_LR0R_MP<-chisq.test(bi[,c(4,5)], simulate.p.value = TRUE)
E_LR0R_MP<-chisq.test(bi[,c(6,7)], simulate.p.value = TRUE)
G_LR0R_MP<-chisq.test(bi[,c(8,9)], simulate.p.value = TRUE)
Mo_LR0R_MP <-chisq.test(bi[,c(10,11)], simulate.p.value = TRUE)
Mu_LR0R_MP <-chisq.test(bi[,c(12,13)], simulate.p.value = TRUE)

A_LR0R_MA<-chisq.test(bi[,c(2,3)],simulate.p.value = TRUE)
B_LR0R_MA<-chisq.test(bi[,c(4,5)],simulate.p.value = TRUE)
E_LR0R_MA<-chisq.test(bi[,c(6,7)],simulate.p.value = TRUE)
G_LR0R_MA<-chisq.test(bi[,c(8,9)],simulate.p.value = TRUE)
Mo_LR0R_MA<-chisq.test(bi[,c(10,11)],simulate.p.value = TRUE)


chisqPostHoc(A_LR0R_MA, digits=20)

all_combine<-as.data.frame(table_sum)

Xsq<-chisq.test(all_combine)

chisqPostHoc(Xsq, digits=20)

Ch<-chisq.test(all_combine,  simulate.p.value = TRUE)

Ph_LR3_MP_P<-Ch$p.value
Ph_LR3_MP_sq<-Ch$statistic

Ph_LR3_MA_P<-Ch$p.value
Ph_LR3_MA_sq<-Ch$statistic



Ph_LR0R_MA_P<-Ch$p.value
Ph_LR0R_MA_sq<-Ch$statistic

Ph_LR0R_MP_P<-Ch$p.value
Ph_LR0R_MP_sq<-Ch$statistic


chi_phy_result<-data.frame()




##Class_level
class_as<- list()
class_sum<- list()
class_as[["total"]] <- read_N # class_as$total <- read_N #t
class_names<- trimws( levels(raw_read$Class) )

for(i in 1:length(class_names) ){
  index<- grep(class_names[i], raw_read$Class)
  sub_data<- raw_read_data[index, ]
  as_sum <- apply(sub_data, 2, sum)
  
  class_as[[class_names[i]]] <- as_sum
  class_as[[paste0("non_",class_names[i])]] <- class_as$total - as_sum
  
  class_sum[[ class_names[i] ]] <- sapply(tmbList, function(x){ 
    sum(class_as[[class_names[i] ]][x])
  })
}

all_combine_class<-as.data.frame(class_sum)

Xsq<-chisq.test(all_combine_class)

Ch_c<-chisq.test(all_combine_class,  simulate.p.value = TRUE)
Cl_LR0R_MP_P<-Ch_c$p.value
Cl_LR0R_MP_sq<-Ch_c$statistic




### mapped read##
LR3_MP<-c(384,892,674,248,1238,1104,1164,1242,2640)
LR3_MP<-as.data.frame(LR3_MP)
LR3_MP$LR3_MP_un <- c(27043936,41394964,33443752,39876044,59551578,31508724,49390618,40959986,45524538)
row.names(LR3_MP) <- c("t1","m1","b1","t2","m2","b2","t3","m3","b3")
LR3_MP$layer <- c("top","middle","bottom","top","middle","bottom","top","middle","bottom")
sum_map_t <- sum(LR3_MP$LR3_MP[c(1,4,7)])
sum_map_m <- sum(LR3_MP$LR3_MP[c(2,5,8)])
sum_map_b <- sum(LR3_MP$LR3_MP[c(3,6,9)])
un_map_t <- sum(LR3_MP$LR3_MP_un[c(1,4,7)])
un_map_m <- sum(LR3_MP$LR3_MP_un[c(2,5,8)])
un_map_b <- sum(LR3_MP$LR3_MP_un[c(3,6,9)])
map <- as.data.frame(c(sum_map_t, sum_map_m, sum_map_b))
map$unmap<- c(un_map_t,un_map_m, un_map_b)
row.names(map) <- c("top","middle","bottom")
colnames(map) <- c("map","unmap")
test_LR3MP<-chisq.test(map)
chisqPostHoc(test_LR3MP, digits=20)

LR3_MA<-c(684,1678,1316,408,1986,1380,2104,1958,3422)
LR3_MA<-as.data.frame(LR3_MA)
LR3_MA$LR3_MA_un<-c(27043636,41394178,33443110,39875884,59550830,31508448,49389678,40959270,45523756)
row.names(LR3_MA) <- c("t1","m1","b1","t2","m2","b2","t3","m3","b3")
LR3_MA$layer <- c("top","middle","bottom","top","middle","bottom","top","middle","bottom")
sum_map_t <- sum(LR3_MA$LR3_MA[c(1,4,7)])
sum_map_m <- sum(LR3_MA$LR3_MA[c(2,5,8)])
sum_map_b <- sum(LR3_MA$LR3_MA[c(3,6,9)])
un_map_t <- sum(LR3_MA$LR3_MA_un[c(1,4,7)])
un_map_m <- sum(LR3_MA$LR3_MA_un[c(2,5,8)])
un_map_b <- sum(LR3_MA$LR3_MA_un[c(3,6,9)])
map <- as.data.frame(c(sum_map_t, sum_map_m, sum_map_b))
map$unmap<- c(un_map_t,un_map_m, un_map_b)
row.names(map) <- c("top","middle","bottom")
colnames(map) <- c("map","unmap")
test_LR3MA<-chisq.test(map)
chisqPostHoc(test_LR3MA, digits=20)

LR3_Pr<-c(872,1822,1688,442,2524,1980,1924,2602,5158)
LR3_Pr<-as.data.frame(LR3_Pr)
LR3_Pr$LR3_Pr_un <- c(27043448,41394034,33442738,39875850,59550292,31507848,49389858,40958626,45522020)
row.names(LR3_Pr) <- c("t1","m1","b1","t2","m2","b2","t3","m3","b3")
LR3_Pr$layer <- c("top","middle","bottom","top","middle","bottom","top","middle","bottom")
sum_Prp_t <- sum(LR3_Pr$LR3_Pr[c(1,4,7)])
sum_Prp_m <- sum(LR3_Pr$LR3_Pr[c(2,5,8)])
sum_Prp_b <- sum(LR3_Pr$LR3_Pr[c(3,6,9)])
un_Prp_t <- sum(LR3_Pr$LR3_Pr_un[c(1,4,7)])
un_Prp_m <- sum(LR3_Pr$LR3_Pr_un[c(2,5,8)])
un_Prp_b <- sum(LR3_Pr$LR3_Pr_un[c(3,6,9)])
Prp <- as.data.frame(c(sum_Prp_t, sum_Prp_m, sum_Prp_b))
Prp$unPrp<- c(un_Prp_t,un_Prp_m, un_Prp_b)
row.names(Prp) <- c("top","middle","bottom")
colnames(Prp) <- c("Prp","unPrp")
test_LR3Pr<-chisq.test(Prp)
chisqPostHoc(test_LR3Pr, digits=20)


LR0R_MP<-c(444,718,606,214,926,738,902,1136,2070)
LR0R_MP<-as.data.frame(LR0R_MP)
LR0R_MP$LR0R_MP_un<-c(27043876,41395138,33443820,39876078,59551890,31509090,49390880,40960092,45525108)
row.names(LR0R_MP) <- c("t1","m1","b1","t2","m2","b2","t3","m3","b3")
LR0R_MP$layer <- c("top","middle","bottom","top","middle","bottom","top","middle","bottom")
sum_MPp_t <- sum(LR0R_MP$LR0R_MP[c(1,4,7)])
sum_MPp_m <- sum(LR0R_MP$LR0R_MP[c(2,5,8)])
sum_MPp_b <- sum(LR0R_MP$LR0R_MP[c(3,6,9)])
un_MPp_t <- sum(LR0R_MP$LR0R_MP_un[c(1,4,7)])
un_MPp_m <- sum(LR0R_MP$LR0R_MP_un[c(2,5,8)])
un_MPp_b <- sum(LR0R_MP$LR0R_MP_un[c(3,6,9)])
Prp <- as.data.frame(c(sum_MPp_t, sum_MPp_m, sum_MPp_b))
Prp$unPrp<- c(un_MPp_t,un_MPp_m, un_MPp_b)
row.names(Prp) <- c("top","middle","bottom")
colnames(Prp) <- c("Prp","unPrp")
test_LR0RPr<-chisq.test(Prp)
chisqPostHoc(test_LR0RPr, digits=20)

LR0R_MA<-c(980,2278,1988,574,3210,2134,2820,2640,5002)
LR0R_MA<-as.data.frame(LR0R_MA)
LR0R_MA$LR0R_MA_un<-c(27043340,41393578,33442438,39875718,59549606,31507694,49388962,40958588,45522176)
row.names(LR0R_MA) <- c("t1","m1","b1","t2","m2","b2","t3","m3","b3")
LR0R_MA$layer <- c("top","middle","bottom","top","middle","bottom","top","middle","bottom")
sum_map_t <- sum(LR0R_MA$LR0R_MA[c(1,4,7)])
sum_map_m <- sum(LR0R_MA$LR0R_MA[c(2,5,8)])
sum_map_b <- sum(LR0R_MA$LR0R_MA[c(3,6,9)])
un_map_t <- sum(LR0R_MA$LR0R_MA_un[c(1,4,7)])
un_map_m <- sum(LR0R_MA$LR0R_MA_un[c(2,5,8)])
un_map_b <- sum(LR0R_MA$LR0R_MA_un[c(3,6,9)])
map <- as.data.frame(c(sum_map_t, sum_map_m, sum_map_b))
map$unmap<- c(un_map_t,un_map_m, un_map_b)
row.names(map) <- c("top","middle","bottom")
colnames(map) <- c("map","unmap")
test_LR0RMA<-chisq.test(map)
chisqPostHoc(test_LR0RMA, digits=20)

LR0R_Pr<-c(710,1294,1210,386,1848,1254,1590,1782,3462)
LR0R_Pr<-as.data.frame(LR0R_Pr)
LR0R_Pr$LR0R_Pr_un<-c(27043610,41394562,33443216,39875906,59550968,31508574,49390192,40959446,45523716)
row.names(LR0R_Pr) <- c("t1","m1","b1","t2","m2","b2","t3","m3","b3")
LR0R_Pr$layer <- c("top","middle","bottom","top","middle","bottom","top","middle","bottom")
sum_map_t <- sum(LR0R_Pr$LR0R_Pr[c(1,4,7)])
sum_map_m <- sum(LR0R_Pr$LR0R_Pr[c(2,5,8)])
sum_map_b <- sum(LR0R_Pr$LR0R_Pr[c(3,6,9)])
un_map_t <- sum(LR0R_Pr$LR0R_Pr_un[c(1,4,7)])
un_map_m <- sum(LR0R_Pr$LR0R_Pr_un[c(2,5,8)])
un_map_b <- sum(LR0R_Pr$LR0R_Pr_un[c(3,6,9)])
map <- as.data.frame(c(sum_map_t, sum_map_m, sum_map_b))
map$unmap<- c(un_map_t,un_map_m, un_map_b)
row.names(map) <- c("top","middle","bottom")
colnames(map) <- c("map","unmap")
test_LR0RPr<-chisq.test(map)
chisqPostHoc(test_LR0RPr, digits=20)


#### ANOVA map ####
LR3_MP_map<- c(1.41989E-05,2.15481E-05,2.01528E-05,6.21923E-06,2.07883E-05,3.50367E-05,2.35667E-05,3.03214E-05,5.79873E-05)
LR3_MP_map <- as.data.frame(LR3_MP_map)
row.names(LR3_MP_map) <- c("t1","m1","b1","t2","m2","b2","t3","m3","b3")
LR3_MP_map$layer <- c("top","middle","bottom","top","middle","bottom","top","middle","bottom")
result<-lm(LR3_MP_map~layercode,LR3_MP_map)

summary(aov(result))
TukeyHSD(aov(mo))
mo_P <- result_mo[[1]][["Pr(>F)"]][[1]]



####IGNORE######
#################################################

t <- sum(table_as$as_sum[c(1,4,7)])
m <- sum(table_as$as_sum[c(2,5,8)])
b <- sum(table_as$as_sum[c(3,6,9)])
A <- c(t,m,b)

allT<-sum(table_as$total[c(1,4,7)])
allM<-sum(table_as$total[c(2,5,8)])
allB<-sum(table_as$total[c(3,6,9)])
TOTAL<-c(allT,allM,allB)


nt<-sum(table_as$nonA[c(1,4,7)])
nm<-sum(table_as$nonA[c(2,5,8)])
nb<-sum(table_as$nonA[c(3,6,9)])
nonA<-c(nt,nm,nb)

M <- data.frame(A)
M$nonA<-nonA
row.names(M) <- c("top","middle","bottom")

Xsq<-chisq.test(all_combine)

prop.test(A,TOTAL)
pairwise.prop.test(A,TOTAL)

## test Basidiomycota ###

Basi_ori<- grep("Basidiomycota", raw_read$Phylum)
Basi<-raw_read[Basi_ori,c(2:10)]
ba_sum <- apply(Basi, 2, function(x){
  sum(x)
})
table_as$B<-ba_sum
t <- sum(table_as$B[c(1,4,7)])
m <- sum(table_as$B[c(2,5,8)])
b <- sum(table_as$B[c(3,6,9)])
B <- c(t,m,b)

prop.test(B,TOTAL)
pairwise.prop.test(B,TOTAL)

### Test for Mortierellomycotina ###
Mo_ori<- grep("Mortierellomycotina", raw_read$Phylum)
Mo<-raw_read[Mo_ori,c(2:10)]
Mo_sum <- apply(Mo, 2, function(x){
  sum(x)
})
table_as$Mo<-ba_sum
t <- sum(table_as$Mo[c(1,4,7)])
m <- sum(table_as$Mo[c(2,5,8)])
b <- sum(table_as$Mo[c(3,6,9)])
Mo <- c(t,m,b)
prop.test(Mo,TOTAL)
pairwise.prop.test(Mo,TOTAL)

## Glomeromycota
G_ori<- grep("Glomeromycota", raw_read$Phylum)
G<-raw_read[Mo_ori,c(2:10)]
G_sum <- apply(G, 2, function(x){
  sum(x)
})
table_as$G<-G_sum
t <- sum(table_as$G[c(1,4,7)])
m <- sum(table_as$G[c(2,5,8)])
b <- sum(table_as$G[c(3,6,9)])
G <- c(t,m,b)
prop.test(G, TOTAL)
pairwise.prop.test(G,TOTAL)


## Mucoromycotina
Mu_ori<- grep("Mucoromycotina", raw_read$Phylum)
Mu<-raw_read[Mu_ori,c(2:10)]
Mu_sum <- apply(Mu, 2, function(x){
  sum(x)
})
table_as$Mu<-Mu_sum
t <- sum(table_as$Mu[c(1,4,7)])
m <- sum(table_as$Mu[c(2,5,8)])
b <- sum(table_as$Mu[c(3,6,9)])
Mu <- c(t,m,b)
prop.test(Mu, TOTAL)
pairwise.prop.test(Mu,TOTAL)


## Entomophromycota
En_ori<- grep("Entomophthoromycota", raw_read$Phylum)
En<-raw_read[En_ori,c(2:10)]
En_sum <- apply(En, 2, function(x){
  sum(x)
})
table_as$En<-En_sum
t <- sum(table_as$En[c(1,4,7)])
m <- sum(table_as$En[c(2,5,8)])
b <- sum(table_as$En[c(3,6,9)])
En <- c(t,m,b)
prop.test(En, TOTAL)
pairwise.prop.test(En,TOTAL)


## chisqure test on the entire table##
all_phy <- data.frame(A)
all_phy <- data.frame(cbind(A,B,Mo,Mu,En, G))
all_phy$B <- B
all_phy$Mo <- Mo
all_phy$Mu <- Mu
all_phy$En <- En
all_phy$G <- G
row.names(all_phy) <- c("top","middle","bottom")
Ch<-chisq.test(all_phy,  simulate.p.value = TRUE)
chisqPostHoc(Ch, digits=10)

chisqPostHoc2<- function(chi, popsInRows = TRUE, control = stats::p.adjust.methods, 
                         digits = 4) 
{
  control <- match.arg(control)
  tbl <- chi$observed
  if (!popsInRows) 
    tbl <- t(tbl)
  popsNames <- rownames(tbl)
  prs <- utils::combn(1:nrow(tbl), 2)
  tests <- ncol(prs)
  pvals <- numeric(tests)
  lbls <- character(tests)
  for (i in 1:tests) {
    pvals[i] <- stats::chisq.test(tbl[prs[, i], ], simulate.p.value = TRUE, B=1e6)$p.value
    lbls[i] <- paste(popsNames[prs[, i]], collapse = " vs. ")
  }
  adj.pvals <- stats::p.adjust(pvals, method = control)
  cat("Adjusted p-values used the", control, "method.\n\n")
  data.frame(comparison = lbls, raw.p = round(pvals, digits), 
             adj.p = round(adj.pvals, digits))
}




################IGNORE###########

##############################################
## chisaure test on percentage normalized table##
per_o<-read.delim("/Users/kc178/Documents/HiSeq1269/newtaxa_combine_0119/taxa_sum/phyla/LR0R_MP_taxa_L2.txt", sep="\t", header=TRUE)
ID<-per_o[1:6,1]
per<-t(per_o[,-1])
colnames(per)<- ID
per_P<-all_phy
colnames(per_P)<-colnames(per)

for (i in 1:NCOL(per)) {
  per_P[1,i]<-(sum(per[c(1,4,7),i]))/3
  per_P[2,i]<-(sum(per[c(2,5,8),i]))/3
  per_P[3,i]<-(sum(per[c(3,6,9),i]))/3  
}
chisq.test(per_P)




######### Working on ratio #######
##########        ################
##phylum
setwd("/Users/kc178/Documents/HiSeq1269/newtaxa_combine_0119/taxa_sum/phyla")
number<- read.delim("moss_culture_list_otutable_tax_L2.txt", sep="\t", header=TRUE)
#tnumber <- t(number)
tnumber<-t(number[,-1])
colnames(tnumber)<-number$ID
tnumber <- data.frame(tnumber)
layercode <- c("top", "middle","bottom","top","middle","bottom", "top", "middle", "bottom")
tnumber <- cbind(tnumber,layercode)
tnumber[tnumber==0.0000000000] <- NA


name <- colnames(tnumber)
pvalue <- list()

for (i in 1:(NCOL(tnumber)-1)  ) {
  testT <- is.na( tnumber[c(1,4,7), i])
  testT<-sum(testT==FALSE)
  testM <- is.na( tnumber[c(2,5,8), i])
  testM<-sum(testM==FALSE)
  testB <- is.na( tnumber[c(3,6,9), i])
  testB<-sum(testB==FALSE)
  if (testT>0 & testM>0 & testB>0) {
    vari<-fligner.test(tnumber[,i]~layercode, data=tnumber)$p.value
    nor<-shapiro.test(tnumber[,i])$p.value
    if (vari>0.05 & nor>0.05) {
      test <- lm(tnumber[,i]~layercode,tnumber)
      result<-summary(aov(test))
      pvalue[[name[i] ]] <- result[[1]][["Pr(>F)"]][[1]]
  
      if (result[[1]][["Pr(>F)"]][[1]] <0.05) {
  
        print (names(pvalue)[i])
  }
    
    }
  }
}

data_pvalue <- as.numeric(pvalue)
FDR<- p.adjust(data_pvalue, "fdr")
FDR<-as.data.frame(FDR)
row.names(FDR)<-names(pvalue)

for (j in 1:NROW(FDR)) {
  if (FDR[j,] <0.05) {
    print (row.names(FDR)[j])
  } 
}




shapiro.test(tnumber$Chytridiomycota)
fligner.test(Chytridiomycota~layercode, data=tnumber)

shapiro.test(tnumber$Basidiomycota)
fligner.test(Basidiomycota~layercode, data=tnumber)

shapiro.test(tnumber$Entomophthoromycota)
fligner.test(Entomophthoromycota~layercode, data=tnumber)

shapiro.test(tnumber$Glomeromycota)
fligner.test(Glomeromycota~layercode, data=tnumber)

shapiro.test(tnumber$Mortierellomycotina)
fligner.test(Mortierellomycotina~layercode, data=tnumber)

shapiro.test(tnumber$Mucoromycotina)
fligner.test(Mucoromycotina~layercode, data=tnumber)


asco <- lm(Ascomycota~layercode,tnumber)
result_a <- summary(aov(asco))
(TukeyHSD(aov(asco))$layercode)[,4]
as_P <- result[[1]][["Pr(>F)"]][[1]]


basi <- lm(Basidiomycota~layercode,tnumber)
result_b <- summary(aov(basi))
(TukeyHSD(aov(basi))$layercode)[,4]
basi_P <- result_b[[1]][["Pr(>F)"]][[1]]

en <- lm(Entomophthoromycota~layercode,tnumber)
result_en <- summary(aov(en))
TukeyHSD(aov(en))
en_P <- result_en[[1]][["Pr(>F)"]][[1]]

g <- lm(Glomeromycota~layercode,tnumber)
result_g <- summary(aov(g))
TukeyHSD(aov(g))
g_P <- result_g[[1]][["Pr(>F)"]][[1]]


mo <- lm(Mortierellomycotina~layercode,tnumber)
result_mo <- summary(aov(mo))
TukeyHSD(aov(mo))
mo_P <- result_mo[[1]][["Pr(>F)"]][[1]]

p.adjust(p, method = p.adjust.methods, n = length(p))

p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",






### classes###
setwd("/Users/kc178/Documents/HiSeq1269/newtaxa_combine_0119/taxa_sum/class")
number<- read.delim("LR0R_MA_taxa_L3.txt", sep="\t", header=TRUE)
tnumber<-t(number[,-1])
colnames(tnumber)<-number$ID
tnumber <- data.frame(tnumber)
layercode <- c("top", "middle","bottom","top","middle","bottom", "top", "middle", "bottom")
tnumber <- cbind(tnumber,layercode)
tnumber[tnumber==0.0000000000] <- NA


name <- colnames(tnumber)
pvalue <- list()




for (i in 1:(NCOL(tnumber)-1)  ) {
  testT <- is.na( tnumber[c(1,4,7), i])
  testT<-sum(testT==FALSE)
  testM <- is.na( tnumber[c(2,5,8), i])
  testM<-sum(testM==FALSE)
  testB <- is.na( tnumber[c(3,6,9), i])
  testB<-sum(testB==FALSE)
  if (testT>0 & testM>0 & testB>0) {
    vari<-fligner.test(tnumber[,i]~layercode, data=tnumber)$p.value
    nor<-shapiro.test(tnumber[,i])$p.value
    if (vari>0.05 & nor>0.05) {
      test <- lm(tnumber[,i]~layercode,tnumber)
      result<-summary(aov(test))
      pvalue[[name[i] ]] <- result[[1]][["Pr(>F)"]][[1]]
      
      if (result[[1]][["Pr(>F)"]][[1]] <0.05) {
        
        print (names(pvalue)[i])
      }
      
    }
  }
}





data_pvalue <- as.numeric(pvalue)
FDR<- p.adjust(data_pvalue, "fdr")
FDR<-as.data.frame(FDR)
row.names(FDR)<-names(pvalue)

for (j in 1:NROW(FDR)) {
  if (FDR[j,] <0.05) {
    print (row.names(FDR)[j])
  } 
}


##LR0R_MP
do <- lm(Dothideomycetes~layercode,tnumber)
result_do <- summary(aov(do))
TukeyHSD(aov(do))

##LR3_MP
tr <- lm(Tremellomycetes~layercode,tnumber)
result_tr <- summary(aov(tr))
TukeyHSD(aov(tr))

ch <- lm(Chytridiomycetes~layercode,tnumber)
result_ch <- summary(aov(ch))
TukeyHSD(aov(ch))


shapiro.test(tnumber$Dothideomycetes)
fligner.test(Dothideomycetes~layercode, data=tnumber)

  
  index<- grep(group_names[i], raw_read$Phylum)
  sub_data<- raw_read_data[index, ]
  as_sum <- apply(sub_data, 2, sum)
  
  table_as[[group_names[i]]] <- as_sum
  table_as[[paste0("non_",group_names[i])]] <- table_as$total - as_sum
  table_sum[[ group_names[i] ]] <- sapply(tmbList, function(x){ 
    sum(table_as[[group_names[i] ]][x])
  })
}






test_LR0R_MA_A <- lm(Ascomycota~layercode,tnumber)
result <- summary(aov(test_LR0R_MA_A))
TukeyHSD(aov(test_logP))


test_LR0R_MA_B <- lm(Basidiomycota~layercode,tnumber)
summary(aov(test_LR0R_MA_B))

test_logA <- lm(log_MA~layer,Rnumber)
summary(aov(test_logA))
TukeyHSD(aov(test_logA))

test_logMP <- lm(log_MP~layer,Rnumber)
summary(aov(test_logMP))
TukeyHSD(aov(test_logMP))

test_p <- lm(primary~layer, Rnumber)
summary(aov(test_p))
TukeyHSD(aov(test_p))


test_MA <- lm(MA~layer,Rnumber)
summary(aov(test_MA))
TukeyHSD(aov(test_MA))
a_test_MA <- anova(test_MA)
a_test_MA


test_MP <- lm(MP~layer, Rnumber)
a_test_MP <- anova(test_MP)
summary(aov(test_MP))
TukeyHSD(aov(test_MP))

asco<-read.delim("asco_LRORMA.txt", header=TRUE, sep="\t")
test_OMA <-lm(LROR_MA~layer, asco)
summary(aov(test_OMA))
TukeyHSD(aov(test_OMA))

asco<-read.delim("asco_LRORMA.txt", header=TRUE, sep="\t")
test_OMP <-lm(LROR_MP~layer, asco)
summary(aov(test_OMP))
TukeyHSD(aov(test_OMP))

asco<-read.delim("asco_LRORMA.txt", header=TRUE, sep="\t")
test_3MA <-lm(LR3_MA~layer, asco)
summary(aov(test_3MA))
TukeyHSD(aov(test_3MA))

asco<-read.delim("asco_LRORMA.txt", header=TRUE, sep="\t")
test_3MP <-lm(LR3_MP~layer, asco)
summary(aov(test_3MP))
TukeyHSD(aov(test_3MP))

test_cul <-lm(culture~layer, asco)
summary(aov(test_cul))
TukeyHSD(aov(test_cul))

eu <- read.delim("euro_percentage.txt")
eu_oma <- lm(LROR_MA~layer, eu)
summary(aov(eu_oma))
TukeyHSD(aov(eu_oma))

eu_omp <- lm(LROR_MP~layer, eu)
summary(aov(eu_omp))
TukeyHSD(aov(eu_omp))

eu_3ma <- lm(LR3_MA~layer, eu)
summary(aov(eu_3ma))
TukeyHSD(aov(eu_3ma))

eu_3mp <- lm(LR3_MP~layer, eu)
summary(aov(eu_3mp))
TukeyHSD(aov(eu_3mp))

eu_cu <- lm(culture~layer, eu)
summary(aov(eu_cu))
TukeyHSD(aov(eu_cu))

ba<-read.delim("basidiotxt", header=TRUE, sep="\t")
ba_oma <- lm(LROR_MA~layer, ba)
summary(aov(ba_oma))
TukeyHSD(aov(ba_oma))

ba_omp <- lm(LROR_MP~layer, ba)
summary(aov(ba_omp))
TukeyHSD(aov(ba_omp))

ba_3ma <- lm(LR3_MA~layer, ba)
summary(aov(ba_3ma))
TukeyHSD(aov(ba_3ma))

ba_3mp <- lm(LR3_MP~layer, ba)
summary(aov(ba_3mp))
TukeyHSD(aov(ba_3mp))

culture <- lm(culture~layer, ba)
summary(aov(culture))
TukeyHSD(aov(culture))

ea<-read.delim("early.txt", header=TRUE, sep="\t")
ea_oma <- lm(LROR_MA~layer,ea)
summary(aov(ea_oma))
TukeyHSD(aov(ea_oma))

ea_omp <- lm(LROR_MP~layer,ea)
summary(aov(ea_omp))
TukeyHSD(aov(ea_omp))

ea_3ma <- lm(LR3_MA~layer, ea)
summary(aov(ea_3ma))
TukeyHSD(aov(ea_3ma))

ea_3mp <- lm(LR3_MP~layer, ea)
summary(aov(ea_3mp))
TukeyHSD(aov(ea_3mp))

ea_cul<-lm(culture~layer,ea)
summary(aov(ea_cul))
TukeyHSD(aov(ea_cul))


glo<-read.delim("glo.txt", header=TRUE, sep="\t")
glo_oma<-lm(LROR_MA~layer,glo)
summary(aov(glo_oma))
TukeyHSD(aov(glo_oma))

glo_omp<-lm(LROR_MP~layer,glo)
summary(aov(glo_omp))
TukeyHSD(aov(glo_omp))

glo_3ma<-lm(LR3_MA~layer,glo)
summary(aov(glo_3ma))
TukeyHSD(aov(glo_3ma))

glo_3mp<-lm(LR3_MP~layer,glo)
summary(aov(glo_3mp))
TukeyHSD(aov(glo_3mp))


ch<-read.delim("chtri.txt", header=TRUE, sep="\t")
ch_3ma <- lm(LR3_MA~layer,ch)
summary(aov(ch_3ma))
TukeyHSD(aov(ch_3ma))

ch_3mp<-lm(LR3_MP~layer,ch)
summary(aov(ch_3mp))
TukeyHSD(aov(ch_3mp))


en<-read.delim("ento.txt", header=TRUE, sep="\t")
en_oma <- lm(LROR_MA~layer,en)
summary(aov(en_oma))
TukeyHSD(aov(en_oma))


en_omp <- lm(LROR_MP~layer,en)
summary(aov(en_omp))
TukeyHSD(aov(en_omp))


