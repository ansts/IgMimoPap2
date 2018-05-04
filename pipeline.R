
require(stringi)
require(matrixStats)
require(pepStat)
require(limma)

require(sva)
#require(energy)
require(mixtools)
require(reshape2)
require(d3heatmap)
require(Rtsne)
cpl=colorRampPalette(c("#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))

#
# Collect the data from the .gpr files in a folder /gpr containing also a key.csv file descibing the slides
#

alldata=blogpr("gpr\\")
Coor=alldata[[1]] # 1-chip/2-zone/3-diag/4-channel/5-row/6-column
iD=alldata[[2]]   # ID of spots
Res=alldata[[3]]  # Results - locally normalized data (the background is reconstituted using the duplicate diff. and subtracted using normexp)
rownames(Res)=iD  # 
Fl=alldata[[4]]   # Flags data
Nam=alldata[[5]]  # Names data (not used)

coix=as.data.frame(t(sapply(alldata[[6]], function(x){as.character(x)})), stringsAsFactors=FALSE) # 1.Chip;2.Zone;3.Channel;4.Patient ID _ Diagnosis;5.Positive Column; 6.Negative Column
FLdf=as.data.frame(Fl[1:ncol(Res)])
FL=rowSums(FLdf)==0
Cop=Coor[as.double(coix[,5])]

WD=Res[FL,]
WDa=aggregate.data.frame(WD, by=list(rownames(WD)), FUN=mean)

WDam=data.frame(log10(WDa[,2:ncol(WDa)]), stringsAsFactors = FALSE, row.names = WDa[,1])

WDcolmeans=colMeans(as.matrix(WDam))
plot(colMaxs(as.matrix(WDam)),colMeans(as.matrix(WDam)),cex=0.3)
text(colMaxs(as.matrix(WDam)),colMeans(as.matrix(WDam)),labels = colnames(WDam))

# fix column names
cnD=colnames(WDam)
cnD=gsub("\\.","",cnD)
cnD[cnD=="RKD1_C"]="RKD2_C"
cnD[cnD=="RKD_C"]="RKD1_C"
colnames(WDam)=cnD

ccnD=stri_split(cnD, regex = "(?=[\\d_])")
ccnD1=lengths(ccnD)==2
tcnD=table(sapply(ccnD[!ccnD1],function(x){x[1]}))

# discard low signal arrays
final=sapply(seq_along(ccnD), function(i){if (ccnD1[i]) return(TRUE) else if (ccnD[[i]][2]==tcnD[ccnD[[i]][1]]) return(TRUE) else return(FALSE)})
final=final&(WDcolmeans>=2)
WDam=WDam[,final]
cnD=cnD[final]
ch=key$ch[1:ncol(Res)]
ch=ch[final]
dg=sapply(Cop,function(x){x[[3]]})
dg=dg[final]

# normalize for amino acid composition dependent binding (e.g. non-specific stickiness of charged amino acids)
D=pepnorm(WDam)

# fuse data on peppos from previous format chips with peptotest chips
load("D7")
idall=intersect(rownames(D), rownames(D7))
Dall=cbind(D[idall,], D7[idall,])
dgall=as.factor(c(as.character(dg),c("GBM", "ML","GBM", "C","C", "GBM", "ML","ML","ML","GBM","C")))
dgalless=dgall
dgalless[dgalless=="GBMo"]="GBM"
dgalless=factor(dgalless, levels = unique(dgalless))

# label batches: red channel, green channel, old arrays (which are green too)
batches=as.factor(c(as.character(ch),rep("O",11)))
Dnall=normalizeCyclicLoess(Dall, method="affy", span=.15, iterations=4)
Dncall=ComBat(Dnall,batches, model.matrix(~dgalless))
blncDnall=!(dgalless %in% c("MB", "G")) & !((batches == "R") & (dgalless=="GBM") & !(colnames(Dnall) %in% the5))
Dncll=ComBat(Dnall[,blncDnall],batches[blncDnall])
dgis=dgi[blncDnall]
batchs=batches[blncDnall]
bgis=as.integer(batchs)
dgll=factor(dgall[blncDnall])
dgllC=dgll=="C"
dgllGBM=dgll=="GBM"
dgllML=dgll=="ML"
peps=rownames(Dnall)

plotMDS(Dnall, col=dgi)
plotMDS(Dnall, col=as.integer(batches))
plotMDS(Dncall, col=dgi)
plotMDS(Dncall, col=as.integer(batches))
plotMDS(Dncll, col=dgis)
plotMDS(Dncll, col=bgis)

Dnll3=Dnall[,batches=="R" & dgall %in% c("C","GBM","ML", "MB")]
meanall=rowMeans(Dncall)
meanall3=rowMeans(Dnll3)
sdall=rowSds(Dncall)
sdall3=rowSds(Dnll3)


#find background spots

tsDncall=Rtsne(Dncall, initial_dims = 10, perplexity = 5)
plot(tsDncall$Y, pch=16, col=cpl(586)[rank(meanall)], xlab="F1", ylab="F2", main="Rtsne Plot of Data")
tsnDncll=Rtsne(Dncll, initial_dims = 10, perplexity = 5)
plot(tsnDncll$Y, pch=16, col=cpl(586)[rank(meanall)], xlab="F1", ylab="F2", main="Rtsne Plot of Data")
tsnDnll3=Rtsne(Dnll3, initial_dims = 3, perplexity = 20)
plot(tsnDnll3$Y, pch=16, col=cpl(586)[rank(meanall3)], xlab="F1", ylab="F2", main="Rtsne Plot of Data")


bckg=tsDncall$Y[,1]<(-25) # the exact position of the clusters is random so this should be adjusted accordingly
mnbg=mean(Dncall[bckg,])
bckfg=tsnDncll$Y[,2]<(-30) # the exact position of the clusters is random so this should be adjusted accordingly
mnfbg=mean(Dncll[bckfg,])
Dncc=Dncll-mnfbg
bck3g=tsnDnll3$Y[,2]>(33)
mn3bg=mean(Dnll3[bck3g,])
Dn3=Dnll3-mn3bg
dgl3=as.factor(dgalless[batches=="R" & dgall %in% c("C","GBM","ML","MB")])
dg3=as.character(dgl3)
dg3=as.factor(dg3)
dg3GBM=dg3 == "GBM"

#limma

vsC05_=limfit(dgalless,Dncall-mnbg, cmf=4, p=0.05)
vbC01_=limfit(dgll,Dncc, cmf=4, p=0.01)
vbc=unique(unlist(vbC01_))
v3c01_=limfit(dg3,Dn3, cmf=4, p=0.01)
v3c=unique(unlist(v3c01_))
vbdiff=limfit(dgll,Dncc, cmf=1, p=0.05)

plot(meanall, sdall, col=(idall %in% vbC05_all)*1+1)
plotMDS(Dncll[vbC05_all,], col=dgi)


d3New=(dg3GBM&(!(colnames(Dn3) %in% the5)))|(dg3=="MB")
d3Ntrue=dg3GBM[d3New]

Dncllnms=colnames(Dncll)
Dncc=Dncll-mnfbg


DLOOX3GBM_=lapply(seq_along(dgll),function(p){
  scan3crit(Dncc[vbc,-p], dgllGBM[-p])
})

DLOOX3C_=lapply(seq_along(dgll),function(p){
  scan3crit(Dncc[vbc,-p], dgllC[-p])
})

DLOOX3ML_=list()
for (p in 1:28) {
  x=scan3crit(Dncc[vbc,-p], dgllML[-p])
  DLOOX3ML_[[p]]=x
  flnm=paste(c("DLOOX3ML_",p), collapse="")
  save(DLOOX3ML_, file=flnm)
  print(paste(c("Case",p),collapse=" "))
}


DLOOX3GBMe_=list()
for (p in 1:28) {
  x=scan3crit(Dncc[vbc,-p], dgllGBM[-p], d=2)
  DLOOX3GBMe_[[p]]=x
  flnm=paste(c("DLOOX3GBMe_",p), collapse="")
  save(DLOOX3GBMe_, file=flnm)
  print(paste(c("Case",p),collapse=" "))
}
DLOOX3Ce_=list()
for (p in 1:28) {
  x=scan3crit(Dncc[vbc,-p], dgllC[-p], d=2)
  DLOOX3Ce_[[p]]=x
  flnm=paste(c("DLOOX3Ce_",p), collapse="")
  save(DLOOX3Ce_, file=flnm)
  print(paste(c("Case",p),collapse=" "))
}
DLOOX3MLe_=list()
for (p in 1:28) {
  x=scan3crit(Dncc[vbc,-p], dgllML[-p], d=2)
  DLOOX3MLe_[[p]]=x
  flnm=paste(c("DLOOX3MLe_",p), collapse="")
  save(DLOOX3MLe_, file=flnm)
  print(paste(c("Case",p),collapse=" "))
}



daloox3GBM_=lapply(DLOOX3GBM_, function(i){rownames(i[1][[1]])})
names(daloox3GBM_)=Dncllnms
daloox3GBM_t=table(unlist(daloox3GBM_))
daloox3GBM_tco=sapply(0:27, function(i){x=names(daloox3GBM_t)[daloox3GBM_t>i]; x[order(daloox3GBM_t[x])]})
daloox3GBM_tco3=lapply(daloox3GBM_tco,function(x){x[x %in% v3c]})

daloox3C_=lapply(DLOOX3C_, function(i){rownames(i[1][[1]])})
names(daloox3C_)=Dncllnms
daloox3C_t=table(unlist(daloox3C_))
daloox3C_tco=sapply(0:27, function(i){x=names(daloox3C_t)[daloox3C_t>i]; x[order(daloox3C_t[x])]})
daloox3C_tco3=lapply(daloox3C_tco,function(x){x[x %in% v3c]})

daloox3ML_=lapply(DLOOX3ML_, function(i){rownames(i[1][[1]])})
names(daloox3ML_)=Dncllnms
daloox3ML_t=table(unlist(daloox3ML_))
daloox3ML_tco=sapply(0:27, function(i){x=names(daloox3ML_t)[daloox3ML_t>i]; x[order(daloox3ML_t[x])]})
daloox3ML_tco3=lapply(daloox3ML_tco,function(x){x[x %in% v3c]})


daloox3GBMe_=lapply(DLOOX3GBMe_, function(i){rownames(i[1][[1]])})
names(daloox3GBMe_)=Dncllnms
daloox3GBMe_t=table(unlist(daloox3GBMe_))
daloox3GBMe_tco=sapply(0:27, function(i){x=names(daloox3GBMe_t)[daloox3GBMe_t>i]; x[order(daloox3GBMe_t[x])]})
daloox3GBMe_tco3=lapply(daloox3GBMe_tco,function(x){x[x %in% v3c]})

daloox3Ce_=lapply(DLOOX3Ce_, function(i){rownames(i[1][[1]])})
names(daloox3Ce_)=Dncllnms
daloox3Ce_t=table(unlist(daloox3Ce_))
daloox3Ce_tco=sapply(0:27, function(i){x=names(daloox3Ce_t)[daloox3Ce_t>i]; x[order(daloox3Ce_t[x])]})
daloox3Ce_tco3=lapply(daloox3Ce_tco,function(x){x[x %in% v3c]})

daloox3MLe_=lapply(DLOOX3MLe_, function(i){rownames(i[1][[1]])})
names(daloox3MLe_)=Dncllnms
daloox3MLe_t=table(unlist(daloox3MLe_))
daloox3MLe_tco=sapply(0:27, function(i){x=names(daloox3MLe_t)[daloox3MLe_t>i]; x[order(daloox3MLe_t[x])]})
daloox3MLe_tco3=lapply(daloox3MLe_tco,function(x){x[x %in% v3c]})

plotMDS(Dncll[names(daloox3GBM_t)[daloox3GBM_t>23],], col=as.integer(dgll))
for (i in 1:28) plot(cmdscale(fraD(t(Dncll[names(daloox3GBM_t[daloox3GBM_t==i]),]),.25)), main=i, pch=16, cex=1.2, col=as.integer(dgll))

daltco_hc=hclust(dist((Dn3[daloox3GBM_tco[[19]],])))
daltco_dg=as.dendrogram(daltco_hc)
d3heatmap(t(Dncc[daloox3GBM_tco[[19]],]), Colv = daltco_dg, distfun = fraD, hclustfun = hclwrd, col=cpl(128), dendrogram = "row")
d3heatmap(t(Dncc[daloox3GBM_tco[[19]],]),distfun = fraD, hclustfun = hclwrd, col=cpl(128))

# plots of MDS D1xD2 of the data for the features common for 1 to 28 of the LOOV bootstrap 
# derived sets (e.g. "Set 14" means the features comon for 14 or more bootstrap sets,
# 1 means all possible unique features found in the bootstrap, 2 - the non-specific features 
# and 28 - those that are found in each and every bootstrap set)
for (i in 1:28) plotMDS(Dncc[daloox3Ce_tco[[i]],], col=as.integer(dgll), main=i)
for (i in 1:28) plotMDS(Dncc[daloox3GBMe_tco[[i]],], col=as.integer(dgll), main=i)
for (i in 1:28) plotMDS(Dncc[daloox3MLe_tco[[i]],], col=as.integer(dgll), main=i)



for (i in 1:28) {
  plot(cmdscale(fraD(t(Dn3[daloox3GBM_tco[[i]],]))), cex=0, main=i)
  text(cmdscale(fraD(t(Dn3[daloox3GBM_tco[[i]],]))), col=as.integer(dg3), labels = colnames(Dn3))
}

for (i in 1:28) {
  plot(cmdscale(fraD(t(Dn3[daloox3GBMe_tco[[i]],]),2)), cex=0, main=i)
  text(cmdscale(fraD(t(Dn3[daloox3GBMe_tco[[i]],]),2)), col=as.integer(dg3), labels = colnames(Dn3))
}

for (i in 1:28) {
  plot(cmdscale(fraD(t(Dn3[daloox3GBM_tco3[[i]],]))), cex=0, main=i)
  text(cmdscale(fraD(t(Dn3[daloox3GBM_tco3[[i]],]))), col=as.integer(dg3), labels = colnames(Dn3))
}

for (i in 1:28) {
  plotMDS(Dn3[daloox3GBM_tco3[[i]],], col=as.integer(dg3), main=paste(c("E",i)))
}



################### Plot of best separation of batch "R" using the generalizing feature set common for 19 
# or more bootstrap sets derived from the batch compensated data set. Note the proper classification of the 
# labeled cases - those were not included in the batch compensated set

plot(cmdscale(fraD(t(Dn3[daloox3GBM_tco[[19]],]),.25)), col=as.integer(dg3),cex=1.5, pch=16, main="Fractional Distance")
par(new=TRUE)
plot(cmdscale(fraD(t(Dn3[daloox3GBM_tco[[19]],]),.25)), col=as.integer(!(colnames(Dn3) %in% Dncllnms)),cex=2.5)
par(new=FALSE)

plot(cmdscale(fraD(t(Dn3[daloox3GBM_tco[[19]],]),2)), col=as.integer(dg3),cex=1.5, pch=16, main="Euclidean Distance")
par(new=TRUE)
plot(cmdscale(fraD(t(Dn3[daloox3GBM_tco[[19]],]),2)), col=as.integer(!(colnames(Dn3) %in% Dncllnms)),cex=2.5)
par(new=FALSE)


daltco_hc=hclust(dist((Dn3[daloox3GBM_tco[[19]],])))
daltco_dg=as.dendrogram(daltco_hc)
##############
d3heatmap(t(Dncc[daloox3GBM_tco[[19]],]), Colv = daltco_dg, distfun = fraD, hclustfun = hclwrd, col=cpl(128), dendrogram = "row")
##############
d3heatmap(t(Dncc[daloox3GBM_tco[[19]],]),distfun = fraD, hclustfun = hclwrd, col=cpl(128), dendrogram = "row")
d3mds=cmdscale(fraD(t(Dn3[daloox3GBM_tco[[19]],])), k=2)
d3MDS_2cor=apply(Dn3[daloox3GBM_tco[[19]],],1,function(x){cor(x,d3mds[,2])})
d3MDS_2cor[order(d3MDS_2cor)]
d3MDS2=names(d3MDS_2cor[abs(d3MDS_2cor)>0.3])
dalmds2_hc=hclust(dist((Dn3[d3MDS2,])))
dalmds2_dg=as.dendrogram(dalmds2_hc)
d3heatmap(t(Dn3[d3MDS2,]), distfun = fraD, Colv =, hclustfun = hclwrd, col=cpl(128))

# Original script for MCCtr_val.pdf - see figmaker.R

MCCesvmt=unlist(MCCesvmt)
MCCesvmv=unlist(MCCesvmv)
names(MCCesvmt)=1:28
names(MCCesvmv)=1:28
op=par(mfrow=c(2,1), omi=c(0.85,0.85,0,0), mai=c(.5,0,0.75,0))
barplot(MCCesvmt, ylim = c(-0.5,1), main="Training")
barplot(MCCesvmv, ylim = c(-0.5,1), main="Validation")
title(xlab = "Commonality", ylab="MCC", outer=TRUE)
par(op)


## svm
MCCesvmv=sapply(1:28,function(i){
  y=cmdscale(fraD(t(Dn3[daloox3GBMe_tco[[i]],]),2))
  x=svm(y[!d3New,], as.factor(dg3GBM[!d3New]), gamma = 0.1, cost=1000)
  t=table(d3Ntrue,predict(x, newdata=y[d3New,]))
  TP=t[4]; TN=t[1];FP=t[2];FN=t[3]
  z=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (is.na(z)) 0 else z
})
print(c(max(MCCesvmv),which.max(MCCesvmv)))
MCCesvmt=sapply(1:28,function(i){
  y=cmdscale(fraD(t(Dn3[daloox3GBMe_tco[[i]],]),2))
  x=svm(y, as.factor(dg3GBM), gamma = 0.1, cost=1000)
  t=table(dg3GBM,predict(x, newdata=y))
  TP=t[4]; TN=t[1];FP=t[2];FN=t[3]
  (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
})
print(c(max(MCCesvmt),which.max(MCCesvmt)))
xp=unlist(MCCesvmt)+runif(28,-0.01,0.01)
yp=unlist(MCCesvmv)+runif(28,-0.01,0.01)
plot(xp,yp,cex=0, xlim=c(0,1),ylim=c(-1,1))
lines(xp,yp,cex=0, xlim=c(0,1),ylim=c(-1,1))
text(xp,yp,1:28, cex=0.8, xlim=c(0,1),ylim=c(-1,1))

