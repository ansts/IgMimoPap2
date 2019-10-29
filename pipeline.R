
require(stringi)
require(matrixStats)
require(pepStat)
require(limma)
require(prallel)
require(sva)
require(mixtools)
require(reshape2)
require(d3heatmap)
require(Rtsne)
require(umap)
require(clusterCrit)
require(rgl)
require(Rfast)
require(e1071)
require(plot3D)

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

# selecting 5 GBM cases from batch "R" to include to balance the batch/group distribution

G5m=combn(colnames(Dall)[dgall=="GBM" & batches=="R"],5)
GBM5means=apply(G5m,2,function(s){rowMeans(Dnall[,s])})
x=rowMeans(Dnall[,colnames(Dnall)[dgalless=="GBM" & batches=="R"]])
GBM5means=cbind(x,GBM5means)
colnames(GBM5means)=c(0,1:126)
GBM5CVs=apply(G5m,2,function(s){rowSds(Dnall[,s])})
x=rowSds(Dnall[,colnames(Dnall)[dgalless=="GBM" & batches=="R"]])
GBM5CVs=cbind(x,GBM5CVs)
colnames(GBM5CVs)=c(0,1:126)
GBM5CVs=GBM5CVs/GBM5means

# find the distance from the mean/CV of the rows of all R GBMs to the mean/var of the rows of each combination of 4

dG5v=as.matrix(dist(t(GBM5CVs)))[2:127,1]
dG5m=as.matrix(dist(t(GBM5means)))[2:127,1]

# and the distance from the origin to the point defined by the ranks of mean/var

dG5mv=sqrt(rank(dG5m)^2+rank(dG5v)^2)

# the 5 cases are
the5=G5m[,which.min(dG5mv)]

# ComBat batch compensation

blncDnall=!(dgalless %in% c("MB", "G")) & !((batches == "R") & (dgalless=="GBM") & !(colnames(Dnall) %in% the5))
batchs=batches[blncDnall]
Dncll=ComBat(Dnall[,blncDnall],batchs)
dgis=dgi[blncDnall]
bgis=as.integer(batchs)
dgll=factor(dgall[blncDnall])
dgllGBM=dgll=="GBM"
meanall=rowMeans(Dncll)
sdall=rowSds(Dncll)
peps=rownames(Dnall)

Dnll3=Dnall[,batches=="R"]
meanall3=rowMeans(Dnll3)
sdall3=rowSds(Dnll3)

#find background spots
# for batch 3 alone
tsnDnll3=Rtsne(Dnll3, initial_dims = 3, perplexity = 20)
plot(tsnDnll3$Y, pch=16, col=cpl(586)[rank(meanall3)], xlab="F1", ylab="F2", main="Rtsne Plot of Data")
bck3g=tsnDnll3$Y[,2]<(-35)
mn3bg=mean(Dnll3[bck3g,])
Dn3=Dnll3-mn3bg

dgl3=as.factor(dgalless[batches=="R"])
dg3=as.character(dgl3)
dg3=as.factor(dg3)
dg3GBM=dg3 == "GBM"
d3New=(dg3GBM&(!(colnames(Dn3) %in% the5)))|(dg3=="MB")
d3Ntrue=dg3GBM[d3New]
d3old=!d3New&dg3!="G"

# ... and for the batch compensated sample pool
tsnDncll=Rtsne(Dncll, initial_dims = 3, perplexity = 20, theta=0)
plot(tsnDncll$Y, pch=16, col=cpl(586)[rank(meanall)], xlab="F1", ylab="F2", main="Rtsne Plot of Data")
bckfg=tsnDncll$Y[,2]>(63) # the exact position of the clusters is random so this should be adjusted accordingly
plot(tsnDncll$Y, pch=16, col=bckfg+1, xlab="F1", ylab="F2", main="Rtsne Plot of Data")

mnfbg=mean(Dncll[bckfg,])
Dncllnms=colnames(Dncll)
Dncc=Dncll-mnfbg



# limma
# for the comBat-ed set

vbc_=lapply(1:28,function(x){
  f=rep(0,28)
  f[x]=1
  f=factor(f)
  r=limfit_(f,Dncc, p=0.05)
  return(r[[1]])
})

vbc0=unique(unlist(vbc_)) # the mimotopes with significant reactivity with at least one pateint's IgM
lv=length(vbc0)

vbcranks=t(sapply(vbc0,function(i){     # convert the table of mimotope reactivities (rows - the mimotopes in descending order columns - patients)
  sapply(vbc_, function(j){             # ordered by the limma algorithm by degree of expression to ranks for each column (patient)
    if (i %in% j) which(j==i) else length(j)+1
  })
}))
rownames(vbcranks)=vbc0
vbcrmd=rowMedians(vbcranks)
vbcranks=vbcranks[order(vbcrmd),]       # the rows of the table of ranks are ordered by the median rank for each mimotope reactivity

bf=rownames(vbcranks)[1:294]  
FSboots=lapply(1:28, function(i){                                    # bootstrap based selection of profiles using a leave one out scheme for the patients and starting with the upper half of the mimotopes 
  proct=proc.time()                                                  # ordered by median rank which organizes the profiles around highly expressed reactivities,
  xb=profileSearch(Dncc[,-i],bf,dgllGBM[-i])                          # the profilesearch algorithm uses them only for the initial recursive elimination search
  x=xb[[2]]                                                          # further adding features from the whole set in the subsequent forward search. The two steps are repeated 
  plotMDS(Dncc[x,], col=as.integer(dgll), main=rownames(Dncc)[i])    # 5 times which yields slightly different sets due to the a stochastic factor in the criterion, 
  fnm=paste("FSboots_",sample(10000,1),".txt", sep="", collapse="")  # then the consensus features from the 5 runs are selected and a final recursive elimination step is performed on them. 
  write(x, file=fnm)
  beep(3)
  print(proc.time()-proct)
  return(xb)                             # xb[[1]] is a list of the consensus profiles after the first and the second step
})                                       # xb[[2]] is the final profile after the recursive elimination applied on the consensus of the second step

# summarize
FSbootsi=sapply(FSboots, function(l){l[[2]]})
FSboot=table(unlist(FSbootsi))
FSbootcom=sapply(0:27, function(i){x=names(FSboot)[FSboot>i]; x[order(FSboot[x])]}) # profiles ordered by the number of bootstrap sets the features 
                                         # are common for starting from none (all features are concatenated) to 
                                         # 28 - the features found in all bootstrap sets (from a leave one out scheme relative to the patients)
                                         # the 28th feature set contains only one peptide common for all sets and is ommitted from now on

for (i in 1:27) plotMDS(Dncc[FSbootcom[[i]],], col=as.integer(dgll), main=i)

for (i in 1:27) {
  plot(cmdscale(dist(t(Dn3[FSbootcom[[i]],]))), cex=0, main=i)
  text(cmdscale(dist(t(Dn3[FSbootcom[[i]],]))), col=as.integer(dg3), labels = colnames(Dn3))
}

MCCscrv=t(sapply(1:1000, function(j){     # A thousand iterations of scrambling the R batch data matrix Dn3 (before batch compensation) by row
  scrD3=apply(Dn3,2,sample)               # and training svm models on the patients common with the comBat-ed set Dncc, followed by predicting the     
  rownames(scrD3)=rownames(Dn3)           # the rest of the patients. At each step 27 models are trained using one of the 27 feature sets selected 
  colnames(scrD3)=colnames(Dn3)           # in the 'summarize' step of the bootstrap feature selection algorithm described above.
  m=sapply(1:27,function(i){              # the Mathiew's Correlation Coefficient is calculated to assessthe quality of each model 
    y=cmdscale(dist(t(scrD3[FSbootcom[[i]],])))   #and is saved for further processing i MCCscrv, the atter having dimensions 1000,27
    x=svm(y[d3old,], as.factor(dg3GBM[d3old]), gamma = 0.003, cost=1000)
    t=table(d3Ntrue,predict(x, newdata=y[d3New,]))
    TP=t[4]; TN=t[1];FP=t[2];FN=t[3]
    z=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if (is.nan(z)) return(0) else return(z)
  })
  return(m)
}))
dim(MCCscrv)
xvpru=apply(MCCscrv,2,quantile,0.95)    # The 0.05 and the 0.95 percentile's of MCC values of the 1000 models for each column of MCCscrv are taken as a confidence interval
xvprl=apply(MCCscrv,2,quantile,0.05)    # for MCC of the actual models that will use the feature set corresponding to that column.

MCCesvmv=sapply(1:27,function(i){       #  MCC of svm models trained and tested on the same patient sets as above (for MCCscrv) and the 27 feature sets 
    y=cmdscale(dist(t(Dn3[FSbootcom[[i]],])))    # but using the actual data matrix Dn3
    x=svm(y[d3old,], as.factor(dg3GBM[d3old]), gamma = 0.003, cost=1000)
    t=table(d3Ntrue,predict(x, newdata=y[d3New,]))
    TP=t[4]; TN=t[1];FP=t[2];FN=t[3]
    z=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if (is.nan(z)) return(0) else return(z)
})

MCCesvmt=sapply(1:27,function(i){       #  MCC of svm models trained on the same patient sets as above and the 27 feature sets but tested back on the training set
  y=cmdscale(dist(t(Dn3[FSbootcom[[i]],])))
  x=svm(y[d3old,], as.factor(dg3GBM[d3old]), gamma = 0.003, cost=1000)
  t=table(dg3GBM[d3old],predict(x, newdata=y[d3old,]))
  TP=t[4]; TN=t[1];FP=t[2];FN=t[3]
  z=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (is.nan(z)) return(0) else return(z)
})
names(MCCesvmt)=1:27
names(MCCesvmv)=1:27

# Figures

pdf(file="MCC_models.pdf", width=6, height=9.7)
op=par(mfrow=c(2,1), omi=c(1,0.85,0,0), mai=c(0.5,0,1,0))
barplot(t(MCCesvmv), ylim=c(-1,1), xlim=c(0,28), width=0.9, space=0.11111,  col=1, las=3)
title(main="                                                          Validation", line=0)
par(new=T)
barplot(xvpru,ylim=c(-1,1), xlim=c(0,28), width=1, space=0, border=rgb(0.8,0.8,0.8,1), col=rgb(0.5,0.5,0.5,0.5))
par(new=T)
barplot(xvprl,ylim=c(-1,1), xlim=c(0,28), width=1, space=0, border=rgb(0.8,0.8,0.8,1), col=rgb(0.5,0.5,0.5,0.5))
par(new=F)
barplot(t(MCCesvmt), ylim=c(-1,1), xlim=c(0,28), width=0.9, space=0.11111,  col=1, las=3)
title(main="                                                          Training", line=0)
par(new=T)
barplot(xvpru,ylim=c(-1,1), xlim=c(0,28), width=1, space=0, border=rgb(0.8,0.8,0.8,0.5), col=rgb(0.5,0.5,0.5,0.5))
par(new=T)
barplot(xvprl,ylim=c(-1,1), xlim=c(0,28), width=1, space=0, border=rgb(0.8,0.8,0.8,0.5), col=rgb(0.5,0.5,0.5,0.5))
title(xlab = "Commonality", ylab="MCC", outer=TRUE)
par(op)
dev.off()



y=cmdscale(dist(t(Dn3[FSbootcom[[15]],])))
r=apply(y,2,range)
x=svm(y[d3old,], as.factor(dg3GBM[d3old]), probability=T, gamma = 0.003, cost=1000)
pry=sapply(seq(r[,1][1], r[,1][2], length.out = 10), function(x1){         # preparing the smooth color coded probability field for the svm plot
    sapply(seq(r[,2][1], r[,2][2], length.out = 10), function(x2){
      attr(predict(x, newdata=t(c(x1,x2)), probability = TRUE), "probabilities")[1]
  })
})

pdf(file="SVMpred.pdf", width=7, height=6)
dg3x=stri_extract_all(colnames(Dn3), regex="(?<=_)\\w+")
coldg3=as.integer(dg3)
coldg3[coldg3==3]=2
coldg3[coldg3==5]=3
d3=d3New|d3old
op=par(mai=c(1,1,1,1))
image2D(t(pry),resfac = 10, col=cm.colors(1000), xaxt="n", yaxt="n", xlab="", ylab="", clab="GBM Probability")
op=par(new=TRUE,mai=c(1,1,1,1.575))
par(new=T)
plot(y[,], cex=0, xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2")
points(y[d3New,], cex=5, xlim=r[,1], ylim=r[,2], pch=1, xlim=r[,1], ylim=r[,2])
text(y[d3,], labels=dg3x[d3], col=coldg3[d3], xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2")
dev.off()
