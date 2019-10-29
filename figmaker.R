
    
pdf(file="NseqperDiagE.pdf", width=4, height=4)
plot(c(1,2,3),c(lengths(dlooe3C),lengths(dlooe3GBM),lengths(dlooe3ML),lengths(DLOOE3C),lengths(DLOOE3GBM),lengths(DLOOE3ML)), 
     cex=0.75, xlab = "", ylab="N of Sequences", pch=c(rep(16,32),rep(1,32)),xlim=c(0,4), xaxt="n", cex.axis=0.75, cex.lab=0.75)
axis(1,at=c(1,2,3),labels = c("Control","GBM","ML"), cex.axis=0.75)
legend("topright", legend=c("P=0.25","Euclidean" ), bty="n",pch=c(16,1),cex=0.5)
dev.off()

pdf(file="MCCpred.pdf", width=7, height=7)
par(mfrow=c(2,2), mai=c(0.2,0.2,0.2,0.2), cex.axis=0.75)
par(omi=c(1,1,1,1))
clrs=c(1,3,6)
lapply(1:3, function(i){plot(nms[[i]],MCCs$Eu[[i]], col=clrs[i], pch=1+(MCCs$Eu[[i]]>MCCsig)*15, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", cex.main=1, main="Model Minkowski Distance p=2")
    lines(nms[[i]],MCCs$Eu[[i]], col=clrs[i], lwd=1, xlim=c(0,1), ylim=c(0,1));grid(nx=NULL, lwd=0.3);par(new=TRUE)})
par(new=FALSE)
lapply(4:6, function(i){plot(nms[[i-3]],MCCs$Eu[[i]], col=clrs[i-3], pch=1+(MCCs$Eu[[i]]>MCCsig)*15, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", cex.main=1, main="Model Minkowski Distance p=0.25")
    lines(nms[[i-3]],MCCs$Eu[[i]], col=clrs[i-3], lwd=1, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="");grid(nx=NULL, lwd=0.3);par(new=TRUE)})
par(new=FALSE)
lapply(1:3, function(i){plot(nms[[i]],MCCs$Fra[[i]], col=clrs[i], pch=1+(MCCs$Fra[[i]]>MCCsig)*15, xlim=c(0,1), ylim=c(0,1), xlab="",ylab="")
    lines(nms[[i]],MCCs$Fra[[i]], col=clrs[i], lwd=1, xlim=c(0,1), ylim=c(0,1));grid(nx=NULL, lwd=0.3);par(new=TRUE)})
par(new=FALSE)
lapply(4:6, function(i){plot(nms[[i-3]],MCCs$Fra[[i]], col=clrs[i-3], pch=1+(MCCs$Fra[[i]]>MCCsig)*15,xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    lines(nms[[i-3]],MCCs$Fra[[i]], col=clrs[i-3], lwd=1, xlim=c(0,1), ylim=c(0,1, xlab="", ylab=""));grid(nx=NULL, lwd=0.3);par(new=TRUE)})
par(new=FALSE)
mtext(text=c("Feature Selection p=2","Feature Selection p=0.25"), side=4, outer=TRUE, srt=90, at=c(0.75,0.25))
title(ylab = "MCC", xlab="Feature Subset Commonality",outer=TRUE)
dev.off()

dlooC=unique(unlist(c(dlooe3C,dloox3C)))
dlooGBM=unique(unlist(c(dlooe3GBM,dloox3GBM)))
dlooML=unique(unlist(c(dlooe3ML,dloox3ML)))

pdf(file="CommonalitiesC.pdf", width=4.5, height=4)
dlooc=sapply(1:7, function(i){
c(length(setdiff(come3C[[i]], comx3C[[i]]))/463,
length(intersect(come3C[[i]], comx3C[[i]]))/463,
length(setdiff(comx3C[[i]], come3C[[i]]))/463)
})
dlooxcunq=rownames(xC[xC==1])
dlooecunq=rownames(eC[eC==1])
x=c(length(setdiff(dlooecunq, dlooxcunq))/463,
length(intersect(dlooecunq, dlooxcunq))/463,
length(setdiff(dlooxcunq, dlooecunq))/463)
dlooc=cbind(dlooc,x)
colnames(dlooc)=c(8:1)/8
rownames(dlooc)=c("P=2", "Both","P=0.25")
par0=par(mai=c(1,1,0.5,1))
palette(gray(seq(0,.9,len=3)))
barplot(dlooc1, beside = FALSE, col=(1:3), ylim=c(0,1.25),cex.axis = 0.5, las=2, cex=0.5, main = "Control")
legend("top", legend = rownames(dlooc1), fill=(1:3), bty = "n", cex=0.5)
title(xlab = "Proportion of Sets the Subset is Common for", ylab = "Proportion of Features Selected", cex.lab=0.5)
palette("default")
par(par0)
dev.off()

pdf(file="CommonalitiesGBM.pdf", width=4, height=4)
dlooGBM=sapply(1:14, function(i){
  c(length(setdiff(come3GBM[[i]], comx3GBM[[i]]))/463,
    length(intersect(come3GBM[[i]], comx3GBM[[i]]))/463,
    length(setdiff(comx3GBM[[i]], come3GBM[[i]]))/463)
})
dlooxgbmunq=rownames(xGBM[xGBM==1])
dlooegbmunq=rownames(eGBM[eGBM==1])
x=c(length(setdiff(dlooegbmunq, dlooxgbmunq))/463,
    length(intersect(dlooegbmunq, dlooxgbmunq))/463,
    length(setdiff(dlooxgbmunq, dlooegbmunq))/463)
dlooGBM=cbind(dlooGBM,x)
colnames(dlooGBM)=round(c(15:1)/15,2)
rownames(dlooGBM)=c("P=2", "Both","P=0.25")
par0=par(mai=c(1,1,0.5,1))
palette(gray(seq(0,.9,len=3)))
barplot(dlooGBM1, beside = FALSE, col=(1:3), cex.axis = 0.5, ylim=c(0,1.2), las=2, cex=0.5, main = "GBM")
legend("top", legend = rownames(dlooc), fill=(1:3), bty = "n", cex=0.5)
title(xlab = "Proportion of Sets the Subset is Common for", ylab = "Proportion of Features Selected", cex.lab=0.5)
palette("default")
par(par0)
dev.off()

pdf(file="CommonalitiesML.pdf", width=4, height=4)
dlooML=sapply(1:8, function(i){
  c(length(setdiff(come3ML[[i]], comx3ML[[i]]))/463,
    length(intersect(come3ML[[i]], comx3ML[[i]]))/463,
    length(setdiff(comx3ML[[i]], come3ML[[i]]))/463)
  
})
dlooxMLunq=rownames(xML[xML==1])
dlooeMLunq=rownames(eML[eML==1])
x=c(length(setdiff(dlooeMLunq, dlooxMLunq))/463,
    length(intersect(dlooeMLunq, dlooxMLunq))/463,
    length(setdiff(dlooxMLunq, dlooeMLunq))/463)
dlooML=cbind(dlooML,x)
colnames(dlooML)=round(c(9:1)/9,2)
rownames(dlooML)=c("P=2", "Both","P=0.25")
par0=par(mai=c(1,1,0.5,1))
palette(gray(seq(0,.9,len=3)))
barplot(dlooML1, beside = FALSE, col=(1:3), cex.axis = 0.5, ylim=c(0,1.2), las=2, cex=0.5, main = "ML")
legend("top", legend = rownames(dlooc), fill=(1:3), bty = "n", cex=0.5)
palette("default")
title(xlab = "Proportion of Sets the Subset is Common for", ylab = "Proportion of Features Selected", cex.lab=0.5)
par(par0)
dev.off()


pdf(file="VennsGBM_df.pdf",onefile = FALSE,  width=8, height=24)
require(eulerr)
featbdg=unique(unlist(bydg05))
L=length(featbdg)
i=1:7
j=list(1:7,c(1,3,5,7,9,12,14), c(1,2,3,4,6,7,8))
X=lapply(i, function(i){list("Control"=comx3C[[j[[1]][i]]],"GBM"=comx3GBM[[j[[2]][i]]],"ML"=comx3ML[[j[[3]][i]]],"GBM - all"=dx3GBM) })

mx_dfxcomx0=lapply(X,function(x){sapply(x,function(s){featbdg %in% s}, USE.NAMES = TRUE)})
mx_dfxcomx=do.call(rbind,mx_dfxcomx0)
j0=sapply(j,function(j){j})
mxnm=c(apply(j0,1,function(j){rep(paste(c("C","GBM","ML","GBMall"), round((c(9,16,10)-j)/c(8,15,9),2),collapse = " ", sep="_"),L)}))

venn_comdf=euler(mx_dfxcomx,by=mxnm)
plot(venn_comdf, layout=c(2,4),  counts=TRUE, auto.key=TRUE, fontsize=5)
dev.off()
###############
pdf(file="metaVennsGBM_df.pdf", width=4, height=4)
x=comx3GBM
names(x)=c(round((16-1:14)/15,2))
barplot(sapply(x,function(s){c(length(setdiff(s,dx3GBM)),length(intersect(s,dx3GBM)),length(setdiff(dx3GBM,s)))}), 
        xlab="Proportion of overlapping sets", ylab="N of features", cex.names=0.75,  cex.lab=0.75, cex.axis=0.75,las=2, 
        legend.text = c("Only in bootstrap overlap","In both","Only in 'all'"), args.legend = list(x="top", cex=0.5))

dev.off()
#############
pdf(file="MDS_bydg5C.pdf",  width=6, height=6)
plot(cmdscale(fraD(t(Dncall[bydg05$C,]))), xlim=c(-1,1), ylim=c(-1,1),pch=16, cex=.75, cex.axis=0.75, cex.lab=0.75, cex.main=0.75, col=dgi, xlab="D1", ylab="D2", main="MDS on significant spots for Control")
legend("bottomleft",levels(dgall), col = 1:length(levels(dgall)), pch = 16, bty="n", cex=0.75)
dev.off()

pdf(file="loopredcomx3C4.pdf", width=6, height=6)
LOOpred(S=comx3C[[4]], Fl=dgallC, P=pats,  g=0.0005, C=1000, dm=2, tit="C", figs = TRUE)
dev.off()

pdf(file="loopredcomx3ML4.pdf", width=6, height=6)
LOOpred(S=comx3ML[[4]], Fl=dgallML, P=pats,  g=0.001, C=1000, dm=2, tit="ML", figs = TRUE)
dev.off()

pdf(file="loopredcomx3GBM11.pdf", width=6, height=6)
LOOpred(S=comx3GBM[[11]], Fl=dgallGBM, P=pats,  g=0.5, C=1000, dm=2, tit="GBM", figs = TRUE)
dev.off()

pdf(file="loointersectionsGBM.pdf", width=4, height=4)
require(plot3D)
scatter3D(rsprimgGBM[,1],rsprimgGBM[,2],rsprimgGBM[,3], type="h", ticktype="detailed", xlab="Left limit", ylab="Right limit", zlab="% Correctly predicted", cex=0.3, cex.axis=0.5,
    xlim=c(1,15), ylim=c(1,15), zlim=c(0,1))
dev.off()

pdf(file="loointersectionsC.pdf", width=4, height=4)
require(plot3D)
scatter3D(rsprimgC[,1],rsprimgC[,2],rsprimgC[,3], type="h", ticktype="detailed", xlab="Left limit", ylab="Right limit", zlab="% Correctly predicted", cex=0.3, cex.axis=0.5,
          xlim=c(1,8), ylim=c(1,8), zlim=c(0,1))
dev.off()

pdf(file="loointersectionsML.pdf", width=4, height=4)
require(plot3D)
scatter3D(rsprimgML[,1],rsprimgML[,2],rsprimgML[,3], type="h", ticktype="detailed", xlab="Left limit", ylab="Right limit", zlab="% Correctly predicted", cex=0.3, cex.axis=0.5,
          xlim=c(1,8), ylim=c(1,8), zlim=c(0,1))
dev.off()

pdf(file="loocomx3GBMall.pdf", width=6, height=6)
#par(mai=c(0.75,0.75,0.5,0.5))
#plot(cmdscale(fraD(t(Dncall[rownames(loocomx3GBMall),]))),col=dgi,cex.axis=0.5, cex.lab=0.5,pch=16, xlab="D1", ylab="D2")
x=rownames(loocomx3GBMall)
LOOpred(S=x, Fl=dgallGBM, P=pats,  g=0.001, C=100, dm=4, tit="GBM", figs = TRUE)
dev.off()

pdf(file="loocomx3Call.pdf", width=6, height=6)
x=rownames(loocomx3Call)
LOOpred(S=x, Fl=dgallC, P=pats,  g=0.0008, C=100, dm=4, tit="Control", figs = TRUE)
dev.off()

pdf(file="loocomx3MLall.pdf", width=6, height=6)
x=rownames(loocomx3MLall)
LOOpred(S=x, Fl=dgallML, P=pats,  g=0.5, C=100, dm=4, tit="ML", figs = TRUE)
dev.off()

pdf(file="loocomx3ALL_GBM.pdf", width=6, height=6)
x=loocomx3ALL
LOOpred(S=x, Fl=dgallGBM, P=pats,  g=0.3, C=100, dm=3, tit="ALL - GBM", figs = TRUE)
dev.off()

pdf(file="loocomx3ALL_C.pdf", width=6, height=6)
x=loocomx3ALL
LOOpred(S=x, Fl=dgallC, P=pats,  g=0.3, C=100, dm=3, tit="ALL - GBM", figs = TRUE)
dev.off()

pdf(file="loocomx3ALL_ML.pdf", width=6, height=6)
x=loocomx3ALL
LOOpred(S=x, Fl=dgallML, P=pats,  g=.001, C=100, dm=3, tit="ALL - GBM", figs = TRUE)
dev.off()

###################
x=sapply(x,function(s){c(length(setdiff(s,dx3GBM)),length(intersect(s,dx3GBM)),length(setdiff(dx3GBM,s)))})
i=which.min(x[2,]/(rowMins(cbind(x[2,]+x[1,],x[2,]+x[3,]))))

pdf(file="loodx3comx3GBM.pdf", width=6, height=6)
x=intersect(comx3GBM[[5]],dx3GBM)
LOOpred(S=x, Fl=dgallGBM, P=pats,  g=0.006, C=100, dm=2, tit="GBM recursive intersection", figs = TRUE)
dev.off()

##################
cmds=t(Dn0[rownames(dhx3GBMisctbl[dhx3GBMisctbl>4]),])
d3heatmap(cmds,distfun = fraD, hclustfun=hclwrd, scale="none",colors=csch, show_grid = FALSE, cexRow = 0.75)
pdf(file="dhx3comisctbl_5up_MDS.pdf", width=6, height=6)
cmds=cmdscale(fraD(t(Dn0[rownames(dhx3GBMisctbl[dhx3GBMisctbl>4]),])),k=3)
par(mai=c(1,1,1,1))
plot(cmds, col=dgall, pch=16, cex=0.8, xlab="D1", ylab="D2", main="MDS on bootstrap derived features for GBM")
legend("bottomright",levels(dgall), col = 1:length(levels(dgall)), pch = 16, bty="n", cex=0.75)
dev.off()

##################

pdf(file="VennsGBM_dhx.pdf",onefile = FALSE,  width=4, height=12)
require(eulerr)

featbdg=unique(unlist(bydg05))
X1=list(comx3C[[1]],comx3GBM[[1]],comhx3GBM)
X2=list(comx3C[[4]],comx3GBM[[7]],comhx3GBM)
X3=list(comx3C[[7]],comx3GBM[[14]],comhx3GBM)
names(X1)=c("Control","GBM","GBM-NOB")
names(X2)=c("Control","GBM","GBM-NOB")
names(X3)=c("Control","GBM","GBM-NOB")
mx_dfxcomx1=sapply(X1,function(s){featbdg %in% s}, USE.NAMES = TRUE)
mx_dfxcomx2=sapply(X2,function(s){featbdg %in% s}, USE.NAMES = TRUE)
mx_dfxcomx3=sapply(X3,function(s){featbdg %in% s}, USE.NAMES = TRUE)
mx_dfxcomx=rbind(mx_dfxcomx2,mx_dfxcomx1,mx_dfxcomx3)
mx1=rep("Full overlap",nrow(mx_dfxcomx1))
mx2=rep("Overlap of half",nrow(mx_dfxcomx2))
mx3=rep("Overlap in 2",nrow(mx_dfxcomx3))
mxx=c(mx2,mx1,mx3)
venn_comdf=euler(mx_dfxcomx, by=mxx)
plot(venn_comdf, layout=c(1,3), counts=TRUE, auto.key=FALSE, cex=0.6)

#venn_comdf2=euler(mx_dfxcomx2[,1:4], by=mx_dfxcomx2[,5])
#plot(venn_comdf2)
#venn_comdf3=euler(mx_dfxcomx3[,1:4], by=mx_dfxcomx3[,5])
#plot(venn_comdf3)
dev.off()

pdf(file="MinkoD_REA.pdf", width=4, height=4)
x=cbind(DLOOX3GBM[[1]][[3]][,1],DLOOX3GBM[[2]][[3]][,1],DLOOX3GBM[[3]][[3]][,1],DLOOX3GBM[[4]][[3]][,1],DLOOX3GBM[[5]][[3]][,1],DLOOX3GBM[[6]][[3]][,1])
y1=cbind(DLOOX3GBM[[1]][[3]][,2],DLOOX3GBM[[2]][[3]][,2],DLOOX3GBM[[3]][[3]][,2],DLOOX3GBM[[4]][[3]][,2],DLOOX3GBM[[5]][[3]][,2],DLOOX3GBM[[6]][[3]][,2])
y2=cbind(DLOOE3GBM[[1]][[3]][,2],DLOOE3GBM[[2]][[3]][,2],DLOOE3GBM[[3]][[3]][,2],DLOOE3GBM[[4]][[3]][,2],DLOOE3GBM[[5]][[3]][,2],DLOOE3GBM[[6]][[3]][,2])
matplot(x,y1,type="l", lty=1, lwd=1.5, ylim=c(0,max(y1)),xlab="Number of features", ylab="Clustering criterion", cex.axis=0.5, cex.lab=0.75)
par(new=TRUE)
matplot(x,y2, type="l", lty=1, lwd=0.75,  ylim=c(0, max(y1)),cex=0.1, xlab="", ylab="", axes=FALSE)
dev.off()


#################

pdf(file="VennsGBM_15bootstr.pdf",onefile = FALSE,  width=4, height=12)
require(eulerr)

featbdg=unique(unlist(bydg05))
X=dloox3GBM
mx_dloox3GBM=sapply(X,function(s){featbdg %in% s})
venn_dloox3GBM=euler(mx_dloox3GBM[,1:10])
plot(venn_dloox3GBM, counts=TRUE, auto.key=FALSE, cex=0.6)
dev.off()

###########
pdf(file="randomdraw.pdf",onefile = FALSE,  width=12, height=12)
rndraw=matrix(rep(0,10000),1000)
rndraw0=sapply(1:1000,function(y){unlist(table(table(unlist(sapply(dloox3GBM,function(x){sample(1:468,length(x))})))))})
for ( i in 1:1000) rndraw[i,1:length(rndraw0[[i]])]=c(rndraw0[[i]])
boxplot(rndraw, ylim=c(0,160), xlim=c(0,16), ylab="Number of features in the overlap", xlab="Number of overlapping sets")
par(new=TRUE)
plot(1:15,table(table(unlist(dloox3GBM))), ylim=c(0,160), xlim=c(0,16), xlab="", ylab="")
lines(1:15,table(table(unlist(dloox3GBM))), ylim=c(0,160), xlim=c(0,16))


dev.off()
#################

pdf(file="hmGBM11.pdf",width=6, height=12)

cpl1=colorRampPalette(c("black","black","#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))(n=128)
heatmap.2(cmds,distfun = fraD, hclustfun=hclwrd, scale="none",col=cpl1,  trace="none", cexRow = 0.5, cexCol = 0.75)


dev.off()
#################
