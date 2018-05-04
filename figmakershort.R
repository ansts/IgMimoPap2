
require(gplots)
pdf(file="hmGBM15.pdf",width=6, height=10)
cpl1=colorRampPalette(c("black","black","#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))(n=128)
heatmap.2(Dncc[daloox3GBMe_tco[[15]],],distfun = dist,  scale="none",col=cpl1,  trace="none", cexRow = 0.5, cexCol = 0.75)
dev.off()
#################
#################
require(gplots)
pdf(file="hmC15.pdf",width=6, height=10)
cpl1=colorRampPalette(c("black","black","#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))(n=128)
heatmap.2(Dncc[daloox3Ce_tco[[15]],],distfun = dist,  scale="none",col=cpl1,  trace="none", cexRow = 0.5, cexCol = 0.75)
dev.off()
#################
require(gplots)
pdf(file="selpats.pdf",width=4, height=4)
cpl1=colorRampPalette(c("black","black","#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))(n=128)
heatmap.2(Dncall[peps3,selpats], key=FALSE,dendrogram="none",  Rowv=FALSE, Colv=FALSE,
scale="none",col=cpl1,  trace="none", cexRow = 0.5, cexCol = 0.75)
dev.off()
#################
require(gplots)
pdf(file="hmML15.pdf",width=6, height=10)
cpl1=colorRampPalette(c("black","black","#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))(n=128)
heatmap.2(Dncc[daloox3MLe_tco[[15]],],distfun = dist,  scale="none",col=cpl1,  trace="none", cexRow = 0.5, cexCol = 0.75)
dev.off()
#################
#################
pdf(file="MCCtr_val.pdf",width=6, height=9)
op=par(mfrow=c(2,1), omi=c(0.85,0.85,0,0), mai=c(.5,0,0.75,0))
barplot(MCCesvmt, ylim = c(-0.5,1), main="Training")
barplot(MCCesvmv, ylim = c(-0.5,1), main="Validation")
title(xlab = "Commonality", ylab="MCC", outer=TRUE)
par(op)
dev.off()
######################
pdf(file="cmdsvalGBM.pdf",  width=7, height=6)
require(plot3D)
y=cmdscale(dist(t(Dn3[daloox3GBMe_tco[[15]],])))
ym=svm(y[!d3New,], as.factor(dg3GBM[!d3New]), gamma = 0.1, cost=1000, probability = TRUE)
pry=sapply(seq(range(y[,1])[1], range(y[,1])[2], length.out = 10), function(x1){
      sapply(seq(range(y[,2])[1], range(y[,2])[2], length.out = 10), function(x2){
        attr(predict(ym, newdata=t(c(x1,x2)), probability = TRUE), "probabilities")[1]
      })
})
op=par(mai=c(1,1,1,1))
image2D(t(pry),resfac = 10, col=cm.colors(1000), xaxt="n", yaxt="n", xlab="", ylab="", colkey = NULL)
op=par(new=TRUE,mai=c(1,1,1,1.56))
plot(cmdscale(dist(t(Dn3[daloox3GBMe_tco[[15]],]))), col=as.integer(dg3),cex=1.5, pch=16, xlab="D1", ylab="D2")
par(new=TRUE)
plot(cmdscale(dist(t(Dn3[daloox3GBMe_tco[[15]],]))), col=as.integer(!(colnames(Dn3) %in% Dncllnms)),cex=2.5, xlab="", ylab="")
legend("bottomright",levels(dg3), col = 1:length(d2), pch = 16, bty="n", cex=0.75)
par(op)
dev.off()
#######################

pdf(file="Venn14.pdf",onefile = FALSE,  width=5, height=5)
require(eulerr)
daloox14e_=unique(c(daloox3Ce_tco[[14]], daloox3GBMe_tco[[14]],daloox3MLe_tco[[14]]))
venndaloox14=sapply(list(C=daloox3Ce_tco[[14]], GBM=daloox3GBMe_tco[[14]], ML=daloox3MLe_tco[[14]]),function(l){daloox14e_ %in% l})
venndl14=euler(venndaloox14)

plot(venndl14, counts=TRUE, auto.key=FALSE)
dev.off()
######################
pdf(file="MinComm.pdf",onefile = FALSE,  width=6, height=4)
barplot(lengths(daloox3GBMe_tco[2:28]), names.arg=2:28, cex.axis = 0.5, xlab="Minimal Commonality", ylab="N of Features")
dev.off()
######################
pdf(file="MDSDncll.pdf",onefile = FALSE,  width=3, height=3)
op=par(mai=c(0.8,0.8,0.2,0.2))
plotMDS(Dncll[vbC05_$C,], col=as.integer(dgll), labels = NULL, pch=16, cex=1, cex.axis=0.3, cex.lab=0.5, xlab="D1", ylab = "D2")
legend("bottomright",levels(dgll), col = 1:length(d2), pch = 16, bty="n", cex=0.5)
par(op)
dev.off()

######################
pdf(file="Traces.pdf",onefile = FALSE,  width=4, height=4)
par(mai=c(1,1,0.5,0.5))
for (i in sample(1:28,3)) {
  plot(DLOOXe_traces[[1]][[i]], ylim=c(0,8), cex=0, xlab="", ylab="", xaxt="n", yaxt="n")
  lines(DLOOXe_traces[[1]][[i]], ylim=c(0,8), col=1, xlab="", ylab="", xaxt="n", yaxt="n")
  par(new=TRUE)
  plot(DLOOXe_traces[[2]][[i]], ylim=c(0,8), cex=0, xlab="", ylab="", xaxt="n", yaxt="n")
  lines(DLOOXe_traces[[2]][[i]], ylim=c(0,8), col=2, xlab="", ylab="", xaxt="n", yaxt="n")
  par(new=TRUE)
  plot(DLOOXe_traces[[3]][[i]], ylim=c(0,8), cex=0, cex.lab=0.5, cex.axis=0.5, xlab="Number of Features", ylab="Clustering Criterion")
  lines(DLOOXe_traces[[3]][[i]], ylim=c(0,8), col=3, xlab="",ylab="", xaxt="n", yaxt="n")
  legend("topright",levels(dgll), col = 1:3, pch = 16, bty="n", cex=0.75)
  par(new=TRUE)
}
par(new=FALSE)
dev.off()





