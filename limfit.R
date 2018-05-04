# A wrapper for the limma linear model based extraction of calls

limfit<-function(fdsn,vals, n=1000, p=0.05, ip=FALSE, cmf=1){
  require(limma)
  require(sva)
  require(parallel)
  m=length(unique(fdsn))
  ncl=ncol(vals)
  if (ncl!=length(fdsn)) {return("Error: ncol of design and matrix differ!")}
  clin=character(0)
  #fdsn=factor(dsn)
  #fsex=factor(sex)
  
  if (cmf==1) {
    design <- model.matrix(~0+fdsn)
    cmx=rep(0,m)
    ix=combn(1:m,2)
    cmx=apply(ix,2,function(i){a=cmx; a[i[1]]=1; a[i[2]]=-1; a})
    fit=lmFit(vals,design)
    fit2=contrasts.fit(fit,cmx)
    fit2=eBayes(fit2, trend=FALSE, proportion = 0.5)
    tptb=topTable(fit2,number=n) #,adjust.method="holm"
    clin=rownames(tptb[tptb$adj.P.Val<p,])
    res=unique(unlist(clin))
    #print(length(clin))
    #print(head(tptb))
    }
  if (cmf==2) {
    tptb=list()
    clin=list()
   
    for (i in 1:m){
      fdsnr=relevel(fdsn,ref=i)
      for (j in 1:(m-1)){
      design = model.matrix(~fdsnr)
      cn=paste("f",0:(m-1),sep ="")
      colnames(design)=cn
      cntr=paste("f",j,"-f0", sep="")
      cmx=makeContrasts(contrasts=cntr, levels=cn)
      fit=lmFit(vals,design)
      fit2=contrasts.fit(fit,cmx)
      fit2=eBayes(fit2, trend=FALSE, proportion = 0.5)
      tptb=topTable(fit2,number=n) #,adjust.method="holm"
      #print(tptb[1:20,])
      #if (nrow(tptb[tptb$adj.P.Val<p,])>=nr) clin[[i]]=rownames(tptb[1:nr,]) else clin[[i]]=NULL
      nr=nrow(tptb[tptb$adj.P.Val<p,])
      print(c(i,j))
      nm=paste(fdsnr[1],fdsnr[j+1])
      clin[[(i-1)*3+j]]=list(rownames(tptb[1:nr,]))
      names(clin)[[(i-1)*3+j]]=nm
    }
    }
 
    l=lengths(clin)
    print(l)
    #print(head(tptb))  
    res=clin
  }
  if (cmf==4) {
    tptb=list()
    clin=list()
  
    for (i in 1:m){
      fdsnr=relevel(fdsn,ref=i)
      #fdsnr=factor(fdsn, levels = levels(fdsn)[c(i:m,(i>1):(i-1))])
      design = model.matrix(~fdsnr)
      cn=paste("f",0:(m-1),sep ="")

      colnames(design)=cn
      #cntr="f1-f0"
      cntr=sapply(1:(m-1),function(i){paste("f",i,"-f0", sep="")})
      #cntr="f0"
      cmx=makeContrasts(contrasts=cntr, levels=cn)
      fit=lmFit(vals,design)
      fit2=contrasts.fit(fit,cmx)
      fit2=eBayes(fit2, trend=FALSE, proportion = 0.5)
      tptb=topTable(fit2,number=n) #,adjust.method="holm"
      #print(tptb[1:20,])
      #if (nrow(tptb[tptb$adj.P.Val<p,])>=nr) clin[[i]]=rownames(tptb[1:nr,]) else clin[[i]]=NULL
      nr=nrow(tptb[tptb$adj.P.Val<p,])
      clin[[i]]=rownames(tptb[1:nr,])

    }
    names(clin)=levels(fdsn)
    l=lengths(clin)
    print(l)
    #print(head(tptb))  
    res=clin
  }
  if (cmf==5) {
    tptb=list()
    clin=list()
    i=1
    for (i in 1:m){
      fdsnr=relevel(fdsn,ref=i)
      design = model.matrix(~fdsnr+batches)
      cn=c("f0", "f1", "f2","f3","f4","b1","b2")
      colnames(design)=cn
      cmx=makeContrasts(f1-f0, f2-f0, f3-f0, f4-f0, levels=cn)
      fit=lmFit(vals,design)
      fit2=contrasts.fit(fit,cmx)
      fit2=eBayes(fit2, trend=TRUE)
      tptb=topTable(fit2,number=n) #,adjust.method="holm"
      #print(tptb[1:20,])
      #if (nrow(tptb[tptb$adj.P.Val<p,])>=nr) clin[[i]]=rownames(tptb[1:nr,]) else clin[[i]]=NULL
      nr=nrow(tptb[tptb$adj.P.Val<p,])
      clin[[i]]=rownames(tptb[1:nr,])
      i=i+1
    }
    names(clin)=levels(fdsn)
    l=lengths(clin)
    print(l)
    #print(head(tptb))  
    res=clin
  }   
  return(res)
  
}  

