# A funciton which iteratively scans the set of features and at each cycle finds the feature  
# the removal of which imrpoves most the separation between the classes of interest c, then removes it 
# leaving one less features for the next cycle. The separation is measured using a composite clustering 
# criterion multiplying Dunn's and BHgamma criteria and the inverse of Connectivity.

scan3crit=function(m0,c,d=2,random=FALSE){
  require(clValid)
  require(parallel)
  require(reshape2)
  require(matrixStats)
  if (is.logical(c)) cc=as.integer(c)+1 else cc=as.integer(c)
  pt=proc.time()
  D0=c()
  bad=list()
  m=m0
  nc=unique(c)
  repeat {
      D=c()
      N=nrow(m)
      print(N)
      if (N<3) break
     # d=fraD(t(m),f=.1)

      cl <- makeCluster(6)
      ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
      clusterEvalQ(cl, {library(clValid)})
      clusterExport(cl, ex)
      clusterExport(cl, list("c","cc","m"), envir=environment())
      if (random==TRUE) d=sample(c(0.25,0.5,1,2),1)
      D=parSapply(cl,1:N,  function(i) {
          d1=fraD(t(m[-i,]),f=d)
          cri3fix(N)*(dunn(d1, cc)*BHgamma(d1,c))/(connectivity(d1,cc)+1)
        })
      
      stopCluster(cl)
      md=max(D)
      D0=rbind(D0,c(N,mean(D), range(D)))
      bad=c(bad,rownames(m)[D==md])
      ix=sapply(D,function(i){i!=md})
      m=m[ix,]
      #print(md)
  }
  bad=unlist(bad)
  plot(D0[,c(1,2)], ylim=range(D0[,2:4]), ylab="D0")
  par(new=TRUE)
  plot(D0[,c(1,3)], ylim=range(D0[,2:4]), cex=0.3, ylab=NULL)
  lines(D0[,c(1,3)], ylim=range(D0[,2:4]), ylab=NULL)
  par(new=TRUE)
  plot(D0[,c(1,4)], ylim=range(D0[,2:4]), cex=0.3, ylab=NULL)
  lines(D0[,c(1,4)], ylim=range(D0[,2:4]), ylab=NULL)

  crit=which(D0[,2]==max(D0[,2]))
  m=m0[!(rownames(m0) %in% bad[1:(crit-1)]),]
  pt=proc.time()-pt
  print(pt)
  return(list(m,bad,D0))
}