#
# Function for recursive elimination feature selection based on a composite clustering criterion
# using Dunn's criterion, connectivity and Bake-Hubert gamma index from the clValid package
#

scan3crit=function(m0,c,d=0.25,random=FALSE, cntr=FALSE){
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
      if (cntr) print(N)
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
          100*(dunn(d1, cc)*BHgamma(d1,c))/connectivity(d1,cc)
        })
      
      stopCluster(cl)
      md=max(D)
      D0=rbind(D0,c(N,mean(D), range(D)))
      bad=c(bad,rownames(m)[D==md])
      ix=sapply(D,function(i){i!=md})
      m=m[ix,]
      if (cntr) print(md)
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