# A script converting the libraries peppos (SYM), peprnd (RND) and
# pepnegrnd (NGR) to fasta format, saving themas txt files and
# processing the XML files from the results of the on-line blastp
# search. XML is converted to list followed by blastQstrip() to parse 
# the XML derived lists. 


require(XML)
require(multcomp)
require(reshape2)
require(RVAideMemoire)

load(file="C10MnSD")
pepcl=rownames(c10MnSD)
pepcl=substr(pepcl,1,7)
pepposnoneg=read.csv("pepposnoneg.csv",header = F)
pepnrnd=pepcl[c10MnSD[,3]=="pepnegrnd"]
peprnd=pepcl[c10MnSD[,3]=="peprnd"]
peppos=substr(pepposnoneg[,1],1,7)
pps=list(peppos,pepnrnd,peprnd)
fnpos="peppos.txt"
fnnrnd="pepnrnd.txt"
fnrnd="peprnd.txt"
fns=list(fnpos,fnnrnd,fnrnd)
for (n in 1:3){
sink(fns[[n]], append = T)
for(w in seq_along(pps[[n]])){
  p=pps[[n]][w]
  cat(paste(">",w,"\n",sep=""))
  cat(paste(p,"\n",sep=""))
}
sink()
}

posalnm=xmlToList(xmlParse("peppos.xml"))
rndalnm=xmlToList(xmlParse("peprnd.xml"))
nrndalnm=xmlToList(xmlParse("pepnrnd.xml"))

tbpos=blastQStrip(posalnm)
names(tbpos)=peppos
tbrnd=blastQStrip(rndalnm)
names(tbrnd)=peprnd
tbnrnd=blastQStrip(nrndalnm)
names(tbnrnd)=pepnrnd


# Extracting data about the hits that fall within the
# preferred limits of alignment parameters

tbrnd1=lapply(tbrnd,function(t){return(t[t$Positive>5 & t$Identity>5 & t$AlLength==t$Positive,])})
tbpos1=lapply(tbpos,function(t){return(t[t$Positive>5 & t$Identity>5 & t$AlLength==t$Positive,])})
tbnrnd1=lapply(tbnrnd,function(t){return(t[t$Positive>5 & t$Identity>5 & t$AlLength==t$Positive,])})

imjp=sum(unlist(lapply(tbpos1,function(j){
x=any(grepl("immunoglobulin",rownames(j)))
})))

imjr=sum(unlist(lapply(tbrnd1,function(j){
  x=any(grepl("immunoglobulin",rownames(j)))
})))

imjn=sum(unlist(lapply(tbnrnd1,function(j){
  x=any(grepl("immunoglobulin",rownames(j)))
})))

# Summarize data and plot relative and absolute frequences of 
# immunoglobulin and non-immunoglobulin hits for each library

z1=unlist(lapply(tbpos1,nrow))
z2=unlist(lapply(tbrnd1,nrow))
z3=unlist(lapply(tbnrnd1,nrow))

print(imjp/length(z1[z1>0]))
print(imjr/length(z2[z2>0]))
print(imjn/length(z3[z3>0]))

print(imjp/length(z1))
print(imjr/length(z2))
print(imjn/length(z3))

imjpl=unlist(lapply(tbpos1,function(j){
  x=any(grepl("immunoglobulin",rownames(j)))*2
  y=(nrow(j)>0)*1
  return(x+y)
}))

imjrl=unlist(lapply(tbrnd1,function(j){
  x=any(grepl("immunoglobulin",rownames(j)))*2
  y=(nrow(j)>0)*1
  return(x+y)
}))

imjnl=1*(unlist(lapply(tbnrnd1,function(j){
  x=any(grepl("immunoglobulin",rownames(j)))*2
  y=(nrow(j)>0)*1
  return(x+y)
})))

imjl=rbind(cbind(imjpl,rep("p",length(imjpl))),cbind(imjrl,rep("r",length(imjrl))),cbind(imjnl,rep("n",length(imjnl))))
timjl=table(imjl[,1],imjl[,2])
stimjl=sweep(timjl,2,colSums(timjl),FUN="/")
colnames(stimjl)=c("NGR","SYM","RND")
rownames(stimjl)=c("No Hit","Non-Ig Hit","Hit in Ig J Region")
par(mar=c(5,4,5,9),xpd=T)
barplot(stimjl[3:1,], col=c(1,rgb(0.5,0.5,0.5,0.5),0), main="Proportion of Number of Hits",xlab="Library", bty="L")
legend("topright",inset=c(-0.5,0),legend=rownames(stimjl), fill=c(0,rgb(0.5,0.5,0.5,0.5),1))
text(c(0.75,0.75,0.75,1.9,1.9,1.9,3.15,3.15,3.15),c(0.023,0.18,0.67,0.085,0.3,0.72,0.058,0.3,0.75), c(timjl[3:1,]), col=c(0,1,1,0,1,1,0,1,1))
par(mar=c(5.1,4.1,4.1,2.1),xpd=F)
                                    
GBM14_55pp=colnames(GBM14df)
GBM14_55pp=GBM14_55pp[1:55]
GBMimjl=imjl[rownames(imjl) %in% GBM14_55pp,]
tGBM=table(GBMimjl[,1],GBMimjl[,2])

x=pcaGBM14$var$contrib[,1:5]
GBMbest=apply(x,2,function(j){rownames(x)[j>6]})
GBMbest=unlist(GBMbest)

x=imjl[rownames(imjl) %in% GBMbest,]
table(x[,1],x[,2])

timjl=cbind(timjl, tGBM[,2])
stimjl=cbind(stimjl,tGBM[,2]/sum(tGBM[,2]))
colnames(stimjl)=c("NGR","SYM","RND", "GBM")
rownames(stimjl)=c("No Hit","Non-Ig Hit","Hit in Ig J Region")
par(mar=c(5,4,5,9),xpd=T)
barplot(stimjl[3:1,], col=c(1,rgb(0.5,0.5,0.5,0.5),0), main="Proportion of Number of Hits",xlab="Library", bty="L")
legend("topright",inset=c(-0.5,0),legend=rownames(stimjl), fill=c(0,rgb(0.5,0.5,0.5,0.5),1))
text(c(0.75,0.75,0.75,1.9,1.9,1.9,3.15,3.15,3.15,4.3,4.3,4.3),c(0.023,0.18,0.67,0.085,0.3,0.72,0.058,0.3,0.75,0.12,0.36,0.75), c(timjl[3:1,]), col=c(0,1,1,0,1,1,0,1,1))
par(mar=c(5.1,4.1,4.1,2.1),xpd=F)

timjl=apply(timjl,2,as.double)
chisq.test(timjl)
multinomial.theo.multcomp(timjl[,4],timjl[,2]/sum(timjl[,2]))
