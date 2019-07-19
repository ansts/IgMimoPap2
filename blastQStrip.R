# A function extracting information from XML derived list Q
# representing results from batch NCBI/BLASTp as a list of tables 
# with columns: 
# bit score/E value/Length Positive/Length Identity/Alignment Length/Query Sequence/Hit Sequence
# - one table per query sequence. 

blastQStrip=function(Q){
  z=Q$BlastOutput_iterations
  z=lapply(z, function(u){
    u=u$Iteration_hits
    ti=c()
    for (uj in u){
      un=strsplit(uj$Hit_def,split = "\\[")[[1]][1]
      ul=c(uj$Hit_hsps$Hsp$`Hsp_bit-score`,uj$Hit_hsps$Hsp$`Hsp_evalue`,uj$Hit_hsps$Hsp$`Hsp_positive`,uj$Hit_hsps$Hsp$`Hsp_identity`,uj$Hit_hsps$Hsp$`Hsp_align-len`,uj$Hit_hsps$Hsp$`Hsp_qseq`,uj$Hit_hsps$Hsp$`Hsp_hseq`)
      ul=data.frame(Bits=as.double(ul[1]),E=as.double(ul[2]),Positive=as.double(ul[3]),Identity=as.double(ul[4]),AlLength=as.double(ul[5]),Query=ul[6],Target=ul[7], row.names=un, stringsAsFactors = F)
      ti=rbind(ti,ul)
    }
    return(ti)
  })
  names(z)=NULL
  return(z)
}