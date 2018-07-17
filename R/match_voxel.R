match_voxel<-function(l,VSeq,InVSeq,wmNonzeroId)
{
  subVSeq = VSeq[l,1:InVSeq[l]];
  loc = match(subVSeq,wmNonzeroId);
  loc[is.na(loc)]=0
  loc=c(loc,rep(0,dim(VSeq)[2]-InVSeq[l]))
  loc
  
}
