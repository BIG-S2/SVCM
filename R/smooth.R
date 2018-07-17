smooth<-function(l,bb,RImgdata,cord)
{
  Tid=which(bb==l)
  TRImgdata = matrix(RImgdata[,Tid],ncol=length(Tid));
  XYZCoord= matrix(cord[Tid,],nrow=length(Tid));
  m=siv(TRImgdata,XYZCoord,HSeq);
  if(m[[2]]==0)
  {
    a = TRImgdata;
  }
  else a = m[[1]];
  return(a)
  
}