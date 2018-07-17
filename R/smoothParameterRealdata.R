smoothParameterRealdata<-function(mxBeta,covEvector,vxlNberSeq,vxlNbySeqid,vxlDistSeq,xMatrix,covEvalue,sigmaError)
{
  result=list(fnlBeta=0,test=0);
  
  inlBeta = t(mxBeta);
  xDesign = xMatrix;

  ch = 1.1;
  ss = 1:10;
  hSeq = ch^ss;
  
  chiSeq =qchisq(0.8/(ss[4:10]-2),1);
  
  n = dim(xDesign)[1];
  Cn =log(n)*qchisq(0.95,1);
  
  IXX=solve(t(xDesign)%*%xDesign);
  result1=aspV2NYU(inlBeta,covEvector,hSeq,vxlNberSeq,vxlNbySeqid,vxlDistSeq,xDesign,chiSeq,sigmaError,Cn,IXX,covEvalue);
  result[[1]]=result1[[1]];
  result[[2]]=result1[[3]];
  result
}