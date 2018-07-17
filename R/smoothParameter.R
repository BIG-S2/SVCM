smoothParameter<-function(vxlNbySeqid,vxlDistSeq, vxlNberSeq, mxR,smxR,mxBeta,xMatrix,hh,f1,f2)
{
  result=list(fnlBeta=0,test=0);
  
  sResidual = t(smxR);
  inlBeta = t(mxBeta);
  xDesign = xMatrix;
  
  mError = mxR - smxR;
  sigmaError = apply(mError,2,var);
  

  ch = 1.1;
  ss = 1:hh;
  hSeq = ch^ss;
  
  chiSeq = f1*qchisq(0.8/(ss[4:hh]-2),1);
  
  n = dim(xDesign)[1];
  Cn = f2*log(n)*qchisq(0.95,1);
  
  result1=aspV2(inlBeta,sResidual,hSeq,vxlNberSeq,vxlNbySeqid,vxlDistSeq,xDesign,chiSeq,sigmaError,Cn);
  result[[1]]=result1[[1]];
  result[[2]]=result1[[3]];
  result;
}