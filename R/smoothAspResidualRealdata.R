smoothAspResidualRealdata<-function(Xdesign,Imgdata,aspBeta,wmNonzeroId,dimn,part)
{
  results=list(aspEta=0,sAspEta=0,covEvalue=0,covEvector=0);
  wmImgsMask=array(rep(0,prod(dimn)),dimn);
  wmImgsMask[wmNonzeroId]=1;
  Mask=wmImgsMask;
  
  aspBeta=t(aspBeta);
  aspEta = Imgdata-Xdesign%*%aspBeta;
  
  
  sAspEta= sr(aspEta,Mask,part);
  
  n=dim(Xdesign)[1];
  p=dim(Xdesign)[2];
  meanSAspEta=colMeans(sAspEta);
  covSRID = (sAspEta-kronecker(rep(1,n),t(meanSAspEta)))%*%t(sAspEta-kronecker(rep(1,n),t(meanSAspEta)))/(n-p);
  covEvector= eigen(covSRID)$vectors;
  covEvalue=eigen(covSRID)$values;
  covEvector = t(sAspEta-kronecker(rep(1,n),t(meanSAspEta)))%*%covEvector;
  
  results[[1]]=aspEta;
  results[[2]]=sAspEta;
  results[[3]]=covEvalue;
  results[[4]]=covEvector;
}
