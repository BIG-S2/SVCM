smoothResidualRealdata_para<-function(newXdesign,newImgData,wmNonzeroId,dim,part,cores=2)
{
  result1=mlrWresd(newXdesign,newImgData);
  mxBeta=result1[[1]];
  mxR=result1[[3]];
  SRImgdata=RImgdata=mxR;
  id=wmNonzeroId;
  Dimn=dim;
  wmImgsMask=array(rep(0,prod(Dimn)),Dimn);
  wmImgsMask[id]=1;
  Mask=wmImgsMask;
  cord=which(wmImgsMask !=0, arr.ind = T)
  Lx = seq(from=1,to=Dimn[1],by=part[1]);
  Ux = seq(from=part[1],to=Dimn[1],by=part[1]);
  Ly = seq(from=1,to=Dimn[2],by=part[2]);
  Uy = seq(from=part[2],to=Dimn[2],by=part[2]);
  Lz = seq(from=1,to=Dimn[3],by=part[3]);
  Uz = seq(from=part[3],to=Dimn[3],by=part[3]);
  Nx = Ux-Lx+1;
  Ny = Uy-Ly+1;
  Nz = Uz-Lz+1;
  NH = floor(max(15,max(c(Nx,Ny,Nz)/2)));
  HSeq = logspace(log10(0.99),log10(max(c(Nx,Ny,Nz))),NH)/4; 
  sequ=0
  bb=wmNonzeroId
  
  for(ii in 1:length(Lx))
  {
    XStartInd = Lx[ii];
    XEndInd = Ux[ii];
    for(jj in 1:length(Ly))
    {
      YStartInd = Ly[jj];
      YEndInd = Uy[jj];
      for(kk in 1:length(Lz))
      {
        ZStartInd = Lz[kk];
        ZEndInd = Uz[kk];
        XYZCoord = matrix(rep(0,Nx[ii]*Ny[jj]*Nz[kk]*3),ncol=3);
        NonZeroCount = 0;
        for (ZIndk in ZStartInd:ZEndInd)
          for (YIndj in YStartInd:YEndInd)
            for (XIndi in XStartInd:XEndInd)
              if(Mask[XIndi,YIndj,ZIndk]!=0)
              {
                NonZeroCount = NonZeroCount + 1; 
                XYZCoord[NonZeroCount,] = c(XIndi,YIndj,ZIndk);
              }
        if(NonZeroCount>0)
        {
          sequ=sequ+1;
          XYZCoord = matrix(XYZCoord[1:NonZeroCount,],ncol=3);
          XYZInd = XYZCoord[,1] + (XYZCoord[,2]-1)*Dimn[1] + (XYZCoord[,3]-1)*Dimn[1]*Dimn[2];
          Tid =rep(0,length(XYZInd));
          mm=1;
          for(ll in 1:length(XYZInd))
          {
            while(id[mm]!= XYZInd[ll])
              mm = mm+1;
            bb[mm]=sequ 
          }
        }
      }
    }
  }
  
  
  
  len=max(bb)
  library(parallel)
  res1.p <- mclapply(1:len, FUN = function(x) { smooth(x,bb,RImgdata,cord)}, mc.cores = cores)
  tt=matrix(unlist(res1.p),nrow = dim(RImgdata)[1])
  order=c()
  
  for (i in 1:length(res1.p))
  {
    order=c(order,which(bb==i))
    
  }
  SRImgdata[,order]=tt
  mError = mxR -  SRImgdata;
  sigmaError = apply(mError,2,var);
  Cov=SRImgdata%*%t(SRImgdata);
  covEvector= (eigen(Cov)$vectors);
  s1=eigen(Cov)$values;
  aa=1
  for(i in 1:length(s1))
    if(sum(s1[1:aa])<0.7*sum(s1))
      aa=aa+1
  covEvalue=s1[1:aa];
  covEvector=(t(SRImgdata)%*%covEvector)[,1:aa];
  for(i in 1:dim(covEvector)[2])
    covEvector[,i]=covEvector[,i]/sqrt(sum(covEvector[,i]^2))
  result2=list(mxBeta=0,sigmaError=0,covEvalue=0,covEvector=0);
  result2[[1]]=mxBeta;
  result2[[2]]=sigmaError;
  result2[[3]]=covEvalue;
  result2[[4]]=covEvector;
  result2;
}

