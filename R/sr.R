sr<-function(RImgdata,Mask)
{
  SRImgdata=matrix(rep(0,dim(RImgdata)[1]*dim(RImgdata)[2]),dim(RImgdata)[1]);
  Dimn=dim(Mask);
  id=which(c(Mask)!=0);
  Lx = seq(from=1,to=64*2,by=8);
  Ux = seq(from=8,to=64*2,by=8);
  Ly = seq(from=1,to=64*2,by=8);
  Uy = seq(from=8,to=64*2,by=8);
  Lz = seq(from=1,to=48*2,by=6);
  Uz = seq(from=6,to=48*2,by=6);
  Nx = Ux-Lx+1;
  Ny = Uy-Ly+1;
  Nz = Uz-Lz+1;
  NH = floor(max(15,max(c(Nx,Ny,Nz)/2)));
  HSeq = logspace(log10(0.999999),log10(max(c(Nx,Ny,Nz))),NH)/4; 
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
          XYZCoord = XYZCoord[1:NonZeroCount,];
          if(NonZeroCount==1)
          {
            XYZInd = XYZCoord[1] + (XYZCoord[2]-1)*Dimn[1] + (XYZCoord[3]-1)*Dimn[1]*Dimn[2];
            XYZCoord=matrix(c(XYZCoord),1,length(XYZCoord))
          }
          if(NonZeroCount>1)
            XYZInd = XYZCoord[,1] + (XYZCoord[,2]-1)*Dimn[1] + (XYZCoord[,3]-1)*Dimn[1]*Dimn[2];
          Tid =rep(0,length(XYZInd));
          mm=1;
          for(ll in 1:length(XYZInd))
          {
            while(id[mm]!= XYZInd[ll])
              mm = mm+1;
            Tid[ll] = mm; 
          }
          TRImgdata = as.matrix(RImgdata[,Tid]);
          if(length(XYZInd)>0)
          {
            m=siv(TRImgdata,XYZCoord,HSeq);
            if(m[[2]]==0)
            {
              SRImgdata[,Tid] = TRImgdata;
            }
            else SRImgdata[,Tid] = m[[1]];
          }
        }
      }
    }
  }
  SRImgdata;
}