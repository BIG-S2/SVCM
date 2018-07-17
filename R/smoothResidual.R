smoothResidual<-function(RImgdata)
{
  Mask=array(rep(1,64*64*8),c(64,64,8));
  SRImgdata=RImgdata
  Dimn=dim(Mask);
  id=which(c(Mask)!=0);
  Lx = seq(from=1,to=64,by=8);
  Ux = seq(from=8,to=64,by=8);
  Ly = seq(from=1,to=64,by=8);
  Uy = seq(from=8,to=64,by=8);
  Lz = 1;
  Uz = 8;
  Nx = Ux-Lx+1;
  Ny = Uy-Ly+1;
  Nz = Uz-Lz+1;
  NH = floor(max(10,max(c(Nx,Ny,Nz)/2)));
  HSeq = logspace(log10(1),log10(max(c(Nx,Ny,Nz))),NH); 
  NonZeroCount=c();
  XYZCoord=matrix(rep(1,length(id)*3),length(id))
  t=0;
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
        pp=Mask[XStartInd:XEndInd,YStartInd:YEndInd,ZStartInd:ZEndInd]
        nn=sum(pp)
        NonZeroCount=c(NonZeroCount,nn)
        for (ZIndk in ZStartInd:ZEndInd)
          for (YIndj in YStartInd:YEndInd)
            for (XIndi in XStartInd:XEndInd)
              if(Mask[XIndi,YIndj,ZIndk]!=0)
              {
                t=t+1; 
                XYZCoord[t,] = c(XIndi,YIndj,ZIndk);
              }
      }
    }
  }
  NonZeroCount=NonZeroCount[NonZeroCount!=0]
  
  a=length(NonZeroCount)
  t1=1:a
  t2=cumsum(NonZeroCount)
  t1[2:a]=(t2+1)[1:(a-1)]
  
  for(i in 1:length(NonZeroCount))
  {
    Coord = XYZCoord[t1[i]:t2[i],];
    if(NonZeroCount[i]==1)
    {
      XYZInd = Coord[1] + (Coord[2]-1)*Dimn[1] + (Coord[3]-1)*Dimn[1]*Dimn[2];
      Coord=matrix(c(Coord),1,length(Coord))
    }
    if(NonZeroCount[i]>1)
      XYZInd = Coord[,1] + (Coord[,2]-1)*Dimn[1] + (Coord[,3]-1)*Dimn[1]*Dimn[2];
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
      m=siv(TRImgdata,Coord,HSeq);
      if(m[[2]]==0)
      {
        SRImgdata[,Tid] = TRImgdata;
      }
      else SRImgdata[,Tid] = m[[1]];
    }
  }
  SRImgdata;
}