findVoxelSequenceRealdata_para<-function(wmNonzeroId,dimn,cores=2)
{
  wmImgsMask=array(rep(0,prod(dimn)),dimn);
  wmImgsMask[wmNonzeroId]=1;
  maskVector = c(wmImgsMask);
  nId = length(wmNonzeroId); 
  
  ch=1.1;
  S=10;
  
  yy=xx=zz=array(rep(1,prod(dimn)),c(dimn[1],dimn[2],dimn[3]));
  
  for(i in 1:dimn[3])
  {
    xx[,,i]=meshgrid(1:dimn[2], 1:dimn[1])$y;
    yy[,,i]=meshgrid(1:dimn[2], 1:dimn[1])$x;
    zz[,,i]=i*zz[,,i];
  }
  
  xx = xx[wmNonzeroId]; 
  yy = yy[wmNonzeroId];
  zz = zz[wmNonzeroId]; 
  Lx = min(xx);
  Ux = max(xx);
  Ly = min(yy);
  Uy = max(yy);
  Lz = min(zz);
  Uz = max(zz);
  hS = ch^S;
  hS2 = hS^2;
  hSx = floor(hS -0.00001);
  hSy = floor(hS -0.00001);
  hSz = floor(hS -0.00001);
  MnVSeq = 50;
  VSeq =VSeqDist= matrix(rep(0,nId*MnVSeq),ncol=50);
  InVSeq = rep(0,nId);
  
  for(Nidi in 1:nId)
  {
    TVSeq = c();
    TVSeqDist = c();
    tx = xx[Nidi];
    ty = yy[Nidi];
    tz = zz[Nidi];
    
    Ltx = max(tx-hSx,Lx);
    Utx = min(tx+hSx,Ux);
    Lty = max(ty-hSy,Ly);
    Uty = min(ty+hSy,Uy);
    Ltz = max(tz-hSz,Lz);
    Utz = min(tz+hSz,Uz);
    
    Ntx = Utx-Ltx+1;
    Nty = Uty-Lty+1;
    Ntz = Utz-Ltz+1;
    
    Newty=Newtx=Newtz=array(rep(1,Ntz*Ntx*Nty),c(Nty,Ntx,Ntz));
    for(i in 1:Ntz)
    {
      Newtx[,,i]=meshgrid(Lty:Uty,Ltx:Utx)$y;
      Newty[,,i]=meshgrid(Lty:Uty,Ltx:Utx)$x;
      Newtz[,,i]=(i+Ltz-1)*Newtz[,,i];
    }
    Newtx = c(Newtx);
    Newty = c(Newty);
    Newtz = c(Newtz);
    tdist = (tx-Newtx)^2+(ty-Newty)^2+(tz-Newtz)^2;
    tdistid = which(tdist<hS2);
    tdist = tdist[tdistid];
    SeqID = Newtx[tdistid] + (Newty[tdistid]-1)*dimn[1] + (Newtz[tdistid]-1)*dimn[1]*dimn[2];
    tMask = maskVector[SeqID];
    
    for (ii in 1:length(tdistid))
    {
      if (length(tMask)==0)
        break;
      if (tMask[ii]!=0)
      {
        TVSeq = c(TVSeq,SeqID[ii]);
        TVSeqDist = c(TVSeqDist,tdist[ii]);
      }
    }
    TnVSeq = length(TVSeq);
    
    if(MnVSeq < TnVSeq)
    {
      a=matrix(rep(0,nId*(TnVSeq-MnVSeq)),nId)
      VSeq = cbind(VSeq,a);
      VSeqDist = cbind(VSeqDist,a);
      MnVSeq = TnVSeq;
    }
    
    VSeq[Nidi,1:TnVSeq] = TVSeq;
    VSeqDist[Nidi,1:TnVSeq] = TVSeqDist; 
    InVSeq[Nidi] = TnVSeq;
  }
  
  
  
  
  m1=VSeq = VSeq[,1:MnVSeq];
  m2=VSeqDist = VSeqDist[,1:MnVSeq];
  
  
  for (iId in 1:nId)
  {
    a=sort(VSeqDist[iId,1:InVSeq[iId]]);
    b=order(VSeqDist[iId,1:InVSeq[iId]]);
    VSeqDist[iId,1:InVSeq[iId]] = t(a);
    VSeq[iId,1:InVSeq[iId]] = VSeq[iId,b];
  }
  
  VSeqId = matrix(rep(0,dim(VSeq)[1]*dim(VSeq)[2]),dim(VSeq)[1]);

  library(parallel)
  # multicores on Linux
  system.time(
    res1.p <- mclapply(1:nId, FUN = function(x) { match_voxel(x,VSeq,InVSeq,wmNonzeroId)}, mc.cores = cores)
  )
  
  VSeqId=matrix(unlist(res1.p),nrow = nId,byrow = TRUE)
  
  VSeqDist = VSeqDist^0.5;
  
  result=list(VSeqId=0,VSeqDist=0,InVSeq=0,VSeq=0);
  result[[1]]= VSeqId;
  result[[2]]= VSeqDist;
  result[[3]]= InVSeq;
  result[[4]]= VSeq;
  result;
  
}