Simtest<-function(n, pattern, pr, nf1,ef, f1,f2)
{
  nSim = 100;
  ch = 1.1;
  S = 10;
  dimn = c(64,64,8);
  S1 = 5;
  
  a=findVoxelSequence(ch, S, dimn);
  vxlNbySeqid=a[[1]];
  vxlDistSeq=a[[2]];
  vxlNberSeq=a[[3]];
  
  b=findVoxelSequence(ch, S1, dimn);
  vxlNbySeqid1=b[[1]];
  vxlDistSeq1=b[[2]];
  vxlNberSeq1=b[[3]];
  
  simResults =array(rep(0,180),c(3,4,5,3));
  simResults1 =array(rep(0,60),c(2,2,5,3));
  tempBias = tempRms2=tempRms1 =tempSd=array(rep(0,64*64*8*3*3),c(64*64*8,3,3));
  tempTrue = matrix(rep(0,64*64*8*3),ncol=3)
  tempEs =array(rep(0,64*64*8*3*2),c(64*64*8,3,2));
  m= generateDesign(pattern,pr);
  beta=m[[6]];

  for (ii in 1:3)
  {
  temp1 = beta[,,ii];
  temp1 = c(temp1);
  temp1 = kronecker(rep(1,8),temp1);
  tempTrue[,ii] = temp1;
  } 
  
  Mask = array(rep(1,64*64*8),c(64,64,8));
  lt=list(a=0,b=0,c=0,d=0,e=0,f=0);
  
  for (iSum in 1:nSim) 
  {
    xMatrix =matrix(rep(1,n*3),n);
    xMatrix[,2]=rnorm(n,1,0.5);
    xMatrix[,3]=runif(n,1,2);
    xMatrix[,2] = (xMatrix[,2]-mean(xMatrix[,2]))/sd(xMatrix[,2]);
    xMatrix[,3] = (xMatrix[,3]-mean(xMatrix[,3]))/sd(xMatrix[,3]);
    #invOmegaX =solve(t(xMatrix)%*%xMatrix);
    
    
    imgData=generateImageData(beta,xMatrix,nf1,ef);
    
    
    s=mlrWresd(xMatrix,t(imgData));
    mxBeta=s[[1]];
    mxCovb=s[[2]];
    mxR=s[[3]];
    
    smxR=smoothResidual(mxR);
    
    h=smoothParameter(vxlNbySeqid1,vxlDistSeq1, vxlNberSeq1, mxR,smxR,mxBeta,xMatrix,S1,f1,f2);
    mxAspBeta1=h[[1]];
    mxAspCovb1=h[[2]];
    
    
    m1=smoothParameter(vxlNbySeqid,vxlDistSeq, vxlNberSeq, mxR,smxR,mxBeta,xMatrix,S,f1,f2);
    mxAspBeta=m1[[1]];
    mxAspCovb=m1[[2]];
    
    pSn = 2*(1-pt(abs(mxAspBeta/mxAspCovb^(0.5)),dim(xMatrix)[1]-dim(xMatrix)[2]))<0.05;
    pSn0=2*(1-pt(abs(mxBeta/mxCovb^(0.5)),dim(xMatrix)[1]-dim(xMatrix)[2]))<0.05;
    pSn0=t(pSn0);
    
    
    
    
    
    lt[[1]]= mxBeta;
    lt[[2]]= mxCovb;
    lt[[3]]= mxAspBeta1;
    lt[[4]]= mxAspCovb1;
    lt[[5]]= mxAspBeta;
    lt[[6]]= mxAspCovb;
    
    for(ii in 1:3)
    {
      temp=lt[[2*ii-1]];
      if(dim(temp)[1]<dim(temp)[2])
        temp=t(temp);
      tempBias[,,ii] = tempBias[,,ii] + temp - tempTrue;
      tempRms1[,,ii] = tempRms1[,,ii] + temp;
      tempRms2[,,ii] = tempRms2[,,ii] + temp^2; 
      temp1=lt[[2*ii]];
      if(dim(temp1)[1]<dim(temp1)[2])
        temp1=t(temp1);
      tempSd[,,ii] = tempSd[,,ii] + temp1^(0.5);
      
    }
    tempEs[,,1] = tempEs[,,1] + pSn0;
    tempEs[,,2] = tempEs[,,2] + pSn;
    
  }
  
  
  tempBias = tempBias/nSim;
  tempRms = ((tempRms2-tempRms1^2/nSim)/nSim)^(0.5);
  tempSd = tempSd/nSim;
  tempRe = tempRms/tempSd;
  
  tempEs = tempEs/nSim;
  
  for(ii in 1:5)
  {
    
    tempId = m[[ii]];
    if(pattern==7)
    {
      jj=1;
      simResults[,2,ii,jj] = apply(drop(tempRms[,jj,]),2,mean);
      simResults[,3,ii,jj] = apply(drop(tempSd[,jj,]),2,mean);
      simResults[,4,ii,jj] = apply(drop(tempRe[,jj,]),2,mean);
      
      simResults1[,1,ii,jj] = apply(drop(tempEs[,jj,]),2,mean);
      simResults1[,2,ii,jj] = apply(drop(tempEs[,jj,]),2,sd);
      for(j in 2:3)
      {
        tempId1=c();
        for(kk in 1:8)
          tempId1=cbind(tempId1,tempId[,jj-1]+64*64*(kk-1));
        tempId1 = c(tempId1);
        simResults[,1,ii,jj] = apply(drop(tempBias[tempId1,jj,]),2,mean);
        simResults[,2,ii,jj] = apply(drop(tempRms[tempId1,jj,]),2,mean);
        simResults[,3,ii,jj] = apply(drop(tempSd[tempId1,jj,]),2,mean);
        simResults[,4,ii,jj] = apply(drop(tempRe[tempId1,jj,]),2,mean);
        
        simResults1[,1,ii,jj] = apply(drop(tempEs[tempId1,jj,]),2,mean);
        simResults1[,2,ii,jj] = apply(drop(tempEs[tempId1,jj,]),2,sd);
        
      }
    }
    else
    {
      for(jj in 1:3)
      {
        tempId1 =c();
        for(kk in 1:8)
          tempId1=cbind(tempId1,tempId[,jj]+64*64*(kk-1));
        tempId1 = c(tempId1);
        simResults[,1,ii,jj] = apply(drop(tempBias[tempId1,jj,]),2,mean);
        simResults[,2,ii,jj] = apply(drop(tempRms[tempId1,jj,]),2,mean);
        simResults[,3,ii,jj] = apply(drop(tempSd[tempId1,jj,]),2,mean);
        simResults[,4,ii,jj] = apply(drop(tempRe[tempId1,jj,]),2,mean);
        
        simResults1[,1,ii,jj] = apply(drop(tempEs[tempId1,jj,]),2,mean);
        simResults1[,2,ii,jj] = apply(drop(tempEs[tempId1,jj,]),2,sd);
        
      }
      
    }
          
  }
  write.table(tempBias,file="Bias");
  write.table(tempRms,file="Rms");
  write.table(tempSd,file="Sd");
  write.table(tempRe,file="Re");
  write.table(tempEs,file="Es");
  write.table(simResults,file="Results");
  write.table(simResults1,file="Results1");
  
}