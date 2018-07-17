generateImageData<-function(beta, xMatrix,nf1,ef)
{
  
  n = dim(xMatrix)[1];
  imgMatrixT = array(rep(0,64*64*n),c(n,64,64));
  for (i in 1:n)
  {
  imgMatrixT[i,,] = xMatrix[i,1]*beta[,,1] + xMatrix[i,2]*beta[,,2] + xMatrix[i,3]*beta[,,3];
  }
  d1 = 1:64;
  d2 = 1:64;
  d3 = 1:8;
  phi1 = sqrt(1/4)*sin(2*pi*d1/64);
  phi2 = sqrt(1/4)*cos(2*pi*d2/64);
  phi3 = sqrt(1/2.625)*(9/8-d3/4);
  lambda1 = sqrt(.3*nf1);
  lambda2 = sqrt(.15*nf1);
  lambda3 = sqrt(.05*nf1);
  
  sigma = sqrt(1);
  imgData = matrix(rep(0,64*64*8*n),nrow=64*64*8);
  
  if(ef==1)
  {
  xiMatrix = matrix(rnorm(n*3),n,3);
  }
  if (ef==2)
  {
  xiMatrix = matrix(rnorm(n*3),n,3)^2 + matrix(rnorm(n*3),n,3)^2 + matrix(rnorm(n*3),n,3)^2 -3;
  }
  EE=svd(cov(xiMatrix))$d;
  aa=svd(cov(xiMatrix))$u
  for (i in 1:3)
  {
  xiMatrix[,i]=xiMatrix[,i]-mean(xiMatrix[,i]);
  }
  DD=xiMatrix%*%aa;
  for (i in 1:3)
  {
  xiMatrix[,i] = DD[,i]/sqrt(EE[i]);
  }
  for (i in 1:n)
  {
    imgDataT = array(rnorm(64*64*8,0,sigma),c(64,64,8));
    
    eta1 = xiMatrix[i,1]*phi1*lambda1;
    eta2 = xiMatrix[i,2]*phi2*lambda2;
    eta3 = xiMatrix[i,3]*phi3*lambda3;
    for (j in 1:64)
    {
      for (k in 1:64)
        for (l in 1:8)
        {
          imgDataT[j,k,l] = imgMatrixT[i,j,k]+eta1[j]+eta2[k]+eta3[l]+imgDataT[j,k,l];
        }
    }
    imgData[,i] = c(imgDataT); 
  }
  imgData ;
}