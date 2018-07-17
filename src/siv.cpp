#include <Rcpp.h>
#include <math.h>
#define TINY 1.0e-30 
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List siv(NumericMatrix RImg,NumericMatrix XYZ,NumericVector H){
  
  double minh;
  int n, N, NH;
  int i, j, k, l, ii,i1,j1,k1,i2,j2,imax=0,ip,ii1=0;
  double trS, canh, tkzh,minGCV, currentGCV,big,dum,sum,temp,sum1;
  
  
  n = RImg.nrow(); 
  N = RImg.ncol(); 
  NH = H.size();  
  
  NumericMatrix SRImg(n,N);
  NumericVector GCV(NH);
  NumericMatrix Zh(4,N);
  NumericMatrix KZh(4,N);
  NumericMatrix KZh2(4,4);
  NumericMatrix S(N,N);
  NumericMatrix ISR(N,n);
  NumericMatrix InvKZh2(4,4);
  NumericVector col(4);
  IntegerVector indx(4);
  List result(3);
  NumericVector vv(4);
  
  
  minGCV=10000000000000000.0;
  minh=0.0;
  for (i=0;i<NH;i++){
    canh = H[i];
    
    for (j=0;j<N;j++){
      for (k=0;k<N;k++){
        Zh(0,k)=1;
        Zh(1,k)=(XYZ(k,0)-XYZ(j,0))/canh;
        Zh(2,k)=(XYZ(k,1)-XYZ(j,1))/canh;
        Zh(3,k)=(XYZ(k,2)-XYZ(j,2))/canh;
        tkzh=exp(-0.5*(Zh(1,k)*Zh(1,k)+Zh(2,k)*Zh(2,k)+Zh(3,k)*Zh(3,k)))/canh;
        KZh(0,k)=tkzh;
        KZh(1,k)=tkzh*Zh(1,k);
        KZh(2,k)=tkzh*Zh(2,k);
        KZh(3,k)=tkzh*Zh(3,k);                 
      }         
      
      
      
      for (k=0;k<4;k++){
        for (l=0;l<4;l++){
          KZh2(k,l)=0.0;
          for (ii=0;ii<N;ii++){
            KZh2(k,l)+=KZh(k,ii)*Zh(l,ii);
          }    
        }
      }
      
      for (i1=0;i1<4;i1++) {
        big=0.0; 
        for (j1=0;j1<4;j1++)
          if ((temp=fabs(KZh2(i1,j1))) > big) 
          {
            big=temp;
          }
          vv[i1]=1.0/big; 
      }
      for (j1=0;j1<4;j1++) { 
        for (i1=0;i1<j1;i1++) { 
          sum=KZh2(i1,j1);
          for (k1=0;k1<i1;k1++) sum -= KZh2(i1,k1)*KZh2(k1,j1);
          KZh2(i1,j1)=sum;
        }
        big=0.0; 
        for (i1=j1;i1<4;i1++) { 
          sum=KZh2(i1,j1); 
          for (k1=0;k1<j1;k1++)
            sum -= KZh2(i1,k1)*KZh2(k1,j1);
          KZh2(i1,j1)=sum;
          if ( (dum=vv[i1]*fabs(sum)) >= big) {
            big=dum;
            imax=i1;
          }
        }
        if (j1 != imax) { 
          for (k1=0;k1<4;k1++) { 
            dum=KZh2(imax,k1);
            KZh2(imax,k1)=KZh2(j1,k1);
            KZh2(j1,k1)=dum;
          }
          vv[imax]=vv[j1]; 
        }
        indx[j1]=imax;
        if (KZh2(j1,j1) == 0.0) KZh2(j1,j1)=TINY;
        if (j1 != 3) { 
          dum=1.0/(KZh2(j1,j1));
          for (i1=j1+1;i1<4;i1++) KZh2(i1,j1) *= dum;
        }
      } 
      
      for (k=0;k<4;k++) { 
        for (l=0;l<4;l++) {
          col[l]=0.0;
        }
        col[k]=1.0;
        ii1=0;
        for (i2=0;i2<4;i2++) { 
          ip=indx[i2];
          sum1=col[ip];
          col[ip]=col[i2];
          if (ii1>0)
            for (j2=ii1-1;j2<=i2-1;j2++) sum1 -= KZh2(i2,j2)*col[j2];
          else if (sum1>0) ii1=i2+1; 
          col[i2]=sum1; 
        }
        for (i2=3;i2>=0;i2--) { 
          sum1=col[i2];
          for (j2=i2+1;j2<4;j2++) sum1 -= KZh2(i2,j2)*col[j2];
          col[i2]=sum1/KZh2(i2,i2); 
        } 
        for(l=0;l<4;l++) {
          InvKZh2(l,k)=col[l];
        }
      }
      
      
      for (k=0;k<N;k++) {
        S(j,k)=0.0;
        for (l=0;l<4;l++) { 
          S(j,k)+=InvKZh2(0,l)*KZh(l,k);
        } 
      }
    }       
    
    
    
    trS=0.0;
    for (j=0;j<N;j++){
      trS+=S(j,j);
    }
    
    for (j=0;j<N;j++){
      for (k=0;k<n;k++) {
        ISR(j,k)=RImg(k,j);
        for (l=0;l<N;l++) {
          ISR(j,k)-=S(j,l)*RImg(k,l);
        }    
      }
    }         
    currentGCV=0.0;
    for (j=0;j<n;j++){
      for (k=0;k<N;k++){
        currentGCV+=ISR(k,j)*ISR(k,j);
      }
    } 
    currentGCV=currentGCV/(1-trS/N)/(1-trS/N);
    
    
    GCV[i]=currentGCV;
    if (currentGCV<minGCV) {
      minGCV=currentGCV;
      minh=canh;
    }
  }
  
  
  
  canh=minh;
  for (j=0;j<N;j++){
    for (k=0;k<N;k++){
      Zh(0,k)=1;
      Zh(1,k)=(XYZ(k,0)-XYZ(j,0))/canh;
      Zh(2,k)=(XYZ(k,1)-XYZ(j,1))/canh;
      Zh(3,k)=(XYZ(k,2)-XYZ(j,2))/canh;
      tkzh=exp(-0.5*(Zh(1,k)*Zh(1,k)+Zh(2,k)*Zh(2,k)+Zh(3,k)*Zh(3,k)))/canh;
      KZh(0,k)=tkzh;
      KZh(1,k)=tkzh*Zh(1,k);
      KZh(2,k)=tkzh*Zh(2,k);
      KZh(3,k)=tkzh*Zh(3,k);                 
    }         
    
    for (k=0;k<4;k++){
      for (l=0;l<4;l++){
        KZh2(k,l)=0.0;
        for (ii=0;ii<N;ii++){
          KZh2(k,l)+=KZh(k,ii)*Zh(l,ii);
        }    
      }
    }
    
    
    for (i=0;i<4;i++) {
      big=0.0; 
      for (j1=0;j1<4;j1++)
        if ((temp=fabs(KZh2(i,j1))) > big) 
        {
          big=temp;
        }
        vv[i]=1.0/big; 
    }
    
    for (j1=0;j1<4;j1++) { 
      for (i=0;i<j1;i++) { 
        sum=KZh2(i,j1);
        for (k=0;k<i;k++) sum -= KZh2(i,k)*KZh2(k,j1);
        KZh2(i,j1)=sum;
      }
      big=0.0; 
      for (i=j1;i<4;i++) { 
        sum=KZh2(i,j1); 
        for (k=0;k<j1;k++)
          sum -= KZh2(i,k)*KZh2(k,j1);
        KZh2(i,j1)=sum;
        if ( (dum=vv[i]*fabs(sum)) >= big) {
          big=dum;
          imax=i;
        }
      }
      if (j1 != imax) { 
        for (k=0;k<4;k++) { 
          dum=KZh2(imax,k);
          KZh2(imax,k)=KZh2(j1,k);
          KZh2(j1,k)=dum;
        }
        vv[imax]=vv[j1]; 
      }
      indx[j1]=imax;
      if (KZh2(j1,j1) == 0.0) KZh2(j1,j1)=TINY;
      if (j1 != 3) { 
        dum=1.0/(KZh2(j1,j1));
        for (i=j1+1;i<4;i++) KZh2(i,j1) *= dum;
      }
    }
    
    
    for (k=0;k<4;k++) { 
      for (l=0;l<4;l++) {
        col[l]=0.0;
      }
      col[k]=1.0;
      ii1=0;
      for (i=0;i<4;i++) { 
        ip=indx[i];
        sum1=col[ip];
        col[ip]=col[i];
        if (ii1>0)
          for (j1=ii1-1;j1<=i-1;j1++) 
            sum1 -= KZh2(i,j1)*col[j1];
        else if (sum1>0) ii1=i+1; 
        col[i]=sum1; 
      }
      for (i=3;i>=0;i--) { 
        sum1=col[i];
        for (j1=i+1;j1<4;j1++) sum1 -= KZh2(i,j1)*col[j1];
        col[i]=sum1/KZh2(i,i); 
      } 
      for(l=0;l<4;l++) {
        InvKZh2(l,k)=col[l];
      }
    }
    
    
    for (k=0;k<N;k++) {
      S(j,k)=0.0;
      for (l=0;l<4;l++) { 
        S(j,k)+=InvKZh2(0,l)*KZh(l,k);
      } 
    }
    
  }       
  
  for (j=0;j<n;j++){
    for (k=0;k<N;k++) {
      SRImg(j,k)=0.0;
      for (l=0;l<N;l++) {
        SRImg(j,k)+=S(k,l)*RImg(j,l);
      }    
    }
  }
  result[0]=SRImg;
  result[1]=minh;
  result[2]=GCV;
  return result;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

