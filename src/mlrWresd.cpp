#include <Rcpp.h>
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

// [[Rcpp::export]]
List mlrWresd(NumericMatrix X,NumericMatrix Y){
  
  
  int n, p, N,ii=0,ip;
  int i, j, k, l,imax=0,i1,j1;
  double sigma2,tempR,big,dum,sum,temp,sum1;
  
  n = X.nrow(); 
  p = X.ncol(); 
  N = Y.ncol(); 
  NumericMatrix Beta(p,N);
  NumericMatrix Covb(p,N);
  NumericMatrix R(n,N);
  NumericMatrix XX(p,p);
  NumericMatrix IXX(p,p);
  IntegerVector indx(p);
  NumericVector col(p);
  NumericMatrix PX(p,n);
  NumericMatrix HX(n,n);
  NumericVector vv(p);
  List result(3);
  
  
  
  for (i=0;i<p;i++){
    for (j=0;j<p;j++){
      XX(i,j) = 0.0; 
      for (k=0;k<n;k++){
        XX(i,j) += X(k,i)*X(k,j);
      }
    }
  } 
  
  for (i=0;i<p;i++) {
    big=0.0; 
    for (j=0;j<p;j++)
      if ((temp=fabs(XX(i,j))) > big) 
      {
        big=temp;
      }
    vv[i]=1.0/big; 
  }
  for (j=0;j<p;j++) { 
    for (i=0;i<j;i++) { 
      sum=XX(i,j);
      for (k=0;k<i;k++) sum -= XX(i,k)*XX(k,j);
      XX(i,j)=sum;
    }
    big=0.0; 
    for (i=j;i<p;i++) { 
      sum=XX(i,j); 
      for (k=0;k<j;k++)
        sum -= XX(i,k)*XX(k,j);
      XX(i,j)=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) { 
      for (k=0;k<p;k++) { 
        dum=XX(imax,k);
        XX(imax,k)=XX(j,k);
        XX(j,k)=dum;
      }
      vv[imax]=vv[j]; 
    }
    indx[j]=imax;
    if (XX(j,j) == 0.0) XX(j,j)=TINY;
    if (j != p-1) { 
      dum=1.0/(XX(j,j));
      for (i=j+1;i<p;i++) XX(i,j) *= dum;
    }
  } 
  
  
  for (j=0;j<p;j++) { 
    for (i=0;i<p;i++) 
      col[i]=0.0;
    col[j]=1.0;
    ii=0;
    for (i1=0;i1<p;i1++) { 
      ip=indx[i1];
      sum1=col[ip];
      col[ip]=col[i1];
      if (ii>0)
        for (j1=ii-1;j1<=i1-1;j1++) sum1 -= XX(i1,j1)*col[j1];
      else if (sum1>0) ii=i1+1; 
      col[i1]=sum1; 
    }
    for (i1=p-1;i1>=0;i1--) { 
      sum1=col[i1];
      for (j1=i1+1;j1<p;j1++) sum1 -= XX(i1,j1)*col[j1];
      col[i1]=sum1/XX(i1,i1); 
    } 
    for(i=0;i<p;i++) 
      IXX(i,j)=col[i];
  }
  
  
  
  
  for (i=0;i<p;i++){
    for (j=0;j<n;j++){
      PX(i,j) = 0.0;
      for (k=0;k<p;k++) PX(i,j) +=IXX(i,k)*X(j,k); 
    }
  }
  
  
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      HX(i,j) = 0.0;
      for (k=0;k<p;k++) HX(i,j) -=X(i,k)*PX(k,j); 
      if (i == j) HX(i,j)++;                                                        
    }
  }
  
  
  for (i=0;i<p;i++){
    for (j=0;j<N;j++){
      Beta(i,j) = 0.0;
      for (k=0;k<n;k++) Beta(i,j) += PX(i,k)*Y(k,j);
    }
  }
  
  for (i=0;i<N;i++){
    sigma2 = 0.0;
    for (j=0;j<n;j++){
      for (k=0;k<n;k++){   
        sigma2 += Y(j,i)*HX(j,k)*Y(k,i);
      }
    }
    sigma2 = sigma2/(n-p);
    for (l=0;l<p;l++) Covb(l,i) = sigma2*IXX(l,l);           
  }
  
  for (i=0;i<N;i++){
    for (j=0;j<n;j++){
      tempR = 0.0;
      for (k=0;k<n;k++){   
        tempR += HX(j,k)*Y(k,i);
      }
      R(j,i) = tempR;
    }
  }
  result[0]=Beta;
  result[1]=Covb;
  result[2]=R;
  
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


