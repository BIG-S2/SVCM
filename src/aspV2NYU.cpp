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
List aspV2NYU(NumericMatrix BetaIn,NumericMatrix EtaIn,NumericVector HsIn,NumericVector VnsIn,NumericMatrix VnsidIn,NumericMatrix VdsIn,NumericMatrix XIn, NumericVector ChiseqIn,NumericVector SeIn,double CnIn,NumericMatrix IXX)
  
  
  
{
  int n, N, p, nh, m, nh1,dh,imax=0;
  int i, j, k, l, ii, jj,i1,j1,ip,ii1=0;
  double tSigma2, dBeta, Kloc, Kst,tDist,big,dum,sum,temp,sum1;
  int nbycount, tId1, tId2;
  
  n = EtaIn.ncol(); 
  N = EtaIn.nrow(); 
  nh = HsIn.size(); 
  nh1 = ChiseqIn.size(); 
  p = BetaIn.ncol(); 
  m = VnsidIn.ncol();         
  NumericMatrix BetaOut(N,p);
  NumericMatrix HsOut(N,p);
  NumericVector col(p);
  IntegerVector indx(p);
  NumericMatrix BetaFix(N,p);
  NumericMatrix Sigma2Fix(N,p);
  NumericMatrix BetaLast(N,p);
  NumericMatrix Sigma2Last(N,p);
  NumericMatrix Sigma2Out(N,p);
  NumericMatrix wBeta(m,p);
  NumericVector sumwBeta(p);
  IntegerMatrix stopUpdate(N,p);
  NumericVector vv(p);
  List result(3);
  
 
  dh = nh-nh1;
  
  
  for (i=0;i<dh;i++){
    
    if (i==0)
    {
      for (j=0;j<N;j++)
      {
        tSigma2 = 0.0;
        for (k=0;k<n;k++)
        {
          tSigma2 += EtaIn(j,k)*EtaIn(j,k);
        }
        tSigma2 = tSigma2/(n-p) + SeIn[j];
        for (k=0;k<p;k++)
        {
          BetaLast(j,k) = BetaIn(j,k);
          Sigma2Last(j,k) = tSigma2*IXX(k,k);                     
        }
      }
    } 
    else 
    {
      for (j=0;j<N;j++)
      {
        for (k=0;k<p;k++)
        {
          BetaLast(j,k) = BetaOut(j,k);
          Sigma2Last(j,k) = Sigma2Out(j,k);                     
        }
      }
    }
    
    for (j=0;j<N;j++){
      
      nbycount = 0;
      for (k=0;k<VnsIn[j];k++)
      {                 
        if (VdsIn(j,nbycount) <= HsIn[i])
        {
          nbycount++;
        }
      }
      
      for (k=0;k<p;k++)
      {
        sumwBeta[k] = 0.0;
      }
      
      for (k=0;k<nbycount;k++)
      {
        Kloc = 1-VdsIn(j,k)/HsIn[i]; 
        tId1 = VnsidIn(j,k);
        for (l=0;l<p;l++)
        {
          dBeta = (BetaLast(j,l)-BetaLast(tId1-1,l))*(BetaLast(j,l)-BetaLast(tId1-1,l))/Sigma2Last(j,l)/(CnIn/(10-i));  
          Kst = exp(-dBeta);
          wBeta(k,l) = Kloc*Kst;  
          sumwBeta[l] = sumwBeta[l]+wBeta(k,l);  
        }
      }             
      
      for (k=0;k<p;k++)
      {
        BetaOut(j,k) = 0.0;
        Sigma2Out(j,k) = 0.0;
        for (l=0;l<nbycount;l++)
        {
          tId1 = VnsidIn(j,l);
          BetaOut(j,k) = BetaOut(j,k) + wBeta(l,k)*BetaLast(tId1-1,k);
          for (ii=0;ii<nbycount;ii++)
          {
            tId2 = VnsidIn(j,ii);
            tSigma2 = 0.0;
            for (jj=0;jj<n;jj++)
            {
              tSigma2 = tSigma2 + EtaIn(tId1-1,jj)*EtaIn(tId2-1,jj);
            }
            tSigma2 = tSigma2/(n-p)*IXX(k,k);
            if (l==ii)
            {
              tSigma2 = tSigma2 + SeIn[l]*IXX(k,k);
            }    
            Sigma2Out(j,k) = Sigma2Out(j,k) + wBeta(l,k)*wBeta(ii,k)*tSigma2;
          }
        }
        BetaOut(j,k) = BetaOut(j,k)/sumwBeta[k];
        Sigma2Out(j,k) = Sigma2Out(j,k)/(sumwBeta[k]*sumwBeta[k]);
      } 
      
    }
    
  }  
  
  
  for (i=0;i<N;i++)
  {
    for (j=0;j<p;j++)
    {
      BetaFix(i,j) = BetaOut(i,j);
      Sigma2Fix(i,j)= Sigma2Out(i,j);
      stopUpdate(i,j) = 0;
    }
  }   
  
  
  for (i=dh;i<nh;i++)
  {
    for (j=0;j<N;j++)
    {
      for (k=0;k<p;k++)
      {
        if (stopUpdate(j,k)==0) 
        {    
          BetaLast(j,k) = BetaOut(j,k);
          Sigma2Last(j,k) = Sigma2Out(j,k);     
        }
      }
    }
    
    for (j=0;j<N;j++){  
      nbycount = 0;
      for (k=0;k<VnsIn[j];k++)
      {                 
        if (VdsIn(j,nbycount) <= HsIn[i])
        {
          nbycount++;
        }
      }
      
      for (k=0;k<p;k++)
      {
        sumwBeta[k] = 0.0;
      }
      
      for (k=0;k<nbycount;k++)
      {
        Kloc = 1-VdsIn(j,k)/HsIn[i];
        tId1 = VnsidIn(j,k);
        for (l=0;l<p;l++)
        {
          dBeta = (BetaLast(j,l)-BetaLast(tId1-1,l))*(BetaLast(j,l)-BetaLast(tId1-1,l))/Sigma2Last(j,l)/(CnIn/(10-i));  
          Kst = exp(-dBeta);                      
          
          wBeta(k,l) = Kloc*Kst;  
          sumwBeta[l] = sumwBeta[l]+wBeta(k,l);  
        }
      }             
      
      for (k=0;k<p;k++)
      {
        if (stopUpdate(j,k)==0)
        {
          BetaOut(j,k) = 0.0;
          Sigma2Out(j,k) = 0.0;
          for (l=0;l<nbycount;l++)
          {
            tId1 = VnsidIn(j,l);
            BetaOut(j,k) = BetaOut(j,k) + wBeta(l,k)*BetaLast(tId1-1,k);
            for (ii=0;ii<nbycount;ii++)
            {
              tId2 = VnsidIn(j,ii);
              tSigma2 = 0.0;
              for (jj=0;jj<n;jj++)
              {
                tSigma2 = tSigma2 + EtaIn(tId1-1,jj)*EtaIn(tId2-1,jj);
              }
              tSigma2 = tSigma2/(n-p)*IXX(k,k);
              if (l==ii)
              {
                tSigma2 = tSigma2 + SeIn[l]*IXX(k,k);
              } 
              Sigma2Out(j,k) = Sigma2Out(j,k) + wBeta(l,k)*wBeta(ii,k)*tSigma2;
            }
          }
          BetaOut(j,k) = BetaOut(j,k)/sumwBeta[k];
          Sigma2Out(j,k) = Sigma2Out(j,k)/(sumwBeta[k]*sumwBeta[k]);
        }
      }
      
    }
      
      for (j=0;j<N;j++)
      {
        for (k=0;k<p;k++)
        {
          if (stopUpdate(j,k)==0)
          {
            tDist = (BetaFix(j,k)-BetaOut(j,k))*(BetaFix(j,k)-BetaOut(j,k))/Sigma2Fix(j,k);
            if (tDist > ChiseqIn[i-3])
            {
              stopUpdate(j,k) = 1;
              HsOut(j,k) = HsIn[i];
            }
          }                     
        }
      }
      
      
    }
    
    
    
    
    for (j=0;j<N;j++)
    {
      for (k=0;k<p;k++)
      {
        if (stopUpdate(j,k)==0)
        {
          stopUpdate(j,k) = 1;
          HsOut(j,k) = HsIn[nh];
        }
      }                     
    }
    
    for (j=0;j<N;j++)
    {
      for (k=0;k<p;k++)
      {
        if (stopUpdate(j,k)==1) 
        {    
          BetaOut(j,k) = BetaLast(j,k);
          Sigma2Out(j,k) = Sigma2Last(j,k);     
        }
      }
    }
  
  result(0)=BetaOut;
  result(1)=HsOut;
  result(2)=Sigma2Out;
  return result;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


