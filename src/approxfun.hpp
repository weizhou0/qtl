#ifndef APPROXFUN_HPP
#define APPROXFUN_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// 2021-01-27 (by Wenjian Bi): Most of the below codes are from R::stats::approxfun
namespace approxfun{

class approxfunClass
{
private:
  
  arma::vec m_xVec, m_yVec;
  double m_ylow, m_yhigh;
  int m_n;
  arma::vec m_slopeVec;
  
public:
  
  void setApproxFun(arma::vec t_xVec,
                    arma::vec t_yVec)
  {
    m_xVec = t_xVec;
    m_yVec = t_yVec;
    m_n = t_xVec.size();
    m_ylow = t_yVec(0);
    m_yhigh = t_yVec(m_n - 1);
    m_slopeVec.zeros(m_n - 1);
    
    for(int i = 0; i < m_n - 1; i ++)
      if(t_xVec(i+1) <= t_xVec(i)) Rcpp::stop("xVec(i+1) should be greater than xVec(i).");
    
    for(int i = 0; i < m_n - 1; i ++)
      m_slopeVec(i) = (t_yVec(i+1) - t_yVec(i)) / (t_xVec(i+1) - t_xVec(i));
  }
  
  double getValue(double t_v)
  {
    int i, j, ij;
    i = 0;
    j = m_n - 1;
    
    // handle out-of-domain points
    if(t_v < m_xVec(i)) return m_ylow;
    if(t_v > m_xVec(j)) return m_yhigh;
    
    // find the correct interval by bisection
    while(i < j - 1) { /* x[i] <= v <= x[j] */
      ij = (i + j)/2; /* i+1 <= ij <= j-1 */
      if(t_v < m_xVec(ij)) j = ij; else i = ij;
      /* still i < j */
    }
    
    // interpolation
    if(t_v == m_xVec(j)) return m_yVec(j);
    if(t_v == m_xVec(i)) return m_yVec(i);
    
    // linear interpolation
    return m_yVec(i) + (t_v - m_xVec(i)) * m_slopeVec(i);
  }
  
  arma::vec getVector(arma::vec t_vVec)
  {
    int p = t_vVec.size();
    arma::vec outVec(p);
    for(int i = 0; i < p; i++)
    {
      outVec(i) = getValue(t_vVec(i));
    }
    return outVec;
  }
  
};

}

#endif

// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// 
// #include "approxfun.hpp"
// 
// // [[Rcpp::export]]
// double tempFun(double t_v,
//                arma::vec t_xVec,
//                arma::vec t_yVec) {
//   
//   static approxfun::approxfunClass* ptr_approxfun = NULL;
//   ptr_approxfun = new approxfun::approxfunClass(t_xVec,
//                                                 t_yVec);
//   
//   double out = ptr_approxfun->getValue(t_v);
//   return out;
// }
// 
// /*** R
// xVec = 1:10
// yVec = rnorm(10)
// v = 2.5
// yVec[2:3]
// tempFun(v, xVec, yVec)
// */

