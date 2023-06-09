// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <cmath>
#include <R.h>
#include "SKATO.hpp"
#include <Rmath.h>
#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>

//double SKAT_Optimal_Integrate_Func_single_Davies_Rcpp(double x, arma::vec& pmin_q, arma::vec& tau, double MuQ, arma::vec& lambda, double VarQ, double VarRemain, arma::vec& r_all);

using namespace Numer;

class IntDavies: public Func
{
private:
    arma::vec & pmin_q;
    arma::vec & tau;
    double MuQ;
    arma::vec & lambda;
    double VarQ;
    double VarRemain;
    arma::vec & r_all;

public:
    IntDavies(arma::vec & pmin_q_, arma::vec & tau_, double MuQ_, arma::vec & lambda_, double  VarQ_, double  VarRemain_, arma::vec & r_all_) : pmin_q(pmin_q_), tau(tau_), MuQ(MuQ_), lambda(lambda_), VarQ(VarQ_), VarRemain(VarRemain_), r_all(r_all_) {} 

    double operator()(const double& x) const
    {
        return SKAT_Optimal_Integrate_Func_single_Davies_Rcpp(x, pmin_q, tau, MuQ, lambda, VarQ, VarRemain, r_all);
    }
};



// [[Rcpp::export]]
Rcpp::List SKAT_META_Optimal_Get_Q_Rcpp(const arma::vec& Score, const arma::vec& r_all) {
  int n_r = r_all.n_elem;
  arma::vec Q_r(n_r, arma::fill::zeros);
  
  for (int i = 0; i < n_r; i++) {
    double r_corr = r_all(i);
    Q_r(i) = (1 - r_corr) * arma::accu(Score % Score) + r_corr * std::pow(arma::accu(Score), 2);
  }
  
  Q_r /= 2.0;
  
  //Rcpp::List re;
  //re["Qr"] = Q_r;
  //return re;
  return Rcpp::List::create(Rcpp::Named("Qr") = Q_r);
}



// [[Rcpp::export]]
Rcpp::List SKAT_META_Optimal_Get_Q_Res_Rcpp(const arma::mat& Score_res, const arma::vec& r_all) {
  int n_r = r_all.n_elem;
  int p = Score_res.n_rows;
  arma::mat Q_r(p, n_r, arma::fill::zeros);
  
  for (int i = 0; i < n_r; i++) {
    double r_corr = r_all(i);
    Q_r.col(i) = (1 - r_corr) * arma::sum(Score_res % Score_res, 1) + r_corr * arma::pow(arma::sum(Score_res, 1), 2);
  }
  
  Q_r /= 2.0;
  
  //Rcpp::List re;
  //re["Qr"] = Q_r;
  return Rcpp::List::create(Rcpp::Named("Qr") = Q_r);
  //return re;
}


// [[Rcpp::export]]
arma::vec Get_Lambda_Rcpp(const arma::mat& K, bool isFast, int maxK) {
  arma::vec lambda;
  
  arma::vec eigenvalues;
  arma::eig_sym(eigenvalues, K);
  
  arma::uvec IDX1 = arma::find(eigenvalues >= 0);
  
  // eigenvalue bigger than sum(eigenvalues)/1000
  arma::uvec IDX2 = arma::find(eigenvalues > arma::mean(eigenvalues(IDX1)) / 100000);
  
  if (IDX2.n_elem == 0) {
    Rcpp::stop("No Eigenvalue is bigger than 0!!");
  }
  
  lambda = eigenvalues(IDX2);
  
  return lambda;
}




// [[Rcpp::export]]
Rcpp::List SKAT_META_Optimal_Param_Rcpp(const arma::mat& Phi, const arma::vec& r_all) {
  int p_m = Phi.n_cols;
  int r_n = r_all.n_elem;
  
  // ZMZ
  arma::mat Z_item1_1 = Phi * arma::ones<arma::vec>(p_m);
  arma::mat ZZ = Phi;
  arma::mat ZMZ = (Z_item1_1 * Z_item1_1.t()) / arma::accu(ZZ);
  
  // W3.2 Term: mixture chisq
  arma::mat W3_2_t = ZZ - ZMZ;
  bool iffast=false;
  int maxK = 100;
  arma::vec lambda = Get_Lambda_Rcpp(W3_2_t, iffast, maxK);
  
  // W3.3 Term: variance of remaining ...
  double W3_3_item = arma::accu(ZMZ % (ZZ - ZMZ)) * 4.0;
  
  // tau term
  double z_mean_2 = arma::accu(ZZ) / std::pow(p_m, 2);
  double tau1 = arma::accu(ZZ * ZZ) / std::pow(p_m, 2) / z_mean_2;
  
  // Mixture Parameters
  double MuQ = arma::accu(lambda);
  double VarQ = arma::accu(arma::square(lambda)) * 2.0 + W3_3_item;
  double KerQ = arma::accu(arma::square(lambda)) / std::pow(arma::accu(lambda), 2) * 12.0;
  double Df = 12.0 / KerQ;
  
  // W3.1 Term: tau1 * chisq_1
  arma::vec tau(r_n, arma::fill::zeros);
  for (int i = 0; i < r_n; i++) {
    double r_corr = r_all(i);
    double term1 = p_m * p_m * r_corr * z_mean_2 + tau1 * (1 - r_corr);
    tau(i) = term1;
  }
  
  Rcpp::List out;
  out["MuQ"] = MuQ;
  out["VarQ"] = VarQ;
  out["KerQ"] = KerQ;
  out["lambda"] = lambda;
  out["VarRemain"] = W3_3_item;
  out["Df"] = Df;
  out["tau"] = tau;
  out["z_mean_2"] = z_mean_2;
  out["p.m"] = p_m;
  out["tau.1"] = tau1;
  out["tau.2"] = p_m * z_mean_2;
  
  return out;
}


// [[Rcpp::export]]
Rcpp::List SKAT_Optimal_Each_Q_Rcpp(const arma::mat& Q_all, const arma::vec& r_all, const Rcpp::List& lambda_all, std::string method) {
  int n_r = r_all.n_elem;
  int n_q = Q_all.n_rows;
  
  arma::vec c1(4, arma::fill::zeros);
  
  arma::mat pval(n_q, n_r, arma::fill::zeros);
  arma::mat pmin_q(n_q, n_r, arma::fill::zeros);
  
  arma::mat param_mat;
  //std::cout << "SKAT_Optimal_Each_Q_Rcpp 1" << std::endl; 
  for (int i = 0; i < n_r; i++) {
    arma::vec Q = Q_all.col(i);
    double r_corr = r_all(i);
    arma::vec lambda_temp = lambda_all[i];
    lambda_temp.print("lambda_temp") ;
    c1[0] = arma::accu(lambda_temp);
    c1[1] = arma::accu(arma::pow(lambda_temp, 2));
    c1[2] = arma::accu(arma::pow(lambda_temp, 3));
    c1[3] = arma::accu(arma::pow(lambda_temp, 4));
    c1.print("c1");
    Rcpp::List param_temp = Get_Liu_Params_Mod_Rcpp(c1);
    
    double muQ = param_temp["muQ"];
    double sigmaQ = param_temp["sigmaQ"];
    double varQ = std::pow(sigmaQ, 2);
    double df = param_temp["l"];
    
    // get p-value
    arma::vec Q_Norm = (Q - muQ) / std::sqrt(varQ) * std::sqrt(2 * df) + df;


  //std::cout << "SKAT_Optimal_Each_Q_Rcpp 2" << std::endl; 
  //Q_Norm.print("Q_Norm");
  //std::cout << "df " << df << std::endl;
    double pval0;
    boost::math::chi_squared_distribution<> dist(df);
    for (int j = 0; j < Q_Norm.n_elem; j++) {
        pval0 = boost::math::cdf(complement(dist, Q_Norm[j])); //lower.tail=FALSE
	pval(j, i) = 1 - pval0;
    }
    arma::vec df1; 
     if (method == "optimal.mod" || method == "optimal.adj" || method == "optimal.moment.adj") {
        Rcpp::List pvalue_lambda = Get_PValue_Lambda_Rcpp(lambda_temp, Q, df1);
	arma::vec pvalvec = pvalue_lambda["p.value"];
        pval.col(i) = pvalvec;
      }
    
    param_mat = arma::join_vert(param_mat, arma::rowvec({muQ, varQ, df}));
  }
  //pval.print("pval");
  //param_mat.print("param_mat");
  //std::cout << "SKAT_Optimal_Each_Q_Rcpp 3" << std::endl; 
  arma::vec pmin = arma::min(pval, 1);
  //pmin.print("pmin");
  for (int i = 0; i < n_r; i++) {
    double muQ = param_mat(i, 0);
    double varQ = param_mat(i, 1);
    double df = param_mat(i, 2);
    
    //std::cout << "SKAT_Optimal_Each_Q_Rcpp 4" << std::endl; 
     //df, nc 
     arma::vec q_org(pmin.n_elem);
     double nc = 0;
     //int df0 = static_cast<int>(df);
     boost::math::chi_squared_distribution<double> chiSqDist(df);
     for (int j = 0; j < pmin.n_elem; j++) {
           q_org[j] = boost::math::quantile(chiSqDist, 1-pmin[j]);
     } 

    //std::cout << "SKAT_Optimal_Each_Q_Rcpp 4a " << i << std::endl;	
    //q_org.print("q_org");

     arma::vec q_q = (q_org - df) / std::sqrt(2 * df) * std::sqrt(varQ) + muQ;
     pmin_q.col(i) = q_q;
  }
  
  //std::cout << "SKAT_Optimal_Each_Q_Rcpp 5" << std::endl; 
   return Rcpp::List::create(
    Rcpp::Named("pmin") = pmin,
    Rcpp::Named("pval") = pval,
    Rcpp::Named("pmin_q") = pmin_q
  ); 
}

// [[Rcpp::export]]
Rcpp::List Get_Liu_Params_Mod_Rcpp(const arma::vec& c1) {
  double muQ = c1(0);
  double sigmaQ = std::sqrt(2 * c1(1));
  double s1 = c1(2) / std::pow(c1(1), 1.5);
  double s2 = c1(3) / std::pow(c1(1), 2.0);
  double beta1 = std::sqrt(8) * s1;
  double beta2 = 12 * s2;
  int type1 = 0;

  double l, a, d, muX, sigmaX;

  if (std::pow(s1, 2) > s2) {
    a = 1 / (s1 - std::sqrt(std::pow(s1, 2) - s2));
    d = s1 * std::pow(a, 3) - std::pow(a, 2);
    l = std::pow(a, 2) - 2 * d;
  } else {
    type1 = 1;
    l = 1 / s2;
    a = std::sqrt(l);
    d = 0;
  }

  muX = l + d;
  sigmaX = std::sqrt(2) * a;

  return Rcpp::List::create(
    Rcpp::Named("l") = l,
    Rcpp::Named("d") = d,
    Rcpp::Named("muQ") = muQ,
    Rcpp::Named("muX") = muX,
    Rcpp::Named("sigmaQ") = sigmaQ,
    Rcpp::Named("sigmaX") = sigmaX
  );

}


// SKAT_davies function implementation
Rcpp::List SKAT_davies_Rcpp(double q, arma::vec& lambda, arma::ivec& h, arma::vec& delta, double sigma, int lim, double acc){


  int r = lambda.n_elem;
  if (h.n_elem != r)
    Rcpp::stop("lambda and h should have the same length!");
  if (delta.n_elem != r)
    Rcpp::stop("lambda and delta should have the same length!");

  arma::vec trace(7, arma::fill::zeros);
  int ifault = 0;
  double result = 0;
  std::vector<int> hi = arma::conv_to < std::vector<int> >::from(h);

  //lim = 10000000000000; 
//double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res
  qfc_1(&lambda[0], &delta[0], &hi[0], &r, &sigma, &q, &lim, &acc, &trace[0], &ifault, &result);
  result = 1 - result;
  //std::cout << "ifault " << ifault << std::endl;
  return Rcpp::List::create(
    Rcpp::Named("trace") = trace,
    Rcpp::Named("ifault") = ifault,
    Rcpp::Named("Qq") = result);

}


// Function to calculate SKAT optimal p-value using Davies method (DONE)
// [[Rcpp::export]]
arma::vec SKAT_Optimal_Integrate_Func_Davies_Rcpp(const arma::vec& x, const arma::vec& pmin_q, arma::vec& tau, double MuQ, arma::vec& lambda, double VarQ, double VarRemain, const arma::vec& r_all){ 

  int n_r = r_all.n_elem;
  int n_x = x.n_elem;
  
  arma::mat temp1 = arma::kron(tau, x.t());
  
  arma::mat temp = (pmin_q - temp1) / (1 - r_all);
  arma::rowvec temp_min = arma::min(temp, 1);
 
  arma::ivec h(lambda.n_elem);
  h.ones();
  arma::vec delta(lambda.n_elem);
  delta.zeros();
  double sigma = 0;
  int lim = 10000;
  double acc = 0.000001; 
  
 
  arma::vec result(n_x);
  for (int i = 0; i < n_x; ++i) {
    double min1 = temp_min(i);
    double temp_value = 0.0;
    if (min1 <= arma::accu(lambda) * 10e4) {
      double min1_temp = min1 - MuQ;
      double sd1 = std::sqrt(VarQ - VarRemain) / std::sqrt(VarQ);
      double min1_st = min1_temp * sd1 + MuQ;
      
      Rcpp::List dav_re = SKAT_davies_Rcpp(min1_st, lambda, h, delta, sigma, lim, acc);
      temp_value = dav_re["Qq"];
      
      //int ifault0 = dav_re["ifault"];
      //std::cout << "ifault0 " << ifault0 << std::endl;
      /*if (dav_re["ifault"] != 0) {
        Rcpp::stop("dav_re$ifault is not 0");
      }
      */
    }
    
    if (temp_value > 1.0) {
      temp_value = 1.0;
    }
    
    result(i) = (1 - temp_value) *  R::dchisq(x(i), 1, false);
  }
  
  return result;
}



// Function to calculate SKAT optimal p-value using Davies method (DONE)
// [[Rcpp::export]]
double SKAT_Optimal_Integrate_Func_single_Davies_Rcpp(double x, arma::vec& pmin_q, arma::vec& tau, double MuQ, arma::vec& lambda, double VarQ, double VarRemain, arma::vec& r_all){


 // int n_r = r_all.n_elem;
  int n_x = 1;
  arma::vec xvec(1);
  xvec[0] = x;
  arma::mat temp1 = arma::kron(tau, xvec);
  //temp1.print("temp1");
  arma::vec temp = (pmin_q - temp1.col(0)) / (1 - r_all);
  //temp.print("temp");
  double temp_min = temp.min();
  //std::cout << "temp_min " << temp_min << std::endl;
  
  arma::ivec h(lambda.n_elem);
  h.ones();
  arma::vec delta(lambda.n_elem);
  delta.zeros();
  double sigma = 0;
  int lim = 10000;
  double acc = 0.000001;
  arma::vec result(n_x);
  for (int i = 0; i < n_x; ++i) {
    double min1 = temp_min;
    double temp_value = 0.0;
    if (min1 <= arma::accu(lambda) * 10e4) {
      double min1_temp = min1 - MuQ;
      double sd1 = std::sqrt(VarQ - VarRemain) / std::sqrt(VarQ);
      double min1_st = min1_temp * sd1 + MuQ;

      Rcpp::List dav_re = SKAT_davies_Rcpp(min1_st, lambda, h, delta, sigma, lim, acc);
      temp_value = dav_re["Qq"];
      //int ifault0 = dav_re["ifault"];
      //std::cout << "ifault0 " << ifault0 << std::endl;
      //std::cout << "temp " << temp_value << std::endl;
      /*if (ifault0 != 0) {
        Rcpp::stop("dav_re$ifault is not 0");
      }*/
    }else{
         temp_value = 1.0;
    }

    if (temp_value > 1.0) {
      temp_value = 1.0;
    }
    double xi = xvec(i);
    result(i) = (1 - temp_value) *  R::dchisq(xi, 1, false);
  }
  double result_s = result(0);
  return result_s;
}



// Function for calculating SKAT p-value using Davies method
double SKAT_Optimal_PValue_Davies_Rcpp(arma::vec& pmin_q, arma::vec& tau, double MuQ, arma::vec& lambda, double VarQ, double VarRemain, arma::vec& r_all, double pmin) {

  //using namespace boost::math::quadrature;
  int n_r = r_all.n_elem;
  using namespace Numer;

  // Numerical integration using Davies method
  //std::cout << "SKAT_Optimal_PValue_Davies_Rcpp 1" << std::endl;
  //double xtest = 0.5;
  //double testx = SKAT_Optimal_Integrate_Func_single_Davies_Rcpp(xtest,  pmin_q, tau, MuQ, lambda, VarQ, VarRemain, r_all);
  //std::cout << "testx " << testx << std::endl;

  IntDavies f(pmin_q, tau, MuQ, lambda, VarQ, VarRemain, r_all);
  const double lower = 0.0, upper = 40.0;
  double err_est;
  int err_code;
  double integral = Numer::integrate(f, lower, upper, err_est, err_code);

  //auto integrand = [pmin_q, tau, MuQ, lambda, VarQ, VarRemain, r_all](double x) mutable { return SKAT_Optimal_Integrate_Func_single_Davies_Rcpp(x, pmin_q, tau, MuQ, lambda, VarQ, VarRemain, r_all); };  // Lambda function with additional parameters
  //std::cout << "SKAT_Optimal_PValue_Davies_Rcpp 2" << std::endl;

  //double result, error;
  //double integral = gauss<double, 10000>::integrate(integrand, 0.0, 40.0);
  //std::cout << "SKAT_Optimal_PValue_Davies_Rcpp 3" << std::endl;
  //std::cout << "integral " << integral << std::endl;

  double p_value = 1.0 - integral;
  if (pmin != R_NegInf) {
    double threshold = pmin * n_r;
    if (threshold < p_value) {
      p_value = threshold;
    }
  }
  return p_value;
}

/*
// Function for numerical integration using Liu method
arma::vec SKAT_Optimal_Integrate_Func_Liu_Rcpp(const arma::vec& x, const arma::vec& pmin_q, const arma::vec& tau, double MuQ, double VarQ, const arma::vec& Df, const arma::vec& r_all) {
  int n_r = r_all.n_elem;
  int n_x = x.n_elem;

  arma::mat temp1 = arma::kron(tau, x.t());
  arma::mat temp = (pmin_q - temp1) / (1 - r_all);
  arma::vec temp_min = arma::min(temp, 1);

  arma::vec temp_q = (temp_min - MuQ) / std::sqrt(VarQ) * std::sqrt(2 * Df) + Df;
  arma::vec result(n_x);
   for(unsigned int k = 0; k < n_x; k++){
	result[k] = R::pchisq(temp_q[k], Df[k], 0, 1) * R::dchisq(x[k], 1, false);
   }

  return result;
}

// Function for calculating SKAT p-value using Liu method
double SKAT_Optimal_PValue_Liu_Rcpp(const arma::mat& pmin_q, const arma::mat& tau, const arma::vec& MuQ, double VarQ, const arma::vec& Df, const arma::vec& r_all, double pmin = R_NegInf) {
  // Numerical integration using Liu method
  arma::vec integrand = SKAT_Optimal_Integrate_Func_Liu_Rcpp(arma::regspace(0.0, 40.0, 2000), pmin_q, tau, MuQ, VarQ, Df, r_all);
  double integral = arma::accu(integrand) * 40.0 / 2000.0;

  double p_value = 1.0 - integral;
  if (pmin != R_NegInf) {
    double threshold = pmin * r_all.n_elem;
    if (threshold < p_value) {
      p_value = threshold;
    }
  }
  return p_value;
}


// Main function to calculate p-values using SKAT-O method
Rcpp::List SKAT_Optimal_Get_Pvalue_Rcpp(const arma::mat& Q_all, const arma::mat& Z1, const arma::vec& r_all, const std::string& method) {
  int n_r = r_all.n_elem;
  int n_q = Q_all.n_rows;
  int p_m = Z1.n_cols;

  Rcpp::List lambda_all(n_r);
  for (int i = 0; i < n_r; i++) {
    double r_corr = r_all(i);
    arma::mat R_M = arma::diagmat(1 - r_corr, p_m) + arma::repmat(arma::vec(p_m, arma::fill::ones) * r_corr, 1, p_m);
    arma::mat L;
    arma::chol(L, R_M, "lower");
    arma::mat Z2 = Z1 * L.t();
    arma::mat K1 = Z2.t() * Z2;
    lambda_all(i) = Get_Lambda_Rcpp(K1);
  }

  Rcpp::List param_m = SKAT_Optimal_Param_Rcpp(Z1, r_all);
  Rcpp::List each_info = SKAT_Optimal_Each_Q_Rcpp(param_m, Q_all, r_all, lambda_all, method);
  arma::mat pmin_q = each_info["pmin.q"];
  arma::vec pmin = each_info["pmin"];
  arma::vec pval(n_q);

  if (method == "davies" || method == "optimal" || method == "optimal.mod" || method == "optimal.adj") {
    for (int i = 0; i < n_q; i++) {
      pval(i) = SKAT_Optimal_PValue_Davies_Rcpp(pmin_q.row(i), param_m, r_all, pmin(i));
    }
  } else if (method == "liu" || method == "liu.mod" || method == "optimal.moment" || method == "optimal.moment.adj") {
    for (int i = 0; i < n_q; i++) {
      pval(i) = SKAT_Optimal_PValue_Liu_Rcpp(pmin_q.row(i), param_m, r_all, pmin(i));
    }
  } else {
    Rcpp::stop("Invalid Method!");
  }

  int multi = 3;
  if (n_r < 3) {
    multi = 2;
  }

  for (int i = 0; i < n_q; i++) {
    arma::vec pval_each = each_info["pval"].row(i).t();
    arma::uvec idx = arma::find(pval_each > 0);

    double pval1 = arma::min(pval_each) * multi;
    if (pval(i) <= 0 || idx.n_elem < n_r) {
      pval(i) = pval1;
    }

    if (pval(i) == 0) {
      if (idx.n_elem > 0) {
        pval(i) = arma::min(pval_each.elem(idx));
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("p.value") = pval, Rcpp::Named("p.val.each") = each_info["pval"]);
}
*/

// Main function to perform SKAT-META optimal analysis
Rcpp::List SKAT_META_Optimal_Rcpp(const arma::mat& Score, const arma::mat& Phi, arma::vec & r_all, const std::string& method, const arma::mat& Score_Resampling) {
  arma::uvec idx = arma::find(r_all >= 0.999);
  if (idx.n_elem > 0) {
    r_all.elem(idx).fill(0.999);
  }

  int p_m = Phi.n_cols;
  int n_r = r_all.n_elem;
 //std::cout << "Score_Resampling_t1_1" << std::endl;
  Rcpp::List out_Q = SKAT_META_Optimal_Get_Q_Rcpp(Score, r_all);
arma::vec Q_r = out_Q["Qr"];
//Q_r.print("Q_r");
 //std::cout << "Score_Resampling_t1_a" << std::endl;
  arma::mat Q_all;
  arma::mat Q_res;
  if (Score_Resampling.size() > 0) {
    Rcpp::List out_Q_res = SKAT_META_Optimal_Get_Q_Res_Rcpp(Score_Resampling, r_all);
    Rcpp::NumericVector Q_res_vec0 = out_Q_res["Qr"];
    arma::vec Q_res_vec1 = Rcpp::as< arma::vec >(Q_res_vec0);
    //arma::vec Q_res_vec1 = arma::conv_to < arma::vec >::from(Q_res_vec0);
    Q_res.set_size(1, Q_res_vec1.n_elem);
    Q_res.row(0) = Q_res_vec1;	
    Q_all = arma::join_vert(Q_r, Q_res);
  }else{
// std::cout << "Score_Resampling_t1_b" << std::endl;
    Q_all.set_size(1, Q_r.n_elem);
    Q_all.row(0) = Q_r.t();   
  }

 //std::cout << "Score_Resampling_t1_2" << std::endl;
  bool isFast = false;
  arma::mat Phi2 = Phi / 2;
  Rcpp::List out = SKAT_META_Optimal_Get_Pvalue_Rcpp(Q_all, Phi2, r_all, method, isFast);

 //std::cout << "Score_Resampling_t1_3" << std::endl;
  Rcpp::List param;
  arma::mat p_val_each_mat = out["p_val_each"];
  param["p_val_each"] = p_val_each_mat.row(0);
  param["q_val_each"] = Q_all.row(0);
  param["rho"] = r_all;
  arma::vec p_val_each_vec = param["p_val_each"];
  
  param["minp"] = arma::min(p_val_each_vec);

  arma::uvec id_temp = arma::find(p_val_each_vec == param["minp"]);

  arma::vec rho = param["rho"];
  arma::uvec id_temp1 = arma::find(rho >= 0.999);
  if (id_temp1.n_elem > 0) {
    rho.elem(id_temp1).fill(1);
  }

  arma::vec rho_est = rho.elem(id_temp);
  param["rho_est"] = rho_est;
  arma::vec p_value_vec = out["p_value"];
  double p_value = p_value_vec(0);
  arma::vec p_value_resampling;
  if (!Q_res.is_empty()) {
    p_value_resampling = p_value_vec.subvec(1, p_value_vec.n_elem - 1);   
  }

  //std::cout << "p_value " << p_value << std::endl; 
  //p_val_each_vec.print("p_val_each_vec");
  //arma::vec q_val_each_vec = param["q_val_each"];
  //q_val_each_vec.print("q_val_each_vec");
  //r_all.print("rho");

Rcpp::List re;
  re["p_value"] = p_value;
  re["param"] = param;
  re["p_value_resampling"] = p_value_resampling;
   
  //std::cout << "Score_Resampling_t1_4" << std::endl;

  return re;
}


// Function to compute p-values for each variant (DONE)
Rcpp::List SKAT_META_Optimal_Get_Pvalue_Rcpp(const arma::mat& Q_all, const arma::mat& Phi, arma::vec& r_all, const std::string& method, bool isFast) {
  int n_r = r_all.n_elem;
  int n_q = Q_all.n_rows;
  int p_m = Phi.n_cols;

  //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp 0" << std::endl;
  Rcpp::List lambda_all(n_r);
  arma::vec lambdavec;
  //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp " << std::endl;
  //r_all.print("r_all");
for (int i = 0; i < n_r; ++i) {
    double r_corr = r_all[i];
    arma::mat R_M = arma::diagmat(arma::vec(p_m, arma::fill::ones) * (1 - r_corr)) + arma::mat(p_m, p_m, arma::fill::ones) * r_corr;
    R_M.print("R_M");
    arma::mat L = arma::chol(R_M, "upper");
    arma::mat Phi_rho = L * (Phi * L.t());
    int maxK = 100;
    lambdavec.clear();
    //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp ite " << i << std::endl;
    //Phi_rho.print("Phi_rho");
    lambdavec = Get_Lambda_Rcpp(Phi_rho, isFast, maxK);
    Rcpp::NumericVector lambda_vec(lambdavec.begin(), lambdavec.end()); // Convert arma::vec to NumericVector
    Rcpp::CharacterVector element_name = Rcpp::CharacterVector::create("element_" + std::to_string(i + 1));
    lambda_vec.attr("names") = element_name; // Set the names attribute of the NumericVector
    lambda_all[i] = lambda_vec;
  }

  //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp 2" << std::endl;
  Rcpp::List param_m = SKAT_META_Optimal_Param_Rcpp(Phi, r_all);
  Rcpp::List Each_Info = SKAT_Optimal_Each_Q_Rcpp(Q_all, r_all, lambda_all, method);
  //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp 2a" << std::endl;
  arma::mat pmin_q = Each_Info["pmin_q"];
  //pmin_q.print("pmin_q");
  //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp 2b" << std::endl;
  arma::mat pval(n_q, 1, arma::fill::zeros);

  // added
  arma::vec pmin = Each_Info["pmin"];
  //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp 2c" << std::endl;

  if (method == "davies" || method == "optimal" || method == "optimal.mod" || method=="optimal.adj") {
    for (int i = 0; i < n_q; ++i) {
	arma::vec pmin_q_vec = pmin_q.row(i).t();
	arma::vec tau =  param_m["tau"];
	double muQ = param_m["MuQ"];
	arma::vec lambda = param_m["lambda"];
	double VarQ = param_m["VarQ"];
	double VarRemain = param_m["VarRemain"];
	double pminval =  pmin(i);
  //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp 2d" << std::endl;
	pval(i) = SKAT_Optimal_PValue_Davies_Rcpp(pmin_q_vec, tau, muQ, lambda, VarQ, VarRemain, r_all, pminval);

  //std::cout << "SKAT_META_Optimal_Get_Pvalue_Rcpp 2e" << std::endl;

	//pval(i) = SKAT_Optimal_PValue_Davies_Rcpp(pmin_q_vec, param_m["tau"], param_m["muQ"], param_m["lambda"], param_m["VarQ"], param_m["VarRemain"],  r_all, pmin(i));
    }
  }else {
    Rcpp::stop("Invalid Method:", method);
  }
  /*
  else if (method == "liu" || method == "liu.mod") {
    for (int i = 0; i < n_q; ++i) {
      pval(i) = SKAT_Optimal_PValue_Liu_Rcpp(pmin_q.row(i).t(), param_m, r_all, pmin(i));
    }
  }
  */
  

//pval.print("pval");


  // Check the pval
  // Since SKAT-O is between burden and SKAT, SKAT-O p-value should be <= min(p-values) * 2
  // To correct conservatively, we use min(p-values) * 3
  int multi = 3;
  if (n_r < 3) {
    multi = 2;
  }
  arma::mat pvalmat = Each_Info["pval"];
  for (int i = 0; i < n_q; ++i) {
    arma::rowvec pval_each = pvalmat.row(i);
    arma::uvec idx = arma::find(pval_each > 0);
    double pval_1 = arma::min(pval_each) * multi;
    if (pval(i) <= 0 || idx.n_elem < n_r) {
      pval(i) = pval_1;
    }

    // if pval == 0, use nonzero min each.pval as p-value
    if (pval(i) == 0) {
      if (idx.n_elem > 0) {
        pval(i) = arma::min(pval_each(idx));
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("p_value") = pval, Rcpp::Named("p_val_each") = Each_Info["pval"]);
}




// [[Rcpp::export]]
Rcpp::List Get_Liu_Params_Mod_Lambda_Rcpp(const arma::vec& lambda, arma::vec& df1) {
  int n = lambda.n_elem;;

  if (df1.n_elem == 0) {
    df1.set_size(n);
    df1.ones();
  }

  arma::vec c1(4, arma::fill::zeros);
  for (int i = 0; i < 4; i++) {
    c1[i] = arma::accu(arma::pow(lambda, i + 1) % df1);
  }

  double muQ = c1[0];
  double sigmaQ = std::sqrt(2 * c1[1]);
  double s1 = c1[2] / std::pow(c1[1], 1.5);
  double s2 = c1[3] / std::pow(c1[1], 2);

  double beta1 = std::sqrt(8) * s1;
  double beta2 = 12 * s2;
  int type1 = 0;

  double a, d, l;
  if (std::pow(s1, 2) > s2) {
    a = 1 / (s1 - std::sqrt(std::pow(s1, 2) - s2));
    d = s1 * std::pow(a, 3) - std::pow(a, 2);
    l = std::pow(a, 2) - 2 * d;
  } else {
    type1 = 1;
    l = 1 / s2;
    a = std::sqrt(l);
    d = 0;
  }

  double muX = l + d;
  double sigmaX = std::sqrt(2) * a;



return Rcpp::List::create(Rcpp::Named("l") = l,
                      Rcpp::Named("d") = d,
                      Rcpp::Named("muQ") = muQ,
                      Rcpp::Named("muX") = muX,
                      Rcpp::Named("sigmaQ") = sigmaQ,
                      Rcpp::Named("sigmaX") = sigmaX);
}




// Compute p-value for modified Liu method with known parameters
// Q: test statistic
// muQ, muX, sigmaQ, sigmaX, l, d: parameters
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
std::string Get_Liu_PVal_MOD_Lambda_Zero_Rcpp(double Q, double muQ, double muX, double sigmaQ, double sigmaX, double l, double d) {
  double QNorm = (Q - muQ) / sigmaQ;
  double QNorm1 = QNorm * sigmaX + muX;

  arma::vec temp = {0.05, pow(10, -10), pow(10, -20), pow(10, -30), pow(10, -40), pow(10, -50), pow(10, -60), pow(10, -70), pow(10, -80), pow(10, -90), pow(10, -100)};

   double quantile;
   arma::vec out;
   boost::math::non_central_chi_squared_distribution<double> chiSqDist(l, d);
   for (int i = 0; i < temp.n_elem; i++) {
     quantile = boost::math::quantile(chiSqDist, 1-temp[i]);
     //lower.tail=F
     //    double lowerTailProbability = boost::math::cdf(boost::math::chi_squared(degreesOfFreedom), x);
     out[i] = quantile;
   }

  int idx = arma::max(arma::find(out < QNorm1));

  std::string pvalMsg = "Pvalue < " + std::to_string(temp[idx]);
  return pvalMsg;
}



// Function to calculate the p-value using arma::vec and pow()
// [[Rcpp::export]]
arma::vec Get_Liu_PVal_MOD_Lambda_Rcpp(const arma::vec& Q_all, const arma::vec& lambda, arma::vec& df1, bool log_p) {
  int n = Q_all.n_elem;
  arma::vec p_value(n);

    // Calculate parameters
    
	lambda.print("lambda");
	df1.print("df1");
    Rcpp::List  param = Get_Liu_Params_Mod_Lambda_Rcpp(lambda, df1);
    double muQ = param["muQ"];
    double muX = param["muX"];
    double sigmaQ = param["sigmaQ"];
    double sigmaX = param["sigmaX"];
    double l = param["l"];
    double d = param["d"];
    //std::cout << muQ << " " << muX << " " << sigmaQ << " " << sigmaX << " " << l << " " << d << std::endl;
  //Q_all.print("Q_all");
  for (int i = 0; i < n; i++) {
    // Compute normalized values
    double Q_Norm = (Q_all(i) - muQ) / sigmaQ;
    double Q_Norm1 = Q_Norm * sigmaX + muX;

    // Calculate p-value
    double pval;
    //boost::math::chi_squared_distribution<> dist(l, d);

    boost::math::non_central_chi_squared_distribution<double> dist(l, d);


    //std::cout << "Q_Norm1 " << Q_Norm1 << std::endl;
  if(Q_Norm1 > 0){  
    if(!log_p){
        pval = boost::math::cdf(complement(dist, Q_Norm1)); //lower.tail=FALSE
    }else{
        pval = boost::math::log1p( boost::math::cdf(complement(dist, Q_Norm1)));
    }
  }else{
    pval = 1.0;  
  }
     p_value(i) = pval;
  }

  return p_value;
}




// Main function to calculate p-values
// lambda - eigenvalues
// Q - test statistics
// df1 - optional degrees of freedom vector
// Returns a list containing p-values, Liu p-values, convergence indicators, log-transformed p-values, and zero p-value messages
// [[Rcpp::export]]
Rcpp::List Get_PValue_Lambda_Rcpp(arma::vec lambda, arma::vec Q, arma::vec df1) {
  int n1 = Q.n_elem;

  arma::vec pValue(n1, arma::fill::zeros);
  arma::vec pValueLiu(n1, arma::fill::zeros);
  arma::ivec isConverge(n1, arma::fill::zeros);
  bool logp = false;
  pValueLiu = Get_Liu_PVal_MOD_Lambda_Rcpp(Q, lambda, df1, logp);

  Rcpp::List out;
  double q, sigma, acc;
  int lim = 10000;
  sigma=0;
  acc=std::pow(10.0,-6.0);
  arma::ivec h(lambda.n_elem);
  h.ones();
  arma::vec delta(lambda.n_elem);
  delta.zeros();
  for (int i = 0; i < n1; ++i) {
    q = Q[i];
    if (df1.n_elem == 0) {
	out=SKAT_davies_Rcpp(q,lambda,h,delta, sigma, lim, acc);
    } else {
        h = arma::conv_to < arma::ivec >::from(df1);
	out=SKAT_davies_Rcpp(q,lambda,h,delta, sigma, lim, acc);
    }
    pValue(i) = out["Qq"];
    // pValueLiu(i) = SKAT_liu(Q[i], lambda);

    isConverge(i) = 1;
    int ifault0 = out["ifault"];
    if (lambda.n_elem == 1) {
      pValue(i) = pValueLiu(i);
    } else if (ifault0 != 0) {
      isConverge(i) = 0;
    }

    if (pValue(i) > 1 || pValue(i) <= 0) {
      isConverge(i) = 0;
      pValue(i) = pValueLiu(i);
    }
  }

  std::string pValueMsg;
  arma::vec pValueLog;

  if (pValue(0) == 0) {
    Rcpp::List param = Get_Liu_Params_Mod_Lambda_Rcpp(lambda, df1);
    pValueMsg = Get_Liu_PVal_MOD_Lambda_Zero_Rcpp(Q(0), param["muQ"], param["muX"], param["sigmaQ"], param["sigmaX"], param["l"], param["d"]);
    arma::vec Q0(1);
    Q0(0) = Q(0);
    df1.clear();
    bool logp=false;
    pValueLog = Get_Liu_PVal_MOD_Lambda_Rcpp(Q0, lambda, df1, logp);
  }

  return Rcpp::List::create(Rcpp::Named("p.value") = pValue,
                      Rcpp::Named("p.val.liu") = pValueLiu,
                      Rcpp::Named("is_converge") = isConverge,
                      Rcpp::Named("p.val.log") = pValueLog,
                      Rcpp::Named("pval.zero.msg") = pValueMsg);
}


arma::vec Get_Lambda_Approx_Rcpp(const arma::mat& K, int maxK = 100) {

  int p1 = K.n_cols;
  arma::vec lambda(p1, arma::fill::zeros);
  
  double lambda_sum = arma::accu(K.diag());
  double lambda2_sum = arma::accu(K % K);
    
  int k1 = std::min(p1 / 10, maxK);
  
  arma::vec lambda1(k1);
  
  if (k1 > maxK) {
    k1 = maxK;
  }
  
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, K);
  
  arma::uvec IDX1 = find(eigval >= 0);
  arma::uvec IDX2 = find(eigval > arma::mean(eigval(IDX1)) / 100000);
  
  if (IDX2.is_empty()) {
    Rcpp::stop("No Eigenvalue is bigger than 0!!");
  } 
  
  if (IDX2.n_elem < k1) {
    k1 = IDX2.n_elem;
  }  
      
  lambda1 = eigval.subvec(IDX2(0), IDX2(k1 - 1));
  
  lambda_sum -= arma::accu(lambda1);
  lambda2_sum -= arma::accu(lambda1 % lambda1);
  
  int df1 = std::round(std::pow(lambda_sum, 2) / lambda2_sum);
  
  arma::vec df1_vec(lambda1.n_elem + 1, arma::fill::ones);
  arma::vec lambda_last = arma::vec(1).fill(lambda_sum / df1);
  
  arma::vec lambda_out(lambda1.n_elem + 1);
  lambda_out.subvec(0, lambda1.n_elem - 1) = lambda1;
  lambda_out(lambda1.n_elem) = lambda_last[0];
  
  df1_vec(lambda1.n_elem) = df1;
  
  return lambda_out;
}


Rcpp::List Get_PValue_Rcpp(const arma::mat& K, const arma::vec& Q, bool isFast = false) {
  
  arma::vec lambda;
  arma::vec df1;
  
  int p1 = K.n_cols;
  int maxK = 100; 
  if (!isFast || p1 < 2000) {
    lambda = Get_Lambda_Rcpp(K, isFast, maxK);
    return Get_PValue_Lambda_Rcpp(lambda, Q, df1);
  } else {
    arma::vec lambda_approx = Get_Lambda_Approx_Rcpp(K);
    lambda = lambda_approx.head(lambda_approx.n_elem - 1);
    df1 = lambda_approx.tail(1);
    return Get_PValue_Lambda_Rcpp(lambda, Q, df1);
  }
}



Rcpp::List Get_Davies_PVal_Rcpp(const arma::mat& Q_mat, const arma::mat& W_mat, const arma::mat& Q_resampling_mat, bool isFast = false) {
  
  
  arma::mat K = W_mat / 2.0;
  
  arma::vec Q_all;
  if (!Q_resampling_mat.is_empty()) {
    arma::mat Q_mat2  = arma::join_cols(Q_mat, Q_resampling_mat);
    Q_all = arma::vectorise(Q_mat2);
  } else {
    Q_all = arma::vectorise(Q_mat);
  }
  //std::cout << "Get_Davies_PVal_Rcpp 1" << std::endl; 
  //Q_all.print("Q_all");
  Rcpp::List re = Get_PValue_Rcpp(K, Q_all, isFast);
  //std::cout << "Get_Davies_PVal_Rcpp 2" << std::endl; 
  Rcpp::List param;
 
  double liu_pval;
  arma::vec pvalliu = re["p.val.liu"];
  liu_pval = pvalliu[0];
  arma::ivec isconverge = re["is_converge"];
  int is_converge = isconverge[0];


  param["liu_pval"] = liu_pval ;

  //std::cout << "liu_pval " << liu_pval << std::endl;

  param["Is_Converged"] = is_converge;
  
  arma::vec p_value_resampling;
  arma::vec liu_pval_resampling;
  arma::ivec is_converged_resampling;
  
  if (!Q_resampling_mat.is_empty()) {
    
    arma::vec p_value_resampling0 = re["p.value"];
    arma::vec liu_pval_resampling0 = pvalliu;
    arma::ivec is_converged_resampling0 = re["is_converge"];
    p_value_resampling = p_value_resampling0.subvec(1, p_value_resampling0.n_elem - 1);
    liu_pval_resampling = liu_pval_resampling0.subvec(1, liu_pval_resampling0.n_elem - 1);
    is_converged_resampling = is_converged_resampling0.subvec(1, is_converged_resampling0.n_elem - 1);
  }  
  arma::vec pvaluevec = re["p.value"];
  re["p_value"] = pvaluevec[0];
  re["param"] = param;
  re["p.value.resampling"] = p_value_resampling;
  re["liu_pval.resampling"] = liu_pval_resampling;
  re["Is_Converged.resampling"] = is_converged_resampling;
  
  return re;
}


Rcpp::List Met_SKAT_Get_Pvalue_Rcpp(const arma::vec& Score, const arma::mat& Phi, arma::vec& r_corr, const std::string& method, const arma::mat& Score_Resampling) {

  arma::uword p_m = Phi.n_rows;

  // If Phi == 0
  if (arma::accu(arma::abs(Phi)) == 0) {
    Rcpp::warning("No polymorphic SNPs!");
    return Rcpp::List::create(
      Rcpp::Named("p_value") = 1,
      Rcpp::Named("p.value.resampling") = R_NilValue,
      Rcpp::Named("pval.zero.msg") = R_NilValue
    );
  }

  // Transpose Score_Resampling if it is provided
  arma::mat Score_Resampling_t;
  if (!Score_Resampling.is_empty()) {
    Score_Resampling_t = Score_Resampling.t();
  }
  //std::cout << "Score_Resampling_t1" << std::endl;
  if (Phi.n_cols <= 1 && Phi.n_rows <= 1) {
    r_corr.set_size(1);
    r_corr[0] = 0;
  } else {
    if (Phi.n_cols <= 10) {
      if (arma::rank(Phi) <= 1) {
    	r_corr.set_size(1);
    	r_corr[0] = 0;
      }
    }
  }

  if (r_corr.n_elem > 1) {
    Rcpp::List re = SKAT_META_Optimal_Rcpp(Score, Phi, r_corr, method, Score_Resampling);
    return re;
  }
  //std::cout << "Score_Resampling_t2" << std::endl;

  Rcpp::List QL;
  arma::vec Q_1;
  arma::mat Q_res_1;
  arma::mat Phi_in = Phi;
  if (r_corr[0] == 0) {
    Q_1.set_size(1);
    Q_1(0) = arma::accu(arma::pow(Score, 2.0)) / 2.0;
    if (!Score_Resampling.is_empty()) {
      arma::mat Q_res = arma::sum(arma::square(Score_Resampling_t), 1) / 2.0;
            Q_res_1 = Q_res;
    }
  //std::cout << "Score_Resampling_t2 1a SKAT" << std::endl;
  } else if (r_corr[0] == 1) {
    QL = SKAT_META_Optimal_Get_Q_Rcpp(Score, r_corr);
    arma::vec Q = QL["Qr"];
    Q_1 = Q;
    
  //std::cout << "Score_Resampling_t2 1a" << std::endl;

    if (!Score_Resampling.is_empty()) {

      QL = SKAT_META_Optimal_Get_Q_Res_Rcpp(Score_Resampling_t, r_corr);
      arma::mat Q_res = QL["Qr"];
            Q_res_1 = Q_res;
    }

    arma::mat a(1,1);
    a(0,0) = arma::accu(Phi);
   //std::cout << "Score_Resampling_t2 1a1" << std::endl;
    Rcpp::List re = Get_Liu_PVal_Rcpp(Q_1, a, Q_res_1);
   //std::cout << "Score_Resampling_t2 1b" << std::endl;
    return re;
  } else {
   //std::cout << "Score_Resampling_t2 1c" << std::endl;
    QL = SKAT_META_Optimal_Get_Q_Rcpp(Score, r_corr);
    arma::vec Q = QL["Qr"];
    Q_1 = Q;

    if (!Score_Resampling.is_empty()) {
      QL = SKAT_META_Optimal_Get_Q_Res_Rcpp(Score_Resampling_t, r_corr);
      arma::mat Q_res = QL["Qr"];
      Q_res_1 = Q_res;
    }

    arma::mat R_M = arma::diagmat(arma::ones(p_m) * (1 - r_corr)) + arma::repmat(r_corr, p_m, p_m);
    arma::mat L = arma::chol(R_M, "lower");
    Phi_in = L * (Phi * L.t());
  }
   //std::cout << "Score_Resampling_t3" << std::endl;
   //Q_1.print("Q");
   //Phi_in.print("Phi");


  Rcpp::List re = Get_Davies_PVal_Rcpp(Q_1, Phi_in, Q_res_1);
   //std::cout << "Score_Resampling_t4" << std::endl; 
  if (r_corr.n_elem == 1) {
    re["Q"] = Q_1;
  }



  return re;
}


// RcppArmadillo implementation of Get_Liu_PVal function
// [[Rcpp::export]]
Rcpp::List Get_Liu_PVal_Rcpp(const arma::mat& Q, const arma::mat& W, const arma::mat& Q_resampling) {
  arma::mat Q_all_mat = arma::join_cols(Q, Q_resampling);
  arma::vec Q_all = arma::vectorise(Q_all_mat);
  arma::mat A1 = W / 2;
  arma::mat A2 = A1 * A1;

  arma::vec c1(4);
  c1(0) = arma::accu(arma::diagmat(A1));
  c1(1) = arma::accu(arma::diagmat(A2));
  c1(2) = arma::accu(A1 % arma::trans(A2));
  c1(3) = arma::accu(A2 % arma::trans(A2));
  //std::cout << "Get_Liu_PVal_Rcpp 0" << std::endl;

  Rcpp::List param = Get_Liu_Params_Rcpp(c1);
  double muQ = param["muQ"];
  double muX = param["muX"];
  double sigmaQ = param["sigmaQ"];
  double sigmaX = param["sigmaX"];

  arma::vec Q_Norm = (Q_all - muQ)/sigmaQ;
  arma::vec Q_Norm1 = Q_Norm * sigmaX + muX;
 
  double l = param["l"];
  double d = param["d"];
  //std::cout << "Get_Liu_PVal_Rcpp 1" << std::endl;
  boost::math::non_central_chi_squared_distribution<double> dist(l, d);

   arma::vec p_value(Q_Norm1.n_elem);
  for (int i = 0; i < Q_Norm1.n_elem; ++i) {
  double Q_Norm1_in = Q_Norm1[i];
  if(Q_Norm1_in > 0){
  //std::cout << "Q_Norm1_in " << Q_Norm1_in << std::endl;
      p_value[i] = 1 - boost::math::cdf(dist, Q_Norm1_in);
  }else{
     p_value[i] = 1.0; 	

  }
}
  std::cout << "Get_Liu_PVal_Rcpp 2" << std::endl;


  arma::vec p_value_resampling;
  if (Q_resampling.n_elem > 0) {
    p_value_resampling = p_value.subvec(1, p_value.n_elem - 1);
  }

  return Rcpp::List::create(
  			Rcpp::Named("p_value") = p_value(0),
  			Rcpp::Named("p.value") = p_value(0),
                      Rcpp::Named("param") = param,
                      Rcpp::Named("p.value.resampling") = p_value_resampling);
}

// RcppArmadillo implementation of Get_Liu_Params function
// [[Rcpp::export]]
Rcpp::List Get_Liu_Params_Rcpp(const arma::vec& c1) {
  double muQ = c1(0);
  double sigmaQ = std::sqrt(2 * c1(1));
  double s1 = c1(2) / std::pow(c1(1), 1.5);
  double s2 = c1(3) / std::pow(c1(1), 2.0);

  double beta1 = std::sqrt(8) * s1;
  double beta2 = 12 * s2;
  int type1 = 0;

  double a, d, l, muX, sigmaX;
  if (std::pow(s1, 2.0) > s2) {
    a = 1 / (s1 - std::sqrt(std::pow(s1, 2.0) - s2));
    d = s1 * std::pow(a, 3.0) - std::pow(a, 2.0);
    l = std::pow(a, 2.0) - 2 * d;
  } else {
    type1 = 1;
    a = 1 / s1;
    d = 0;
    l = 1 / std::pow(s1, 2.0);
  }

  muX = l + d;
  sigmaX = std::sqrt(2) * a;

  return Rcpp::List::create(Rcpp::Named("l") = l,
                      Rcpp::Named("d") = d,
                      Rcpp::Named("muQ") = muQ,
                      Rcpp::Named("muX") = muX,
                      Rcpp::Named("sigmaQ") = sigmaQ,
                      Rcpp::Named("sigmaX") = sigmaX);
}

