#ifndef SKATO_HPP
#define SKATO_HPP


#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include "qfc.hpp"



Rcpp::List SKAT_META_Optimal_Get_Q_Rcpp(const arma::vec& Score, const arma::vec& r_all);
Rcpp::List SKAT_META_Optimal_Get_Q_Res_Rcpp(const arma::mat& Score_res, const arma::vec& r_all);
arma::vec Get_Lambda_Rcpp(const arma::mat& K, bool isFast, int maxK);
Rcpp::List SKAT_META_Optimal_Param_Rcpp(const arma::mat& Phi, const arma::vec& r_all);
Rcpp::List SKAT_Optimal_Each_Q_Rcpp(const arma::mat& Q_all, const arma::vec& r_all, const Rcpp::List& lambda_all, std::string method);
Rcpp::List SKAT_davies_Rcpp(double q, arma::vec& lambda, arma::ivec& h, arma::vec& delta, double sigma, int lim, double acc);




arma::vec SKAT_Optimal_Integrate_Func_Davies_Rcpp(const arma::vec& x, const arma::vec& pmin_q, arma::vec& tau, double MuQ, arma::vec& lambda, double VarQ, double VarRemain, const arma::vec& r_all);

double SKAT_Optimal_Integrate_Func_single_Davies_Rcpp(double x, arma::vec& pmin_q, arma::vec& tau, double MuQ, arma::vec& lambda, double VarQ, double VarRemain, arma::vec& r_all);


double SKAT_Optimal_PValue_Davies_Rcpp(arma::vec& pmin_q, arma::vec& tau, double MuQ, arma::vec& lambda, double VarQ, double VarRemain, arma::vec& r_all, double pmin);


Rcpp::List SKAT_META_Optimal_Get_Pvalue_Rcpp(const arma::mat& Q_all, const arma::mat& Phi, arma::vec& r_all, const std::string& method, bool isFast);

Rcpp::List SKAT_Optimal_Get_Pvalue_Rcpp(const arma::mat& Q_all, const arma::mat& Z1, const arma::vec& r_all, const std::string& method);

Rcpp::List SKAT_META_Optimal_Rcpp(const arma::mat& Score, const arma::mat& Phi, arma::vec & r_all, const std::string& method, const arma::mat& Score_Resampling);


Rcpp::List Met_SKAT_Get_Pvalue_Rcpp(const arma::vec& Score, const arma::mat& Phi, arma::vec& r_corr, const std::string& method, const arma::mat& Score_Resampling);


Rcpp::List Get_PValue_Lambda_Rcpp(arma::vec lambda, arma::vec Q, arma::vec df1);
Rcpp::List Get_Liu_Params_Mod_Lambda_Rcpp(const arma::vec& lambda, arma::vec& df1);



std::string Get_Liu_PVal_MOD_Lambda_Zero_Rcpp(double Q, double muQ, double muX, double sigmaQ, double sigmaX, double l, double d);

//arma::vec SKAT_Optimal_Integrate_Func_Liu_Rcpp(const arma::vec& x, const arma::mat& pmin_q, const arma::mat& tau, const arma::vec& MuQ, double VarQ, const arma::vec& Df, const arma::vec& r_all);

//double SKAT_Optimal_PValue_Liu_Rcpp(const arma::mat& pmin_q, const arma::mat& tau, const arma::vec& MuQ, double VarQ, const arma::vec& Df, const arma::vec& r_all, double pmin = R_NegInf);

Rcpp::List Get_Liu_Params_Mod_Rcpp(const arma::vec& c1);

Rcpp::List Get_Liu_PVal_Rcpp(const arma::mat& Q, const arma::mat& W, const arma::mat& Q_resampling);
Rcpp::List Get_Liu_Params_Rcpp(const arma::vec& c1);


arma::vec Get_Liu_PVal_MOD_Lambda_Rcpp(const arma::vec& Q_all, const arma::vec& lambda, arma::vec& df1, bool log_p);

#endif
