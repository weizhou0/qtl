#define ARMA_USE_SUPERLU 1

// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



#include "SAIGE_test.hpp"
#include "SPA.hpp"
#include "ER_binary_func.hpp"
#include "UTIL.hpp"
#include "approxfun.hpp"
#include "getMem.hpp"
#include "getMem.hpp"
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/date_time.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

namespace SAIGE {

SAIGEClass::SAIGEClass(
	arma::mat & t_XVX,
	arma::mat t_XXVX_inv,
	arma::mat & t_XV,
	arma::mat & t_XVX_inv_XV,
	arma::mat & t_Sigma_iXXSigma_iX,
	arma::mat & t_X,
	arma::vec &  t_S_a,
	arma::vec & t_res,
	arma::vec & t_mu2,
	arma::vec & t_mu,
	arma::vec & t_varRatio_sparse,
	arma::vec & t_varRatio_null,
	arma::vec & t_varRatio_null_noXadj,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
	double t_SPA_Cutoff,
	arma::vec & t_tauvec,
	arma::vec & t_varWeightsvec,
	std::string t_traitType,
	arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	bool t_isnoadjCov,
	double t_pval_cutoff_for_fastTest,
	bool t_isCondition,
        std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
        double t_pCutoffforFirth,
        arma::vec & t_offset,
	arma::vec & t_resout, 
	arma::sp_mat & t_SigmaMat_sp,
	float t_tauVal_sp,
	 arma::sp_mat & t_Ilongmat,
        arma::vec & t_I_longl_vec,
        arma::sp_mat & t_Tlongmat,
        arma::vec & t_T_longl_vec,
	bool t_is_EmpSPA,
        arma::mat & t_cumul){


    m_XVX = t_XVX;
    m_XV = t_XV;
    m_XXVX_inv = t_XXVX_inv;
    m_XVX_inv_XV = t_XVX_inv_XV;
    m_Sigma_iXXSigma_iX = t_Sigma_iXXSigma_iX;
    m_isVarPsadj = false;
    if(m_Sigma_iXXSigma_iX.n_cols == 1 && m_Sigma_iXXSigma_iX.n_rows == 1){
	m_isVarPsadj = false;
    }else{
	m_isVarPsadj = true;
    }
    m_X = t_X;
    m_S_a = t_S_a;
    m_res = t_res;
    m_resout = t_resout;
    m_mu2 = t_mu2 % t_varWeightsvec;
    m_mu = t_mu;
    m_varRatio_sparse = t_varRatio_sparse;
    m_varRatio_null = t_varRatio_null;
    m_varRatio_null_noXadj = t_varRatio_null_noXadj;
    m_cateVarRatioMinMACVecExclude = t_cateVarRatioMinMACVecExclude;
    m_cateVarRatioMaxMACVecInclude = t_cateVarRatioMaxMACVecInclude;
    m_tauvec = t_tauvec;
    m_varWeightsvec = t_varWeightsvec;
    m_traitType = t_traitType;
    m_y = t_y;

    m_case_indices = arma::find(m_y == 1);
    m_ctrl_indices = arma::find(m_y == 0);

    m_n = t_y.size();
    m_p = t_XV.n_rows;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_impute_method =  t_impute_method;
    m_isCondition = t_isCondition;
    m_condition_genoIndex = t_condition_genoIndex;
    if(m_isCondition){	
	        m_numMarker_cond = t_condition_genoIndex.size();      
    }else{
		m_numMarker_cond = 0;
    }	    
    //m_p = t_X.nrow();

    if(m_traitType == "binary"){
      //if(m_isOutputAFinCaseCtrl){
        m_case_indices = arma::find(m_y == 1);
        m_ctrl_indices = arma::find(m_y == 0);
	m_n_case = m_case_indices.n_elem;
	m_n_ctrl = m_ctrl_indices.n_elem;
	m_is_Firth_beta = t_is_Firth_beta;
	m_pCutoffforFirth = t_pCutoffforFirth;
	m_offset = t_offset;
      //}
    }
    //m_dimNum = t_dimNum;
    m_flagSparseGRM = t_flagSparseGRM;
    //m_isFastTest = t_isFastTest;
    m_isnoadjCov = t_isnoadjCov;
    m_pval_cutoff_for_fastTest = t_pval_cutoff_for_fastTest;
    //if(m_dimNum != 0){
    //	m_locationMat = t_locationMat;
    //	m_valueVec = t_valueVec;
    //}
    //
    //if(t_SigmaMat_sp.n_rows > 2){
	//std::cout << "here m_SigmaMat_sp" << std::endl;
	m_SigmaMat_sp = t_SigmaMat_sp;
    //}	   
    m_tauVal_sp = t_tauVal_sp;
   g_I_longl_mat = t_Ilongmat;
   arma::uvec t_I_longl_vec_new = arma::conv_to< arma::uvec >::from(t_I_longl_vec);
   g_I_longl_vec = t_I_longl_vec_new;
   g_T_longl_mat = t_Tlongmat;
   g_T_longl_vec = t_T_longl_vec;

   m_res_sample = (g_I_longl_mat.t()) * m_res;
   m_mu2_sample = (g_I_longl_mat.t()) * m_mu2;


   m_is_EmpSPA = false;
   if(t_is_EmpSPA){
      m_K_0_emp.setApproxFun(t_cumul.col(0), t_cumul.col(1));
      m_K_1_emp.setApproxFun(t_cumul.col(0), t_cumul.col(2));
      m_K_2_emp.setApproxFun(t_cumul.col(0), t_cumul.col(3));
      m_is_EmpSPA = t_is_EmpSPA;
      m_cumul = t_cumul;
      m_varResid = arma::var(m_res);
   }
}    


// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void SAIGEClass::set_seed(unsigned int seed){
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

void SAIGEClass::scoreTest(arma::vec & t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde,
		     arma::vec & t_P2Vec,
		     double& t_gy, 
		     bool t_is_region,
		     arma::uvec & t_indexForNonZero){
    arma::vec Sm, var2m;
    double S, var2;
    getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);


    if(t_is_region && m_traitType == "binary"){
      t_gy = dot(t_gtilde, m_y);
    }

    S = dot(t_gtilde, m_res % m_varWeightsvec);

    //std::cout << "S " << S << std::endl;

    S = S/m_tauvec[0]; 

    double varRatioVal_var2 = m_varRatioVal;
    if(!m_flagSparseGRM_cur){
      t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0]; 
      //t_P2Vec = t_gtilde % m_mu2; 
      var2m = dot(t_P2Vec, t_gtilde);
    }else{
      if(m_SigmaMat_sp.n_rows > 2){	
      //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      t_P2Vec = m_SigmaMat_sp * t_gtilde;
      var2m = dot(t_P2Vec, t_gtilde);
      if(m_isVarPsadj){
	varRatioVal_var2 = 1;     
	var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;	
      }
     }else{
	//t_P2Vec = m_sigmainvG_noV;
	t_P2Vec = getSigma_G_V(t_gtilde, 500, 1e-5);
	var2m = dot(t_P2Vec, t_gtilde);
      if(m_isVarPsadj){
	varRatioVal_var2 = 1;      
	var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;	
      }
     }
    }

    var2 = var2m(0,0);
    //std::cout << "var2 " << var2 << std::endl;
    //std::cout << "m_varRatioVal " << m_varRatioVal << std::endl;
    //double var1 = var2 * m_varRatioVal;
    double var1 = var2 * varRatioVal_var2;
    //std::cout << "var1 " << var1 << std::endl;
    double stat = S*S/var1;
    double t_pval;
    //std::cout << "S " << S << std::endl;    
    //std::cout << "var1 " << var1 << std::endl;    

    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }

    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "t_pval_str " << t_pval_str << std::endl;
}


void SAIGEClass::scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){

    arma::vec g1 = t_GVec.elem(t_indexForNonZero);

    arma::vec m_varWeightsvec_new = m_varWeightsvec.elem(t_indexForNonZero);

    arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    res1 = res1 % m_varWeightsvec_new;
    arma::vec Z = A1.t() * g1;
    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
    arma::vec S_a2;
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
    if(m_traitType == "binary" || m_traitType == "count"){
      mu21  = m_mu2.elem(t_indexForNonZero);
      g1tildemu2 = dot(square(g1_tilde), mu21);
      Bmu2 = arma::dot(square(B),  mu21);
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative" || m_traitType == "count_nb"){
      Bmu2 = dot(g1, B % m_varWeightsvec_new);
      //Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1 % m_varWeightsvec_new) - 2*Bmu2;
      //var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }
    
    var1 = var2 * m_varRatioVal;
    S1 = dot(res1, g1_tilde);
    arma::mat res1X1_temp = (res1.t()) * X1;
    arma::vec res1X1 = res1X1_temp.t();
    S_a2 = m_S_a - res1X1;
    S2 = - arma::dot(S_a2,  Z);
    S = S1 + S2;
    S = S/m_tauvec[0];

    double stat = S*S/var1;
    double t_pval;

    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1);
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }
    

    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "t_pval_str scoreTestFast " << t_pval_str << std::endl;
    //std::cout << "end of scoreTestFast" << std::endl;
}



void SAIGEClass::scoreTestFast_noadjCov(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){

    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::vec m_varWeightsvec_new = m_varWeightsvec.elem(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res_sample.elem(t_indexForNonZero);
    res1 = res1 % m_varWeightsvec_new;
    double g1mu2, mu2f, mu2fg, altFreq2;
    altFreq2 = t_altFreq * 2;
    mu2f = arma::sum(m_mu2_sample * pow(altFreq2,2));
    mu21  = m_mu2_sample.elem(t_indexForNonZero);
    g1mu2 = dot(square(g1), mu21);

    mu2fg = dot(g1, mu21 * altFreq2);
    double var2 = g1mu2 - 2*mu2fg + mu2f;

    //if(m_traitType == "binary" || m_traitType == "count"){
    if(m_traitType == "quantitative" || m_traitType == "count_nb"){
      var2 = var2 * m_tauvec[0];
    }


    double var1 = var2 * m_varRatioVal;
    double S1 = dot(res1, g1);
    double S2 = arma::sum(m_res * altFreq2);
    double S = S1 - S2;
    S = S/m_tauvec[0];


    double stat = S*S/var1;
    double t_pval;

    //if (var1 <= std::pow(std::numeric_limits<double>::min(), 2)){
    if (var1 <= std::numeric_limits<double>::min()){
        t_pval = 1;
    } else{
      boost::math::chi_squared chisq_dist(1);
      t_pval = boost::math::cdf(complement(chisq_dist, stat));
    }


    char pValueBuf[100];
    if (t_pval != 0)
        sprintf(pValueBuf, "%.6E", t_pval);
    else {
        double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    //std::cout << "t_pval_str scoreTestFast " << t_pval_str << std::endl;
    //std::cout << "end of scoreTestFast" << std::endl;
}


void SAIGEClass::getadjG(arma::vec & t_GVec, arma::vec & g){
   g = m_XV * t_GVec;
   //t_GVec.t().print("t_GVec");
   //std::cout << "sum(t_GVec) " << arma::accu(t_GVec) << std::endl;
   //      m_XV.print("m_XV");
    g = t_GVec - m_XXVX_inv * g;
    //m_XXVX_inv.print("m_XXVX_inv");
    //g.t().print("g");
}


void SAIGEClass::getadjGFast(arma::vec & t_GVec, arma::vec & g, arma::uvec & iIndex)
{

   arma::vec t_GVec0;	
   arma::uvec indexNonZeroVec0_arma;	
   if(g_I_longl_mat.n_cols != g_I_longl_mat.n_rows && t_GVec.n_elem != m_XV.n_cols){
      t_GVec0 = g_I_longl_mat * t_GVec;
      indexNonZeroVec0_arma = arma::find(t_GVec0 > 0.0);
   }else{
      t_GVec0 = t_GVec;
      indexNonZeroVec0_arma = iIndex;
   }	   


  // To increase computational efficiency when lots of GVec elements are 0
 arma::vec m_XVG(m_p, arma::fill::zeros);
  for(int i = 0; i < indexNonZeroVec0_arma.n_elem; i++){
      m_XVG += m_XV.col(indexNonZeroVec0_arma(i)) * t_GVec0(indexNonZeroVec0_arma(i));
  }
  g = t_GVec0 - m_XXVX_inv * m_XVG; 
}


void SAIGEClass::get_mu(arma::vec & t_mu){
    t_mu = m_mu;
}

void SAIGEClass::getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices){
     t_case_indices = m_case_indices;
     t_ctrl_indices = m_ctrl_indices;
  }


void SAIGEClass::setupSparseMat(int r, arma::umat & locationMatinR, arma::vec & valueVecinR) {
    m_locationMat = locationMatinR;
    m_valueVec = valueVecinR;
    m_dimNum = r;
}



arma::sp_mat SAIGEClass::gen_sp_SigmaMat() {
    //std::cout << "gen_sp_SigmaMat " << std::endl;
    arma::sp_mat resultMat(m_locationMat, m_valueVec, m_dimNum, m_dimNum);
    //std::cout << "gen_sp_SigmaMat2 " << std::endl;
    return resultMat;
}

// revised
// need to add sparse Sigma version 
// This function only uses variance ratio and does not use sparse GRM
void SAIGEClass::getMarkerPval(arma::vec & t_GVec,
			       arma::uvec & iIndex,
			       arma::uvec & iIndexComVec,
                               double& t_Beta, 
                               double& t_seBeta, 
                               double& t_pval,
			       double& t_pval_noSPA, 
                               double t_altFreq,
                               double& t_Tstat,
			       double& t_gy,
			       double& t_var1, 
			       bool & t_isSPAConverge, 
			       arma::vec & t_gtilde, 
			       bool & is_gtilde,
			       bool  is_region,
                               arma::vec & t_P2Vec,
			       bool t_isCondition,
			       double& t_Beta_c,
                           	double& t_seBeta_c,
                           	double& t_pval_c,
                          	double& t_pval_noSPA_c,
                           	double& t_Tstat_c,
                           	double& t_varT_c,
			   	arma::rowvec & t_G1tilde_P_G2tilde, 
				bool & t_isFirth,
				bool & t_isFirthConverge, 
				bool t_isER, 
				bool t_isnoadjCov,
                                bool t_isSparseGRM)
{



  t_isFirth = false;
  std::string t_pval_str;
  double t_var2, t_SPApval;
   double altFreq0;
   arma::vec t_GVec0;
   arma::uvec indexNonZeroVec0_arma, indexZeroVec0_arma;
   if(g_I_longl_mat.n_cols != g_I_longl_mat.n_rows){
      t_GVec0 = g_I_longl_mat * t_GVec;
      altFreq0 = arma::mean(t_GVec0) /2;
      indexNonZeroVec0_arma = arma::find(t_GVec0 > 0.0);
      indexZeroVec0_arma = arma::find(t_GVec0 == 0.0);
   }else{
      t_GVec0 = t_GVec;
      altFreq0 = t_altFreq;
      indexNonZeroVec0_arma = iIndex; 	
      indexZeroVec0_arma = iIndexComVec;
   }

 //arma::vec timeoutput3 = getTime();
 //
 //
 if(t_isSparseGRM){
 	t_isnoadjCov = false;
 }
  if(!t_isnoadjCov){
	//std::cout << "scoreTest " << std::endl;  
	unsigned int nonzero = indexNonZeroVec0_arma.n_elem;
	unsigned int nGvec = t_GVec0.n_elem;
 	//std::cout << "nonzero " << nonzero << std::endl;
	//std::cout << "nGvec " << nGvec << std::endl;
	//
  	//scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, iIndex);
	if(!t_isSparseGRM){
  	  is_gtilde = false;
	  //std::cout << "t_isSparseGRM " << t_isSparseGRM << std::endl;	
	  scoreTestFast(t_GVec0, indexNonZeroVec0_arma, t_Beta, t_seBeta, t_pval_str, altFreq0, t_Tstat, t_var1, t_var2);
	}else{
  	  is_gtilde = true;
	  scoreTest(t_GVec0, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, indexNonZeroVec0_arma);	
	}	
  }else{
  	is_gtilde = false;
	unsigned int nonzero = iIndex.n_elem;
	unsigned int nGvec = t_GVec.n_elem;
	scoreTestFast_noadjCov(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_str, altFreq0,t_Tstat, t_var1, t_var2);
  }





  double StdStat = std::abs(t_Tstat) / sqrt(t_var1);
  //std::cout << "before SPA t_Tstat " << t_Tstat << std::endl;
  //std::cout << "before SPA t_var1 " << t_var1 << std::endl;
  //std::cout << "t_altFreq " << t_altFreq << std::endl;
  t_isSPAConverge = false;

  double pval_noadj;
  try {
        pval_noadj = std::stod(t_pval_str);
  } catch (const std::invalid_argument&) {
        pval_noadj = 0;
        std::cerr << "Argument is invalid\n";
        //throw;
  } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        //throw;
        pval_noadj = 0;
  }

  
 //arma::vec timeoutput3_a = getTime();
  double q, qinv, m1, NAmu, NAsigma, tol1, p_iIndexComVecSize;

  //arma::uvec iIndexComVec = arma::find(t_GVec == 0);
  //arma::uvec iIndexVec = arma::find(t_GVec != 0);
 
  
  unsigned int iIndexComVecSize = iIndexComVec.n_elem;
  unsigned int iIndexSize = iIndex.n_elem; 
  arma::vec gNB(iIndexSize, arma::fill::none);
  arma::vec gNA(iIndexComVecSize, arma::fill::none);
  arma::vec muNB(iIndexSize, arma::fill::none);
  arma::vec muNA(iIndexComVecSize, arma::fill::none);

  /*
    std::cout << "iIndexComVecSize " << iIndexComVecSize << std::endl;
    std::cout << "iIndexSize " << iIndexSize << std::endl;
    std::cout << "gNA.n_elem 1 " << gNA.n_elem << std::endl;
    std::cout << "gNB.n_elem 1 " << gNB.n_elem << std::endl;
    std::cout << "muNA.n_elem 1 " << muNA.n_elem << std::endl;
    std::cout << "muNB.n_elem 1 " << muNB.n_elem << std::endl;
  */


  double gmuNB;

if(!t_isER){


  //if(StdStat > m_SPA_Cutoff && m_traitType != "quantitative"){
  if(StdStat > m_SPA_Cutoff && m_traitType == "binary"){

       if(!is_gtilde){
          t_gtilde.resize(m_n);
          getadjGFast(t_GVec, t_gtilde, iIndex);
	  is_gtilde = true;
       }

    if(!m_is_EmpSPA){

	//int t_gtilden = t_gtilde.n_elem;
        p_iIndexComVecSize = double(iIndexComVecSize)/m_n;
   	m1 = dot(m_mu, t_gtilde);

	if(p_iIndexComVecSize >= 0.5){
		unsigned int j1 = 0;
		unsigned int j2 = 0;
/*
		for(unsigned int j = 0; j < m_n ; j++){	
			//std::cout << "j " << j << std::endl;
			if(t_GVec(j) != 0){
			//std::cout << "j1 " << j1 << std::endl;
				gNB(j1) = t_gtilde(j);
				muNB(j1) = m_mu(j);
				j1 = j1 + 1;	
			}else{
			//std::cout << "j2 " << j2 << std::endl;
				gNA(j2) = t_gtilde(j);
				muNA(j2) = m_mu(j);
				j2 = j2 + 1;
          //process_mem_usage(mem1, mem2);
//   std::cout << "VM 4 a 1.3c: " << mem1/1000 << "; RSS 4 a 1.3: " << mem2/1000 << std::endl;
			}	
		}
		*/
//	std::cout << "gNB.n_elem " <<  gNB.n_elem << std::endl;	
//	std::cout << "gNA.n_elem " <<  gNA.n_elem << std::endl;	
	gNB = t_gtilde(iIndex);
	gNA = t_gtilde(iIndexComVec);
   	muNB = m_mu(iIndex);
   	muNA = m_mu(iIndexComVec);

	/*
	    std::cout << "gNA.n_elem 2 " << gNA.n_elem << std::endl;
        std::cout << "gNB.n_elem 2 " << gNB.n_elem << std::endl;
        std::cout << "muNA.n_elem 2 " << muNA.n_elem << std::endl;
        std::cout << "muNB.n_elem 2 " << muNB.n_elem << std::endl;
*/


  	gmuNB = dot(gNB,muNB);	 
   	NAmu= m1-gmuNB;
	}
	/*else{
		gNA.clear();
		gNB.clear();
		muNA.clear();
		gNB.clear();

	}*/

   	if(m_traitType == "binary"){
                q = t_Tstat/sqrt(t_var1/t_var2) + m1;

                if((q-m1) > 0){
                        qinv = -1 * std::abs(q-m1) + m1;
                }else if ((q-m1) == 0){
                        qinv =  m1;
                }else{
                        qinv = std::abs(q-m1) + m1;
                }
		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % (1-muNB) % arma::pow(gNB,2));
		}
        }else if(m_traitType == "survival" || m_traitType == "count"){
                q = t_Tstat/sqrt(t_var1/t_var2);
                qinv = -q;
		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % arma::pow(gNB,2));
		}
        }
    	bool logp=false;
	double tol0 = std::numeric_limits<double>::epsilon();
	tol1 = std::pow(tol0, 0.25);
	if(p_iIndexComVecSize >= 0.5){
		//std::cout << "SPA_fast" << std::endl;
        	SPA_fast(m_mu, t_gtilde, q, qinv, pval_noadj, false, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge);
	}else{
		//std::cout << "SPA" << std::endl;
		SPA(m_mu, t_gtilde, q, qinv, pval_noadj, tol1, logp, m_traitType, t_SPApval, t_isSPAConverge);	
	}

	//std::cout << "t_SPApval " << t_SPApval << std::endl;
   }else{ //if(!m_is_EmpSPA){
	double g_altFreq_new = arma::mean(t_GVec)/2;
        int g_N = t_GVec.n_elem;
        //double t_Tstat_G = dot(t_GVec, m_res % m_varWeightsvec); 
        //double t_Tstat_G = dot(t_gtilde, m_res % m_varWeightsvec); 
        //t_Tstat_G = t_Tstat_G/m_tauvec[0];
	//t_SPApval = EmpSPA_getMarkerPval(t_gtilde, t_Tstat_G, g_altFreq_new, g_N);
	t_SPApval = EmpSPA_getMarkerPval(t_GVec, StdStat, g_altFreq_new, g_N);
	//std::cout << "Empirical SPA p value " << t_SPApval << std::endl;   
	//std::cout << "g_altFreq_new " << g_altFreq_new << std::endl;   
   }	   



    	boost::math::normal ns;
	t_pval = t_SPApval;
    	double t_qval;
        try {
           t_qval = boost::math::quantile(ns, t_pval/2);
           t_qval = fabs(t_qval);
           t_seBeta = fabs(t_Beta)/t_qval;
        }catch (const std::overflow_error&) {
          t_qval = std::numeric_limits<double>::infinity();
          t_seBeta = 0;
        } 
  }
   t_pval_noSPA = pval_noadj; 
   if(m_traitType!="quantitative"){
   //if(m_traitType=="binary"){
        if(t_isSPAConverge){
                t_pval = t_SPApval;
                //t_pval = pval_noadj;
        }else{
                t_pval = pval_noadj;
        }
   }else{
        t_pval = t_pval_noSPA;
   }


}else{ //if(!t_isER){

    t_pval_noSPA = pval_noadj; 
    arma::mat Z_er(t_GVec.n_elem, 1);
    Z_er.col(0) = t_GVec;
    arma::vec res_er = m_res;
    arma::vec pi1_er = m_mu;
    arma::vec resout_er = m_resout;
    t_pval =  SKATExactBin_Work(Z_er, res_er, pi1_er, m_n_case, iIndex, iIndexComVec, resout_er, 2e+6, 1e+4, 1e-6, 1);
	
    boost::math::normal ns;
    double t_qval;
    try{
      t_qval = boost::math::quantile(ns, t_pval/2);
      t_qval = fabs(t_qval);
      t_seBeta = fabs(t_Beta)/t_qval;
    }catch (const std::overflow_error&) {
      t_qval = std::numeric_limits<double>::infinity();
      t_seBeta = 0;
    }
}




   if(m_traitType == "binary" & m_is_Firth_beta & t_pval <= m_pCutoffforFirth){
	t_isFirth = true;

	if(!is_gtilde){
                getadjGFast(t_GVec, t_gtilde, iIndex);
                is_gtilde = true;
        }
	arma::mat x(t_GVec.n_elem, 2, arma::fill::ones);	
	x.col(1) = t_gtilde;
	//x.col(1) = t_GVec;
	arma::vec init(2, arma::fill::zeros);
	//std::cout << "t_Beta " << t_Beta << std::endl;
	//std::cout << "t_seBeta " << t_seBeta << std::endl;
	fast_logistf_fit_simple(x, m_y, m_offset, true, init, 50, 15, 15, 1e-5, 1e-5, 1e-5, t_Beta ,t_seBeta, t_isFirthConverge);	
	//std::cout << "t_Beta after " << t_Beta << std::endl;
	//std::cout << "t_seBeta after " << t_seBeta << std::endl;
   }
   
//arma::vec timeoutput4 = getTime();
//printTime(timeoutput3, timeoutput3_a, "Test Marker  ScoreTest");
//printTime(timeoutput3, timeoutput4, "Test Marker 3 to 4");
//printTime(timeoutput3_a, timeoutput4, "Test Marker SPA");







   //condition
   if(t_isCondition){
	if(!is_gtilde){
        	getadjGFast(t_GVec, t_gtilde, iIndex);
        	is_gtilde = true;
        }
        t_G1tilde_P_G2tilde = sqrt(m_varRatioVal) * t_gtilde.t() * m_P2Mat_cond;
        arma::vec t_Tstat_ctemp =  t_G1tilde_P_G2tilde * m_VarInvMat_cond * m_Tstat_cond;
	arma::mat tempgP2 = t_gtilde.t() * m_P2Mat_cond;

    	t_Tstat_c = t_Tstat - t_Tstat_ctemp(0);
    	arma::vec t_varT_ctemp = t_G1tilde_P_G2tilde * m_VarInvMat_cond * (t_G1tilde_P_G2tilde.t());
    	t_varT_c = t_var1 - t_varT_ctemp(0);

    double S_c = t_Tstat_c;

    double stat_c = S_c*S_c/t_varT_c;


     //if (t_varT_c <= std::pow(std::numeric_limits<double>::min(), 2)){
     if (t_varT_c <= std::numeric_limits<double>::min()){
        t_pval_noSPA_c = 1;
	stat_c = 0;
     }else{
        boost::math::chi_squared chisq_dist(1);
        t_pval_noSPA_c = boost::math::cdf(complement(chisq_dist, stat_c));
     }

    char pValueBuf_c[100];
    if (t_pval_noSPA_c != 0)
        sprintf(pValueBuf_c, "%.6E", t_pval_noSPA_c);
    else {
        double log10p_c = log10(2.0) - M_LOG10E*stat_c/2 - 0.5*log10(stat_c*2*M_PI);
        int exponent_c = floor(log10p_c);
        double fraction_c = pow(10.0, log10p_c - exponent_c);
        if (fraction_c >= 9.95) {
          fraction_c = 1;
           exponent_c++;
         }
        sprintf(pValueBuf_c, "%.1fE%d", fraction_c, exponent_c);
    }
    std::string buffAsStdStr_c = pValueBuf_c;
    std::string& t_pval_noSPA_str_c = buffAsStdStr_c;

    t_Beta_c = S_c/t_varT_c;
    t_seBeta_c = fabs(t_Beta_c) / sqrt(stat_c);
    t_Tstat_c = S_c;
/*
*/

  double pval_noSPA_c;  
  try {
        pval_noSPA_c = std::stod(t_pval_noSPA_str_c);
  } catch (const std::invalid_argument&) {
        pval_noSPA_c = 0;
        std::cerr << "Argument is invalid\n";
        //throw;
  } catch (const std::out_of_range&) {
        std::cerr << "Argument is out of range for a double\n";
        //throw;
        pval_noSPA_c = 0;
  }
  t_pval_noSPA_c = pval_noSPA_c; 


    //if(m_traitType != "quantitative" && stat_c > std::pow(m_SPA_Cutoff,2)){
    if(m_traitType == "binary" && stat_c > std::pow(m_SPA_Cutoff,2)){
	bool t_isSPAConverge_c;
	double q_c, qinv_c, pval_noadj_c, SPApval_c;    
	if(m_traitType == "binary"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2) + m1;

                if((q_c-m1) > 0){
                        qinv_c = -1 * std::abs(q_c-m1) + m1;
                }else if ((q_c-m1) == 0){
                        qinv_c =  m1;
                }else{
                        qinv_c = std::abs(q_c-m1) + m1;
                }
        }else if(m_traitType == "count"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2);
                qinv = -q_c;
        }

        bool logp=false;

        if(p_iIndexComVecSize >= 0.5){
		//std::cout << "SPA_fast " << std::endl;
                SPA_fast(m_mu, t_gtilde, q_c, qinv_c, pval_noadj_c, false, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, SPApval_c, t_isSPAConverge_c);
        }else{
		//std::cout << "SPA " << std::endl;
                SPA(m_mu, t_gtilde, q_c, qinv_c, pval_noadj_c, tol1, logp, m_traitType, SPApval_c, t_isSPAConverge_c);
        }

        boost::math::normal ns;
        t_pval_c = SPApval_c;
        double t_qval_c;
        try {
           t_qval_c = boost::math::quantile(ns, t_pval_c/2);
           t_qval_c = fabs(t_qval_c);
           t_seBeta_c = fabs(t_Beta_c)/t_qval_c;
        }catch (const std::overflow_error&) {
          t_qval_c = std::numeric_limits<double>::infinity();
          t_seBeta_c = 0;
        }
    }else{
    	t_pval_c = t_pval_noSPA_c;	    
    }	    
   }

    gNA.clear();
    gNB.clear();
    muNA.clear();
    gNB.clear();



    if(is_region && !is_gtilde){
	    arma::vec t_GVec0 = g_I_longl_mat * t_GVec;
	    arma::uvec iIndex0 = arma::find(t_GVec0 > 0.0);


	getadjGFast(t_GVec, t_gtilde, iIndex);
	//getadjGFast(t_GVec0, t_gtilde, iIndex0);
	is_gtilde = true; 
    }

    //if(is_region && isScoreFast){
    if(is_region && (t_isnoadjCov || (!t_isnoadjCov && !t_isSparseGRM))){

      t_gy = dot(t_gtilde, m_y);


      if(!m_flagSparseGRM_cur){
        t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0];
      }else{

	 if(m_SigmaMat_sp.n_rows > 2){


      t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      //var2m = dot(t_P2Vec , t_gtilde);
      //if(m_isVarPsadj){
      //  var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;
      //}
     }else{
        t_P2Vec = m_sigmainvG_noV;
        //var2m = dot(t_P2Vec , t_GVec);
     }



        //arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat();
        //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      }
    }
    //std::cout << "end of getPval" << std::endl;
}


bool SAIGEClass::assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
	m_varRatio = m_varRatio_sparse;
    }else{
	if(!isnoXadj){    
	    m_varRatio = m_varRatio_null;
	}else{
	    m_varRatio = m_varRatio_null_noXadj;
	}	
    }
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){    	    
		m_varRatioVal = m_varRatio(i);
		hasVarRatio = true;
	}	
    }

    if(!hasVarRatio){	
	if(MAC < m_cateVarRatioMinMACVecExclude(0)){
		m_varRatioVal = m_varRatio(0);
		hasVarRatio = true;
	}	
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
		//m_varRatioVal = m_varRatio(a-1);
		m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }

    return(hasVarRatio);    
}

void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj){ 
    arma::vec m_varRatio;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse;
    }else{
	if(isnoXadj){    
            m_varRatio = m_varRatio_null_noXadj;
	}else{
	    m_varRatio = m_varRatio_null;	
	}	
    }	
    m_varRatioVal = m_varRatio(0);
}


void SAIGEClass::assignSingleVarianceRatio_withinput(double t_varRatioVal){
        m_varRatioVal = t_varRatioVal;
}


void SAIGEClass::assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::vec & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      double t_qsum_cond,
      arma::vec & t_gsum_cond,
      arma::vec & t_p_cond
      ){
	m_P2Mat_cond = t_P2Mat_cond;
	m_VarInvMat_cond = t_VarInvMat_cond;
	m_VarMat_cond = t_VarMat_cond;
	m_Tstat_cond = t_Tstat_cond;
	m_MAF_cond = t_MAF_cond;
	m_qsum_cond = t_qsum_cond;
	m_gsum_cond = t_gsum_cond;
	m_G2_Weight_cond = t_G2_Weight_cond;
	m_p_cond = t_p_cond;
}

void SAIGEClass::assignConditionFactors_scalefactor(
	arma::vec & t_scalefactor_G2_cond	
		){
	m_scalefactor_G2_cond = t_scalefactor_G2_cond;
	arma::mat scalefactor_G2_cond_Mat = arma::diagmat(arma::sqrt(m_scalefactor_G2_cond));
	arma::mat weightMat_G2_G2 = m_G2_Weight_cond * m_G2_Weight_cond.t(); 
	arma::mat VarMat_cond_scaled = scalefactor_G2_cond_Mat * m_VarMat_cond * scalefactor_G2_cond_Mat;
	arma::mat VarMat_cond_scaled_weighted = VarMat_cond_scaled % weightMat_G2_G2;
	m_VarInvMat_cond_scaled_weighted = VarMat_cond_scaled_weighted.i();
	//m_VarInvMat_cond_region_binary = (1/scalefactor_G2_cond_Mat) * m_VarInvMat_cond	* (1/scalefactor_G2_cond_Mat);
	
}

void SAIGEClass::extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv){
	t_XV = m_XV;
	t_XXVX_inv = m_XXVX_inv;	
}



void SAIGEClass::fast_logistf_fit_simple(arma::mat & x,
                arma::vec & y,
                arma::vec & offset,
                bool firth,
        arma::vec init,
        int maxit,
        int maxstep,
        int maxhs,
        double lconv,
        double gconv,
        double xconv,
        double & beta_G,
        double & sebeta_G, 
	bool & isfirthconverge){
  isfirthconverge = false;	
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  arma::vec pi_0 = -x * beta - offset;
  pi_0 = arma::exp(pi_0) + 1;
  arma::vec pi = 1/pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  while(iter <= maxit){
        beta_old = beta;
        arma::vec wpi = pi % (1 - pi);
        arma::vec W2 = arma::sqrt(wpi);
        //arma::vec wpi_sqrt = arma::sqrt(wpi);
        //arma::vec W2 = weight % wpi_sqrt;
        arma::mat XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                XW2.col(j) = x.col(j) % W2;
        }

        arma::mat Q;
        arma::mat R;
        arma::qr_econ(Q, R, XW2);
        arma::vec h = Q % Q * oneVec;
        arma::vec U_star(2, arma::fill::zeros);
        arma::vec ypih;
        if(firth){
                ypih = (y - pi) + (h % (0.5 - pi));
        }else{
                ypih = (y - pi);
        }
        //ypih.print();
        arma::vec xcol(n, arma::fill::zeros);
        U_star = x.t() * ypih;

        arma::mat XX_XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                xcol = x.col(j);
                XX_XW2.col(j) = xcol % W2;
        }
        arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
        bool isinv = arma::inv_sympd (XX_covs, XX_Fisher);
        
	if(!isinv){
                break;
        }
        //}
        arma::vec delta = XX_covs * U_star;
        //delta.replace(arma::datum::nan, 0);

        double mx = arma::max(arma::abs(delta))/maxstep;
        if(mx > 1){
                delta = delta/mx;
        }
        evals = evals + 1;
        iter = iter + 1;
        beta = beta + delta;
        pi_0 = -x * beta - offset;
        pi_0 = arma::exp(pi_0) + 1;
        pi = 1/pi_0;
        if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
		isfirthconverge = true;
                break;
        }
  }
        arma::mat var;
        if(XX_covs.has_nan()){
                var = XX_covs;
                beta_G = arma::datum::nan;
                sebeta_G = arma::datum::nan;
        }else{
                beta_G = beta(1);
                sebeta_G = sqrt(XX_covs(1,1));
        }

	//std::cout << "beta_G " << beta_G << std::endl;
	//std::cout << "sebeta_G " << sebeta_G << std::endl;
        //return beta;
}

void SAIGEClass::set_flagSparseGRM_cur(bool t_flagSparseGRM_cur){
	m_flagSparseGRM_cur = t_flagSparseGRM_cur;
}




arma::vec SAIGEClass::getDiagOfSigma_V(){
	double tauVal0 = m_tauvec(0);

        int Nnomissing = m_mu2.n_elem;
        arma::vec diagVec(Nnomissing);
        arma::sp_vec diagVecV0;
        arma::vec diagVecG, diagVecV, diagVecG_I, diagVecG_T, diagVecG_IT,diagVecV_I, diagVecV_T, diagVecV_IT;
        unsigned int tauind = 0;

        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){

             //diagVec = tauVal0/m_mu2;
             diagVec = 1/m_mu2;
             tauind = tauind + 1;
             //diagVecG = spV.diag();
             diagVec = diagVec + m_tauVal_sp;
             tauind = tauind + 1;


        }else{ //if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
                //diagVec = tauVal0/m_mu2;
                diagVec = 1/m_mu2;
                //tauind = tauind + 1;
                //diagVecG = spV.diag();
                //diagVecG_I = diagVecG.elem(g_I_longl_vec);
                //diagVec = diagVec + tauVec(indforV) * diagVecG_I;
                diagVec = diagVec + m_tauVal_sp;

                /*
                tauind = tauind + 1 + 3 * (indforV - 1);
                tauind = tauind + 1;

                if(g_T_longl_mat.n_rows > 0){
                  diagVecG_I.ones(Nnomissing);
                  diagVecG_IT = diagVecG_I % g_T_longl_vec;
                  diagVecG_T = diagVecG_IT % g_T_longl_vec;
                  diagVecG_IT = 2 * diagVecG_IT;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_IT;
                  tauind = tauind + 1;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_T;
                  tauind = tauind + 1;
                }
                */
        }

        for(unsigned int i=0; i< Nnomissing; i++){
                if(diagVec(i) < 1e-4){
                        diagVec(i) = 1e-4 ;
                }
        }

        return(diagVec);
}







arma::vec SAIGEClass::getCrossprod_V(arma::vec& bVec){
        //indforV starts from 0
        arma::vec crossProdVec;
        arma::vec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;
	double tauVal0 = m_tauvec(0); 

        unsigned int tau_ind = 0;
        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){ //it must have specified GRM
                  crossProd1 = bVec;
                  //crossProdVec = tauVal0*(bVec % (1/m_mu2)) + m_tauVal_sp*crossProd1;
                  crossProdVec = (bVec % (1/m_mu2)) + m_tauVal_sp*crossProd1;
                  tau_ind = tau_ind + 2;
        }else{

                  Ibvec = g_I_longl_mat.t() * bVec;

                  GRM_I_bvec = Ibvec;
                  crossProd1 = g_I_longl_mat * GRM_I_bvec;
                  //crossProdVec = tauVal0*(bVec % (1/m_mu2)) + m_tauVal_sp*crossProd1;
                  crossProdVec = (bVec % (1/m_mu2)) + m_tauVal_sp*crossProd1;
                /*
                  tau_ind = tau_ind + 1 + 3 * (indforV - 1);
                  tau_ind = tau_ind + 1;

                if(g_T_longl_mat.n_rows > 0){
                        Tbvec = g_T_longl_mat.t() * bVec;
                        GRM_T_bvec = spV * Tbvec;
                        crossProdGRM_TGIb = g_T_longl_mat * GRM_I_bvec;
                        crossProdGRM_IGTb = g_I_longl_mat * GRM_T_bvec;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (crossProdGRM_TGIb + crossProdGRM_IGTb);
                        tau_ind = tau_ind + 1;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (g_T_longl_mat * GRM_T_bvec);
                        tau_ind = tau_ind + 1;
                }
                */

        }

        return(crossProdVec);
}



arma::vec SAIGEClass::getPCG1ofSigmaAndVector_V(arma::vec& bVec, int maxiterPCG, double tolPCG){
    // Start Timers
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    int Nnomissing = m_mu2.n_elem;
    arma::vec xVec(Nnomissing);
    xVec.zeros();

    //if(g_isStoreSigma){
    //    std::cout << " arma::spsolve(g_spSigma, bVec) 0" << std::endl;
    //    //xVec = arma::spsolve(g_spSigma, bVec);
    //    xVec = arma::spsolve(g_spSigma_V, bVec);
    //    std::cout << " arma::spsolve(g_spSigma, bVec) 1" << std::endl;
    //}else{
        arma::vec rVec = bVec;
        arma::vec r1Vec;
        arma::vec crossProdVec(Nnomissing);
        arma::vec zVec(Nnomissing);
        arma::vec minvVec(Nnomissing);
        minvVec = 1/getDiagOfSigma_V();
        zVec = minvVec % rVec;

        float sumr2 = sum(rVec % rVec);
        arma::vec z1Vec(Nnomissing);
        arma::vec pVec = zVec;

        int iter = 0;

        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                arma::colvec ApVec = getCrossprod_V(pVec);
                arma::vec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                double a = preA(0);
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
                z1Vec = minvVec % r1Vec;

                arma::vec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                double bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
        }

        if (iter >= maxiterPCG){
		std::cout << "pcg did not converge. You may increase maxiter number." << std::endl;

        }
	std::cout << "iter from getPCG1ofSigmaAndVector " << iter << std::endl;
//}
        return(xVec);
}


arma::vec SAIGEClass::getSigma_G_V(arma::vec & bVec, int maxiterPCG, double tolPCG){
	arma::vec Sigma_iG;
        Sigma_iG = getPCG1ofSigmaAndVector_V(bVec, maxiterPCG, tolPCG);
        return(Sigma_iG);
}


double SAIGEClass::EmpSPA_K_0(double t, 
             int N0, 
             double adjG0, 
             arma::vec & adjG1)        // adjusted Genotype 
{
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * m_K_0_emp.getValue(t_adjG0) + arma::sum(m_K_0_emp.getVector(t_adjG1));
    return out;
}

double SAIGEClass::EmpSPA_K_1(double t,
             int N0, 
             double adjG0, 
             arma::vec & adjG1,        // adjusted Genotype
             double q2)
{
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * adjG0 * m_K_1_emp.getValue(t_adjG0) + arma::sum(adjG1 % m_K_1_emp.getVector(t_adjG1)) - q2;
    return out;
}


double SAIGEClass::EmpSPA_K_2(double t, 
             int N0, 
             double adjG0, 
             arma::vec & adjG1)       // adjusted Genotype
{
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * pow(adjG0, 2) * m_K_2_emp.getValue(t_adjG0) + arma::sum(pow(adjG1, 2) % m_K_2_emp.getVector(t_adjG1));
    return out;
}

Rcpp::List SAIGEClass::EmpSPA_fastgetroot_K1(double t_initX,
                            int N0, 
                            double adjG0, 
                            arma::vec & adjG1,        // adjusted Genotype
                            double q2)
  {
    double x = t_initX, oldX;
    double K1 = 0, K2 = 0, oldK1;
    double diffX = arma::datum::inf, oldDiffX;
    bool converge = true;
    double tol = 0.001;
    int maxiter = 100;
    int iter = 0;
    
    for(iter = 0; iter < maxiter; iter ++){
      
      oldX = x;
      oldDiffX = diffX;
      oldK1 = K1;
      
      K1 = EmpSPA_K_1(x, N0, adjG0, adjG1, q2);
      K2 = EmpSPA_K_2(x, N0, adjG0, adjG1);
      
      diffX = -1 * K1 / K2;
      
      // Checked on 03/25/2021: Expected!!!!
      // std::cout << "x:\t" << x << std::endl;
      // std::cout << "K1:\t" << K1 << std::endl;
      // std::cout << "arma::sign(K1):\t" << arma::sign(K1) << std::endl;
      // std::cout << "arma::sign(oldK1):\t" << arma::sign(oldK1) << std::endl;
      // std::cout << "K2:\t" << K2 << std::endl;
      // std::cout << "diffX:\t" << diffX << std::endl;
      
      if(!std::isfinite(K1)){
        // checked it on 07/05:
        // if the solution 'x' tends to infinity, 'K2' tends to 0, and 'K1' tends to 0 very slowly.
        // then we can set the one sided p value as 0 (instead of setting converge = F)
        x = arma::datum::inf;
        K2 = 0;
        break;
      }
      
      if(arma::sign(K1) != arma::sign(oldK1)){
        while(std::abs(diffX) > std::abs(oldDiffX) - tol){
          diffX = diffX / 2;
        }
      }
      
      if(std::abs(diffX) < tol) break;
      
      x = oldX + diffX;
    }
    
    if(iter == maxiter) 
      converge = false;

    Rcpp::List yList = Rcpp::List::create(Rcpp::Named("root") = x,
                                          Rcpp::Named("iter") = iter,
                                          Rcpp::Named("converge") = converge,
                                          Rcpp::Named("K2") = K2);
    return yList;
  }

double SAIGEClass::EmpSPA_GetProb_SPA(double adjG0, 
                     arma::vec & adjG1,
                     int N0, 
                     double q2, 
                     bool lowerTail)
  {
    double initX = 0;
    
    // The following initial values are validated on 03/25/2021
    if(q2 > 0) initX = 3;
    if(q2 <= 0) initX = -3;
    
    Rcpp::List rootList = EmpSPA_fastgetroot_K1(initX, N0, adjG0, adjG1, q2);
    double zeta = rootList["root"];
    
    double k1 = EmpSPA_K_0(zeta,  N0, adjG0, adjG1);
    double k2 = EmpSPA_K_2(zeta,  N0, adjG0, adjG1);
    double temp1 = zeta * q2 - k1;
    
    double w = arma::sign(zeta) * sqrt(2 * temp1);
    double v = zeta * sqrt(k2);
    
    double pval = arma::normcdf(arma::sign(lowerTail-0.5) * (w + 1/w * log(v/w)));

    return pval;
  }

double SAIGEClass::EmpSPA_getMarkerPval(arma::vec & t_g,
                       double t_Tstat,
		       double g_altFreq_new, 
		       int g_N)
{

    // estimated variance after adjusting for covariates

    arma::vec adjGVec2 = pow(t_g, 2);
    double VarS = m_varResid * sum(adjGVec2) * m_varRatioVal;
    double VarS_rawG = m_varResid*2*g_altFreq_new*(1-g_altFreq_new)*g_N * m_varRatioVal; 
    //double t_zScore = t_Tstat / sqrt(VarS);
    double t_zScore = t_Tstat;
    //arma::vec adjGVecNorm = t_g / sqrt(VarS);
    arma::vec adjGVecNorm = t_g / sqrt(sum(adjGVec2));
    double pval_noEmpSPA = arma::normcdf(-1*std::abs(t_zScore))*2;
    std::cout << "pval_noEmpSPA " << pval_noEmpSPA << std::endl;
    std::cout << "t_Tstat " << t_Tstat << std::endl;
    std::cout << "VarS " << VarS << std::endl;
    std::cout << "VarS_rawG " << VarS_rawG << std::endl;
    std::cout << "m_varResid " << m_varResid << std::endl;
    std::cout << "m_varRatioVal " << m_varRatioVal << std::endl;

    double varg = arma::mean(t_g);
    std::cout << "varg " << varg << std::endl;

    int N0 = 0;
    double adjG0 = 0;   // since N0=0, this value actually does not matter

    double pval1 = EmpSPA_GetProb_SPA(adjG0, adjGVecNorm, N0, std::abs(t_zScore), false);
    double pval2 = EmpSPA_GetProb_SPA(adjG0, adjGVecNorm, N0, -1*std::abs(t_zScore), true);
    double pval = pval1 + pval2;

    return pval;
  }


}
