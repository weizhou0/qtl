
#ifndef SAIGE_HPP
#define SAIGE_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "approxfun.hpp"

namespace SAIGE{

class SAIGEClass
{
    private:
      arma::mat m_XVX;
      arma::mat m_XVX_inv_XV;
      arma::mat m_X;
      arma::mat m_Sigma_iXXSigma_iX;
      arma::vec m_res;
      arma::vec m_res_sample;
      arma::vec m_resout;
      arma::vec m_mu;
      arma::vec  m_S_a;
      std::string m_traitType; 
      std::string m_impute_method;
      std::vector<uint32_t> m_condition_genoIndex;

 


    public:
      arma::vec m_mu2;
      arma::vec m_mu2_sample;
      arma::vec m_tauvec;
      arma::vec m_varWeightsvec; 
      arma::mat m_XXVX_inv;
      arma::mat m_XV;
      int m_n, m_p; //MAIN Dimensions: sample size, number of covariates
      double m_varRatioVal;
      arma::vec m_varRatio_sparse;
      arma::vec m_varRatio_null;
      arma::vec m_varRatio_null_noXadj;
      arma::vec m_varRatio_null_eg;
      arma::vec m_varRatio_sparse_eg;
      arma::vec m_y;

      bool m_isOutputAFinCaseCtrl;
      bool m_isOutputNinCaseCtrl;
      bool m_isOutputHetHomCountsinCaseCtrl;
      arma::uvec m_case_indices;
      arma::uvec m_ctrl_indices;
      arma::uvec m_case_hom_indices;
      arma::uvec m_case_het_indices;
      arma::uvec m_ctrl_hom_indices;
      arma::uvec m_ctrl_het_indices;
      int m_n_case;
      int m_n_ctrl;
      //arma::sp_mat m_SigmaMat_sp;
      bool m_flagSparseGRM;
      bool m_flagSparseGRM_cur;
      bool m_isFastTest;
      bool m_isnoadjCov;      
      bool m_is_EmpSPA;

      double m_pval_cutoff_for_fastTest;
      double m_SPA_Cutoff;
      arma::umat m_locationMat;
      arma::vec m_valueVec;
      int m_dimNum;	
      arma::vec m_cateVarRatioMinMACVecExclude; 
      arma::vec m_cateVarRatioMaxMACVecInclude;
      arma::mat m_P2Mat_cond;
      int m_numMarker_cond;
      arma::mat m_VarInvMat_cond;
      arma::mat m_VarMat_cond;
      arma::vec m_Tstat_cond;
      arma::vec m_G2_Weight_cond;
      arma::vec m_MAF_cond;
      double  m_qsum_cond;
      arma::vec m_gsum_cond;
      arma::vec m_p_cond;
      arma::vec m_scalefactor_G2_cond;
      arma::mat m_VarInvMat_cond_scaled_weighted;
      //arma::mat m_VarInvMat_cond_region_binary;
      bool m_isCondition;
      bool m_is_Firth_beta;
      double m_pCutoffforFirth;
     arma::vec  m_offset;	
      bool m_isVarPsadj;

      arma::vec m_sigmainvG_noV;	
      arma::sp_mat m_SigmaMat_sp;
      arma::sp_mat g_I_longl_mat;
      arma::sp_mat g_T_longl_mat;
      arma::uvec g_I_longl_vec;
      arma::vec g_T_longl_vec;
      double m_tauVal_sp;	
      //arma::m_var2m;
  ////////////////////// -------------------- functions ---------------------------------- //////////////////////

      approxfun::approxfunClass m_K_0_emp;
      approxfun::approxfunClass m_K_1_emp;
      approxfun::approxfunClass m_K_2_emp;
      arma::mat m_cumul;
      double m_varResid;

SAIGEClass(
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
	arma::vec & t_varRatio_null_eg,
	arma::vec & t_varRatio_sparse_eg,
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
        arma::mat & t_cumul);


   void set_seed(unsigned int seed);

   void scoreTest(arma::vec & t_GVec,
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
		     arma::uvec & t_indexForNonZero);

    void scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2);

    void scoreTestFast_noadjCov(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2);	

     void set_flagSparseGRM_cur(bool t_flagSparseGRM_cur);

     void get_mu(arma::vec & t_mu);

     void getadjG(arma::vec & t_GVec, arma::vec & g);
     void getadjGFast(arma::vec & t_GVec, arma::vec & g,  arma::uvec & iIndex);

     void getMarkerPval(arma::vec & t_GVec,
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
				bool t_isSparseGRM);


    void getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices);


    void setupSparseMat(int r, arma::umat & locationMatinR, arma::vec & valueVecinR);

    arma::sp_mat gen_sp_SigmaMat();


    bool assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj);

    void assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj);


    void assignSingleVarianceRatio_withinput(double t_varRatioVal);


    void assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
            arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
       arma::vec & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      double t_qsum_cond,
      arma::vec & t_gsum_cond,
      arma::vec & t_p_cond);

     void assignConditionFactors_scalefactor(
        arma::vec & t_scalefactor_G2_cond);	


    void extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv);

    void fast_logistf_fit_simple(arma::mat & x,
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
	bool & isfirthconverge);	

     arma::vec getSigma_G_V(arma::vec & bVec, int maxiterPCG, double tolPCG);
     arma::vec getPCG1ofSigmaAndVector_V(arma::vec & bVec, int maxiterPCG, double tolPCG);

     arma::vec getCrossprod_V(arma::vec& bVec);

     arma::vec getDiagOfSigma_V();

     double EmpSPA_K_0(double t,
             int N0,
             double adjG0,
             arma::vec & adjG1); 	


     double EmpSPA_K_1(double t,
             int N0,
             double adjG0,
             arma::vec & adjG1,        // adjusted Genotype
             double q2);


     double EmpSPA_K_2(double t,
             int N0,
             double adjG0,
             arma::vec & adjG1);


     Rcpp::List EmpSPA_fastgetroot_K1(double t_initX,
                            int N0,
                            double adjG0,
                            arma::vec & adjG1,        // adjusted Genotype
                            double q2);

     double EmpSPA_GetProb_SPA(double adjG0,
                     arma::vec & adjG1,
                     int N0,
                     double q2,
                     bool lowerTail);


     double EmpSPA_getMarkerPval(arma::vec & t_g,
                       double t_Tstat,
		       double g_altFreq_new,
                       int g_N);

     void getMarkerPval_multi(arma::mat & t_GMat_center,
                               std::vector <double> & t_Beta,
                               std::vector <double> & t_seBeta,
                               std::vector <double> & t_pval,
                               std::vector <double> & t_pval_noSPA,
                               std::vector <double> & t_Tstat,
                               std::vector <double> & t_var1,
                               std::vector <std::string> & t_pval_str_vec);


     void scoreTestFast_noadjCov_multi(arma::mat & t_GMat_centered,
                     arma::vec & t_Beta,
                     arma::vec & t_seBeta,
                     std::vector <std::string> & t_pval_str,
                     arma::vec &t_Tstat,
                     arma::vec &t_var1,
                     arma::vec &t_var2,
                     arma::uvec & t_skipSPAind_vec,
                     arma::vec & t_pval_noadj_vec); 



};
}
#endif
