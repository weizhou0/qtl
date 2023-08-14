void Unified_getMarkerPval_multi(
                           arma::mat & t_GMat,
                           bool t_isOnlyOutputNonZero,
                           arma::uvec & t_indexForNonZero_vec,
                           arma::uvec & t_indexForZero_vec,
                           std::vector <double> & t_Beta,
                           std::vector <double> & t_seBeta,
                           std::vector <double> & t_pval,
                           std::vector <double> & t_pval_noSPA,
                           std::vector <double> & t_Tstat,
                           std::vector <double> & t_gy,
                           std::vector <double> & t_varT,
                           std::vector <double> t_altFreq,
                           std::vector <bool>  & t_isSPAConverge,
                           arma::vec & t_gtilde,
                           std::vector <bool>& is_gtilde,
                           bool  is_region,
                           arma::vec & t_P2Vec,
                           bool  t_isCondition,
                           std::vector <double> & t_Beta_c,
                           std::vector <double> & t_seBeta_c,
                           std::vector <double> & t_pval_c,
                           std::vector <double> & t_pval_noSPA_c,
                           std::vector <double> & t_Tstat_c,
                           std::vector <double> & t_varT_c,
                           arma::rowvec & t_G1tilde_P_G2tilde_Vec,
                           bool & t_isFirth,
                           bool & t_isFirthConverge,
                           bool t_isER,
                           bool t_isnoadjCov,
                           bool t_isSparseGRM)
{
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SAIGE method to calculate marker-level p-values, 't_isOnlyOutputNonZero' should be false.");


    ptr_gSAIGEobj->getMarkerPval_multi(t_GMat, t_indexForNonZero_vec, t_indexForZero_vec, t_Beta, t_seBeta, t_pval, t_pval_noSPA, t_altFreq, t_Tstat, t_gy, t_varT, t_isSPAConverge, t_gtilde, is_gtilde, is_region, t_P2Vec, t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c, t_pval_noSPA_c, t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde_Vec, t_isFirth, t_isFirthConverge, t_isER, t_isnoadjCov, t_isSparseGRM); //SAIGE_new.cpp

    //t_indexForNonZero_vec.clear();

}
