void SAIGEClass::getMarkerPval_multi(arma::mat & t_GMat,
                               arma::uvec & iIndex,
                               arma::uvec & iIndexComVec,
			       std::vector <double> & t_Beta,
                               std::vector <double> & t_seBeta,
                               std::vector <double> & t_pval,
                               std::vector <double> & t_pval_noSPA,
                               arma::vec & t_altFreq_vec,
                               std::vector <double> & t_Tstat,
                               std::vector <double> & t_gy,
                               std::vector <double> & t_var1,
			       std::vector <bool> & t_isSPAConverge,
			       arma::mat & t_gtilde,
			       std::vector <bool> & is_gtilde,
                               bool  is_region,
                               arma::vec & t_P2Vec,
                               bool t_isCondition,
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
  t_isFirth = false;
  std::vector <std::string>  t_pval_str_vec;
  arma::vec t_var2_vec; 
  double t_pval, t_SPApval;
  arma::vec altFreq0_vec, t_pval_noadj_vec;
  arma::ivec t_skipSPAind_vec;
  std::cout << "t_isnoadjCov " << t_isnoadjCov << std::endl;
  is_gtilde = false;
  scoreTestFast_noadjCov_multi(t_GMat, t_Beta_vec, t_seBeta_vec, t_pval_str_vec, altFreq0_vec, t_Tstat_vec, t_var1_vec, t_var2_vec, t_skipSPAind_vec, t_pval_noadj_vec);

  double StdStat = std::abs(t_Tstat) / sqrt(t_var1);
  t_isSPAConverge = false;
  double altFreq0, pval_noadj;
  arma::vec t_GVec0;
  arma::uvec indexNonZeroVec0_arma, indexZeroVec0_arma;
  int numMarkers = t_var1.n_elem; 
  
  for(int i = 0; i < t_skipSPAind_vec.n_elem; i++){
      uint j =  t_skipSPAind_vec(i);
      if((g_I_longl_mat.n_cols != g_I_longl_mat.n_rows) && (t_GVec.n_elem < m_y.n_elem)){
          t_GVec0 = g_I_longl_mat * t_GMat.col(j);
          altFreq0  = arma::mean(t_GVec0, 0) /2; //column mean
      }else{
          t_GVec0 = t_GMat.col(j);
          altFreq0 = altFreq0_vec(j);
      }
      indexNonZeroVec0_arma = arma::find(t_GVec0 > 0.0);
      indexZeroVec0_arma = arma::find(t_GVec0 == 0.0);
      pval_noadj = t_pval_noadj_vec(j);

      double q, qinv, m1, NAmu, NAsigma, tol1, p_iIndexComVecSize;
      unsigned int iIndexComVecSize = iIndexComVec.n_elem;
      unsigned int iIndexSize = iIndex.n_elem;
      arma::vec gNB(iIndexSize, arma::fill::none);
      arma::vec gNA(iIndexComVecSize, arma::fill::none);
      arma::vec muNB(iIndexSize, arma::fill::none);
      arma::vec muNA(iIndexComVecSize, arma::fill::none);
      double gmuNB;

      if(!is_gtilde){
          t_gtilde.resize(m_n);
          getadjGFast(t_GVec0, t_gtilde, iIndex);
          is_gtilde = true;
      }

        p_iIndexComVecSize = double(iIndexComVecSize)/m_n;
        m1 = dot(m_mu, t_gtilde);

        if(p_iIndexComVecSize >= 0.5){
                unsigned int j1 = 0;
                unsigned int j2 = 0;
        	gNB = t_gtilde(iIndex);
        	gNA = t_gtilde(iIndexComVec);
        	muNB = m_mu(iIndex);
        	muNA = m_mu(iIndexComVec);
        	gmuNB = dot(gNB,muNB);
        	NAmu= m1-gmuNB;
        }

        if(m_traitType == "survival" || m_traitType == "count"){
                q = (t_Tstat(j)) /sqrt((t_var1(j))/(t_var2(j)));
                qinv = -q;
                if(p_iIndexComVecSize >= 0.5){
                        NAsigma = t_var2(j) - arma::sum(muNB % arma::pow(gNB,2));
                }
        }
        bool logp=false;
        double tol0 = std::numeric_limits<double>::epsilon();
        tol1 = std::pow(tol0, 0.25);
        if(p_iIndexComVecSize >= 0.5){
                SPA_fast(m_mu, t_gtilde, q, qinv, pval_noadj, false, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge);
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
   t_pval_noSPA = pval_noadj;
   if(m_traitType!="quantitative"){
        if(t_isSPAConverge){
                t_pval = t_SPApval;
        }else{
                t_pval = pval_noadj;
        }
   }else{
        t_pval = t_pval_noSPA;
   }
}

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

