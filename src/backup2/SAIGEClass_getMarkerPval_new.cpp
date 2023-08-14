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
  std::cout << "herehrere" << std::endl;
  t_isFirth = false;
  std::string t_pval_str;
  double t_var2, t_SPApval;
  arma::vec altFreq0_vec;
  

  arma::mat t_GMat0;
  arma::uvec indexNonZeroVec0_arma, indexZeroVec0_arma;
  if((g_I_longl_mat.n_cols != g_I_longl_mat.n_rows) && (t_GVec.n_elem < m_y.n_elem)){
      t_GMat0 = g_I_longl_mat * t_GMat;
      altFreq0_vec = arma::mean(t_GMat0, 0) /2; //column mean
      //indexNonZeroVec0_arma = arma::find(t_GVec0 > 0.0);
      //indexZeroVec0_arma = arma::find(t_GVec0 == 0.0);
   }else{
      t_GMat0 = t_GMat;
      altFreq0_vec = t_altFreq_vec;
      //indexNonZeroVec0_arma = iIndex;
      //indexZeroVec0_arma = iIndexComVec;
   }

 //arma::vec timeoutput3 = getTime();
 //
 //
 if(t_isSparseGRM){
        t_isnoadjCov = false;
 }
std::cout << "t_isnoadjCov " << t_isnoadjCov << std::endl;

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
          scoreTestFast(t_GVec0, indexNonZeroVec0_arma, t_Beta, t_seBeta, t_pval_str, altFreq0, t_Tstat, t_var1, t_var2);
        }else{
          is_gtilde = true;
          scoreTest(t_GVec0, t_Beta, t_seBeta, t_pval_str, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, indexNonZeroVec0_arma);
        }
  }else{
        is_gtilde = false;
        //unsigned int nonzero = iIndex.n_elem;
        //unsigned int nGvec = t_GVec.n_elem;
       scoreTestFast_noadjCov_multi(t_GMat, t_Beta, t_seBeta, t_pval_str, altFreq0,t_Tstat, t_var1, t_var2);
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
//      std::cout << "gNB.n_elem " <<  gNB.n_elem << std::endl;
//      std::cout << "gNA.n_elem " <<  gNA.n_elem << std::endl;
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

