// [[Rcpp::export]]
void mainMarkerInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
                           std::string & t_traitType,
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           bool & t_isMoreOutput,
                           bool & t_isImputation,
                           bool & t_isFirth)
{

  std::cout << "Here1 mainMarkerInCPP" << std::endl;
  //std::cout << "ptr_gSAIGEobj->m_flagSparseGRM_cur " << ptr_gSAIGEobj->m_flagSparseGRM_cur << std::endl;
  //std::cout << "ptr_gSAIGEobj->m_flagSparseGRM " << ptr_gSAIGEobj->m_flagSparseGRM << std::endl;


  int q = t_genoIndex.size();  // number of markers
  // set up output
  std::vector<std::string> markerVec(q);  // marker IDs
  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs

  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q);
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);
  std::vector<double> pvalVec(q, arma::datum::nan);
  std::vector<double> TstatVec(q, arma::datum::nan);
  std::vector<double> varTVec(q, arma::datum::nan);
  std::vector<double> pvalNAVec(q, arma::datum::nan);

  bool isCondition = ptr_gSAIGEobj->m_isCondition;
  //if(isCondition){
  std::vector<double> Beta_cVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBeta_cVec(q, arma::datum::nan);
  std::vector<double> pval_cVec(q, arma::datum::nan);
  std::vector<double> Tstat_cVec(q, arma::datum::nan);
  std::vector<double> varT_cVec(q, arma::datum::nan);
  std::vector<double> pvalNA_cVec(q, arma::datum::nan);
  //}

  std::vector<std::string> Beta_ge_cStrVec(q, "NA");         // beta value for ALT allele
  std::vector<std::string> seBeta_ge_cStrVec(q, "NA");
  std::vector<std::string> pval_ge_cStrVec(q, "NA");
  std::vector<std::string> pval_noSPA_ge_cStrVec(q, "NA");
  std::vector<double> pval_SKATO_ge_cVec(q, arma::datum::nan);
  std::vector<double> pval_SKAT_ge_cVec(q, arma::datum::nan);
  //std::vector<std::string> Tstat_ge_cStrVec(q, arma::datum::nan);
  //std::vector<std::string> varT_ge_cStrVec(q, arma::datum::nan);



  arma::rowvec G1tilde_P_G2tilde_Vec(ptr_gSAIGEobj->m_numMarker_cond);

  std::vector<bool>  isSPAConvergeVec(q);
  std::vector<double>  AF_caseVec(q);
  std::vector<double>  AF_ctrlVec(q);
  std::vector<uint32_t>  N_caseVec(q);
  std::vector<uint32_t>  N_ctrlVec(q);
    //if(t_isMoreOutput){
  std::vector<double>  N_case_homVec(q);
  std::vector<double>  N_ctrl_hetVec(q);
  std::vector<double>  N_case_hetVec(q);
  std::vector<double>  N_ctrl_homVec(q);
  std::vector<uint32_t>  N_Vec(q);
  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;
  std::vector<uint> indexForMissing;

  //int n = ptr_gSAIGEobj->m_n;
  int n = ptr_gSAIGEobj->m_n;

  if(g_n_unique > 0){
    n = g_n_unique;
  }

  arma::vec t_GVec(n);
  arma::vec gtildeVec(n);
 arma::vec t_P2Vec;
  if(ptr_gSAIGEobj->m_isFastTest){
    ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
  }else{
    ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
  }

  bool hasVarRatio = true;;
  bool isSingleVarianceRatio = true;
  if((ptr_gSAIGEobj->m_varRatio_null).n_elem == 1){
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
        //ptr_gSAIGEobj->assignSingleVarianceRatio(false);
  }else{
        isSingleVarianceRatio = false;
  }


  int mFirth = 0;
  int mFirthConverge = 0;

  for(int i = 0; i < q; i++){
    if((i+1) % g_marker_chunksize == 0){
      std::cout << "Completed " << (i+1) << "/" << q << " markers in the chunk." << std::endl;
    }
    // information of marker
    double altFreq, altCounts, missingRate, imputeInfo, AF_case, AF_ctrl, N_case_hom, N_ctrl_het, N_case_het, N_ctrl_hom;
    std::string chr, ref, alt, marker;
    uint32_t pd, N_case, N_ctrl, N;

    //free(end);

    bool flip = false;
    std::string t_genoIndex_str = t_genoIndex.at(i);
    char* end;
    uint64_t gIndex = std::strtoull( t_genoIndex_str.c_str(), &end,10 );
    std::remove(end);

    uint64_t gIndex_prev = 0;
    if(i == 0){
        gIndex_prev = 0;
    }else{
        char* end_prev;
        std::string t_genoIndex_prev_str;
        if(t_genoType == "bgen"){
            t_genoIndex_prev_str = t_genoIndex_prev.at(i-1);
        }else if(t_genoType == "plink"){
            t_genoIndex_prev_str = t_genoIndex.at(i-1);
        }
        gIndex_prev = std::strtoull( t_genoIndex_prev_str.c_str(), &end_prev,10 );
        std::remove(end_prev);
    }




    //Main.cpp
    //PLINK or BGEN
    //uint32_t gIndex_temp = gIndex;
    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false;

   //clear vectors
   indexZeroVec.clear();
   indexNonZeroVec.clear();
   indexForMissing.clear();

   bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec, t_GVec, t_isImputation);
   //arma::vec timeoutput2 = getTime();
   //printTime(timeoutput1, timeoutput2, "Unified_getOneMarker");
//
   if(!isReadMarker){
      //std::cout << "isReadMarker " << isReadMarker << std::endl;
      g_markerTestEnd = true;
      bool isEndFile = check_Vcf_end();
      break;
    }

    std::string pds = std::to_string(pd);
    std::string info = chr+":"+pds+":"+ref+":"+alt;
    //std::cout << "info " << info << std::endl;

    chrVec.at(i) = chr;
    posVec.at(i) = pds;
    refVec.at(i) = ref;
    altVec.at(i) = alt;
    // record basic information for the marker
    markerVec.at(i) = marker;               // marker IDs
    infoVec.at(i) = info;    // marker information: CHR:POS:REF:ALT
   altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    //altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    imputationInfoVec.at(i) = imputeInfo;



    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    int nG = t_GVec.n_elem;
    double MAC = MAF * n * (1 - missingRate) *2;

    // Quality Control (QC) based on missing rate, MAF, and MAC
    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff || imputeInfo < g_marker_minINFO_cutoff)){
      continue;
    }else{

    // Check UTIL.cpp
    //
    //
    //arma::vec timeoutput3 = getTime();
    indexZeroVec.clear();
    indexNonZeroVec.clear();

    flip = imputeGenoAndFlip(t_GVec, altFreq, altCounts,indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);

//arma::vec timeoutput4 = getTime();
//printTime(timeoutput3, timeoutput4, "imputeGenoAndFlip");
    altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    MAC = std::min(altCounts, 2*n-altCounts);

    // analysis results for single-marker
    double Beta, seBeta, pval, pval_noSPA, Tstat, varT, gy;
    double Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c;

    bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
    //arma::vec t_P2Vec;
    //arma::vec t_P2Vec;
    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;


    //std::cout << "g_n_unique " << g_n_unique << std::endl;
    indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
    indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);

    indexZeroVec.clear();
    indexNonZeroVec.clear();
    t_P2Vec.clear();
    G1tilde_P_G2tilde_Vec.clear();
    //arma::vec timeoutput5 = getTime();
    //std::cout << "Here1a" << std::endl;
    //set_varianceRatio(MAC, isSingleVarianceRatio);

    //ptr_gSAIGEobj->set_flagSparseGRM_cur(m_isSparseGRM);

    //if(ptr_gSAIGEobj->m_isnoadjCov){
    //}



    ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
    if(isSingleVarianceRatio){
       ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
    }else{
       hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
    }
    //check 'Main.cpp'
    bool is_region = false;

