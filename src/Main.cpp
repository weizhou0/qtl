#define ARMA_USE_SUPERLU 1
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <vector>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
// std::this_thread::sleep_for (std::chrono::seconds(1));
#include <cstdio>         // std::remove
#include <fstream>
#include <string.h>
// Currently, omp does not work well, will check it later
// error: SET_VECTOR_ELT() can only be applied to a 'list', not a 'character'
// remove all Rcpp::List to check if it works
// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]]

#include "Main.hpp"
#include "PLINK.hpp"
#include "BGEN.hpp"
#include "VCF.hpp"
#include "SAIGE_test.hpp"
#include "UTIL.hpp"
#include "CCT.hpp"
#include "GENO_null.hpp"
#include "SKATO.hpp"


#include <Rcpp.h>
#include "getMem.hpp"

#include <boost/math/distributions/beta.hpp>

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;
static VCF::VcfClass* ptr_gVCFobj = NULL;
static SAIGE::SAIGEClass* ptr_gSAIGEobj = NULL;
static NullGENO::NullGenoClass* ptr_gNULLGENOobj = NULL;
//single, SAIGE
//Region, SAIGE-GENE+
//Step 1: NULL model
//


//Step 1
bool isUseSparseSigmaforModelFitting = false;
std::vector<arma::sp_fmat> Kmat_vec;
arma::sp_fmat g_I_longl_mat;
arma::sp_fmat g_T_longl_mat;
arma::uvec g_I_longl_vec;
arma::fvec g_T_longl_vec;
arma::sp_fmat g_spGRM;
arma::sp_fmat g_spSigma;
arma::sp_fmat g_spSigma_noV;
arma::sp_fmat g_spSigma_V;
bool g_isSparseGRM;
bool g_isStoreSigma;
int g_num_Kmat;
bool g_isGRM;
arma::umat g_covarianceidxMat;
arma::uvec g_covarianceidxMat_col1;
arma::uvec g_covarianceidxMat_col2;
arma::uvec g_covarianceidxMat_col3;
arma::uvec g_covarianceidxMat_notcol1;
//arma::fvec g_var_weights;



//Step 2
// global variables for analysis
std::string g_impute_method;      // "mean", "minor", or "drop", //drop is not allowed
double g_missingRate_cutoff;
//unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_minMAC_cutoff;    // for Rare Variants (RVs) whose MAC < this value, we aggregate these variants like SAIGE-GENE+ 
double g_marker_minINFO_cutoff;
arma::vec g_region_maxMAF_cutoff;
double g_min_gourpmac_for_burdenonly;
double g_maxMAFLimit;
unsigned int g_region_maxMarkers_cutoff;   // maximal number of markers in one chunk, only used for region-based analysis to reduce memory usage
bool g_isOutputMoreDetails;
int g_marker_chunksize;

std::string g_method_to_CollapseUltraRare;
double g_DosageCutoff_for_UltraRarePresence;

double g_dosage_zerod_MAC_cutoff;
double g_dosage_zerod_cutoff;
bool g_markerTestEnd = false;

//arma::vec g_weights_beta(2);

bool  g_is_Firth_beta;
double g_pCutoffforFirth;
double g_MACCutoffforER;

arma::uvec g_dup_sample_Index;
int g_n_unique = 0;

//std::ofstream OutFile;
//std::ofstream OutFile_singleInGroup;
//std::ofstream OutFile_single;
//std::ofstream OutFile_singleInGroup_temp;
std::vector<std::ofstream> OutFile_single_vec;
std::vector<std::ofstream> OutFile_vec;
std::vector<std::ofstream> OutFile_singleInGroup_vec;
std::vector<std::ofstream> OutFile_singleInGroup_temp_vec;


std::string g_outputFilePrefix0;
std::string g_outputFilePrefixGroup;
std::string g_outputFilePrefixSingleInGroup;
std::string g_outputFilePrefixSingleInGroup_temp;
std::string g_outputFilePrefixSingle;


arma::mat g_emat;
arma::fmat g_emat_f;
bool g_isgxe;
double g_pval_cutoff_for_gxe;

// [[Rcpp::export]]
void setAssocTest_GlobalVarsInCPP(std::string t_impute_method,
		double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
			       double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,
			       std::string t_outputFilePrefix,
			       double t_MACCutoffforER)
{
			       //arma::vec & t_weights_beta, 
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
  g_marker_minINFO_cutoff = t_min_info_marker;
  g_dosage_zerod_cutoff = t_dosage_zerod_cutoff;
  g_dosage_zerod_MAC_cutoff = t_dosage_zerod_MAC_cutoff;
  //g_weights_beta = t_weights_beta;
  g_outputFilePrefixGroup = t_outputFilePrefix;
  g_outputFilePrefixSingleInGroup = t_outputFilePrefix + ".singleAssoc.txt";
  g_outputFilePrefixSingleInGroup_temp = t_outputFilePrefix + ".singleAssoc.txt_temp";
  g_outputFilePrefixSingle = t_outputFilePrefix;
  g_MACCutoffforER = t_MACCutoffforER;
}


// [[Rcpp::export]]
void setAssocTest_GlobalVarsInCPP_GbyE(
                                arma::mat & t_emat,
                                bool t_isgxe, 
				double t_pval_cutoff_for_gxe,
                                arma::mat & t_XV_gxe,
                                arma::mat & t_XXVX_inv_gxe,
                                arma::mat & t_y_gxe,
                                arma::mat & t_res_gxe,
                                arma::mat & t_mu2_gxe,
                                arma::mat & t_mu_gxe,
                                arma::mat & t_varWeights_gxe)
{
	//g_emat_f=t_emat;
	g_emat=t_emat;
	//g_emat = arma::conv_to< arma::mat >::from(g_emat_f);
	g_isgxe=t_isgxe;
	g_pval_cutoff_for_gxe = t_pval_cutoff_for_gxe;
	//ptr_gSAIGEobj->m_XV_gxe_mt = t_XV_gxe;
	//ptr_gSAIGEobj->m_XXVX_inv_gxe_mt = t_XXVX_inv_gxe;
	//ptr_gSAIGEobj->m_y_gxe_mt = t_y_gxe;
	//ptr_gSAIGEobj->m_res_gxe_mt = t_res_gxe;
	//ptr_gSAIGEobj->m_mu2_gxe_mt = t_mu2_gxe;
	//ptr_gSAIGEobj->m_mu_gxe_mt = t_mu_gxe;
	//ptr_gSAIGEobj->m_varWeights_gxe_mt = t_varWeights_gxe;
}


// [[Rcpp::export]]
void setMarker_GlobalVarsInCPP(
			       bool t_isOutputMoreDetails,
			       int t_marker_chunksize
			       )

{
  g_isOutputMoreDetails = t_isOutputMoreDetails;
  g_marker_chunksize = t_marker_chunksize;
}


//double t_DosageCutoff_for_UltraRarePresence,
			       //std::string t_method_to_CollapseUltraRare,

  //g_method_to_CollapseUltraRare = t_method_to_CollapseUltraRare;
  //g_DosageCutoff_for_UltraRarePresence = t_DosageCutoff_for_UltraRarePresence;
// [[Rcpp::export]]
void setRegion_GlobalVarsInCPP(
                               arma::vec t_max_maf_region,
                               unsigned int t_max_markers_region,
			       double t_MACCutoff_to_CollapseUltraRare, 
			       double t_min_gourpmac_for_burdenonly)
{
  g_region_maxMAF_cutoff = t_max_maf_region;
  g_maxMAFLimit = g_region_maxMAF_cutoff.max();
  g_region_maxMarkers_cutoff = t_max_markers_region;
  g_region_minMAC_cutoff = t_MACCutoff_to_CollapseUltraRare;
  g_min_gourpmac_for_burdenonly = t_min_gourpmac_for_burdenonly;
}



//////// ---------- Main function for marker-level analysis --------- ////////////
                           //std::vector<uint32_t> & t_genoIndex,

// [[Rcpp::export]]
void mainMarkerInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
			   std::vector<std::string> & t_traitType,
			   std::vector<std::string> & t_genoIndex_prev,
			   std::vector<std::string> & t_genoIndex,
			   bool & t_isMoreOutput,
			   bool & t_isImputation,
			   bool & t_isFirth) 
{

  //std::cout << "Here1 mainMarkerInCPP" << std::endl;
  //std::cout << "t_traitType.size() here " << t_traitType.size() << std::endl;
  //std::cout << "ptr_gSAIGEobj->m_flagSparseGRM_cur " << ptr_gSAIGEobj->m_flagSparseGRM_cur << std::endl;
  //std::cout << "ptr_gSAIGEobj->m_flagSparseGRM " << ptr_gSAIGEobj->m_flagSparseGRM << std::endl;


  int q = (t_genoIndex.size()) * (t_traitType.size());  // number of markers
 unsigned int ne = (g_emat.n_cols)/(t_traitType.size());
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

 //std::cout << "n " << n << std::endl;

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

  std::cout << "ptr_gSAIGEobj->m_varRatio_null_mt.n_rows " << ptr_gSAIGEobj->m_varRatio_null_mt.n_rows << std::endl;
  if(ptr_gSAIGEobj->m_varRatio_null_mt.n_rows == 1){
        //ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
        //ptr_gSAIGEobj->assignSingleVarianceRatio(false);
	isSingleVarianceRatio = true;
  }else{		
	isSingleVarianceRatio = false;
  }


  int mFirth = 0;
  int mFirthConverge = 0;

  for(int i = 0; i < t_genoIndex.size(); i++){

//std::cout << "i " << i << std::endl;

    if((i+1) % g_marker_chunksize == 0){
      //std::cout << "Completed " << (i+1) << "/" << q << " markers in the chunk." << std::endl;
      std::cout << "Completed " << (i+1) << "/" << t_genoIndex.size() << " markers in the chunk." << std::endl;
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
   //t_GVec0.clear();
   //t_GVec.clear();
   std::cout << "Unified_getOneMarker " << std::endl;
   std::cout << "n " << n << std::endl;
//   std::cout << "t_GVec.n_elem " << t_GVec.n_elem << std::endl; 

   //t_GVec.set_size(n);
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


   //std::cout << "t_GVec0.size()) " << t_GVec0.size() << std::endl;
   //arma::vec t_GVec(t_GVec0.size());
   //arma::vec t_GVec = arma::conv_to< arma::colvec >::from(t_GVec0);

   //arma::vec t_GVec(t_GVec0);
   //t_GVec0.clear(); 

   //for(uint j = 0; j < n; j++){
   //	t_GVec(j) = t_GVec0.at(j);	
   //}

    //for(int indi = 0; indi < indexForNonZero.size(); indi++){
    //  std::cout << indexForNonZero[indi] << std::endl;
    //}
//   std::cout << "marker " << marker << std::endl;
//   std::cout << "indexForMissing.size() " << indexForMissing.size() << std::endl;
//   std::cout << "indexNonZeroVec.size() " << indexNonZeroVec.size() << std::endl;
    //int n = t_GVec.size();
    //arma::vec gtildeVec(n);

  


    std::string pds = std::to_string(pd); 
    std::string info = chr+":"+pds+":"+ref+":"+alt;
//    std::cout << "info " << info << std::endl;

  //std::cout << "Here2 mainMarkerInCPP" << std::endl;

 for(int i_mt0 = 0; i_mt0 < t_traitType.size(); i_mt0++){
    int j_mt0 = i_mt0*t_genoIndex.size()+i;
    chrVec.at(j_mt0) = chr;
    posVec.at(j_mt0) = pds;
    refVec.at(j_mt0) = ref;
    altVec.at(j_mt0) = alt; 
    // record basic information for the marker
    markerVec.at(j_mt0) = marker;               // marker IDs
    infoVec.at(j_mt0) = info;    // marker information: CHR:POS:REF:ALT
    altFreqVec.at(j_mt0) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    //altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(j_mt0) = missingRate;
    imputationInfoVec.at(j_mt0) = imputeInfo;
 }
//  std::cout << "Here2b mainMarkerInCPP" << std::endl;


    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    int nG = t_GVec.n_elem;
    double MAC = MAF * n * (1 - missingRate) *2;
   
   
   std::cout << "MAC " << MAC << std::endl;
   std::cout << "MAF " << MAF << std::endl;
   std::cout << "n " << n << std::endl;
   std::cout << "nG " << nG << std::endl;
 
    
   std::cout << "missingRate " << missingRate << std::endl;
   std::cout << "MAF " << MAF << std::endl;
   std::cout << "MAC " << MAC << std::endl;
   std::cout << "altFreq " << altFreq << std::endl;
   std::cout << "n " << n << std::endl;
   


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
   
  std::cout << "Here2c mainMarkerInCPP" << std::endl;
//arma::vec timeoutput4 = getTime();
//printTime(timeoutput3, timeoutput4, "imputeGenoAndFlip");
for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
    int j_mt = i_mt*t_genoIndex.size()+i;
    altFreqVec.at(j_mt) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(j_mt) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
//}    
    MAC = std::min(altCounts, 2*n-altCounts);
    
   /* 
    std::cout << "MAC " << MAC << std::endl; 
    std::cout << "altCounts " << altCounts << std::endl; 
   //std::cout << "altFreq after flip " << altFreq << std::endl; 
   //std::cout << "info " << info << std::endl;
   std::cout << "flip " << flip << std::endl; 
  */

    // analysis results for single-marker
    double Beta, seBeta, pval, pval_noSPA, Tstat, varT, gy;
    double Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c;

    bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
    //arma::vec t_P2Vec;
    //arma::vec t_P2Vec;
    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;


    std::cout << "g_n_unique " << g_n_unique << std::endl;
    indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
    indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);
/*
    arma::vec t_GVec0, t_GVec_cell;
    double altFreq_new, altCounts_new; 
    
    if(g_n_unique == 0){
      t_GVec0 = t_GVec;
    }else{
      //t_GVec0 = t_GVec.elem(g_dup_sample_Index);
      //t_GVec.set_size(t_GVec0.n_elem);
      //t_GVec = t_GVec0;
      t_GVec_cell = g_I_longl_mat * t_GVec;
      t_GVec0 = t_GVec;
      //indexZeroVec_arma = arma::find(t_GVec0 == 0.0);
      //indexNonZeroVec_arma = arma::find(t_GVec0 > 0.0);
      altCounts_new = arma::sum(t_GVec_cell);
      altFreq_new = arma::mean(t_GVec_cell) /2 ;
      //std::cout << "altCounts_new " << altCounts_new << std::endl;
      //
      //
      //get_indexinAnotherVector(indexNonZeroVec, g_dup_sample_Index,indexNonZeroVec_arma);
      //get_indexinAnotherVector(indexZeroVec, g_dup_sample_Index,indexZeroVec_arma);
      //std::cout << "indexNonZeroVec_arma.n_elem " << indexNonZeroVec_arma.n_elem << std::endl;
      //std::cout << "indexZeroVec_arma.n_elem " << indexZeroVec_arma.n_elem << std::endl;
      //gtildeVec.set_size(t_GVec_cell.n_elem);
    }
  */ 

    indexZeroVec.clear();
    indexNonZeroVec.clear();
    t_P2Vec.clear();
    G1tilde_P_G2tilde_Vec.clear();    
    //arma::vec timeoutput5 = getTime(); 
    std::cout << "Here1a" << std::endl;
    //set_varianceRatio(MAC, isSingleVarianceRatio);

    //ptr_gSAIGEobj->set_flagSparseGRM_cur(m_isSparseGRM);

    //if(ptr_gSAIGEobj->m_isnoadjCov){    
    //}	    

  std::cout << "Here2d mainMarkerInCPP" << std::endl;
  //
  //

std::cout << "g_isgxe " << g_isgxe << std::endl;
std::cout << "ptr_gSAIGEobj->m_isCondition " << ptr_gSAIGEobj->m_isCondition << std::endl;

if(g_isgxe && ptr_gSAIGEobj->m_isCondition){
    if(ptr_gSAIGEobj->m_isCondition){
                ptr_gSAIGEobj->m_numMarker_cond = (ptr_gSAIGEobj->m_condition_genoIndex).size();
    }else{
                ptr_gSAIGEobj->m_numMarker_cond = 0;
    }
}
  ptr_gSAIGEobj->assign_for_itrait(i_mt);
    //ptr_gSAIGEobj->assign_for_trait_i(i_mt);
  std::cout << "Here2e mainMarkerInCPP" << std::endl;

    ptr_gSAIGEobj->set_flagSparseGRM_cur(false);		
    if(isSingleVarianceRatio){
  std::cout << "Here2f mainMarkerInCPP" << std::endl;
       ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
       std::cout << "Here2g mainMarkerInCPP" << std::endl;
    }else{	
  std::cout << "Here2gb mainMarkerInCPP" << std::endl;
       hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
  std::cout << "Here2gc mainMarkerInCPP" << std::endl;
    }
    //check 'Main.cpp'
    bool is_region = false;


  std::cout << "Here3 mainMarkerInCPP" << std::endl;

  std::cout << "Here1" << std::endl;
  //std::cout << "ptr_gSAIGEobj->m_flagSparseGRM_cur " << ptr_gSAIGEobj->m_flagSparseGRM_cur << std::endl;

  /*
    if(ptr_gSAIGEobj->m_flagSparseGRM_cur && ptr_gSAIGEobj->m_SigmaMat_sp.n_rows == 2){
	ptr_gSAIGEobj->getadjGFast(t_GVec0, gtildeVec, indexNonZeroVec_arma);    
	arma::fvec tauvec_f = arma::conv_to< arma::fvec >::from(ptr_gSAIGEobj->m_tauvec);
	arma::fvec m_mu2_f = arma::conv_to< arma::fvec >::from(ptr_gSAIGEobj->m_mu2);
	arma::fvec t_GVec_f = arma::conv_to< arma::fvec >::from(gtildeVec);
	arma::fvec sigmainvG_noV_vec = getSigma_G_noV(m_mu2_f, tauvec_f, t_GVec_f, 500, 1e-5, false);
	ptr_gSAIGEobj->m_sigmainvG_noV = arma::conv_to< arma::vec >::from(sigmainvG_noV_vec);
    }
*/
//for(int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
//    int j_mt = i_mt*t_genoIndex.size()+i;



    if(MAC > g_MACCutoffforER){
      std::cout << "Here" << std::endl;
      Unified_getMarkerPval( 
		    t_GVec, 
                          false, // bool t_isOnlyOutputNonZero, 
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT,   
			  altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, true, false);
    }else{
      Unified_getMarkerPval( 
		    t_GVec, 
		    false, // bool t_isOnlyOutputNonZero, 
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT,   
			  altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, false, false);
    }


//    std::cout << "pval " << pval << std::endl;
    //std::cout << "ptr_gSAIGEobj->m_pval_cutoff_for_fastTest " << ptr_gSAIGEobj->m_pval_cutoff_for_fastTest << std::endl;
     
    if(pval < (ptr_gSAIGEobj->m_pval_cutoff_for_fastTest)){
       //ptr_gSAIGEobj->set_flagSparseGRM_cur(true);

       ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
       if(isSingleVarianceRatio){
         ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
       }else{
         hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
       }

       if(ptr_gSAIGEobj->m_flagSparseGRM_cur && ptr_gSAIGEobj->m_SigmaMat_sp.n_rows == 2){
         arma::vec t_GVec0 = g_I_longl_mat * t_GVec;	       
         ptr_gSAIGEobj->getadjGFast(t_GVec0, gtildeVec, indexNonZeroVec_arma);
         arma::fvec tauvec_f = arma::conv_to< arma::fvec >::from(ptr_gSAIGEobj->m_tauvec);
         arma::fvec m_mu2_f = arma::conv_to< arma::fvec >::from(ptr_gSAIGEobj->m_mu2);
         //arma::fvec t_GVec_f = arma::conv_to< arma::fvec >::from(t_GVec);
         arma::fvec t_GVec_f = arma::conv_to< arma::fvec >::from(gtildeVec);
         arma::fvec sigmainvG_noV_vec = getSigma_G_noV(m_mu2_f, tauvec_f, t_GVec_f, 500, 1e-5, false);
         ptr_gSAIGEobj->m_sigmainvG_noV = arma::conv_to< arma::vec >::from(sigmainvG_noV_vec);	
      }



     if(MAC > g_MACCutoffforER || t_traitType.at(i_mt) != "binary"){

      Unified_getMarkerPval(
                    t_GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT,
                          altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, false, ptr_gSAIGEobj->m_flagSparseGRM_cur);



     }else{
      Unified_getMarkerPval(
                    t_GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT,
                          altFreq, isSPAConverge, gtildeVec, is_gtilde, is_region, t_P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, false, false);
 
 
     }     
     }


   if(t_traitType.at(i_mt) == "binary"){
     if(is_Firth){
       mFirth = mFirth + 1;
       if(is_FirthConverge){
		mFirthConverge = mFirthConverge + 1;
       }
     }
   }
   //arma::vec timeoutput6 = getTime();
   //printTime(timeoutput5, timeoutput6, "Unified_getMarkerPval");

   indexNonZeroVec_arma.clear();
   indexZeroVec_arma.clear();
   //std::cout << "isSPAConverge " << isSPAConverge << std::endl;
    BetaVec.at(j_mt) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
    seBetaVec.at(j_mt) = seBeta;       
    pvalVec.at(j_mt) = pval;
    pvalNAVec.at(j_mt) = pval_noSPA;
    TstatVec.at(j_mt) = Tstat * (1 - 2*flip);
    varTVec.at(j_mt) = varT;

    if(isCondition){
    	Beta_cVec.at(j_mt) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true
    	seBeta_cVec.at(j_mt) = seBeta_c;
    	pval_cVec.at(j_mt) = pval_c;
    	pvalNA_cVec.at(j_mt) = pval_noSPA_c;
    	Tstat_cVec.at(j_mt) = Tstat_c * (1 - 2*flip);
    	varT_cVec.at(j_mt) = varT_c;
    }
	
    if(t_traitType.at(i_mt) == "binary"){
       if(t_traitType.size() > 1){
               ptr_gSAIGEobj->assign_for_itrait_binaryindices(i_mt);
	  }
      arma::vec dosage_case = t_GVec.elem(ptr_gSAIGEobj->m_case_indices);
      arma::vec dosage_ctrl = t_GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
      AF_case = arma::mean(dosage_case) /2;
      AF_ctrl = arma::mean(dosage_ctrl) /2;
      N_case = dosage_case.n_elem;
      N_ctrl = dosage_ctrl.n_elem;
      if(flip){
         AF_case = 1-AF_case;
         AF_ctrl = 1-AF_ctrl;
      }
      AF_caseVec.at(j_mt) = AF_case;
      AF_ctrlVec.at(j_mt) = AF_ctrl;

      N_caseVec.at(j_mt) = N_case;
      N_ctrlVec.at(j_mt) = N_ctrl;

      arma::uvec N_case_ctrl_het_hom0;
      if(t_isMoreOutput){	
   	N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5); 
    	N_case_homVec.at(j_mt)  = N_case_ctrl_het_hom0.n_elem;
    	N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
    	N_case_hetVec.at(j_mt) = N_case_ctrl_het_hom0.n_elem;
    	N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
    	N_ctrl_homVec.at(j_mt) = N_case_ctrl_het_hom0.n_elem;
    	N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
    	N_ctrl_hetVec.at(j_mt) = N_case_ctrl_het_hom0.n_elem;
	if(flip){
		N_case_homVec.at(j_mt) = N_case - N_case_hetVec.at(j_mt) -  N_case_homVec.at(j_mt);
		N_ctrl_homVec.at(j_mt) = N_ctrl - N_ctrl_hetVec.at(j_mt) - N_ctrl_homVec.at(j_mt);
	}		
      }	
    }else if(t_traitType.at(i_mt) == "quantitative" || t_traitType.at(i_mt) == "count"){
      N_Vec.at(j_mt) = n;
    }

  

   if(t_traitType.at(i_mt) == "binary" || t_traitType.at(i_mt) == "count"){
      isSPAConvergeVec.at(j_mt) = isSPAConverge;
   }   

//G x E
    //unsigned int ne = g_emat.n_cols;
    std::vector<std::string> Beta_c_ge_vec(ne);
    std::vector<std::string> seBeta_c_ge_vec(ne);
    std::vector<std::string> pval_noSPA_c_ge_vec(ne);
    std::vector<std::string> pval_c_ge_vec(ne);
    arma::vec evec;
    if(g_isgxe && pval < g_pval_cutoff_for_gxe){

  	std::cout << "Here3 G x E" << std::endl;
    	unsigned int qg = 1;
    	unsigned int nrow_e = g_emat.n_rows;
    	arma::mat P2Mat_g(nrow_e, qg);
    	arma::mat VarInvMat_g(qg, qg);
    	arma::mat VarMat_g(qg, qg);
    	arma::vec TstatVec_g(qg);
    	arma::vec pVec_g(qg);
    	arma::vec MAFVec_g(qg);
    	arma::vec w0G2_cond_Vec_g(qg);
    	arma::vec gsumVec(nrow_e, arma::fill::zeros);
    	arma::vec gyVec_g(qg);
    	gyVec_g(0) = gy;
    	double qsum_g = arma::accu(gyVec_g);
    	arma::vec qsumVec_g(qg);
    	qsumVec_g(0) = qsum_g;
    	std::string Beta_c_ge_str;
    	std::string seBeta_c_ge_str;
    	std::string pval_noSPA_c_ge_str;
    	std::string pval_c_ge_str;
    	arma::vec t_GVec1;
    	arma::mat P1Matge(3, nrow_e);
    	arma::mat P2Matge(nrow_e, 3);	
    	arma::mat VarMatge;
    	arma::vec Score_c_ge_vec(ne);
	std::cout << "HEREa0" << std::endl;

	w0G2_cond_Vec_g(0) = 1;
     	TstatVec_g(0) = Tstat;
     	pVec_g(0) = pval;
	MAFVec_g(0) = MAF;
	if((ptr_gSAIGEobj->g_I_longl_mat.n_cols != ptr_gSAIGEobj->g_I_longl_mat.n_rows) && (t_GVec.n_elem < ptr_gSAIGEobj->m_y_gxe_mt.n_rows)){
		t_GVec1 = ptr_gSAIGEobj->g_I_longl_mat * t_GVec;
		indexNonZeroVec_arma = arma::find(t_GVec1 > 0.0);
		indexZeroVec_arma = arma::find(t_GVec1 == 0.0);
	}else{

		t_GVec1 = t_GVec;
	}
	//t_P2Vec.print("t_P2Vec");
	if((ptr_gSAIGEobj->g_I_longl_mat.n_cols != ptr_gSAIGEobj->g_I_longl_mat.n_rows) && (t_GVec.n_elem < ptr_gSAIGEobj->m_y_gxe_mt.n_rows)){
	//if(t_GVec.n_elem < ptr_gSAIGEobj->m_y.n_elem){
	        //t_GVec1 = ptr_gSAIGEobj->g_I_longl_mat * t_GVec;
  		//if(!is_gtilde){
			std::cout << "HEREa0a_0" << std::endl;
                ptr_gSAIGEobj->getadjGFast_gxe(t_GVec1, gtildeVec, indexNonZeroVec_arma);
			std::cout << "HEREa0a_1" << std::endl;
        	//}
	 }else{
		t_GVec1 = t_GVec;
			std::cout << "HEREa0a_2" << std::endl;
  		if(!is_gtilde){
                	ptr_gSAIGEobj->getadjGFast_gxe(t_GVec1, gtildeVec, indexNonZeroVec_arma);
        	}
			std::cout << "HEREa0a_3" << std::endl;

	 }



	std::cout << "HEREa1" << std::endl;

        if(t_P2Vec.n_elem == 0){
                if(!ptr_gSAIGEobj->m_flagSparseGRM_cur){
                        t_P2Vec = gtildeVec % ((ptr_gSAIGEobj->m_mu2_gxe_mt).col(i_mt)) *((ptr_gSAIGEobj->m_tauvec_mt)(0,i_mt));
                }else{
                        t_P2Vec = ptr_gSAIGEobj->getSigma_G_V(gtildeVec, 500, 1e-5);
                }
        }


	std::cout << "HEREa1a" << std::endl;
	P2Mat_g.col(0) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*t_P2Vec;
	std::cout << "HEREa1b" << std::endl;
	VarMat_g = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t() * P2Mat_g;
	std::cout << "HEREa1c" << std::endl;
	VarInvMat_g = VarMat_g.i();	
	std::cout << "HEREa2" << std::endl;

  	ptr_gSAIGEobj->assignConditionFactors(
                                        P2Mat_g,
                                        VarInvMat_g,
                                        VarMat_g,
                                        TstatVec_g,
                                        w0G2_cond_Vec_g,
                                        MAFVec_g,
                                        qsumVec_g,
                                        gtildeVec,
                                        pVec_g);
	ptr_gSAIGEobj->assign_for_itrait(i_mt);				
	std::cout << "HEREa3" << std::endl;
	double Beta_ge, seBeta_ge, pval_ge, Tstat_ge, varT_ge, pval_noSPA_ge, Beta_c_ge, seBeta_c_ge, pval_c_ge, Tstat_c_ge, varT_c_ge, pval_noSPA_c_ge, gy_ge, altFreq_ge;
	arma::vec gtildeVec_ge, t_P2Vec_ge, t_GEVec;
	arma::rowvec G1tilde_P_G2tilde_Vec_ge;
	bool isSPAConverge_ge, is_gtilde_ge, is_region_ge, is_Firth_ge, is_FirthConverge_ge;
	//(ptr_gSAIGEobj->m_varRatio_null_eg).print("ptr_gSAIGEobj->m_varRatio_null_eg");
	unsigned int ksub;
    	for(unsigned int k = 0; k < ne; k++){
      	    ksub = i_mt * ne + k;
	    std::cout << "ksub " << ksub << std::endl;
	    evec = g_emat.col(ksub);
	    std::cout << "ksub b " << ksub << std::endl;
	    t_GEVec = t_GVec1 % evec;
	    altFreq_ge = arma::mean(t_GEVec)/2;
	    is_gtilde_ge = false;
	    gtildeVec_ge.clear();
	    t_P2Vec_ge.clear();
	    G1tilde_P_G2tilde_Vec_ge.clear();

		
	    if(!ptr_gSAIGEobj->m_flagSparseGRM_cur){
	    std::cout << "ksub c " << ksub << std::endl;
	    	ptr_gSAIGEobj->m_varRatioVal = ptr_gSAIGEobj->m_varRatio_null_eg_mt(k,i_mt);		
	//	ptr_gSAIGEobj->m_varRatio_null_eg(k);
	    }else{
	    	ptr_gSAIGEobj->m_varRatioVal =  ptr_gSAIGEobj->m_varRatio_sparse_eg_mt(k,i_mt);
	//	ptr_gSAIGEobj->m_varRatioVal = ptr_gSAIGEobj->m_varRatio_sparse_eg(k);
	    }
	std::cout << "HEREa4" << std::endl;
	std::cout << "ptr_gSAIGEobj->m_varRatioVa " << ptr_gSAIGEobj->m_varRatioVal << std::endl;
	    Unified_getMarkerPval_gxe(
                    t_GEVec,
                          false, // bool t_isOnlyOutputNonZero,
                          indexNonZeroVec_arma, indexZeroVec_arma, Beta_ge, seBeta_ge, pval_ge, pval_noSPA_ge, Tstat_ge, gy_ge, varT_ge, altFreq_ge, isSPAConverge_ge, gtildeVec_ge, is_gtilde_ge, is_region_ge, t_P2Vec_ge, true, Beta_c_ge, seBeta_c_ge, pval_c_ge, pval_noSPA_c_ge, Tstat_c_ge, varT_c_ge, G1tilde_P_G2tilde_Vec_ge, is_Firth_ge, is_FirthConverge_ge, false, true, false);
	std::cout << "HEREa5" << std::endl;

			Beta_c_ge_vec.at(k) = std::to_string(Beta_c_ge);
			seBeta_c_ge_vec.at(k) = std::to_string(seBeta_c_ge);
			pval_c_ge_vec.at(k) = std::to_string(pval_c_ge);
			pval_noSPA_c_ge_vec.at(k) = std::to_string(pval_noSPA_c_ge);
			Score_c_ge_vec(k) = Tstat_c_ge;

	     		P1Matge.row(k) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*(gtildeVec_ge.t());
	std::cout << "HEREa5a" << std::endl;

      if(t_P2Vec_ge.n_elem == 0){
                if(!ptr_gSAIGEobj->m_flagSparseGRM_cur){
                        t_P2Vec_ge = gtildeVec_ge % ((ptr_gSAIGEobj->m_mu2_gxe_mt).col(i_mt)) *((ptr_gSAIGEobj->m_tauvec_mt)(0,i_mt));
                }else{
                        t_P2Vec_ge = ptr_gSAIGEobj->getSigma_G_V(gtildeVec_ge, 500, 1e-5);
                }
        }
	std::cout << "HEREa5b" << std::endl;



	     P2Matge.col(k) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*t_P2Vec_ge;
        }
        std::cout << "HEREa5c" << std::endl;
	VarMatge = P1Matge * P2Matge;
	//arma::vec r_corr = {0.0000, 0.0100, 0.0400, 0.0900, 0.2500, 0.5000, 1.0000};
	arma::vec r_corr = {0.000}; //SKAT test
	//arma::vec r_corr = {1.000};
        std::cout << "HEREa6" << std::endl;
	//Score_c_ge_vec.save("/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/score_example.txt", arma::csv_ascii);
	//VarMatge.save("/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/VarMat_example.txt", arma::csv_ascii);
	Rcpp::List groupOutListge = get_SKAT_pvalue_Rcpp(Score_c_ge_vec, VarMatge, r_corr);
        std::cout << "HEREa7" << std::endl;
	if(ne > 1){     
            Beta_c_ge_str = join(Beta_c_ge_vec, ",");
            seBeta_c_ge_str = join(seBeta_c_ge_vec, ",");
            pval_c_ge_str = join(pval_c_ge_vec, ",");
            pval_noSPA_c_ge_str = join(pval_noSPA_c_ge_vec, ",");
	}else{
	    Beta_c_ge_str = Beta_c_ge_vec.at(0);
	    seBeta_c_ge_str = seBeta_c_ge_vec.at(0);
	    pval_c_ge_str = pval_c_ge_vec.at(0);
	    pval_noSPA_c_ge_str = pval_noSPA_c_ge_vec.at(0);
	}

	Beta_ge_cStrVec.at(j_mt) = Beta_c_ge_str;
	seBeta_ge_cStrVec.at(j_mt) = seBeta_c_ge_str;
	pval_ge_cStrVec.at(j_mt) = pval_c_ge_str;
	pval_noSPA_ge_cStrVec.at(j_mt) = pval_noSPA_c_ge_str;
	//pval_SKATO_ge_cVec.at(i) = Rcpp::as<double>(groupOutListge["Pvalue_SKATO"]);
	pval_SKAT_ge_cVec.at(j_mt) = Rcpp::as<double>(groupOutListge["Pvalue_SKAT"]);


    }//if(g_isgxe){


}//loop through j_mt

   } //    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff || imputeInfo < g_marker_minINFO_cutoff)){
 
  
  
    //t_GVec.clear();
  }

for(unsigned int j_mt = 0; j_mt < t_traitType.size(); j_mt++){

  //output
  writeOutfile_single(t_isMoreOutput,
      t_isImputation,
      isCondition,
      t_isFirth,
  mFirth,
  mFirthConverge,
  t_traitType,
  chrVec,
  posVec,
  markerVec,
  refVec,
  altVec,
  altCountsVec,
  altFreqVec,
  imputationInfoVec,
  missingRateVec,
  BetaVec,
  seBetaVec,
  TstatVec,
  varTVec,
  pvalVec,
  pvalNAVec,
  isSPAConvergeVec,
  Beta_cVec,
  seBeta_cVec,
  Tstat_cVec,
  varT_cVec,
  pval_cVec,
  pvalNA_cVec,
  AF_caseVec,
  AF_ctrlVec,
  N_caseVec,
  N_ctrlVec,
  N_case_homVec,
  N_ctrl_hetVec,
  N_case_hetVec,
  N_ctrl_homVec,
  N_Vec,
  g_isgxe,
Beta_ge_cStrVec,
seBeta_ge_cStrVec,
pval_ge_cStrVec,
pval_noSPA_ge_cStrVec,
pval_SKAT_ge_cVec,
j_mt,
t_genoIndex.size()
);

} //for(unsigned int j_mt = 0; j_mt < t_traitType.size(); j_mt++){



}




// a unified function to get single marker from genotype file
bool Unified_getOneMarker(std::string & t_genoType,   // "PLINK", "BGEN", "Vcf"
                               uint64_t & t_gIndex_prev,        // different meanings for different genoType
                               uint64_t & t_gIndex,        // different meanings for different genoType
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string& t_marker,    // marker ID extracted from genotype file
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool & t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                               std::vector<uint>& t_indexForMissing,     // index of missing genotype data
                               bool & t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                               std::vector<uint>& t_indexForNonZero, //
			       arma::vec & t_GVec,
			       bool t_isImputation 
			       )     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
{
  //arma::vec GVec(ptr_gSAIGEobj->m_n);
  bool isBoolRead = true;
  if(t_genoType == "plink"){
   bool isTrueGenotype = true;
   //t_gIndex_prev is after reading the last marker

   //arma::vec timeoutput1 = getTime();
   ptr_gPLINKobj->getOneMarker(t_gIndex_prev, t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                       t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                       isTrueGenotype, t_GVec);   // t_isTrueGenotype, only used for PLINK format.
   //arma::vec timeoutput2 = getTime();
     //   printTime(timeoutput1, timeoutput2, "Unified_getOneMarker a");
  }
  
  if(t_genoType == "bgen"){
    //bool isBoolRead = true;
    ptr_gBGENobj->getOneMarker(t_gIndex_prev, t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead, t_GVec, t_isImputation);
  }

  if(t_genoType == "vcf"){
    ptr_gVCFobj->getOneMarker(t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero, isBoolRead, t_GVec, t_isImputation);
    ptr_gVCFobj->move_forward_iterator(1);
  }	  
  
  return isBoolRead;
}


// [[Rcpp::export]]
uint32_t Unified_getSampleSizeinGeno(std::string & t_genoType){
    uint32_t N0;
    if(t_genoType == "plink"){
	N0 = ptr_gPLINKobj->getN0();
    }
    if(t_genoType == "bgen"){
	N0 = ptr_gBGENobj->getN0();
    }	    

    if(t_genoType == "vcf"){
       N0 = ptr_gVCFobj->getN0();
    }
    return(N0);
}	

// [[Rcpp::export]]
uint32_t Unified_getSampleSizeinAnalysis(std::string & t_genoType){
    uint32_t N;
    if(t_genoType == "plink"){
        N = ptr_gPLINKobj->getN();
    }
    if(t_genoType == "bgen"){
        N = ptr_gBGENobj->getN();
    }

    if(t_genoType == "vcf"){
       N = ptr_gVCFobj->getN();
    }
    return(N);
}




// a unified function to get marker-level p-value
void Unified_getMarkerPval(
                           arma::vec & t_GVec,
                           bool t_isOnlyOutputNonZero,
			   arma::uvec & t_indexForNonZero_vec,
			   arma::uvec & t_indexForZero_vec,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
			   double& t_pval_noSPA,
                           double& t_Tstat,
			   double& t_gy,
                           double& t_varT,
                           double t_altFreq, 
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
			   arma::rowvec & t_G1tilde_P_G2tilde_Vec, 
			   bool & t_isFirth,
			   bool & t_isFirthConverge, 
			   bool t_isER, 
			   bool t_isnoadjCov,
                           bool t_isSparseGRM)
{
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SAIGE method to calculate marker-level p-values, 't_isOnlyOutputNonZero' should be false.");   


    ptr_gSAIGEobj->getMarkerPval(t_GVec, t_indexForNonZero_vec, t_indexForZero_vec, t_Beta, t_seBeta, t_pval, t_pval_noSPA, t_altFreq, t_Tstat, t_gy, t_varT, t_isSPAConverge, t_gtilde, is_gtilde, is_region, t_P2Vec, t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c, t_pval_noSPA_c, t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde_Vec, t_isFirth, t_isFirthConverge, t_isER, t_isnoadjCov, t_isSparseGRM); //SAIGE_new.cpp
    
    //t_indexForNonZero_vec.clear();
  
}


// a unified function to get marker-level p-value
void Unified_getMarkerPval_gxe(
                           arma::vec & t_GVec,
                           bool t_isOnlyOutputNonZero,
                           arma::uvec & t_indexForNonZero_vec,
                           arma::uvec & t_indexForZero_vec,
                           double& t_Beta,
                           double& t_seBeta,
                           double& t_pval,
                           double& t_pval_noSPA,
                           double& t_Tstat,
                           double& t_gy,
                           double& t_varT,
                           double t_altFreq,
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
                           arma::rowvec & t_G1tilde_P_G2tilde_Vec,
                           bool & t_isFirth,
                           bool & t_isFirthConverge,
                           bool t_isER,
                           bool t_isnoadjCov,
                           bool t_isSparseGRM)
{
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SAIGE method to calculate marker-level p-values, 't_isOnlyOutputNonZero' should be false.");


    ptr_gSAIGEobj->getMarkerPval_gxe(t_GVec, t_indexForNonZero_vec, t_indexForZero_vec, t_Beta, t_seBeta, t_pval, t_pval_noSPA, t_altFreq, t_Tstat, t_gy, t_varT, t_isSPAConverge, t_gtilde, is_gtilde, is_region, t_P2Vec, t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c, t_pval_noSPA_c, t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde_Vec, t_isFirth, t_isFirthConverge, t_isER, t_isnoadjCov, t_isSparseGRM); //SAIGE_new.cpp

    //t_indexForNonZero_vec.clear();

}



//////// ---------- Main functions to set objects for different genotype format --------- ////////////

// [[Rcpp::export]]
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> & t_SampleInModel,
                      std::string t_AlleleOrder)
{
  ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                        t_famFile,
                                        t_bedFile,
                                        //t_SampleInModel,
                                        t_AlleleOrder);
  ptr_gPLINKobj->setPosSampleInPlink(t_SampleInModel);
  
}

// [[Rcpp::export]]
void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> & t_SampleInBgen,
                     std::vector<std::string> & t_SampleInModel,
                     std::string t_AlleleOrder)
{
  std::cout << "t_SampleInBgen " << t_SampleInBgen.size() << std::endl;
  ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                     t_bgenFileIndex,
                                     t_SampleInBgen,
                                     t_SampleInModel,
				     false,
				     false,
                                     t_AlleleOrder);
  //int n = ptr_gBGENobj->getN();
}


// [[Rcpp::export]]
void setVCFobjInCPP(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            std::vector<std::string> & t_SampleInModel)
{
  ptr_gVCFobj = new VCF::VcfClass(t_vcfFileName,
		  		t_vcfFileIndex,
				t_vcfField,
				false,
				t_SampleInModel);

bool isEnd = ptr_gVCFobj->check_iterator_end();

}



//////// ---------- Main functions to set objects for different analysis methods --------- ////////////

// [[Rcpp::export]]
void setSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
	arma::mat & t_Sigma_iXXSigma_iX,
        arma::mat & t_X,
        arma::mat &  t_S_a,
        arma::mat & t_res,
        arma::mat & t_mu2,
        arma::mat & t_mu,
        arma::mat & t_varRatio_sparse,
        arma::mat & t_varRatio_null,
        arma::mat & t_varRatio_null_noXadj,
        arma::mat & t_varRatio_null_eg,
        arma::mat & t_varRatio_sparse_eg,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::mat & t_tauvec,
	arma::mat & t_varWeightsvec,
        std::vector<std::string> & t_traitType,
        arma::mat & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	bool t_isnoadjCov,
	double t_pval_cutoff_for_fastTest,
	bool t_isCondition,
	std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
	double t_pCutoffforFirth,
	arma::mat & t_offset,
	arma::mat & t_resout, 
	arma::sp_mat & t_SigmaMat_sp,
	float t_tauVal_sp,
	arma::sp_mat & t_Ilongmat,
	arma::vec & t_I_longl_vec,
	arma::sp_mat & t_Tlongmat, 
	arma::vec & t_T_longl_vec,
	bool t_is_EmpSPA,
	arma::mat & t_cumul,
	        bool t_is_gxe,
        arma::mat &  t_XV_gxe,
	arma::mat & t_XXVX_inv_gxe,
	arma::mat & t_y_gxe,
	arma::mat & t_res_gxe,
	arma::mat & t_mu2_gxe,
	arma::mat & t_mu_gxe,
	arma::mat & t_varWeights_gxe
	)
{
	//t_SigmaMat_sp.print("t_SigmaMat_sp");
  // check SAIGE.cpp
  ptr_gSAIGEobj = new SAIGE::SAIGEClass(
	t_XVX,
        t_XXVX_inv,
        t_XV,
        t_XVX_inv_XV,
	t_Sigma_iXXSigma_iX,
        t_X,
        t_S_a,
        t_res,
        t_mu2,
        t_mu,
	t_varRatio_sparse,
        t_varRatio_null,
	t_varRatio_null_noXadj,
	t_varRatio_null_eg,
	t_varRatio_sparse_eg,
	t_cateVarRatioMinMACVecExclude,
	t_cateVarRatioMaxMACVecInclude,
        t_SPA_Cutoff,
        t_tauvec,
	t_varWeightsvec,
        t_traitType,
        t_y,
	t_impute_method,
	t_flagSparseGRM,
	t_isnoadjCov,
	t_pval_cutoff_for_fastTest,
	t_isCondition,
	t_condition_genoIndex,
	t_is_Firth_beta,
        t_pCutoffforFirth,
	t_offset, 
	t_resout, 
	t_SigmaMat_sp,
	t_tauVal_sp,
	t_Ilongmat,
        t_I_longl_vec,
        t_Tlongmat,
        t_T_longl_vec,
	t_is_EmpSPA,
	t_cumul,
	t_is_gxe,
        t_XV_gxe,
        t_XXVX_inv_gxe,
        t_y_gxe,
        t_res_gxe,
        t_mu2_gxe,
        t_mu_gxe,
        t_varWeights_gxe);
  //ptr_gSAIGEobj->m_flagSparseGRM = false;
}


/*
// [[Rcpp::export]]
void setSparseSigmaInCPP(int r, arma::umat & t_locationMatinR, arma::vec & t_valueVecinR)
{
  ptr_gSAIGEobj->setupSparseMat(r, t_locationMatinR, t_valueVecinR);
  ptr_gSAIGEobj->m_flagSparseGRM = true;
}
*/

// [[Rcpp::export]]
Rcpp::List RegionSetUpConditional_binary_InCPP(arma::vec & t_weight_cond){

	unsigned int q_cond = (ptr_gSAIGEobj->m_VarInvMat_cond).n_rows;
  	//boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);
  	arma::vec w0G2Vec_cond(q_cond);
  	double w0G2_cond, MAFG2_cond;
        for(unsigned int ci = 0; ci < q_cond; ci++){
		//if(!(t_weight_cond.is_zero())){
			w0G2_cond = t_weight_cond(ci);
		//}else{
                //	MAFG2_cond = (ptr_gSAIGEobj->m_MAF_cond)[ci];
                //	w0G2_cond = boost::math::pdf(beta_dist, MAFG2_cond);
		//}
                	w0G2Vec_cond.at(ci) = w0G2_cond;
        }
	unsigned int nt = ptr_gSAIGEobj->m_traitType_vec.size();
	arma::mat VarMat_weighted_cond_sub;
	arma::mat m_VarMat_weighted_cond(q_cond, q_cond*nt);
	unsigned int startimt, endimt;
	arma::mat woG2Mat_cond = w0G2Vec_cond * (w0G2Vec_cond.t());
	for(unsigned int i_mt = 0; i_mt < nt; i_mt++){
		startimt = i_mt*q_cond;
		endimt = (i_mt+1)*q_cond - 1;
		VarMat_weighted_cond_sub = woG2Mat_cond % ((ptr_gSAIGEobj->m_VarMat_cond).cols(startimt, endimt));
		m_VarMat_weighted_cond.cols(startimt, endimt) = VarMat_weighted_cond_sub;
	}

	 Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("VarMat_G2_cond") = m_VarMat_weighted_cond,
                                          Rcpp::Named("Score_G2_cond") = ptr_gSAIGEobj->m_Tstat_cond,
                                          Rcpp::Named("pval_G2_cond") = ptr_gSAIGEobj->m_p_cond,
                                          Rcpp::Named("gsum_G2_cond") = ptr_gSAIGEobj->m_gsum_cond_Mat,
                                          Rcpp::Named("qsum_G2_cond") = ptr_gSAIGEobj->m_qsum_cond_Vec
					  );
	 return(OutList);

}


//////// ---------- Main function for region-level analysis --------- ////////////
// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
			   arma::mat & annoIndicatorMat,
			   std::vector<std::string> & t_weightlistvec,
			   arma::vec & maxMAFVec, 
			   arma::vec & minMAFVec, 
                           std::string t_outputFile,
			   std::vector<std::string> & t_traitType,
                           unsigned int t_n,           // sample size  
                           arma::mat & P1Mat,            // edited on 2021-08-19: to avoid repeated memory allocation of P1Mat and P2Mat
                           arma::mat & P2Mat, 
			   std::string t_regionTestType, 
			   bool t_isImputation,
                           arma::mat & t_Beta_param,            // edited on 2021-08-19: to avoid repeated memory allocation of P1Mat and P2Mat
			   arma::mat & t_weight,
			   arma::vec & t_weight_cond,
			   bool t_isIncludeNoWeights,
			   bool t_isSingleinGroupTest,
			   bool t_isOutputMarkerList, 
			   std::vector<std::string> & annoStringVec,
			   std::string regionName, 
			   bool t_isFastTest,
			   bool t_isMoreOutput) 
{

  //create the output list
  Rcpp::List OutList = Rcpp::List::create();
  //arma::vec timeoutput1 = getTime();

  //set up weights
  bool isWeightCustomized = false;
  unsigned int q0 = t_genoIndex.size();                 // number of markers (before QC) in one region
 unsigned int q_weight_customize = 0; 
 if(!(t_weight.is_zero()) && t_weight.n_rows == q0){
        isWeightCustomized = true;
        q_weight_customize = t_weight.n_cols;
 }

  unsigned int q_weight_beta = 0;
  bool is_Beta_weight = false;
  if(!t_Beta_param.is_zero()){
     	is_Beta_weight = true;
	q_weight_beta = t_Beta_param.n_rows;
  }

  unsigned int q_anno = annoIndicatorMat.n_cols;
  unsigned int q_maf = maxMAFVec.n_elem;
  unsigned int q_weight = q_weight_customize + q_weight_beta;
  if(t_isIncludeNoWeights){
	q_weight = q_weight + 1;
  }
  arma::vec w0_vec(q_weight);

 
  unsigned int q_anno_maf = q_anno*q_maf;
  unsigned int q_anno_maf_weight = q_anno*q_maf*q_weight;
  arma::mat genoURMat(t_n, q_anno_maf_weight, arma::fill::zeros);
  unsigned int q = q0 + q_anno_maf_weight;
  unsigned int q_multTrait = q * t_traitType.size(); 

  arma::imat annoMAFIndicatorMat(q, q_anno_maf_weight, arma::fill::zeros);
  arma::ivec annoMAFIndicatorVec(q_anno_maf_weight);
  annoMAFIndicatorVec.zeros();
  arma::vec maxMAFperAnno(q_anno, arma::fill::zeros);
  arma::vec minMAFperAnno(q_anno, arma::fill::ones);
  arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
  MAFIndicatorVec.zeros();


  //setting up conditional markers  
  unsigned int q_cond = (ptr_gSAIGEobj->m_VarInvMat_cond).n_rows;	
  arma::rowvec G1tilde_P_G2tilde_Vec(q_cond);
  //boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);

  bool isCondition = ptr_gSAIGEobj->m_isCondition;
  arma::vec w0G2Vec_cond(q_cond);
  double w0G2_cond, MAFG2_cond;
  arma::mat P1Matsub, P2Matsub;


  if(isCondition){
	for(unsigned int ci = 0; ci < q_cond; ci++){
		//if(!(t_weight_cond.is_zero())){
			w0G2_cond = t_weight_cond(ci);
		//}else{	
		//	MAFG2_cond = (ptr_gSAIGEobj->m_MAF_cond)[ci];	
  		//	w0G2_cond = boost::math::pdf(beta_dist, MAFG2_cond);
		//}	
		w0G2Vec_cond.at(ci) = w0G2_cond;
	}
  }

  arma::mat w0G2Mat_cond(q_cond, q_cond);
  w0G2Mat_cond = w0G2Vec_cond * (w0G2Vec_cond.t());


  arma::mat genoSumMat(t_n, q_anno_maf_weight, arma::fill::zeros); //for Phi_cc for binary traits and BURDEN test
  arma::vec genoSumcount_noweight(q_anno_maf_weight, arma::fill::zeros);
  //arma::sp_mat genoSumMat_sp(t_n, q_anno_maf); //for Phi_cc for binary traits and BURDEN test
  std::vector<double> Beta_cVec(q_multTrait, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBeta_cVec(q_multTrait, arma::datum::nan);
  std::vector<double> pval_cVec(q_multTrait, arma::datum::nan);
  std::vector<double> Tstat_cVec(q_multTrait, arma::datum::nan);
  std::vector<double> varT_cVec(q_multTrait, arma::datum::nan);
  std::vector<double> pvalNA_cVec(q_multTrait, arma::datum::nan);
  arma::mat G1tilde_P_G2tilde_Weighted_Mat(q_multTrait, q_cond);
  //group test output
  //arma::vec MAC_GroupVec = arma::zeros<vec>(q_anno_maf);

  arma::vec MAC_GroupVec(q_anno_maf_weight);
  MAC_GroupVec.zeros();
  arma::vec MACCase_GroupVec(q_anno_maf_weight);
  MACCase_GroupVec.zeros();
  arma::vec MACControl_GroupVec(q_anno_maf_weight);
  MACControl_GroupVec.zeros();
  arma::vec NumRare_GroupVec(q_anno_maf_weight);
  NumRare_GroupVec.zeros();
  arma::vec NumUltraRare_GroupVec(q_anno_maf_weight);
  NumUltraRare_GroupVec.zeros();
  arma::vec gtildeVec;
  double MACgroup, MACcasegroup, MACcontrolgroup, AF_case, AF_ctrl, N_case_hom, N_ctrl_het, N_case_het, N_ctrl_hom;
  uint32_t N_case, N_ctrl, N; 

  //single-variant assoc output
  arma::uvec indicatorVec(q, arma::fill::zeros);       // 0: does not pass QC, 1: non-URV, 2: URV
  std::vector<std::string> markerVec(q);
  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altCountsVec(q, arma::datum::nan);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q, arma::datum::nan);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q, arma::datum::nan);
  std::vector<double> altFreqVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MACVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MAFVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_caseVec(q_multTrait, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_ctrlVec(q_multTrait, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<uint32_t> N_caseVec(q_multTrait, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<uint32_t> N_ctrlVec(q_multTrait, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double>  N_case_homVec(q_multTrait, arma::datum::nan);
  std::vector<double>  N_ctrl_hetVec(q_multTrait, arma::datum::nan);
  std::vector<double>  N_case_hetVec(q_multTrait, arma::datum::nan);
  std::vector<double>  N_ctrl_homVec(q_multTrait, arma::datum::nan);
  std::vector<uint32_t> N_Vec(q_multTrait, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> BetaVec(q_multTrait, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q_multTrait, arma::datum::nan);
  std::vector<double> pvalVec(q_multTrait, arma::datum::nan);
  std::vector<double> TstatVec(q_multTrait, arma::datum::nan);
  std::vector<double> TstatVec_flip(q_multTrait, arma::datum::nan);
  std::vector<double> gyVec(q_multTrait, arma::datum::nan);
  std::vector<double> varTVec(q_multTrait, arma::datum::nan);
  std::vector<double> pvalNAVec(q_multTrait, arma::datum::nan);  
  std::vector<bool>  isSPAConvergeVec(q_multTrait);

  unsigned int m1 = g_region_maxMarkers_cutoff;     // number of markers in all chunks expect for the last chunk
  // P2Mat should be of dimension: t_n * m1 
  // Suppose that 
  // n is the sample size in analysis 
  // m (<q) is the number of markers that pass the marker-level QC (e.g., g_missingRate_cutoff and g_region_maxMAF_cutoff)
  // VarMat (m x m) is the variance matrix of these m markers
  // VarMat = P1Mat %*% P2Mat, where P1Mat is of (m x n) and P2Mat is of (n x m)
  // Added on 09-17-2021: we collapse all ultra-rare variants (URV) to get one "fake" marker. 
  // That part has been moved to function mainRegionURVInCPP()
  
  std::vector<unsigned int> mPassCVVec;
    
  // conduct marker-level analysis
  double Beta, seBeta, pval, pval_noSPA, Tstat, varT, gy;
  double Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c;
  bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
  arma::vec P1Vec(t_n), P2Vec(t_n);
  //arma::vec GVec(t_n);
  //arma::vec GZeroVec(t_n);

  int n = ptr_gSAIGEobj->m_n;

  if(g_n_unique > 0){
    n = g_n_unique;
  }

  arma::vec GVec(n);


  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;

  //variance ratio
  bool hasVarRatio = true;;
  bool isSingleVarianceRatio = true;

  if((ptr_gSAIGEobj->m_varRatio_null_mt).n_rows > 1){
    isSingleVarianceRatio = false;
  }else{
    isSingleVarianceRatio = true;
    //ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
  }  
  // initiate chunk information
  unsigned int nchunks = 0; //number of chunks
  unsigned int ichunk = 0; //ith chunk
  unsigned int i1InChunk = 0; //i1th marker in ith chunk
  unsigned int i1 = 0;    // index of Markers (non-URV)
  unsigned int i2 = 0;    // index of Markers (Ultra-Rare Variants, URV)
  unsigned int jm, jmr; 

  double cctpval;
  double cctpval_cond;


  arma::ivec ichunk_vec(t_traitType.size());
  ichunk_vec.zeros();
  arma::ivec i1InChunk_vec(t_traitType.size());
  i1InChunk_vec.zeros();
  arma::ivec nchunks_vec(t_traitType.size());
  nchunks_vec.zeros();
  arma::ivec i1_vec(t_traitType.size());
  i1_vec.zeros();
  arma::ivec i2_vec(t_traitType.size());
  i2_vec.zeros();




  // cycle for q0 markers
  for(unsigned int i = 0; i < q0; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false;

    //GVec.resize(t_n);
    GVec.zeros();

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


    bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec,
					  GVec,
					  t_isImputation);

   //arma::vec timeoutput2a = getTime();
   //printTime(timeoutput1a, timeoutput2a, "Unified_getOneMarker");

   if(!isReadMarker){
      std::cout << "ERROR: Reading " <<  i << "th marker failed." << std::endl;
      break;
    }	    
    std::string pds = std::to_string(pd);
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;

    double MAF = std::min(altFreq, 1 - altFreq);
    double w0;
    //double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
    double MAC = MAF * 2 * n * (1 - missingRate);   // checked on 08-10-2021
    flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);
    MAF = std::min(altFreq, 1 - altFreq);
    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma, indexNonZeroVec_arma0;
    arma::vec t_GVec0;
    //std::cout << "g_n_unique " << g_n_unique << std::endl;
    
    
    indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
    indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);


/*
    if(g_n_unique == 0){
      t_GVec0 = GVec;
      indexNonZeroVec_arma0 = indexNonZeroVec_arma;
    }else{
      //t_GVec0 = GVec.elem(g_dup_sample_Index);
      //t_GVec.set_size(t_GVec0.n_elem);
      //t_GVec = t_GVec0;
      t_GVec0 = g_I_longl_mat * GVec;
      //indexZeroVec_arma0 = arma::find(t_GVec0 == 0.0);
      indexNonZeroVec_arma0 = arma::find(t_GVec0 > 0.0);
      double altCounts_new = arma::sum(t_GVec0);
      gtildeVec.set_size(t_GVec0.n_elem);
      //GVec.resize(t_n);
      //GVec.zeros();
      //GVec = t_GVec0;
      altCounts = arma::accu(t_GVec0);
      altFreq = altCounts/(2*double(t_n));
    }
*/

    //int GVecn = GVec.n_elem;
    //std::cout << "GVecn " << GVecn << std::endl;
    //MAF = std::min(altFreq, 1 - altFreq);
    MAC = std::min(altCounts, t_n *2 - altCounts);
    chrVec.at(i) = chr;
    posVec.at(i) = pds;
    refVec.at(i) = ref;
    altVec.at(i) = alt;
    markerVec.at(i) = marker;             // marker IDs
    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;

    altCountsVec.at(i) = altCounts;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;
    imputationInfoVec.at(i) = imputeInfo;

  for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
    unsigned int j_mt = i_mt*q+i;
    //arma::vec timeoutput3a = getTime();
    //printTime(timeoutput2a, timeoutput3a, "Unified_getOneMarker 2");
   if((missingRate > g_missingRate_cutoff) || (MAF > g_maxMAFLimit) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff) || (imputeInfo < g_marker_minINFO_cutoff)){
      continue;
   }else{
     if(i_mt == 0){	
	 for(unsigned int w = 0; w < q_weight_customize; w++){
		w0_vec(w) = t_weight(i,w); 
	 }
   	for(unsigned int w = 0; w < q_weight_beta; w++){
		boost::math::beta_distribution<> beta_dist(t_Beta_param(w,0), t_Beta_param(w,1));
	        w0_vec(w+q_weight_customize) = boost::math::pdf(beta_dist, MAF);;
	}
	
	if(t_isIncludeNoWeights){
		w0_vec(q_weight_customize+q_weight_beta) = 1;
	}
     } //if(i_mt == 0){

    uint nNonZero = indexNonZeroVec_arma.n_elem;
    ptr_gSAIGEobj->assign_for_itrait(i_mt);

    if(MAC > g_region_minMAC_cutoff){  // not Ultra-Rare Variants

      indicatorVec.at(i) = 1;      
      if(i1InChunk_vec(i_mt) == 0 && i_mt == 0){
        std::cout << "Start analyzing chunk " << ichunk << "....." << std::endl;
      }


      ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
      if(isSingleVarianceRatio){
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur,  ptr_gSAIGEobj->m_isnoadjCov);
      }else{
        hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur,  ptr_gSAIGEobj->m_isnoadjCov);
      }

      if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){ //perform single-variant assoc tests 
 
        //set_varianceRatio(MAC, isSingleVarianceRatio);
        if(MAC > g_MACCutoffforER || t_traitType.at(i_mt) != "binary"){
	  //is region
	  //no ER
	  //
	  //std::cout << "ptr_gSAIGEobj->m_isnoadjCov " << ptr_gSAIGEobj->m_isnoadjCov << std::endl;
          Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, ptr_gSAIGEobj->m_isnoadjCov, ptr_gSAIGEobj->m_flagSparseGRM_cur);



/*
	if(pval < (ptr_gSAIGEobj->m_pval_cutoff_for_fastTest)){
	  ptr_gSAIGEobj->set_flagSparseGRM_cur(true);
	  if(!is_gtilde){
             ptr_gSAIGEobj->getadjGFast(GVec, gtildeVec, indexNonZeroVec_arma);
	  }
          arma::fvec tauvec_f = arma::conv_to< arma::fvec >::from(ptr_gSAIGEobj->m_tauvec);
          arma::fvec m_mu2_f = arma::conv_to< arma::fvec >::from(ptr_gSAIGEobj->m_mu2);
          arma::fvec t_GVec_f = arma::conv_to< arma::fvec >::from(gtildeVec);
	  if(!isSingleVarianceRatio){
            hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
          }else{
            ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
          }	

	  Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, false, true);
	}
*/

	}else{	
          Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
          indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, false, false);
	}

	BetaVec.at(j_mt) = Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true       
        seBetaVec.at(j_mt) = seBeta;       
        pvalVec.at(j_mt) = pval;
        pvalNAVec.at(j_mt) = pval_noSPA;
        TstatVec.at(j_mt) = Tstat * (1 - 2*flip);
        TstatVec_flip.at(j_mt) = Tstat;
        gyVec.at(j_mt) = gy;
        varTVec.at(j_mt) = varT;
        isSPAConvergeVec.at(j_mt) = isSPAConverge;

        if(isCondition){ 	
      	  Beta_cVec.at(j_mt) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
      	  seBeta_cVec.at(j_mt) = seBeta_c;
      	  pval_cVec.at(j_mt) = pval_c;
      	  pvalNA_cVec.at(j_mt) = pval_noSPA_c;
      	  Tstat_cVec.at(j_mt) = Tstat_c * (1 - 2*flip);
      	  varT_cVec.at(j_mt) = varT_c;
	  G1tilde_P_G2tilde_Weighted_Mat.row(j_mt) = G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t() * w0;	
        }

        if(t_regionTestType != "BURDEN"){

	  //std::cout << "P1Mat.n_rows " << P1Mat.n_rows << " P1Mat.n_cols " << P1Mat.n_cols << std::endl;
	  //std::cout << "i_mt*q + i1InChunk " << i_mt*q + i1InChunk << std::endl;
	  //std::cout << "i_mt " << i_mt << std::endl;
	  //std::cout << "q " << q << std::endl;
	  //std::cout << "i1InChunk " << i1InChunk << std::endl;
          //P1Mat.row(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
          P1Mat.row(i_mt*m1 + i1InChunk_vec(i_mt)) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
	  //std::cout << "P2Mat.n_rows " << P2Mat.n_rows << " P2Mat.n_cols " << P2Mat.n_cols << std::endl;
          //P2Mat.col(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
          P2Mat.col(i_mt*m1 + i1InChunk_vec(i_mt)) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
	}
     }//if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){ 

    //if(i_mt == (t_traitType.size()-1)){
    //i1 += 1;
    //}
    //i1InChunk += 1;
    //
    i1_vec(i_mt) = i1_vec(i_mt) + 1;
    i1InChunk_vec(i_mt) = i1InChunk_vec(i_mt) + 1;
     //arma::vec timeoutput3aa = getTime();
     //printTime(timeoutput3a, timeoutput3aa, "Unified_getOneMarker 3a");
     arma::vec dosage_case, dosage_ctrl;
     if(t_traitType.at(i_mt) == "binary"){
       if(t_traitType.size() > 1){
        ptr_gSAIGEobj->assign_for_itrait_binaryindices(i_mt);
       }	
        dosage_case = GVec.elem(ptr_gSAIGEobj->m_case_indices);
        dosage_ctrl = GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
        MACcasegroup = arma::accu(dosage_case);
        MACcontrolgroup = arma::accu(dosage_ctrl);	
    }
      //arma::vec timeoutput3ab = getTime();
      //printTime(timeoutput3aa, timeoutput3ab, "Unified_getOneMarker 3b");
      //arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
   //std::cout << "i_mt first loop a " << i_mt << std::endl; 
   if(i_mt == 0){
      MAFIndicatorVec.zeros();
      MAFIndicatorVec.elem( find(maxMAFVec >= MAF && minMAFVec < MAF) ).ones();	
      annoMAFIndicatorVec.zeros();
    for(unsigned int r = 0; r < q_weight; r++){
      for(unsigned int j = 0; j < q_anno; j++){
        if(annoIndicatorMat(i,j) == 1){
		maxMAFperAnno(j) = std::max(maxMAFperAnno(j), MAF);
		minMAFperAnno(j) = std::min(minMAFperAnno(j), MAF);
		for(unsigned int m = 0; m < q_maf; m++){
			if(MAFIndicatorVec(m) == 1){
			
			   //arma::vec timeoutput3ab0 = getTime();
					jm = j*q_maf + m;	
					jmr = r*q_maf*q_anno + j*q_maf + m;
					annoMAFIndicatorVec(jmr) = 1;
					MAC_GroupVec(jmr) = MAC_GroupVec(jmr) + MAC;
					if(t_traitType.at(i_mt) == "binary"){
						MACCase_GroupVec(jmr) = MACCase_GroupVec(jmr) + MACcasegroup;
						MACControl_GroupVec(jmr) = MACControl_GroupVec(jmr) + MACcontrolgroup;
						//genoSumMat.col(jm) = genoSumMat.col(jm) + w0*GVec;
					}
					w0 = w0_vec(r);
					for(unsigned int k = 0; k < nNonZero; k++){

					//std::cout << "genoSumMat.n_rows " << genoSumMat.n_rows << std::endl;
					//std::cout << "genoSumMat.n_cols " << genoSumMat.n_cols << std::endl;
					//indexNonZeroVec_arma0.print("indexNonZeroVec_arma0");
					//std::cout << "t_GVec0.n_elem " << t_GVec0.n_elem << std::endl;
					genoSumcount_noweight(jmr) = genoSumcount_noweight(jmr) + GVec(indexNonZeroVec_arma(k));

					
					//genoSumMat(indexNonZeroVec_arma(k), jm) = genoSumMat(indexNonZeroVec_arma(k), jm) + w0*t_GVec0(indexNonZeroVec_arma0(k));
					genoSumMat(indexNonZeroVec_arma(k), jmr) = genoSumMat(indexNonZeroVec_arma(k), jmr) + w0*GVec(indexNonZeroVec_arma(k));
					//std::cout << "genoSumcount_noweight.n_elem " << genoSumcount_noweight.n_elem << std::endl;

					//genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + t_GVec0(indexNonZeroVec_arma0(k));
					}
					NumRare_GroupVec(jmr) = NumRare_GroupVec(jmr) + 1;					
			}//if(MAFIndicatorVec(m) == 1){ 

			//std::cout << "jm " << jm << std::endl;
	
  //arma::vec timeoutput3ab2 = getTime();
  //printTime(timeoutput3ab1, timeoutput3ab2, "Unified_getOneMarker 3b2");
			}

		}	
        }
      }
  //arma::vec timeoutput3ac = getTime();
   //    printTime(timeoutput3ab, timeoutput3ac, "Unified_getOneMarker 3c");
     annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();

   }//   if(i_mt == 0){

     if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){
      if(t_traitType.at(i_mt) == "binary"){
        AF_case = arma::mean(dosage_case) /2;
        AF_ctrl = arma::mean(dosage_ctrl) /2;
        if(flip){
          AF_case = 1-AF_case;
           AF_ctrl = 1-AF_ctrl;
        }
        AF_caseVec.at(j_mt) = AF_case;
        AF_ctrlVec.at(j_mt) = AF_ctrl;
	N_case = dosage_case.n_elem;
	N_ctrl = dosage_ctrl.n_elem;
        N_caseVec.at(j_mt) = N_case;
        N_ctrlVec.at(j_mt) = N_ctrl;


        arma::uvec N_case_ctrl_het_hom0;
        if(t_isMoreOutput){
          N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
          N_case_homVec.at(j_mt)  = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
          N_case_hetVec.at(j_mt) = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
          N_ctrl_homVec.at(j_mt) = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
          N_ctrl_hetVec.at(j_mt) = N_case_ctrl_het_hom0.n_elem;
          if(flip){
                N_case_homVec.at(j_mt) = N_case - N_case_hetVec.at(j_mt) -  N_case_homVec.at(j_mt);
                N_ctrl_homVec.at(j_mt) = N_ctrl - N_ctrl_hetVec.at(j_mt) - N_ctrl_homVec.at(j_mt);
          }
        }
      }else if(t_traitType.at(i_mt) == "quantitative" || t_traitType.at(i_mt) == "count"){
        N_Vec.at(j_mt) = t_n;
      }      
     } //if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){

      //arma::vec timeoutput4a = getTime();
      //printTime(timeoutput3ab, timeoutput4a, "Unified_getOneMarker 3");

    }else{   // Ultra-Rare Variants (URV)
     if(i_mt == 0){

      indicatorVec.at(i) = 2;
      arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
      MAFIndicatorVec.zeros();
      MAFIndicatorVec.elem( find(maxMAFVec >= MAF && minMAFVec < MAF) ).ones();
      annoMAFIndicatorVec.zeros();
for(unsigned int r = 0; r < q_weight; r++){
      for(unsigned int j = 0; j < q_anno; j++){
        if(annoIndicatorMat(i,j) == 1){
		maxMAFperAnno(j) = std::max(maxMAFperAnno(j), MAF);
		minMAFperAnno(j) = std::min(minMAFperAnno(j), MAF);
		for(unsigned int m = 0; m < q_maf; m++){
                        if(MAFIndicatorVec(m) == 1){
					jm = j*q_maf + m;	
					jmr = r*q_maf*q_anno + j*q_maf + m;
					annoMAFIndicatorVec(jmr) = 2;
					w0 = w0_vec(r);
					if(q_weight_beta > 0 && r >= q_weight_customize){
						w0 = 1;
					}
				//if(!isWeightCustomized){
				//  for(unsigned int k = 0; k < nNonZero; k++){
				//	genoURMat(indexNonZeroVec_arma(k), jm) = std::max(genoURMat(indexNonZeroVec_arma(k), jm), t_GVec0(indexNonZeroVec_arma(k)));
				//  }
					//genoURMat.col(jm) = arma::max(genoURMat.col(jm), GVec);
				//}else{
                                  	for(unsigned int k = 0; k < nNonZero; k++){

						genoURMat(indexNonZeroVec_arma(k), jmr) = std::max(genoURMat(indexNonZeroVec_arma(k), jmr) , w0 * (GVec(indexNonZeroVec_arma(k))));
					//weightURMat_cnt(indexNonZeroVec_arma(k), jm) = weightURMat_cnt(indexNonZeroVec_arma(k), jm) + 1;
				  	}
				  	NumUltraRare_GroupVec(jmr) = NumUltraRare_GroupVec(jmr) + 1;
				//}
				}
			}
		}
	}
      }
      annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();

      i2 += 1;
     }//if(i_mt == 0){
    }//else{   // Ultra-Rare Variants (URV)


 //}//else if((missingRate > g_missingRate_cutoff) || (MAF > g_maxMAFLimit) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff) || (imputeInfo < g_marker_minINFO_cutoff)){
//
   
    //std::cout << "i_mt first loop " << i_mt << std::endl; 

    if(i1InChunk_vec(i_mt) == m1){
      //std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      if(i_mt == 0){
      std::cout << "In chunks 0-" << ichunk_vec(i_mt) << ", " << i2 << " markers are ultra-rare and " << i1_vec(i_mt) << " markers are not ultra-rare." << std::endl;
	}
	if(t_regionTestType != "BURDEN"){
        P1Matsub = P1Mat.rows(i_mt*m1, i_mt*m1 + i1InChunk_vec(i_mt) -1);
	P2Matsub = P2Mat.cols(i_mt*m1, i_mt*m1 + i1InChunk_vec(i_mt) -1);
        //std::cout << "ichunk " << ichunk << "i_mt " << i_mt << std::endl;
	P1Matsub.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk_vec(i_mt)) + "_trait_"+std::to_string(i_mt)+".bin");
	P2Matsub.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk_vec(i_mt))  + "_trait_"+std::to_string(i_mt) + ".bin");
        //P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        //P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        //P1Mat.rows(i_mt*q_multTrait, i_mt*q_multTrait + i1InChunk -1).save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + "_trait_"+std::to_string(i_mt)+".bin");
        //P2Mat.cols(i_mt*q_multTrait, i_mt*q_multTrait + i1InChunk -1).save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk)  + "_trait_"+std::to_string(i_mt) + ".bin");
      }
      if(i_mt == 0){
      	mPassCVVec.push_back(m1);
      }
      //ichunk += 1;
	ichunk_vec(i_mt) = ichunk_vec(i_mt) + 1;
	//i1InChunk = 0;
	i1InChunk_vec(i_mt) = 0;
      	//nchunks = nchunks + 1;
        nchunks_vec(i_mt) = nchunks_vec(i_mt) + 1;
    }
    Rcpp::checkUserInterrupt();

 }//if((missingRate > g_missingRate_cutoff) || (MAF > g_maxMAFLimit) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff) || (imputeInfo < g_marker_minINFO_cutoff)){, else{
  }//for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
  //    int j_mt = i_mt*t_genoIndex.size()+i;

} //  for(unsigned int i = 0; i < q0; i++)



//the second last chunk
 for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
  if(i1InChunk_vec(i_mt) != 0){
    //std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
    if(t_regionTestType != "BURDEN"){	  
      //P1Mat = P1Mat.rows(0, i1InChunk - 1);
      //P2Mat = P2Mat.cols(0, i1InChunk - 1);
      if(i_mt == 0){	
	std::cout << "In chunks 0-" << ichunk_vec(i_mt) << ", " << i2 << " markers are ultra-rare and " << i1_vec(i_mt) << " markers are not ultra-rare." << std::endl;
	}

	P1Matsub = P1Mat.rows(i_mt*m1, i_mt*m1 + i1InChunk_vec(i_mt) -1);
      //std::cout << "ichunk " << ichunk_vec(i_mt) << "i_mt " << i_mt << std::endl;
      P1Matsub.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk_vec(i_mt)) + "_trait_"+std::to_string(i_mt) + ".bin");
      //P1Mat.rows(i_mt*q_multTrait, i_mt*q_multTrait + i1InChunk -1).save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + "_trait_"+std::to_string(i_mt) + ".bin");
	P2Matsub = P2Mat.cols(i_mt*m1, i_mt*m1 + i1InChunk_vec(i_mt) -1);
	P2Matsub.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk_vec(i_mt)) + "_trait_"+std::to_string(i_mt) + ".bin");
	//P2Mat.cols(i_mt*q_multTrait, i_mt*q_multTrait + i1InChunk -1).save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + "_trait_"+std::to_string(i_mt) + ".bin");
    if(i_mt == 0){
    	mPassCVVec.push_back(i1InChunk_vec(0));
    }
        ichunk_vec(i_mt) = ichunk_vec(i_mt) + 1;
        nchunks_vec(i_mt) = nchunks_vec(i_mt) + 1; 
        i1InChunk_vec(i_mt) = 0;
     //} 

    }else{
    if(i_mt == 0){
      //for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
        std::cout << "In chunks 0-" << ichunk_vec(i_mt) << ", " << i2 << " markers are ultra-rare and " << i1_vec(i_mt) << " markers are not ultra-rare." << std::endl;
    	mPassCVVec.push_back(i1InChunk_vec(0));
    }
	ichunk_vec(i_mt) = ichunk_vec(i_mt) + 1;
	nchunks_vec(i_mt) = nchunks_vec(i_mt) + 1;
	i1InChunk_vec(i_mt) = 0;
      //}	

    }
    //}
  }
}

//for all UR variants
if(i2 > 0){
  int m1new = std::max(m1, q_anno_maf);

  if(t_regionTestType != "BURDEN"){
    P1Mat.resize(m1new * t_traitType.size(), P1Mat.n_cols);
    P2Mat.resize(P2Mat.n_rows,m1new * t_traitType.size());
  }
  arma::mat XV, XXVX_inv;
  ptr_gSAIGEobj->extract_XV_XXVX_inv(XV, XXVX_inv);
  //the last chunk for UR
  unsigned int i_ur, i_i_mt;
 for(unsigned int r = 0; r < q_weight; r++){
  for(unsigned int j = 0; j < q_anno; j++){
   for(unsigned int m = 0; m < q_maf; m++){
   	jm = j*q_maf + m;
	jmr = r*q_anno*q_maf + j*q_maf + m;
	arma::vec genoURVec = genoURMat.col(jmr);
	int n = genoURVec.size();
	arma::uvec indexForNonZero = arma::find(genoURVec != 0); 
	i_ur = q0 + jmr;
	markerVec.at(i_ur) = "UR";             // marker IDs
	if(indexForNonZero.n_elem > 0){
	  double altFreq = arma::mean(genoURVec)/2;
	  double altCounts = arma::accu(genoURVec);
	  double missingRate = 0;
	  double imputeInfo = 1;
    	  std::string chr, ref, alt, marker;
    	  bool flip = false;
	  std::string info = "UR";	
    	  double MAF = std::min(altFreq, 1 - altFreq);
	  double w0 = w0_vec(r);
	  double MAC = MAF * 2 * t_n * (1 - missingRate);
	  std::vector<uint32_t> indexForMissing;
    	  flip = imputeGenoAndFlip(genoURVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);
       	  for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jmr) = genoSumMat(indexForNonZero(k), jmr) + genoURVec(indexForNonZero(k)) * w0;
                genoSumcount_noweight(jmr) = genoSumcount_noweight(jmr) + genoURVec(indexForNonZero(k));
          }

  	 if(t_regionTestType != "BURDEN"){
	    arma::vec genoSumMatvec1 = genoSumMat.col(jmr);
	    arma::vec genoSumMatvec2 = XV * genoSumMatvec1;
	    arma::vec genoSumMatvec3 = genoSumMatvec1 - XXVX_inv * genoSumMatvec2;
	    genoSumMat.col(jmr) = genoSumMatvec3;
          }//if(t_regionTestType != "BURDEN"){

	  MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
          arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
    for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
	 if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){

	if(i_mt == 0){
	    annoMAFIndicatorVec.zeros();
	    annoMAFIndicatorVec(jmr) = 1;
	    annoMAFIndicatorMat.row(i_ur) = annoMAFIndicatorVec.t();
    	    infoVec.at(i_ur) = info;                 // marker information: CHR:POS:REF:ALT
    	    altFreqVec.at(i_ur) = altFreq;	// allele frequencies of ALT allele, this is not always < 0.5.
	    altCountsVec.at(i_ur) = altCounts;
	    missingRateVec.at(i_ur) = missingRate;
    	    MACVec.at(i_ur) = MAC;
    	    MAFVec.at(i_ur) = MAF;
            indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
            indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);
	}


         i_i_mt = q*i_mt + q0 + jmr;
         ptr_gSAIGEobj->assign_for_itrait(i_mt);
         if(!isSingleVarianceRatio){
            hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
         }else{
            ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
         }

	    if(MAC <= g_MACCutoffforER && t_traitType.at(i_mt) == "binary"){	

              ptr_gSAIGEobj->getMarkerPval(genoURVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, false, false);
	    }else{
              ptr_gSAIGEobj->getMarkerPval(genoURVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, true, ptr_gSAIGEobj->m_flagSparseGRM_cur);

	    }

            BetaVec.at(i_i_mt) = Beta* (1 - 2*flip);
            seBetaVec.at(i_i_mt) = seBeta;
            pvalVec.at(i_i_mt) = pval;
            pvalNAVec.at(i_i_mt) = pval_noSPA;
            TstatVec.at(i_i_mt) = Tstat * (1 - 2*flip);
            TstatVec_flip.at(i_i_mt) = Tstat;
            gyVec.at(i_i_mt) = gy;
            varTVec.at(i_i_mt) = varT;
            isSPAConvergeVec.at(i_i_mt) = isSPAConverge;
            if(isCondition){
              Beta_cVec.at(i_i_mt) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
              seBeta_cVec.at(i_i_mt) = seBeta_c;
              pval_cVec.at(i_i_mt) = pval_c;
              pvalNA_cVec.at(i_i_mt) = pval_noSPA_c;
              Tstat_cVec.at(i_i_mt) = Tstat_c * (1 - 2*flip);
              varT_cVec.at(i_i_mt) = varT_c;
              G1tilde_P_G2tilde_Weighted_Mat.row(i_i_mt) = G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t() * w0;
	      //check cond
            }
        arma::vec dosage_case, dosage_ctrl;

	if(i_mt == 0){
    	    chrVec.at(i_ur) = "UR";
    	    posVec.at(i_ur) = "UR";
    	    refVec.at(i_ur) = "UR";
    	    altVec.at(i_ur) = "UR";

	    std::string str = std::to_string(maxMAFVec.at(m));
	    str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
	    std::string strmin = std::to_string(minMAFVec.at(m));
	    strmin.erase ( strmin.find_last_not_of('0') + 1, std::string::npos );
    	    markerVec.at(i_ur) = regionName + ":" + annoStringVec.at(j) + ":" + strmin+":"+str ;
            MAC_GroupVec(jmr) = MAC_GroupVec(jmr) + MAC;
	} //if(i_mt == 0){
            if(t_traitType.at(i_mt) == "binary"){
	           if(t_traitType.size() > 1){
		              ptr_gSAIGEobj->assign_for_itrait_binaryindices(i_mt);
		    }
                        dosage_case = genoURVec.elem(ptr_gSAIGEobj->m_case_indices);
                        dosage_ctrl = genoURVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                        MACcasegroup = arma::accu(dosage_case);
                        MACcontrolgroup = arma::accu(dosage_ctrl);
                        MACCase_GroupVec(jmr) = MACCase_GroupVec(jmr) + MACcasegroup;
                        MACControl_GroupVec(jmr) = MACControl_GroupVec(jmr) + MACcontrolgroup;
            //}
          //if(t_traitType.at(i_mt) == "binary"){
            AF_case = arma::mean(dosage_case) /2;
            AF_ctrl = arma::mean(dosage_ctrl) /2;
            if(flip){
              AF_case = 1-AF_case;
              AF_ctrl = 1-AF_ctrl;
            }
            AF_caseVec.at(i_i_mt) = AF_case;
            AF_ctrlVec.at(i_i_mt) = AF_ctrl;
            N_caseVec.at(i_i_mt) = dosage_case.n_elem;
            N_ctrlVec.at(i_i_mt) = dosage_ctrl.n_elem;
            arma::uvec N_case_ctrl_het_hom0;
            if(t_isMoreOutput){
          	N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
          	N_case_homVec.at(i_i_mt)  = N_case_ctrl_het_hom0.n_elem;
          	N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
          	N_case_hetVec.at(i_i_mt) = N_case_ctrl_het_hom0.n_elem;
          	N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
          	N_ctrl_homVec.at(i_i_mt) = N_case_ctrl_het_hom0.n_elem;
          	N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
          	N_ctrl_hetVec.at(i_i_mt) = N_case_ctrl_het_hom0.n_elem;
          	if(flip){
                	N_case_homVec.at(i_i_mt) = N_case - N_case_hetVec.at(i_i_mt) -  N_case_homVec.at(i_i_mt);
                	N_ctrl_homVec.at(i_i_mt) = N_ctrl - N_ctrl_hetVec.at(i_i_mt) - N_ctrl_homVec.at(i_i_mt);
          	}
            }//if(t_isMoreOutput){

        }else if(t_traitType.at(i_mt) == "quantitative" || t_traitType.at(i_mt) == "count"){
          	N_Vec.at(i_i_mt) = n;
        }

      if(t_regionTestType != "BURDEN"){	

	  
	  P1Mat.row(i_mt*m1new + i1InChunk_vec(i_mt)) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
          //P2Mat.col(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
          P2Mat.col(i_mt*m1new + i1InChunk_vec(i_mt)) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;


        //P1Mat.row(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
	//gtildeVec.print("gtildeVec");
        //P2Mat.col(i1InChunk) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
	//P2Vec.print("P2Vec");

      } //if(t_regionTestType != "BURDEN"){


    }else{//if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){	
	if(i_mt == 0){
            MAC_GroupVec(jmr) = MAC_GroupVec(jmr) + MAC;
	}

	    arma::vec dosage_case, dosage_ctrl;
            if(t_traitType.at(i_mt) == "binary"){
		if(t_traitType.size() > 1){
	     		ptr_gSAIGEobj->assign_for_itrait_binaryindices(i_mt);
		}

                        dosage_case = genoURVec.elem(ptr_gSAIGEobj->m_case_indices);
                        dosage_ctrl = genoURVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                        MACcasegroup = arma::accu(dosage_case);
                        MACcontrolgroup = arma::accu(dosage_ctrl);
                        MACCase_GroupVec(jmr) = MACCase_GroupVec(jmr) + MACcasegroup;
                        MACControl_GroupVec(jmr) = MACControl_GroupVec(jmr) + MACcontrolgroup;
            }
    }
  //}
   	//if(i_mt == (t_traitType.size()-1)){
    		i1InChunk_vec(i_mt) = i1InChunk_vec(i_mt) + 1;
    		i1_vec(i_mt) = i1_vec(i_mt) + 1;
   	//}

    }//    }//for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){

    }
   }  
  }
 }

  if(i1InChunk_vec(0) != 0){
    //nchunks = nchunks + 1;
    if(t_regionTestType != "BURDEN"){
      //P1Mat = P1Mat.rows(0, i1InChunk - 1);
      //P2Mat = P2Mat.cols(0, i1InChunk - 1);
//if(nchunks != 1){
      //P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      //P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");

     for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
      P1Matsub = P1Mat.rows(i_mt*m1new, i_mt*m1new + i1InChunk_vec(i_mt) -1);
      P2Matsub = P2Mat.cols(i_mt*m1new, i_mt*m1new + i1InChunk_vec(i_mt) -1);
      //std::cout << "ichunk " << ichunk << "i_mt " << i_mt << std::endl;
      P1Matsub.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk_vec(i_mt)) + "_trait_"+std::to_string(i_mt) + ".bin");
      P2Matsub.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk_vec(i_mt)) + "_trait_"+std::to_string(i_mt) + ".bin");
      //P1Mat.rows(i_mt*m1new, i_mt*m1new + i1InChunk -1).save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + "_trait_"+std::to_string(i_mt) + ".bin");
      //P2Mat.col(i_mt*m1new, i_mt*m1new + i1InChunk -1).save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + "_trait_"+std::to_string(i_mt) + ".bin");
	nchunks_vec(i_mt) = nchunks_vec(i_mt) + 1;
	ichunk_vec(i_mt) = ichunk_vec(i_mt) + 1;
	}
    }else{
	for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
		nchunks_vec(i_mt) = nchunks_vec(i_mt) + 1;
		ichunk_vec(i_mt) = ichunk_vec(i_mt) + 1;	
	}	
    }
      //ichunk = ichunk + 1;
//    }
    mPassCVVec.push_back(i1InChunk_vec(0));
  }
   //std::cout << "P1Mat.n_rows ok2 " << P1Mat.n_rows << std::endl; 

  }// if(i2 > 0)    
  int mPassCVVecsize = mPassCVVec.size();
  for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
   nchunks_vec(i_mt) = mPassCVVecsize;
  }

arma::mat VarMat, VarMatsub;
//(i1, i1);

i1=i1_vec(0);

if(t_regionTestType != "BURDEN"){
  VarMat.resize(i1, i1*t_traitType.size());
  VarMatsub.resize(i1, i1);
  for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){

  //VarMatsub.clear();
  if(nchunks_vec(i_mt) == 1){
    arma::vec VarMat_sub = P1Mat.rows(i_mt*i1, (i_mt+1)*i1-1) * P2Mat.cols(i_mt*i1, (i_mt+1)*i1-1);
    VarMat.submat(0, i_mt*i1, i1-1, ((i_mt+1)*i1-1)) = VarMat_sub; 
  }

  // the region includes more markers than limitation, so P1Mat and P2Mat have been put in hard drive
  //
  //
  //std::cout << "nchunks_vec(i_mt) " << nchunks_vec(i_mt) << std::endl;
  //nchunks_vec.print("nchunks_vec");
  
  //for(unsigned int mp = 0; mp <   mPassCVVec.size(); mp++){
  //	 std::cout << "mp " << mp << "   mPassCVVec.at(mp) " << mPassCVVec.at(mp) << std::endl; 
  //}
  

  if(nchunks_vec(i_mt) > 1)
  {
    int first_row = 0, first_col = 0, last_row = 0, last_col = 0;
    
    for(unsigned int index1 = 0; index1 < nchunks_vec(i_mt); index1++)
    {
      last_row = first_row + mPassCVVec.at(index1) - 1;
      
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + "_trait_"+std::to_string(i_mt) + ".bin";
      
      P1Mat.load(P1MatFile);
   
       //P1Mat.print("P1Mat");

      if(P1Mat.n_cols == 0) continue;
      
      // off-diagonal sub-matrix
      for(unsigned int index2 = 0; index2 < index1; index2++)
      {
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks_vec(i_mt) - 1 << ", " << index2 << "/" << nchunks_vec(i_mt) - 1 << ")........" << std::endl;
       
        P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index2) + "_trait_"+std::to_string(i_mt) + ".bin");
       //P2Mat.print("P2Mat");
        
        if(P2Mat.n_cols == 0) continue;
        arma::mat offVarMat = P1Mat * P2Mat;
        last_col = first_col + mPassCVVec.at(index2) - 1;
        
        VarMatsub.submat(first_row, first_col, last_row, last_col) = offVarMat;
        VarMatsub.submat(first_col, first_row, last_col, last_row) = offVarMat.t();
        first_col = last_col + 1;
      }
      
      // diagonal sub-matrix
      last_col = first_col + mPassCVVec.at(index1) - 1;
      std::cout << "Analyzing chunks (" << index1 << "/" << nchunks_vec(i_mt) - 1 << ", " << index1 << "/" << nchunks_vec(i_mt) - 1 << ")........" << std::endl;
      //std::cout << "P2Mat.n_cols " << P2Mat.n_cols << std::endl;
      P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1)  + "_trait_"+std::to_string(i_mt) + ".bin");
	//P2Mat.print("P2Mat");
	//P1Mat.print("P1Mat");

      //std::cout << "P1Mat.n_rows " << P1Mat.n_rows << " P1Mat.n_cols " << P1Mat.n_cols << std::endl;
      //std::cout << "P2Mat.n_rows " << P2Mat.n_rows << " P2Mat.n_cols " << P2Mat.n_cols << std::endl;
	
      
      arma::mat diagVarMat = P1Mat * P2Mat;
    
      VarMatsub.submat(first_row, first_col, last_row, last_col) = diagVarMat;
      //diagVarMat.print("diagVarMat") ;
      
      first_row = last_row + 1;
      first_col = 0;
      Rcpp::checkUserInterrupt();
    }
     
    for(unsigned int index1 = 0; index1 < nchunks_vec(i_mt); index1++)
    {
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + "_trait_"+std::to_string(i_mt) + ".bin";
      std::string P2MatFile = t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + "_trait_"+std::to_string(i_mt) + ".bin";
      const char* File1 = P1MatFile.c_str();
      const char* File2 = P2MatFile.c_str();
      std::remove(File1);
      std::remove(File2);
    }
    
    VarMat.submat(0, i_mt*i1, i1-1, ((i_mt+1)*i1-1)) = VarMatsub;	
 
	//VarMatsub.print("VarMatsub");


    }//if(nchunks > 1)

  } //    for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){

}  //if(t_regionTestType != "BURDEN"){
//}
//std::cout << "VarMat.n_rows " << VarMat.n_rows << " " << VarMat.n_cols << std::endl;
//VarMat.print("VarMat");
//read and test single markers done
//arma::vec timeoutput2 = getTime();
//printTime(timeoutput1, timeoutput2, "read and test single markers done");
//check the max MAF of all markers
//if there are sets containing the same markers
//arma::uvec q_maf_for_anno(q_anno);
arma::umat q_maf_for_anno(q_anno, q_maf);
q_maf_for_anno.zeros();
arma::uvec uniqmax_ind;
arma::uvec uniqmin_ind;
arma::uvec uniq_ind, uniq_ct;
arma::vec minMAFsub, maxMAFsub, uniqmax, uniqmin;
arma::uvec minMAFsubind, maxMAFsubind;
for(unsigned int j = 0; j < q_anno; j++){
	//maxMAFperAnno.print("maxMAFperAnno");
	//minMAFperAnno.print("minMAFperAnno");
	//maxMAFVec.print("maxMAFVec");
	//minMAFVec.print("minMAFVec");
	arma::uvec jtemp = find(maxMAFVec >= maxMAFperAnno(j) && minMAFVec < minMAFperAnno(j));
	if(jtemp.n_elem > 1){
		for(unsigned int jt = 1; jt < jtemp.n_elem; jt++){
			q_maf_for_anno(j,jtemp(jt)) = 1;
			//q_maf_for_anno(j) = jtemp.min();
			//q_maf_for_anno.rows() = jtemp.min();
		}
	}

	uniqmax = arma::unique(maxMAFVec);

	for(unsigned int id = 0; id < uniqmax.n_elem; id++){
		uniq_ct = arma::find(maxMAFVec == uniqmax(id));
		//uniq_ct = arma::find(maxMAFVec == maxMAFVec(uniqmax_ind(id)));
		if(uniq_ct.n_elem > 1){
			minMAFsub = minMAFVec(uniq_ct);
			if(minMAFsub(0) <  minMAFperAnno(j)){
				for(unsigned int im = 1; im < uniq_ct.n_elem; im++){
					q_maf_for_anno(j,uniq_ct(im)) = 1;
				}
			}
		}
	}

	uniqmin = arma::unique(minMAFVec);
	for(unsigned int id = 0; id < uniqmin_ind.n_elem; id++){
		uniq_ct = arma::find(minMAFVec == uniqmin(id));
		if(uniq_ct.n_elem > 1){
			maxMAFsub = maxMAFVec(uniq_ct);
			if(maxMAFsub(0) >  maxMAFperAnno(j)){
				for(unsigned int im = 1; im < uniq_ct.n_elem; im++){
					q_maf_for_anno(j,uniq_ct(im)) = 1;
				}
			}
		}	
	}
}


//If only conduct Burden test
arma::vec BURDEN_pval_Vec(q_anno_maf_weight * (t_traitType.size()));
BURDEN_pval_Vec.fill(-1.0);
arma::vec BURDEN_pval_cVec(q_anno_maf_weight * (t_traitType.size()));
BURDEN_pval_cVec.fill(-1.0);
std::vector<std::string> BURDEN_AnnoName_Vec(q_anno_maf_weight * (t_traitType.size()));
std::vector<std::string> BURDEN_maxMAFName_Vec(q_anno_maf_weight * (t_traitType.size()));
std::vector<double> BURDEN_Beta_Vec(q_anno_maf_weight * (t_traitType.size()));
std::vector<double> BURDEN_seBeta_Vec(q_anno_maf_weight * (t_traitType.size()));
//std::vector<double> BURDEN_pval_cVec(q_anno_maf);
std::vector<double> BURDEN_Beta_cVec(q_anno_maf_weight * (t_traitType.size()));
std::vector<double> BURDEN_seBeta_cVec(q_anno_maf_weight * (t_traitType.size()));



bool iswriteOutput = false;
bool isregion = true;
if(!ptr_gSAIGEobj->m_flagSparseGRM){
	isregion = false;
}
//Rcpp::DataFrame OUT_BURDEN = Rcpp::DataFrame::create();
unsigned int i_b= 0;
unsigned int q_maf_m;
bool isPolyMarker = true;
std::string AnnoName;
double maxMAFName;
unsigned int startt, endt;
//std::cout << "If only conduct Burden test b " << std::endl;

if(t_regionTestType == "BURDEN"){
for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
  ptr_gSAIGEobj->assign_for_itrait(i_mt);
  startt = i_mt*q_anno*q_maf*q_weight;
  endt = (i_mt+1)*q_anno*q_maf*q_weight - 1;

     for(unsigned int r = 0; r < q_weight; r++){
     for(unsigned int j = 0; j < q_anno; j++){
       AnnoName = annoStringVec[j];
       isPolyMarker = true;	
       for(unsigned int m = 0; m < q_maf; m++){
        q_maf_m = q_maf_for_anno(j,m);
	maxMAFName = maxMAFVec(m); 
	jm = j*q_maf+m;
	//i = jm;
       //if(m <= q_maf_m){
       if(q_maf_m != 1){

	jmr = r*q_anno*q_maf + j*q_maf + m;
	i_b = i_mt*q_anno*q_maf*q_weight + jmr;

        arma::vec genoSumVec = genoSumMat.col(jmr);
        int n = genoSumVec.size();
        arma::uvec indexNonZeroVec_arma = arma::find(genoSumVec != 0);
	arma::uvec indexZeroVec_arma = arma::find(genoSumVec == 0);
	double altCounts =  genoSumcount_noweight(jmr);
	double altFreq = altCounts/(2*t_n);
	double MAC, MAF;
	if(altFreq > 1){
		MAF = 1;
		MAC = t_n;
	}else{
		MAF = std::min(altFreq, 1 - altFreq);
		MAC = MAF*2*t_n;
	}
	//std::cout << "altCounts " << altCounts << std::endl;
	//std::cout << "altFreq " << altCounts << std::endl;
        double missingRate = 0;
        double imputeInfo = 1;
        std::string chr, ref, alt, marker;
        bool flip = false;
        std::string info = "UR";
        //double MAF = std::min(altFreq, 1 - altFreq);
        double w0;
        //double MAC = MAF * 2 * t_n * (1 - missingRate);
        if(indexNonZeroVec_arma.n_elem > 0 && MAC >= g_min_gourpmac_for_burdenonly){
	  //if(MAC >= g_min_gourpmac_for_burdenonly){
          isPolyMarker = true;   
	  std::vector<uint32_t> indexForMissing;
			
    	if(isSingleVarianceRatio){
       		ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
    	}else{
       		hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
    	}	
	  //arma::vec timeoutput_getp = getTime();
	  if(MAC <= g_MACCutoffforER && t_traitType.at(i_mt) == "binary"){
          ptr_gSAIGEobj->getMarkerPval(genoSumVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, isregion, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, false, ptr_gSAIGEobj->m_flagSparseGRM_cur);
	  }else{
          ptr_gSAIGEobj->getMarkerPval(genoSumVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, isregion, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, false, ptr_gSAIGEobj->m_flagSparseGRM_cur);
	  }	  
	  //arma::vec timeoutput_getp2 = getTime();
	  //printTime(timeoutput_getp, timeoutput_getp2, "get p  done");
	  
	  if(isCondition){
	    BURDEN_pval_cVec(i_b) = pval_c;
	    BURDEN_Beta_cVec.at(i_b) = Beta_c;
	    BURDEN_seBeta_cVec.at(i_b) = seBeta_c;
          }

	   BURDEN_AnnoName_Vec.at(i_b) = AnnoName;
	   std::string str = std::to_string(maxMAFName);
	   str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
	   BURDEN_maxMAFName_Vec.at(i_b) = str;
	   BURDEN_pval_Vec(i_b) = pval;
	   BURDEN_Beta_Vec.at(i_b) = Beta;
	   BURDEN_seBeta_Vec.at(i_b) = seBeta;

        }else{ //if(indexNonZeroVec_arma.n_elem > 0 && MAC >= g_min_gourpmac_for_burdenonly){
	   isPolyMarker = false;	 
	}
     }else{
	if(isPolyMarker){
	  if(isCondition){
            BURDEN_pval_cVec(i_b) = pval_c;
            BURDEN_Beta_cVec.at(i_b) = Beta_c;
            BURDEN_seBeta_cVec.at(i_b) = seBeta_c;
          }

           BURDEN_AnnoName_Vec.at(i_b) = AnnoName;
           BURDEN_maxMAFName_Vec.at(i_b) = std::to_string(maxMAFName);
           BURDEN_pval_Vec(i_b) = pval;
           BURDEN_Beta_Vec.at(i_b) = Beta;
           BURDEN_seBeta_Vec.at(i_b) = seBeta; 
	 }
     }//else{  //if(q_maf_m != 1){	     
    

      }//for(unsigned int r = 0; r < q_weight; r++){	
     }
   }

//std::cout << "If only conduct Burden test c " << std::endl;
//}//for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
//}
	arma::vec BURDEN_pval_Vec_onetrait = BURDEN_pval_Vec.subvec(startt, endt);
	   arma::uvec nonMissingPvalVecInd = arma::find(BURDEN_pval_Vec_onetrait >= 0);
	   arma::vec nonMissingPvalVec = BURDEN_pval_Vec_onetrait.elem(nonMissingPvalVecInd);

	   cctpval = CCT_cpp(nonMissingPvalVec);
           if(isCondition){
	   arma::vec BURDEN_pval_cVec_onetrait = BURDEN_pval_cVec.subvec(startt, endt);
	   arma::vec nonMissingPvalVec_cond = BURDEN_pval_cVec_onetrait.elem(nonMissingPvalVecInd);
	   cctpval_cond = CCT_cpp(nonMissingPvalVec_cond);	   
           }

	   if(!t_isFastTest){
		iswriteOutput = true;
	   }else{
		if(ptr_gSAIGEobj->m_flagSparseGRM_cur){
			iswriteOutput = true;
		}else{
		  if(cctpval >= 0.1){
			iswriteOutput = true;
		  } 	
		}
	   }
	
	   if(iswriteOutput){

//for(unsigned int j_mt = 0; j_mt < t_traitType.size(); j_mt++){

		writeOutfile_BURDEN(regionName,
			BURDEN_AnnoName_Vec,
			BURDEN_maxMAFName_Vec,
			BURDEN_pval_Vec,
			BURDEN_Beta_Vec,
			BURDEN_seBeta_Vec,
			BURDEN_pval_cVec,
			BURDEN_Beta_cVec,
			BURDEN_seBeta_cVec,
			MAC_GroupVec,
			MACCase_GroupVec,
			MACControl_GroupVec,
			NumRare_GroupVec,
			NumUltraRare_GroupVec,
			cctpval,
			cctpval_cond,
			q_weight,
			q_anno,
			q_maf,
			isCondition,
			i_mt);
	   }else{
	     OutList.push_back(iswriteOutput, "iswriteOutput");
	   }

 //arma::vec timeoutput3 = getTime();
 //printTime(timeoutput2, timeoutput3, "burden test done");
 //
 }//for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
 }else{//if(t_regionTestType == "BURDEN"){
  //q_maf_for_anno = q_maf_for_anno + 1;
  OutList.push_back(MAC_GroupVec, "MAC_GroupVec");
  //OutList.push_back(q_maf_for_anno, "q_maf_for_annoVec");
  OutList.push_back(q_maf_for_anno, "q_maf_for_annoMat");
  if(t_traitType.at(0) == "binary"){
    OutList.push_back(MACCase_GroupVec, "MACCase_GroupVec");
    OutList.push_back(MACControl_GroupVec, "MACCtrl_GroupVec");
    OutList.push_back(genoSumMat, "genoSumMat");
    OutList.push_back(gyVec, "gyVec");
  }

    OutList.push_back(VarMat, "VarMat");	
    OutList.push_back(MAFVec, "MAFVec");	
    OutList.push_back(TstatVec_flip, "TstatVec_flip");	
  //arma::mat scaled_m_VarInvMat_cond;
    if(isCondition){
  //std::cout << "okk5" << std::endl;
   std::cout << "G1tilde_P_G2tilde_Weighted_Mat.n_rows " << G1tilde_P_G2tilde_Weighted_Mat.n_rows << std::endl;	
   std::cout << "G1tilde_P_G2tilde_Weighted_Mat.n_cols " << G1tilde_P_G2tilde_Weighted_Mat.n_cols << std::endl;	
   std::cout << "(ptr_gSAIGEobj->m_VarInvMat_cond).n_rows " << (ptr_gSAIGEobj->m_VarInvMat_cond).n_rows  << std::endl;
   std::cout << "(ptr_gSAIGEobj->m_VarInvMat_cond).n_cols " << (ptr_gSAIGEobj->m_VarInvMat_cond).n_cols  << std::endl;

   arma::mat AdjCondMat(q_multTrait, q_cond);
   arma::mat VarMatAdjCond(q_multTrait, q);
   arma::vec TstatAdjCond(q_multTrait);
   arma::mat AdjCondMat_sub(q, q_cond);
   arma::mat VarInvMat_cond_sub(q_cond, q_cond);
   arma::mat VarMatAdjCond_sub(q,q);
   arma::vec TstatAdjCond_sub(q);
   arma::vec Tstat_cond_sub(q_cond);
   arma::mat G1tilde_P_G2tilde_Weighted_Mat_sub(q, q_cond);
unsigned int startt_qcond, endt_qcond;
   for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
  	startt = i_mt*q;
  	endt = (i_mt+1)*q - 1;
	startt_qcond = i_mt*q_cond;
	endt_qcond = (i_mt+1)*q_cond - 1;
	G1tilde_P_G2tilde_Weighted_Mat_sub = G1tilde_P_G2tilde_Weighted_Mat.rows(startt,endt);
	VarInvMat_cond_sub = (ptr_gSAIGEobj->m_VarInvMat_cond).cols(startt_qcond, endt_qcond);
	AdjCondMat_sub = G1tilde_P_G2tilde_Weighted_Mat_sub * (VarInvMat_cond_sub / (w0G2Mat_cond));
	//w0G2Mat_cond.print("w0G2Mat_cond");
	//(ptr_gSAIGEobj->m_VarInvMat_cond).print("ptr_gSAIGEobj->m_VarInvMat_cond");
	//AdjCondMat_sub.print("AdjCondMat_sub");
	//G1tilde_P_G2tilde_Weighted_Mat_sub.print("G1tilde_P_G2tilde_Weighted_Mat_sub");
	AdjCondMat.rows(startt,endt) = AdjCondMat_sub;
   std::cout << "G1tilde_P_G2tilde_Weighted_Mat.n_rows " << G1tilde_P_G2tilde_Weighted_Mat.n_rows << std::endl;	
	VarMatAdjCond_sub = AdjCondMat_sub * (G1tilde_P_G2tilde_Weighted_Mat_sub.t());
	//VarMatAdjCond_sub.print("VarMatAdjCond_sub");
	VarMatAdjCond.rows(startt,endt) = VarMatAdjCond_sub;
   std::cout << "G1tilde_P_G2tilde_Weighted_Mat.n_rows " << G1tilde_P_G2tilde_Weighted_Mat.n_rows << std::endl;
   	Tstat_cond_sub = (ptr_gSAIGEobj->m_Tstat_cond).subvec(startt_qcond,endt_qcond);
	TstatAdjCond_sub = AdjCondMat_sub * (Tstat_cond_sub % w0G2Vec_cond );
	TstatAdjCond.subvec(startt, endt) = TstatAdjCond_sub;
   std::cout << "G1tilde_P_G2tilde_Weighted_Mat.n_rows " << G1tilde_P_G2tilde_Weighted_Mat.n_rows << std::endl;	
   }
      //arma::mat AdjCondMat = G1tilde_P_G2tilde_Weighted_Mat * (ptr_gSAIGEobj->m_VarInvMat_cond / (w0G2Mat_cond));
      //arma::mat VarMatAdjCond = AdjCondMat * (G1tilde_P_G2tilde_Weighted_Mat.t());
      //arma::vec TstatAdjCond = AdjCondMat * (ptr_gSAIGEobj->m_Tstat_cond % w0G2Vec_cond ); 
      OutList.push_back(G1tilde_P_G2tilde_Weighted_Mat, "G1tilde_P_G2tilde_Weighted_Mat"); 
      OutList.push_back(ptr_gSAIGEobj->m_scalefactor_G2_cond, "scalefactor_G2_cond");
      OutList.push_back(ptr_gSAIGEobj->m_VarInvMat_cond_scaled_weighted, "VarInvMat_G2_cond_scaled"); 
      OutList.push_back(ptr_gSAIGEobj->m_Tstat_cond, "Tstat_G2_cond"); //m_Tstat_cond is weighted
      OutList.push_back(ptr_gSAIGEobj->m_G2_Weight_cond, "G2_Weight_cond");
      OutList.push_back(TstatAdjCond, "TstatAdjCond");
      OutList.push_back(VarMatAdjCond, "VarMatAdjCond"); 

     

    }

    if(!t_isFastTest){
            iswriteOutput = true;
     }else{
        if(ptr_gSAIGEobj->m_flagSparseGRM_cur){
            iswriteOutput = true;
        }
     }

  //}//for(unsigned int i_mt = 0; i_mt < t_traitType.size(); i_mt++){
   
}//if(t_regionTestType == "BURDEN"){


//}


//std::cout << "If only conduct Burden test d " << std::endl;
 int numofUR = q_anno_maf_weight;
 int numofUR0;
 int mFirth = 0;
if(t_isSingleinGroupTest){
  OutList.push_back(pvalVec, "pvalVec");
if(iswriteOutput){

for(unsigned int j_mt = 0; j_mt < t_traitType.size(); j_mt++){  
  std::ofstream & t_OutFile_singleInGroup = OutFile_singleInGroup_vec.at(j_mt);
  //std::cout << "If only conduct Burden test e " << std::endl;
  numofUR0 = writeOutfile_singleInGroup(t_isMoreOutput,
      t_isImputation,
      isCondition,
      is_Firth,
      mFirth,
      is_FirthConverge,
      t_traitType,
      chrVec,
      posVec,
      markerVec,
      refVec,
      altVec,
      altCountsVec,
      altFreqVec,
      imputationInfoVec,
      missingRateVec,
      BetaVec,
      seBetaVec,
      TstatVec,
      varTVec,
      pvalVec,
      pvalNAVec,
      isSPAConvergeVec,
      Beta_cVec,
      seBeta_cVec,
      Tstat_cVec,
      varT_cVec,
      pval_cVec,
      pvalNA_cVec,
      AF_caseVec,
      AF_ctrlVec,
      N_caseVec,
      N_ctrlVec,
      N_case_homVec,
      N_ctrl_hetVec,
      N_case_hetVec,
      N_ctrl_homVec,
      N_Vec,
      t_OutFile_singleInGroup,
      j_mt,
      q);

}//for(unsigned int j_mt = 0; j_mt < t_traitType.size(); j_mt++){

}else{
  //OutFile_singleInGroup_temp.open(g_outputFilePrefixSingleInGroup_temp.c_str(), std::ofstream::out);

for(unsigned int j_mt = 0; j_mt < t_traitType.size(); j_mt++){
    std::ofstream & t_OutFile_singleInGroup_temp = OutFile_singleInGroup_temp_vec.at(j_mt);
  numofUR0 = writeOutfile_singleInGroup(t_isMoreOutput,
      t_isImputation,
      isCondition,
      is_Firth,
      mFirth,
      is_FirthConverge,
      t_traitType,
      chrVec,
      posVec,
      markerVec,
      refVec,
      altVec,
      altCountsVec,
      altFreqVec,
      imputationInfoVec,
      missingRateVec,
      BetaVec,
      seBetaVec,
      TstatVec,
      varTVec,
      pvalVec,
      pvalNAVec,
      isSPAConvergeVec,
      Beta_cVec,
      seBeta_cVec,
      Tstat_cVec,
      varT_cVec,
      pval_cVec,
      pvalNA_cVec,
      AF_caseVec,
      AF_ctrlVec,
      N_caseVec,
      N_ctrlVec,
      N_case_homVec,
      N_ctrl_hetVec,
      N_case_hetVec,
      N_ctrl_homVec,
      N_Vec,
      t_OutFile_singleInGroup_temp,
      j_mt,
      q);
      t_OutFile_singleInGroup_temp.close();
   }//for(unsigned int j_mt = 0; j_mt < t_traitType.size(); j_mt++){
 }//iswriteOutput
 OutList.push_back(numofUR0, "numofUR");
}//if(t_isSingleinGroupTest){

 OutList.push_back(iswriteOutput, "iswriteOutput");

 if(t_isOutputMarkerList){
	OutList.push_back(indicatorVec, "markerIndcatorVec");
 }

 
  OutList.push_back(NumRare_GroupVec, "NumRare_GroupVec");
  OutList.push_back(NumUltraRare_GroupVec, "NumUltraRare_GroupVec");


 if(t_regionTestType != "BURDEN" || t_isOutputMarkerList){
  OutList.push_back(annoMAFIndicatorMat, "annoMAFIndicatorMat");
 }

  //std::cout << "If only conduct Burden test g " << std::endl;


return OutList;
}



// [[Rcpp::export]]
void assign_conditionMarkers_factors(
                           std::string t_genoType,     // "plink", "bgen", "vcf"
			   std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           unsigned int t_n, 
			   arma::vec & t_weight_cond
			   )           // sample size
{
  ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
  bool isImpute = false;	
  unsigned int q = t_genoIndex.size();

  std::cout << "q " << q << std::endl;


  unsigned int nt = ptr_gSAIGEobj->m_traitType_vec.size();
  arma::mat P1Mat(q*nt, t_n);
  arma::mat P2Mat(t_n, q*nt);
  arma::mat VarInvMat(q, q*nt);
  arma::vec TstatVec(q*nt);
  arma::vec pVec(q*nt);
  arma::vec MAFVec(q*nt);
  arma::vec gyVec(q*nt);
  arma::vec w0G2_cond_Vec(q*nt);
  arma::vec gsumVec(t_n, arma::fill::zeros);
  arma::mat gsumMat(t_n, nt, arma::fill::zeros);
  //boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);
  arma::vec GVec(t_n);
  double Beta, seBeta, pval, pval_noSPA, Tstat, varT, gy, w0G2_cond;
  bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
  arma::vec P2Vec(t_n);

  //std::vector<uint> indexZeroVec;
  //std::vector<uint> indexNonZeroVec;
  double Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c;
  arma::rowvec G1tilde_P_G2tilde_Vec;
  bool isCondition = false;
  for(unsigned int i = 0; i < q; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::vector<uint> indexZeroVec;
    std::vector<uint> indexNonZeroVec;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;

    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false; 




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
            gIndex_prev = std::strtoull( t_genoIndex_prev_str.c_str(), &end_prev,10 );
        }else if(t_genoType == "plink"){
            t_genoIndex_prev_str = t_genoIndex.at(i-1);
            gIndex_prev = std::strtoull( t_genoIndex_prev_str.c_str(), &end_prev,10 );
        }
	//else if(t_genoType == "vcf"){
         //   gIndex_prev = 0;	
	//}
	std::remove(end_prev);
    }

    bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec, GVec, isImpute);
     //arma::vec GVec(GVec0);
     //GVec0.clear();
    if(!isReadMarker){
      break;
    }

    std::string info = chr+":"+std::to_string(pd)+"_"+ref+"/"+alt;

  double MAF = std::min(altFreq, 1 - altFreq);
  double MAC = MAF * 2 * t_n * (1 - missingRate);

  bool hasVarRatio;
 
  flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);


 arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
       indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
       indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);


  MAF = std::min(altFreq, 1 - altFreq);
  MAC = std::min(altCounts, 2*t_n-altCounts);

  arma::vec gtildeVec;

  std::cout << "herebefore " << std::endl;
  bool isSingleVarianceRatio;
  if(ptr_gSAIGEobj->m_varRatio_null_mt.n_rows == 1){
        //ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
        //ptr_gSAIGEobj->assignSingleVarianceRatio(false);
        isSingleVarianceRatio = true;
  }else{
        isSingleVarianceRatio = false;
  }


  for(unsigned int i_mt = 0; i_mt < nt; i_mt++){

     ptr_gSAIGEobj->assign_for_itrait(i_mt);
     unsigned int j_mt = i_mt*q+i;

    if(isSingleVarianceRatio){
  //std::cout << "Here2f mainMarkerInCPP" << std::endl;
       ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
       //std::cout << "Here2g mainMarkerInCPP" << std::endl;
    }else{
 // std::cout << "Here2gb mainMarkerInCPP" << std::endl;
       hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
 // std::cout << "Here2gc mainMarkerInCPP" << std::endl;
    }


  std::cout << "j_mt " << j_mt << std::endl;


  if(MAC > g_MACCutoffforER){
     Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
                    indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, true, ptr_gSAIGEobj->m_flagSparseGRM);
  }else{
     Unified_getMarkerPval(
                    GVec,
                    false, // bool t_isOnlyOutputNonZero,
                    indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, true, ptr_gSAIGEobj->m_flagSparseGRM);

  }

      std::cout << "after p value" << std::endl;
      ptr_gSAIGEobj->getadjG(GVec, gtildeVec);
      P1Mat.row(j_mt) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*gtildeVec.t();
      P2Mat.col(j_mt) = sqrt(ptr_gSAIGEobj->m_varRatioVal)*P2Vec;
      //P1Mat.row(i) = gtildeVec.t();
      //P2Mat.col(i) = P2Vec;
      MAFVec(j_mt) = MAF;
      //w0G2_cond = boost::math::pdf(beta_dist, MAF);
      //t_weight_cond.print();
      if(MAF == 0.0){
	std::cerr << "ERROR: Conditioning marker is monomorphic\n";
      }	      

      std::cout << "after p value 2" << std::endl;
     //if(!t_weight_cond.is_zero()){
     w0G2_cond = t_weight_cond(i);
    //}else{
//	 w0G2_cond = boost::math::pdf(beta_dist, MAF);
  //  }
     w0G2_cond_Vec(j_mt) = w0G2_cond;
     gyVec(j_mt) = gy * w0G2_cond;
     if(i_mt == 0){
     	gsumVec = gsumVec + GVec * w0G2_cond;
     }
           std::cout << "after p value 3" << std::endl;
     //gsumMat.col(i_mt) = gsumMat.col(i_mt) + GVec * w0G2_cond;
     TstatVec(j_mt) = Tstat;
     pVec(j_mt) = pval;

     }
  }


 arma::mat VarMatsub, P1Matsub, P2Matsub;
 arma::vec qsumVec(nt*q);
 arma::mat gsumtildeMat(gsumVec.n_elem, nt); 
 arma::vec gsumtildeVec;
 arma::mat VarMat(q, q*nt);

 gyVec.print("gyVec");
 for(unsigned int i_mt = 0; i_mt < nt; i_mt++){
     std::cout << "after p value 4" << std::endl;
     ptr_gSAIGEobj->assign_for_itrait(i_mt);
     P1Matsub = P1Mat.rows(i_mt*q, (i_mt+1)*q - 1);
     P2Matsub = P2Mat.cols(i_mt*q, (i_mt+1)*q - 1);
	//P1Matsub.print("P1Matsub");
	//P2Matsub.print("P2Matsub");	
	std::cout << "P1Matsub.n_rows " << P1Matsub.n_rows << std::endl;
	std::cout << "P1Matsub.n_cols " << P1Matsub.n_cols << std::endl;
	std::cout << "P2Matsub.n_rows " << P2Matsub.n_rows << std::endl;
	std::cout << "P2Matsub.n_cols " << P2Matsub.n_cols << std::endl;

     std::cout << "after p value 5" << std::endl;
     VarMatsub = P1Matsub * P2Matsub;
     VarMat.cols(i_mt*q, (i_mt+1)*q - 1) = VarMatsub;
     VarInvMat.cols(i_mt*q, (i_mt+1)*q - 1) = VarMatsub.i();
     std::cout << "after p value 6" << std::endl;
     


     qsumVec(i_mt) = arma::accu(gyVec.subvec(i_mt*q, (i_mt+1)*q - 1));
     ptr_gSAIGEobj->getadjG(gsumVec, gsumtildeVec);
     std::cout << "after p value 7" << std::endl;
     gsumtildeMat.col(i_mt) = gsumtildeVec;
     std::cout << "after p value 8" << std::endl;
}
  //double qsum = arma::accu(gyVec);


  ptr_gSAIGEobj->assignConditionFactors(
		   			P2Mat,
					VarInvMat,
					VarMat,
					TstatVec,
				        w0G2_cond_Vec,	
					MAFVec,
					qsumVec,
					gsumtildeMat,
					pVec);
  ptr_gSAIGEobj->m_VarInvMat_cond_scaled_weighted.resize(q, q*nt);					
}

// [[Rcpp::export]]
void assign_conditionMarkers_factors_binary_region(
			   arma::vec & scalefactor_G2_cond){
	//std::cout << "assign_conditionMarkers_factors_binary_region" << std::endl;
	ptr_gSAIGEobj->assignConditionFactors_scalefactor(scalefactor_G2_cond);
}

// [[Rcpp::export]]
void assign_conditionMarkers_factors_binary_region_multiTrait(
                           arma::mat & scalefactor_G2_cond,
			   unsigned int oml){
        //std::cout << "assign_conditionMarkers_factors_binary_region" << std::endl;
        ptr_gSAIGEobj->assignConditionFactors_scalefactor_multiTrait(scalefactor_G2_cond, oml);
}



// [[Rcpp::export]]
void set_iterator_inVcf(std::string & variantList, std::string & chrom, int & beg_pd, int & end_pd){
   if(!variantList.empty()){
	ptr_gVCFobj->set_iterator(variantList);	
   }else{
	ptr_gVCFobj->set_iterator(chrom, beg_pd, end_pd);
   }	   
}	

// [[Rcpp::export]]
bool check_Vcf_end(){
	bool isEnd = false;
	isEnd = ptr_gVCFobj->check_iterator_end();
	return(isEnd);
}


// [[Rcpp::export]]
void move_forward_iterator_Vcf(int i){
	ptr_gVCFobj->move_forward_iterator(i);
}



// [[Rcpp::export]]
arma::vec fast_logistf_fit(arma::mat & x,
		arma::vec & y,
		arma::vec & weight,
		arma::vec & offset,
		bool firth,
		arma::uvec & col_fit,
    	arma::vec init, 
	int maxit, 
	int maxstep, 
	int maxhs, 
	double lconv, 
	double gconv, 
	double xconv, 
	bool & isfirthconverge){
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  arma::vec pi_0 = -x * beta - offset; 
  pi_0 = arma::exp(pi_0) + 1;
  arma::vec pi_1 = 1/pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  isfirthconverge = false;
  while(iter <= maxit){
	beta_old = beta;
	arma::vec wpi = weight % pi_1 % (1 - pi_1);
	arma::vec wpi_sqrt = arma::sqrt(wpi);
	arma::vec W2 = weight % wpi_sqrt;
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
		ypih = (weight % (y - pi_1)) + (h % (0.5 - pi_1));
	}else{
		ypih = (weight % (y - pi_1));
	}
	//ypih.print();
	arma::vec xcol(n, arma::fill::zeros);
	U_star = x.t() * ypih;
	
	arma::mat XX_XW2(n, k, arma::fill::zeros);
	for(int j = 0; j < k; j++){
		xcol = x.col(j);
		XX_XW2.col(j) = xcol % wpi_sqrt; 	
	}
	arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
	bool isinv = arma::inv_sympd (XX_covs, XX_Fisher); 
	if(!isinv){
		break;
	}	
	//}
	arma::vec delta = XX_covs * U_star;
	delta.replace(arma::datum::nan, 0);	

	double mx = arma::max(arma::abs(delta))/maxstep;
	if(mx > 1){
		delta = delta/mx;
	}
	evals = evals + 1;
	iter = iter + 1;
	beta = beta + delta;
	pi_0 = -x * beta - offset;
  	pi_0 = arma::exp(pi_0) + 1;
  	pi_1 = 1/pi_0;
	if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
		isfirthconverge = true;
		break;
	}
  }
	arma::mat var;
	if(XX_covs.has_nan()){
		var = XX_covs;
		beta = arma::datum::nan;
	}
	return beta;
}

// [[Rcpp::export]]
void closeGenoFile(std::string & t_genoType)
{
  if(t_genoType == "bgen"){
    ptr_gBGENobj->closegenofile();
  }else if(t_genoType == "vcf"){
    ptr_gVCFobj->closegenofile();
  }else if(t_genoType == "plink"){
    ptr_gPLINKobj->closegenofile();
  }	  
}

// [[Rcpp::export]]
bool openOutfile(std::string t_traitType, bool isappend){
	bool isopen;
	std::ofstream OutFile;
	if(!isappend){
	OutFile.open(g_outputFilePrefixGroup.c_str());
	isopen = OutFile.is_open();
	if(isopen){
		OutFile << "Region\tGroup\tmax_MAF\tPvalue_Burden\tBETA_Burden\tSE_Burden\t";
		if(ptr_gSAIGEobj->m_isCondition){
			OutFile << "Pvalue_Burden_c\tBeta_Burden_c\tseBeta_Burden_c\t";
		}
		OutFile << "MAC\t";
		if(t_traitType == "binary"){	
			OutFile << "MAC_case\tMAC_control\t";
		}
		OutFile << "Number_rare\tNumber_ultra_rare\n";
	}
     }else{
	OutFile.open(g_outputFilePrefixGroup.c_str(), std::ofstream::out | std::ofstream::app) ;				
	isopen = OutFile.is_open();
     }

     OutFile_vec.push_back(std::move(OutFile));	

     return(isopen);
}

// [[Rcpp::export]]
bool openOutfile_singleinGroup(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput){
     bool isopen;
     std::ofstream OutFile_singleInGroup;
     std::ofstream OutFile_singleInGroup_temp;
     if(!isappend){
        OutFile_singleInGroup.open(g_outputFilePrefixSingleInGroup.c_str());
        OutFile_singleInGroup_temp.open(g_outputFilePrefixSingleInGroup_temp.c_str());
	isopen = OutFile_singleInGroup.is_open();
        if(isopen){
                OutFile_singleInGroup << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
		if(t_isImputation){
			OutFile_singleInGroup << "imputationInfo\t";
		}else{
			
			OutFile_singleInGroup << "MissingRate\t";
		}
		OutFile_singleInGroup << "BETA\tSE\tTstat\tvar\tp.value\t";
	        if(t_traitType == "binary" || t_traitType == "count"){
                        OutFile_singleInGroup << "p.value.NA\tIs.SPA\t";
                }	

                if(ptr_gSAIGEobj->m_isCondition){
                        OutFile_singleInGroup << "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
			if(t_traitType == "binary" || t_traitType == "count"){
				OutFile_singleInGroup << "p.value.NA_c\t";		
			}
                }
		
	        if(t_traitType == "binary"){
                        OutFile_singleInGroup << "AF_case\tAF_ctrl\tN_case\tN_ctrl";
 			if(t_isMoreOutput){
                                OutFile_singleInGroup << "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
                        }
			OutFile_singleInGroup << "\n";
                }else if(t_traitType == "quantitative" || t_traitType == "count"){
			OutFile_singleInGroup << "N\n";	
			
		}	

        }
    }else{
       OutFile_singleInGroup.open(g_outputFilePrefixSingleInGroup.c_str(), std::ofstream::out | std::ofstream::app);
       isopen = OutFile_singleInGroup.is_open();
    }
    OutFile_singleInGroup_vec.push_back(std::move(OutFile_singleInGroup));
    OutFile_singleInGroup_temp_vec.push_back(std::move(OutFile_singleInGroup_temp));	

    //if(t_testType != "BURDEN" && t_isfastTest){
    //  OutFile_singleInGroup_temp.open(g_outputFilePrefixSingleInGroup_temp.c_str(), std::ofstream::out | std::ofstream::app);
    //}
    return(isopen);
}


// [[Rcpp::export]]
void removeOutfile_singleinGroup_temp(){
	const char* File1 = g_outputFilePrefixSingleInGroup_temp.c_str();
	std::remove(File1);
	std::cout << "remove temp file " << g_outputFilePrefixSingleInGroup_temp << std::endl;
}


// [[Rcpp::export]]
void removeOutfile_inGroup(){
	const char* File1 = g_outputFilePrefixGroup.c_str();
	std::remove(File1);
}

// [[Rcpp::export]]
void removeOutfile_inSingle(){
	 std::cout << "remove empty file " << g_outputFilePrefix0 << std::endl;
	const char* File1 = g_outputFilePrefix0.c_str();
	std::remove(File1);
}



// [[Rcpp::export]]
bool openOutfile_single(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput,  bool t_isGbyE){
      bool isopen;
      std::ofstream OutFile_single;
      if(!isappend){
        OutFile_single.open(g_outputFilePrefixSingle.c_str());
        isopen = OutFile_single.is_open();
        if(isopen){
                OutFile_single << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
                if(t_isImputation){
                        OutFile_single << "imputationInfo\t";
                }else{

                        OutFile_single << "MissingRate\t";
                }
                OutFile_single << "BETA\tSE\tTstat\tvar\tp.value\t";
                if(t_traitType == "binary" || t_traitType == "count"){
                        OutFile_single << "p.value.NA\tIs.SPA\t";
                }

                if(ptr_gSAIGEobj->m_isCondition){
                        OutFile_single << "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
			if(t_traitType == "binary"){
				OutFile_single << "p.value.NA_c\t";
			}
                }
		

                if(t_traitType == "binary"){
                        OutFile_single << "AF_case\tAF_ctrl\tN_case\tN_ctrl";

			if(t_isMoreOutput){	
				OutFile_single << "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
			}
			//OutFile_single << "\n";

                }else if(t_traitType == "quantitative" || t_traitType == "count" || t_traitType == "count_nb"){
                        OutFile_single << "N";

                }

		if(t_isGbyE){

                        OutFile_single << "\tBeta_ge\tseBeta_ge\tpval_ge\tpval_noSPA_ge\tpval_ge_SKAT";
                }

                OutFile_single << "\n";

	
        }
     }else{
        OutFile_single.open(g_outputFilePrefixSingle.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = OutFile_single.is_open();
     }
     OutFile_single_vec.push_back(std::move(OutFile_single));
     return(isopen);
}


void writeOutfile_single(bool t_isMoreOutput,
			bool t_isImputation,
			bool t_isCondition,
			bool t_isFirth,
                         int mFirth,
                         int mFirthConverge,
                        std::vector<std::string> & t_traitType,
                        std::vector<std::string> & chrVec,
                        std::vector<std::string> & posVec,
                        std::vector<std::string> & markerVec,
                        std::vector<std::string> & refVec,
                        std::vector<std::string> & altVec,
                        std::vector<double> & altCountsVec,
                        std::vector<double> & altFreqVec,
                        std::vector<double> & imputationInfoVec,
                        std::vector<double> & missingRateVec,
                        std::vector<double> & BetaVec,
                        std::vector<double> & seBetaVec,
                        std::vector<double> & TstatVec,
                        std::vector<double> & varTVec,
                        std::vector<double> & pvalVec,
                        std::vector<double> & pvalNAVec,
                        std::vector<bool>  & isSPAConvergeVec,
                        std::vector<double> & Beta_cVec,
                        std::vector<double> & seBeta_cVec,
                        std::vector<double> & Tstat_cVec,
                        std::vector<double> & varT_cVec,
                        std::vector<double> & pval_cVec,
                        std::vector<double> & pvalNA_cVec,
                        std::vector<double> & AF_caseVec,
                        std::vector<double> & AF_ctrlVec,
                        std::vector<uint32_t> & N_caseVec,
                        std::vector<uint32_t> & N_ctrlVec,
                        std::vector<double>  & N_case_homVec,
                        std::vector<double>  & N_ctrl_hetVec,
                        std::vector<double>  & N_case_hetVec,
                        std::vector<double>  & N_ctrl_homVec,
                        std::vector<uint32_t> & N_Vec,
			bool is_GbyE,
			std::vector<std::string> & Beta_ge_cStrVec,
			std::vector<std::string> & seBeta_ge_cStrVec,
			std::vector<std::string> & pval_ge_cStrVec,
			std::vector<std::string> & pval_noSPA_ge_cStrVec,
			std::vector<double> & pval_SKAT_ge_cVec,
			unsigned int itt,
			unsigned int nmarkers
){
  int numtest = 0;
  //unsigned int nmarkers = pvalVec.size()/(t_traitType.size());
  unsigned int k;
  //std::cout << "nmarkers " << nmarkers << std::endl; 
  //std::cout << "t_traitType.size() " << t_traitType.size() << std::endl; 
  //for(unsigned int k = 0; k < pvalVec.size(); k++){
//for(unsigned int i = 0; i < t_traitType.size(); i++){

//std::cout << "i " << i << std::endl;
  std::ofstream & OutFile_single = OutFile_single_vec.at(itt);
  numtest = 0;
     for(unsigned int j = 0; j < nmarkers; j++){
  	k = itt*nmarkers + j;
        if(!std::isnan(pvalVec.at(k))){
                numtest = numtest + 1;
		//OutFile_single << std::to_string(i);
		//OutFile_single << "\t";
                OutFile_single << chrVec.at(k);
                OutFile_single << "\t";
                OutFile_single << posVec.at(k);
                OutFile_single << "\t";
                OutFile_single << markerVec.at(k);
                OutFile_single << "\t";
                OutFile_single << refVec.at(k);
                OutFile_single << "\t";
                OutFile_single << altVec.at(k);
                OutFile_single << "\t";
                OutFile_single << altCountsVec.at(k);
                OutFile_single << "\t";
                OutFile_single << altFreqVec.at(k);
                OutFile_single << "\t";

                if(t_isImputation){
                        OutFile_single << imputationInfoVec.at(k);
                        OutFile_single << "\t";
                }else{
                        OutFile_single << missingRateVec.at(k);
                        OutFile_single << "\t";

                }
                OutFile_single << BetaVec.at(k);
                OutFile_single << "\t";
                OutFile_single << seBetaVec.at(k);
                OutFile_single << "\t";
                OutFile_single << TstatVec.at(k);
                OutFile_single << "\t";
                OutFile_single << varTVec.at(k);
                OutFile_single << "\t";
                OutFile_single << pvalVec.at(k);
                OutFile_single << "\t";

                if(t_traitType.at(itt) == "binary" || t_traitType.at(itt) == "count"){
                        OutFile_single << pvalNAVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << std::boolalpha << isSPAConvergeVec.at(k);
                        OutFile_single << "\t";
                }
                if(t_isCondition){
                        OutFile_single << Beta_cVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << seBeta_cVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << Tstat_cVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << varT_cVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << pval_cVec.at(k);
                        OutFile_single << "\t";
			if(t_traitType.at(itt) == "binary" || t_traitType.at(itt) == "count"){	
                        	OutFile_single << pvalNA_cVec.at(k);
                        	OutFile_single << "\t";
			}	
                }
                if(t_traitType.at(itt) == "binary"){
                        OutFile_single << AF_caseVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << AF_ctrlVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << N_caseVec.at(k);
                        OutFile_single << "\t";
                        OutFile_single << N_ctrlVec.at(k);

                        if(t_isMoreOutput){
                                OutFile_single << "\t";
                                OutFile_single << N_case_homVec.at(k);
                                OutFile_single << "\t";
                                OutFile_single << N_case_hetVec.at(k);
                                OutFile_single << "\t";
                                OutFile_single << N_ctrl_homVec.at(k);
                                OutFile_single << "\t";
                                OutFile_single << N_ctrl_hetVec.at(k);
                        }
                        //OutFile_single << "\n";
                }else if(t_traitType.at(itt) == "quantitative" || t_traitType.at(itt) == "count" || t_traitType.at(itt) == "count_nb"){
                        OutFile_single << N_Vec.at(k);
                        //OutFile_single << "\n";

                }


		if(is_GbyE){

			OutFile_single << "\t";
			OutFile_single << Beta_ge_cStrVec.at(k);
			OutFile_single << "\t";
			OutFile_single << seBeta_ge_cStrVec.at(k);
			OutFile_single << "\t";
			OutFile_single << pval_ge_cStrVec.at(k);
			OutFile_single << "\t";
			OutFile_single << pval_noSPA_ge_cStrVec.at(k);
			OutFile_single << "\t";
			OutFile_single << pval_SKAT_ge_cVec.at(k);
		}

		OutFile_single << "\n";

        }
  }

  if(itt == 0){
  std::cout << numtest << " markers were tested." << std::endl;
  }
  if(t_traitType.at(itt) == "binary"){
      if(t_isFirth){
        std::cout << "Firth approx was applied to " << mFirth << " markers. " << mFirthConverge << " sucessfully converged." <<std::endl;
      }
   }


  //}  
}


// [[Rcpp::export]]
void set_flagSparseGRM_cur_SAIGE(bool t_flagSparseGRM_cur){
	ptr_gSAIGEobj->set_flagSparseGRM_cur(t_flagSparseGRM_cur);
}

// [[Rcpp::export]]
void set_flagSparseGRM_cur_SAIGE_org(){
	ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
}



void set_varianceRatio(double MAC, bool isSingleVarianceRatio, bool isnoXadj){
    bool hasVarRatio;
    if(!ptr_gSAIGEobj->m_isFastTest){
       ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
       if(!isSingleVarianceRatio){
         hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM, isnoXadj);
       }else{
         ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM, isnoXadj);
       }
     }else{
       ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
       if(!isSingleVarianceRatio){
         hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, false, isnoXadj);
       }else{
         ptr_gSAIGEobj->assignSingleVarianceRatio(false, isnoXadj);
       }
     }
}




void writeOutfile_BURDEN(std::string regionName,
			std::vector<std::string>  & BURDEN_AnnoName_Vec,
			std::vector<std::string> & BURDEN_maxMAFName_Vec,
			arma::vec & BURDEN_pval_Vec,
			std::vector<double> & BURDEN_Beta_Vec,
			std::vector<double> & BURDEN_seBeta_Vec,
			arma::vec & BURDEN_pval_cVec,
			std::vector<double> & BURDEN_Beta_cVec,
			std::vector<double> & BURDEN_seBeta_cVec,
			arma::vec & MAC_GroupVec,
			arma::vec & MACCase_GroupVec,
			arma::vec & MACControl_GroupVec,
			arma::vec & NumRare_GroupVec,
			arma::vec & NumUltraRare_GroupVec,
			double cctpval,
			double cctpval_cond,
			unsigned int  q_weight,
			unsigned int q_anno,
			unsigned int q_maf,
			bool isCondition,
			unsigned int itt){
     std::ofstream & OutFile = OutFile_vec.at(itt); 			
     unsigned int i, q;
   for(unsigned int r = 0; r < q_weight; r++){  
     for(unsigned int j = 0; j < q_anno; j++){
       for(unsigned int m = 0; m < q_maf; m++){
           //i = j*q_maf+m;
	   q = j*q_anno + m; 
	   i = itt*q_anno*q_maf*q_weight + j*q_weight*q_anno+m;
	   if(BURDEN_pval_Vec(i) != -1){
           OutFile << regionName;
           OutFile << "\t";
           OutFile << BURDEN_AnnoName_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_maxMAFName_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_pval_Vec(i);
           OutFile << "\t";
           OutFile << BURDEN_Beta_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_seBeta_Vec.at(i);
           OutFile << "\t";
           if(isCondition){
               OutFile << BURDEN_pval_cVec(i);
               OutFile << "\t";
               OutFile << BURDEN_Beta_cVec.at(i);
               OutFile << "\t";
               OutFile << BURDEN_seBeta_cVec.at(i);
               OutFile << "\t";
	   }
	   OutFile << MAC_GroupVec(i);
           OutFile << "\t";
           if(ptr_gSAIGEobj->m_traitType_vec.at(itt) == "binary"){
               OutFile << MACCase_GroupVec(i);
               OutFile << "\t";
               OutFile << MACControl_GroupVec(i);
               OutFile << "\t";
           }
           OutFile << NumRare_GroupVec(i);
           OutFile << "\t";
           OutFile << NumUltraRare_GroupVec(i);
           OutFile << "\n";
	}
     }
   }
  } 
     OutFile << regionName;
     OutFile << "\tCauchy\tNA\t";
     OutFile << cctpval;
     OutFile << "\tNA\tNA\t";	
     if(isCondition){
	OutFile << cctpval_cond;
	OutFile << "\tNA\tNA\t";
     }
     OutFile << "NA\t";
     if(ptr_gSAIGEobj->m_traitType_vec.at(itt) == "binary"){
        OutFile << "NA\t";
        OutFile << "NA\t";
     }
     OutFile << "NA\t";
     OutFile << "NA\n";	
}


int writeOutfile_singleInGroup(bool t_isMoreOutput,
                        bool t_isImputation,
                        bool t_isCondition,
                        bool t_isFirth,
                         int mFirth,
                         int mFirthConverge,
                        std::vector<std::string> & t_traitType,
                        std::vector<std::string> & chrVec,
                        std::vector<std::string> & posVec,
                        std::vector<std::string> & markerVec,
                        std::vector<std::string> & refVec,
                        std::vector<std::string> & altVec,
                        std::vector<double> & altCountsVec,
                        std::vector<double> & altFreqVec,
                        std::vector<double> & imputationInfoVec,
                        std::vector<double> & missingRateVec,
                        std::vector<double> & BetaVec,
                        std::vector<double> & seBetaVec,
                        std::vector<double> & TstatVec,
                        std::vector<double> & varTVec,
                        std::vector<double> & pvalVec,
                        std::vector<double> & pvalNAVec,
                        std::vector<bool>  & isSPAConvergeVec,
                        std::vector<double> & Beta_cVec,
                        std::vector<double> & seBeta_cVec,
                        std::vector<double> & Tstat_cVec,
                        std::vector<double> & varT_cVec,
                        std::vector<double> & pval_cVec,
                        std::vector<double> & pvalNA_cVec,
                        std::vector<double> & AF_caseVec,
                        std::vector<double> & AF_ctrlVec,
                        std::vector<uint32_t> & N_caseVec,
                        std::vector<uint32_t> & N_ctrlVec,
                        std::vector<double>  & N_case_homVec,
                        std::vector<double>  & N_ctrl_hetVec,
                        std::vector<double>  & N_case_hetVec,
                        std::vector<double>  & N_ctrl_homVec,
                        std::vector<uint32_t> & N_Vec,
                        std::ofstream & t_OutFile_singleInGroup, 
			unsigned int itt,
                        unsigned int nmarkers){
  int numofUR = 0;
  unsigned int k;
  //for(unsigned int k = 0; k < pvalVec.size(); k++){
  //for(unsigned int k = 0; k < nmarkers; k++){
  //std::ofstream & OutFile_single = OutFile_single_vec.at(itt);
     //numtest = 0;
     //std::cout << "If only conduct Burden test f " << std::endl;
     for(unsigned int j = 0; j < nmarkers; j++){
        k = itt*nmarkers + j;    
	//std::cout << "k " << k << std::endl;
    //if(std::isfinite(pvalVec.at(k))){
        if(!std::isnan(pvalVec.at(k))){
		if(chrVec.at(j) == "UR"){
                        numofUR = numofUR + 1;
                }
                t_OutFile_singleInGroup << chrVec.at(j);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << posVec.at(j);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << markerVec.at(j);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << refVec.at(j);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altVec.at(j);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altCountsVec.at(j);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altFreqVec.at(j);
                t_OutFile_singleInGroup << "\t";

                if(t_isImputation){
                        t_OutFile_singleInGroup << imputationInfoVec.at(j);
                        t_OutFile_singleInGroup << "\t";
                }else{
                        t_OutFile_singleInGroup << missingRateVec.at(j);
                        t_OutFile_singleInGroup << "\t";

                }
                t_OutFile_singleInGroup << BetaVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << seBetaVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << TstatVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << varTVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << pvalVec.at(k);
                t_OutFile_singleInGroup << "\t";
	//std::cout << "k c " << k << std::endl;

                if(t_traitType.at(itt) == "binary" || t_traitType.at(itt) == "count"){
                        t_OutFile_singleInGroup << pvalNAVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << std::boolalpha << isSPAConvergeVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                }
	//std::cout << "k b " << k << std::endl;
                if(t_isCondition){
                        t_OutFile_singleInGroup << Beta_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << seBeta_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << Tstat_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << varT_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << pval_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
			if(t_traitType.at(itt) == "binary" || t_traitType.at(itt) == "count"){
                        	t_OutFile_singleInGroup << pvalNA_cVec.at(k);
                        	t_OutFile_singleInGroup << "\t";
			}
                }
                if(t_traitType.at(itt) == "binary"){
                        t_OutFile_singleInGroup << AF_caseVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << AF_ctrlVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << N_caseVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << N_ctrlVec.at(k);

                        if(t_isMoreOutput){
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_case_homVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_case_hetVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_ctrl_homVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_ctrl_hetVec.at(k);
                        }
                        t_OutFile_singleInGroup << "\n";
                }else if(t_traitType.at(itt) == "quantitative" || t_traitType.at(itt) == "count"){
                        t_OutFile_singleInGroup << N_Vec.at(k);
                        t_OutFile_singleInGroup << "\n";

                }
        }
  }
  return(numofUR);

}

// [[Rcpp::export]]
void copy_singleInGroup(){

  for(unsigned int j_mt = 0; j_mt < ptr_gSAIGEobj->m_traitType_vec.size(); j_mt++){

  std::ifstream ini_file;
  std::ofstream & t_OutFile_singleInGroup = OutFile_singleInGroup_vec.at(j_mt);
  std::string t_outputFilePrefixSingleInGroup_temp = g_outputFilePrefix0+"_"+std::to_string(j_mt)+".singleAssoc.txt_temp";
  ini_file.open(t_outputFilePrefixSingleInGroup_temp.c_str());
  if (!ini_file)
  {
    std::cout << "Error in Opening the temp file!" << std::endl;
  }
  std::string str;
  while (getline(ini_file, str))
  {
    t_OutFile_singleInGroup << str << "\n"; 	
  }
  ini_file.close();

  }
}



// [[Rcpp::export]]
void set_dup_sample_index(arma::uvec & t_dup_sample_Index){
        g_dup_sample_Index = t_dup_sample_Index;
        arma::uvec uniq_sample_Inddex  = arma::find_unique(g_dup_sample_Index);
        g_n_unique = uniq_sample_Inddex.n_elem;
}


// [[Rcpp::export]]
void setupSparseGRM_new(arma::sp_mat & t_spGRM){
        arma::sp_mat t_spGRM_1 = t_spGRM;
        arma::sp_fmat g_spGRM_f = arma::conv_to<arma::sp_fmat>::from(t_spGRM_1);
        g_spGRM = g_spGRM_f;
}


// [[Rcpp::export]]
arma::sp_mat getSparseSigma_new(){
	arma::sp_mat g_spGRM_double = arma::conv_to<arma::sp_mat>::from(g_spGRM);
	return(g_spGRM_double);
}	


// [[Rcpp::export]]
void set_I_longl_mat(arma::sp_mat & t_Ilongmat, arma::vec & t_I_longl_vec){
        arma::sp_fmat t_Kmat_new = arma::conv_to< arma::sp_fmat >::from(t_Ilongmat);
        g_I_longl_mat = t_Kmat_new;
        arma::uvec t_I_longl_vec_new = arma::conv_to< arma::uvec >::from(t_I_longl_vec);
        g_I_longl_vec = t_I_longl_vec_new;

}

// [[Rcpp::export]]
void set_I_longl_mat_SAIGEtest(arma::sp_mat & t_Ilongmat, arma::vec & t_I_longl_vec){
        //arma::sp_fmat t_Kmat_new = arma::conv_to< arma::sp_fmat >::from(t_Ilongmat);
	std::cout << "ok1 " << std::endl;
	int n = (ptr_gSAIGEobj->g_I_longl_mat).n_rows;
	std::cout << "n " << n << std::endl;
	ptr_gSAIGEobj->g_I_longl_mat = t_Ilongmat;
	std::cout << "ok2 " << std::endl;
        arma::uvec t_I_longl_vec_new = arma::conv_to< arma::uvec >::from(t_I_longl_vec);
        ptr_gSAIGEobj->g_I_longl_vec = t_I_longl_vec_new;
	std::cout << "ok3 " << std::endl;

}


// [[Rcpp::export]]
void set_T_longl_mat(arma::sp_mat & t_Tlongmat, arma::vec & t_T_longl_vec){
        arma::sp_fmat t_Kmat_new = arma::conv_to< arma::sp_fmat >::from(t_Tlongmat);
        g_T_longl_mat = t_Kmat_new;
        arma::fvec t_T_longl_vec_new = arma::conv_to< arma::fvec >::from(t_T_longl_vec);
        g_T_longl_vec = t_T_longl_vec_new;
}


// [[Rcpp::export]]
void set_T_longl_mat_SAIGEtest(arma::sp_mat & t_Tlongmat, arma::vec & t_T_longl_vec){
        //arma::sp_fmat t_Kmat_new = arma::conv_to< arma::sp_fmat >::from(t_Tlongmat);
        ptr_gSAIGEobj->g_T_longl_mat = t_Tlongmat;
        //arma::fvec t_T_longl_vec_new = arma::conv_to< arma::fvec >::from(t_T_longl_vec);
        ptr_gSAIGEobj->g_T_longl_vec = t_T_longl_vec;
}



// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin(arma::fcolvec& bVec, bool LOCO){
    arma::fvec crossProdVec;
    if(g_isSparseGRM){
        crossProdVec = g_spGRM * bVec;
    }else{
	if(!LOCO){    
          crossProdVec = parallelCrossProd(bVec);
	}else{
	  crossProdVec = parallelCrossProd_LOCO(bVec);
	}	
    }
    return(crossProdVec);
}




// [[Rcpp::export]]
arma::fcolvec getCrossprod_multiV(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, bool LOCO){

        arma::fcolvec crossProdVec;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0 && tauVec.n_elem == 2){
                crossProdVec = tauVec(0)*(bVec % (1/wVec));
                return(crossProdVec);
        }
        //
        arma::fvec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;

        unsigned int tau_ind = 0;
        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){ //it must have specified GRM
		if(g_isGRM){
                  crossProd1 = getCrossprodMatAndKin(bVec, LOCO);
                  crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;
                  tau_ind = tau_ind + 2;
		}else{
                  crossProdVec = tauVec(0)*(bVec % (1/wVec));
		  tau_ind = tau_ind + 1;
	        }		
        }else{
                Ibvec = g_I_longl_mat.t() * bVec;
		if(g_isGRM){
                  GRM_I_bvec = getCrossprodMatAndKin(Ibvec, LOCO);
                  crossProd1 = g_I_longl_mat * GRM_I_bvec;
                  crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;
                  tau_ind = tau_ind + 2;
		}else{
                  crossProdVec = tauVec(0)*(bVec % (1/wVec));
		  tau_ind = tau_ind + 1;
		}	


                if(g_T_longl_mat.n_rows > 0){
                        Tbvec = g_T_longl_mat.t() * bVec;
			if(g_isGRM){
                          GRM_T_bvec = getCrossprodMatAndKin(Tbvec, LOCO);
                          crossProdGRM_TGIb = g_T_longl_mat * GRM_I_bvec;
                          crossProdGRM_IGTb = g_I_longl_mat * GRM_T_bvec;
                          crossProdVec = crossProdVec + tauVec(tau_ind) * (crossProdGRM_TGIb + crossProdGRM_IGTb);
                          tau_ind = tau_ind + 1;
                          crossProdVec = crossProdVec + tauVec(tau_ind) * (g_T_longl_mat * GRM_T_bvec);
                          tau_ind = tau_ind + 1;
			}
                }

        }


        if(g_T_longl_mat.n_rows == 0 && g_I_longl_mat.n_rows == 0){

                if(Kmat_vec.size() > 0){
                        for(unsigned int i = 0; i < Kmat_vec.size(); i++){
                                crossProdVec  = crossProdVec + tauVec(tau_ind)*(Kmat_vec[i] * bVec);
                                tau_ind = tau_ind + 1;
                        }
                }

        }else{ //if(g_T_longl_mat.n_rows == 0 && g_I_longl_mat.n_rows == 0){
                crossProdVec  = crossProdVec + tauVec(tau_ind) * (g_I_longl_mat * Ibvec);
                tau_ind = tau_ind + 1;

                if(g_T_longl_mat.n_rows > 0){
                        crossProdGRM_TIb = g_T_longl_mat * Ibvec;
                        crossProdGRM_ITb = g_I_longl_mat * Tbvec;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (crossProdGRM_TIb + crossProdGRM_ITb);
                        tau_ind = tau_ind + 1;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (g_T_longl_mat * Tbvec);
                        tau_ind = tau_ind + 1;
                        //crossProdGRM_TGIb = g_T_longl_mat * GRM_I_bvec;
                        //crossProdGRM_IGTb = g_I_longl_mat * GRM_T_bvec;
                }


                if(Kmat_vec.size() > 0){
                        for(unsigned int i = 0; i < Kmat_vec.size(); i++){
                                V_I_bvec = Kmat_vec[i] * Ibvec;
                                crossProdVec  = crossProdVec + tauVec(tau_ind) * (g_I_longl_mat * V_I_bvec);
                                tau_ind = tau_ind + 1;
                                if(g_T_longl_mat.n_rows > 0){
                                        V_T_bvec = Kmat_vec[i] * Tbvec;
                                        crossProdV_TGIb = g_T_longl_mat * V_I_bvec;
                                        crossProdV_IGTb = g_I_longl_mat * V_T_bvec;
                                        crossProdVec = crossProdVec + tauVec(tau_ind) * (crossProdV_TGIb + crossProdV_IGTb);
                                        tau_ind = tau_ind + 1;
                                        crossProdVec = crossProdVec + tauVec(tau_ind) * (g_T_longl_mat * V_T_bvec);
                                        tau_ind = tau_ind + 1;
                                }
                        }
                }

        }
        return(crossProdVec);
}





// [[Rcpp::export]]
arma::fvec getDiagOfSigma_multiV(arma::fvec& wVec, arma::fvec& tauVec, bool LOCO){
        int Nnomissing = wVec.n_elem;
        arma::fvec diagVec(Nnomissing);
        arma::sp_fvec diagVecV0;
        arma::fvec diagVecG, diagVecV, diagVecG_I, diagVecG_T, diagVecG_IT,diagVecV_I, diagVecV_T, diagVecV_IT;
        unsigned int tauind = 0;
        
	if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){

         if(g_isGRM){		
          if(!(ptr_gNULLGENOobj->setKinDiagtoOne)){
           if(!g_isSparseGRM){		  
            if(!LOCO){
              int MminMAF = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
              diagVec = tauVec(1)* (*ptr_gNULLGENOobj->Get_Diagof_StdGeno()) /MminMAF + tauVec(0)/wVec;
            }else{
              diagVec = tauVec(1)* (*ptr_gNULLGENOobj->Get_Diagof_StdGeno_LOCO());
              int Msub_MAFge_minMAFtoConstructGRM_in_b = ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_in();
              int Msub_MAFge_minMAFtoConstructGRM_singleVar_b = ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_singleChr_in();
              diagVec = diagVec/(Msub_MAFge_minMAFtoConstructGRM_in_b - Msub_MAFge_minMAFtoConstructGRM_singleVar_b) + tauVec(0)/wVec;
            }
	    tauind = tauind + 2;
	   }else{
	     diagVec = tauVec(0)/wVec;
	     tauind = tauind + 1;
	     diagVecG = g_spGRM.diag();
	     diagVec = diagVec + tauVec(tauind) * diagVecG;
             tauind = tauind + 1;
	
	   }		   
          }else{ //if(!(ptr_gNULLGENOobj->setKinDiagtoOne)){
            diagVec = tauVec(1) + tauVec(0)/wVec;
	    tauind = tauind + 2;
          }
	}else{//if(g_isGRM){
          diagVec = tauVec(0)/wVec;
	  tauind = tauind + 1;
	}	


          if(Kmat_vec.size() > 0){
            for(unsigned int i = 0; i < Kmat_vec.size(); i++){
              diagVec = diagVec + (Kmat_vec[i]).diag() * tauVec(tauind);
            }
          }

        }else{ //if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
                diagVec = tauVec(0)/wVec;
                tauind = tauind + 1;
	  if(g_isGRM && g_isSparseGRM){
                diagVecG = g_spGRM.diag();
                diagVecG_I = diagVecG.elem(g_I_longl_vec);
                diagVec = diagVec + tauVec(tauind) * diagVecG_I;
                tauind = tauind + 1;

                if(g_T_longl_mat.n_rows > 0){
		//std::cout << "getDiagOfSigma_multiV Here1" << std::endl;
		//diagVecG_I.print("diagVecG_I");
                  diagVecG_IT = diagVecG_I % g_T_longl_vec;
		//std::cout << "getDiagOfSigma_multiV Here2" << std::endl;
                  diagVecG_T = diagVecG_IT % g_T_longl_vec;
		//std::cout << "getDiagOfSigma_multiV Here3" << std::endl;
                  diagVecG_IT = 2 * diagVecG_IT;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_IT;
                  tauind = tauind + 1;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_T;
                  tauind = tauind + 1;
                }
	    }else{
		diagVecG.ones(g_n_unique);
		diagVecG_I = diagVecG.elem(g_I_longl_vec);
                diagVec = diagVec + tauVec(tauind) * diagVecG_I;
		tauind = tauind + 1;
		if(g_T_longl_mat.n_rows > 0){
                  //std::cout << "getDiagOfSigma_multiV Here1" << std::endl;
                  diagVecG_IT = diagVecG_I % g_T_longl_vec;
                  //std::cout << "getDiagOfSigma_multiV Here2" << std::endl;
                  diagVecG_T = diagVecG_IT % g_T_longl_vec;
                  //std::cout << "getDiagOfSigma_multiV Here3" << std::endl;
                  diagVecG_IT = 2 * diagVecG_IT;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_IT;
                  tauind = tauind + 1;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_T;
                  tauind = tauind + 1;
                  //std::cout << "getDiagOfSigma_multiV Here4" << std::endl;
                }

           } 		    

	

                if(Kmat_vec.size() > 0){
                  for(unsigned int i = 0; i < Kmat_vec.size(); i++){
                    diagVecV0 = Kmat_vec[i].diag();
                    arma::fvec diagVecVtemp(diagVecV0);
                    diagVecV = diagVecVtemp;
                    diagVecV_I = diagVecV.elem(g_I_longl_vec);
                    diagVec = diagVec + tauVec(tauind) * diagVecV_I;
                    tauind = tauind + 1;
                    if(g_T_longl_mat.n_rows > 0){
                      diagVecV_IT = diagVecV_I % g_T_longl_vec;
                      diagVecV_T = diagVecV_IT % g_T_longl_vec;
                      diagVecV_IT = 2 * diagVecV_IT;
                      diagVec = diagVec + tauVec(tauind) * diagVecV_IT;
                      tauind = tauind + 1;
                      diagVec = diagVec + tauVec(tauind) * diagVecV_T;
                      tauind = tauind + 1;
                    }
                  }
                }

        }

        for(unsigned int i=0; i< Nnomissing; i++){
                if(diagVec(i) < 1e-4){
                        diagVec(i) = 1e-4 ;
                }
        }

        return(diagVec);
}



// [[Rcpp::export]]
void gen_sp_Sigma_multiV(arma::fvec& wVec,  arma::fvec& tauVec){
   arma::fvec dtVec = (1/wVec) * (tauVec(0));
   arma::sp_fmat GRM_Imat, GRM_Tmat;
   arma::fvec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;
   unsigned int tau_ind = 0;
   //arma::sp_fmat g_spGRM_f = arma::conv_to< arma::sp_fmat >::from(g_spGRM);

   if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){

       if(g_isGRM && g_isSparseGRM){ 
         g_spSigma = g_spGRM * tauVec(1);
         g_spSigma.diag() = g_spSigma.diag() + dtVec;
         tau_ind = tau_ind + 2;
	}else{
          g_spSigma.zeros(wVec.n_elem, wVec.n_elem);
          g_spSigma.diag() = g_spSigma.diag() + dtVec;
	  tau_ind = tau_ind + 1;	  
	}	
   }else{
       if(g_isGRM && g_isSparseGRM){	   
         GRM_Imat = g_spGRM * (g_I_longl_mat.t());
         g_spSigma = g_I_longl_mat * GRM_Imat;
         g_spSigma = g_spSigma * tauVec(1);
	 g_spSigma.diag() = g_spSigma.diag() + dtVec;
	 tau_ind = tau_ind + 2;
         if(g_T_longl_mat.n_rows > 0){
           GRM_Tmat = g_spGRM * (g_T_longl_mat.t());
           g_spSigma = g_spSigma + tauVec(tau_ind) * (g_T_longl_mat * GRM_Imat + g_I_longl_mat * GRM_Tmat);
           tau_ind = tau_ind + 1;
           g_spSigma = g_spSigma + tauVec(tau_ind) * (g_T_longl_mat * GRM_Tmat);
           tau_ind = tau_ind + 1;
         }
       }else{
	 g_spSigma.zeros(wVec.n_elem, wVec.n_elem);      
         g_spSigma.diag() = g_spSigma.diag() + dtVec;
         tau_ind = tau_ind + 1;
	}
   }


   if(g_T_longl_mat.n_rows == 0 && g_I_longl_mat.n_rows == 0){
        if(Kmat_vec.size() > 0){
             for(unsigned int i = 0; i < Kmat_vec.size(); i++){
                     g_spSigma = g_spSigma + tauVec(tau_ind)*(Kmat_vec[i]);
                     tau_ind = tau_ind + 1;
             }
         }
   }else{
       g_spSigma = g_spSigma + tauVec(tau_ind) * g_I_longl_mat * (g_I_longl_mat.t());
       tau_ind = tau_ind + 1;

       if(g_T_longl_mat.n_rows > 0){

           g_spSigma = g_spSigma + tauVec(tau_ind) * (g_T_longl_mat * (g_I_longl_mat.t()) + g_I_longl_mat * (g_T_longl_mat.t()));
           tau_ind = tau_ind + 1;
           g_spSigma = g_spSigma + tauVec(tau_ind) * (g_T_longl_mat * (g_T_longl_mat.t()));
           tau_ind = tau_ind + 1;
                        //crossProdGRM_TGIb = g_T_longl_mat * GRM_I_bvec;
                        //crossProdGRM_IGTb = g_I_longl_mat * GRM_T_bvec;
       }

       if(Kmat_vec.size() > 0){
           for(unsigned int i = 0; i < Kmat_vec.size(); i++){
                GRM_Imat = Kmat_vec[i] * (g_I_longl_mat.t());
                g_spSigma = g_spSigma + tauVec(tau_ind) * (g_I_longl_mat * GRM_Imat);
                tau_ind = tau_ind + 1;
                if(g_T_longl_mat.n_rows > 0){
                        GRM_Tmat = Kmat_vec[i] * (g_T_longl_mat.t());
                        g_spSigma = g_spSigma + tauVec(tau_ind) * ((g_T_longl_mat * GRM_Imat) + (g_I_longl_mat * GRM_Tmat));
                        tau_ind = tau_ind + 1;
                        g_spSigma = g_spSigma + tauVec(tau_ind) * (g_T_longl_mat * GRM_Tmat);
                        tau_ind = tau_ind + 1;
                }
           }
        }

   }

    //return g_spSigma;
}



// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_multiV(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG, bool LOCO){
    // Start Timers
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    int Nnomissing = wVec.n_elem;
    arma::fvec xVec(Nnomissing);
    xVec.zeros();

    if(g_isStoreSigma){
        std::cout << " arma::spsolve(g_spSigma, bVec) 0" << std::endl;
        xVec = arma::spsolve(g_spSigma, bVec);
        std::cout << " arma::spsolve(g_spSigma, bVec) 1" << std::endl;
    }else{
        arma::fvec rVec = bVec;
        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
	//std::cout << "getPCG1ofSigmaAndVector_multiV Here1" << std::endl;
        minvVec = 1/getDiagOfSigma_multiV(wVec, tauVec, LOCO);
        zVec = minvVec % rVec;

        float sumr2 = sum(rVec % rVec);
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;

        int iter = 0;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
		iter = iter + 1;
	//	        std::cout << "getPCG1ofSigmaAndVector_multiV Here2" << std::endl;
                arma::fcolvec ApVec = getCrossprod_multiV(pVec, wVec, tauVec, LOCO);
	//	        std::cout << "getPCG1ofSigmaAndVector_multiV Here3" << std::endl;
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
                z1Vec = minvVec % r1Vec;

                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} 
        return(xVec);
}



// [[Rcpp::export]]
void set_seed(unsigned int seed) {
        Rcpp::Environment base_env("package:base");
        Rcpp::Function set_seed_r = base_env["set.seed"];
        set_seed_r(seed);
}

// [[Rcpp::export]]
Rcpp::NumericVector nb(int n) {
        return(rbinom(n,1,0.5));
}

// [[Rcpp::export]]
void setStartEndIndex(int startIndex, int endIndex, int chromIndex){
  ptr_gNULLGENOobj->startIndex = startIndex;
  ptr_gNULLGENOobj->endIndex = endIndex;
  ptr_gNULLGENOobj->Msub = 0;
  ptr_gNULLGENOobj->chromIndex = chromIndex;

  for(size_t i=0; i< ptr_gNULLGENOobj->M; i++){
        if(i < startIndex || i > endIndex){
                if(ptr_gNULLGENOobj->alleleFreqVec[i] >= ptr_gNULLGENOobj->minMAFtoConstructGRM && ptr_gNULLGENOobj->alleleFreqVec[i] <= 1-ptr_gNULLGENOobj->minMAFtoConstructGRM){

                        ptr_gNULLGENOobj->Msub = ptr_gNULLGENOobj->Msub + 1;
                }
        }
  }
}

// [[Rcpp::export]]
void setStartEndIndexVec( arma::ivec & startIndex_vec,  arma::ivec & endIndex_vec){

  startIndex_vec.print("startIndex_vec");

  ptr_gNULLGENOobj->setstartendIndexVec(startIndex_vec, endIndex_vec);
  //ptr_gNULLGENOobj->startIndexVec = startIndex_vec;
  //ptr_gNULLGENOobj->endIndexVec = endIndex_vec;
  //ptr_gNULLGENOobj->Msub = ptr_gNULLGENOobj->M - (endIndex - startIndex + 1);
}

//This function calculates the coefficients of variation for mean of a vector
// [[Rcpp::export]]
float calCV(arma::fvec& xVec){
  int veclen = xVec.n_elem;
  float vecMean = arma::mean(xVec);
  float vecSd = arma::stddev(xVec);
  float vecCV = (vecSd/vecMean)/veclen;
  return(vecCV);
}

// [[Rcpp::export]]
arma::fmat getSigma_X_multiV(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, int maxiterPCG, float tolPCG, bool LOCO){


        int Nnomissing = Xmat.n_rows;
        int colNumX = Xmat.n_cols;
        arma::fmat Sigma_iX1(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;

        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);
                Sigma_iX1.col(i) = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG, LOCO);
        }
        return(Sigma_iX1);
}

// [[Rcpp::export]]
arma::fvec  getSigma_G_multiV(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG, bool LOCO){
        arma::fvec Sigma_iG;
        Sigma_iG = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, Gvec, maxiterPCG, tolPCG, LOCO);
        return(Sigma_iG);
}


// [[Rcpp::export]]
Rcpp::List fitglmmaiRPCG_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec &wVec,  arma::fvec & tauVec, arma::ivec & fixtauVec, arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff, bool LOCO){
        Function warning("warning");

        unsigned int k1 = g_num_Kmat;

        int q2 = arma::sum(fixtauVec==0);
        arma::uvec idxtau = arma::find(fixtauVec==0);
        arma::fvec tau0 = tauVec;

        std::cout << "check 1" << std::endl;
        Rcpp::List re = getAIScore_multiV(Yvec, Xmat,wVec,  tauVec, fixtauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO);

        std::cout << "check 2" << std::endl;

        arma::fvec YPAPY = re["YPAPY"];
        arma::fvec Trace = re["Trace"];
        arma::fvec score1 = YPAPY - Trace;
        score1.print("score1");
        YPAPY.print("YPAPY");
        Trace.print("Trace");
        arma::fmat AI1 = re["AI"];
        arma::fvec Dtau(AI1.n_cols);
        score1.print("score");
        try{
                Dtau = arma::solve(AI1, score1, arma::solve_opts::allow_ugly);
        }
        catch(std::runtime_error){
                std::cout << "arma::solve(AI, score): AI seems singular, using less variant components matrix is suggested." << std::endl;
                Dtau.zeros();
        }
        //std::cout << "check 3" << std::endl;
        //Dtau.print("Dtau");
        //score1.print("score1");
        //AI1.print("AI1");
        //
        //
        arma::fvec Dtau_k1(k1);
        Dtau_k1.zeros();

        //std::cout << "k1 " << k1 << std::endl;
        //fixtauVec.print("fixtauVec");

        // fill dtau using dtau_pre, padding 0
        int i2 = 0;
        for(int i=0; i<k1; i++){
                std::cout << "i " << i << std::endl;
                if(fixtauVec(i)==0){ // not fixed
                        Dtau_k1(i) = Dtau(i2);
                        i2++;
                }
        } // end for i
       // Dtau_k1.print("Dtau_k1");
        tau0 = tauVec;
        tauVec = tauVec + Dtau_k1;
        arma::fvec tauVecabs, tauVec_sub1;
        arma::ivec fixrhoidx0, fixrhoidx, tauupdateidx;
        //arma::ivec fixrhoidx0, fixrhoidx, covarianceidxVec1, covarianceidxVec_sub1, covarianceidxVec2, covarianceidxVec_sub2,covarianceidxVec3, covarianceidxVec_sub3, tauupdateidx;
        arma::uvec covarianceidxVec1, covarianceidxVec_sub1, covarianceidxVec2, covarianceidxVec_sub2,covarianceidxVec3, covarianceidxVec_sub3, zeroindVec1;

        if(g_covarianceidxMat.n_rows == 0){
                tauVec.elem( arma::find(tauVec < tol && tau0 < tol) ).zeros();
                float step = 1.0;
                while ( arma::any(tauVec < 0.0) ){
                        step = step*0.5;
                        tauVec = tau0 + step*Dtau_k1;
                        tauVec.elem( arma::find(tauVec < tol && tau0 < tol) ).zeros();
                } // end while
                tauVec.elem( arma::find(tauVec < tol) ).zeros();
        }else{
                fixrhoidx0 = updatefixrhoidx0(tau0, tol);
                tauVec.elem( arma::find(tauVec < tol && tau0 < tol) ).zeros();
                tauVec.elem(g_covarianceidxMat_col1) = tauVec.elem(g_covarianceidxMat_col1);
                fixrhoidx = updatefixrhoidx0(tauVec, tol);
                covarianceidxVec_sub1 = g_covarianceidxMat_col1.elem(arma::find(fixrhoidx0 == 1 && fixrhoidx == 1));
                covarianceidxVec_sub2 = g_covarianceidxMat_col2.elem(arma::find(fixrhoidx0 == 1 && fixrhoidx == 1));
                covarianceidxVec_sub3 = g_covarianceidxMat_col3.elem(arma::find(fixrhoidx0 == 1 && fixrhoidx == 1));

		tauVec_sub1 = tauVec.elem(covarianceidxVec_sub1);
                zeroindVec1 = arma::find(tauVec_sub1 != 0);
		if(zeroindVec1.n_elem == tauVec_sub1.n_elem){
                        tauVecabs = tauVec.elem(covarianceidxVec_sub1) / arma::abs(tauVec.elem(covarianceidxVec_sub1));
                        tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
                }else{
                        tauVecabs = tauVec_sub1;
                        tauVecabs.zeros();
                        tauVecabs.elem(zeroindVec1) = tauVecabs.elem(zeroindVec1)/arma::abs(tauVecabs.elem(zeroindVec1));
                        tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
                }



                //tauVecabs = tauVec.elem(covarianceidxVec_sub1) / arma::abs(tauVec.elem(covarianceidxVec_sub1));
                //tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
                tauupdateidx = tauUpdateValue(tauVec);
                float step = 1.0;
                 while ( arma::any(tauVec.elem(g_covarianceidxMat_notcol1) < 0.0) || tauVec(0) < 0.0 || arma::any(tauupdateidx == 0)){
                 tauVec.print("tauVec");
                 tauupdateidx.print("tauupdateidx");
                        step = step*0.5;
                        tauVec = tau0 + step*Dtau_k1;
                        tauVec.elem( arma::find(tauVec < tol && tau0 < tol) ).zeros();
                        tauVec.elem(g_covarianceidxMat_col1) = tauVec.elem(g_covarianceidxMat_col1);
                        fixrhoidx = updatefixrhoidx0(tauVec, tol);
                        covarianceidxVec_sub1 = g_covarianceidxMat_col1.elem(arma::find(fixrhoidx0 == 1 && fixrhoidx == 1));
                        covarianceidxVec_sub2 = g_covarianceidxMat_col2.elem(arma::find(fixrhoidx0 == 1 && fixrhoidx == 1));
                        covarianceidxVec_sub3 = g_covarianceidxMat_col3.elem(arma::find(fixrhoidx0 == 1 && fixrhoidx == 1));
                        tauVecabs = tauVec.elem(covarianceidxVec_sub1) / arma::abs(tauVec.elem(covarianceidxVec_sub1));
                        tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
                        tauupdateidx = tauUpdateValue(tauVec);
                        Rcpp::checkUserInterrupt();
                } // end while
                tauVec.elem( arma::find(tauVec < tol) ).zeros();
                tauVec.elem(g_covarianceidxMat_col1) = tauVec.elem(g_covarianceidxMat_col1);

		tauVec_sub1 = tauVec.elem(covarianceidxVec_sub1);
		zeroindVec1 = arma::find(tauVec_sub1 != 0);	

		if(zeroindVec1.n_elem == tauVec_sub1.n_elem){
                	tauVecabs = tauVec.elem(covarianceidxVec_sub1) / arma::abs(tauVec.elem(covarianceidxVec_sub1));
                	tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
		}else{
			tauVecabs = tauVec_sub1;
			tauVecabs.zeros();
			tauVecabs.elem(zeroindVec1) = tauVecabs.elem(zeroindVec1)/arma::abs(tauVecabs.elem(zeroindVec1));
			tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));			
		}	
        }
         return List::create(Named("tau") = tauVec, Named("AI") = AI1, Named("score") = score1);
}



// [[Rcpp::export]]
arma::fvec getMeanDiagofKmat(bool LOCO){

		                  std::cout << "Here3a" << std::endl;	
	arma::fvec mean_diag_kins_vec(g_num_Kmat - 1);
		                  std::cout << "Here3b" << std::endl;	

        arma::sp_vec diagVecG0;
        arma::sp_fvec diagVecV0;
        arma::fvec diagVecG, diagVecV, diagVecG_I, diagVecG_T, diagVecG_IT,diagVecV_I, diagVecV_T, diagVecV_IT;
        arma::fvec diagVec;
        unsigned int tauind = 0;
        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){

        if(g_isGRM){
           if(!LOCO){
              if(!g_isSparseGRM){
                int MminMAF = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
                diagVec = (*ptr_gNULLGENOobj->Get_Diagof_StdGeno()) /MminMAF;
              }else{
                diagVec = arma::diagvec(g_spGRM);
              }
            }else{
              diagVec = (*ptr_gNULLGENOobj->Get_Diagof_StdGeno_LOCO());
              int Msub_MAFge_minMAFtoConstructGRM_in_b = ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_in();
              int Msub_MAFge_minMAFtoConstructGRM_singleVar_b = ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_singleChr_in();
              diagVec = diagVec/(Msub_MAFge_minMAFtoConstructGRM_in_b - Msub_MAFge_minMAFtoConstructGRM_singleVar_b);
            }
	   mean_diag_kins_vec(tauind) = arma::mean(diagVec);
	   tauind = tauind + 1;
	  }

          if(Kmat_vec.size() > 0){
            for(unsigned int i = 0; i < Kmat_vec.size(); i++){
              diagVec = (Kmat_vec[i]).diag();
              mean_diag_kins_vec(i+tauind) = arma::mean(diagVec);
	      tauind = tauind + 1;
            }
          }

        }else{ //if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
	    /*if(g_isGRM && g_isSparseGRM){	
                diagVecG = arma::diagvec(g_spGRM);
                diagVecG_I = diagVecG.elem(g_I_longl_vec);
                diagVec = diagVecG_I;
                mean_diag_kins_vec(0) = arma::mean(diagVec);
                tauind = tauind + 1;
            }*/
		                  std::cout << "Here3c" << std::endl;	
		std::cout << "g_T_longl_mat.n_rows " << g_T_longl_mat.n_rows << std::endl;
            if(g_T_longl_mat.n_rows > 0){
		std::cout << "g_isGRM " << g_isGRM << std::endl;
		std::cout << "g_isSparseGRM " << g_isSparseGRM << std::endl;

	        if(g_isGRM && g_isSparseGRM){
			//g_spGRM.print("g_spGRM");		
		  diagVecG = arma::diagvec(g_spGRM);

                  diagVecG_I = diagVecG.elem(g_I_longl_vec);
                  diagVec = diagVecG_I;
                  mean_diag_kins_vec(0) = arma::mean(diagVec);
                  tauind = tauind + 1;

                  diagVecG_IT = diagVecG_I % g_T_longl_vec;
                  diagVecG_T = diagVecG_IT % g_T_longl_vec;
                  diagVecG_IT = 2 * diagVecG_IT;
                  diagVec = diagVecG_IT;
                  std::cout << "Here1" << std::endl;
                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  tauind = tauind + 1;
                  std::cout << "Here2" << std::endl;
                  diagVec = diagVecG_T;
		  std::cout << "tauind " << tauind << std::endl;
                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  std::cout << "Here2" << std::endl;
                  tauind = tauind + 1;
		}  

                  //diagVecV = diagVecG;
		std::cout << "g_n_unique " << g_n_unique << std::endl;
                  diagVecV.ones(g_n_unique);
		  //diagVecV.print("diagVecV");
                  diagVecV_I = diagVecV.elem(g_I_longl_vec);
                  diagVec = diagVecV_I;
                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
		  std::cout << "Here2a" << std::endl;
                  tauind = tauind + 1;
                  diagVecV_IT = diagVecV_I % g_T_longl_vec;
                  diagVecV_T = diagVecV_IT % g_T_longl_vec;
                  diagVecV_IT = 2 * diagVecV_IT;
                  diagVec = diagVecV_IT;
		  std::cout << "Here2b" << std::endl;

                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  tauind = tauind + 1;
                  std::cout << "Here2" << std::endl;
                  diagVec = diagVecV_T;
                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  tauind = tauind + 1;


           }else{
		                  std::cout << "Here3d" << std::endl;	
		if(g_isGRM && g_isSparseGRM){	
		                  std::cout << "Here3e" << std::endl;	
                  diagVecG = arma::diagvec(g_spGRM);
                  diagVecG_I = diagVecG.elem(g_I_longl_vec);
                  diagVec = diagVecG_I;
                  mean_diag_kins_vec(0) = arma::mean(diagVec);
                  tauind = tauind + 1;
	  	}
		                  std::cout << "Here3" << std::endl;	
		  diagVecV = diagVecG;
                  diagVecV.ones();
                  diagVecV_I = diagVecV.elem(g_I_longl_vec);
                  diagVec = diagVecV_I;
                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  tauind = tauind + 1;		
	   }	   

                if(Kmat_vec.size() > 0){
                  for(unsigned int i = 0; i < Kmat_vec.size(); i++){
                    diagVecV0 = Kmat_vec[i].diag();
                    arma::fvec diagVecVtemp(diagVecV0);
                    diagVecV = diagVecVtemp;
                    diagVecV_I = diagVecV.elem(g_I_longl_vec);
                    diagVec = diagVecV_I;
                    mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                    tauind = tauind + 1;
                    if(g_T_longl_mat.n_rows > 0){
                      diagVecV_IT = diagVecV_I % g_T_longl_vec;
                      diagVecV_T = diagVecV_IT % g_T_longl_vec;
                      diagVecV_IT = 2 * diagVecV_IT;
                      diagVec = diagVecV_IT;
                      mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                      tauind = tauind + 1;
                      diagVec = diagVecV_T;
                      mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                      tauind = tauind + 1;
                    }
                  }
                }

        }

        return(mean_diag_kins_vec);
}


arma::ivec updatefixrhoidx0(arma::fvec & t_tau0Vec, float tol){
        arma::ivec fixrhoidx0Vec(g_covarianceidxMat.n_rows);
        arma::uvec covarianceidxVec;
        float tau0_x1, tau0_x2, tau0_x3, tau0_x2_x3, tau0temp;
        for(int i=0; i<g_covarianceidxMat.n_rows; i++){
                covarianceidxVec = (g_covarianceidxMat.row(i) - 1).t();
                tau0_x1 = t_tau0Vec((covarianceidxVec(0)));
                tau0_x2 = t_tau0Vec((covarianceidxVec(1)));
                tau0_x3 = t_tau0Vec((covarianceidxVec(2)));
                tau0_x2_x3 = std::sqrt(tau0_x2 * tau0_x3);
                tau0temp = (1 - 1.01 * tol) * tau0_x2_x3;
                if(std::abs(tau0_x1) > tau0temp){
                        fixrhoidx0Vec(i) = 1;
                }else{
                        fixrhoidx0Vec(i) = 0;
                }
        }
        return(fixrhoidx0Vec);
}


arma::ivec tauUpdateValue(arma::fvec & t_tau0Vec){
        arma::ivec fixrhoidx0Vec(g_covarianceidxMat.n_rows);
        arma::uvec covarianceidxVec;
        float tau0_x1, tau0_x2, tau0_x3, tau0_x2_x3, tau0temp;
        for(int i=0; i<g_covarianceidxMat.n_rows; i++){
                covarianceidxVec = (g_covarianceidxMat.row(i) - 1).t();
                tau0_x1 = t_tau0Vec((covarianceidxVec(0)));
                tau0_x2 = t_tau0Vec((covarianceidxVec(1)));
                tau0_x3 = t_tau0Vec((covarianceidxVec(2)));
                tau0_x2_x3 = std::sqrt(tau0_x2 * tau0_x3);
                if(std::abs(tau0_x1) > tau0_x2_x3){
                        fixrhoidx0Vec(i) = 0;
                }else{
                        fixrhoidx0Vec(i) = 1;
                }
        }
        return(fixrhoidx0Vec);
}


// [[Rcpp::export]]
Rcpp::List getAIScore_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, arma::ivec & fixtauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO){
	fixtauVec.print("fixtauVec");

        int q2 = arma::sum(fixtauVec==0);
        arma::uvec idxtau = arma::find(fixtauVec==0);
        arma::fvec tau0;
        unsigned int k1 = g_num_Kmat;
        arma::fmat AI(k1,k1);
        arma::fvec YPAPY(k1);
        YPAPY.zeros();
        arma::fvec Trace(k1);
        Trace.zeros();
        //std::cout << "k1 " << k1 << std::endl;
        arma::fmat Sigma_iXt = Sigma_iX.t();
        arma::fmat Xmatt = Xmat.t();

        //Sigma_iY.print("Sigma_iY");
        //Sigma_iX.print("Sigma_iX");
        //cov.print("cov");
        //Sigma_iXt.print("Sigma_iXt");
        //Yvec.print("Yvec");

        arma::fvec PY1 = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Yvec));

        //PY1.print("PY1");
        //arma::fvec APY = getCrossprodMatAndKin(PY1);
        //float YPAPY = dot(PY1, APY);
        //arma::fvec A0PY = PY1; ////Quantitative
        //float YPA0PY = dot(PY1, A0PY); ////Quantitative
        //arma::fvec Trace = GetTrace_q(Sigma_iX, Xmat, wVec, tauVec, cov1, nrun, maxiterPCG, tolPCG, traceCVcutoff);
        unsigned int n = PY1.n_elem;
        //arma::fvec PA0PY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, A0PY, maxiterPCG, tolPCG);
        //arma::fvec PA0PY = PA0PY_1 - Sigma_iX * (cov1 * (Sigma_iXt * PA0PY_1));
        arma::fvec PAPY_1, PAPY,  APY;
        arma::fmat APYmat(n, k1);

        //if(q2>0){
        //
	arma::fvec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;


        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
                for(int i=0; i<k1; i++){
                    if(fixtauVec(i) == 0){
                        if(i==0){
                                APY = PY1;
                        }else if(i==1){
				if(g_isGRM){
                                        APY = getCrossprodMatAndKin(PY1, LOCO);
                                }else{
					if(Kmat_vec.size() > 0){
                                        	APY = Kmat_vec[0]*PY1;
                                 	}			
				}	
                        }else{
				if(g_isGRM){
                                	if(Kmat_vec.size() > 0){
                                        	APY = Kmat_vec[i-2]*PY1;
                                	}
				}else{
					if(Kmat_vec.size() > 0){
                                                APY = Kmat_vec[i-1]*PY1;
                                        }	
				}	
                        }
                        APYmat.col(i) = APY;
                        PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG, LOCO);
                        PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
                        for(int j=0; j<=i; j++){
                                AI(i,j) = arma::dot(APYmat.col(j), PAPY);
                                if(j !=i){
                                        AI(j,i) = AI(i,j);
                                }
                        }
                        YPAPY(i) = arma::dot(PY1, APYmat.col(i));
                     }
                }
        }else{
               Ibvec = g_I_longl_mat.t() * PY1;
                std::cout << "Ibvec.n_elem " << Ibvec.n_elem << std::endl;

                if(g_isGRM){
                        GRM_I_bvec = getCrossprodMatAndKin(Ibvec, LOCO);
                }

                std::cout << "g_I_longl_mat.n_rows " << g_I_longl_mat.n_rows << std::endl;
               if(g_T_longl_mat.n_rows == 0){
                  for(int i=0; i<k1; i++){
                    if(fixtauVec(i) == 0){
                        if(i==0){
                                APY = PY1;
                        }else if(i==1){
				if(g_isGRM){
                                  std::cout << "GRM_I_bvec.n_elem " << GRM_I_bvec.n_elem << std::endl;
                                  APY = g_I_longl_mat * GRM_I_bvec;
				}else{
				  APY = g_I_longl_mat * Ibvec;
				}	
                        }else if (i == 2){
				if(g_isGRM){
                                  APY = g_I_longl_mat * Ibvec;
				}else{
				  if(Kmat_vec.size() > 0){
					V_I_bvec = Kmat_vec[0]*Ibvec;
					APY = g_I_longl_mat * V_I_bvec;
				  }	  
				}	
                                //APY = Kmat_vec[i-2]*PY1;
                        }else{
			     if(g_isGRM){	
                                if(Kmat_vec.size() > 0){
                                        V_I_bvec = Kmat_vec[i-3]*Ibvec;
                                        APY = g_I_longl_mat * V_I_bvec;
                                }
			     }else{
				if(Kmat_vec.size() > 0){
                                        V_I_bvec = Kmat_vec[i-2]*Ibvec;
                                        APY = g_I_longl_mat * V_I_bvec;
                                }
			     }		     

                        }
                        APYmat.col(i) = APY;
			PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG, LOCO);
                        PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
                        for(int j=0; j<=i; j++){
                                AI(i,j) = arma::dot(APYmat.col(j), PAPY);
                                if(j !=i){
                                        AI(j,i) = AI(i,j);
                                }
                        }
                        YPAPY(i) = arma::dot(PY1, APYmat.col(i));
                     }
                  }

                }else{
                  Tbvec = g_T_longl_mat.t() * PY1;
                  if(g_isGRM){
                    GRM_T_bvec = getCrossprodMatAndKin(Tbvec, LOCO);
                  }

                  unsigned int kmatind = 0;
                  for(int i=0; i<k1; i++){
                    if(fixtauVec(i) == 0){
		     if(g_isGRM){	    
                        if(i==0){
                                APY = PY1;
                        }else if(i==1){
                                APY = g_I_longl_mat * GRM_I_bvec;
                        }else if(i == 2){
                                APY = (g_T_longl_mat * GRM_I_bvec) + (g_I_longl_mat * GRM_T_bvec);
                        }else if(i == 3){
                                APY = g_T_longl_mat * GRM_T_bvec;

                        }else if(i == 4){
                                APY = g_I_longl_mat * Ibvec;
                        }else if(i == 5){

                                APY = (g_T_longl_mat * Ibvec) + (g_I_longl_mat * Tbvec);

                        }else if(i == 6){
                                APY = g_T_longl_mat * Tbvec;

                        }else{
                                if(i % 3 == 1){
                                        V_I_bvec = Kmat_vec[kmatind] * Ibvec;
                                        APY = g_I_longl_mat * V_I_bvec;
                                }else if(i % 3 == 2){
                                        V_T_bvec = Kmat_vec[kmatind] * Tbvec;
                                        APY = (g_T_longl_mat * V_I_bvec) + (g_I_longl_mat * V_T_bvec);
                                }else if(i % 3 == 0){
                                        APY = g_T_longl_mat * V_T_bvec;
                                        kmatind = kmatind + 1;

                                }
                        }
		     }else{ //if(g_isGRM){
			if(i==0){
                                APY = PY1;
                        }else if(i==1){
                                APY = g_I_longl_mat * Ibvec;
                        }else if(i == 2){
                                APY = (g_T_longl_mat * Ibvec) + (g_I_longl_mat * Tbvec);
                        }else if(i == 3){
                                APY = g_T_longl_mat * Tbvec;

                        }else{
                                if(i % 3 == 1){
                                        V_I_bvec = Kmat_vec[kmatind] * Ibvec;
                                        APY = g_I_longl_mat * V_I_bvec;
                                }else if(i % 3 == 2){
                                        V_T_bvec = Kmat_vec[kmatind] * Tbvec;
                                        APY = (g_T_longl_mat * V_I_bvec) + (g_I_longl_mat * V_T_bvec);
                                }else if(i % 3 == 0){
                                        APY = g_T_longl_mat * V_T_bvec;
                                        kmatind = kmatind + 1;

                                }
                        }     



		     }//else for if(g_isGRM){	     
                        APYmat.col(i) = APY;
                        PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG, LOCO);
                        PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
                        for(int j=0; j<=i; j++){
                                AI(i,j) = arma::dot(APYmat.col(j), PAPY);
                                if(j !=i){
                                        AI(j,i) = AI(i,j);
                                }
                        }
			                        YPAPY(i) = arma::dot(PY1, APYmat.col(i));
                     }
                  }

                }
        }

        AI.print("AI");
	YPAPY.print("YPAPY");
        arma::fmat AI_update = AI.submat(idxtau, idxtau);
        arma::fvec YPAPY_update = YPAPY.elem(idxtau);
	YPAPY.print("YPAPY");

        //vector with length=q2
        Trace = GetTrace_multiV(Sigma_iX, Xmat, wVec, tauVec, fixtauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO);
        //YPAPY_update.print("YPAPY_update");
        Trace.print("Trace");
        //arma::fvec PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG);
        //arma::fvec PAPY = PAPY_1 - Sigma_iX * (cov1 * (Sigma_iXt * PAPY_1));
        return Rcpp::List::create(Named("YPAPY") = YPAPY_update, Named("Trace") = Trace,Named("PY") = PY1,Named("AI") = AI_update);
}



// [[Rcpp::export]]
arma::fvec GetTrace_multiV(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::ivec & fixtauVec, arma::fmat& cov1,  int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO){

        set_seed(200);

        int q2 = arma::sum(fixtauVec==0);
        arma::uvec idxtau = arma::find(fixtauVec==0);

        idxtau.print("idxtau");

        arma::fmat Sigma_iXt = Sigma_iX.t();
        int Nnomissing = wVec.n_elem;;
        //unsigned int k = Kmat_vec.size();
        //unsigned int k1 = k + 2;
        unsigned int  k1 = g_num_Kmat;
        arma::fmat temp_mat(nrun, k1);
        arma::fmat temp_mat_update(nrun, q2);
        arma::fvec temp_vec(nrun);
        arma::fvec temp_vec_double(Nnomissing);
        temp_mat.zeros();
        temp_mat_update.zeros();

        arma::fvec Sigma_iu;
        arma::fvec Pu;
        arma::fmat Au_mat(Nnomissing, k1);
        arma::fvec uVec;
	Rcpp::NumericVector uVec0;

        int nrun_trace_start = 0;
        int nrun_trace_end = nrun;
        arma::fvec traceCV(q2);
        traceCV.fill(traceCVcutoff + 0.1);

        arma::uvec covarianceidxVec;
        arma::fvec traceCVsub;
        arma::uvec indexsubvec =  { 1, 2 };
        if(g_covarianceidxMat.n_cols > 0){
             covarianceidxVec = arma::vectorise(g_covarianceidxMat.cols(indexsubvec));
             covarianceidxVec = covarianceidxVec - 1;
        }
        bool isConverge = false;
        //while((traceCV > cutoff_trace) | (traceCV0 > cutoff_trace)){
        //while( arma::any(traceCV > traceCVcutoff) ){
        //
        //
	arma::fvec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;
        while( !isConverge ){
                for(int i = nrun_trace_start; i < nrun_trace_end; i++){

                        uVec0 = nb(Nnomissing);
                        uVec = as<arma::fvec>(uVec0);
                        uVec = uVec*2 - 1;

			//std::cout << "GetTrace_multiV Here 1" << std::endl;
                        Sigma_iu = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, uVec, maxiterPCG, tolPCG, LOCO);
                        Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));

                        if(fixtauVec(0) == 0)   {
                                Au_mat.col(0) = uVec;
                                temp_mat(i,0) = dot(Au_mat.col(0), Pu);
                        }
                        // conversion for ops with sp_mat
                     if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
                        if(fixtauVec(1) == 0){
                                if(g_isGRM){
                                        temp_vec_double = getCrossprodMatAndKin(uVec, LOCO);
                                }else{
					if(Kmat_vec.size() > 0){
						temp_vec_double = 0.0+Kmat_vec[0] * uVec;
					}	
                                }
                                Au_mat.col(1) = temp_vec_double;
                                temp_mat(i,1) = dot(temp_vec_double, Pu);
                        }


			if(g_isGRM){
                          for(int j=2; j<k1;j++){
                                if(fixtauVec(j) == 0){
                                        Au_mat.col(j) = 0.0+Kmat_vec[j-2] * uVec;
                                        temp_mat(i,j) = dot(Au_mat.col(j), Pu);
                                }
                          } // end for j in 2:k1
			}else{
			  for(int j=2; j<k1;j++){
                                if(fixtauVec(j) == 0){
                                        Au_mat.col(j) = 0.0+Kmat_vec[j-1] * uVec;
                                        temp_mat(i,j) = dot(Au_mat.col(j), Pu);
                                }
                          } // end for j in 1:k1	

			}	
                   }else{
                        Ibvec = g_I_longl_mat.t() * uVec;
			//std::cout << "GetTrace_multiV Here 2" << std::endl;
                        if(g_isGRM){
                                GRM_I_bvec = getCrossprodMatAndKin(Ibvec, LOCO);
                        }else{
				GRM_I_bvec = Ibvec;
			}	
                        if(g_T_longl_mat.n_rows == 0){
			   if(g_isGRM){	
                                if(fixtauVec(1) == 0)   {
                                        temp_vec_double = g_I_longl_mat * GRM_I_bvec;
                                        Au_mat.col(1) = temp_vec_double;
                                        temp_mat(i,1) = dot(temp_vec_double, Pu);
                                }

                                if(fixtauVec(2) == 0){
                                        temp_vec_double = g_I_longl_mat * Ibvec;
                                        Au_mat.col(2) = temp_vec_double;
                                        temp_mat(i,2) = dot(temp_vec_double, Pu);
                                }

                                for(int j=3; j<k1;j++){
                                        if(fixtauVec(j) == 0){
                                                V_I_bvec = Kmat_vec[j-3] * Ibvec;
                                                temp_vec_double = g_I_longl_mat * V_I_bvec;
                                                Au_mat.col(j) = temp_vec_double;
                                                temp_mat(i,j) = dot(temp_vec_double, Pu);
                                        }
                                }
			   }else{

                                if(fixtauVec(1) == 0){
                                        temp_vec_double = g_I_longl_mat * Ibvec;
                                        Au_mat.col(1) = temp_vec_double;
                                        temp_mat(i,1) = dot(temp_vec_double, Pu);
                                }

                                for(int j=2; j<k1;j++){
                                        if(fixtauVec(j) == 0){
                                                V_I_bvec = Kmat_vec[j-2] * Ibvec;
                                                temp_vec_double = g_I_longl_mat * V_I_bvec;
                                                Au_mat.col(j) = temp_vec_double;
                                                temp_mat(i,j) = dot(temp_vec_double, Pu);
                                        }
                                }		
			   }	   
                        }else{
                             Tbvec = g_T_longl_mat.t() * uVec;
                             if(g_isGRM){
                                        GRM_T_bvec = getCrossprodMatAndKin(Tbvec, LOCO);
                                if(fixtauVec(1) == 0)   {
                                        temp_vec_double = g_I_longl_mat * GRM_I_bvec;
                                        Au_mat.col(1) = temp_vec_double;
                                        temp_mat(i,1) = dot(temp_vec_double, Pu);
                                }

                                if(fixtauVec(2) == 0)   {
                                        temp_vec_double = (g_I_longl_mat * GRM_T_bvec) + (g_T_longl_mat * GRM_I_bvec);
                                        Au_mat.col(2) = temp_vec_double;
                                        temp_mat(i,2) = dot(temp_vec_double, Pu);
                                }

                                if(fixtauVec(3) == 0)   {
                                        temp_vec_double = g_T_longl_mat * GRM_T_bvec;
                                        Au_mat.col(3) = temp_vec_double;
                                        temp_mat(i,3) = dot(temp_vec_double, Pu);
                                }

                                 if(fixtauVec(4) == 0)   {
                                        temp_vec_double = g_I_longl_mat * Ibvec;
                                        Au_mat.col(4) = temp_vec_double;
                                        temp_mat(i,4) = dot(temp_vec_double, Pu);
                                }

                                if(fixtauVec(5) == 0)   {
                                        temp_vec_double = (g_I_longl_mat * Tbvec) + (g_T_longl_mat * Ibvec);
                                        Au_mat.col(5) = temp_vec_double;
                                        temp_mat(i,5) = dot(temp_vec_double, Pu);
                                }

                                if(fixtauVec(6) == 0)   {
                                        temp_vec_double = g_T_longl_mat * Tbvec;
                                        Au_mat.col(6) = temp_vec_double;
                                        temp_mat(i,6) = dot(temp_vec_double, Pu);
                                }


                                int j = 7;
                                while(j < k1){
                                        V_I_bvec = Kmat_vec[j-7] * Ibvec;
                                        V_T_bvec = Kmat_vec[j-7] * Tbvec;
                                //for(int j=7; j<k1;j++){
                                        if(fixtauVec(j) == 0){
                                                temp_vec_double = g_I_longl_mat * V_I_bvec;
                                                Au_mat.col(j) = temp_vec_double;
                                                temp_mat(i,j) = dot(temp_vec_double, Pu);
                                        }
                                        j = j + 1;
                                        if(fixtauVec(j) == 0){
                                                temp_vec_double = (g_I_longl_mat * V_T_bvec) + (g_T_longl_mat * V_I_bvec);
                                                Au_mat.col(j) = temp_vec_double;
                                                temp_mat(i,j) = dot(temp_vec_double, Pu);
                                        }
                                        j = j + 1;

                                        if(fixtauVec(j) == 0){
                                                temp_vec_double = g_T_longl_mat * V_T_bvec;
                                                Au_mat.col(j) = temp_vec_double;
                                                temp_mat(i,j) = dot(temp_vec_double, Pu);

                                        }
                                        j = j + 1;
                                }

                        }else{// if(g_isGRM){


                                 if(fixtauVec(1) == 0)   {
                                        temp_vec_double = g_I_longl_mat * Ibvec;
                                        Au_mat.col(1) = temp_vec_double;
                                        temp_mat(i,1) = dot(temp_vec_double, Pu);
                                }

                                if(fixtauVec(2) == 0)   {
                                        temp_vec_double = (g_I_longl_mat * Tbvec) + (g_T_longl_mat * Ibvec);
                                        Au_mat.col(2) = temp_vec_double;
                                        temp_mat(i,2) = dot(temp_vec_double, Pu);
                                }

                                if(fixtauVec(3) == 0)   {
                                        temp_vec_double = g_T_longl_mat * Tbvec;
                                        Au_mat.col(3) = temp_vec_double;
                                        temp_mat(i,3) = dot(temp_vec_double, Pu);
                                }


                                int j = 4;
                                while(j < k1){
                                        V_I_bvec = Kmat_vec[j-4] * Ibvec;
                                        V_T_bvec = Kmat_vec[j-4] * Tbvec;
                                //for(int j=7; j<k1;j++){
                                        if(fixtauVec(j) == 0){
                                                temp_vec_double = g_I_longl_mat * V_I_bvec;
                                                Au_mat.col(j) = temp_vec_double;
                                                temp_mat(i,j) = dot(temp_vec_double, Pu);
                                        }
                                        j = j + 1;
                                        if(fixtauVec(j) == 0){
                                                temp_vec_double = (g_I_longl_mat * V_T_bvec) + (g_T_longl_mat * V_I_bvec);
                                                Au_mat.col(j) = temp_vec_double;
                                                temp_mat(i,j) = dot(temp_vec_double, Pu);
                                        }
                                        j = j + 1;

                                        if(fixtauVec(j) == 0){
                                                temp_vec_double = g_T_longl_mat * V_T_bvec;
                                                Au_mat.col(j) = temp_vec_double;
                                                temp_mat(i,j) = dot(temp_vec_double, Pu);

                                        }
                                        j = j + 1;
                                }



			}//}else{// if(g_isGRM){	
		}		

	}



                        uVec.clear();
                        Pu.clear();
                        Sigma_iu.clear();

                } // end for i
                temp_mat_update = temp_mat.cols(idxtau);

                std::cout << "dim temp_mat_update" << temp_mat_update.n_rows << " " << temp_mat_update.n_cols << std::endl;;
                // update trace cv vector
                for(int k=0; k<q2; k++){
                        temp_vec = temp_mat_update.col(k);
                        traceCV(k) = calCV(temp_vec);
                }
                //traceCV.print("traceCV");
                //traceCVcutoff = 1.0;
                // if not converge, increase nrun_trace and rerun
                //temp_mat.print("temp_mat");
                //

                if(g_covarianceidxMat.n_cols > 0){
                        traceCVsub = traceCV.elem(covarianceidxVec);
                        //if(arma::any(traceCVsub > traceCVcutoff)){
                        if(arma::any(traceCV > traceCVcutoff)){
                                isConverge = false;
                        }else{
                                isConverge = true;
                        }
                }else{
                        if( arma::any(traceCV > traceCVcutoff) ){
                                isConverge = false;
                        }else{
                                isConverge = true;
                        }
                }



                if( !isConverge){
                        nrun_trace_start = nrun_trace_end;
                        nrun_trace_end = nrun_trace_end + 10;
                        temp_mat.resize(nrun_trace_end,k1);
                        temp_mat_update.resize(nrun_trace_end,q2);
                        //std::cout << "arma::mean(temp_mat0): " << arma::mean(temp_mat0) << std::endl;
                        Rcout << "CV for trace random estimator using "<< nrun_trace_start << " runs is " << traceCV <<  " > " << traceCVcutoff << std::endl;
                        Rcout << "try " << nrun_trace_end << "runs" << std::endl;
                } // end if arma::any(traceCV > traceCVcutoff)

        } // end while  arma::any(traceCV > traceCVcutoff)
        Au_mat.clear();
        Pu.clear();
        Sigma_iu.clear();
        uVec.clear();
        temp_vec.clear();

        arma::fvec traVec(q2);
        for(int i=0; i<q2; i++){
                traVec(i) = arma::mean(temp_mat_update.col(i));
        }
        temp_mat.clear();
        temp_mat_update.clear();
        return(traVec);
}



// [[Rcpp::export]]
Rcpp::List getCoefficients_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG, bool LOCO){

        int Nnomissing = wVec.n_elem;
        arma::fvec Sigma_iY;
	std::cout << "before Sigma_iY" << std::endl;
	std::cout << "Yvec.n_elem " << Yvec.n_elem << std::endl;
        Sigma_iY = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, Yvec, maxiterPCG, tolPCG, LOCO);
        std::cout << "after Sigma_iY" << std::endl;
        int colNumX = Xmat.n_cols;
        arma::fmat Sigma_iX(Nnomissing,colNumX);
        arma::fvec XmatVecTemp;
        for(int i = 0; i < colNumX; i++){
                XmatVecTemp = Xmat.col(i);
                Sigma_iX.col(i) = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG, LOCO);

        }
        arma::fmat Xmatt = Xmat.t();
        //arma::fmat cov = inv_sympd(Xmatt * Sigma_iX);
        arma::fmat cov;
        try {
          cov = arma::inv_sympd(arma::symmatu(Xmatt * Sigma_iX));
        } catch (const std::exception& e) {
          cov = arma::pinv(arma::symmatu(Xmatt * Sigma_iX));
          cout << "inv_sympd failed, inverted with pinv" << endl;
        }
	cov.print("cov"); 	
        arma::fmat Sigma_iXt = Sigma_iX.t();
	//Sigma_iXt.print("Sigma_iXt");
	//Yvec.print("Yvec");
        arma::fvec SigmaiXtY = Sigma_iXt * Yvec;
        arma::fvec alpha = cov * SigmaiXtY;
	
        arma::fvec eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
        return Rcpp::List::create(Named("Sigma_iY") = Sigma_iY, Named("Sigma_iX") = Sigma_iX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
}


// [[Rcpp::export]]
void setminMAC_VarianceRatio(float t_minMACVarRatio, float t_maxMACVarRatio, bool t_isVarianceRatioinGeno){
        ptr_gNULLGENOobj->g_minMACVarRatio = t_minMACVarRatio;
        ptr_gNULLGENOobj->g_maxMACVarRatio = t_maxMACVarRatio;
        ptr_gNULLGENOobj->isVarRatio = t_isVarianceRatioinGeno;
        std::cout << "ptr_gNULLGENOobj->g_minMACVarRatio " << ptr_gNULLGENOobj->g_minMACVarRatio << std::endl;
        std::cout << "ptr_gNULLGENOobj->g_maxMACVarRatio " << ptr_gNULLGENOobj->g_maxMACVarRatio << std::endl;
}


// [[Rcpp::export]]
arma::fvec get_GRMdiagVec(){
  int mMarker = gettotalMarker();
  int MminMAF = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        //cout << "MminMAF=" << MminMAF << endl;

  arma::fvec diagGRMVec = (*ptr_gNULLGENOobj->Get_Diagof_StdGeno())/MminMAF;
  return(diagGRMVec);
}

// [[Rcpp::export]]
void setminMAFforGRM(float minMAFforGRM){
  ptr_gNULLGENOobj->minMAFtoConstructGRM = minMAFforGRM;
}

// [[Rcpp::export]]
void setmaxMissingRateforGRM(float maxMissingforGRM){
  ptr_gNULLGENOobj->maxMissingRate = maxMissingforGRM;
}


// [[Rcpp::export]]
void set_Diagof_StdGeno_LOCO(){


  int Nnomissing = ptr_gNULLGENOobj->getNnomissing();
  int chrlength = ptr_gNULLGENOobj->startIndexVec.n_elem;
  (ptr_gNULLGENOobj->mtx_DiagStd_LOCO).zeros(Nnomissing, chrlength);
  (ptr_gNULLGENOobj->Msub_MAFge_minMAFtoConstructGRM_byChr).zeros(chrlength);
//  std::cout << "debug1" << std::endl;
    int starti, endi;
    arma::fvec * temp = &ptr_gNULLGENOobj->m_OneSNP_StdGeno;
for(size_t k=0; k< chrlength; k++){
   starti = ptr_gNULLGENOobj->startIndexVec[k];
   endi = ptr_gNULLGENOobj->endIndexVec[k];
//  std::cout << "debug2" << std::endl;
  if((starti != -1) && (endi != -1)){
        for(int i=starti; i<= endi; i++){
                        ptr_gNULLGENOobj->Get_OneSNP_StdGeno(i, temp);
                        (ptr_gNULLGENOobj->mtx_DiagStd_LOCO).col(k) = (ptr_gNULLGENOobj->mtx_DiagStd_LOCO).col(k) + (*temp) % (*temp);
                        ptr_gNULLGENOobj->Msub_MAFge_minMAFtoConstructGRM_byChr[k] = ptr_gNULLGENOobj->Msub_MAFge_minMAFtoConstructGRM_byChr[k] + 1;

        }
  (ptr_gNULLGENOobj->mtx_DiagStd_LOCO).col(k) = *ptr_gNULLGENOobj->Get_Diagof_StdGeno() -  (ptr_gNULLGENOobj->mtx_DiagStd_LOCO).col(k);
  }
}
}



// [[Rcpp::export]]
arma::fvec get_DiagofKin(){
    int Nnomissing = ptr_gNULLGENOobj->getNnomissing();
    arma::fvec x(Nnomissing);

    if(!(ptr_gNULLGENOobj->setKinDiagtoOne)){
           x  = (*ptr_gNULLGENOobj->Get_Diagof_StdGeno());
           int MminMAF = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
           x = x/MminMAF;
    }else{
           x  = arma::ones<arma::fvec>(Nnomissing);
    }
    return(x);
}





//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_usingSubMarker : public Worker
{
        // source vectors
        arma::fcolvec & m_bVec;
        unsigned int m_N;
        unsigned int m_M_Submarker;
        unsigned int m_M;
        arma::ivec subMarkerIndex ;

        // product that I have accumulated
        arma::fvec m_bout;


        // constructors
        CorssProd_usingSubMarker(arma::fcolvec & y)
                : m_bVec(y) {

                //m_Msub = ptr_gNULLGENOobj->getMsub();
                subMarkerIndex = getSubMarkerIndex();
                m_M_Submarker = subMarkerIndex.n_elem;
                m_N = ptr_gNULLGENOobj->getNnomissing();
                m_bout.zeros(m_N);
        }
        CorssProd_usingSubMarker(const CorssProd_usingSubMarker& CorssProd_usingSubMarker, Split)
                : m_bVec(CorssProd_usingSubMarker.m_bVec)
        {

                m_N = CorssProd_usingSubMarker.m_N;
                //m_M = CorssProd_usingSubMarker.m_M;
                m_M_Submarker = CorssProd_usingSubMarker.m_M_Submarker;
                subMarkerIndex = CorssProd_usingSubMarker.subMarkerIndex;
                m_bout.zeros(m_N);

        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                float val1;
                int j;
                for(unsigned int i = begin; i < end; i++){
                        j = subMarkerIndex[i];
//                      std::cout << "j: " << j << std::endl;
                        ptr_gNULLGENOobj->Get_OneSNP_StdGeno(j, &vec);
                        val1 = dot(vec,  m_bVec);
                        m_bout += val1 * (vec);
                }
        }

        // join my value with that of another InnerProduct
        void join(const  CorssProd_usingSubMarker & rhs) {
        m_bout += rhs.m_bout;
        }
};



// [[Rcpp::export]]
arma::fvec parallelCrossProd_usingSubMarker(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        int m_M_Submarker = getSubMarkerNum();
        CorssProd_usingSubMarker CorssProd_usingSubMarker(bVec);
        parallelReduce(0, m_M_Submarker, CorssProd_usingSubMarker);
        return CorssProd_usingSubMarker.m_bout/m_M_Submarker;
}



// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin_usingSubMarker(arma::fcolvec& bVec){

        arma::fvec crossProdVec = parallelCrossProd_usingSubMarker(bVec) ;

        return(crossProdVec);
}



//The code below is from http://gallery.rcpp.org/articles/parallel-inner-product/
struct InnerProduct : public Worker
{
   // source vectors
   std::vector<float> x;
   std::vector<float> y;

   // product that I have accumulated
   float product;

   // constructors
   InnerProduct(const std::vector<float> x, const std::vector<float> y)
      : x(x), y(y), product(0) {}
   InnerProduct(const InnerProduct& innerProduct, Split)
      : x(innerProduct.x), y(innerProduct.y), product(0) {}

   // process just the elements of the range I've been asked to
   void operator()(std::size_t begin, std::size_t end) {
      product += std::inner_product(x.begin() + begin,
                                    x.begin() + end,
                                    y.begin() + begin,
                                    0.0);
   }

   // join my value with that of another InnerProduct
   void join(const InnerProduct& rhs) {
     product += rhs.product;
   }
};



// [[Rcpp::export]]
float parallelInnerProduct(std::vector<float> &x, std::vector<float> &y) {

   int xsize = x.size();
   // declare the InnerProduct instance that takes a pointer to the vector data
   InnerProduct innerProduct(x, y);

   // call paralleReduce to start the work
   parallelReduce(0, x.size(), innerProduct);

   // return the computed product
   return innerProduct.product/xsize;
}

// [[Rcpp::export]]
Rcpp::List createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){

        int nSubMarker = markerIndexVec.n_elem;
        int Ntotal = wVec.n_elem;
        std::vector<unsigned int>     iIndexVec;
        std::vector<unsigned int>     iIndexVec2;
        std::vector<unsigned int>     jIndexVec;
        std::vector<unsigned int>     jIndexVec2;
        std::vector<unsigned int>     allIndexVec;
        std::vector<float>     kinValueVec;
        std::vector<float>     kinValueVec2;
        std::vector<float> stdGenoMultiMarkers;
        stdGenoMultiMarkers.resize(Ntotal*nSubMarker);

        //std::cout << "createSparseKin1" << std::endl;
        size_t sizeTemp;
        float kinValue;
        float kinValueTemp;

        Get_MultiMarkersBySample_StdGeno(markerIndexVec, stdGenoMultiMarkers);
        std::cout << "createSparseKin2" << std::endl;
        arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), nSubMarker, Ntotal);



        for(unsigned int i=0; i< Ntotal; i++){
              for(unsigned int j = i; j < Ntotal; j++){
                        //kinValueTemp = arma::dot(stdGenoMultiMarkersMat.row(i), stdGenoMultiMarkersMat.row(j));
                        if(j > i){
                                kinValueTemp = arma::dot(stdGenoMultiMarkersMat.col(i), stdGenoMultiMarkersMat.col(j));
                                kinValueTemp = kinValueTemp/nSubMarker;
                                if(kinValueTemp >= relatednessCutoff){
                                        iIndexVec.push_back(i);
                                        jIndexVec.push_back(j);

                                }
                        }else{
                                iIndexVec.push_back(i);
                                jIndexVec.push_back(j);
                        }
                }
        }

        arma::fvec * temp = &(ptr_gNULLGENOobj->m_OneSNP_StdGeno);
        size_t ni = iIndexVec.size();
        kinValueVec.resize(ni);
        std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

        int Mmarker = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
                //ptr_gNULLGENOobj->getM();
        for(size_t i=0; i< Mmarker; i++){
                ptr_gNULLGENOobj->Get_OneSNP_StdGeno(i, temp);
                for(size_t j=0; j < ni; j++){
                        kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
                }

        }


        for(size_t j=0; j < ni; j++){
                if(kinValueVec[j] >= relatednessCutoff){
        //      std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
                        kinValueVec[j] = tauVec(1)*kinValueVec[j];
                        iIndexVec2.push_back(iIndexVec[j]+1);
                        jIndexVec2.push_back(jIndexVec[j]+1);
                        if(iIndexVec[j] == jIndexVec[j]){
                                kinValueVec[j] = kinValueVec[j] + tauVec(0)/(wVec(iIndexVec[j]));
                        }
                        kinValueVec2.push_back(kinValueVec[j]);
                }

        }

        return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);
}




// [[Rcpp::export]]
Rcpp::List refineKin(float relatednessCutoff){
        std::vector<unsigned int>     iIndexVec2;
        std::vector<unsigned int>     jIndexVec2;
//      std::vector<float>     kinValueVec;
        std::vector<float>     kinValueVec2;
 //       std::vector<float>     kinValueVec_orig; //for test original kinship

        arma::fvec * temp = &(ptr_gNULLGENOobj->m_OneSNP_StdGeno);
        (*temp).clear();
        size_t ni = ptr_gNULLGENOobj->indiceVec.size();
        std::cout << "ni: " << ni << std::endl;

        initKinValueVecFinal(ni);

        int Mmarker = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();

        arma::fvec kinValueVecTemp2;
        arma::fvec GRMvec;
        GRMvec.set_size(ni);
        for(size_t i=0; i< Mmarker; i++){
                float freqv = ptr_gNULLGENOobj->alleleFreqVec[i];

                ptr_gNULLGENOobj->Get_OneSNP_Geno(i);
                float invstdv = ptr_gNULLGENOobj->invstdvVec[i];
                ptr_gNULLGENOobj->setSparseKinLookUpArr(freqv, invstdv);


                parallelcalsparseGRM(GRMvec);
                parallelsumTwoVec(GRMvec);
                (*temp).clear();
        }



        int a1;
        int a2;
        for(size_t j=0; j < ni; j++){
                ptr_gNULLGENOobj->kinValueVecFinal[j] = (ptr_gNULLGENOobj->kinValueVecFinal[j]) /(Mmarker);
                if((ptr_gNULLGENOobj->kinValueVecFinal[j]) >= relatednessCutoff){
                                 a1 = (ptr_gNULLGENOobj->indiceVec)[j].first + 1;
                                 a2 = (ptr_gNULLGENOobj->indiceVec)[j].second + 1;
                                 iIndexVec2.push_back(a1);
                                 jIndexVec2.push_back(a2);

                        kinValueVec2.push_back((ptr_gNULLGENOobj->kinValueVecFinal)[j]);
                }

        }



        std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;
        return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);
}

// [[Rcpp::export]]
arma::fmat getColfromStdGenoMultiMarkersMat(arma::uvec & a){
        return((ptr_gNULLGENOobj->stdGenoMultiMarkersMat).cols(a));
}

// [[Rcpp::export]]
int getNColStdGenoMultiMarkersMat(){
        return((ptr_gNULLGENOobj->stdGenoMultiMarkersMat).n_cols);
}

// [[Rcpp::export]]
int getNRowStdGenoMultiMarkersMat(){
        return((ptr_gNULLGENOobj->stdGenoMultiMarkersMat).n_rows);
}


// [[Rcpp::export]]
void setSubMarkerIndex(arma::ivec &subMarkerIndexRandom){
        ptr_gNULLGENOobj->subMarkerIndex = subMarkerIndexRandom;
        int Nnomissing = ptr_gNULLGENOobj->getNnomissing();
        (ptr_gNULLGENOobj->stdGenoMultiMarkersMat).set_size(subMarkerIndexRandom.n_elem, Nnomissing);
}

// [[Rcpp::export]]
void setRelatednessCutoff(float a){
        ptr_gNULLGENOobj->relatednessCutoff = a;
}


// [[Rcpp::export]]
double innerProduct(Rcpp::NumericVector x, Rcpp::NumericVector y) {
   return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}


// [[Rcpp::export]]
arma::fvec getDiagOfSigma_noV(arma::fvec& wVec, arma::fvec& tauVec, bool LOCO){
        int Nnomissing = wVec.n_elem;
        arma::fvec diagVec(Nnomissing);
        arma::sp_fvec diagVecV0;
        arma::fvec diagVecG, diagVecV, diagVecG_I, diagVecG_T, diagVecG_IT,diagVecV_I, diagVecV_T, diagVecV_IT;
        unsigned int tauind = 0;

        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){

          if(!(ptr_gNULLGENOobj->setKinDiagtoOne)){
             diagVec = tauVec(0)/wVec;
             tauind = tauind + 1;
             diagVecG = g_spGRM.diag();
             diagVec = diagVec + tauVec(tauind) * diagVecG;
             tauind = tauind + 1;
          }else{ //if(!(ptr_gNULLGENOobj->setKinDiagtoOne)){
            diagVec = tauVec(1) + tauVec(0)/wVec;
            tauind = tauind + 2;
          }


        }else{ //if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
                diagVec = tauVec(0)/wVec;
                tauind = tauind + 1;
                diagVecG = g_spGRM.diag();
                diagVecG_I = diagVecG.elem(g_I_longl_vec);
                diagVec = diagVec + tauVec(tauind) * diagVecG_I;
                tauind = tauind + 1;

                if(g_T_longl_mat.n_rows > 0){
                  diagVecG_IT = diagVecG_I % g_T_longl_vec;
                  diagVecG_T = diagVecG_IT % g_T_longl_vec;
                  diagVecG_IT = 2 * diagVecG_IT;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_IT;
                  tauind = tauind + 1;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_T;
                  tauind = tauind + 1;
                }

        }

        for(unsigned int i=0; i< Nnomissing; i++){
                if(diagVec(i) < 1e-4){
                        diagVec(i) = 1e-4 ;
                }
        }

        return(diagVec);
}



// [[Rcpp::export]]
arma::fvec getDiagOfSigma_V(arma::fvec& wVec, float tauVal, float tauVal0){
        int Nnomissing = wVec.n_elem;
        arma::fvec diagVec(Nnomissing);
        arma::sp_fvec diagVecV0;
        arma::fvec diagVecG, diagVecV, diagVecG_I, diagVecG_T, diagVecG_IT,diagVecV_I, diagVecV_T, diagVecV_IT;
        unsigned int tauind = 0;

        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){

             diagVec = tauVal0/wVec;
             tauind = tauind + 1;
             //diagVecG = spV.diag();
             diagVec = diagVec + tauVal;
             tauind = tauind + 1;


        }else{ //if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
                diagVec = tauVal0/wVec;
                //tauind = tauind + 1;
                //diagVecG = spV.diag();
                //diagVecG_I = diagVecG.elem(g_I_longl_vec);
                //diagVec = diagVec + tauVec(indforV) * diagVecG_I;
		diagVec = diagVec + tauVal;

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




// [[Rcpp::export]]
arma::fcolvec getCrossprod_noV(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, bool LOCO){

        arma::fcolvec crossProdVec;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0 && tauVec.n_elem == 2){
                crossProdVec = tauVec(0)*(bVec % (1/wVec));
                return(crossProdVec);
        }
        //
        arma::fvec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;

        unsigned int tau_ind = 0;
        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){ //it must have specified GRM
                  crossProd1 = getCrossprodMatAndKin(bVec, LOCO);
                  crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;
                  tau_ind = tau_ind + 2;
        }else{

                  Ibvec = g_I_longl_mat.t() * bVec;

                  GRM_I_bvec = getCrossprodMatAndKin(Ibvec, LOCO);
                  crossProd1 = g_I_longl_mat * GRM_I_bvec;
                  crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;
                  tau_ind = tau_ind + 2;


                if(g_T_longl_mat.n_rows > 0){
                        Tbvec = g_T_longl_mat.t() * bVec;
                        GRM_T_bvec = getCrossprodMatAndKin(Tbvec, LOCO);
                        crossProdGRM_TGIb = g_T_longl_mat * GRM_I_bvec;
                        crossProdGRM_IGTb = g_I_longl_mat * GRM_T_bvec;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (crossProdGRM_TGIb + crossProdGRM_IGTb);
                        tau_ind = tau_ind + 1;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (g_T_longl_mat * GRM_T_bvec);
                        tau_ind = tau_ind + 1;
                }

        }

        return(crossProdVec);
}




// [[Rcpp::export]]
arma::fcolvec getCrossprod_V(arma::fcolvec& bVec, arma::fvec& wVec, float tauVal, float tauVal0){
	//indforV starts from 0
        arma::fcolvec crossProdVec;
        arma::fvec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;

        unsigned int tau_ind = 0;
        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){ //it must have specified GRM
                  crossProd1 = bVec;
                  crossProdVec = tauVal0*(bVec % (1/wVec)) + tauVal*crossProd1;
                  tau_ind = tau_ind + 2;
        }else{

                  Ibvec = g_I_longl_mat.t() * bVec;

                  GRM_I_bvec = Ibvec;
                  crossProd1 = g_I_longl_mat * GRM_I_bvec;
                  crossProdVec = tauVal0*(bVec % (1/wVec)) + tauVal*crossProd1;
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


// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_noV(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG, bool LOCO){
    // Start Timers
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    int Nnomissing = wVec.n_elem;
    arma::fvec xVec(Nnomissing);
    xVec.zeros();

    if(g_isStoreSigma){
        std::cout << " arma::spsolve(g_spSigma, bVec) 0" << std::endl;
        //xVec = arma::spsolve(g_spSigma, bVec);
        xVec = arma::spsolve(g_spSigma_noV, bVec);
        std::cout << " arma::spsolve(g_spSigma, bVec) 1" << std::endl;
    }else{
        arma::fvec rVec = bVec;
        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        minvVec = 1/getDiagOfSigma_noV(wVec, tauVec, LOCO);
        zVec = minvVec % rVec;

        float sumr2 = sum(rVec % rVec);
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;

        int iter = 0;

        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                arma::fcolvec ApVec = getCrossprod_noV(pVec, wVec, tauVec, LOCO);
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
                z1Vec = minvVec % r1Vec;

                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
}
        return(xVec);
}




// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_V(arma::fvec& wVec,  float tauVal, float tauVal0, arma::fvec& bVec, int maxiterPCG, float tolPCG){
    // Start Timers
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    int Nnomissing = wVec.n_elem;
    arma::fvec xVec(Nnomissing);
    xVec.zeros();

    if(g_isStoreSigma){
        std::cout << " arma::spsolve(g_spSigma, bVec) 0" << std::endl;
        //xVec = arma::spsolve(g_spSigma, bVec);
        xVec = arma::spsolve(g_spSigma_V, bVec);
        std::cout << " arma::spsolve(g_spSigma, bVec) 1" << std::endl;
    }else{
        arma::fvec rVec = bVec;
        arma::fvec r1Vec;
        arma::fvec crossProdVec(Nnomissing);
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
        minvVec = 1/getDiagOfSigma_V(wVec, tauVal, tauVal0);
        zVec = minvVec % rVec;

        float sumr2 = sum(rVec % rVec);
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;

        int iter = 0;

        while (sumr2 > tolPCG && iter < maxiterPCG) {
                iter = iter + 1;
                arma::fcolvec ApVec = getCrossprod_V(pVec, wVec, tauVal, tauVal0);
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

                float a = preA(0);
                xVec = xVec + a * pVec;
                r1Vec = rVec - a * ApVec;
                z1Vec = minvVec % r1Vec;

                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;
                sumr2 = sum(rVec % rVec);
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
}
        return(xVec);
}





// [[Rcpp::export]]
arma::fvec  getSigma_G_V(arma::fvec& wVec, float tauVal, float tauVal0, arma::fvec& Gvec, int maxiterPCG, float tolPCG){
        arma::fvec Sigma_iG;
        Sigma_iG = getPCG1ofSigmaAndVector_V(wVec, tauVal, tauVal0, Gvec, maxiterPCG, tolPCG);
        return(Sigma_iG);
}


// [[Rcpp::export]]
arma::fvec  getSigma_G_noV(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG, bool LOCO){
        arma::fvec Sigma_iG;
        Sigma_iG = getPCG1ofSigmaAndVector_noV(wVec, tauVec, Gvec, maxiterPCG, tolPCG, LOCO);
        return(Sigma_iG);
}




// [[Rcpp::export]]
void set_useGRMtoFitNULL(bool useGRMtoFitNULL){
        g_isGRM = useGRMtoFitNULL;
}


// [[Rcpp::export]]
void set_isSparseGRM(bool t_isSparseGRM){
        g_isSparseGRM = t_isSparseGRM;
}


// [[Rcpp::export]]
void set_store_sigma(bool isstoreSigma){
        //g_longl_vec = arma::conv_to< arma::fvec >::from(longlVec);
        g_isStoreSigma = isstoreSigma;
}

// [[Rcpp::export]]
void set_num_Kmat(int t_num_Kmat){
        //g_longl_vec = arma::conv_to< arma::fvec >::from(longlVec);
        g_num_Kmat = t_num_Kmat;
}


// [[Rcpp::export]]
int get_numofV(){
        int k = Kmat_vec.size();
        return(k);
}


// [[Rcpp::export]]
arma::umat set_covarianceidx_Mat(){
        //unsigned int k = Kmat_vec.size();
        unsigned int k = g_num_Kmat;

        unsigned int q = (k+1)/3;
        g_covarianceidxMat.set_size(q, 3);
        g_covarianceidxMat.zeros();
        arma::uvec g_covarianceidxVec(3);
        for(unsigned int j=0; j < q; j++){
                g_covarianceidxVec = {j*3+3, j*3+2, j*3+4};
                g_covarianceidxMat.row(j) = g_covarianceidxVec.t();
        }
        g_covarianceidxMat_col1 = g_covarianceidxMat.col(0) - 1;
        g_covarianceidxMat_col2 = g_covarianceidxMat.col(1) - 1;
        g_covarianceidxMat_col3 = g_covarianceidxMat.col(2) - 1;

        arma::uvec indexsubvec =  { 1, 2 };
        g_covarianceidxMat_notcol1 = arma::vectorise(g_covarianceidxMat.cols(indexsubvec));

        g_covarianceidxMat_notcol1 = g_covarianceidxMat_notcol1 - 1;

        return(g_covarianceidxMat);
}

/*
// [[Rcpp::export]]
void set_Vmat_vec_longlVar(){
        std::cout << "here" << std::endl;
        int k  = Kmat_vec.size();
        arma::sp_fmat spGRM_longl_0;
        arma::sp_fmat spGRM_longl;
        for(int j=0; j < 1+k; j++){
          if(j == 0){
            spGRM_longl_0 = g_spGRM;
          }else{
            //spGRM_longl = Kmat_vec[j-1];
            //spGRM_longl_0 = arma::conv_to< arma::sp_mat >::from(Kmat_vec[j-1]);
            spGRM_longl_0 = Kmat_vec[j-1];
          }
          //spGRM_longl_0 = g_longl_vec.t() % (spGRM_longl_0.each_row());
          for(int q=0; q < spGRM_longl_0.n_rows; q++){
          //for(int q=0; q < spGRM_longl.n_rows; q++){
                 spGRM_longl_0.row(q) = g_longl_vec.t() % spGRM_longl_0.row(q);
          }
          //spGRM_longl = arma::conv_to< arma::sp_fmat >::from(spGRM_longl_0);
          //Kmat_vec.push_back(spGRM_longl);
          Kmat_vec.push_back(spGRM_longl_0);

          //spGRM_longl_0 = spGRM_longl_0 % g_longl_vec;
          //spGRM_longl_0 = g_longl_vec.t() % (spGRM_longl_0.t().each_row());
          //spGRM_longl_0 = spGRM_longl_0.t();
          for(int q=0; q < spGRM_longl_0.n_cols; q++){
                 spGRM_longl_0.col(q) = spGRM_longl_0.col(q) % g_longl_vec;
                 //spGRM_longl.col(q) = spGRM_longl.col(q) % g_longl_vec;
          }
          //spGRM_longl = arma::conv_to< arma::sp_fmat >::from(spGRM_longl_0);
          //Kmat_vec.push_back(spGRM_longl);
          Kmat_vec.push_back(spGRM_longl_0);
          //spGRM_longl_0.clear();
          spGRM_longl.clear();
        }
}
i

*/

// [[Rcpp::export]]
void closeGenoFile_plink()
{
    if(g_isGRM && !g_isSparseGRM){	
  //genoToTest_plainDosage.test_genoGZfile.close();
        for (int i = 0; i < ptr_gNULLGENOobj->numofGenoArray; i++){
                (*ptr_gNULLGENOobj->genoVecofPointers[i]).clear();
                delete ptr_gNULLGENOobj->genoVecofPointers[i];
        }

        ptr_gNULLGENOobj->genoVecofPointers.clear();

        //ptr_gNULLGENOobj->genoVec.clear();
        ptr_gNULLGENOobj->invstdvVec.clear();
        ptr_gNULLGENOobj->ptrsubSampleInGeno.clear();
        ptr_gNULLGENOobj->alleleFreqVec.clear();
        ptr_gNULLGENOobj->m_OneSNP_Geno.clear();
        ptr_gNULLGENOobj->m_OneSNP_StdGeno.clear();
        ptr_gNULLGENOobj->m_DiagStd.clear();
        printf("closed the plinkFile!\n");
    }	
}


// [[Rcpp::export]]
int gettotalMarker(){
        int numMarker = ptr_gNULLGENOobj->getM();
        return(numMarker);
}

// [[Rcpp::export]]
arma::fvec getAlleleFreqVec(){
        return(ptr_gNULLGENOobj->alleleFreqVec);
}

// [[Rcpp::export]]
arma::ivec getMACVec(){
        return(ptr_gNULLGENOobj->MACVec);
}


// [[Rcpp::export]]
arma::ivec getMACVec_forVarRatio(){
        return(ptr_gNULLGENOobj->MACVec_forVarRatio);
}

// [[Rcpp::export]]
arma::ivec getIndexVec_forVarRatio(){
        return(ptr_gNULLGENOobj->markerIndexVec_forVarRatio);
}

// [[Rcpp::export]]
bool getIsVarRatioGeno(){
        return(ptr_gNULLGENOobj->isVarRatio);
}
// [[Rcpp::export]]
arma::ivec getSubMarkerIndex(){
        return(ptr_gNULLGENOobj->subMarkerIndex);
}

// [[Rcpp::export]]
std::vector<bool> getQCdMarkerIndex(){
        return(ptr_gNULLGENOobj->MarkerswithMAFge_minMAFtoConstructGRM_indVec);
}


// [[Rcpp::export]]
int getSubMarkerNum(){
        return(ptr_gNULLGENOobj->subMarkerIndex.n_elem);
}


void initKinValueVecFinal(int ni){
        ptr_gNULLGENOobj->kinValueVecFinal.resize(ni);
        std::fill(ptr_gNULLGENOobj->kinValueVecFinal.begin(), ptr_gNULLGENOobj->kinValueVecFinal.end(), 0);
};

// [[Rcpp::export]]
int getNnomissingOut(){
        return(ptr_gNULLGENOobj->getNnomissing());
}


// [[Rcpp::export]]
int getMsub_MAFge_minMAFtoConstructGRM(){
        return(ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_in());
}

// [[Rcpp::export]]
int getMsub_MAFge_minMAFtoConstructGRM_singleChr(){
        return(ptr_gNULLGENOobj->getMsub_MAFge_minMAFtoConstructGRM_singleChr_in());
}


// [[Rcpp::export]]
void Get_MultiMarkersBySample_StdGeno_Mat(){
        //ptr_gNULLGENOobj->subMarkerIndex
        //int m_M_Submarker = markerIndexVec.n_elem;
        int m_M_Submarker = getSubMarkerNum();
        //arma::fvec stdGenoMultiMarkers;
        int Nnomissing = ptr_gNULLGENOobj->getNnomissing();
          //int nSubMarker = markerIndexVec.n_elem;
          //int Ntotal = ptr_gNULLGENOobj->getNnomissing();
        //std::vector<float> stdGenoMultiMarkers;
        //stdGenoMultiMarkers.resize(Nnomissing*m_M_Submarker);

        int indexOfVectorPointer;
        int SNPIdxinVec;
        size_t Start_idx;
        size_t ind= 0;
        size_t indtotal = 0;
        unsigned char geno1;
        float freq;
        float invStd;
        int flag;
        int SNPIdx;

//      std::cout << "createSparseKin1d" << std::endl;
        for(size_t k=0; k< m_M_Submarker; k++){
                ind = 0;
                flag = 0;
                //SNPIdx = markerIndexVec[k];
                SNPIdx = (ptr_gNULLGENOobj->subMarkerIndex)[k];
                indexOfVectorPointer = SNPIdx/(ptr_gNULLGENOobj->numMarkersofEachArray);
                SNPIdxinVec = SNPIdx % (ptr_gNULLGENOobj->numMarkersofEachArray);
                Start_idx = (ptr_gNULLGENOobj->m_size_of_esi) * SNPIdxinVec;
                freq = (ptr_gNULLGENOobj->alleleFreqVec)[SNPIdx];
                invStd = (ptr_gNULLGENOobj->invstdvVec)[SNPIdx];
                if(k == 0){
                        std::cout << "freq: " << freq << " invStd: " << invStd << "  SNPIdx: " << SNPIdx << std::endl;
                }

                while(flag == 0){
//              std::cout << "createSparseKin1e" << std::endl;
                for(size_t i=Start_idx; i< Start_idx+(ptr_gNULLGENOobj->m_size_of_esi); i++){
                        geno1 = (ptr_gNULLGENOobj->genoVecofPointers)[indexOfVectorPointer]->at(i);
                        //std::cout << "createSparseKin1f" << std::endl;

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        (ptr_gNULLGENOobj->stdGenoMultiMarkersMat)(k, ind) = ((2-(a+b)) - 2*freq)* invStd;
//                      std::cout << "k,ind " << k << " " << ind << std::endl;
//                      std::cout << "(ptr_gNULLGENOobj->stdGenoMultiMarkersMat)(k, ind) " << (ptr_gNULLGENOobj->stdGenoMultiMarkersMat)(k, ind) << std::endl;

//                        stdGenoMultiMarkers[ind*m_M_Submarker+k] = ((2-(a+b)) - 2*freq)* invStd;;
//                      if(k == 0){
    //                    std::cout << "ind*m_M_Submarker+k: " << ind*m_M_Submarker+k << " stdGenoMultiMarkers[ind*m_M_Submarker+k]: " << stdGenoMultiMarkers[ind*m_M_Submarker+k] <<  std::endl;
  //              }


                        indtotal++;
                        ind++;
                        geno1 = geno1 >> 1;

                                if(ind == Nnomissing){
                                        flag = 1;
                                        break;
                                }
                        }// end of for(int j=0; j<4; j++){
                    }// end of for(size_t i=Start_idx
                } //end of while(flag == 0){

        }

        std::cout << "stdGenoMultiMarkersMat.n_rows: " << ptr_gNULLGENOobj->stdGenoMultiMarkersMat.n_rows << std::endl;
        std::cout << "stdGenoMultiMarkersMat.n_cols: " << ptr_gNULLGENOobj->stdGenoMultiMarkersMat.n_cols << std::endl;
//      arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), m_M_Submarker, Nnomissing);

//      return(stdGenoMultiMarkersMat);
        //std::cout << "stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] " << stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] << std::endl;

}


// [[Rcpp::export]]
void Get_MultiMarkersBySample_StdGeno(arma::fvec& markerIndexVec, std::vector<float> &stdGenoMultiMarkers){

//      std::cout << "createSparseKin1c" << std::endl;
        int indexOfVectorPointer;
        int SNPIdxinVec;
        size_t Start_idx;
        size_t ind= 0;
        size_t indtotal = 0;
        unsigned char geno1;
        float freq;
        float invStd;
        int flag;
        int SNPIdx;

        int m_M_Submarker = markerIndexVec.n_elem;
        //arma::fvec stdGenoMultiMarkers;
        int Nnomissing = ptr_gNULLGENOobj->getNnomissing();


//      std::cout << "createSparseKin1d" << std::endl;
        for(size_t k=0; k< m_M_Submarker; k++){
                ind = 0;
                flag = 0;
                SNPIdx = markerIndexVec[k];
                indexOfVectorPointer = SNPIdx/(ptr_gNULLGENOobj->numMarkersofEachArray);
                SNPIdxinVec = SNPIdx % (ptr_gNULLGENOobj->numMarkersofEachArray);
                Start_idx = (ptr_gNULLGENOobj->m_size_of_esi) * SNPIdxinVec;
                freq = (ptr_gNULLGENOobj->alleleFreqVec)[SNPIdx];
                invStd = (ptr_gNULLGENOobj->invstdvVec)[SNPIdx];
                //if(k == 0){
                //      std::cout << "freq: " << freq << " invStd: " << invStd << "  SNPIdx: " << SNPIdx << std::endl;
                //}

                while(flag == 0){
//              std::cout << "createSparseKin1e" << std::endl;
                for(size_t i=Start_idx; i< Start_idx+(ptr_gNULLGENOobj->m_size_of_esi); i++){
                        geno1 = (ptr_gNULLGENOobj->genoVecofPointers)[indexOfVectorPointer]->at(i);
                        //std::cout << "createSparseKin1f" << std::endl;

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        stdGenoMultiMarkers[ind*m_M_Submarker+k] = ((2-(a+b)) - 2*freq)* invStd;;
//                      stdGenoMultiMarkers[ind*m_M_Submarker+k] = 2-(a+b);
//                      if(k == 0){
    //                    std::cout << "ind*m_M_Submarker+k: " << ind*m_M_Submarker+k << " stdGenoMultiMarkers[ind*m_M_Submarker+k]: " << stdGenoMultiMarkers[ind*m_M_Submarker+k] <<  std::endl;
  //              }


                        indtotal++;
                        ind++;
                        geno1 = geno1 >> 1;

                                if(ind == Nnomissing){
                                        flag = 1;
                                        break;
                                }
                        }// end of for(int j=0; j<4; j++){
                    }// end of for(size_t i=Start_idx
                } //end of while(flag == 0){

        }

        //std::cout << "stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] " << stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] << std::endl;

}




//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd : public Worker
{
        // source vectors
        arma::fcolvec & m_bVec;
        unsigned int m_N;
        unsigned int m_M;

        // product that I have accumulated
        arma::fvec m_bout;
        int Msub_mafge1perc;

        // constructors
        CorssProd(arma::fcolvec & y)
                : m_bVec(y) {

                m_M = ptr_gNULLGENOobj->getM();
                m_N = ptr_gNULLGENOobj->getNnomissing();
                m_bout.zeros(m_N);
                Msub_mafge1perc=0;
                //ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        }
        CorssProd(const CorssProd& CorssProd, Split)
                : m_bVec(CorssProd.m_bVec)
        {

                m_N = CorssProd.m_N;
                m_M = CorssProd.m_M;
                m_bout.zeros(m_N);
                Msub_mafge1perc=0;
                //CorssProd.Msub_mafge1perc;

        }
        // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                for(unsigned int i = begin; i < end; i++){
                        //if(ptr_gNULLGENOobj->alleleFreqVec[i] >= minMAFtoConstructGRM && ptr_gNULLGENOobj->alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
                                ptr_gNULLGENOobj->Get_OneSNP_StdGeno(i, &vec);
                                float val1 = dot(vec,  m_bVec);
                                m_bout += val1 * (vec) ;
                                Msub_mafge1perc += 1;
                        //}
                        //std::cout << "i: " << i << std::endl;
                        //for(unsigned int j = 0; j < 10; j++){
                        //      std::cout << "m_bVec[j] " << m_bVec[j] << std::endl;
                        //      std::cout << "vec[j] " << vec[j] << std::endl;
                        //}

                        //m_bout += val1 * (vec) / m_M;
                }
        }

        // join my value with that of another InnerProduct
        void join(const CorssProd & rhs) {
                m_bout += rhs.m_bout;
                Msub_mafge1perc += rhs.Msub_mafge1perc;
        }
};



//http://gallery.rcpp.org/articles/parallel-inner-product/
struct CorssProd_LOCO : public Worker
{
        // source vectors
        arma::fcolvec & m_bVec;
        unsigned int m_N;
        unsigned int m_Msub;
        unsigned int m_M;
        int startIndex;
        int endIndex;
        // product that I have accumulated
        arma::fvec m_bout;
        unsigned int m_Msub_mafge1perc;

        // constructors
        CorssProd_LOCO(arma::fcolvec & y)
                : m_bVec(y) {

                m_Msub = ptr_gNULLGENOobj->getMsub(); //LOCO
                startIndex = ptr_gNULLGENOobj->getStartIndex();
                endIndex = ptr_gNULLGENOobj->getEndIndex();
                m_M = ptr_gNULLGENOobj->getM(); //LOCO
                m_N = ptr_gNULLGENOobj->getNnomissing();
                m_bout.zeros(m_N);
                m_Msub_mafge1perc=0;
                //ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        }
        CorssProd_LOCO(const CorssProd_LOCO& CorssProd_LOCO, Split)
                : m_bVec(CorssProd_LOCO.m_bVec)
        {

                m_N = CorssProd_LOCO.m_N;
                m_M = CorssProd_LOCO.m_M;
                m_Msub = CorssProd_LOCO.m_Msub;
                startIndex = ptr_gNULLGENOobj->getStartIndex();
                endIndex = ptr_gNULLGENOobj->getEndIndex();
                m_bout.zeros(m_N);
                m_Msub_mafge1perc=0;
                //ptr_gNULLGENOobj->getnumberofMarkers_byChr(uint chr);
        }

           // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                float val1;
                for(unsigned int i = begin; i < end; i++){
                        ptr_gNULLGENOobj->Get_OneSNP_StdGeno(i, &vec);
                //      if(i >= startIndex && i <= endIndex){
                //              val1 = 0;
                                        //if(endIndex == 4){
                                        //              cout << "i: " << i << endl;
                                        //}
                //      }else{
                        val1 = dot(vec,  m_bVec);
                        m_Msub_mafge1perc += 1;
                //      }
                        m_bout += val1 * (vec);
                }
        }

        // join my value with that of another InnerProduct
        void join(const CorssProd_LOCO & rhs) {
        m_bout += rhs.m_bout;
        m_Msub_mafge1perc += rhs.m_Msub_mafge1perc;
        }
};


// [[Rcpp::export]]
arma::fvec parallelCrossProd(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        //int M = ptr_gNULLGENOobj->getM();
        //int Msub_mafge1perc = ptr_gNULLGENOobj->getMmafge1perc();
        int Msub_mafge1perc = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        CorssProd CorssProd(bVec);

  // call paralleReduce to start the work
        parallelReduce(0, Msub_mafge1perc, CorssProd);

        //cout << "print test; M: " << M << endl;
        //for(int i=0; i<10; ++i)
        //{
        //        cout << (CorssProd.m_bout)[i] << ' ' << endl;
        //        cout << bVec[i] << ' ' << endl;
        //        cout << (CorssProd.m_bout/M)[i] << ' ' << endl;
        //}
        ////cout << endl;
  // return the computed product
        //std::cout << "number of markers with maf ge " << minMAFtoConstructGRM << " is " << CorssProd.Msub_mafge1perc << std::endl;
        return CorssProd.m_bout/(CorssProd.Msub_mafge1perc);
        //return CorssProd.m_bout;
}

// [[Rcpp::export]]
float innerProductFun(std::vector<float> &x, std::vector<float> & y) {
   return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}




// [[Rcpp::export]]
arma::fvec parallelCrossProd_full(arma::fcolvec & bVec, int & markerNum) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        //int M = ptr_gNULLGENOobj->getM();
        //
        int Msub_mafge1perc = ptr_gNULLGENOobj->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        CorssProd CorssProd(bVec);

        //std::cout << "Msub_mafge1perc ok  " << Msub_mafge1perc << std::endl;
  // call paralleReduce to start the work
        parallelReduce(0, Msub_mafge1perc, CorssProd);
        markerNum = CorssProd.Msub_mafge1perc;
        //std::cout << "markerNum " << markerNum << std::endl;

        //cout << "print test; M: " << M << endl;
        //for(int i=0; i<10; ++i)
        //{
        //        cout << (CorssProd.m_bout)[i] << ' ' << endl;
        //        cout << bVec[i] << ' ' << endl;
        //        cout << (CorssProd.m_bout/M)[i] << ' ' << endl;
        //}
        ////cout << endl;
  // return the computed product
        //std::cout << "number of markers with maf ge " << minMAFtoConstructGRM << " is " << CorssProd.Msub_mafge1perc << std::endl;
        return CorssProd.m_bout;
        //return CorssProd.m_bout;
}


// [[Rcpp::export]]
arma::fvec parallelCrossProd_LOCO(arma::fcolvec & bVec) {

  // declare the InnerProduct instance that takes a pointer to the vector data
        //int Msub = ptr_gNULLGENOobj->getMsub();
        //int M = ptr_gNULLGENOobj->getM();
        int numberMarker_full = 0;
        arma::fvec outvec = parallelCrossProd_full(bVec, numberMarker_full);

        //CorssProd_LOCO CorssProd_LOCO(bVec);
        CorssProd CorssProd(bVec);
  // call paralleReduce to start the work
        int startIndex = ptr_gNULLGENOobj->getStartIndex();
        int endIndex = ptr_gNULLGENOobj->getEndIndex();

        parallelReduce(startIndex, endIndex+1, CorssProd);



        outvec = outvec - CorssProd.m_bout;

        /*
        for(int i=0; i<10; ++i)
        {
                std::cout << (outvec)[i] << ' ';
        }
        std::cout << std::endl;
*/

        int markerNum = numberMarker_full - CorssProd.Msub_mafge1perc;
        //std::cout << "markerNum: " << markerNum << std::endl;
        // return the computed product
        //cout << "Msub: " << Msub << endl;
        //for(int i=0; i<100; ++i)
        //{
        //      cout << (CorssProd_LOCO.m_bout/Msub)[i] << ' ';
        //}
        //cout << endl;
        //return CorssProd_LOCO.m_bout/Msub;
        //return CorssProd_LOCO.m_bout/(CorssProd_LOCO.m_Msub_mafge1perc);
        return outvec/markerNum;
}

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
struct indicesRelatedSamples : public RcppParallel::Worker {

  int  Ntotal;
  tbb::concurrent_vector< std::pair<int, int> > &output;

  indicesRelatedSamples(int Ntotal, tbb::concurrent_vector< std::pair<int, int> > &output) :
    Ntotal(Ntotal), output(output) {}


  void operator()(std::size_t begin, size_t end) {
    int m_M_Submarker = getSubMarkerNum();
    for(std::size_t k=begin; k < end; k++) {
      int i = (int)(k / Ntotal);
      int j = (int)(k % Ntotal);
      if((j <= i)){
                        i = Ntotal - i - 2;
                        j = Ntotal - j - 1;
      }
      //std::cout << "i,j,k debug: " << i << " " << j << " " << k << std::endl;
      float kinValueTemp = arma::dot((ptr_gNULLGENOobj->stdGenoMultiMarkersMat).col(i), (ptr_gNULLGENOobj->stdGenoMultiMarkersMat).col(j));
      kinValueTemp = kinValueTemp/m_M_Submarker;
      if(kinValueTemp >=  ptr_gNULLGENOobj->relatednessCutoff) {
        output.push_back( std::pair<int, int>(i, j) );
      }
    }
  }

};


// [[Rcpp::export]]
void printComb(int N){
  int x = N*(N-1)/2 - 1;
  for(std::size_t k=0; k < x; k++) {
      int i = k / N;
      int j = k % N;
      if((j < i)){
                        i = N - i - 2;
                        j = N - j - 1;
      }
     std::cout << "i,j " << i << "," << j << std::endl;
  }

}


// [[Rcpp::export]]
void findIndiceRelatedSample(){

  int Ntotal = ptr_gNULLGENOobj->getNnomissing();
//  tbb::concurrent_vector< std::pair<float, float> > output;

//  indicesRelatedSamples indicesRelatedSamples(Ntotal,output);
  indicesRelatedSamples indicesRelatedSamples(Ntotal,ptr_gNULLGENOobj->indiceVec);

  long int Ntotal2 = (long int)Ntotal;

  long int totalCombination = Ntotal2*(Ntotal2-1)/2 - 1;
  std::cout << "Ntotal: " << Ntotal << std::endl;
  std::cout << std::numeric_limits<int>::max() << std::endl;
  std::cout << std::numeric_limits<long int>::max() << std::endl;
  std::cout << std::numeric_limits<long long int>::max() << std::endl;
  std::cout << "totalCombination: " << totalCombination << std::endl;
  long int x = 1000001;
  int b = (int)(x / Ntotal);
  int a = (int)(x % Ntotal);
  std::cout << "a " << a << std::endl;
  std::cout << "b " << b << std::endl;

  parallelFor(0, totalCombination, indicesRelatedSamples);

//  arma::fmat xout(output.size()+Ntotal,2);

//  for(int i=0; i<output.size(); i++) {
//    xout(i,0) = output[i].first;
//    xout(i,1) = output[i].second;
//  }
//  for(int i=output.size(); i < output.size()+Ntotal; i++) {
//    xout(i,0) = i - output.size();
//    xout(i,1) = xout(i,0);
//  }

/*
  for(int i=0; i < Ntotal; i++){
    (ptr_gNULLGENOobj->indiceVec).push_back( std::pair<int, int>(i, i) );
  }
*/

//  return(xout);
}


struct sparseGRMUsingOneMarker : public Worker {
   // input matrix to read from
  // arma::imat & iMat;
   // output matrix to write to
   arma::fvec & GRMvec;

   //int M = ptr_gNULLGENOobj->getM();
   // initialize from Rcpp input and output matrixes (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
//   sparseGRMUsingOneMarker(arma::imat & iMat, arma::fvec &GRMvec)
//      : iMat(iMat), GRMvec(GRMvec) {}


  sparseGRMUsingOneMarker(arma::fvec &GRMvec)
      : GRMvec(GRMvec) {}


   // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
            // rows we will operate on
//            int iint = iMat(i,0);
//            int jint = iMat(i,1);
           int iint = (ptr_gNULLGENOobj->indiceVec)[i].first;
           int jint = (ptr_gNULLGENOobj->indiceVec)[i].second;
/*
            float ival = ptr_gNULLGENOobj->m_OneSNP_StdGeno(iint);
            float jval = ptr_gNULLGENOobj->m_OneSNP_StdGeno(jint);
            // write to output matrix
            //rmat(i,j) = sqrt(.5 * (d1 + d2));
            GRMvec(i) = ival*jval/M;
*/
        //use Look-Up table for calucate GRMvec(i)
            int ival = ptr_gNULLGENOobj->m_OneSNP_Geno(iint);
            int jval = ptr_gNULLGENOobj->m_OneSNP_Geno(jint);
            GRMvec(i) = ptr_gNULLGENOobj->sKinLookUpArr[ival][jval];

      }
   }
};


// [[Rcpp::export]]
void parallelcalsparseGRM(arma::fvec &GRMvec) {

//  int n1 = ptr_gNULLGENOobj->indiceVec.size();
  // allocate the output matrix
  //GRMvec.set_size(n1);
//  std::cout << "OKKK3: "  << std::endl;
//  sparseGRMUsingOneMarker sparseGRMUsingOneMarker(iMat, GRMvec);
  sparseGRMUsingOneMarker sparseGRMUsingOneMarker(GRMvec);
//  std::cout << "OKKK4: "  << std::endl;

//  std::cout << "n1 " << n1 << std::endl;
//  std::cout << "iMat.n_cols " << iMat.n_cols << std::endl;
  // call parallelFor to do the work
//  parallelFor(0, iMat.n_rows, sparseGRMUsingOneMarker);
  parallelFor(0, (ptr_gNULLGENOobj->indiceVec).size(), sparseGRMUsingOneMarker);

  // return the output matrix
  // return GRMvec;
}

struct sumTwoVec : public Worker
{
   // source vectors
   arma::fvec &x;

   arma::fvec &sumVec;

   //int M = ptr_gNULLGENOobj->getM();
   // constructors
   sumTwoVec(arma::fvec &x,arma::fvec &sumVec)
      : x(x), sumVec(sumVec) {}

     // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
            // rows we will operate on
            sumVec(i) = x(i)+(ptr_gNULLGENOobj->kinValueVecFinal)[i];
            (ptr_gNULLGENOobj->kinValueVecFinal)[i] = sumVec(i);
      }
   }

};

// [[Rcpp::export]]
void  parallelsumTwoVec(arma::fvec &x) {
  int n1 = x.n_elem;
  // allocate the output matrix
  arma::fvec sumVec;
  sumVec.set_size(n1);

  sumTwoVec sumTwoVec(x, sumVec);

  // call parallelFor to do the work
  parallelFor(0, x.n_elem, sumTwoVec);

}



// [[Rcpp::export]]
void setgenoNULL(){
	ptr_gNULLGENOobj = new NullGENO::NullGenoClass();
}
// [[Rcpp::export]]
void setgeno(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool isDiagofKinSetAsOne)
{
        int start_s=clock();

	//ptr_gNULLGENOobj = new NullGENO::NullGenoClass(bedfile, bimfile, famfile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne
//			);


        ptr_gNULLGENOobj->setGenoObj(bedfile, bimfile, famfile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne);
        //ptr_gNULLGENOobj->printAlleleFreqVec();
        //ptr_gNULLGENOobj->printGenoVec();
        int stop_s=clock();
        cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
}




// [[Rcpp::export]]
arma::ivec Get_OneSNP_Geno(int SNPIdx)
{

        arma::ivec temp = * ptr_gNULLGENOobj->Get_OneSNP_Geno(SNPIdx);
        return(temp);

}

// [[Rcpp::export]]
arma::ivec Get_OneSNP_Geno_forVarRatio(int SNPIdx)
{

        arma::ivec temp = * ptr_gNULLGENOobj->Get_OneSNP_Geno_forVarRatio(SNPIdx);
        return(temp);

}
// [[Rcpp::export]]
arma::fvec Get_OneSNP_StdGeno(int SNPIdx)
{
        arma::fvec temp;
        ptr_gNULLGENOobj->Get_OneSNP_StdGeno(SNPIdx, & temp);
        return(temp);

}

// [[Rcpp::export]]
void print_g_n_unique(){
	std::cout << "g_n_unique print " << g_n_unique << std::endl;
}


/*
// [[Rcpp::export]]
void set_var_weights(arma::vec & t_var_weights){
	arma::fvec t_var_weights_f = arma::conv_to< arma::fvec >::from(t_var_weights); 
	g_var_weights = t_var_weights_f;
}
*/

void  get_indexinAnotherVector(std::vector<uint> & nonzeroInd_orig, arma::uvec & dupInd, arma::uvec &  nonzeroInd_new_arma){
   std::vector<uint> nonzeroInd_new;
   arma::uvec newIndfori;
   for (unsigned int i = 0; i < nonzeroInd_orig.size(); i++) {
	newIndfori = arma::find(dupInd == nonzeroInd_orig.at(i));
	if(newIndfori.n_elem > 0){
		for (unsigned int j = 0; j < newIndfori.n_elem; j++){
			nonzeroInd_new.push_back(newIndfori(j));
		}
	}
	newIndfori.clear();	
   }
   nonzeroInd_new_arma =  arma::conv_to< arma::uvec >::from(nonzeroInd_new);
  nonzeroInd_new.clear(); 
}

// [[Rcpp::export]]
arma::sp_fmat get_sp_Sigma_to_R(){
	return(g_spSigma);
}



std::string join(std::vector<std::string> const &strings, std::string delim)
{
    std::stringstream ss;
    std::copy(strings.begin(), strings.end(),
        std::ostream_iterator<std::string>(ss, delim.c_str()));
    return ss.str();
}



// [[Rcpp::export]]
Rcpp::List get_SKAT_pvalue_Rcpp(arma::vec & Score, arma::mat &  Phi, arma::vec & r_corr) {
  //std::cout << "get_SKAT_pvalue_Rcpp0" << std::endl;
  
  // Call SKAT:::Met_SKAT_Get_Pvalue with specified arguments
  Rcpp::List out_SKAT_List;
  arma::vec Pvalue(3); 


std::string method = "optimal.adj";
//std::cout << "get_SKAT_pvalue_Rcpp2a0" << std::endl;
arma::mat Score_Resampling;
out_SKAT_List = Met_SKAT_Get_Pvalue_Rcpp(Score, Phi, r_corr, method, Score_Resampling);


//    out_SKAT_List = skatGetPvalue(_["Score"] = Score,
//                                  _["Phi"] = Phi,
 //                                 _["r.corr"] = r_corr,
 //                                 _["method"] = "optimal.adj",
 //                                 _["Score.Resampling"] = R_NilValue);
//std::cout << "get_SKAT_pvalue_Rcpp2a" << std::endl;
//  } catch (std::exception& e) {

//std::cout << "get_SKAT_pvalue_Rcpp2b" << std::endl;

//    Pvalue[0] = arma::datum::nan;
//    Pvalue[1] = arma::datum::nan;
//    Pvalue[2] = arma::datum::nan;
//    return List::create(Named("Pvalue_SKATO") = Pvalue[0],
//                      Named("Pvalue_Burden") = Pvalue[2],
//                      Named("Pvalue_SKAT") = Pvalue[1]);
//  }

  Rcpp::List paramList = out_SKAT_List["param"];
  //std::cout << "get_SKAT_pvalue_Rcpp2b" << std::endl;
  arma::vec rho;
  if(paramList.containsElementNamed("rho")){
      rho = Rcpp::as<arma::vec>(paramList["rho"]);
  }else{
      rho = r_corr;
  }

  //rho.print("rho");
  //std::cout << "get_SKAT_pvalue_Rcpp2c" << std::endl;
  arma::vec pvalueout = Rcpp::as<arma::vec>(out_SKAT_List["p_value"]); 
  //pvalueout.print("pvalueout");
  Rcpp::List pvaleachList;
  if(paramList.containsElementNamed("p_val_each")){
	pvaleachList = paramList["p_val_each"];
  }
  
  if(!(any(rho == 0.0)) && !(any(rho == 1.0)) && !arma::is_finite(pvalueout)){
      Pvalue[0] = pvalueout[0];
      Pvalue[1] = pvalueout[0];
      Pvalue[2] = pvalueout[0];
      int error_code = 0;
  }else{    
      arma::uvec pos00 = arma::find(rho == 0.0);
      arma::uvec pos01 = arma::find(rho == 1.0);
      Pvalue[0] = out_SKAT_List["p_value"];
     if(pos00.n_elem > 0 && paramList.containsElementNamed("p_val_each")){
      Pvalue[1] = pvaleachList[pos00[0]];
     }else{
	Pvalue[1] = out_SKAT_List["p_value"];
     }
     if(pos01.n_elem > 0 && paramList.containsElementNamed("p_val_each")){
      Pvalue[2] = pvaleachList[pos01[0]];
     }else{
	Pvalue[2]  = out_SKAT_List["p_value"];
     }
  }

  return Rcpp::List::create(Named("Pvalue_SKATO") = Pvalue[0],
                      Named("Pvalue_Burden") = Pvalue[2],
                      Named("Pvalue_SKAT") = Pvalue[1]);
}

//////// ---------- Main function for marker-level analysis - vectorize--------- ////////////


// [[Rcpp::export]]
void mainMarkerInCPP_multi(
                           std::string & t_genoType,     // "PLINK", "BGEN"
                           std::vector<std::string> & t_traitType,
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           bool & t_isMoreOutput,
                           bool & t_isImputation,
                           bool & t_isFirth)
{

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
  arma::vec mafVec(q);
  arma::vec macVec(q);
  arma::vec meanGVec(q);

  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q);

/*
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);
  std::vector<double> pvalVec(q, arma::datum::nan);
  std::vector<double> TstatVec(q, arma::datum::nan);
  std::vector<double> varTVec(q, arma::datum::nan);
  std::vector<double> pvalNAVec(q, arma::datum::nan);
*/
  std::vector<double> BetaVec;
  std::vector<double> seBetaVec;
  std::vector<double> pvalVec;
  std::vector<double> TstatVec;
  std::vector<double> varTVec;
  std::vector<double> pvalNAVec;



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
  std::vector<uint> indexNonZeroVec_multi;
  std::vector<uint> indexForMissing;

  //int n = ptr_gSAIGEobj->m_n;
  int n = ptr_gSAIGEobj->m_n;
  arma::vec t_GVec_cell(n);

  

  if(g_n_unique > 0){
    n = g_n_unique;
  }
  arma::vec t_GVec(n);
  arma::mat t_GMat(n, q);

  arma::ivec t_skipScoreIdx(q);
  t_skipScoreIdx.zeros();
  


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
  int passj = 0;


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
   if(!isReadMarker){
      //std::cout << "isReadMarker " << isReadMarker << std::endl;
      g_markerTestEnd = true;
      bool isEndFile = check_Vcf_end();
      break;
    }

    std::string pds = std::to_string(pd);
    std::string info = chr+":"+pds+":"+ref+":"+alt;
    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    int nG = t_GVec.n_elem;
    double MAC = MAF * n * (1 - missingRate) *2;
    // Quality Control (QC) based on missing rate, MAF, and MAC
    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff || imputeInfo < g_marker_minINFO_cutoff)){
      t_skipScoreIdx(i) = 1;
      continue;
    }else{

    // Check UTIL.cpp
    //
    //
    //arma::vec timeoutput3 = getTime();
    indexZeroVec.clear();
    indexNonZeroVec.clear();
    
    
    chrVec.at(passj) = chr;
    posVec.at(passj) = pds;
    refVec.at(passj) = ref;
    altVec.at(passj) = alt;
    markerVec.at(passj) = marker;               // marker IDs
    infoVec.at(passj) = info;    // marker information: CHR:POS:REF:ALT
    altFreqVec.at(passj) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(passj) = missingRate;
    imputationInfoVec.at(passj) = imputeInfo;


    flip = imputeGenoAndFlip(t_GVec, altFreq, altCounts,indexForMissing, g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);

    altFreqVec.at(passj) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(passj) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    MAC = std::min(altCounts, 2*n-altCounts);

    // analysis results for single-marker
    double Beta, seBeta, pval, pval_noSPA, Tstat, varT, gy;
    double Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c;

    bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
    //t_GMat.col(passj) = t_GVec;

    //if(g_n_unique >  0){
//	t_GVec_cell = g_I_longl_mat * t_GVec;
//	meanGVec(passj) = arma::mean(t_GVec_cell);
//    }else{
	meanGVec(passj) = arma::mean(t_GVec);
//    }
   
    //std::cout << "hererere" << t_GVec_cell.n_elem << " " << t_GMat.n_rows << std::endl;

    t_GMat.col(passj) = t_GVec;

    mafVec(passj) = MAF;
    macVec(passj) = MAC;
    //meanGVec(passj) = arma::mean(t_GVec);
    //std::cout << "t_GMat here1 " << std::endl;
    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
    indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
    indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);
    indexZeroVec.clear();
    indexNonZeroVec.clear();
    passj = passj + 1;
   }

}

    t_GMat.resize(n, passj);
    mafVec.resize(passj);     
    macVec.resize(passj);     
    BetaVec.resize(passj);     
    seBetaVec.resize(passj);     
    pvalVec.resize(passj);     
    pvalNAVec.resize(passj);     
    TstatVec.resize(passj);     
    varTVec.resize(passj);     
    meanGVec.resize(passj);

    //std::cout << "t_GMat here2 " << std::endl;
    //std::cout << "t_GMat(0,0) " << t_GMat(0,0) << std::endl;
    ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
    //double MACforvr = 20.5;
    //if(isSingleVarianceRatio){
    //   ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
    //}else{
    //   hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MACforvr, ptr_gSAIGEobj->m_flagSparseGRM_cur, true);
    //}

    std::vector <std::string> t_pval_str_vec;
    
    ptr_gSAIGEobj->getMarkerPval_multi(t_GMat, meanGVec, macVec, BetaVec, seBetaVec, pvalVec, pvalNAVec, TstatVec, varTVec, t_pval_str_vec);

    
    //std::cout << "t_GMat here3 " << std::endl;

    //output
unsigned int itt = 0;

writeOutfile_single(t_isMoreOutput,
      t_isImputation,
      isCondition,
      t_isFirth,
  mFirth,
  mFirthConverge,
  t_traitType,
  chrVec,
  posVec,
  markerVec,
  refVec,
  altVec,
  altCountsVec,
  altFreqVec,
  imputationInfoVec,
  missingRateVec,
  BetaVec,
  seBetaVec,
  TstatVec,
  varTVec,
  pvalVec,
  pvalNAVec,
  isSPAConvergeVec,
  Beta_cVec,
  seBeta_cVec,
  Tstat_cVec,
  varT_cVec,
  pval_cVec,
  pvalNA_cVec,
  AF_caseVec,
  AF_ctrlVec,
  N_caseVec,
  N_ctrlVec,
  N_case_homVec,
  N_ctrl_hetVec,
  N_case_hetVec,
  N_ctrl_homVec,
  N_Vec,
  g_isgxe,
  Beta_ge_cStrVec,
  seBeta_ge_cStrVec,
  pval_ge_cStrVec,
  pval_noSPA_ge_cStrVec,
  pval_SKAT_ge_cVec,
  itt,
  q
);

}





// [[Rcpp::export]]
Rcpp::List  getOneMarkerID_VCF(
                               std::string & t_ref,       // REF allele
                               std::string & t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string & t_marker,    // marker ID extracted from genotype file
                               uint32_t & t_pd,           // base position
                               std::string & t_chr       // chromosome
                               )     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
{
     
    bool isReadVariant = ptr_gVCFobj->getOneMarker_ID(t_ref, t_alt, t_marker, t_pd, t_chr);
    ptr_gVCFobj->move_forward_iterator(1);
    Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("REF") = t_ref,
                                          Rcpp::Named("ALT") = t_alt,
                                          Rcpp::Named("ID") = t_marker, 
                                          Rcpp::Named("POS") = t_pd,
                                          Rcpp::Named("CHROM") = t_chr,
					  Rcpp::Named("isReadVariant") = isReadVariant
                                          );


    return(OutList);
}

// [[Rcpp::export]]
void assign_g_outputFilePrefixSingle( std::string t_outputFilePrefixSingle){
	g_outputFilePrefixSingle = t_outputFilePrefixSingle;
}


// [[Rcpp::export]]
void assign_g_outputFilePrefix( std::string t_outputFilePrefix){
	g_outputFilePrefixGroup = t_outputFilePrefix;
	g_outputFilePrefixSingleInGroup = t_outputFilePrefix + ".singleAssoc.txt";
	g_outputFilePrefixSingleInGroup_temp = t_outputFilePrefix + ".singleAssoc.txt_temp";
}

// [[Rcpp::export]]
void assign_g_outputFilePrefix0( std::string t_outputFilePrefix){
	g_outputFilePrefix0 = t_outputFilePrefix;
}
