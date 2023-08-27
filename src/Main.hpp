#ifndef MAIN_HPP
#define MAIN_HPP

#define ARMA_USE_SUPERLU 1
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


void set_dup_sample_index(arma::uvec & t_dup_sample_Index);

void setAssocTest_GlobalVarsInCPP(std::string t_impute_method,
                double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
                               double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,
                               arma::vec & t_weights_beta,
			       std::string t_outputFilePrefix,
			        double t_MACCutoffforER);



void setAssocTest_GlobalVarsInCPP_GbyE(
				arma::fmat & t_emat,
				bool t_isgxe
					);


void setMarker_GlobalVarsInCPP(
                               bool t_isOutputMoreDetails,
                               int t_marker_chunksize,
			       arma::mat & t_emat,	                                bool t_isgxe
                               );


void setRegion_GlobalVarsInCPP(
                               arma::vec t_max_maf_region,
                               unsigned int t_max_markers_region,
                               double t_MACCutoff_to_CollapseUltraRare,
			       double t_min_gourpmac_for_burdenonly);



void mainMarkerInCPP(
                           std::string & t_genoType,     // "PLINK", "BGEN"
                           std::string & t_traitType,
			   std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string>  & t_genoIndex,
                           bool & t_isMoreOutput,
                           bool & t_isImputation,
			   bool & t_isFirth);

bool Unified_getOneMarker(std::string & t_genoType,   // "PLINK", "BGEN"
				uint64_t & t_gIndex_prev,
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
                               std::vector<uint>& t_indexForNonZero,
			       arma::vec & t_GVec,
			       bool t_isImputation);

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
                                bool t_isSparseGRM); 



Rcpp::List mainRegionInCPP(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           arma::mat & annoIndicatorMat,
                           arma::vec & maxMAFVec,
                           std::string t_outputFile,
                           std::string t_traitType,
                           unsigned int t_n,           // sample size
                           arma::mat & P1Mat,            // edited on 2021-08-19: to avoid repeated memory allocation of P1Mat and P2Mat
                           arma::mat & P2Mat,
                           std::string t_regionTestType,
                           bool t_isImputation,
                           arma::vec & t_weight,
                           arma::vec & t_weight_cond,
                           bool t_isSingleinGroupTest,
                           bool t_isOutputMarkerList,
                           std::vector<std::string> & annoStringVec,
                           std::string regionName,
			   bool t_isFastTest,
                           bool t_isMoreOutput);

void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> & t_SampleInModel,
                      std::string t_AlleleOrder);



void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> & t_SampleInBgen,
                     std::vector<std::string> & t_SampleInModel,
                     std::string t_AlleleOrder);

void setVCFobjInCPP(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            std::vector<std::string> & t_SampleInModel);

void setSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
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


void assign_conditionMarkers_factors(
                           std::string t_genoType,     // "plink", "bgen", "vcf"
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           unsigned int t_n,
                           arma::vec & t_weight_cond
                           );


void assign_conditionMarkers_factors_binary_region(
                           arma::vec & scalefactor_G2_cond);

void set_iterator_inVcf(std::string & variantList);

bool check_Vcf_end();

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
        double xconv);

Rcpp::List RegionSetUpConditional_binary_InCPP(arma::vec & t_weight_cond);

void closeGenoFile(std::string & t_genoType);

bool openOutfile(std::string t_traitType, bool isappend);


bool openOutfile_singleinGroup(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput);

bool openOutfile_single(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput, bool t_isGbyE);

void writeOutfile_single(bool t_isMoreOutput,
      bool t_isImputation,
                        bool t_isCondition,
			bool t_isFirth,
                         int mFirth,
                         int mFirthConverge,
                        std::string t_traitType,
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
                        std::vector<double>  & pval_SKAT_ge_cVec);




int writeOutfile_singleinGroup(bool t_isMoreOutput,
      bool t_isImputation,
                        bool t_isCondition,
                        bool t_isFirth,
                         int mFirth,
                         int mFirthConverge,
                        std::string t_traitType,
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
                        std::vector<uint32_t> & N_Vec);


void set_flagSparseGRM_cur_SAIGE(bool t_flagSparseGRM_cur);

void set_flagSparseGRM_cur_SAIGE_org();


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
                        unsigned int q_anno,
                        unsigned int q_maf,
                        bool isCondition,
                        std::string t_traitType);

void copy_singleInGroup();

void set_varianceRatio(double MAC, bool isSingleVarianceRatio, bool isnoXadj);

int writeOutfile_singleInGroup(bool t_isMoreOutput,
                        bool t_isImputation,
                        bool t_isCondition,
                        bool t_isFirth,
                         int mFirth,
                         int mFirthConverge,
                        std::string t_traitType,
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
                        std::ofstream & t_OutFile_singleInGroup);

uint32_t Unified_getSampleSizeinGeno(std::string & t_genoType);
uint32_t Unified_getSampleSizeinAnalysis(std::string & t_genoType);


void setupSparseGRM_new(arma::sp_mat & t_spGRM);

void set_I_longl_mat(arma::sp_mat & t_Ilongmat, arma::vec & t_I_longl_vec);

void set_T_longl_mat(arma::sp_mat & t_Tlongmat, arma::vec & t_T_longl_vec);

arma::fvec getCrossprodMatAndKin(arma::fcolvec& bVec, bool LOCO);

arma::fcolvec getCrossprod_multiV(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, bool LOCO);

arma::fvec getDiagOfSigma_multiV(arma::fvec& wVec, arma::fvec& tauVec, bool LOCO);

void gen_sp_Sigma_multiV(arma::fvec& wVec,  arma::fvec& tauVec);

arma::fvec getPCG1ofSigmaAndVector_multiV(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG, bool LOCO);


void set_seed(unsigned int seed);

Rcpp::NumericVector nb(int n);

void setStartEndIndex(int startIndex, int endIndex, int chromIndex);

void setStartEndIndexVec( arma::ivec & startIndex_vec,  arma::ivec & endIndex_vec);

float calCV(arma::fvec& xVec);

arma::fmat getSigma_X_multiV(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, int maxiterPCG, float tolPCG, bool LOCO);

arma::fvec  getSigma_G_multiV(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG, bool LOCO);

Rcpp::List fitglmmaiRPCG_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec &wVec,  arma::fvec & tauVec, arma::ivec & fixtauVec, arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff, bool LOCO);

arma::fvec getMeanDiagofKmat(bool LOCO);

arma::ivec updatefixrhoidx0(arma::fvec & t_tau0Vec, float tol);

arma::ivec tauUpdateValue(arma::fvec & t_tau0Vec);

Rcpp::List getAIScore_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, arma::ivec & fixtauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO);


arma::fvec GetTrace_multiV(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::ivec & fixtauVec, arma::fmat& cov1,  int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO);


Rcpp::List getCoefficients_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG, bool LOCO);


void setminMAC_VarianceRatio(float t_minMACVarRatio, float t_maxMACVarRatio, bool t_isVarianceRatioinGeno);


arma::fvec get_GRMdiagVec();


void setminMAFforGRM(float minMAFforGRM);


void setmaxMissingRateforGRM(float maxMissingforGRM);


void set_Diagof_StdGeno_LOCO();


arma::fvec get_DiagofKin();

arma::fvec parallelCrossProd_usingSubMarker(arma::fcolvec & bVec);

arma::fvec getCrossprodMatAndKin_usingSubMarker(arma::fcolvec& bVec);

float parallelInnerProduct(std::vector<float> &x, std::vector<float> &y);

Rcpp::List createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec);


Rcpp::List refineKin(float relatednessCutoff);

arma::fmat getColfromStdGenoMultiMarkersMat(arma::uvec & a);

int getNColStdGenoMultiMarkersMat();

int getNRowStdGenoMultiMarkersMat();

void setSubMarkerIndex(arma::ivec &subMarkerIndexRandom);

void setRelatednessCutoff(float a);

double innerProduct(Rcpp::NumericVector x, Rcpp::NumericVector y) ;

arma::fvec getDiagOfSigma_noV(arma::fvec& wVec, arma::fvec& tauVec, bool LOCO);

arma::fcolvec getCrossprod_noV(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, bool LOCO);

arma::fvec getPCG1ofSigmaAndVector_noV(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG, bool LOCO);

arma::fvec  getSigma_G_noV(arma::fvec& wVec, arma::fvec& tauVec,arma::fvec& Gvec, int maxiterPCG, float tolPCG, bool LOCO);


void set_useGRMtoFitNULL(bool useGRMtoFitNULL);

void set_isSparseGRM(bool t_isSparseGRM);

void set_store_sigma(bool isstoreSigma);

void set_num_Kmat(int t_num_Kmat);

int get_numofV();

arma::umat set_covarianceidx_Mat();

//void set_Vmat_vec_longlVar();

void closeGenoFile_plink();

int gettotalMarker();

arma::fvec getAlleleFreqVec();

arma::ivec getMACVec();

arma::ivec getMACVec_forVarRatio();

arma::ivec getIndexVec_forVarRatio();

bool getIsVarRatioGeno();

arma::ivec getSubMarkerIndex();

std::vector<bool> getQCdMarkerIndex();

int getSubMarkerNum();

void initKinValueVecFinal(int ni);

int getNnomissingOut();

int getMsub_MAFge_minMAFtoConstructGRM();

int getMsub_MAFge_minMAFtoConstructGRM_singleChr();

void Get_MultiMarkersBySample_StdGeno_Mat();

void Get_MultiMarkersBySample_StdGeno(arma::fvec& markerIndexVec, std::vector<float> &stdGenoMultiMarkers);


arma::fvec parallelCrossProd(arma::fcolvec & bVec);

float innerProductFun(std::vector<float> &x, std::vector<float> & y);

arma::fvec parallelCrossProd_full(arma::fcolvec & bVec, int & markerNum);

arma::fvec parallelCrossProd_LOCO(arma::fcolvec & bVec);

void printComb(int N);

void findIndiceRelatedSample();

void parallelcalsparseGRM(arma::fvec &GRMvec);

void  parallelsumTwoVec(arma::fvec &x);

void setgenoNULL();

void setgeno(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool isDiagofKinSetAsOne);

arma::ivec Get_OneSNP_Geno(int SNPIdx);

arma::ivec Get_OneSNP_Geno_forVarRatio(int SNPIdx);

arma::fvec Get_OneSNP_StdGeno(int SNPIdx);

arma::fvec  getSigma_G_V(arma::fvec& wVec, float tauVal, float tauVal0, arma::fvec& Gvec, int maxiterPCG, float tolPCG);

arma::fvec getPCG1ofSigmaAndVector_V(arma::fvec& wVec,  float tauVal, float tauVal0, arma::fvec& bVec, int maxiterPCG, float tolPCG);

arma::fcolvec getCrossprod_V(arma::fcolvec& bVec, arma::fvec& wVec, float tauVal, float tauVal0);

arma::fvec getDiagOfSigma_V(arma::fvec& wVec, float tauVal, float tauVal0);

void set_T_longl_mat_SAIGEtest(arma::sp_mat & t_Tlongmat, arma::vec & t_T_longl_vec);

void set_I_longl_mat_SAIGEtest(arma::sp_mat & t_Ilongmat, arma::vec & t_I_longl_vec);

void  get_indexinAnotherVector(std::vector<uint> & nonzeroInd_orig, arma::uvec & dupInd, arma::uvec &  nonzeroInd_new_arma);

arma::sp_fmat get_sp_Sigma_to_R();

std::string join(std::vector<std::string> const &strings, std::string delim);

Rcpp::List get_SKAT_pvalue_Rcpp(arma::vec & Score, arma::mat &  Phi, arma::vec & r_corr);



void mainMarkerInCPP_multi(
                           std::string & t_genoType,     // "PLINK", "BGEN"
                           std::string & t_traitType,
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           bool & t_isMoreOutput,
                           bool & t_isImputation,
                           bool & t_isFirth);

Rcpp::List  getOneMarkerID_VCF(
                               std::string & t_ref,       // REF allele
                               std::string & t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string & t_marker,    // marker ID extracted from genotype file
                               uint32_t & t_pd,           // base position
                               std::string & t_chr       // chromosome
                               );

#endif
