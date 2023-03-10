
#ifndef NullGENO_HPP
#define NullGENO_HPP


#define ARMA_USE_SUPERLU 1
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <omp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>// include this header for calculating execution time
#include <cassert>
#include <boost/date_time.hpp> // for gettimeofday and timeval
#include "getMem.hpp"
using namespace Rcpp;
using namespace std;
using namespace RcppParallel;

namespace NullGENO {

class NullGenoClass{
private:
        //COPY from RVTEST:
        // we reverse the two bits as defined in PLINK format
        const static unsigned char HOM_REF = 0x0;  // 0b00 ;
        const static unsigned char HET = 0x2;      // 0b10 ;
        const static unsigned char HOM_ALT = 0x3;  // 0b11 ;
        const static unsigned char MISSING = 0x1;  // 0b01 ;

	arma::sp_fmat g_I_longl_mat;
	arma::sp_fmat g_T_longl_mat;
	arma::uvec g_I_longl_vec;
	arma::fvec g_T_longl_vec;

public:
	//to chunk the geno vector to avoid large continuous memory usage
        int numMarkersofEachArray;
        int numofGenoArray;
        int numMarkersofLastArray;
        std::vector< std::vector<unsigned char>* > genoVecofPointers;
        ///////////
        std::vector< std::vector<unsigned char>* > genoVecofPointers_forVarRatio;
        //arma::fvec g_cateVarRatioMinMACVecExclude;
        //arma::fvec g_cateVarRatioMaxMACVecInclude;
        float g_minMACVarRatio;
        float g_maxMACVarRatio;
        bool isVarRatio = false;
        int numberofMarkers_varRatio = 0;
        int numberofMarkers_varRatio_common = 0;
        arma::ivec g_randMarkerIndforVR;
        std::vector<float>      invstdvVec0_forVarRatio;
        arma::fvec      invstdvVec_forVarRatio;
         std::vector<float>      alleleFreqVec0_forVarRatio;
        arma::fvec      alleleFreqVec_forVarRatio;
        std::vector<int>      MACVec0_forVarRatio;
        std::vector<int>      markerIndexVec0_forVarRatio;
        arma::ivec MACVec_forVarRatio;
        arma::ivec markerIndexVec_forVarRatio;
	        //vector<unsigned char> genoVec;
        size_t M;
        size_t N;
        size_t Nnomissing;
        std::vector<float>      invstdvVec0;
        arma::fvec      invstdvVec;
        vector<int>     ptrsubSampleInGeno;
        std::vector<bool> indicatorGenoSamplesWithPheno_in;
        std::vector<float>      alleleFreqVec0;
        arma::fvec      alleleFreqVec;
        arma::ivec      m_OneSNP_Geno;
        arma::fvec      m_OneSNP_StdGeno;
        arma::fvec      m_DiagStd;
        arma::fvec      m_DiagStd_LOCO;
        arma::fmat      mtx_DiagStd_LOCO;


        std::vector<int>        MACVec0; //for variance ratio based on different MAC categories
        arma::ivec      MACVec;
        arma::ivec      subMarkerIndex; //for sparse GRM
        arma::fmat      stdGenoMultiMarkersMat;
        std::vector<float> stdGenoforSamples; //for sparse GRM
        std::vector<float>     kinValueVecFinal;
        float relatednessCutoff;
        float maxMissingRate;
	        float minMAFtoConstructGRM;


        tbb::concurrent_vector< std::pair<int, int> > indiceVec;
        arma::ivec xout;
        arma::ivec yout;
        bool setKinDiagtoOne;
        int numberofMarkerswithMAFge_minMAFtoConstructGRM = 0;
        std::vector<bool> MarkerswithMAFge_minMAFtoConstructGRM_indVec;


        size_t Msub;
        int startIndex;
	int endIndex;
        int chromIndex;


        arma::ivec startIndexVec;
        arma::ivec endIndexVec;
        arma::ivec startIndexVec_forvr;
        arma::ivec endIndexVec_forvr;


        int Msub_MAFge_minMAFtoConstructGRM;

        int Msub_MAFge_minMAFtoConstructGRM_singleChr;
        arma::ivec Msub_MAFge_minMAFtoConstructGRM_byChr;

        unsigned char m_genotype_buffer[4];
        int geno_idx;
        int m_size_of_esi;
        unsigned char m_bits_val[8];

	//look-up table in a 2D array for sparseKin
        float sKinLookUpArr[3][3] = {{0}};

	//look-up table for std geno
        //float stdGenoLookUpArr[3] = {0};
	//
	//
	//NullGenoClass(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool  isDiagofKinSetAsOne);
	NullGenoClass();

	void setStdGenoLookUpArr(float mafVal, float invsdVal, arma::fvec & stdGenoLookUpArr);


        void setSparseKinLookUpArr(float mafVal, float invsdVal);
	void setBit(unsigned char & ch, int ii, int aVal, int bVal);
	void setGenotype(unsigned char* c, const int pos, const int geno);
	void getGenotype(unsigned char* c, const int pos, int& geno);

	void setstartendIndexVec(arma::ivec & t_startIndexVec, arma::ivec & t_endIndexVec);



	void Init_OneSNP_Geno();

	arma::ivec * Get_OneSNP_Geno(size_t SNPIdx);

	arma::ivec * Get_OneSNP_Geno_forVarRatio(size_t SNPIdx);

	void Get_OneSNP_Geno_atBeginning(size_t SNPIdx, vector<int> & indexNA, vector<unsigned char> & genoVecOneMarkerOld, float & altFreq, float & missingRate, int & mac,  int & alleleCount, bool & passQC, size_t SNPIdx_new, bool & passVarRatio , size_t SNPIdx_vr);

	int Get_OneSNP_StdGeno(size_t SNPIdx, arma::fvec * out );

	arma::fvec * Get_Diagof_StdGeno();

	arma::fvec * Get_Diagof_StdGeno_LOCO();

	void setGenoObj(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool  isDiagofKinSetAsOne);

	void printFromgenoVec(unsigned char genoBinary0);

	int getM() const;
	int getnumberofMarkerswithMAFge_minMAFtoConstructGRM() const;

	int getMsub() const;
	int getStartIndex() const;
	int getEndIndex() const;
        int getN() const;
        int getNnomissing() const;
        float getAC(int m);
	float getMAC(int m);
	int getMsub_MAFge_minMAFtoConstructGRM_in() const;
	int getMsub_MAFge_minMAFtoConstructGRM_singleChr_in() const;
	void Get_Samples_StdGeno(arma::ivec SampleIdsVec);
 
	
};

}

#endif	
