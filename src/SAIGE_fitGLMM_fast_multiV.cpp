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


float minMAFtoConstructGRM = 0;
//This is a class with attritbutes about the genotype informaiton 
class genoClass{
private:
        //COPY from RVTEST:
        // we reverse the two bits as defined in PLINK format
        const static unsigned char HOM_REF = 0x0;  // 0b00 ;
        const static unsigned char HET = 0x2;      // 0b10 ;
        const static unsigned char HOM_ALT = 0x3;  // 0b11 ;
        const static unsigned char MISSING = 0x1;  // 0b01 ;


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
	std::vector<float>	invstdvVec0;
	arma::fvec	invstdvVec;
	vector<int>	ptrsubSampleInGeno;
	std::vector<bool> indicatorGenoSamplesWithPheno_in;	
	

  	std::vector<float> 	alleleFreqVec0;
  	arma::fvec 	alleleFreqVec;
  	arma::ivec	m_OneSNP_Geno;
  	arma::fvec	m_OneSNP_StdGeno;
  	arma::fvec	m_DiagStd;
	arma::fvec	m_DiagStd_LOCO;
  	arma::fmat	mtx_DiagStd_LOCO;


	std::vector<int>	MACVec0; //for variance ratio based on different MAC categories
	arma::ivec	MACVec;
	arma::ivec	subMarkerIndex; //for sparse GRM
	arma::fmat      stdGenoMultiMarkersMat;	
	std::vector<float> stdGenoforSamples; //for sparse GRM
	std::vector<float>     kinValueVecFinal;
        float relatednessCutoff;
	float maxMissingRate;

	tbb::concurrent_vector< std::pair<int, int> > indiceVec;
	arma::ivec xout;
        arma::ivec yout;
	//int Mmafge1perc;
	bool setKinDiagtoOne;
	int numberofMarkerswithMAFge_minMAFtoConstructGRM = 0;
//	arma::SpMat<float> sparseGRMinC(2,2);
	std::vector<bool> MarkerswithMAFge_minMAFtoConstructGRM_indVec;	


        //std::vector<float> stdGenoVec;
	//for LOCO
	//bool LOCO = false;
	//vector<int> chromosomeStartIndex;
	//vector<int> chromosomeEndIndex;
	//vector<int> chromosomeVec;
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
	//end for LOCO

	unsigned char m_genotype_buffer[4];
	int geno_idx;
	int m_size_of_esi;
	unsigned char m_bits_val[8];

	
	//look-up table for std geno
	//float stdGenoLookUpArr[3] = {0};
	void setStdGenoLookUpArr(float mafVal, float invsdVal, arma::fvec & stdGenoLookUpArr){
	//	arma::fvec stdGenoLookUpArr(3);
		float mafVal2 = 2*mafVal;
		stdGenoLookUpArr(0) = (0-mafVal2)*invsdVal;
		stdGenoLookUpArr(1) = (1-mafVal2)*invsdVal;
		stdGenoLookUpArr(2) = (2-mafVal2)*invsdVal;
	//	return(stdGenoLookUpArr)
	}


        //look-up table in a 2D array for sparseKin 
        float sKinLookUpArr[3][3] = {{0}};
	//(g - 2*freq)* invStd;;
        void setSparseKinLookUpArr(float mafVal, float invsdVal){
		float mafVal2 = 2*mafVal;
		float a0 = (0-mafVal2)*invsdVal;
		float a1 = (1-mafVal2)*invsdVal;
		float a2 = (2-mafVal2)*invsdVal;
		
		sKinLookUpArr[0][0] = a0*a0;
		sKinLookUpArr[0][1] = a0*a1;
		sKinLookUpArr[0][2] = a0*a2;
		sKinLookUpArr[1][0] = sKinLookUpArr[0][1];
		sKinLookUpArr[1][1] = a1*a1;
		sKinLookUpArr[1][2] = a1*a2;
		sKinLookUpArr[2][0] = sKinLookUpArr[0][2];
		sKinLookUpArr[2][1] = sKinLookUpArr[1][2];
		sKinLookUpArr[2][2] = a2*a2;

	}




        void setBit(unsigned char & ch, int ii, int aVal, int bVal){

                if (bVal == 1 && aVal == 1){
			ch ^= char(1 << ((ii*2) + 1)); //set a to be 1

                }else if(bVal == 0){
			ch ^= char(1 << (ii*2)); //change b to 0

                        if(aVal == 1){
				ch ^= char(1 << ((ii*2) + 1)); //change a to 1
                        }
                }
        }



	//COPY from RVTEST:
	void setGenotype(unsigned char* c, const int pos, const int geno) {
    		(*c) |= (geno << (pos << 1));
  	}

	void getGenotype(unsigned char* c, const int pos, int& geno) {
    		geno = ((*c) >> (pos << 1)) & 0x3;  // 0b11 = 0x3
  	}



	void Init_OneSNP_Geno(){
		m_size_of_esi = (Nnomissing+3)/4;
		int k = 8;
		while (k > 0){
			-- k;
			m_bits_val[k] = 1 << k;
		}
	}
	

        arma::ivec * Get_OneSNP_Geno(size_t SNPIdx){
                m_OneSNP_Geno.zeros(Nnomissing);

		//avoid large continuous memory usage
		int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
		////////////////

                size_t Start_idx = m_size_of_esi * SNPIdxinVec;
                size_t ind= 0;
                unsigned char geno1;
                int bufferGeno;
                for(size_t i=Start_idx; i< Start_idx+m_size_of_esi - 1; i++){
                        //geno1 = genoVec[i];
			geno1 = genoVecofPointers[indexOfVectorPointer]->at(i); //avoid large continuous memory usage
                        for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
				bufferGeno = 2-(a+b);
				m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                //if(ind >= Nnomissing){

                                ////printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                //        return & m_OneSNP_Geno;
                                //}
                        }
                }

		size_t i = Start_idx+m_size_of_esi - 1;
		geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);
		for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                bufferGeno = 2-(a+b);
                                m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                if(ind >= Nnomissing){

                                ////printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                        return & m_OneSNP_Geno;
                                }
                }

                return & m_OneSNP_Geno;
       }
   
        arma::ivec * Get_OneSNP_Geno_forVarRatio(size_t SNPIdx){
                m_OneSNP_Geno.zeros(Nnomissing);

		//avoid large continuous memory usage
		int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
		////////////////

                size_t Start_idx = m_size_of_esi * SNPIdxinVec;
                size_t ind= 0;
                unsigned char geno1;
                int bufferGeno;
                for(size_t i=Start_idx; i< Start_idx+m_size_of_esi-1; i++){
                        //geno1 = genoVec[i];
			geno1 = genoVecofPointers_forVarRatio[indexOfVectorPointer]->at(i); //avoid large continuous memory usage
                        for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
				bufferGeno = 2-(a+b);
				m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                //if(ind >= Nnomissing){

                                //printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                //        return & m_OneSNP_Geno;
                                //}
                        }
                }

		size_t i = Start_idx+m_size_of_esi-1;
		geno1 = genoVecofPointers_forVarRatio[indexOfVectorPointer]->at(i); //avoid large continuous memory usage
                for(int j=0; j<4; j++){
                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                bufferGeno = 2-(a+b);
                                m_OneSNP_Geno[ind] = bufferGeno;
                                ind++;
                                geno1 = geno1 >> 1;
                                if(ind >= Nnomissing){

                                //printf("%d, %d-%d-%d-%f-%d\n",Start_idx, genoVec[i] ,a ,b , m_OneSNP_Geno[ind-1] , m_size_of_esi);
                                        return & m_OneSNP_Geno;
                                }
                  }

                return & m_OneSNP_Geno;
       }


	void Get_OneSNP_Geno_atBeginning(size_t SNPIdx, vector<int> & indexNA, vector<unsigned char> & genoVecOneMarkerOld, float & altFreq, float & missingRate, int & mac,  int & alleleCount, bool & passQC, size_t SNPIdx_new, bool & passVarRatio , size_t SNPIdx_vr){

		arma::ivec m_OneSNP_GenoTemp;
		m_OneSNP_GenoTemp.zeros(N);
		m_OneSNP_Geno.zeros(Nnomissing);
		int m_size_of_esi_temp = (N+3)/4;
		size_t ind= 0;
		unsigned char geno1;
		int bufferGeno;
		int u;
		alleleCount = 0;
		int numMissing = 0;
		for(int i=0; i< (m_size_of_esi_temp - 1); i++){
			geno1 = genoVecOneMarkerOld[i];
			for(int j=0; j<4; j++){
				u = j & 3;

				int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                if (b == 1 && a == 0){
                                        bufferGeno = 3;
                                }else if(b == 0 && a == 0){
                                        bufferGeno = 2;
                                }else if(b == 0 && a == 1){
                                        bufferGeno = 1;
                                }else if(b == 1 && a == 1){
                                        bufferGeno = 0;
                                }else{
                                        cout << "Error GENO!!\n";
                                        break;
                                }


				//getGenotype(&geno1, u, bufferGeno);
				//printf("%d", geno1);
	//			std::cout << "bufferGeno " << bufferGeno << std::endl;
				/*
				switch(geno1){
    				 case HOM_REF: break;
    				 case HET: sum+=1; break;
    				 case HOM_ALT: sum+=2; break;
    				 case MISSING: numMissing++; break;
    				}
				*/
				m_OneSNP_GenoTemp[ind] = bufferGeno;
				//if(SNPIdx == 0){
				//	std::cout << "[ind] " << ind << std::endl;
				//	std::cout << "indicatorGenoSamplesWithPheno_in[ind] " << indicatorGenoSamplesWithPheno_in[ind] << std::endl;
				//}	
				if(indicatorGenoSamplesWithPheno_in[ind]){
					if(bufferGeno == 3){
	//					std::cout << "SNPIdx " << SNPIdx << std::endl;
	//					std::cout << "ind " << ind << std::endl;
						numMissing++;
						//indexNA.push_back(ptrsubSampleInGeno[ind])
					}else{
						alleleCount = alleleCount + bufferGeno;
					}	
				}	
				ind++;
                                geno1 = geno1 >> 1;	
			  }

	      }	  

		int i = m_size_of_esi_temp - 1;
		geno1 = genoVecOneMarkerOld[i];
		//std::cout << "N " << N << std::endl;
		//while(ind < N){
                        for(int j=0; j<4; j++){
                                u = j & 3;

                                int b = geno1 & 1 ;
                                geno1 = geno1 >> 1;
                                int a = geno1 & 1 ;
                                if (b == 1 && a == 0){
                                        bufferGeno = 3;
                                }else if(b == 0 && a == 0){
                                        bufferGeno = 2;
                                }else if(b == 0 && a == 1){
                                        bufferGeno = 1;
                                }else if(b == 1 && a == 1){
                                        bufferGeno = 0;
                                }else{
                                        cout << "Error GENO!!\n";
                                        break;
                                }

                                m_OneSNP_GenoTemp[ind] = bufferGeno;
                                //if(SNPIdx == 0){
                                //        std::cout << "[ind] " << ind << std::endl;
                                //        std::cout << "indicatorGenoSamplesWithPheno_in[ind] " << indicatorGenoSamplesWithPheno_in[ind] << std::endl;
                                //}
                                if(indicatorGenoSamplesWithPheno_in[ind]){
                                        if(bufferGeno == 3){
        //                                      std::cout << "SNPIdx " << SNPIdx << std::endl;
        //                                      std::cout << "ind " << ind << std::endl;
                                                numMissing++;
                                        }else{
                                                alleleCount = alleleCount + bufferGeno;
                                        }
                                }
                                ind++;
                                geno1 = geno1 >> 1;
				if(ind == (N)){
				break;
				}
                          }
		//}


	      altFreq = alleleCount/float((Nnomissing-numMissing) * 2);
	      //sum = 0;
	      //std::cout << "missingRate " << missingRate << std::endl;
	      //std::cout << "maxMissingRate " << maxMissingRate << std::endl;
	      missingRate = numMissing/float(Nnomissing);	      

              //int indxInOut = 0;
	      //if(minMAFtoConstructGRM > 0){
              //if(altFreq >= minMAFtoConstructGRM && altFreq <= (1-minMAFtoConstructGRM) && missingRate <= maxMissingRate){
	      	int fillinMissingGeno = int(round(2*altFreq)); 
	
		if(numMissing > 0){
		       //for(int indx=0; indx < Nnomissing; indx++){
                                                //cout << "HERE5\n";
			//	u = indx & 3;
			//	bufferGeno = m_OneSNP_GenoTemp[ptrsubSampleInGeno[indx] - 1];
			alleleCount = alleleCount + fillinMissingGeno*numMissing;
			
			//		if(bufferGeno == 3){
			//			bufferGeno = fillinMissingGeno;
			//			alleleCount = alleleCount + bufferGeno;
			//		}
			//	}	
  			/*
				//setGenotype(&geno2, u, bufferGeno);
				if(bufferGeno == 0){
                                        setGenotype(&geno2, u, HOM_ALT);
                                }else if(bufferGeno == 1){
                                        setGenotype(&geno2, u, HET);
                                }else if(bufferGeno == 2){
                                        setGenotype(&geno2, u, HOM_REF);
                                }
				//else{
                                //        setGenotype(&geno1, u, MISSING);
                                        //m_OneSNP_Geno[j] = 0;  //12-18-2017
                                //}	

				if(u == 3 || indx == (Nnomissing-1)){
                                        genoVecofPointers[SNPIdx/numMarkersofEachArray]->push_back(geno2); //avoid large continuous memory usage
                                        geno2 = 0;
               			}
		*/	
		}
			//passQC = true;	
	     //}

	     altFreq = alleleCount/float(Nnomissing * 2);

	     unsigned char geno2;
	     passQC = false;
	     passVarRatio = false;
	     float maf = std::min(altFreq, 1-altFreq);
	     mac = std::min(alleleCount, int(Nnomissing) * 2 - alleleCount);


		if(maf >= minMAFtoConstructGRM && missingRate <= maxMissingRate){
			passQC = true;
		}
		if(isVarRatio){
			if(g_maxMACVarRatio != -1){ //if estimating categorical variance ratios
			   if(mac >= g_minMACVarRatio && mac < g_maxMACVarRatio){
				passVarRatio = true;
				genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
				genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
			   }else if(mac >= g_maxMACVarRatio){
				   //randomly select 200 markers for estimating the variance ratio for the last MAC category	
				   //if(numberofMarkers_varRatio_common < 200){
				   	//if(static_cast<int>(SNPIdx) == 123){
					//	std::cout << "123" << std::endl;
					//	bool isIng_randMarkerIndforVR = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
					//	std::cout << "isIng_randMarkerIndforVR " << isIng_randMarkerIndforVR << std::endl;
					//}
				   	passVarRatio = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
					if(passVarRatio){
						//std::cout << "SNPIdx " << SNPIdx << std::endl;
						genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
						genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
				  		numberofMarkers_varRatio_common = numberofMarkers_varRatio_common + 1;
						//passVarRatio = false;
					}
				  //} 
			//	passVarRatio = true;	
			}
			}else{
				if(mac >= g_minMACVarRatio){
                                   //randomly select 200 markers for estimating the variance ratio for the last MAC category
                                   //if(numberofMarkers_varRatio_common < 200){
                                        //if(static_cast<int>(SNPIdx) == 123){
                                        //      std::cout << "123" << std::endl;
                                        //      bool isIng_randMarkerIndforVR = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
                                        //      std::cout << "isIng_randMarkerIndforVR " << isIng_randMarkerIndforVR << std::endl;
                                        //}
                                        passVarRatio = arma::any(g_randMarkerIndforVR == static_cast<int>(SNPIdx));
                                        if(passVarRatio){
                                                //std::cout << "SNPIdx " << SNPIdx << std::endl;
                                                genoVecofPointers_forVarRatio[SNPIdx_vr] = new vector<unsigned char>;
                                                genoVecofPointers_forVarRatio[SNPIdx_vr]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                                                numberofMarkers_varRatio_common = numberofMarkers_varRatio_common + 1;
                                                //passVarRatio = false;
                                        }
                                  //}
                        //      passVarRatio = true;
                        }	
			

			}
			//avoid the overlap between markers for GRM and markers for variance ratio estimation	   
			if(passVarRatio){
				passQC = false;
			}
		      //}
		}

		if(passQC | passVarRatio){
			for(int indx=0; indx < Nnomissing; indx++){
                                                //cout << "HERE5\n";
                              u = indx & 3;
                              bufferGeno = m_OneSNP_GenoTemp[ptrsubSampleInGeno[indx] - 1];
			      if(bufferGeno == 3){					
				bufferGeno = fillinMissingGeno;
			      }	
			      if(bufferGeno == 0){
                                        setGenotype(&geno2, u, HOM_ALT);
                                }else if(bufferGeno == 1){
                                        setGenotype(&geno2, u, HET);
                                }else if(bufferGeno == 2){
                                        setGenotype(&geno2, u, HOM_REF);
                              }
			      if(u == 3 || indx == (Nnomissing-1)){
				       if(passVarRatio){
						genoVecofPointers_forVarRatio[SNPIdx_vr/numMarkersofEachArray]->push_back(geno2); //avoid large continuous memory usage
					}
				        if(passQC){
						genoVecofPointers[SNPIdx_new/numMarkersofEachArray]->push_back(geno2);
					}	
                                        geno2 = 0;
                              }
			}
		}
	    // altFreq = alleleCount/float(Nnomissing * 2);

   }
		

   int Get_OneSNP_StdGeno(size_t SNPIdx, arma::fvec * out ){
		//avoid large continuous memory usage
                int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
                int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
                ////////////////
		//std::cout << "indexOfVectorPointer " << indexOfVectorPointer << std::endl;
 		out->zeros(Nnomissing);
		//std::cout << "m_size_of_esi " << m_size_of_esi << std::endl;
		//std::cout << "SNPIdxinVec " << SNPIdxinVec << std::endl;
		//std::cout << "genoVecofPointers[indexOfVectorPointer]->size() " << genoVecofPointers[indexOfVectorPointer]->size() << std::endl;



 		size_t Start_idx = m_size_of_esi * SNPIdxinVec;

		//std::cout << "Start_idx " << Start_idx << std::endl;
		size_t ind= 0;
		unsigned char geno1;
		
		float freq = alleleFreqVec[SNPIdx];
		//cout << "Get_OneSNP_StdGeno here" << endl; 
		float invStd = invstdvVec[SNPIdx];

		arma::fvec stdGenoLookUpArr(3);
		setStdGenoLookUpArr(freq, invStd, stdGenoLookUpArr);
//		std::cout << "freq " << freq << endl;
//		std::cout << "invStd " << invStd << endl;

		//setStdGenoLookUpArr(freq, invStd);
		//std::cout << "stdGenoLookUpArr[0]: " << stdGenoLookUpArr[0] << std::endl;
		//std::cout << "stdGenoLookUpArr[1]: " << stdGenoLookUpArr[1] << std::endl;
		//std::cout << "stdGenoLookUpArr[2]: " << stdGenoLookUpArr[2] << std::endl;
//		cout << "Get_OneSNP_StdGeno here2"  << endl;
		for(size_t i=Start_idx; i< Start_idx+m_size_of_esi-1; i++){
//			geno1 = genoVec[i];
			geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

			for(int j=0; j<4; j++){
    			int b = geno1 & 1 ;
    			geno1 = geno1 >> 1;
    			int a = geno1 & 1 ;
    			//(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;;
    			(*out)[ind] = stdGenoLookUpArr(2-(a+b));
//			std::cout << "a " << a << endl;
//			std::cout << "b " << b << endl;
//			std::cout << "(*out)[ind] " << (*out)[ind] << endl;
			ind++;
    			geno1 = geno1 >> 1;
    			
//    			if(ind >= Nnomissing){
//				cout << "Get_OneSNP_StdGeno " << SNPIdx << endl; 
//				cout << "Nnomissing " << Nnomissing << endl; 
//				stdGenoLookUpArr.clear();
//    				return 1;
//    			}
	    		}
		}


		size_t i = Start_idx+m_size_of_esi-1;
                geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

                for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        (*out)[ind] = stdGenoLookUpArr(2-(a+b));
                        ind++;
                        geno1 = geno1 >> 1;

                        if(ind >= Nnomissing){
                                stdGenoLookUpArr.clear();
                                return 1;
                        }
                }
		//cout << "SNPIdx " << SNPIdx << endl; 

		stdGenoLookUpArr.clear();
		return 1;
				
 	}


	arma::fvec * Get_Diagof_StdGeno(){
	
		arma::fvec * temp = &m_OneSNP_StdGeno;
		// Not yet calculated
		//cout << "size(m_DiagStd)[0] " << size(m_DiagStd)[0] << endl;
		if(size(m_DiagStd)[0] != Nnomissing){
			m_DiagStd.zeros(Nnomissing);
			for(size_t i=0; i< numberofMarkerswithMAFge_minMAFtoConstructGRM; i++){
				//if(alleleFreqVec[i] >= minMAFtoConstructGRM && alleleFreqVec[i] <= 1-minMAFtoConstructGRM){


				Get_OneSNP_StdGeno(i, temp);

				/*if(i == 0){
					cout << "setgeno mark7 " << i <<  endl;
					for(int j=0; j<10; ++j)
					{
                				cout << (*temp)[j] << ' ';
                			}
                			cout << endl;
				}
				*/
				m_DiagStd = m_DiagStd + (*temp) % (*temp);

				//}
		//		std::cout << "i " << i << std::endl;
		//		std::cout << "numberofMarkerswithMAFge_minMAFtoConstructGRM " << numberofMarkerswithMAFge_minMAFtoConstructGRM << std::endl;
			}
		
		}
/*
		std::cout << "test\n";
		for(int i=0; i<10; ++i)
        	{
        	  cout << m_DiagStd[i] << ' ';
        	}
		cout << endl;
*/	
		return & m_DiagStd;
	}

	
 	

	arma::fvec * Get_Diagof_StdGeno_LOCO(){
                //if(size(m_DiagStd_LOCO)[0] != Nnomissing){
		//m_DiagStd_LOCO.zeros(Nnomissing);
                  //      for(size_t i=startIndex; i<= endIndex; i++){
				//if(i < startIndex || i > endIndex){
		//			if(alleleFreqVec[i] >= minMAFtoConstructGRM && alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
                  //                		Get_OneSNP_StdGeno(i, temp);
                    //              		m_DiagStd_LOCO = m_DiagStd_LOCO + (*temp) % (*temp);
		//				Msub_MAFge_minMAFtoConstructGRM = Msub_MAFge_minMAFtoConstructGRM + 1;
		//			}
				//}
                 //       }


		//m_DiagStd_LOCO = m_DiagStd - geno.mtx_DiagStd_LOCO.col(geno.chromIndex);
		m_DiagStd_LOCO = mtx_DiagStd_LOCO.col(chromIndex);
                Msub_MAFge_minMAFtoConstructGRM_singleChr  = Msub_MAFge_minMAFtoConstructGRM_byChr(chromIndex); 
                //}

                return & m_DiagStd_LOCO;
        }


 
  	//Function to assign values to all attributes
 
  	//Function to assign values to all attributes
  	//This function is used instead of using a constructor because using constructor can not take
  	//genofile as an argument from runModel.R 
        //genofile is the predix for plink bim, bed, fam, files   
  	void setGenoObj(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool  isDiagofKinSetAsOne){
		//cout << "OK1\n";
		setKinDiagtoOne = isDiagofKinSetAsOne;   
		ptrsubSampleInGeno = subSampleInGeno;
		indicatorGenoSamplesWithPheno_in = indicatorGenoSamplesWithPheno;
		Nnomissing = subSampleInGeno.size(); 
    		// reset
    		//genoVec.clear();
    		alleleFreqVec.clear();
		MACVec.clear();
  		invstdvVec.clear();

   		M=0;
  		N=0;
   	
		//std::string bedfile = genofile+".bed";
		//std::string bimfile = genofile+".bim"; 
		//std::string famfile = genofile+".fam"; 
		std::string junk;
		//cout << "OK2\n";
		//count the number of individuals
		ifstream test_famfile;
		test_famfile.open(famfile.c_str());
        	if (!test_famfile.is_open()){
                	printf("Error! fam file not open!");
                	return ;
        	}
		int indexRow = 0;
		while (std::getline(test_famfile,junk)){
                	indexRow = indexRow + 1;
                	junk.clear();
        	}
		N = indexRow;
		test_famfile.clear();	
		//cout << "OK3\n";
		//count the number of markers
		ifstream test_bimfile;
        	test_bimfile.open(bimfile.c_str());
        	if (!test_bimfile.is_open()){
                	printf("Error! bim file not open!");
                	return ;
        	}
        	indexRow = 0;
        	while (std::getline(test_bimfile,junk)){
                	indexRow = indexRow + 1;
                	junk.clear();
        	}
        	M = indexRow;
        	test_bimfile.clear(); 

    		junk.clear();
		//cout << "OK3b\n";
    		// Init OneSNP Geno
    		Init_OneSNP_Geno();
		//cout << "OK3c\n";
    
    		//std::string junk;
    		indexRow = 0;
    		int buffer;
    		int TotalRead=0;

		std::vector<unsigned char> genoVecOneMarkerOld;
		std::vector<unsigned char> genoVecOneMarkerNew;
		/////////////////////////////
		// Added for reserve for genoVec
		size_t nbyteOld = ceil(float(N)/4);
		size_t nbyteNew = ceil(float(Nnomissing)/4);
		size_t reserve = ceil(float(Nnomissing)/4) * M + M*2;
		cout << "nbyte: " << nbyteOld << endl;
		cout << "nbyte: " << nbyteNew << endl;		
		cout << "reserve: " << reserve << endl;		

    		genoVecOneMarkerOld.reserve(nbyteOld);
    		genoVecOneMarkerOld.resize(nbyteOld);
 		//genoVec.reserve(reserve);
		
		//cout << "OK4\n";

		ifstream test_bedfile;
        	test_bedfile.open(bedfile.c_str(), ios::binary);
        	if (!test_bedfile.is_open()){
                	printf("Error! file open!");
                	return;
        	}		
		//printf("\nM: %zu, N: %zu, Reserve: %zu\n", M, N, reserveTemp);
		printf("\nM: %zu, N: %zu\n", M, N);

        	//test_bedfile.seekg(3);
		//set up the array of vectors for genotype
		//numMarkersofEachArray = floor((memoryChunk*pow (10.0, 9.0))/(ceil(float(N)/4)));
		//cout << "numMarkersofEachArray: " << numMarkersofEachArray << endl;
		numMarkersofEachArray = 1;
		//if(M % numMarkersofEachArray == 0){
                        //numofGenoArray = M / numMarkersofEachArray;
                        numofGenoArray = M;
			genoVecofPointers.resize(numofGenoArray);
			genoVecofPointers_forVarRatio.resize(numofGenoArray);
                        //cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;
                        for (int i = 0; i < numofGenoArray ; i++){
                                genoVecofPointers[i] = new vector<unsigned char>;
                                genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                        }

                //}
		/*else{
                        numofGenoArray = M/numMarkersofEachArray + 1;
                        genoVecofPointers.resize(numofGenoArray);
			numMarkersofLastArray = M - (numofGenoArray-1)*numMarkersofEachArray;
                        cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;
			try{	
                        for (int i = 0; i < numofGenoArray-1; i++){
                                genoVecofPointers[i] = new vector<unsigned char>;
                                genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                        }
			genoVecofPointers[numofGenoArray-1] = new vector<unsigned char>;
			genoVecofPointers[numofGenoArray-1]->reserve(numMarkersofLastArray*ceil(float(Nnomissing)/4));
			}
			catch(std::bad_alloc& ba)
                        {
                                std::cerr << "bad_alloc caught1: " << ba.what() << '\n';
                                exit(EXIT_FAILURE);
                        }
		}*/

		cout << "setgeno mark1" << endl;
		arma::ivec g_randMarkerIndforVR_temp;
		//randomly select common markers for variance ratio
		if(isVarRatio){
			 g_randMarkerIndforVR_temp = arma::randi(1000, arma::distr_param(0,M-1));
			 g_randMarkerIndforVR = arma::unique(g_randMarkerIndforVR_temp);
			 //arma::ivec g_randMarkerIndforVR_sort = arma::sort(g_randMarkerIndforVR);
			 //g_randMarkerIndforVR_sort.print("g_randMarkerIndforVR_sort");
		}
		//alleleFreqVec.zeros(M);
		//invstdvVec.zeros(M);
		//MACVec.zeros(M);
        	float freq, Std, invStd, missingRate;
        	int alleleCount, mac;
		std::vector<int> indexNA;
        	int lengthIndexNA;
        	int indexGeno;
        	int indexBit;
        	int fillinMissingGeno;
        	int b2;
        	int a2;

		size_t ind= 0;
                unsigned char geno1 = 0;
                int bufferGeno;
                int u;
		//std::vector<int> genoVec4Markers(4);
		//test_bedfile.read((char*)(&genoVecTemp[0]),nbyteTemp*M);
		bool isPassQC = false;
		bool isPass_vr = false;
		cout << "setgeno mark2" << endl;
		size_t SNPIdx_new = 0;
		size_t SNPIdx_vr = 0;
		//Mmafge1perc = 0;
		for(int i = 0; i < M; i++){
			genoVecOneMarkerOld.clear();
			genoVecOneMarkerOld.reserve(nbyteOld);
                        genoVecOneMarkerOld.resize(nbyteOld);

			test_bedfile.seekg(3+nbyteOld*i);
			test_bedfile.read((char*)(&genoVecOneMarkerOld[0]),nbyteOld);
 			//printf("\nFile read is done: M: %zu, N: %zu, TotalByte %zu\n", M, N, genoVecTemp.size());
			//cout << "Imputing missing genotypes and extracting the subset of samples with nonmissing genotypes and phenotypes\n";  
	//		cout << "i is " << i << endl;  

      			indexNA.clear();
		//}	
        		Get_OneSNP_Geno_atBeginning(i, indexNA, genoVecOneMarkerOld, freq, missingRate, mac, alleleCount, isPassQC, SNPIdx_new, isPass_vr, SNPIdx_vr);

			//std::cout << "freq " << freq << std::endl;
			//std::cout << "isPassQC " << isPassQC << std::endl;
			if(isPassQC){

      				Std = std::sqrt(2*freq*(1-freq));
      				if(Std == 0){
      					invStd= 0;
      				} else {
      					invStd= 1/Std;
      				}
      				invstdvVec0.push_back(invStd);
				alleleFreqVec0.push_back(freq);
				numberofMarkerswithMAFge_minMAFtoConstructGRM = numberofMarkerswithMAFge_minMAFtoConstructGRM + 1;
		
				MACVec0.push_back(mac);	
				MarkerswithMAFge_minMAFtoConstructGRM_indVec.push_back(true);
				SNPIdx_new = SNPIdx_new + 1;
			}else{
				MarkerswithMAFge_minMAFtoConstructGRM_indVec.push_back(false);
			}

			if(isVarRatio){
				if(isPass_vr){
					Std = std::sqrt(2*freq*(1-freq));
                                	if(Std == 0){
                                        	invStd= 0;
                                	}else {
                                        	invStd= 1/Std;
                                	}
					invstdvVec0_forVarRatio.push_back(invStd);
					alleleFreqVec0_forVarRatio.push_back(freq);
					MACVec0_forVarRatio.push_back(mac);
					markerIndexVec0_forVarRatio.push_back(i);
					SNPIdx_vr = SNPIdx_vr + 1;
					numberofMarkers_varRatio = numberofMarkers_varRatio + 1;
				}
			}	


			//m_OneSNP_Geno.clear();

    		}//end for(int i = 0; i < M; i++){

		if( minMAFtoConstructGRM > 0 | maxMissingRate < 1){
			cout << numberofMarkerswithMAFge_minMAFtoConstructGRM << " markers with MAF >= " << minMAFtoConstructGRM << " and missing rate <= " << maxMissingRate  << endl;
		}
		//else{
		//	cout << M << " markers with MAF >= " << minMAFtoConstructGRM << endl;
		//}

		
		int numofGenoArray_old = numofGenoArray;
		if(numberofMarkerswithMAFge_minMAFtoConstructGRM % numMarkersofEachArray == 0){
                        numofGenoArray = numberofMarkerswithMAFge_minMAFtoConstructGRM / numMarkersofEachArray;
                        //genoVecofPointers.resize(numofGenoArray);
                        //cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;

                }else{
                        numofGenoArray = numberofMarkerswithMAFge_minMAFtoConstructGRM/numMarkersofEachArray + 1;
		}
//		cout << " numofGenoArray "<< numofGenoArray << endl;
//		cout << " numofGenoArray_old "<< numofGenoArray_old << endl;
//		cout << " genoVecofPointers.size() "<< genoVecofPointers.size() << endl;
		if(numofGenoArray > numofGenoArray_old){

                        for (int i = numofGenoArray; i < numofGenoArray_old ; i++){
				delete genoVecofPointers[i];
                                //genoVecofPointers[i] = new vector<unsigned char>;
                                //genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                        }
		}
//		cout << " genoVecofPointers.size() "<< genoVecofPointers.size() << endl;

		//genoVecofPointers.resize(MarkerswithMAFge_minMAFtoConstructGRM_indVec);

		invstdvVec.clear();
		invstdvVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);
		alleleFreqVec.clear();
		alleleFreqVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);
		MACVec.clear();
		MACVec.set_size(numberofMarkerswithMAFge_minMAFtoConstructGRM);

		for(int i = 0; i < numberofMarkerswithMAFge_minMAFtoConstructGRM; i++){

			invstdvVec[i] = invstdvVec0.at(i);
			alleleFreqVec[i] = alleleFreqVec0.at(i);
			MACVec[i] = MACVec0.at(i);

		}
	if(isVarRatio){
		invstdvVec_forVarRatio.clear();
                invstdvVec_forVarRatio.set_size(numberofMarkers_varRatio);
		alleleFreqVec_forVarRatio.clear();
                alleleFreqVec_forVarRatio.set_size(numberofMarkers_varRatio);
		MACVec_forVarRatio.clear();
		MACVec_forVarRatio.set_size(numberofMarkers_varRatio);
	        markerIndexVec_forVarRatio.clear();
	        markerIndexVec_forVarRatio.set_size(numberofMarkers_varRatio);	
		for(int i = 0; i < numberofMarkers_varRatio; i++){
			invstdvVec_forVarRatio[i] = invstdvVec0_forVarRatio.at(i);
			alleleFreqVec_forVarRatio[i] =alleleFreqVec0_forVarRatio.at(i);
			MACVec_forVarRatio[i] = MACVec0_forVarRatio.at(i);
			markerIndexVec_forVarRatio[i] = markerIndexVec0_forVarRatio.at(i);
		}
	}

        	test_bedfile.close();
//		cout << "setgeno mark5" << endl;
//		printAlleleFreqVec();
		//printGenoVec();
   		//Get_Diagof_StdGeno();
//		cout << "setgeno mark6" << endl;
  	}//End Function
 

  	void printFromgenoVec(unsigned char genoBinary0){
		unsigned char genoBinary = genoBinary0;
  		for(int j=0; j<4; j++){
        		int b = genoBinary & 1 ;
                	genoBinary = genoBinary >> 1;
                	int a = genoBinary & 1 ;
			genoBinary = genoBinary >> 1;
			cout << 2-(a+b) << " " << endl;
		}
		cout << endl;
  	}
 
  
  	int getM() const{
    		return(M);
  	}

	int getnumberofMarkerswithMAFge_minMAFtoConstructGRM() const{
		return(numberofMarkerswithMAFge_minMAFtoConstructGRM);
	}
 
	//int getMmafge1perc() const{
	//	return(Mmafge1perc);
 	//}

	int getMsub() const{
                return(Msub);
        }

	int getStartIndex() const{
		return(startIndex);
	}

	int getEndIndex() const{
                return(endIndex);
        }

  	int getN() const{
    		return(N);
  	}
 
  	int getNnomissing() const{
    		return(Nnomissing);
  	}


  	float getAC(int m){
    		return(alleleFreqVec[m]*2*Nnomissing);
  	}

  	float getMAC(int m){
    		if(alleleFreqVec[m] > 0.5){
      			return((1-alleleFreqVec[m])*2*Nnomissing);
    		}else{
      			return(alleleFreqVec[m]*2*Nnomissing);
    		}
  	}

	int getMsub_MAFge_minMAFtoConstructGRM_in() const{
		return(numberofMarkerswithMAFge_minMAFtoConstructGRM);	
	}

	int getMsub_MAFge_minMAFtoConstructGRM_singleChr_in() const{
		return(Msub_MAFge_minMAFtoConstructGRM_singleChr);	
	}


	//int getnumberofMarkerswithMAFge_minMAFtoConstructGRM(){
 	//	return(numberofMarkerswithMAFge_minMAFtoConstructGRM);
	//}
  	//print out the vector of genotypes
  	void printGenoVec(){
    		//for(unsigned int i=0; i<M; ++i)
    		for(unsigned int i=0; i<2; ++i)
    		{
	
    			Get_OneSNP_Geno(i);
    			//for(unsigned int j=0; j< Nnomissing; j++){
    			for(unsigned int j=0; j< 100; j++){
      				cout << m_OneSNP_Geno[j] << ' ';
      			}
      			cout << endl;
    		}
    		//cout << "genoVec.size()" << genoVec.size() << endl;
    		cout << "M = " << M << endl;
    		cout << "N = " << N << endl;
  	}
  
  	//print out the vector of allele frequency
  	void printAlleleFreqVec(){
    		//for(int i=0; i<alleleFreqVec.size(); ++i)
    		for(int i=(M-100); i<M; ++i)
    		{
      			cout << alleleFreqVec[i] << ' ';
    		}
    		cout << endl;
  	}


	void Get_Samples_StdGeno(arma::ivec SampleIdsVec){
        	int indexOfVectorPointer;
        	int SNPIdxinVec;

        	int numSamples = SampleIdsVec.n_elem;
        	//stdGenoVec.zeros(Nnomissing*numSamples);
        	stdGenoforSamples.clear();
        	stdGenoforSamples.resize(M*numSamples);

        	arma::ivec sampleGenoIdxVec;
        	sampleGenoIdxVec.zeros(numSamples);
        	arma::ivec sampleGenoIdxSubVec;
        	sampleGenoIdxSubVec.zeros(numSamples);

        	for(int j=0; j < numSamples; j++){
                	sampleGenoIdxVec[j] = SampleIdsVec[j] / 4;
                	sampleGenoIdxSubVec[j] = SampleIdsVec[j] % 4;
        	}


        	int startidx;
        	unsigned char geno1;

        	for(int i=0; i < M; i++){
                	indexOfVectorPointer = i/numMarkersofEachArray;
                	SNPIdxinVec = i % numMarkersofEachArray;
                	startidx = m_size_of_esi * SNPIdxinVec;

                	float freq = alleleFreqVec[i];
                	float invStd = invstdvVec[i];

                	for(int j=0; j < numSamples; j++){
                        	int k = startidx + sampleGenoIdxVec[j];
                        	geno1 = genoVecofPointers[indexOfVectorPointer]->at(k);
                        	for(int q=0; q<4; q++){
                                	if(q == sampleGenoIdxSubVec[j]){
                                        	int b = geno1 & 1 ;
                                        	geno1 = geno1 >> 1;
                                        	int a = geno1 & 1 ;
                                        	stdGenoforSamples[i*(numSamples)+j] = ((2-(a+b)) - 2*freq)* invStd;
                                //(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;;
                                //ind++;
                                        	geno1 = geno1 >> 1;
                                	}else{
                                        	geno1 = geno1 >> 1;
                                        	geno1 = geno1 >> 1;
                                	}
                        	}
                	}
        	}

        //return(stdGenoVec);
	}



  
};



// //create a geno object as a global variable
genoClass geno;

std::vector<arma::sp_fmat> Kmat_vec;
arma::fvec g_longl_vec;

arma::sp_fmat g_I_longl_mat;
arma::sp_fmat g_T_longl_mat;
arma::uvec g_I_longl_vec;
arma::fvec g_T_longl_vec;

// [[Rcpp::export]]
void set_I_longl_mat(arma::sp_mat & t_Ilongmat, arma::vec & t_I_longl_vec){
	arma::sp_fmat t_Kmat_new = arma::conv_to< arma::sp_fmat >::from(t_Ilongmat);
	g_I_longl_mat = t_Kmat_new;
	arma::uvec t_I_longl_vec_new = arma::conv_to< arma::uvec >::from(t_I_longl_vec);
	g_I_longl_vec = t_I_longl_vec_new;

}


// [[Rcpp::export]]
void set_T_longl_mat(arma::sp_mat & t_Tlongmat, arma::vec & t_T_longl_vec){
        arma::sp_fmat t_Kmat_new = arma::conv_to< arma::sp_fmat >::from(t_Tlongmat);
        g_T_longl_mat = t_Kmat_new;
	arma::fvec t_T_longl_vec_new = arma::conv_to< arma::fvec >::from(t_T_longl_vec);
	g_T_longl_vec = t_T_longl_vec_new;
}


// [[Rcpp::export]]
void addNewKat( arma::sp_mat & t_Kmat){
	arma::sp_fmat t_Kmat_new = arma::conv_to< arma::sp_fmat >::from(t_Kmat);
        //Kmat_vec.push_back(t_Kmat); 	
        Kmat_vec.push_back(t_Kmat_new); 	
        std::cout << "Kmat_vec.size() " << Kmat_vec.size() << std::endl;
}

// [[Rcpp::export]]
arma::sp_fmat getProdTauKmat(arma::fvec & tauVec){
	arma::sp_fmat Kmat_sigma;
	int n = geno.getNnomissing();
	Kmat_sigma.set_size(n, n);
	Kmat_sigma.zeros();
	if(Kmat_vec.size() > 0){
		for (int i = 0; i < Kmat_vec.size(); i++){
			Kmat_sigma = Kmat_vec[i] * tauVec(i);
		}
	}	
}	


arma::sp_fmat g_spGRM;
arma::sp_fmat g_spSigma;
bool g_isStoreSigma;
int g_num_Kmat;


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
arma::fvec getMeanDiagofKmat_largeMem(){
	arma::fvec mean_diag_kins_vec(Kmat_vec.size() + 1);
	arma::fvec diagVecofKmat;
	for (int i = 0; i < Kmat_vec.size(); i++){
		diagVecofKmat = Kmat_vec[i].diag();
		mean_diag_kins_vec(i+1) = arma::mean(diagVecofKmat);
	}
	//arma::vec = g_spGRM.diag();

	//arma::fvec diagVecofKmat;
	diagVecofKmat = arma::diagvec(g_spGRM);
	//arma::vec diagVecofKmat_0;
        //diagVecofKmat_0	= arma::diagvec(g_spGRM);
	//diagVecofKmat = arma::conv_to< arma::fvec >::from(diagVecofKmat_0);
	mean_diag_kins_vec(0) = arma::mean(diagVecofKmat); 
	return(mean_diag_kins_vec);
}	



// [[Rcpp::export]]
int get_numofV(){
	int k = Kmat_vec.size();
	return(k);
}

// [[Rcpp::export]]
void set_longlVar_vec(arma::vec & longlVec){
	g_longl_vec = arma::conv_to< arma::fvec >::from(longlVec);
	//g_longl_vec = longlVec;
}

arma::umat g_covarianceidxMat;
arma::uvec g_covarianceidxMat_col1;
arma::uvec g_covarianceidxMat_col2;
arma::uvec g_covarianceidxMat_col3;
arma::uvec g_covarianceidxMat_notcol1;

//     [,1] [,2] [,3]
//[1,]    4    2    6
//[2,]    5    3    7
//2 GRM			2	
//3 In				5
//4 GRM*T + GRM*t(T)	3
//5 In*T + In*t(In)		6
//6 t(GRM *T) *T	4
//7 t(In * T) * T		7
/*> model5$theta
               dispersion       kins1.var.intercept       kins2.var.intercept
             3.136864e+00              4.582863e+00              0.000000e+00
kins1.cov.intercept.slope kins2.cov.intercept.slope           kins1.var.slope
            6.916919e-323              0.000000e+00              5.120541e-01
          kins2.var.slope
             0.000000e+00

//     [,1] [,2] [,3]
//[1,]  3    2    4 
//[2,]  6    5    7

*/

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
	  std::cout << "x here" << std::endl;
	  //spGRM_longl_0.clear();
	  spGRM_longl.clear();
	}
}

// [[Rcpp::export]]
void closeGenoFile_plink()
{
  //genoToTest_plainDosage.test_genoGZfile.close();
	for (int i = 0; i < geno.numofGenoArray; i++){
		(*geno.genoVecofPointers[i]).clear();	
    		delete geno.genoVecofPointers[i];
  	}

  	geno.genoVecofPointers.clear();

  	//geno.genoVec.clear();
  	geno.invstdvVec.clear();
  	geno.ptrsubSampleInGeno.clear();
  	geno.alleleFreqVec.clear();
  	geno.m_OneSNP_Geno.clear();
  	geno.m_OneSNP_StdGeno.clear();
  	geno.m_DiagStd.clear();
  	printf("closed the plinkFile!\n");
}


// [[Rcpp::export]]
int gettotalMarker(){
  	int numMarker = geno.getM();
  	return(numMarker); 
}

// [[Rcpp::export]]
arma::fvec getAlleleFreqVec(){
  	return(geno.alleleFreqVec);
}

// [[Rcpp::export]]
arma::ivec getMACVec(){
        return(geno.MACVec);
}


// [[Rcpp::export]]
arma::ivec getMACVec_forVarRatio(){
        return(geno.MACVec_forVarRatio);
}

// [[Rcpp::export]] 
arma::ivec getIndexVec_forVarRatio(){
	return(geno.markerIndexVec_forVarRatio);
}	

// [[Rcpp::export]]
bool getIsVarRatioGeno(){
	return(geno.isVarRatio);
}
// [[Rcpp::export]]
arma::ivec getSubMarkerIndex(){
	return(geno.subMarkerIndex);
}

// [[Rcpp::export]]
std::vector<bool> getQCdMarkerIndex(){
	return(geno.MarkerswithMAFge_minMAFtoConstructGRM_indVec);
}


// [[Rcpp::export]]
int getSubMarkerNum(){
        return(geno.subMarkerIndex.n_elem);
}


void initKinValueVecFinal(int ni){
	geno.kinValueVecFinal.resize(ni);
        std::fill(geno.kinValueVecFinal.begin(), geno.kinValueVecFinal.end(), 0);
};

// [[Rcpp::export]]
int getNnomissingOut(){
	return(geno.getNnomissing());
}

// [[Rcpp::export]]
int getMsub_MAFge_minMAFtoConstructGRM(){
	return(geno.getMsub_MAFge_minMAFtoConstructGRM_in());
}

// [[Rcpp::export]]
int getMsub_MAFge_minMAFtoConstructGRM_singleChr(){
        return(geno.getMsub_MAFge_minMAFtoConstructGRM_singleChr_in());
}


//  // [[Rcpp::export]]
//int getMmafge1perc(){
//	return(geno.getMmafge1perc());
//}
//arma::fmat Get_MultiMarkersBySample_StdGeno_Mat(arma::fvec& markerIndexVec){
//arma::fmat Get_MultiMarkersBySample_StdGeno_Mat(){

// [[Rcpp::export]]
void Get_MultiMarkersBySample_StdGeno_Mat(){
	//geno.subMarkerIndex
	//int m_M_Submarker = markerIndexVec.n_elem;
	int m_M_Submarker = getSubMarkerNum();
        //arma::fvec stdGenoMultiMarkers;
        int Nnomissing = geno.getNnomissing();
	  //int nSubMarker = markerIndexVec.n_elem;
          //int Ntotal = geno.getNnomissing();
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
		SNPIdx = (geno.subMarkerIndex)[k];
                indexOfVectorPointer = SNPIdx/(geno.numMarkersofEachArray);
                SNPIdxinVec = SNPIdx % (geno.numMarkersofEachArray);
                Start_idx = (geno.m_size_of_esi) * SNPIdxinVec;
                freq = (geno.alleleFreqVec)[SNPIdx];
                invStd = (geno.invstdvVec)[SNPIdx];
                if(k == 0){
                        std::cout << "freq: " << freq << " invStd: " << invStd << "  SNPIdx: " << SNPIdx << std::endl;
                }

                while(flag == 0){
//              std::cout << "createSparseKin1e" << std::endl;
                for(size_t i=Start_idx; i< Start_idx+(geno.m_size_of_esi); i++){
                        geno1 = (geno.genoVecofPointers)[indexOfVectorPointer]->at(i);
                        //std::cout << "createSparseKin1f" << std::endl;

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
			(geno.stdGenoMultiMarkersMat)(k, ind) = ((2-(a+b)) - 2*freq)* invStd;
//			std::cout << "k,ind " << k << " " << ind << std::endl;
//			std::cout << "(geno.stdGenoMultiMarkersMat)(k, ind) " << (geno.stdGenoMultiMarkersMat)(k, ind) << std::endl;

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

        std::cout << "stdGenoMultiMarkersMat.n_rows: " << geno.stdGenoMultiMarkersMat.n_rows << std::endl;
        std::cout << "stdGenoMultiMarkersMat.n_cols: " << geno.stdGenoMultiMarkersMat.n_cols << std::endl;
//	arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), m_M_Submarker, Nnomissing);

//	return(stdGenoMultiMarkersMat);
        //std::cout << "stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] " << stdGenoMultiMarkers[Nnomissing*m_M_Submarker-1] << std::endl;

}





// [[Rcpp::export]]
void Get_MultiMarkersBySample_StdGeno(arma::fvec& markerIndexVec, std::vector<float> &stdGenoMultiMarkers){

//	std::cout << "createSparseKin1c" << std::endl;
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
	int Nnomissing = geno.getNnomissing();
	

//	std::cout << "createSparseKin1d" << std::endl;
        for(size_t k=0; k< m_M_Submarker; k++){
                ind = 0;
                flag = 0;
                SNPIdx = markerIndexVec[k];
                indexOfVectorPointer = SNPIdx/(geno.numMarkersofEachArray);
                SNPIdxinVec = SNPIdx % (geno.numMarkersofEachArray);
                Start_idx = (geno.m_size_of_esi) * SNPIdxinVec;
		freq = (geno.alleleFreqVec)[SNPIdx];
                invStd = (geno.invstdvVec)[SNPIdx];
		//if(k == 0){
		//	std::cout << "freq: " << freq << " invStd: " << invStd << "  SNPIdx: " << SNPIdx << std::endl;
		//}

                while(flag == 0){
//		std::cout << "createSparseKin1e" << std::endl;
                for(size_t i=Start_idx; i< Start_idx+(geno.m_size_of_esi); i++){
                        geno1 = (geno.genoVecofPointers)[indexOfVectorPointer]->at(i);
			//std::cout << "createSparseKin1f" << std::endl;

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        stdGenoMultiMarkers[ind*m_M_Submarker+k] = ((2-(a+b)) - 2*freq)* invStd;;
//			stdGenoMultiMarkers[ind*m_M_Submarker+k] = 2-(a+b);
//			if(k == 0){
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
  		
  		m_M = geno.getM();
  		m_N = geno.getNnomissing();
  		m_bout.zeros(m_N);
		Msub_mafge1perc=0;
		//geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
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
			//if(geno.alleleFreqVec[i] >= minMAFtoConstructGRM && geno.alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
				geno.Get_OneSNP_StdGeno(i, &vec);
				float val1 = dot(vec,  m_bVec);
				m_bout += val1 * (vec) ;
				Msub_mafge1perc += 1;
			//}
			//std::cout << "i: " << i << std::endl;
			//for(unsigned int j = 0; j < 10; j++){
			//	std::cout << "m_bVec[j] " << m_bVec[j] << std::endl;
			//	std::cout << "vec[j] " << vec[j] << std::endl;
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

                m_Msub = geno.getMsub(); //LOCO
		startIndex = geno.getStartIndex();
		endIndex = geno.getEndIndex();
                m_M = geno.getM(); //LOCO
                m_N = geno.getNnomissing();
                m_bout.zeros(m_N);
		m_Msub_mafge1perc=0;
		//geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        }
        CorssProd_LOCO(const CorssProd_LOCO& CorssProd_LOCO, Split)
                : m_bVec(CorssProd_LOCO.m_bVec)
        {

                m_N = CorssProd_LOCO.m_N;
                m_M = CorssProd_LOCO.m_M;
		m_Msub = CorssProd_LOCO.m_Msub;
		startIndex = geno.getStartIndex();
                endIndex = geno.getEndIndex();
                m_bout.zeros(m_N);
		m_Msub_mafge1perc=0;
		//geno.getnumberofMarkers_byChr(uint chr);
        }
	
	   // process just the elements of the range I've been asked to
        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
		float val1;
                for(unsigned int i = begin; i < end; i++){
                        geno.Get_OneSNP_StdGeno(i, &vec);
		//	if(i >= startIndex && i <= endIndex){
		//		val1 = 0;
					//if(endIndex == 4){
					//		cout << "i: " << i << endl;
					//}
		//	}else{
                        val1 = dot(vec,  m_bVec);
	       		m_Msub_mafge1perc += 1;
		//	}
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
  	//int M = geno.getM();
 	//int Msub_mafge1perc = geno.getMmafge1perc();
	int Msub_mafge1perc = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
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
        //int M = geno.getM();
	//
        int Msub_mafge1perc = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM(); 
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
        //int Msub = geno.getMsub();
	//int M = geno.getM();
        int numberMarker_full = 0;
        arma::fvec outvec = parallelCrossProd_full(bVec, numberMarker_full);

        //CorssProd_LOCO CorssProd_LOCO(bVec);
	CorssProd CorssProd(bVec);
  // call paralleReduce to start the work
	int startIndex = geno.getStartIndex();
        int endIndex = geno.getEndIndex();

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
        //	cout << (CorssProd_LOCO.m_bout/Msub)[i] << ' ';
        //}
	//cout << endl;
        //return CorssProd_LOCO.m_bout/Msub;
        //return CorssProd_LOCO.m_bout/(CorssProd_LOCO.m_Msub_mafge1perc);
	return outvec/markerNum;
}


arma::umat locationMat;
arma::vec valueVec;
int dimNum = 0;

// [[Rcpp::export]]
void setupSparseGRM(int r, arma::umat & locationMatinR, arma::vec & valueVecinR) {
    // sparse x sparse -> sparse
    //arma::sp_mat result(a);
    //int r = a.n_rows;
    locationMat.zeros(2,r);
    valueVec.zeros(r);

    locationMat = locationMatinR;
    valueVec = valueVecinR;
    dimNum = r;

    std::cout << locationMat.n_rows << " locationMat.n_rows " << std::endl;
    std::cout << locationMat.n_cols << " locationMat.n_cols " << std::endl;
    std::cout << valueVec.n_elem << " valueVec.n_elem " << std::endl;
    
    //for(size_t i=0; i< 10; i++){
    //    std::cout << valueVec(i) << std::endl;
    //    std::cout << locationMat(0,i) << std::endl;
    //    std::cout << locationMat(1,i) << std::endl;
    //}

    //arma::vec y = arma::linspace<arma::vec>(0, 5, r);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    //arma::vec x = arma::spsolve( result, y );

    //return x;
}

// [[Rcpp::export]]
void setupSparseGRM_new(arma::sp_mat & t_spGRM){
	arma::sp_mat t_spGRM_1 = t_spGRM;
	arma::sp_fmat g_spGRM_f = arma::conv_to<arma::sp_fmat>::from(t_spGRM_1);
	g_spGRM = g_spGRM_f;
}


bool isUsePrecondM = false;
bool isUseSparseSigmaforInitTau = false;
bool isUseSparseSigmaforModelFitting = false;



// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin(arma::fcolvec& bVec){
    arma::fvec crossProdVec;
    if(isUseSparseSigmaforInitTau | isUseSparseSigmaforModelFitting){
        //arma::dcolvec bVec_new = arma::conv_to<arma::dcolvec>::from(bVec);
        //cout << "use sparse kinship to estimate initial tau and for getCrossprodMatAndKin" <<  endl;
	//arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
	//arma::vec x = result * bVec_new;
	//arma::vec x = g_spGRM * bVec_new;
        // double wall3in = get_wall_time();
        // double cpu3in  = get_cpu_time();
        // cout << "Wall Time in gen_spsolve_v4 = " << wall3in - wall2in << endl;
        // cout << "CPU Time  in gen_spsolve_v4 = " << cpu3in - cpu2in  << endl;
        //crossProdVec = arma::conv_to<arma::fvec>::from(x);
	crossProdVec = g_spGRM * bVec;

    }else{ 
  	crossProdVec = parallelCrossProd(bVec) ;
    }  
    return(crossProdVec);
}






// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin_LOCO(arma::fcolvec& bVec){

        arma::fvec crossProdVec = parallelCrossProd_LOCO(bVec) ;
        //arma::fvec crossProdVec_2 = parallelCrossProd_LOCO_2(bVec) ;

	//for(int k=0; k < 10; k++) {
        //	std::cout << "new crossProdVec " << k << " " << crossProdVec[k] << std::endl;
        //	std::cout << "old crossProdVec " << k << " " << crossProdVec_2[k] << std::endl;

	//}	


        return(crossProdVec);
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
      float kinValueTemp = arma::dot((geno.stdGenoMultiMarkersMat).col(i), (geno.stdGenoMultiMarkersMat).col(j));
      kinValueTemp = kinValueTemp/m_M_Submarker;
      if(kinValueTemp >=  geno.relatednessCutoff) {
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


//arma::fmat findIndiceRelatedSample(){
//arma::fmat findIndiceRelatedSample(){

// [[Rcpp::export]]
void findIndiceRelatedSample(){

  int Ntotal = geno.getNnomissing(); 
//  tbb::concurrent_vector< std::pair<float, float> > output;

//  indicesRelatedSamples indicesRelatedSamples(Ntotal,output);
  indicesRelatedSamples indicesRelatedSamples(Ntotal,geno.indiceVec);

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
    (geno.indiceVec).push_back( std::pair<int, int>(i, i) );
  }
*/

//  return(xout);
}



struct sparseGRMUsingOneMarker : public Worker {
   // input matrix to read from
  // arma::imat & iMat;
   // output matrix to write to
   arma::fvec & GRMvec;

   //int M = geno.getM();
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
	   int iint = (geno.indiceVec)[i].first;	
	   int jint = (geno.indiceVec)[i].second;	
/*
            float ival = geno.m_OneSNP_StdGeno(iint);
            float jval = geno.m_OneSNP_StdGeno(jint);
            // write to output matrix
            //rmat(i,j) = sqrt(.5 * (d1 + d2));
            GRMvec(i) = ival*jval/M;
*/
	//use Look-Up table for calucate GRMvec(i)
	    int ival = geno.m_OneSNP_Geno(iint);	
	    int jval = geno.m_OneSNP_Geno(jint);
	    GRMvec(i) = geno.sKinLookUpArr[ival][jval]; 

      }
   }
};

//void parallelcalsparseGRM(arma::imat & iMat, arma::fvec &GRMvec) {

// [[Rcpp::export]]
void parallelcalsparseGRM(arma::fvec &GRMvec) {

//  int n1 = geno.indiceVec.size();
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
  parallelFor(0, (geno.indiceVec).size(), sparseGRMUsingOneMarker);

  // return the output matrix
  // return GRMvec;
}


struct sumTwoVec : public Worker
{   
   // source vectors
   arma::fvec &x;
   
   arma::fvec &sumVec;
  
   //int M = geno.getM(); 
   // constructors
   sumTwoVec(arma::fvec &x,arma::fvec &sumVec) 
      : x(x), sumVec(sumVec) {}
   
     // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
            // rows we will operate on
            sumVec(i) = x(i)+(geno.kinValueVecFinal)[i];
	    (geno.kinValueVecFinal)[i] = sumVec(i);	
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
void setgeno(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool isDiagofKinSetAsOne)
{
	int start_s=clock();
        geno.setGenoObj(bedfile, bimfile, famfile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne);
	//geno.printAlleleFreqVec();
	//geno.printGenoVec();
	int stop_s=clock();
	cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
}




// [[Rcpp::export]]
arma::ivec Get_OneSNP_Geno(int SNPIdx)
{

	arma::ivec temp = * geno.Get_OneSNP_Geno(SNPIdx);
	return(temp);

}


// [[Rcpp::export]]
arma::ivec Get_OneSNP_Geno_forVarRatio(int SNPIdx)
{
       
        arma::ivec temp = * geno.Get_OneSNP_Geno_forVarRatio(SNPIdx);
        return(temp);

}



// [[Rcpp::export]]
arma::fvec Get_OneSNP_StdGeno(int SNPIdx)
{

	arma::fvec temp; 
	geno.Get_OneSNP_StdGeno(SNPIdx, & temp);
//	for(int j = 0; j < 100; j++){
//                std::cout << "temp(j): " << j << " " << temp(j) << std::endl;

 //       }


	return(temp);

}
  
    
  

//Sigma = tau[1] * diag(1/W) + tau[2] * kins 
// [[Rcpp::export]]
arma::fvec getDiagOfSigma_largeMem__multiV(arma::fvec& wVec, arma::fvec& tauVec, bool LOCO){
  
	int Nnomissing = geno.getNnomissing();
	arma::fvec diagVec(Nnomissing);
	float diagElement;
	float floatBuffer;
	//int M = geno.getM();
	int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
	//cout << "MminMAF=" << MminMAF << endl;
	//cout << "M=" << M << endl;
  	//float minvElement;
  	if(!(geno.setKinDiagtoOne)){ 
	  //diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /M + tauVec(0)/wVec;
          if(!LOCO){ 
	    diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /MminMAF + tauVec(0)/wVec;
	  }else{
	    diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno_LOCO());
            int Msub_MAFge_minMAFtoConstructGRM_in_b = geno.getMsub_MAFge_minMAFtoConstructGRM_in();
            int Msub_MAFge_minMAFtoConstructGRM_singleVar_b = geno.getMsub_MAFge_minMAFtoConstructGRM_singleChr_in();	    
	    diagVec = diagVec/(Msub_MAFge_minMAFtoConstructGRM_in_b - Msub_MAFge_minMAFtoConstructGRM_singleVar_b) + tauVec(0)/wVec;
          }

	}else{
	  diagVec = tauVec(1) + tauVec(0)/wVec;
	}

        if(Kmat_vec.size() > 0){
          for(unsigned int i = 0; i < Kmat_vec.size(); i++){
            diagVec = diagVec + (Kmat_vec[i]).diag() * tauVec(i+2);
          }
        } 


	//std::cout << "M " << M << std::endl;
	//std::cout << "tauVec(0) " << tauVec(0) << std::endl;
	//std::cout << "tauVec(1) " << tauVec(1) << std::endl;
        //for(unsigned int i=0; i< 10; i++){
	//	 std::cout << "diagVec(i) " << diagVec(i) << std::endl;
	//}

	//make diag of kin to be 1 to compare results of emmax and gmmat
	//diagVec = tauVec(1) + tauVec(0)/wVec;
	for(unsigned int i=0; i< Nnomissing; i++){
//	if(i < 100){
//		std::cout << i << "th element of diag of sigma and wVec " << diagVec(i) << " " << wVec(i) << std::endl;
//	}
  		if(diagVec(i) < 1e-4){
  			diagVec(i) = 1e-4 ;
  		}
  	}
 
    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
  	return(diagVec);
}



// [[Rcpp::export]]
arma::fvec getDiagOfSigma_multiV(arma::fvec& wVec, arma::fvec& tauVec, bool LOCO){

        int Nnomissing = geno.getNnomissing();
        arma::fvec diagVec(Nnomissing);
        float diagElement;
        float floatBuffer;
        //int M = geno.getM();
        int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
	arma::sp_vec diagVecG0;
        arma::sp_fvec diagVecV0;
	arma::fvec diagVecG, diagVecV, diagVecG_I, diagVecG_T, diagVecG_IT,diagVecV_I, diagVecV_T, diagVecV_IT;
        //cout << "MminMAF=" << MminMAF << endl;
        //cout << "M=" << M << endl;
        //float minvElement;
	//
	unsigned int tauind = 0;
	if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
          if(!(geno.setKinDiagtoOne)){
          //diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /M + tauVec(0)/wVec;
            if(!LOCO){
              diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno()) /MminMAF + tauVec(0)/wVec;
            }else{
              diagVec = tauVec(1)* (*geno.Get_Diagof_StdGeno_LOCO());
              int Msub_MAFge_minMAFtoConstructGRM_in_b = geno.getMsub_MAFge_minMAFtoConstructGRM_in();
              int Msub_MAFge_minMAFtoConstructGRM_singleVar_b = geno.getMsub_MAFge_minMAFtoConstructGRM_singleChr_in();
              diagVec = diagVec/(Msub_MAFge_minMAFtoConstructGRM_in_b - Msub_MAFge_minMAFtoConstructGRM_singleVar_b) + tauVec(0)/wVec;
            }

          }else{
            diagVec = tauVec(1) + tauVec(0)/wVec;
          }
	  
	  if(Kmat_vec.size() > 0){
            for(unsigned int i = 0; i < Kmat_vec.size(); i++){
              diagVec = diagVec + (Kmat_vec[i]).diag() * tauVec(i+2);
            }
          }
	
	}else{ //if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
		diagVec = tauVec(0)/wVec;
		tauind = tauind + 1;
		//diagVecG0 = g_spGRM.diag();
		//arma::vec diagVecGtemp(diagVecG0);
		//diagVecG = arma::conv_to< arma::fvec >::from(diagVecGtemp);
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


        //std::cout << "M " << M << std::endl;
        //std::cout << "tauVec(0) " << tauVec(0) << std::endl;
        //std::cout << "tauVec(1) " << tauVec(1) << std::endl;
        //for(unsigned int i=0; i< 10; i++){
        //       std::cout << "diagVec(i) " << diagVec(i) << std::endl;
        //}

        //make diag of kin to be 1 to compare results of emmax and gmmat
        //diagVec = tauVec(1) + tauVec(0)/wVec;
        for(unsigned int i=0; i< Nnomissing; i++){
//      if(i < 100){
//              std::cout << i << "th element of diag of sigma and wVec " << diagVec(i) << " " << wVec(i) << std::endl;
//      }
                if(diagVec(i) < 1e-4){
                        diagVec(i) = 1e-4 ;
                }
        }

    //cout << *geno.Get_Diagof_StdGeno() << endl ;
    //cout << diagVec << endl ;
        return(diagVec);
}


// [[Rcpp::export]]
arma::fcolvec getCrossprod_largeMem_multiV(arma::fcolvec& bVec, arma::fvec& wVec, arma::fvec& tauVec, bool LOCO){

        arma::fcolvec crossProdVec;
        // Added by SLEE, 04/16/2017
        if(tauVec(1) == 0){
                crossProdVec = tauVec(0)*(bVec % (1/wVec));
                return(crossProdVec);
        }
        //
        arma::fvec crossProd1;
      
	if(!LOCO){
		crossProd1 = getCrossprodMatAndKin(bVec);
	}else{
		crossProd1 = getCrossprodMatAndKin_LOCO(bVec);
	}	
	//for(int j = 0; j < 100; j++){
        //        std::cout << "bVec(j): " << bVec(j) << std::endl;
        //        std::cout << "crossProd1(j): " << crossProd1(j) << std::endl;

        //}
	
        crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;

        //arma::dcolvec bVec_new = arma::conv_to<arma::dcolvec>::from(bVec);
        // double wall3in = get_wall_time();
        // double cpu3in  = get_cpu_time();
        // cout << "Wall Time in gen_spsolve_v4 = " << wall3in - wall2in << endl;
        // cout << "CPU Time  in gen_spsolve_v4 = " << cpu3in - cpu2in  << endl;
        //crossProdVec = arma::conv_to<arma::fvec>::from(x);
        //arma::dcolvec crossProdVec_temp;
	if(Kmat_vec.size() > 0){
		for(unsigned int i = 0; i < Kmat_vec.size(); i++){
			crossProdVec  = crossProdVec + tauVec(i+2)*(Kmat_vec[i] * bVec); 
			//crossProdVec = crossProdVec + arma::conv_to<arma::fvec>::from(x);
			//arma::conv_to<arma::dcolvec>::from(crossProdVec_temp);
		}	
	}	

	//for(int j = 0; j < 10; j++){
        //        std::cout << "crossProdVec(j): " << j << " " << crossProdVec(j) << std::endl;
        //}
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
	if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
		
        	if(!LOCO){
                	crossProd1 = getCrossprodMatAndKin(bVec);
        	}else{
                	crossProd1 = getCrossprodMatAndKin_LOCO(bVec);
        	}
        	crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;
		tau_ind = tau_ind + 2;


	}else{
		Ibvec = g_I_longl_mat.t() * bVec;
		if(!LOCO){
                        GRM_I_bvec = getCrossprodMatAndKin(Ibvec);
                }else{
                        GRM_I_bvec = getCrossprodMatAndKin_LOCO(Ibvec);
                }
		crossProd1 = g_I_longl_mat * GRM_I_bvec;
		crossProdVec = tauVec(0)*(bVec % (1/wVec)) + tauVec(1)*crossProd1;		
		tau_ind = tau_ind + 2;	


		if(g_T_longl_mat.n_rows > 0){
			Tbvec = g_T_longl_mat.t() * bVec;
			if(!LOCO){
                                GRM_T_bvec = getCrossprodMatAndKin(Tbvec);
                        }else{
                                GRM_T_bvec = getCrossprodMatAndKin_LOCO(Tbvec);
                        }
			crossProdGRM_TGIb = g_T_longl_mat * GRM_I_bvec;
                        crossProdGRM_IGTb = g_I_longl_mat * GRM_T_bvec;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (crossProdGRM_TGIb + crossProdGRM_IGTb);
                        tau_ind = tau_ind + 1;
                        crossProdVec = crossProdVec + tauVec(tau_ind) * (g_T_longl_mat * GRM_T_bvec);
                        tau_ind = tau_ind + 1;

		}	

	}	

        //for(int j = 0; j < 100; j++){
        //        std::cout << "bVec(j): " << bVec(j) << std::endl;
        //        std::cout << "crossProd1(j): " << crossProd1(j) << std::endl;

        //}

	if(g_T_longl_mat.n_rows == 0 && g_I_longl_mat.n_rows == 0){

        //arma::dcolvec bVec_new = arma::conv_to<arma::dcolvec>::from(bVec);
        // double wall3in = get_wall_time();
        // double cpu3in  = get_cpu_time();
        // cout << "Wall Time in gen_spsolve_v4 = " << wall3in - wall2in << endl;
        // cout << "CPU Time  in gen_spsolve_v4 = " << cpu3in - cpu2in  << endl;
        //crossProdVec = arma::conv_to<arma::fvec>::from(x);
        //arma::dcolvec crossProdVec_temp;
	//	
	
        	if(Kmat_vec.size() > 0){
                	for(unsigned int i = 0; i < Kmat_vec.size(); i++){
                        	crossProdVec  = crossProdVec + tauVec(tau_ind)*(Kmat_vec[i] * bVec);
				tau_ind = tau_ind + 1;
                        //crossProdVec = crossProdVec + arma::conv_to<arma::fvec>::from(x);
                        //arma::conv_to<arma::dcolvec>::from(crossProdVec_temp);
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
        //for(int j = 0; j < 10; j++){
        //        std::cout << "crossProdVec(j): " << j << " " << crossProdVec(j) << std::endl;
        //}
        return(crossProdVec);
}




double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}



// [[Rcpp::export]]
arma::sp_mat gen_sp_GRM() {
    // sparse x sparse -> sparse
    arma::sp_mat result(locationMat, valueVec, dimNum, dimNum);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    return result;
}




// [[Rcpp::export]]
arma::sp_fmat gen_sp_Sigma_largeMem_multiV(arma::fvec& wVec,  arma::fvec& tauVec){
   arma::fvec dtVec = (1/wVec) * (tauVec(0));
   std::cout << "tauVec(0) " << tauVec(0) << std::endl;

//   dtVec.print();
/*
   arma::vec valueVecNew = valueVec * tauVec(1);

   int nnonzero = valueVec.n_elem;
   for(size_t i=0; i< nnonzero; i++){
     if(locationMat(0,i) == locationMat(1,i)){
//       std::cout << "i: " << i << " " << valueVecNew(i) << std::endl;
       valueVecNew(i) = valueVecNew(i) + dtVec(locationMat(0,i));
//       std::cout << "i: " << i << " " << valueVecNew(i) << std::endl;
	if(valueVecNew(i) < 1e-4){
  			valueVecNew(i) = 1e-4 ;
  		}


     }
   }

    // sparse x sparse -> sparse
    arma::sp_mat result(locationMat, valueVecNew, dimNum, dimNum);

*/
   arma::sp_fmat result = g_spGRM * tauVec(1);
   //std::cout << "tauVec(1) " << tauVec(1) << std::endl;
   //std::cout << "g_spGRM(1775,1775) " << g_spGRM(1775,1775) << std::endl;
   result.diag() = result.diag() + dtVec;   
    if(Kmat_vec.size() > 0){	
      for(unsigned int i = 0; i < Kmat_vec.size(); i++){
        result = result + Kmat_vec[i] * tauVec(i+2);
	//std::cout << "tauVec(i+2) " << tauVec(i+2) << std::endl;
	//std::cout << "i " << i << std::endl;
	//std::cout << "Kmat_vec[i] " << Kmat_vec[i](1775,1775) << std::endl;
      }	      
    }
//    std::cout << "result.n_rows " << result.n_rows << std::endl;
//    std::cout << "result.n_cols " << result.n_cols << std::endl;
    //result.print();
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    return result;
}



// [[Rcpp::export]]
void gen_sp_Sigma_multiV(arma::fvec& wVec,  arma::fvec& tauVec){
   arma::fvec dtVec = (1/wVec) * (tauVec(0));
   std::cout << "tauVec(0) " << tauVec(0) << std::endl;
   arma::sp_fmat GRM_Imat, GRM_Tmat;

   arma::fvec crossProd1, GRM_I_bvec, Ibvec, Tbvec, GRM_T_bvec, crossProdGRM_TGIb, crossProdGRM_IGTb, V_I_bvec, V_T_bvec, crossProdV_TGIb, crossProdV_IGTb, crossProdGRM_TIb, crossProdGRM_ITb;
   unsigned int tau_ind = 0;
   //arma::sp_fmat g_spGRM_f = arma::conv_to< arma::sp_fmat >::from(g_spGRM); 

   if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
       g_spSigma = g_spGRM * tauVec(1);
       g_spSigma.diag() = g_spSigma.diag() + dtVec;
       tau_ind = tau_ind + 2;


   }else{
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


/*
// [[Rcpp::export]]
arma::fvec gen_spsolve_v4_multiV(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec & yvec){

//    double wall0in = get_wall_time();
// double cpu0in  = get_cpu_time();

    arma::vec yvec2 = arma::conv_to<arma::vec>::from(yvec);
//double wall1in = get_wall_time();
// double cpu1in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall1in - wall0in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu1in - cpu0in  << endl;

    arma::sp_mat result = gen_sp_Sigma_multiV(wVec, tauVec);

    

//double wall2in = get_wall_time();
// double cpu2in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall2in - wall1in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu2in - cpu1in  << endl;

//    std::cout << "yvec.n_elem: " << yvec.n_elem << std::endl;
//    std::cout << "yvec2.n_elem: " << yvec2.n_elem << std::endl;
//    std::cout << "result.n_rows: " << result.n_rows << std::endl;
//    std::cout << "result.n_cols: " << result.n_cols << std::endl;
    arma::vec x = arma::spsolve(result, yvec2);

//double wall3in = get_wall_time();
// double cpu3in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall3in - wall2in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu3in - cpu2in  << endl;


    arma::fvec z = arma::conv_to<arma::fvec>::from(x);

//double wall4in = get_wall_time();
// double cpu4in  = get_cpu_time();
// cout << "Wall Time in gen_spsolve_v4 = " << wall4in - wall3in << endl;
// cout << "CPU Time  in gen_spsolve_v4 = " << cpu4in - cpu3in  << endl;


    return z;
}

*/

//bool isUsePrecondM = false;
//bool isUseSparseSigmaforInitTau = false;


// [[Rcpp::export]]
void setisUsePrecondM(bool isUseSparseSigmaforPCG){
	isUsePrecondM = isUseSparseSigmaforPCG;
}

// [[Rcpp::export]]
void setisUseSparseSigmaforInitTau(bool isUseSparseSigmaforInitTau0){
	isUseSparseSigmaforInitTau = isUseSparseSigmaforInitTau0;
}



// [[Rcpp::export]]
void setisUseSparseSigmaforNullModelFitting(bool isUseSparseSigmaforModelFitting0){
        isUseSparseSigmaforModelFitting = isUseSparseSigmaforModelFitting0;
}


//Modified on 11-28-2018 to allow for a preconditioner for CG (the sparse Sigma)                                                                                                                                     //Sigma = tau[1] * diag(1/W) + tau[2] * kins + tau[3] * V3 + .....
//This function needs the function getDiagOfSigma and function getCrossprod

// [[Rcpp::export]]
arma::fvec getPCG1ofSigmaAndVector_multiV(arma::fvec& wVec,  arma::fvec& tauVec, arma::fvec& bVec, int maxiterPCG, float tolPCG, bool LOCO){       
    // Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    int Nnomissing = geno.getNnomissing();
    arma::fvec xVec(Nnomissing);
    xVec.zeros();
    
    //if(isUseSparseSigmaforInitTau){
    //	cout << "use sparse kinship to estimate initial tau " <<  endl;
    //	xVec = gen_spsolve_v4_multiV(wVec, tauVec, bVec);
    //}else if(isUseSparseSigmaforModelFitting){
    //}else if(!isUseSparseSigmaforModelFitting){
    //	cout << "use sparse kinship to fit the model " << endl;
    //    xVec = gen_spsolve_v4_multiV(wVec, tauVec, bVec);
    if(g_isStoreSigma){
      //arma::vec bVec0 = arma::conv_to<arma::vec>::from(bVec);
      //arma::vec x = arma::spsolve(g_spSigma, bVec0);
      //xVec = arma::conv_to<arma::fvec>::from(x);
      //
      //
	    std::cout << " arma::spsolve(g_spSigma, bVec) 0" << std::endl;
      xVec = arma::spsolve(g_spSigma, bVec);
	    std::cout << " arma::spsolve(g_spSigma, bVec) 1" << std::endl;
    }else{
        arma::fvec rVec = bVec;
        arma::fvec r1Vec;
        //cout << "HELLOb: "  << endl;
        //int Nnomissing = geno.getNnomissing();
        //cout << "HELL1: "  << endl;
        arma::fvec crossProdVec(Nnomissing);
        //arma::SpMat<float> precondM = sparseGRMinC;
        arma::fvec zVec(Nnomissing);
        arma::fvec minvVec(Nnomissing);
 //     double wall1 = get_wall_time();
 //     double cpu1  = get_cpu_time();

//      cout << "Wall Time 1= " << wall1 - wall0 << endl;
//      cout << "CPU Time 1 = " << cpu1  - cpu0  << endl;
        if (!isUsePrecondM){
	    double wall1_gDiag = get_wall_time();
	    double cpu1_gDiag  = get_cpu_time();
/*		for(int i = 0; i < 10; i++){
                	cout << "wVec[i]: " << i << " " << wVec[i] << endl;
                	cout << "bVec[i]: " << i << " " << bVec[i] << endl;
        	}
*/
	    minvVec = 1/getDiagOfSigma_multiV(wVec, tauVec, LOCO);
	    //	for(int i = 0; i < 10; i++){
          //      	cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
        //	}

	    double wall1_gDiag_2 = get_wall_time();
	    double cpu1_gDiag_2  = get_cpu_time();
//		 cout << "Wall Time getDiagOfSigma = " << wall1_gDiag_2 - wall1_gDiag << endl;
// cout << "CPU Time getDiagOfSigma = " << cpu1_gDiag_2 - cpu1_gDiag  << endl;
            zVec = minvVec % rVec;

        }else{
//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
	     //zVec = gen_spsolve_v4_multiV(wVec, tauVec, rVec);
	     //arma::vec rVec0 = arma::conv_to<arma::vec>::from(rVec);	
	      //arma::vec zVec0 = arma::spsolve(g_spSigma, rVec0);
      	      //zVec = arma::conv_to<arma::fvec>::from(zVec0);
	      zVec = arma::spsolve(g_spSigma, rVec);
                //sparseGRMinC = (sparseGRMinC) * (tauVec(1));
                //arma::fvec dtVec = (1/wVec) * (tauVec(0));
                //(sparseGRMinC).diag() = (sparseGRMinC).diag() + dtVec;
                //zVec = gen_spsolve_v4(sparseGRMinC, rVec) ;
        }
// double wall1 = get_wall_time();
// double cpu1  = get_cpu_time();
// cout << "Wall Time 1 = " << wall1 - wall0 << endl;
// cout << "CPU Time 1 = " << cpu1  - cpu0  << endl;


//      cout << "HELL3: "  << endl;
      //for(int i = 0; i < 10; i++){
       //         cout << "zVec[i]: " << zVec[i] << endl;
        //}
        float sumr2 = sum(rVec % rVec);

/*
        if(bVec[0] == 1 && bVec[99] == 1){
        for(int i = 0; i < 100; i++){
                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                cout << "minvVec[i]: " << i << " " << minvVec[i] << endl;
                cout << "wVec[i]: " << i << " " << wVec[i] << endl;
        }
        }
*/
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;
        /*
        if(bVec[0] == 1 && bVec[2] == 1){
        for(int i = 0; i < 10; i++){
                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
        }
        }
*/
        //arma::fvec xVec(Nnomissing);
        //xVec.zeros();

        int iter = 0;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
// double wall2 = get_wall_time();
//        double cpu2  = get_cpu_time();

//    cout << "Wall Time 2= " << wall2 - wall1 << endl;
//    cout << "CPU Time 2 = " << cpu2  - cpu1  << endl;



                iter = iter + 1;
                arma::fcolvec ApVec = getCrossprod_multiV(pVec, wVec, tauVec, LOCO);
//		for(int i = 0; i < 10; i++){
//                	cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
//        	}
                arma::fvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);
// double wall3 = get_wall_time();
//        double cpu3  = get_cpu_time();

//    cout << "Wall Time 3= " << wall3 - wall2 << endl;
//    cout << "CPU Time 3 = " << cpu3  - cpu2  << endl;
                float a = preA(0);

/*           if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "bVec[0] == 1 && bVec[2] == 1: " << endl;
                        for(int i = 0; i < 10; i++){

                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "pVec[i]: " << i << " " << pVec[i] << endl;
                                cout << "zVec[i]: " << i << " " << zVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                        }
                    }
*/

                xVec = xVec + a * pVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        for(int i = 0; i < 10; i++){
                                cout << "xVec[i]: " << i << " " << xVec[i] << endl;
                        }
                }

*/


                r1Vec = rVec - a * ApVec;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        cout << "a: " << a  << endl;
                        for(int i = 0; i < 10; i++){
                                cout << "ApVec[i]: " << i << " " << ApVec[i] << endl;
                                cout << "rVec[i]: " << i << " " << rVec[i] << endl;
                                cout << "r1Vec[i]: " << i << " " << r1Vec[i] << endl;
                        }
                }
*/
//                z1Vec = minvVec % r1Vec;
// double wall3a = get_wall_time();
//       double cpu3a  = get_cpu_time();

        if (!isUsePrecondM){
                z1Vec = minvVec % r1Vec;
        }else{
	      //arma::vec r1Vec0 = arma::conv_to<arma::vec>::from(r1Vec);
              //arma::vec z1Vec0 = arma::spsolve(g_spSigma, r1Vec0);
              //z1Vec = arma::conv_to<arma::fvec>::from(z1Vec0);
		//z1Vec = gen_spsolve_v4_multiV(wVec, tauVec, r1Vec);
                //z1Vec = arma::spsolve(sparseGRMinC, r1Vec) ;
		z1Vec = arma::spsolve(g_spSigma, r1Vec);
        }

//       double wall3b = get_wall_time();
//       double cpu3b  = get_cpu_time();
// cout << "Wall Time 3b = " << wall3b - wall3a << endl;
// cout << "CPU Time 3b = " << cpu3b  - cpu3a  << endl;


                arma::fvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
                float bet = Prebet(0);
                pVec = z1Vec+ bet*pVec;
                zVec = z1Vec;
                rVec = r1Vec;

// double wall4 = get_wall_time();
//        double cpu4  = get_cpu_time();

//    cout << "Wall Time 4= " << wall4 - wall3 << endl;
//    cout << "CPU Time 4 = " << cpu4  - cpu3  << endl;

                sumr2 = sum(rVec % rVec);
                //        std::cout << "sumr2: " << sumr2 << std::endl;
                //        std::cout << "tolPCG: " << tolPCG << std::endl;
/*
                if(bVec[0] == 1 && bVec[2] == 1){
                        std::cout << "sumr2: " << sumr2 << std::endl;
                        std::cout << "tolPCG: " << tolPCG << std::endl;
                }
*/
        }

        if (iter >= maxiterPCG){
                cout << "pcg did not converge. You may increase maxiter number." << endl;

        }
        cout << "iter from getPCG1ofSigmaAndVector " << iter << endl;
} //else if(isUseSparseKinforInitTau){
//  double wall5 = get_wall_time();
//    double cpu5  = get_cpu_time();
//    cout << "Wall Time getPCG1ofSigmaAndVector = " << wall5 - wall0 << endl;
//    cout << "CPU Time  getPCG1ofSigmaAndVector = " << cpu5  - cpu0  << endl;


//      std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
//        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
        return(xVec);
}


//http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
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

/*
// [[Rcpp::export]] 
void setChromosomeIndicesforLOCO(vector<int> chromosomeStartIndexVec, vector<int> chromosomeEndIndexVec, vector<int> chromosomeVecVec){
  LOCO = true;
  chromosomeStartIndex = chromosomeStartIndexVec;
  chromosomeEndIndexV = chromosomeEndIndexVec;
  chromosomeVec = chromosomeVecVec;
}
*/

// [[Rcpp::export]]
void setStartEndIndex(int startIndex, int endIndex, int chromIndex){
  geno.startIndex = startIndex;
  geno.endIndex = endIndex;
  geno.Msub = 0;
  geno.chromIndex = chromIndex;

  for(size_t i=0; i< geno.M; i++){
	if(i < startIndex || i > endIndex){
  		if(geno.alleleFreqVec[i] >= minMAFtoConstructGRM && geno.alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
      
			geno.Msub = geno.Msub + 1;
  		}
	}
  }
  //geno.Msub = geno.M - (endIndex - startIndex + 1);
}



// [[Rcpp::export]]
void setStartEndIndexVec( arma::ivec & startIndex_vec,  arma::ivec & endIndex_vec){	
  geno.startIndexVec = startIndex_vec;
  geno.endIndexVec = endIndex_vec;
  //geno.Msub = geno.M - (endIndex - startIndex + 1);
}

// // [[Rcpp::export]]
//void setStartEndIndexVec_forvr( arma::ivec & startIndex_vec,  arma::ivec & endIndex_vec){
//  geno.startIndexVec_forvr = startIndex_vec;
//  geno.endIndexVec_forvr = endIndex_vec;
  //geno.Msub = geno.M - (endIndex - startIndex + 1);
//}



//This function calculates the coefficients of variation for mean of a vector
// [[Rcpp::export]]
float calCV(arma::fvec& xVec){
  int veclen = xVec.n_elem;
  float vecMean = arma::mean(xVec);
  float vecSd = arma::stddev(xVec);
  float vecCV = (vecSd/vecMean)/veclen;
  return(vecCV);
}






/*add for SPA by Wei 04222017*/
// [[Rcpp::export]]
arma::fmat getSigma_X_multiV(arma::fvec& wVec, arma::fvec& tauVec,arma::fmat& Xmat, int maxiterPCG, float tolPCG, bool LOCO){


  	int Nnomissing = Xmat.n_rows;
  	int colNumX = Xmat.n_cols;

  	//cout << colNumX << endl;
  	//cout << size(wVec) << endl;
  	//cout << size(tauVec) << endl;


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

                //m_Msub = geno.getMsub();
                subMarkerIndex = getSubMarkerIndex();
                m_M_Submarker = subMarkerIndex.n_elem;
                m_N = geno.getNnomissing();
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
//			std::cout << "j: " << j << std::endl;	
                        geno.Get_OneSNP_StdGeno(j, &vec);
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

//	std::cout << "m_M_Submarker: " << m_M_Submarker << std::endl;
        CorssProd_usingSubMarker CorssProd_usingSubMarker(bVec);
//	std::cout << "m_M_Submarker: 2 " << m_M_Submarker << std::endl;
  // call paralleReduce to start the work
        parallelReduce(0, m_M_Submarker, CorssProd_usingSubMarker);
//	std::cout << "m_M_Submarker: 3 " << m_M_Submarker << std::endl;
//	std::cout << "CorssProd_usingSubMarker.m_bout " << CorssProd_usingSubMarker.m_bout << std::endl;
  // return the computed product
        //cout << "Msub: " << Msub << endl;
        //for(int i=0; i<100; ++i)
        //{
        //      cout << (CorssProd_usingSubMarker.m_bout/m_M_Submarker)[i] << ' ';
        //}
//        cout << endl;

//	cout << (CorssProd_usingSubMarker.m_bout).n_elem << endl;
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
float calGRMValueforSamplePair(arma::ivec &sampleidsVec){
        //std::vector<float> stdGenoforSamples = geno.Get_Samples_StdGeno(sampleidsVec);
        geno.Get_Samples_StdGeno(sampleidsVec);
	//std::cout << "here5" << std::endl;
	//for(int i = 0; i < 10; i++){
	//	std::cout << geno.stdGenoforSamples[i] << " ";
	//}
	//std::cout << std::endl;
	//std::cout << geno.stdGenoforSamples.size() << std::endl;
        int Ntotal = geno.getNnomissing();
        float grmValue;
	std::vector<float> stdGenoforSamples2;
	//std::cout << "here5b" << std::endl;
	//std::cout << sampleidsVec.n_elem << std::endl;
	//std::cout << "here5c" << std::endl;
        if(sampleidsVec.n_elem == 2){
                std::vector<float> s1Vec;
                //s1Vec.zeros(Ntotal);

                std::vector<float> s2Vec;
                //arma::fvec s2Vec;
                //s2Vec.zeros(Ntotal);

                for(int i = 0; i < Ntotal; i++){
                        //s1Vec[i] = stdGenoforSamples[i*2+0];
                        s1Vec.push_back(geno.stdGenoforSamples[i*2]);
                        s2Vec.push_back(geno.stdGenoforSamples[i*2+1]);
                }
                grmValue = parallelInnerProduct(s1Vec, s2Vec);
                //grmValue = innerProductFun(s1Vec, s2Vec);
        }else{
	//	std::cout << "here5d" << std::endl;
	//	std::cout << "geno.stdGenoforSamples.size() " << geno.stdGenoforSamples.size() << std::endl;
		stdGenoforSamples2.clear();
		for (int i=0; i< geno.stdGenoforSamples.size(); i++){
			//std::cout << i << " " << geno.stdGenoforSamples[i] << " ";
        		stdGenoforSamples2.push_back(geno.stdGenoforSamples[i]);
		}
	//	std::cout << std::endl;
	//	std::cout << "here6" << std::endl;
                grmValue = parallelInnerProduct(stdGenoforSamples2, geno.stdGenoforSamples);
                //grmValue = innerProductFun(stdGenoforSamples2, geno.stdGenoforSamples);
	//	std::cout << "here7" << std::endl;
        }
        return(grmValue);
}


//Rcpp::List createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
//arma::sp_fmat createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){



// [[Rcpp::export]]
Rcpp::List createSparseKin(arma::fvec& markerIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){

        int nSubMarker = markerIndexVec.n_elem;
        int Ntotal = geno.getNnomissing();
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
	//std::cout << "createSparseKin1b" << std::endl;

	Get_MultiMarkersBySample_StdGeno(markerIndexVec, stdGenoMultiMarkers);
	std::cout << "createSparseKin2" << std::endl;
	//arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), Ntotal, nSubMarker);
	arma::fmat stdGenoMultiMarkersMat(&stdGenoMultiMarkers.front(), nSubMarker, Ntotal);
	//std::cout << "createSparseKin3" << std::endl;
	//std::cout << "stdGenoMultiMarkersMat.n_rows: " << stdGenoMultiMarkersMat.n_rows << std::endl;
	//std::cout << "stdGenoMultiMarkersMat.n_cols: " << stdGenoMultiMarkersMat.n_cols << std::endl;



        for(unsigned int i=0; i< Ntotal; i++){
              for(unsigned int j = i; j < Ntotal; j++){
                        //kinValueTemp = arma::dot(stdGenoMultiMarkersMat.row(i), stdGenoMultiMarkersMat.row(j));
			if(j > i){
                		kinValueTemp = arma::dot(stdGenoMultiMarkersMat.col(i), stdGenoMultiMarkersMat.col(j));
                		kinValueTemp = kinValueTemp/nSubMarker;
                		if(kinValueTemp >= relatednessCutoff){
//                              if(i == 0){
                                //std::cout << "kinValueTemp: " << kinValueTemp << std::endl;
                                //std::cout << "relatednessCutoff: " << relatednessCutoff << std::endl;
                                //std::cout << "i: " << i << std::endl;
//                              std::cout << "j: " << j;
//                              }
                        		iIndexVec.push_back(i);
					jIndexVec.push_back(j);

                		}
			}else{
				iIndexVec.push_back(i);
				jIndexVec.push_back(j);
			}
        	}
	}
	
	arma::fvec * temp = &(geno.m_OneSNP_StdGeno);
        size_t ni = iIndexVec.size();
        kinValueVec.resize(ni);
        std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

        int Mmarker = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
		//geno.getM();
        for(size_t i=0; i< Mmarker; i++){
                geno.Get_OneSNP_StdGeno(i, temp);
                for(size_t j=0; j < ni; j++){
                        kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
                }

        }


	for(size_t j=0; j < ni; j++){
		if(kinValueVec[j] >= relatednessCutoff){
	//	std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
			kinValueVec[j] = tauVec(1)*kinValueVec[j];
			iIndexVec2.push_back(iIndexVec[j]+1);
			jIndexVec2.push_back(jIndexVec[j]+1);
			if(iIndexVec[j] == jIndexVec[j]){
				kinValueVec[j] = kinValueVec[j] + tauVec(0)/(wVec(iIndexVec[j]));	
			}
			kinValueVec2.push_back(kinValueVec[j]);
		}

	}

//	std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;

	//arma::fvec x(kinValueVec2);
	//arma::umat locations(iIndexVec2);
	//arma::uvec jIndexVec2_b(jIndexVec2);
	//locations.insert_cols(locations.n_cols, jIndexVec2_b); 
	//arma::umat locationst = locations.t();
	//locations.clear();
	
	//create a sparse Sigma
//	arma::sp_fmat sparseSigma(locationst, x);
//	arma::sp_fmat sparseSigmab  = arma::symmatu(sparseSigma);
	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);
//	return sparseSigmab;
}



// [[Rcpp::export]]
arma::fmat getColfromStdGenoMultiMarkersMat(arma::uvec & a){
	return((geno.stdGenoMultiMarkersMat).cols(a));
}

// [[Rcpp::export]]
int getNColStdGenoMultiMarkersMat(){
	return((geno.stdGenoMultiMarkersMat).n_cols);
}

// [[Rcpp::export]]
int getNRowStdGenoMultiMarkersMat(){
        return((geno.stdGenoMultiMarkersMat).n_rows);
}


// [[Rcpp::export]]
void setSubMarkerIndex(arma::ivec &subMarkerIndexRandom){
	geno.subMarkerIndex = subMarkerIndexRandom;
//	std::cout << "(geno.subMarkerIndex).n_elem: " << (geno.subMarkerIndex).n_elem << std::endl;
	int Nnomissing = geno.getNnomissing();
	(geno.stdGenoMultiMarkersMat).set_size(subMarkerIndexRandom.n_elem, Nnomissing);
}

// [[Rcpp::export]]
void setRelatednessCutoff(float a){
	geno.relatednessCutoff = a;
}


// [[Rcpp::export]]
double innerProduct(NumericVector x, NumericVector y) {
   return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}


//Rcpp::List refineKin(std::vector<unsigned int> &iIndexVec, std::vector<unsigned int> & jIndexVec, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
//Rcpp::List refineKin(arma::imat &iMat, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){

// [[Rcpp::export]]
Rcpp::List refineKin(float relatednessCutoff){
        std::vector<unsigned int>     iIndexVec2;
        std::vector<unsigned int>     jIndexVec2;
//	std::vector<float>     kinValueVec;
        std::vector<float>     kinValueVec2;
 //       std::vector<float>     kinValueVec_orig; //for test original kinship

	arma::fvec * temp = &(geno.m_OneSNP_StdGeno);
	(*temp).clear();
        //size_t ni = iIndexVec.size();
        //size_t ni = iMat.n_rows;
        size_t ni = geno.indiceVec.size();
	std::cout << "ni: " << ni << std::endl;
 
	initKinValueVecFinal(ni);

//	std::cout << "OKK: "  << std::endl;
//        kinValueVec.resize(ni);
//        std::fill(kinValueVec.begin(), kinValueVec.end(), 0);

        //int Mmarker = geno.getM();
        int Mmarker = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM(); 

        //for(size_t i=0; i< Mmarker; i++){
        //        geno.Get_OneSNP_StdGeno(i, temp);
        //        for(size_t j=0; j < ni; j++){
        //                kinValueVec[j] = kinValueVec[j] + (((*temp)[iIndexVec[j]])*((*temp)[jIndexVec[j]]))/Mmarker;
        //        }
        //}
	//arma::fvec kinValueVecTemp;
	arma::fvec kinValueVecTemp2;
	arma::fvec GRMvec;
	GRMvec.set_size(ni);
	//int Mmarker_mafgr1perc = 0;
  	for(size_t i=0; i< Mmarker; i++){
//		std::cout << "OKKK: "  << std::endl;
//		std::cout << "Mmarker: " << std::endl;

//                geno.Get_OneSNP_StdGeno(i, temp);
		float freqv = geno.alleleFreqVec[i];
		//if(freqv >= minMAFtoConstructGRM && freqv <= 1-minMAFtoConstructGRM){
		//Mmarker_mafgr1perc = Mmarker_mafgr1perc + 1;

                geno.Get_OneSNP_Geno(i);
		float invstdv = geno.invstdvVec[i];
		geno.setSparseKinLookUpArr(freqv, invstdv);			

		//std::cout << "freqv: " << freqv << std::endl;
		//std::cout << "invstdv: " << invstdv << std::endl;
		//for (int j = 0; j < 3; j++){
		//	std::cout << geno.sKinLookUpArr[j][0] << std::endl;	
		//	std::cout << geno.sKinLookUpArr[j][1] << std::endl;	
		//	std::cout << geno.sKinLookUpArr[j][2] << std::endl;	

		//}
		//std::cout << "geno.m_OneSNP_StdGeno(i) " << geno.m_OneSNP_StdGeno(i) <<  std::endl;	
		//kinValueVecTemp = parallelcalsparseGRM(iMat);
//		parallelcalsparseGRM(iMat, GRMvec);

		parallelcalsparseGRM(GRMvec);
		//std::cout << "kinValueVecTemp.n_elem: " << kinValueVecTemp.n_elem << std::endl;
//		std::cout << "OKKK2: "  << std::endl;
		parallelsumTwoVec(GRMvec);
//		for(size_t j=0; j< ni; j++){
//			(geno.kinValueVecFinal)[j] = (geno.kinValueVecFinal)[j] + GRMvec(j);
//		}
		(*temp).clear();
	   //}//if(freqv >= 0.01 && freqv <= 0.99){
		//kinValueVecTemp.clear();
        }



       // for(size_t j=0; j < 100; j++){
       //         std::cout << "iIndexVec[j]: " << iIndexVec[j] << std::endl;
       //         std::cout << "jIndexVec[j]: " << jIndexVec[j] << std::endl;
       //         std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
       // }

	int a1;
	int a2;
        for(size_t j=0; j < ni; j++){
		geno.kinValueVecFinal[j] = (geno.kinValueVecFinal[j]) /(Mmarker);

//		std::cout << "j: " << j << " geno.kinValueVecFinal[j]: " << geno.kinValueVecFinal[j] << std::endl;
            //    if(geno.kinValueVecFinal[j] >= relatednessCutoff){
                if((geno.kinValueVecFinal[j]) >= relatednessCutoff){
        //      std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
			//kinValueVec_orig.push_back((geno.kinValueVecFinal)[j]); //for test	
                        //(geno.kinValueVecFinal)[j] = tauVec(1)*(geno.kinValueVecFinal)[j];
                        //(geno.kinValueVecFinal)[j] = tauVec(1)*(geno.kinValueVecFinal)[j];
 				 a1 = (geno.indiceVec)[j].first + 1;
				 a2 = (geno.indiceVec)[j].second + 1;
				 iIndexVec2.push_back(a1);
				 jIndexVec2.push_back(a2);

                        kinValueVec2.push_back((geno.kinValueVecFinal)[j]);
                }

        }



	std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;
	//return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2,  Named("kinValue_orig") = kinValueVec_orig);	
	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);	
}


// [[Rcpp::export]]
Rcpp::List shortenList(arma::imat &iMat, arma::fvec &kinValueVecTemp, float relatednessCutoff, arma::fvec& wVec,  arma::fvec& tauVec){
	        std::vector<unsigned int>     iIndexVec2;
        std::vector<unsigned int>     jIndexVec2;
	std::vector<float>     kinValueVec2;
	size_t ni = iMat.n_rows;

	for(size_t j=0; j < ni; j++){
                if(kinValueVecTemp(j) >= relatednessCutoff){
        //      std::cout << "kinValueVec[j]: " << kinValueVec[j] << std::endl;
                        kinValueVecTemp(j) = tauVec(1)*(kinValueVecTemp(j));
                        iIndexVec2.push_back(iMat(j,1)+1);
                        //iIndexVec2.push_back(iIndexVec[j]+1);
                        jIndexVec2.push_back(iMat(j,2)+1);
                        //jIndexVec2.push_back(jIndexVec[j]+1);
        //                if(iIndexVec[j] == jIndexVec[j]){
        //                        kinValueVec[j] = kinValueVec[j] + tauVec(0)/(wVec(iIndexVec[j]));
        //                }

                        if(iMat(j,1) == iMat(j,2)){
                                kinValueVecTemp(j) = kinValueVecTemp(j) + tauVec(0)/(wVec(iMat(j,1)));
                        }

                        kinValueVec2.push_back(kinValueVecTemp(j));
                }

        }

        std::cout << "kinValueVec2.size(): " << kinValueVec2.size() << std::endl;
	return Rcpp::List::create(Named("iIndex") = iIndexVec2, Named("jIndex") = jIndexVec2, Named("kinValue") = kinValueVec2);

}

// [[Rcpp::export]]
arma::fvec testTime(int i, arma::fcolvec & m_bVec){
	arma::fvec vec;
	arma::fvec mvec;
	std::cout << "i is " << i << std::endl;
	clock_t t_0;
	t_0 = clock();
        geno.Get_OneSNP_StdGeno(i, &vec);
	clock_t t_1;
	t_1 = clock();
	std::cout << "t_1-t_0 is " << t_1-t_0 << std::endl;
        float val1 = dot(vec,  m_bVec);
	clock_t t_2;
	t_2 = clock();
	std::cout << "t_2-t_1 is " << t_2-t_1 << std::endl;
        mvec = val1 * (vec);
	clock_t t_3;
	t_3 = clock();
	std::cout << "t_3-t_2 is " << t_3-t_2 << std::endl;
	return(mvec);
}


// [[Rcpp::export]]
arma::sp_mat gen_sp_v2(const arma::sp_mat& a) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a);
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;

    return result;
}


// [[Rcpp::export]]
arma::vec gen_spsolve_v2(const arma::sp_mat& a) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a);
    int r = result.n_rows;
    arma::vec y = arma::linspace<arma::vec>(0, 5, r);	
    //arma::sp_fmat A = sprandu<sp_fmat>(100, 200, 0.1);
    //arma::sp_mat result1 = result * A;
    arma::vec x = arma::spsolve( result, y ); 
    	
    return x;
}

// [[Rcpp::export]]
arma::vec gen_spsolve_inR(const arma::sp_mat& a, arma::vec & y) {
    // sparse x sparse -> sparse
    //arma::sp_mat result1 = result * A;
    arma::vec x = arma::spsolve( a, y );

    return x;
}

// [[Rcpp::export]]
arma::fvec get_DiagofKin(){
    //int M = geno.getM();
    int Nnomissing = geno.getNnomissing();
        //cout << "MminMAF=" << MminMAF << endl;
        //cout << "M=" << M << endl; 


    arma::fvec x(Nnomissing);

    if(!(geno.setKinDiagtoOne)){
           x  = (*geno.Get_Diagof_StdGeno());
    	   int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
           x = x/MminMAF; 
    }else{
	   x  = arma::ones<arma::fvec>(Nnomissing);	
    }	
    return(x);
}





//The code below is modified from http://gallery.rcpp.org/articles/parallel-inner-product/
struct stdgenoVectorScalorProduct : public Worker
{
   // source vectors
   arma::fvec & m_bout;
   float  y;
   //unsigned int m_N;
   int jthMarker;

   // constructors
   stdgenoVectorScalorProduct(const int jth, const float y, arma::fvec & prodVec)
      : jthMarker(jth), y(y), m_bout(prodVec) {
        //m_N = geno.getNnomissing();
//      m_bout.zeros(m_N);

  }


   // process just the elements of the range I've been asked to

        void operator()(std::size_t begin, std::size_t end) {
                arma::fvec vec;
                geno.Get_OneSNP_StdGeno(jthMarker, &vec);
                for(unsigned int i = begin; i < end; i++){
                        m_bout[i] = m_bout[i]+vec[i] * y;
                }
        }



};


// [[Rcpp::export]]
void getstdgenoVectorScalorProduct(int jth, float y, arma::fvec & prodVec) {


   stdgenoVectorScalorProduct stdgenoVectorScalorProduct(jth, y, prodVec);

   unsigned int m_N = geno.getNnomissing();

   parallelFor(0, m_N, stdgenoVectorScalorProduct);

   // return the computed product
}





struct getP_mailman : public Worker
{
        // source vectors
        unsigned int ithMarker;
	unsigned int powVal;
	arma::ivec ithGeno;
        // destination vector
        arma::ivec Psubvec;


        // constructors
        getP_mailman(unsigned int ith, unsigned int mmchunksize)
                : ithMarker(ith){
		ithGeno = Get_OneSNP_Geno(ith);			
		//unsigned int k  = pow(3, mmchunksize);
		unsigned m_N = geno.getNnomissing();
		Psubvec.zeros(m_N);
		unsigned int powNumber = mmchunksize - 1 - ith % mmchunksize; 		
		powVal = pow(3, powNumber);
        }


	// take the square root of the range of elements requested
     void operator()(std::size_t begin, std::size_t end) {

	for(unsigned int j = begin; j < end; j++){
		Psubvec[j] = ithGeno[j] * powVal;
        }		 
     }

};


int computePindex(arma::ivec &ithGeno){
	int a = ithGeno.n_elem;
	int q = 0;
	int baseNum;
	for(unsigned int i = 0; i < a; i++){
		baseNum = pow(3, a - i - 1);
		q = q + ithGeno[i] * baseNum;
	}
	return(q);
}


struct getP_mailman_NbyM : public Worker
{
        // source vectors
        unsigned int jthChunk;
        unsigned int mmchunksize;
        //arma::ivec ithGeno;
        // destination vector
        arma::ivec Psubvec;


        // constructors
        getP_mailman_NbyM(unsigned int jthChunk,unsigned int mmchunksize)
                : jthChunk(jthChunk), mmchunksize(mmchunksize){
                //ithGeno = Get_OneSNP_Geno(ith);
                //unsigned int k  = pow(3, mmchunksize);
                unsigned m_M = geno.getM();
                Psubvec.zeros(m_M);
                //powNumber = mmchunksize - 1 - ith % mmchunksize;
        }


        // take the square root of the range of elements requested
     void operator()(std::size_t begin, std::size_t end) {
	arma::ivec ithGeno;	
	arma::ivec ithGenosub;
	unsigned int jthIndvStart = jthChunk * mmchunksize;
	unsigned int jthIndvEnd = (jthChunk+1) * mmchunksize - 1;
	arma::uvec indvIndex = arma::linspace<arma::uvec>(jthIndvStart, jthIndvEnd);
        for(unsigned int i = begin; i < end; i++){
		ithGeno = Get_OneSNP_Geno(i);		
		ithGenosub = ithGeno.elem(indvIndex);
		Psubvec[i] = computePindex(ithGenosub);
        }
     }

};



// // [[Rcpp::export]]
//arma::ivec parallelmmGetP(unsigned int ith, unsigned int mmchunksize) {
  
//  	int M = geno.getM();
//	int N = geno.getNnomissing();	
//  	Pvec.zeros(N);

//  	getP_mailman getP_mailman(ith, mmchunksize);
  
//  	parallelFor(0, N, getP_mailman);
 	
//  	return getP_mailman.Psubvec;
//}


// [[Rcpp::export]]
void sumPz(arma::fvec & Pbvec, arma::fvec & Ubvec, unsigned int mmchunksize){

        for (int i = 0; i < Pbvec.n_elem; i++){
                std::cout << "i: " << i << " " << Pbvec[i] << std::endl;
        }

        unsigned int d = Pbvec.n_elem;;
        Ubvec.zeros(mmchunksize);
        unsigned int i = 0;
        arma::fvec z0;
        arma::fvec z1;
        arma::fvec z2;
        z0.zeros(d/3);
        z1.zeros(d/3);
        z2.zeros(d/3);

        while(i < mmchunksize){
                d = d / 3;
//              std::cout << "d: " << d << std::endl;
                z0.resize(d);
                z1.resize(d);
                z2.resize(d);

//              arma::uvec indexvec = arma::linspace<arma::uvec>(0, d-1);
                z0 = Pbvec.subvec(0, d-1);
/*
                 for (int j = 0; j < z0.n_elem; j++){
                std::cout << "j: " << j << " " << z0[j] << std::endl;
        }
*/
                //indexvec = arma::linspace<arma::uvec>(d, 2*d-1);
                //z1 = Pbvec.elem(indexvec);
                z1 = Pbvec.subvec(d, 2*d-1);
                //indexvec = arma::linspace<arma::uvec>(2*d, 3*d-1);
                //z2 = Pbvec.elem(indexvec);
                z2 = Pbvec.subvec(2*d, 3*d-1);

                Pbvec.resize(d);
                Pbvec = z0 + z1 + z2;
                Ubvec[i] = sum(z1) + 2*sum(z2);
                i = i + 1;
              std::cout << "i: " << i << std::endl;
              std::cout << "Ubvec[i]: " << Ubvec[i] << std::endl;

        }
}



// [[Rcpp::export]]
void mmGetPb_MbyN(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec, arma::fvec & kinbvec) {
	std::cout << "OKKK" << std::endl;
        int M = geno.getM();
        int N = geno.getNnomissing();
	int k = pow(3,mmchunksize);
        arma::ivec Pvec;
	Pvec.zeros(N);
	Pbvec.zeros(k);
	arma::ivec ithGeno;
	ithGeno.ones(N);
	unsigned int Ptemp;
	Ptemp = 1;
	int indL = cthchunk*mmchunksize;
	int indH = (cthchunk+1)*mmchunksize - 1;
	unsigned int j0 = 0;
	//arma::fmat stdGenoMat(mmchunksize, N);
	float ithfreq; 
	float ithinvstd; 
	arma::fvec chunkfreq = geno.alleleFreqVec.subvec(indL, indH);
	arma::fvec chunkinvstd = geno.invstdvVec.subvec(indL, indH); 
	arma::fvec chunkbvec = bvec.subvec(indL, indH); 

	for (int i = indH; i >= indL; i--){
		ithGeno = Get_OneSNP_Geno(i);
		cout << "Ptemp: " << Ptemp << endl;
		//ithfreq = geno.alleleFreqVec(i);
		//ithinvstd = geno.invstdvVec(i);
		Pvec = Pvec + Ptemp * ithGeno; 
		Ptemp = Ptemp * 3;
		//stdGenoMat.row(j) = ithGeno*ithinvstd - 2*ithfreq*ithinvstd;
		//j0 = j0 + 1;

                //unsigned int k  = pow(3, mmchunksize);
                //unsigned m_N = geno.getNnomissing();
                //Psubvec.zeros(m_N);
                //unsigned int powNumber = mmchunksize - 1 - ith % mmchunksize;

	
	//	getP_mailman getP_mailman(i, mmchunksize);
	//	parallelFor(0, N, getP_mailman);
	//	Pvec = Pvec + getP_mailman.Psubvec;
	//	getP_mailman.Psubvec.clear();
  	}
	

	for (int i = 0; i < N; i++){	
//		std::cout << "i: " << i << " " << Pvec[i] << std::endl;	
		Pbvec[Pvec[i]] = Pbvec[Pvec[i]] + bvec[i];
//		std::cout << "Pbvec[Pvec[i]] " << Pbvec[Pvec[i]] << std::endl;
	}

	arma::fvec Gbvectemp;
	sumPz(Pbvec, Gbvectemp, mmchunksize);
	arma::fvec crossKinVec;
	arma::fvec GbvecInvStd = Gbvectemp % chunkinvstd;
        arma::fvec secondTerm = 2*chunkfreq % chunkinvstd * sum(chunkbvec);
        crossKinVec  = GbvecInvStd - secondTerm;

	//getstdgenoVectorScalorProduct(j, crossKinVec[j], kinbvec);
	j0 = 0;
	arma::fvec stdvec;
	for (int i = indL; i <= indH; i++){
		geno.Get_OneSNP_StdGeno(i, &stdvec);
                kinbvec = kinbvec + crossKinVec[j0]*(stdvec);
		j0 = j0 + 1;
	}

//	for (int i = 0; i < k; i++){
//                std::cout << "Pbvec[i]: " << i << " " << Pbvec[i] << std::endl;
//        }

        //return Pbvec;
}

// [[Rcpp::export]]
void mmGetPb_NbyM(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec) {

        int M = geno.getM();
        int N = geno.getNnomissing();
        int k = pow(3,mmchunksize);
        arma::ivec Pvec;
        Pvec.zeros(M);
        Pbvec.zeros(k);
	getP_mailman_NbyM getP_mailman_NbyM(cthchunk,mmchunksize);
	parallelFor(0, M, getP_mailman_NbyM);
	Pvec = getP_mailman_NbyM.Psubvec;
	for (int i = 0; i < M; i++){
		Pbvec[Pvec[i]] = Pbvec[Pvec[i]] + bvec[i];
	}
}



// [[Rcpp::export]]
void muliplyMailman(arma::fvec & bvec, arma::fvec & Gbvec, arma::fvec & kinbvec){
	int M = geno.getM();
        int N = geno.getNnomissing();

        Gbvec.zeros(M);
	std::cout << "Gbvec.n_elem " << Gbvec.n_elem << std::endl;
        unsigned int mmchunksize = ceil(log(N)/log(3));
	std::cout << "mmchunksize " << mmchunksize << std::endl;

        int numchunk = M / mmchunksize; 
	std::cout << "numchunk " << numchunk << std::endl;
        int reschunk = M % mmchunksize;
	std::cout << "reschunk " << reschunk << std::endl;
	//unsigned int indL;
	//unsigned int indH;
	//mmGetPb(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec)
	arma::fvec Pbvec;
	//arma::fvec Gbvectemp;	


	
	//for (unsigned int j = 0; j < 1; j++){
	for (unsigned int j = 0; j < numchunk; j++){
//		std::cout << "j: " << j << std::endl;
		//Pbvec.zeros(M);
		//indL = j*mmchunksize;
		//indH = (j+1)*mmchunksize-1;
//		if(j == 0){
		double wall0ain = get_wall_time();
 		double cpu0ain  = get_cpu_time();
//		}
//		mmGetPb_MbyN(j, mmchunksize, bvec, Pbvec);

//		if(j == 0){

		mmGetPb_MbyN(j, mmchunksize, bvec, Pbvec, kinbvec);



	double wall1ain = get_wall_time();
 double cpu1ain  = get_cpu_time();
 cout << "Wall Time in mmGetPb_MbyN = " << wall1ain - wall0ain << endl;
 cout << "CPU Time  in mmGetPb_MbyN = " << cpu1ain - cpu0ain  << endl;


//}

//		sumPz(Pbvec, Gbvectemp, mmchunksize);

//if(j == 0){
cout << "ith chunk " << j << endl;
//}
		//getstdgenoVectorScalorProduct(int jth, float y, arma::fvec & prodVec)

//		Gbvec.subvec(j*mmchunksize, (j+1)*mmchunksize-1) = Gbvectemp;
  	}

        if(reschunk > 0){
			arma::fvec vec;
		//arma::uvec indexvec = arma::linspace<arma::uvec>(M-reschunk, M-1);
		for (unsigned int j = M-reschunk; j < M; j++){
     		           geno.Get_OneSNP_StdGeno(j, &vec);
			kinbvec = kinbvec + arma::dot(vec, bvec) * vec;
		}	
        }

	kinbvec = kinbvec / M;
}


// [[Rcpp::export]]
void muliplyMailman_NbyM(arma::fvec & bvec, arma::fvec & tGbvec){
        int M = geno.getM();
        int N = geno.getNnomissing();

        tGbvec.zeros(N);

        unsigned int mmchunksize = ceil(log(M)/log(3));

        int numchunk = N / mmchunksize;
        int reschunk = N % mmchunksize;
        unsigned int indL;
        unsigned int indH;
        //mmGetPb(unsigned int cthchunk, unsigned int mmchunksize, arma::fvec & bvec, arma::fvec & Pbvec)
        arma::fvec Pbvec;
        Pbvec.zeros(M);
	arma::fvec tGbvectemp;

        for (unsigned int j = 0; j < numchunk; j++){
                indL = j*mmchunksize;
                indH = (j+1)*mmchunksize-1;
		mmGetPb_NbyM(j, mmchunksize, bvec, Pbvec);
           	sumPz(Pbvec, tGbvectemp, mmchunksize);
                tGbvec.subvec(j*mmchunksize, (j+1)*mmchunksize-1) = tGbvectemp;     
        }

        if(reschunk > 0){
		arma::imat A(reschunk,M);
		A.zeros();
		arma::ivec Gtemp(N);
		Gtemp.zeros();
		arma::ivec Gtemp2(reschunk);
		Gtemp2.zeros();
		arma::uvec indexvec = arma::linspace<arma::uvec>(M - reschunk -1, M);
                for (unsigned int j = 0; j < M; j++){
			Gtemp = Get_OneSNP_Geno(j);
			Gtemp2 = Gtemp.elem(indexvec);
			A.col(j) = Gtemp2;
                }

		Pbvec.elem(indexvec) = Gtemp2 * (bvec.elem(indexvec));
        }
}

// [[Rcpp::export]]
void freqOverStd(arma::fcolvec& freqOverStdVec){
	freqOverStdVec = 2 * (geno.alleleFreqVec) % (geno.invstdvVec);

	 //int M = geno.getM();
/*
	for (unsigned int j = 0; j < M; j++){
		std::cout << "geno.alleleFreqVec " << j << " " << geno.alleleFreqVec[j] << std::endl; 
		std::cout << "geno.invstdvVec " << j << " " << geno.invstdvVec[j] << std::endl; 
		std::cout << "freqOverStdVec " << j << " " << freqOverStdVec[j] << std::endl; 
               }
*/

}
// [[Rcpp::export]]
arma::fvec getCrossprodMatAndKin_mailman(arma::fcolvec& bVec){
	std::cout << "b0: " << std::endl;
	int M = geno.getM();
        int N = geno.getNnomissing();
	arma::fvec Gbvec;


	double wall0in = get_wall_time();
 	double cpu0in  = get_cpu_time();
 	arma::fvec kinbvec;
        kinbvec.zeros(N);

	muliplyMailman(bVec, Gbvec, kinbvec);


double wall1in = get_wall_time();
 double cpu1in  = get_cpu_time();
 cout << "Wall Time in muliplyMailman = " << wall1in - wall0in << endl;
 cout << "CPU Time  in muliplyMailman = " << cpu1in - cpu0in  << endl;



//	for (unsigned int j = 0; j < M; j++){
//                std::cout << "Gbvec " << j << " " << Gbvec[j] << std::endl;
//               }
/*
//	std::cout << "b: " << std::endl;
	arma::fvec freqOverStdVec;
//	std::cout << "a: " << std::endl;
	freqOverStd(freqOverStdVec);
//	std::cout << "c: " << std::endl;
	arma::fvec crossKinVec;
	arma::fvec GbvecInvStd = Gbvec % (geno.invstdvVec);
	arma::fvec secondTerm = freqOverStdVec * sum(bVec);
	crossKinVec  = GbvecInvStd - secondTerm;

double wall2in = get_wall_time();
 double cpu2in  = get_cpu_time();
 cout << "Wall Time in Gtb = " << wall2in - wall1in << endl;
 cout << "CPU Time  in Gtb = " << cpu2in - cpu1in  << endl;


	 for (unsigned int j = 0; j < M; j++){
                std::cout << "GbvecInvStd " << j << " " << GbvecInvStd[j] << std::endl;
                std::cout << "secondTerm " << j << " " << secondTerm[j] << std::endl;
		std::cout << "crossKinVec " << j << " " << crossKinVec[j] << std::endl;
               }
*/
/*	
	arma::fvec kinbvec;
	kinbvec.zeros(N);

	for (unsigned int j = 0; j < M; j++){
		getstdgenoVectorScalorProduct(j, crossKinVec[j], kinbvec);
	}


double wall3in = get_wall_time();
 double cpu3in  = get_cpu_time();
 cout << "Wall Time in getstdgenoVectorScalorProduct = " << wall3in - wall2in << endl;
 cout << "CPU Time  in getstdgenoVectorScalorProduct = " << cpu3in - cpu2in  << endl;



	kinbvec = kinbvec / M;
*/	
        return(kinbvec);

}

// [[Rcpp::export]]
arma::fvec get_GRMdiagVec(){
  int mMarker = gettotalMarker(); 
  int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
        //cout << "MminMAF=" << MminMAF << endl;

  arma::fvec diagGRMVec = (*geno.Get_Diagof_StdGeno())/MminMAF;
  return(diagGRMVec);
}

// [[Rcpp::export]]
void setminMAFforGRM(float minMAFforGRM){
  minMAFtoConstructGRM = minMAFforGRM;
}

// [[Rcpp::export]]
void setmaxMissingRateforGRM(float maxMissingforGRM){
  geno.maxMissingRate = maxMissingforGRM;
}


// [[Rcpp::export]]
void set_Diagof_StdGeno_LOCO(){

      
  int Nnomissing = geno.getNnomissing();
  int chrlength = geno.startIndexVec.n_elem;
  (geno.mtx_DiagStd_LOCO).zeros(Nnomissing, chrlength);
  (geno.Msub_MAFge_minMAFtoConstructGRM_byChr).zeros(chrlength);
//  std::cout << "debug1" << std::endl;
    int starti, endi;
    arma::fvec * temp = &geno.m_OneSNP_StdGeno;
for(size_t k=0; k< chrlength; k++){
   starti = geno.startIndexVec[k];
   endi = geno.endIndexVec[k];
//  std::cout << "debug2" << std::endl;
  if((starti != -1) && (endi != -1)){
  	for(int i=starti; i<= endi; i++){
         		geno.Get_OneSNP_StdGeno(i, temp);
	 		(geno.mtx_DiagStd_LOCO).col(k) = (geno.mtx_DiagStd_LOCO).col(k) + (*temp) % (*temp);
	 		geno.Msub_MAFge_minMAFtoConstructGRM_byChr[k] = geno.Msub_MAFge_minMAFtoConstructGRM_byChr[k] + 1;

  	}
  (geno.mtx_DiagStd_LOCO).col(k) = *geno.Get_Diagof_StdGeno() -  (geno.mtx_DiagStd_LOCO).col(k);
  }
}	
}

/*
// [[Rcpp::export]]
void setminMAC_VarianceRatio(arma::fvec  t_cateVarRatioMinMACVecExclude, arma::fvec  t_cateVarRatioMaxMACVecInclude){
  g_cateVarRatioMinMACVecExclude = t_cateVarRatioMinMACVecExclude;
  g_cateVarRatioMaxMACVecInclude = t_cateVarRatioMaxMACVecInclude;
}
*/


// [[Rcpp::export]]
void setminMAC_VarianceRatio(float t_minMACVarRatio, float t_maxMACVarRatio, bool t_isVarianceRatioinGeno){ 
	geno.g_minMACVarRatio = t_minMACVarRatio;
	geno.g_maxMACVarRatio = t_maxMACVarRatio;
	geno.isVarRatio = t_isVarianceRatioinGeno;
	std::cout << "geno.g_minMACVarRatio " << geno.g_minMACVarRatio << std::endl;
	std::cout << "geno.g_maxMACVarRatio " << geno.g_maxMACVarRatio << std::endl;	
}

// // [[Rcpp::export]] 
//int getNumofMarkersforGRM(){
//  int a = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
//  return(a);
//}
// [[Rcpp::export]]
Rcpp::List getCoefficients_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, int maxiterPCG, float tolPCG, bool LOCO){

        int Nnomissing = geno.getNnomissing();
        arma::fvec Sigma_iY;
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
        arma::fmat Sigma_iXt = Sigma_iX.t();
        arma::fvec SigmaiXtY = Sigma_iXt * Yvec;
        arma::fvec alpha = cov * SigmaiXtY;

        arma::fvec eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
        return Rcpp::List::create(Named("Sigma_iY") = Sigma_iY, Named("Sigma_iX") = Sigma_iX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
}


// [[Rcpp::export]]
arma::fvec GetTrace_largeMem_multiV(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::ivec & fixtauVec, arma::fmat& cov1,  int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO){		
		
  	set_seed(200);
	
	int q2 = arma::sum(fixtauVec==0);
        arma::uvec idxtau = arma::find(fixtauVec==0);

	idxtau.print("idxtau");

  	arma::fmat Sigma_iXt = Sigma_iX.t();
        int Nnomissing = geno.getNnomissing();
	unsigned int k = Kmat_vec.size();
        unsigned int k1 = k + 2;
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
        NumericVector uVec0;

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
	while( !isConverge ){
    		for(int i = nrun_trace_start; i < nrun_trace_end; i++){

    			uVec0 = nb(Nnomissing);
    			uVec = as<arma::fvec>(uVec0);
    			uVec = uVec*2 - 1;
 			Sigma_iu = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, uVec, maxiterPCG, tolPCG, LOCO);
			Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
			
			if(fixtauVec(0) == 0)	{
        			Au_mat.col(0) = uVec;
        			temp_mat(i,0) = dot(Au_mat.col(0), Pu);
			}
			//std::cout << "okkkkk" << std::endl;
        		// conversion for ops with sp_mat
			if(fixtauVec(1) == 0)   {
				if(!LOCO){
					temp_vec_double = getCrossprodMatAndKin(uVec);
				}else{
					temp_vec_double = getCrossprodMatAndKin_LOCO(uVec);	
				}	
				Au_mat.col(1) = temp_vec_double;
        			temp_mat(i,1) = dot(temp_vec_double, Pu);
			}
			//std::cout << "okkkkk2" << std::endl;


    			for(int j=2; j<k1;j++){
				if(fixtauVec(j) == 0){
    					Au_mat.col(j) = 0.0+Kmat_vec[j-2] * uVec;
    					temp_mat(i,j) = dot(Au_mat.col(j), Pu);
				}
    			} // end for j in 2:k1
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
arma::fvec GetTrace_multiV(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, arma::fvec& tauVec, arma::ivec & fixtauVec, arma::fmat& cov1,  int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO){

        set_seed(200);

        int q2 = arma::sum(fixtauVec==0);
        arma::uvec idxtau = arma::find(fixtauVec==0);

        idxtau.print("idxtau");

        arma::fmat Sigma_iXt = Sigma_iX.t();
        int Nnomissing = geno.getNnomissing();
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
        NumericVector uVec0;

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
                        Sigma_iu = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, uVec, maxiterPCG, tolPCG, LOCO);
                        Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));

                        if(fixtauVec(0) == 0)   {
                                Au_mat.col(0) = uVec;
                                temp_mat(i,0) = dot(Au_mat.col(0), Pu);
                        }
                        // conversion for ops with sp_mat
		  if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){		
                        if(fixtauVec(1) == 0)   {
                                if(!LOCO){
                                        temp_vec_double = getCrossprodMatAndKin(uVec);
                                }else{
                                        temp_vec_double = getCrossprodMatAndKin_LOCO(uVec);
                                }
                                Au_mat.col(1) = temp_vec_double;
                                temp_mat(i,1) = dot(temp_vec_double, Pu);
                        }


                        for(int j=2; j<k1;j++){
                                if(fixtauVec(j) == 0){
                                        Au_mat.col(j) = 0.0+Kmat_vec[j-2] * uVec;
                                        temp_mat(i,j) = dot(Au_mat.col(j), Pu);
                                }
                        } // end for j in 2:k1
		   }else{
			Ibvec = g_I_longl_mat.t() * uVec;
                	if(!LOCO){
                        	GRM_I_bvec = getCrossprodMatAndKin(Ibvec);
                	}else{
                        	GRM_I_bvec = getCrossprodMatAndKin_LOCO(Ibvec);
                	}
			if(g_T_longl_mat.n_rows == 0){
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
				Tbvec = g_T_longl_mat.t() * uVec;
				if(!LOCO){
                                	GRM_T_bvec = getCrossprodMatAndKin(Tbvec);
                        	}else{
                                	GRM_T_bvec = getCrossprodMatAndKin_LOCO(Tbvec);
                        	}
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
Rcpp::List getAIScore_largeMem_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, arma::ivec & fixtauVec, 
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO){

	int q2 = arma::sum(fixtauVec==0);
	arma::uvec idxtau = arma::find(fixtauVec==0);
	//int ncol_x = X.n_cols;
	arma::fvec tau0;
	//unsigned int k = Kmat_vec.size();
	//unsigned int k1 = k + 2;	
	unsigned int k1 = idxtau.size();
	arma::fmat AI(k1,k1);
	arma::fvec YPAPY(k1);
	YPAPY.zeros();
        arma::fvec Trace(k1);
	Trace.zeros();
	std::cout << "k1 " << k1 << std::endl;
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
		for(int i=0; i<k1; i++){
		    if(fixtauVec(i) == 0){	
			if(i==0){
				APY = PY1;
			}else if(i==1){
				if(!LOCO){
					APY = getCrossprodMatAndKin(PY1);
				}else{
					APY = getCrossprodMatAndKin_LOCO(PY1);
				}	
			}else{
				APY = Kmat_vec[i-2]*PY1;
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
	AI.print("AI");	
	arma::fmat AI_update = AI.submat(idxtau, idxtau);
	arma::fvec YPAPY_update = YPAPY.elem(idxtau);

	//vector with length=q2
	Trace = GetTrace_multiV(Sigma_iX, Xmat, wVec, tauVec, fixtauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO);
	YPAPY_update.print("YPAPY_update");
	Trace.print("Trace");	
        //arma::fvec PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG);
        //arma::fvec PAPY = PAPY_1 - Sigma_iX * (cov1 * (Sigma_iXt * PAPY_1));
	return Rcpp::List::create(Named("YPAPY") = YPAPY_update, Named("Trace") = Trace,Named("PY") = PY1,Named("AI") = AI_update);
}




// [[Rcpp::export]]
Rcpp::List getAIScore_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec& wVec,  arma::fvec& tauVec, arma::ivec & fixtauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff, bool LOCO){

        std::cout << "getAIScore_multiV " << getAIScore_multiV << std::endl;
        int q2 = arma::sum(fixtauVec==0);
        arma::uvec idxtau = arma::find(fixtauVec==0);
        //int ncol_x = X.n_cols;
        arma::fvec tau0;
        //unsigned int k = Kmat_vec.size();
        //unsigned int k1 = k + 2;
        //unsigned int k1 = idxtau.size();
	unsigned int k1 = g_num_Kmat;
        arma::fmat AI(k1,k1);
        arma::fvec YPAPY(k1);
        YPAPY.zeros();
        arma::fvec Trace(k1);
        Trace.zeros();
        std::cout << "k1 " << k1 << std::endl;
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
                                if(!LOCO){
                                        APY = getCrossprodMatAndKin(PY1);
                                }else{
                                        APY = getCrossprodMatAndKin_LOCO(PY1);
                                }
                        }else{
				if(Kmat_vec.size() > 0){
                                	APY = Kmat_vec[i-2]*PY1;
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

		if(!LOCO){
                        GRM_I_bvec = getCrossprodMatAndKin(Ibvec);
                }else{
                        GRM_I_bvec = getCrossprodMatAndKin_LOCO(Ibvec);
                }

		std::cout << "g_I_longl_mat.n_rows " << g_I_longl_mat.n_rows << std::endl;
	       if(g_T_longl_mat.n_rows == 0){
		  for(int i=0; i<k1; i++){
                    if(fixtauVec(i) == 0){
                        if(i==0){
                                APY = PY1;
                        }else if(i==1){
				std::cout << "GRM_I_bvec.n_elem " << GRM_I_bvec.n_elem << std::endl;
                                //if(!LOCO){
                                //        APY = getCrossprodMatAndKin(PY1);
                                //}else{
                                //        APY = getCrossprodMatAndKin_LOCO(PY1);
                                //}
				APY = g_I_longl_mat * GRM_I_bvec;
                        }else if (i == 2){
				APY = g_I_longl_mat * Ibvec;
                                //APY = Kmat_vec[i-2]*PY1;
                        }else{
				if(Kmat_vec.size() > 0){
					V_I_bvec = Kmat_vec[i-3]*Ibvec;
					APY = g_I_longl_mat * V_I_bvec;				
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
                  if(!LOCO){
                    GRM_T_bvec = getCrossprodMatAndKin(Tbvec);
                  }else{
                    GRM_T_bvec = getCrossprodMatAndKin_LOCO(Tbvec);
                  }

		  unsigned int kmatind = 0; 	
		  for(int i=0; i<k1; i++){
                    if(fixtauVec(i) == 0){
                        if(i==0){
                                APY = PY1;
                        }else if(i==1){
                                //if(!LOCO){
                                //        APY = getCrossprodMatAndKin(PY1);
                                //}else{
                                //        APY = getCrossprodMatAndKin_LOCO(PY1);
                                //}
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
        arma::fmat AI_update = AI.submat(idxtau, idxtau);
        arma::fvec YPAPY_update = YPAPY.elem(idxtau);

        //vector with length=q2
        Trace = GetTrace_multiV(Sigma_iX, Xmat, wVec, tauVec, fixtauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO);
        YPAPY_update.print("YPAPY_update");
        Trace.print("Trace");
        //arma::fvec PAPY_1 = getPCG1ofSigmaAndVector_multiV(wVec, tauVec, APY, maxiterPCG, tolPCG);
        //arma::fvec PAPY = PAPY_1 - Sigma_iX * (cov1 * (Sigma_iXt * PAPY_1));
        return Rcpp::List::create(Named("YPAPY") = YPAPY_update, Named("Trace") = Trace,Named("PY") = PY1,Named("AI") = AI_update);
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
Rcpp::List fitglmmaiRPCG_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec &wVec,  arma::fvec & tauVec, arma::ivec & fixtauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff, bool LOCO){
	Function warning("warning");

	//unsigned int k = Kmat_vec.size();
        //unsigned int k1 = k + 2;
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
	Dtau.print("Dtau");
	score1.print("score1");
	AI1.print("AI1");
	//
	//
	arma::fvec Dtau_k1(k1);
	Dtau_k1.zeros();

	std::cout << "k1 " << k1 << std::endl;
	fixtauVec.print("fixtauVec");

	// fill dtau using dtau_pre, padding 0
	int i2 = 0;
	for(int i=0; i<k1; i++){
		std::cout << "i " << i << std::endl;
		if(fixtauVec(i)==0){ // not fixed
			Dtau_k1(i) = Dtau(i2);
			i2++;
		} 
	} // end for i
	Dtau_k1.print("Dtau_k1");	
	tau0 = tauVec;
	tauVec = tauVec + Dtau_k1;
	arma::fvec tauVecabs;
	arma::ivec fixrhoidx0, fixrhoidx, tauupdateidx;
	//arma::ivec fixrhoidx0, fixrhoidx, covarianceidxVec1, covarianceidxVec_sub1, covarianceidxVec2, covarianceidxVec_sub2,covarianceidxVec3, covarianceidxVec_sub3, tauupdateidx;
	arma::uvec covarianceidxVec1, covarianceidxVec_sub1, covarianceidxVec2, covarianceidxVec_sub2,covarianceidxVec3, covarianceidxVec_sub3;	

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
		tauVecabs = tauVec.elem(covarianceidxVec_sub1) / arma::abs(tauVec.elem(covarianceidxVec_sub1));
		tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
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
		//covarianceidxVec_sub1 = g_covarianceidxMat_col1.elem(arma::find(fixrhoidx == 1));
		//covarianceidxVec_sub2 = g_covarianceidxMat_col2.elem(arma::find(fixrhoidx == 1));
		//covarianceidxVec_sub3 = g_covarianceidxMat_col3.elem(arma::find(fixrhoidx == 1));
		tauVecabs = tauVec.elem(covarianceidxVec_sub1) / arma::abs(tauVec.elem(covarianceidxVec_sub1));
                tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
	}	
	 return List::create(Named("tau") = tauVec, Named("AI") = AI1, Named("score") = score1);
}	




// [[Rcpp::export]]
Rcpp::List fitglmmaiRPCG_largeMem_multiV(arma::fvec& Yvec, arma::fmat& Xmat, arma::fvec &wVec,  arma::fvec & tauVec, arma::ivec & fixtauVec,
arma::fvec& Sigma_iY, arma::fmat & Sigma_iX, arma::fmat & cov,
int nrun, int maxiterPCG, float tolPCG, float tol, float traceCVcutoff, bool LOCO){
        Function warning("warning");

        unsigned int k = Kmat_vec.size();
        unsigned int k1 = k + 2;

        int q2 = arma::sum(fixtauVec==0);
        arma::uvec idxtau = arma::find(fixtauVec==0);
        arma::fvec tau0 = tauVec;


        std::cout << "check 1" << std::endl;
        Rcpp::List re = getAIScore_multiV(Yvec, Xmat,wVec,  tauVec, fixtauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff, LOCO);

        std::cout << "check 2" << std::endl;

        arma::fvec YPAPY = re["YPAPY"];
        arma::fvec Trace = re["Trace"];
        arma::fvec score1 = YPAPY - Trace;
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
        Dtau.print("Dtau");
        score1.print("score1");
        AI1.print("AI1");
        //
        //
        arma::fvec Dtau_k1(k1);
        Dtau_k1.zeros();

        std::cout << "k1 " << k1 << std::endl;
        fixtauVec.print("fixtauVec");

        // fill dtau using dtau_pre, padding 0
        int i2 = 0;
        for(int i=0; i<k1; i++){
                std::cout << "i " << i << std::endl;
                if(fixtauVec(i)==0){ // not fixed
                        Dtau_k1(i) = Dtau(i2);
                        i2++;
                }
        } // end for i
        Dtau_k1.print("Dtau_k1");
        tau0 = tauVec;
        tauVec = tauVec + Dtau_k1;
        arma::fvec tauVecabs;
        arma::ivec fixrhoidx0, fixrhoidx, tauupdateidx;
        //arma::ivec fixrhoidx0, fixrhoidx, covarianceidxVec1, covarianceidxVec_sub1, covarianceidxVec2, covarianceidxVec_sub2,covarianceidxVec3, covarianceidxVec_sub3, tauupdateidx;
        arma::uvec covarianceidxVec1, covarianceidxVec_sub1, covarianceidxVec2, covarianceidxVec_sub2,covarianceidxVec3, covarianceidxVec_sub3;

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
                tauVecabs = tauVec.elem(covarianceidxVec_sub1) / arma::abs(tauVec.elem(covarianceidxVec_sub1));
                tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
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
                //covarianceidxVec_sub1 = g_covarianceidxMat_col1.elem(arma::find(fixrhoidx == 1));
                //covarianceidxVec_sub2 = g_covarianceidxMat_col2.elem(arma::find(fixrhoidx == 1));
                //covarianceidxVec_sub3 = g_covarianceidxMat_col3.elem(arma::find(fixrhoidx == 1));
                tauVecabs = tauVec.elem(covarianceidxVec_sub1) / arma::abs(tauVec.elem(covarianceidxVec_sub1));
                tauVec.elem(covarianceidxVec_sub1) = tauVecabs % (arma::sqrt(tauVec.elem(covarianceidxVec_sub2) % tauVec.elem(covarianceidxVec_sub3)));
        }
         return List::create(Named("tau") = tauVec, Named("AI") = AI1, Named("score") = score1);
}



// [[Rcpp::export]]
arma::fvec getMeanDiagofKmat(bool LOCO){
        arma::fvec mean_diag_kins_vec(g_num_Kmat - 1);


        arma::sp_vec diagVecG0;
        arma::sp_fvec diagVecV0;
        arma::fvec diagVecG, diagVecV, diagVecG_I, diagVecG_T, diagVecG_IT,diagVecV_I, diagVecV_T, diagVecV_IT;
        arma::fvec diagVec;

                unsigned int tauind = 0;
        if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){


           if(!LOCO){
              if(g_spGRM.n_rows == 0){
		        int MminMAF = geno.getnumberofMarkerswithMAFge_minMAFtoConstructGRM();      
                diagVec = (*geno.Get_Diagof_StdGeno()) /MminMAF;
              }else{
                //arma::vec diagVecofKmat_0;
                //diagVecofKmat_0 = arma::diagvec(g_spGRM);
                //diagVec = arma::conv_to< arma::fvec >::from(diagVecofKmat_0);
		diagVec = arma::diagvec(g_spGRM);
              }
            }else{
              diagVec = (*geno.Get_Diagof_StdGeno_LOCO());
              int Msub_MAFge_minMAFtoConstructGRM_in_b = geno.getMsub_MAFge_minMAFtoConstructGRM_in();
              int Msub_MAFge_minMAFtoConstructGRM_singleVar_b = geno.getMsub_MAFge_minMAFtoConstructGRM_singleChr_in();
              diagVec = diagVec/(Msub_MAFge_minMAFtoConstructGRM_in_b - Msub_MAFge_minMAFtoConstructGRM_singleVar_b);
            }
            mean_diag_kins_vec(0) = arma::mean(diagVec);

          if(Kmat_vec.size() > 0){
            for(unsigned int i = 0; i < Kmat_vec.size(); i++){
              diagVec =(Kmat_vec[i]).diag();
              mean_diag_kins_vec(i+1) = arma::mean(diagVec);
            }
          }

        }else{
                //diagVecG0 = g_spGRM.diag();
                //arma::vec diagVecGtemp(diagVecG0);
                //diagVecG = arma::conv_to< arma::fvec >::from(diagVecGtemp);
		diagVec = arma::diagvec(g_spGRM);
                diagVecG_I = diagVecG.elem(g_I_longl_vec);
		diagVec = diagVecG_I;
                mean_diag_kins_vec(0) = arma::mean(diagVec);
                tauind = tauind + 1;


                if(g_T_longl_mat.n_rows > 0){
                  diagVecG_IT = diagVecG_I % g_T_longl_vec;
                  diagVecG_T = diagVecG_IT % g_T_longl_vec;
                  diagVecG_IT = 2 * diagVecG_IT;
                  diagVec = diagVecG_IT;

		  std::cout << "Here1" << std::endl; 
                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  tauind = tauind + 1;
		  std::cout << "Here2" << std::endl; 
                  diagVec = diagVecG_T;
                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  tauind = tauind + 1;
		  
		  diagVecV = diagVecG;
		  diagVecV.ones();
		  diagVecV_I = diagVecV.elem(g_I_longl_vec);
                  diagVec = diagVecV_I;
                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  tauind = tauind + 1; 	
	          diagVecV_IT = diagVecV_I % g_T_longl_vec;
                  diagVecV_T = diagVecV_IT % g_T_longl_vec;
                  diagVecV_IT = 2 * diagVecV_IT;
                  diagVec = diagVecV_IT;

                  mean_diag_kins_vec(tauind) = arma::mean(diagVec);
                  tauind = tauind + 1;
                  std::cout << "Here2" << std::endl;
                  diagVec = diagVecV_T;
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
