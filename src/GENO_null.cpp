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
#include "GENO_null.hpp"


using namespace Rcpp;
using namespace std;
using namespace RcppParallel;

namespace NullGENO {

void NullGenoClass::setStdGenoLookUpArr(float mafVal, float invsdVal, arma::fvec & stdGenoLookUpArr){
                float mafVal2 = 2*mafVal;
                stdGenoLookUpArr(0) = (0-mafVal2)*invsdVal;
                stdGenoLookUpArr(1) = (1-mafVal2)*invsdVal;
                stdGenoLookUpArr(2) = (2-mafVal2)*invsdVal;
}



void NullGenoClass::setSparseKinLookUpArr(float mafVal, float invsdVal){
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

void NullGenoClass::setBit(unsigned char & ch, int ii, int aVal, int bVal){

                if (bVal == 1 && aVal == 1){
                        ch ^= char(1 << ((ii*2) + 1)); //set a to be 1

                }else if(bVal == 0){
                        ch ^= char(1 << (ii*2)); //change b to 0

                        if(aVal == 1){
                                ch ^= char(1 << ((ii*2) + 1)); //change a to 1
                        }
                }
}


void NullGenoClass::setGenotype(unsigned char* c, const int pos, const int geno) {
                (*c) |= (geno << (pos << 1));
}


void NullGenoClass::getGenotype(unsigned char* c, const int pos, int& geno) {
                geno = ((*c) >> (pos << 1)) & 0x3;  // 0b11 = 0x3
}

void NullGenoClass::Init_OneSNP_Geno(){
                m_size_of_esi = (Nnomissing+3)/4;
                int k = 8;
                while (k > 0){
                        -- k;
                        m_bits_val[k] = 1 << k;
                }
}

arma::ivec* NullGenoClass::Get_OneSNP_Geno(size_t SNPIdx){
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



arma::ivec* NullGenoClass::Get_OneSNP_Geno_forVarRatio(size_t SNPIdx){
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



void NullGenoClass::Get_OneSNP_Geno_atBeginning(size_t SNPIdx, vector<int> & indexNA, vector<unsigned char> & genoVecOneMarkerOld, float & altFreq, float & missingRate, int & mac,  int & alleleCount, bool & passQC, size_t SNPIdx_new, bool & passVarRatio , size_t SNPIdx_vr){
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
        //                      std::cout << "bufferGeno " << bufferGeno << std::endl;
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
                                //      std::cout << "[ind] " << ind << std::endl;
                                //      std::cout << "indicatorGenoSamplesWithPheno_in[ind] " << indicatorGenoSamplesWithPheno_in[ind] << std::endl;
                                //}
                                if(indicatorGenoSamplesWithPheno_in[ind]){
                                        if(bufferGeno == 3){
        //                                      std::cout << "SNPIdx " << SNPIdx << std::endl;
        //                                      std::cout << "ind " << ind << std::endl;
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
                        //      u = indx & 3;
                        //      bufferGeno = m_OneSNP_GenoTemp[ptrsubSampleInGeno[indx] - 1];
                        alleleCount = alleleCount + fillinMissingGeno*numMissing;

                        //              if(bufferGeno == 3){
                        //                      bufferGeno = fillinMissingGeno;
                        //                      alleleCount = alleleCount + bufferGeno;
                        //              }
                        //      }
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

int NullGenoClass::Get_OneSNP_StdGeno(size_t SNPIdx, arma::fvec * out ){
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
//              std::cout << "freq " << freq << endl;
//              std::cout << "invStd " << invStd << endl;

                //setStdGenoLookUpArr(freq, invStd);
                //std::cout << "stdGenoLookUpArr[0]: " << stdGenoLookUpArr[0] << std::endl;
                //std::cout << "stdGenoLookUpArr[1]: " << stdGenoLookUpArr[1] << std::endl;
                //std::cout << "stdGenoLookUpArr[2]: " << stdGenoLookUpArr[2] << std::endl;
//              cout << "Get_OneSNP_StdGeno here2"  << endl;
                for(size_t i=Start_idx; i< Start_idx+m_size_of_esi-1; i++){
//                      geno1 = genoVec[i];
                        geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);

                        for(int j=0; j<4; j++){
                        int b = geno1 & 1 ;
                        geno1 = geno1 >> 1;
                        int a = geno1 & 1 ;
                        //(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;;
                        (*out)[ind] = stdGenoLookUpArr(2-(a+b));
//                      std::cout << "a " << a << endl;
//                      std::cout << "b " << b << endl;
//                      std::cout << "(*out)[ind] " << (*out)[ind] << endl;
                        ind++;
                        geno1 = geno1 >> 1;

//                      if(ind >= Nnomissing){
//                              cout << "Get_OneSNP_StdGeno " << SNPIdx << endl;
//                              cout << "Nnomissing " << Nnomissing << endl;
//                              stdGenoLookUpArr.clear();
//                              return 1;
//                      }
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


arma::fvec* NullGenoClass::Get_Diagof_StdGeno(){

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
                //              std::cout << "i " << i << std::endl;
                //              std::cout << "numberofMarkerswithMAFge_minMAFtoConstructGRM " << numberofMarkerswithMAFge_minMAFtoConstructGRM << std::endl;
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


arma::fvec* NullGenoClass::Get_Diagof_StdGeno_LOCO(){
                //if(size(m_DiagStd_LOCO)[0] != Nnomissing){
                //m_DiagStd_LOCO.zeros(Nnomissing);
                  //      for(size_t i=startIndex; i<= endIndex; i++){
                                //if(i < startIndex || i > endIndex){
                //                      if(alleleFreqVec[i] >= minMAFtoConstructGRM && alleleFreqVec[i] <= 1-minMAFtoConstructGRM){
                  //                            Get_OneSNP_StdGeno(i, temp);
                    //                          m_DiagStd_LOCO = m_DiagStd_LOCO + (*temp) % (*temp);
                //                              Msub_MAFge_minMAFtoConstructGRM = Msub_MAFge_minMAFtoConstructGRM + 1;
                //                      }
                                //}
                 //       }


                //m_DiagStd_LOCO = m_DiagStd - geno.mtx_DiagStd_LOCO.col(geno.chromIndex);
                m_DiagStd_LOCO = mtx_DiagStd_LOCO.col(chromIndex);
                Msub_MAFge_minMAFtoConstructGRM_singleChr  = Msub_MAFge_minMAFtoConstructGRM_byChr(chromIndex);
                //}

                return & m_DiagStd_LOCO;
}


void NullGenoClass::setGenoObj(std::string bedfile, std::string bimfile, std::string famfile, std::vector<int> & subSampleInGeno, std::vector<bool> & indicatorGenoSamplesWithPheno, float memoryChunk, bool  isDiagofKinSetAsOne){

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
        //              cout << "i is " << i << endl;

                        indexNA.clear();
                //}
                        Get_OneSNP_Geno_atBeginning(i, indexNA, genoVecOneMarkerOld, freq, missingRate, mac, alleleCount, isPassQC, SNPIdx_new, isPass_vr, SNPIdx_vr);

                        //std::cout << "freq " << freq << std::endl;
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
                //      cout << M << " markers with MAF >= " << minMAFtoConstructGRM << endl;
                //}


                int numofGenoArray_old = numofGenoArray;
                if(numberofMarkerswithMAFge_minMAFtoConstructGRM % numMarkersofEachArray == 0){
                        numofGenoArray = numberofMarkerswithMAFge_minMAFtoConstructGRM / numMarkersofEachArray;
                        //genoVecofPointers.resize(numofGenoArray);
                        //cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;

                }else{
                        numofGenoArray = numberofMarkerswithMAFge_minMAFtoConstructGRM/numMarkersofEachArray + 1;
                }
//              cout << " numofGenoArray "<< numofGenoArray << endl;
//              cout << " numofGenoArray_old "<< numofGenoArray_old << endl;
//              cout << " genoVecofPointers.size() "<< genoVecofPointers.size() << endl;
                if(numofGenoArray > numofGenoArray_old){

                        for (int i = numofGenoArray; i < numofGenoArray_old ; i++){
                                delete genoVecofPointers[i];
                                //genoVecofPointers[i] = new vector<unsigned char>;
                                //genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(Nnomissing)/4));
                        }
                }
//              cout << " genoVecofPointers.size() "<< genoVecofPointers.size() << endl;

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
//              cout << "setgeno mark5" << endl;
//              printAlleleFreqVec();
                //printGenoVec();
                //Get_Diagof_StdGeno();
//              cout << "setgeno mark6" << endl;
}//End Function



void NullGenoClass::printFromgenoVec(unsigned char genoBinary0){
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


int NullGenoClass::getM() const{
                return(M);
}

int NullGenoClass::getnumberofMarkerswithMAFge_minMAFtoConstructGRM() const{
                return(numberofMarkerswithMAFge_minMAFtoConstructGRM);
}


int NullGenoClass::getMsub() const{
                return(Msub);
}

int NullGenoClass::getStartIndex() const{
                return(startIndex);
}

int NullGenoClass::getEndIndex() const{
                return(endIndex);
}

int NullGenoClass::getN() const{
                return(N);
}

int NullGenoClass::getNnomissing() const{
                return(Nnomissing);
}

float NullGenoClass::getAC(int m){
                return(alleleFreqVec[m]*2*Nnomissing);
}

float NullGenoClass::getMAC(int m){
                if(alleleFreqVec[m] > 0.5){
                        return((1-alleleFreqVec[m])*2*Nnomissing);
                }else{
                        return(alleleFreqVec[m]*2*Nnomissing);
                }
}

int NullGenoClass::getMsub_MAFge_minMAFtoConstructGRM_in() const{
                return(numberofMarkerswithMAFge_minMAFtoConstructGRM);
}

int NullGenoClass::getMsub_MAFge_minMAFtoConstructGRM_singleChr_in() const{
                return(Msub_MAFge_minMAFtoConstructGRM_singleChr);
}

void NullGenoClass::Get_Samples_StdGeno(arma::ivec SampleIdsVec){
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



}	
