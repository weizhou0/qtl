


bool isUseSparseSigmaforModelFitting = false;
std::vector<arma::sp_fmat> Kmat_vec;
arma::sp_fmat g_I_longl_mat;
arma::sp_fmat g_T_longl_mat;
arma::uvec g_I_longl_vec;
arma::fvec g_T_longl_vec;
arma::sp_fmat g_spGRM;
arma::sp_fmat g_spSigma;
bool g_isSparseGRM;
bool g_isStoreSigma;
int g_num_Kmat;
bool g_isGRM;

// [[Rcpp::export]]
void setupSparseGRM_new(arma::sp_mat & t_spGRM){
        arma::sp_mat t_spGRM_1 = t_spGRM;
        arma::sp_fmat g_spGRM_f = arma::conv_to<arma::sp_fmat>::from(t_spGRM_1);
        g_spGRM = g_spGRM_f;
}

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
arma::fvec getCrossprodMatAndKin(arma::fcolvec& bVec, bool LOCO){
    arma::fvec crossProdVec;
    if(g_isSparseGRM){
        crossProdVec = g_spGRM * bVec;
    }else{
	if(!LOCO){    
          crossProdVec = parallelCrossProd(bVec);
	}else{
	  crossProdVec = parallelCrossProd_LOCO(bVec)
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
		if(g_isGRM){
                  Ibvec = g_I_longl_mat.t() * bVec;
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
              int MminMAF = geno->getnumberofMarkerswithMAFge_minMAFtoConstructGRM();
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
                  diagVecG_IT = diagVecG_I % g_T_longl_vec;
                  diagVecG_T = diagVecG_IT % g_T_longl_vec;
                  diagVecG_IT = 2 * diagVecG_IT;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_IT;
                  tauind = tauind + 1;
                  diagVec = diagVec + tauVec(tauind) * diagVecG_T;
                  tauind = tauind + 1;
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
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
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
        minvVec = 1/getDiagOfSigma_multiV(wVec, tauVec, LOCO);
        zVec = minvVec % rVec;

        float sumr2 = sum(rVec % rVec);
        arma::fvec z1Vec(Nnomissing);
        arma::fvec pVec = zVec;

        int iter = 0;
        while (sumr2 > tolPCG && iter < maxiterPCG) {
		iter = iter + 1;
                arma::fcolvec ApVec = getCrossprod_multiV(pVec, wVec, tauVec, LOCO);
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
                if(ptr_gNULLGENOobj->alleleFreqVec[i] >= minMAFtoConstructGRM && ptr_gNULLGENOobj->alleleFreqVec[i] <= 1-minMAFtoConstructGRM){

                        ptr_gNULLGENOobj->Msub = ptr_gNULLGENOobj->Msub + 1;
                }
        }
  }
}

// [[Rcpp::export]]
void setStartEndIndexVec( arma::ivec & startIndex_vec,  arma::ivec & endIndex_vec){
  ptr_gNULLGENOobj->startIndexVec = startIndex_vec;
  ptr_gNULLGENOobj->endIndexVec = endIndex_vec;
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
	   tauind = tauind + 1
	  }

          if(Kmat_vec.size() > 0){
            for(unsigned int i = 0; i < Kmat_vec.size(); i++){
              diagVec = (Kmat_vec[i]).diag();
              mean_diag_kins_vec(i+tauind) = arma::mean(diagVec);
            }
          }

        }else{ //if(g_I_longl_mat.n_rows == 0 && g_T_longl_mat.n_rows == 0){
	    if(g_isGRM && g_isSparseGRM){	
                diagVecG = arma::diagvec(g_spGRM);
                diagVecG_I = diagVecG.elem(g_I_longl_vec);
                diagVec = diagVecG_I;
                mean_diag_kins_vec(0) = arma::mean(diagVec);
                tauind = tauind + 1;
            } 

            if(g_T_longl_mat.n_rows > 0){
	        if(g_isGRM && g_isSparseGRM){		
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
		}  

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

        std::cout << "getAIScore_multiV " << getAIScore_multiV << std::endl;
        int q2 = arma::sum(fixtauVec==0);
        arma::uvec idxtau = arma::find(fixtauVec==0);
        arma::fvec tau0;
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
                        if(g_isGRM){
                                GRM_I_bvec = getCrossprodMatAndKin(Ibvec, LOCO);
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
                                        Au_mat.col(2) = temp_vec_double;
                                        temp_mat(i,2) = dot(temp_vec_double, Pu);
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

        int Nnomissing = ptr_gNULLGENOobj->getNnomissing();
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
  minMAFtoConstructGRM = minMAFforGRM;
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
        int Ntotal = ptr_gNULLGENOobj->getNnomissing();
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
//      std::cout << "(ptr_gNULLGENOobj->subMarkerIndex).n_elem: " << (ptr_gNULLGENOobj->subMarkerIndex).n_elem << std::endl;
        int Nnomissing = ptr_gNULLGENOobj->getNnomissing();
        (ptr_gNULLGENOobj->stdGenoMultiMarkersMat).set_size(subMarkerIndexRandom.n_elem, Nnomissing);
}

// [[Rcpp::export]]
void setRelatednessCutoff(float a){
        ptr_gNULLGENOobj->relatednessCutoff = a;
}


// [[Rcpp::export]]
double innerProduct(NumericVector x, NumericVector y) {
   return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

