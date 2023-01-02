
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
