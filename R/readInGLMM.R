getChromNumber = function(chrom = ""){
  if(chrom == ""){
    stop("chrom is not specified\n")
  }else{	  
    chrom_v2 = as.character(chrom)
    chrom_v2 = gsub("CHR", "", chrom_v2, ignore.case = T)
    chrom_v3 = as.numeric(gsub("[^0-9.]", "",chrom_v2))

   if(!is.na(chrom_v3)){

   
    if(chrom_v3 <= 22 &  chrom_v3 >= 1) {
      #stop("chromosome ", chrom, " is out of the range of null model LOCO results\n")
    #}else {
      cat("Leave chromosome ", chrom_v3, " out will be applied\n")
    }else{
      chrom_v3 = chrom	
    }
   }else{
	chrom_v3 = chrom
   }
  }
  
  return(chrom_v3)
}	

removeLOCOResult = function(chromList, obj.glmm.null){
  for (chr in 1:22) {
    if(chr %in% chromList){    
      obj.glmm.null$LOCOResult[chr] = list(NULL)
      cat("chromosome ", chr, " model results are removed to save memory\n")
      gc()
    }
  }
  return(obj.glmm.null)
}	


ReadModel = function(GMMATmodelFile = "", chrom="", LOCO=TRUE, is_Firth_beta=FALSE, is_EmpSPA=FALSE, espa_nt=9999, espa_range=c(-20,20)){
  # Check file existence
  Check_File_Exist(GMMATmodelFile, "GMMATmodelFile")
  if(!LOCO %in% c(TRUE, FALSE))
    stop("LOCO should be TRUE or FALSE.")
  # load GMMATmodelFile
  load(GMMATmodelFile)
  #if(!modglmm$LOCO){
  #     if(LOCO){
  #             stop("LOCO is TRUE but the null model was not fit with LOCO=TRUE\n")
  #             LOCO=FALSE
  #     }
  #}

  modglmm$Y = NULL
  #modglmm$offset = modglmm$linear.predictors - modglmm$coefficients[1]
  modglmm$linear.predictors = NULL
  modglmm$coefficients = NULL
  modglmm$cov = NULL
  #rm(modglmm)
  gc()
  #traitType = modglmm$traitType
  #y = modglmm$y
  #X = modglmm$X
  #N = length(y)
  #tauVec = modglmm$theta
  #indChromCheck = FALSE
  chrom_v3=NULL

  if(!LOCO) {
    print("Leave-one-chromosome-out is not applied")
    if(modglmm$LOCO) {
      for (chr in 1:22) {
        modglmm$LOCOResult[chr] = list(NULL)
        cat("chromosome ", chr, " model results are removed to save memory\n")
        gc()
      }
    }
  }else{
    if (!modglmm$LOCO){
      stop("LOCO is TRUE but the null model file .rda does not contain LOCO results. In order to apply Leave-one-chromosome-out, please run Step 1 using LOCO. Otherwise, please set LOCO=FALSE in this step (Step 2).\n")
    }else{
        if(chrom == ""){
          stop("chrom needs to be specified in order to apply Leave-one-chromosome-out on gene- or region-based tests")
        }else{
          chrom_v3 = getChromNumber(chrom)
        }
   }

 if(is.numeric(chrom_v3)){
  if(chrom_v3 >= 1 & chrom_v3 <= 22){

   chromList = c(1:22)
   chromList = chromList[which(chromList != chrom_v3)]
   modglmm = removeLOCOResult(chromList, modglmm)
   if(!is.null(modglmm$LOCOResult[[chrom_v3]])){
   modglmm$fitted.values = modglmm$LOCOResult[[chrom_v3]]$fitted.values
   modglmm$residuals = modglmm$LOCOResult[[chrom_v3]]$residuals
   modglmm$obj.noK = modglmm$LOCOResult[[chrom_v3]]$obj.noK
   if(is_Firth_beta){
     if(!is.null(modglmm$LOCOResult[[chrom_v3]]$offset)){
       modglmm$offset = modglmm$LOCOResult[[chrom_v3]]$offset
     }
   }
  }
   #modglmm$offset = modglmm$LOCOResult[[chrom_v3]]$linear.predictors -  modglmm$LOCOResult[[chrom_v3]]$coefficients[1]
   modglmm$LOCOResult[chrom_v3] = list(NULL)
  }else{ #if(chrom_v3 >= 1 & chrom_v3 <= 22){
     chromList = c(1:22)
     modglmm = removeLOCOResult(chromList, modglmm)
  }
 }else{
        chromList = c(1:22)
        modglmm = removeLOCOResult(chromList, modglmm)
  }
   gc()

  }

  modglmm$mu = as.vector(modglmm$fitted.values)
  tau = modglmm$theta
  N = length(modglmm$mu)
  if(modglmm$traitType == "binary"){
    modglmm$mu2 = (modglmm$mu) *(1-modglmm$mu)
    modglmm$obj_cc = SKAT::SKAT_Null_Model(modglmm$y ~ modglmm$X-1, out_type="D", Adjustment = FALSE)
    modglmm$obj_cc$mu = modglmm$mu
    modglmm$obj_cc$res = modglmm$res
    modglmm$obj_cc$pi_1 = modglmm$mu2
  }else if(modglmm$traitType == "quantitative"){
    modglmm$mu2 = (1/tau[1])*rep(1,N)
  }else if(modglmm$traitType == "count"){
    modglmm$mu2 = modglmm$mu
  }else if(modglmm$traitType == "count_nb"){
    modglmm$mu2 = (1/tau[1])*rep(1,N)
  }	  
 #if(FALSE){

  modglmm$obj_cc = list(res.out = NULL)
  print("check")
  if(!is.null(modglmm$obj_cc$res.out)){
               modglmm$obj_cc$res.out<-cbind(modglmm$obj_cc$res,
modglmm$obj_cc$res.out)
            } else {
                modglmm$obj_cc$res.out<-modglmm$res
             }



 if(is_Firth_beta){
        if(modglmm$traitType == "binary"){
                if(is.null(modglmm$offset)){
                        #if(FALSE){
                        #print("WARNING. is_Firth_beta = TRUE, but offset was not computed in Step 1. Please re-run Step 1 using the more rencet version of SAIGE/SAIGE-GENE.")
                        cat("Applying is_Firth_beta = TRUE and note the results would be more accurate withe Step 1 results using the more rencet version of SAIGE/SAIGE-GENE. \n")

                        if(ncol(modglmm$X) == 1){
                                covoffset = rep(0,nrow(modglmm$X))
                        }else{
                                formulastr = paste0("y ~ ", paste(colnames(modglmm$X)[-1], collapse="+"))
                                modglmm$X = cbind(modglmm$X, as.vector(modglmm$y))
                                colnames(modglmm$X)[ncol(modglmm$X)] = "y"
                                modglmm$X = as.data.frame(modglmm$X)
                                formula.new = as.formula(formulastr)
                                modwitcov = glm(formula.new, data = modglmm$X, family = binomial)
                                modglmm$X = modglmm$X[,-ncol(modglmm$X)]
                                covoffset = as.matrix(modglmm$X[,-1]) %*%  modwitcov$coefficients[-1]
                                modglmm$X = as.matrix(modglmm$X)
                        }
                        modglmm$offset = covoffset
                        #}
                }
        }
  }


 if(is.null(modglmm$offset)){
         modglmm$offset = rep(0,nrow(modglmm$X))
  }


 if(is.null(modglmm$Sigma_iXXSigma_iX)){
        modglmm$Sigma_iXXSigma_iX = matrix(0, nrow=1, ncol=1)
 }

 if(is.null(modglmm$varWeights)){
	modglmm$varWeights = rep(1, length(modglmm$y))
 }	 


 if(is_EmpSPA){
    idx0 = qcauchy(1:espa_nt/(espa_nt+1))
    idx1 = idx0*max(espa_range)/max(idx0)
    resid = modglmm$residuals
    var_resid = var(resid)
    cat("var_resid ", var_resid, "\n")
    scaled_resid= as.numeric(resid) / as.numeric(sqrt(var_resid))

    cumul<-NULL
    for(id in idx1){
      t<-id
      e_resid<-exp(scaled_resid*t)
      M0<-mean(e_resid)
      M1<-mean(scaled_resid*e_resid)
      M2<-mean(scaled_resid^2*e_resid)
      k0<-log(M0)
      k1<-M1/M0
      k2<-(M0*M2-M1^2)/M0^2
      cumul<-rbind(cumul, c(t, k0, k1, k2))
    }
    
    #K_org_emp<-approxfun(cumul[,1], cumul[,2], rule=2)
    #K_1_emp<-approxfun(cumul[,1], cumul[,3], rule=2)
    #K_2_emp<-approxfun(cumul[,1], cumul[,4], rule=2)    
    #modglmm$K_org_emp=K_org_emp
    #modglmm$K_1_emp=K_1_emp
    #modglmm$K_2_emp=K_2_emp
  }else{
    cumul = matrix(1:8, ncol=4)
  }
  modglmm$cumul = cumul

 return(modglmm)
}




Get_Variance_Ratio<-function(varianceRatioFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest, isSparseGRM, useSparseGRMtoFitNULL){

    iscateVR = FALSE
    # check variance ratio
    if (!file.exists(varianceRatioFile)) {
	if(varianceRatioFile != ""){
	    stop("varianceRatioFile is specified but the file ", varianceRatioFile, " does not exist\n")
	}else{	
            cat("varianceRatioFile is not specified so variance ratio won't be used\n")
	    
	}
	if(isSparseGRM){
          cat("WARNING: Sparse GRM is specified. Please make sure the null model was fit using the same sparse GRM in Step 1.\n")
          ratioVec_sparse = c(1)
	  ratioVec_null = c(-1)
	}else{
          ratioVec_sparse = c(-1)		
	  ratioVec_null = c(1)
	}
    }else{
        varRatioData = data.frame(data.table:::fread(varianceRatioFile, header = F, stringsAsFactors = FALSE))
	if(ncol(varRatioData) == 3){
	    spindex = which(varRatioData[,2] == "sparse")
	    if(length(spindex) > 0){
	        ratioVec_sparse = varRatioData[which(varRatioData[,2] == "sparse"),1]
		ratioVec_sparse = as.numeric(ratioVec_sparse)
		if(!isSparseGRM & sum(ratioVec_sparse > 1.0001 | ratioVec_sparse < 0.9999) > 0){
		       	stop("sparse GRM is not specified but it was used for estimating variance ratios in Step 1. Please specify --sparseGRMFile and --sparseGRMSampleIDFile\n")
		}	
	    }else{    
		ratioVec_sparse = c(-1)
	    	if(isSparseGRM){
			stop("sparse GRM is specified but the variance ratio for sparse GRM was not estimatedin Step 1. Pleae remove --sparseGRMFile and --sparseGRMSampleIDFile\n")
		}	
	    }
	    ratioVec_null = varRatioData[which(varRatioData[,2] == "null"),1]
	    ratioVec_null = as.numeric(ratioVec_null)
            cat("variance Ratio null is ", ratioVec_null, "\n")
	    if(length(ratioVec_null) > 1){
		iscateVR = TRUE
		nrv = length(ratioVec_null)
	    }else{
		nrv = 1
	    }
	}else{
	    cat("Variance ratios were estimated with version < 1.0.6\n")	
	    if(isSparseGRM){
	        ratioVec_sparse = varRatioData[,1]
	        ratioVec_sparse = as.numeric(ratioVec_sparse)
        	cat("variance Ratio is ", ratioVec_sparse, "\n")
		ratioVec_null = rep(-1, length(varRatioData[,1]))
		cat("WARNING: Sparse GRM is specified and the variance ratio(s) were specified. Please make sure the variance ratios were estimated using a full GRM and a sparse GRM.")
	    }else{
	        ratioVec_null = varRatioData[,1]
	        ratioVec_null = as.numeric(ratioVec_null)
                cat("variance Ratio is ", ratioVec_null, "\n")
		cat("WARNING: Sparse GRM is not specified and the variance ratio(s) were specified. Please make sure that in Step 1, 1. the null model was fit using a full GRM (--useSparseGRMtoFitNULL=FALSE) and the variacne ratio was NOT estimated with the sparse GRM (useSparseGRMforVarRatio=FALSE) or 2. the null model was fit using a sparse GRM (--useSparseGRMtoFitNULL=TRUE) and the variacne ratio was estiamted with the sparse GRM and null --skipVarianceRatioEstimation=FALSE\n")
		ratioVec_sparse = c(-1)
	    }
	    if(length(varRatioData[,1]) > 1){
		iscateVR = TRUE
		nrv = length(varRatioData[,1])
	    }else{
		nrv = 1
	    }
	}


	if(iscateVR){
            ln = length(cateVarRatioMinMACVecExclude)
            hn = length(cateVarRatioMaxMACVecInclude)

            if (nrv != ln) {
                stop("ERROR! The number of variance ratios are different from the length of cateVarRatioMinMACVecExclude\n")
            }
            if (ln != (hn + 1)) {
                stop("ERROR! The length of cateVarRatioMaxMACVecInclude does not match with the lenght of cateVarRatioMinMACVecExclude (-1)\n")
            }
        }
	

    }
    return(list(ratioVec_sparse = ratioVec_sparse, ratioVec_null = ratioVec_null))
}


