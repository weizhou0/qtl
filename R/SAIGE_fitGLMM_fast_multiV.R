#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants 
#'
#' @param plinkFile character. Path to plink file to be used for calculating elements of the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF of markers in the plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios
#' @param phenoFile character. Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null model. Please specify the names of the covariates using the argument covarColList and specify categorical covariates using the argument qCovarCol. All categorical covariates must also be included in covarColList.
#' @param phenoCol character. Column name for the phenotype in phenoFile e.g. "CAD"
#' @param traitType character. e.g. "binary" or "quantitative". By default, "binary"
#' @param invNormalize logical. Whether to perform the inverse normalization for the phentoype or not. e.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the null model. e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the null model. All categorical covariates listed in qCovarCol must be also in covarColList,  e,g c("Sex"). 
#' @param sampleIDColinphenoFile character. Column name for the sample IDs in the phenotype file e.g. "IID".  
#' @param tol numeric. The tolerance for fitting the null model to converge. By default, 0.02.
#' @param maxiter integer. The maximum number of iterations used to fit the null GLMMM. By default, 20.
#' @param tolPCG numeric. The tolerance for PCG to converge. By default, 1e-5.
#' @param maxiterPCG integer. The maximum number of iterations for PCG. By default, 500. 
#' @param nThreads integer. Number of threads to be used. By default, 1 
#' @param SPAcutoff numeric. The cutoff for the deviation of score test statistics from the mean in the unit of sd to perform SPA. By default, 2.
#' @param numMarkersForVarRatio integer (>0). Minimum number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed 
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param tauInit vector of numbers. e.g. c(1,1), Initial values for tau. For binary traits, the first element will be always be set to 1. If the tauInit is 0,0, the second element will be 0.5 for binary traits and the initial tau vector for quantitative traits is 1,0 
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out (LOCO) option. By default, TRUE
#' @param traceCVcutoff numeric. The threshold for coefficient of variation (CV) for the trace estimator to increase nrun. By default, 0.0025
#' @param ratioCVcutoff numeric. The threshold for coefficient of variation (CV) for the variance ratio estimate. If ratioCV > ratioCVcutoff. numMarkersForVarRatio will be increased by 10. By default, 0.001 
#' @param outputPrefix character. Path to the output files with prefix.
#' @param outputPrefix_varRatio character. Path to the output variance ratio file with prefix. variace ratios will be output to outputPrefix_varRatio.varianceRatio.txt. If outputPrefix_varRatio is not specified, outputPrefix_varRatio will be the same as the outputPrefix. By default, ""
#' @param IsOverwriteVarianceRatioFile logical. Whether to overwrite the variance ratio file if the file exists. By default, FALSE
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file. By default, "" 
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM. By default, "" 
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples that are included in the sparse GRM. By default, 2000
#' @param relatednessCutoff float. The threshold for coefficient of relatedness to treat two samples as unrelated in the sparse GRM. 
#' @param cateVarRatioIndexVec vector of integer 0 or 1. The length of cateVarRatioIndexVec is the number of MAC categories for variance ratio estimation. 1 indicates variance ratio in the MAC category is to be estimated, otherwise 0. By default, NULL. If NULL, variance ratios corresponding to all specified MAC categories will be estimated. This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(10.5,20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param isCovariateTransform logical. Whether use qr transformation on non-genetic covariates. By default, TRUE
#' @param isDiagofKinSetAsOne logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param useSparseGRMtoFitNULL logical. Whether to use sparse GRM to fit the null model. By default, FALSE
#' @param useSparseGRMforVarRatio logical. Whether to use sparse GRM to estimate the variance Ratios. If TRUE, the variance ratios will be estimated using the full GRM (numerator) and the sparse GRM (denominator). By default, FALSE
#' @param minCovariateCount integer. If binary covariates have a count less than this, they will be excluded from the model to avoid convergence issues. By default, -1 (no covariates will be excluded)
#' @param minMAFforGRM numeric. Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01
#' @param includeNonautoMarkersforVarRatio logical. Whether to allow for non-autosomal markers for variance ratio. By default, FALSE
#' @param FemaleOnly logical. Whether to run Step 1 for females only. If TRUE, sexCol and FemaleCode need to be specified. By default, FALSE
#' @param MaleOnly logical. Whether to run Step 1 for males only. If TRUE, sexCol and MaleCode need to be specified. By default, FALSE
#' @param FemaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for females. By default, '1' 
#' @param MaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for males. By default, '0'
#' @param sexCol character. Coloumn name for sex in the phenotype file, e.g Sex. By default, '' 
#' @param isCovariateOffset logical. Whether to estimate fixed effect coeffciets. By default, FALSE.  
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @export
fitNULLGLMM_multiV = function(plinkFile = "",
		bedFile="",
		bimFile="",
		famFile="",
                phenoFile = "",
                phenoCol = "",
                traitType = "binary",
                invNormalize = FALSE,
                covarColList = NULL,
                qCovarCol = NULL,
		longlCol = "",
                sampleIDColinphenoFile = "",
                tol=0.02,
                maxiter=20,
                tolPCG=1e-5,
                maxiterPCG=500,
                nThreads = 1, 
                SPAcutoff = 2, 
                numMarkersForVarRatio = 30, 
                skipModelFitting = FALSE,
		memoryChunk = 2,
		tauInit = c(0,0),
		LOCO = TRUE,
		isLowMemLOCO = FALSE,
		traceCVcutoff = 0.0025,
		ratioCVcutoff = 0.001, 
                outputPrefix = "",
		outputPrefix_varRatio = "",
		IsOverwriteVarianceRatioFile=FALSE,
		sparseGRMFile="",
                sparseGRMSampleIDFile="",
		numRandomMarkerforSparseKin = 1000,
		relatednessCutoff = 0.125, 
		isCateVarianceRatio = FALSE,
		cateVarRatioIndexVec = NULL,
		cateVarRatioMinMACVecExclude = c(10,20.5),
		cateVarRatioMaxMACVecInclude = c(20.5),
		isCovariateTransform = TRUE,
		isDiagofKinSetAsOne = FALSE,
		minCovariateCount = -1, 
		minMAFforGRM = 0.01,
		maxMissingRateforGRM = 0.15,
		useSparseGRMtoFitNULL=FALSE,
		useSparseGRMforVarRatio=FALSE,
		includeNonautoMarkersforVarRatio = FALSE,
		sexCol = "",
		FemaleCode = 1,
		FemaleOnly = FALSE,
		MaleCode = 0,	
		MaleOnly = FALSE,
		SampleIDIncludeFile = "",
		isCovariateOffset = FALSE,
		skipVarianceRatioEstimation = FALSE, 
		nrun = 30,
		VmatFilelist = "",
		VmatSampleFilelist = "", 
		useGRMtoFitNULL=TRUE)
{

    ##set up output files
    modelOut = paste0(outputPrefix, ".rda")

    if(skipModelFitting){
       if(!file.exists(modelOut)){
          stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
       }
    }else{
       if(LOCO & isLowMemLOCO){	   
	  modelOut = paste(c(outputPrefix,"_noLOCO.rda"), collapse="")	
	}
       file.create(modelOut, showWarnings = TRUE)
    }

    if(plinkFile != ""){
	bimFile = paste0(plinkFile, ".bim")
	bedFile = paste0(plinkFile, ".bed")
	famFile = paste0(plinkFile, ".fam")
    }
    setgenoNULL()

    if(!useGRMtoFitNULL){
        useSparseGRMtoFitNULL = FALSE
        useSparseGRMforVarRatio = FALSE
        LOCO = FALSE
        nThreads = 1
        cat("No GRM will be used to fit the NULL model and nThreads is set to 1\n")
    }


    if (useSparseGRMtoFitNULL & bedFile == ""){
	cat("Sparse GRM is used to fit the null model and plink file is not specified, so variance ratios won't be estimated\n")    
	skipVarianceRatioEstimation = TRUE
    } 

    if(!skipVarianceRatioEstimation){
    	SPAGMMATOut = paste0(outputPrefix, "_", numMarkersForVarRatio, "markers.SAIGE.results.txt")
    	#Check_OutputFile_Create(SPAGMMATOut)

    	if (outputPrefix_varRatio == "") {
            outputPrefix_varRatio = outputPrefix
    	}    
    	varRatioFile = paste0(outputPrefix_varRatio, ".varianceRatio.txt")

    	if (!file.exists(varRatioFile)) {
            file.create(varRatioFile, showWarnings = TRUE)
    	}else {
            if (!IsOverwriteVarianceRatioFile) {
            stop("WARNING: The variance ratio file ", varRatioFile, 
            " already exists. The new variance ratios will be output to ", 
            varRatioFile, ". In order to avoid overwriting the file, please remove the ", 
            varRatioFile, " or use the argument outputPrefix_varRatio to specify a different prefix to output the variance ratio(s). Otherwise, specify --IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten with new variance ratio(s)\n")
            }else{
                cat("The variance ratio file ", varRatioFile, " already exists. IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten\n")
            }
        }
    }else{
	cat("Variance ratio estimation will be skipped\n.")
        useSparseGRMforVarRatio = FALSE
    }


    if (useSparseGRMtoFitNULL){
        #useSparseGRMforVarRatio = FALSE
        LOCO = FALSE
	nThreads = 1
	if(bedFile != ""){
	  cat("sparse GRM will be used to fit the NULL model and nThreads is set to 1\n")
	}
	cat("Leave-one-chromosome-out is not applied\n")
    }



    if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
        if (!file.exists(sparseGRMFile)) {
            stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
        }
        if (!file.exists(sparseGRMSampleIDFile)) {
            stop("sparseGRMSampleIDFile ", sparseGRMSampleIDFile, 
                " does not exist!")
        }
    }


    if (nThreads > 1) {
        RcppParallel:::setThreadOptions(numThreads = nThreads)
        cat(nThreads, " threads will be used ", "\n")
    }

    if (FemaleOnly & MaleOnly) {
        stop("Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job\n")
    }
    
    if (FemaleOnly) {
        outputPrefix = paste0(outputPrefix, "_FemaleOnly")
        cat("Female-specific model will be fitted. Samples coded as ", 
            FemaleCode, " in the column ", sexCol, " in the phenotype file will be included\n")
    }else if (MaleOnly) {
        outputPrefix = paste0(outputPrefix, "_MaleOnly")
        cat("Male-specific model will be fitted. Samples coded as ", 
            MaleCode, " in the column ", sexCol, " in the phenotype file will be included\n")
    }

    if ((!useSparseGRMtoFitNULL & useGRMtoFitNULL)  | !skipVarianceRatioEstimation){
        if (!file.exists(bedFile)) {
            stop("ERROR! bed file does not exsit\n")
        }
        if (!file.exists(bimFile)) {
            stop("ERROR! bim file does not exsit\n")
        }else {
            if (LOCO){
  		chrVec = data.table:::fread(bimFile, header = F, data.table=F , select = 1)
  		updatechrList = updateChrStartEndIndexVec(chrVec)
  		LOCO = updatechrList$LOCO
  		chromosomeStartIndexVec = updatechrList$chromosomeStartIndexVec
  		chromosomeEndIndexVec = updatechrList$chromosomeEndIndexVec
            }	    
	    if(!LOCO) {
                chromosomeStartIndexVec = rep(NA, 22)
                chromosomeEndIndexVec = rep(NA, 22)
            }
        }


        if (!file.exists(famFile)) {
            stop("ERROR! fam file does not exsit\n")
        }else{
            sampleListwithGenov0 = data.table:::fread(famFile, header = F, , colClasses = list(character = 1:4))
            sampleListwithGenov0 = data.frame(sampleListwithGenov0)
            colnames(sampleListwithGenov0) = c("FIDgeno", "IIDgeno", 
                "father", "mother", "sex", "phe")
            sampleListwithGeno = NULL
            sampleListwithGeno$IIDgeno = sampleListwithGenov0$IIDgeno
            sampleListwithGeno = data.frame(sampleListwithGeno)
            sampleListwithGeno$IndexGeno = seq(1, nrow(sampleListwithGeno), 
                by = 1)
            cat(nrow(sampleListwithGeno), " samples have genotypes\n")
        }
    }else{

      if(useSparseGRMtoFitNULL | useSparseGRMforVarRatio){    
        sampleListwithGenov0 = data.table:::fread(sparseGRMSampleIDFile, 
        header = F, , colClasses = c("character"), data.table = F)
        colnames(sampleListwithGenov0) = c("IIDgeno")
        sampleListwithGeno = NULL
        sampleListwithGeno$IIDgeno = sampleListwithGenov0$IIDgeno
        sampleListwithGeno = data.frame(sampleListwithGeno)
        sampleListwithGeno$IndexGeno = seq(1, nrow(sampleListwithGeno), 
            by = 1)
        cat(nrow(sampleListwithGeno), " samples are in the sparse GRM\n")
      }	
    }	    


    if (!file.exists(phenoFile)) {
        stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
    }else{
        if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
            data = data.table:::fread(cmd = paste0("gunzip -c ", 
                phenoFile), header = T, stringsAsFactors = FALSE, 
                colClasses = list(character = sampleIDColinphenoFile), data.table=F)
        }
        else {
            data = data.table:::fread(phenoFile, header = T, 
                stringsAsFactors = FALSE, colClasses = list(character = sampleIDColinphenoFile), data.table=F)
        }
        #data = data.frame(ydat)
	if(longlCol == ""){
		checkColList = c(phenoCol, covarColList, sampleIDColinphenoFile)
	}else{
		checkColList = c(phenoCol, covarColList, sampleIDColinphenoFile, longlCol)
	}

        for (i in checkColList) {
            if (!(i %in% colnames(data))) {
                stop("ERROR! column for ", i, " does not exist in the phenoFile \n")
            }else{
		data = data[,which(colnames(data) %in% checkColList), drop=F]	
		data = data[complete.cases(data),,drop=F]
	    }		    
        }

	if(SampleIDIncludeFile != ""){
		if(!file.exists(SampleIDIncludeFile)){
			stop("ERROR! SampleIDIncludeFile ", SampleIDIncludeFile, " does not exsit\n")
		}else{
			sampleIDInclude = data.table:::fread( SampleIDIncludeFile, header = F,  stringsAsFactors = FALSE, colClasses = c("character"), data.table=F)
			sampleIDInclude = as.vector(sampleIDInclude[!duplicated(sampleIDInclude), ])
			cat(length(sampleIDInclude), " non-duplicated sample IDs were found in SampleIDIncludeFile\n")
			data = data[which(as.vector(data[,which(colnames(data) == sampleIDColinphenoFile)]) %in% sampleIDInclude),,drop=F]
			cat(nrow(data), " samples in sampleIDInclude have non-missing phenotypes and covariates\n")
		}	

	}	


	if(length(qCovarCol) > 0){
	    cat(qCovarCol, "are categorical covariates\n")
	    if(!all(qCovarCol %in% covarColList)){
		stop("ERROR! all covariates in qCovarCol must be in covarColList\n")
	    }else{
		for(q in qCovarCol){
			data[,q] = as.factor(data[,q])
		}
	    }
	}


        if (FemaleOnly | MaleOnly) {
            if (!sexCol %in% colnames(data)) {
                stop("ERROR! column for sex ", sexCol, " does not exist in the phenoFile \n")
            }
            else {
                if (FemaleOnly) {
                  data = data[which(data[, which(colnames(data) == 
                    sexCol)] == FemaleCode), ]
                  if (nrow(data) == 0) {
                    stop("ERROR! no samples in the phenotype are coded as ", 
                      FemaleCode, " in the column ", sexCol, 
                      "\n")
                  }
                }
                else if (MaleOnly) {
                  data = data[which(data[, which(colnames(data) == 
                    sexCol)] == MaleCode), ]
                  if (nrow(data) == 0) {
                    stop("ERROR! no samples in the phenotype are coded as ", 
                      MaleCode, " in the column ", sexCol, "\n")
                  }
                }
            }
        }



        if (length(covarColList) > 0) {
            formula = paste0(phenoCol, "~", paste0(covarColList, 
                collapse = "+"))
            hasCovariate = TRUE
        }else {
            formula = paste0(phenoCol, "~ 1")
            hasCovariate = FALSE
        }

        cat("formula is ", formula, "\n")
        formula.null = as.formula(formula)
        mmat = model.matrix(formula.null, data, na.action = NULL)
	mmat = data.frame(mmat)
        mmat = cbind(mmat,  data[, which(colnames(data) == phenoCol), drop=F])
	colnames(mmat)[ncol(mmat)] = phenoCol


	if (length(covarColList) > 0) {
	    if(length(qCovarCol) > 0){
		covarColList = colnames(mmat)[2:(ncol(mmat)-1)]
		formula = paste0(phenoCol, "~", paste0(covarColList, collapse = "+")) 
		formula.null = as.formula(formula)
            }
	}

        mmat$IID = data[, which(sampleIDColinphenoFile == colnames(data))]
	if(longlCol != ""){
		mmat$longlVar = data[, which(longlCol == colnames(data))]
	}

        mmat_nomissing = mmat[complete.cases(mmat), ]
        mmat_nomissing$IndexPheno = seq(1, nrow(mmat_nomissing), 
            by = 1)
        cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")



        if(sparseGRMSampleIDFile != ""){
	#if((useSparseGRMtoFitNULL & !skipVarianceRatioEstimation) | useSparseGRMforVarRatio){
		sampleListwithGenov0 = data.table:::fread(sparseGRMSampleIDFile,
                header = F, , colClasses = c("character"), data.table = F)
                colnames(sampleListwithGenov0) = c("IIDgeno")
		cat(length(sampleListwithGenov0$IIDgeno), " samples are in the sparse GRM\n")
		mmat_nomissing = mmat_nomissing[which(mmat_nomissing$IID %in% sampleListwithGenov0$IIDgeno), ]
		cat(nrow(mmat_nomissing), " samples who have non-missing phenotypes are also in the sparse GRM\n")
	}


      if(longlCol == ""){	
	if(any(duplicated(mmat_nomissing$IID))){
		cat("Duplicated sample IDs are detected in the phenotype file. Assuming repeated measurements\n")
	}
      }else{
	cat("Longitudinal variable ", longlCol, " is specified\n")      
	if(!any(duplicated(mmat_nomissing$IID))){
		stop("No duplicated sample IDs are detected in the phenotype file\n")
	}	
      }	      
	if(useGRMtoFitNULL){
        	dataMerge = merge(mmat_nomissing, sampleListwithGeno, 
            by.x = "IID", by.y = "IIDgeno")
        	dataMerge_sort = dataMerge[with(dataMerge, order(IndexGeno)),]
        #dataMerge_sort = dataMerge[with(dataMerge, order(IndexPheno)),]
	}else{
		dataMerge_sort = mmat_nomissing
		dataMerge_sort$IIDgeno = dataMerge_sort$IID

	}	

	print("Test")
	print(head(dataMerge_sort))
	#dataMerge_sort = dataMerge

        rm(mmat)
	rm(mmat_nomissing)
	gc()
	if(useGRMtoFitNULL){	
		indicatorGenoSamplesWithPheno = (sampleListwithGeno$IndexGeno %in% dataMerge_sort$IndexGeno)

        	if (length(unique(dataMerge_sort$IIDgeno)) < length(unique(sampleListwithGeno$IIDgeno))) {
            	cat(length(unique(sampleListwithGeno$IIDgeno)) - length(unique(dataMerge_sort$IIDgeno)), 
                " samples in geno file do not have phenotypes\n")
        	}
        	cat(length(unique(dataMerge_sort$IIDgeno)), " samples will be used for analysis\n")
	}else{
		indicatorGenoSamplesWithPheno = rep(TRUE, nrow(dataMerge_sort))

	}	

	if(any(duplicated(dataMerge_sort$IID))){
		cat(nrow(dataMerge_sort), " observations will be used for analysis\n")
		#if(longlCol != ""){
		  #dataMerge_sort = dataMerge_sort[with(dataMerge_sort, order(IndexGeno, longlVar)),]
		#  dataMerge_sort = dataMerge_sort[with(dataMerge_sort, order(longlVar)),]
		#}
		set_I_mat_inR(dataMerge_sort$IID)
		if(longlCol != ""){
			set_T_mat_inR(dataMerge_sort$IID, dataMerge_sort$longlVar)
		}
	
	}else{
		if(!useGRMtoFitNULL){
			stop("No duplicated IDs are observed in the phenotype file, so GRM must be used to fit the null model. Please set useGRMtoFitNULL=TRUE\n")
		}
	}
	set_useGRMtoFitNULL(useGRMtoFitNULL)	
    }



    print("Test3")
    print(head(dataMerge_sort))

    if (traitType == "quantitative" & invNormalize) {
        cat("Perform the inverse nomalization for ", phenoCol, 
            "\n")
        invPheno = qnorm((rank(dataMerge_sort[, which(colnames(dataMerge_sort) == 
            phenoCol)], na.last = "keep") - 0.5)/sum(!is.na(dataMerge_sort[, 
            which(colnames(dataMerge_sort) == phenoCol)])))
        dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] = invPheno
    }
     print("Test4")
    print(head(dataMerge_sort))
    if (traitType == "binary" & (length(covarColList) > 0)) {
        out_checksep = checkPerfectSep(formula.null, data = dataMerge_sort, 
            minCovariateCount)
        covarColList <- covarColList[!(covarColList %in% out_checksep)]
        formula = paste0(phenoCol, "~", paste0(covarColList, 
            collapse = "+"))
        formula.null = as.formula(formula)
        if (length(covarColList) == 1) {
            hasCovariate = FALSE
        }
        else {
            hasCovariate = TRUE
        }
        dataMerge_sort <- dataMerge_sort[, !(names(dataMerge_sort) %in% 
            out_checksep)]
    }
    if (!hasCovariate) {
	print("No covariate is includes so isCovariateOffset = FALSE")
        isCovariateOffset = FALSE
    }

    if (isCovariateTransform & hasCovariate) {
        cat("qr transformation has been performed on covariates\n")
        out.transform <- Covariate_Transform(formula.null, data = dataMerge_sort)
        formulaNewList = c(phenoCol, " ~ ", out.transform$Param.transform$X_name[1])
        if (length(out.transform$Param.transform$X_name) > 1) {
            for (i in c(2:length(out.transform$Param.transform$X_name))) {
                formulaNewList = c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
            }
        }
        formulaNewList = paste0(formulaNewList, collapse = "")
        formulaNewList = paste0(formulaNewList, "-1")
        formula.new = as.formula(paste0(formulaNewList, collapse = ""))
        data.new = data.frame(cbind(out.transform$Y, out.transform$X1))
        colnames(data.new) = c(phenoCol, out.transform$Param.transform$X_name)
        cat("colnames(data.new) is ", colnames(data.new), "\n")
        cat("out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr), 
            "\n")
    }else{
        formula.new = formula.null
        data.new = dataMerge_sort
        out.transform = NULL
    }

    mmat = model.matrix(formula.new, data=data.new, na.action = NULL)
    if (traitType == "binary") {
            modwitcov = glm(formula.new, data = data.new, 
                family = binomial)
    }else if(traitType == "quantitative"){
	    modwitcov = glm(formula.new, data = data.new,
                family = gaussian(link = "identity"))	
    }else if(traitType == "count"){
	     modwitcov = glm(formula.new, data = data.new,
	                     family = "poisson")
    }


    if (isCovariateOffset) {
    	covoffset = mmat[,-1, drop=F] %*%  modwitcov$coefficients[-1]  
        print("isCovariateOffset=TRUE, so fixed effects coefficnets won't be estimated.")
	formula.new.withCov = formula.new
        formula_nocov = paste0(phenoCol, "~ 1")
        formula.new = as.formula(formula_nocov)
        hasCovariate = FALSE
    }else{
	covoffset = rep(0,nrow(data.new))
    }	    
    data.new$covoffset = covoffset


    if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
	getsubGRM_orig(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, dataMerge_sort$IID)
   	#getsubGRM(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, dataMerge_sort$IID, dataMerge_sort$longlVar)    
        #m4 = gen_sp_v2(sparseGRMtest)
        #cat("Setting up sparse GRM using ", sparseGRMFile, " and ", sparseGRMSampleIDFile, "\n")
        #cat("Dimension of the sparse GRM is ", dim(m4), "\n")
        #A = summary(m4)
        #locationMatinR = rbind(A$i - 1, A$j - 1)
        #valueVecinR = A$x
        #setupSparseGRM(dim(m4)[1], locationMatinR, valueVecinR)
	#setupSparseGRM_new(sparseGRMtest)
        #rm(sparseGRMtest)
	gc()
    }

    #allow for multiple variance components
    #set_Vmat_vec(VmatFilelist, VmatSampleFilelist, dataMerge_sort$IID, dataMerge_sort$longlVar)
    set_Vmat_vec_orig(VmatFilelist, VmatSampleFilelist, dataMerge_sort$IID)

    numofV = get_numofV()
    cat("numofV ", numofV, "\n")
    if(any(duplicated(dataMerge_sort$IID))){
        if(longlCol == ""){
	    if(useGRMtoFitNULL){	
	       num_Kmat = numofV + 3
            }else{
	       num_Kmat = numofV + 2	
	    }	    
               #k = num_Kmat + 3
        }else{
	    if(useGRMtoFitNULL){	
	       num_Kmat = 7 + numofV*3
	    }else{
	       num_Kmat = 4 + numofV*3
	    }	    
	       #k = num_Kmat + 2		
        }
    }else{
	if(useGRMtoFitNULL){    
	    num_Kmat = numofV + 2
        }else{
	    num_Kmat = numofV + 1	
	}	
        #k = 2
    }
    k = num_Kmat



    set_num_Kmat(num_Kmat)
    cat("num_Kmat ", num_Kmat, "\n")


    if(longlCol != ""){
       covarianceIdxMat = set_covarianceidx_Mat()
    }else{
       covarianceIdxMat = NULL
    }	    

    if(!skipVarianceRatioEstimation){
            isVarianceRatioinGeno = TRUE
            if(isCateVarianceRatio){
                minMAC_varRatio = min(cateVarRatioMinMACVecExclude)
                maxMAC_varRatio = max(cateVarRatioMaxMACVecInclude)
                cat("Categorical variance ratios will be estimated. Please make sure there are at least 200 markers in each MAC category.\n")
            }else{
                minMAC_varRatio = 20
                maxMAC_varRatio = -1 #will randomly select markers from the plink file and leave them out when constructing GRM
            }
            setminMAC_VarianceRatio(minMAC_varRatio, maxMAC_varRatio, isVarianceRatioinGeno)
    }

    #set up parameters
    if (minMAFforGRM > 0) {
        cat("Markers in the Plink file with MAF < ", minMAFforGRM,
            " will be removed before constructing GRM\n")
    }
    if (maxMissingRateforGRM > 0){
        cat("Markers in the Plink file with missing rate > ", maxMissingRateforGRM, " will be removed before constructing GRM\n")
    }

    setminMAFforGRM(minMAFforGRM)
    setmaxMissingRateforGRM(maxMissingRateforGRM)



    if (traitType == "binary") {
        cat(phenoCol, " is a binary trait\n")
        uniqPheno = sort(unique(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)]))
        if (uniqPheno[1] != 0 | uniqPheno[2] != 1) {
            stop("ERROR! phenotype value needs to be 0 or 1 \n")
        }
	print("formula.new")
	print(formula.new)
	print("head(data.new)")
	print(head(data.new))
        if (!isCovariateOffset) {
            fit0 = glm(formula.new, data = data.new, family = binomial)
	    Xorig = NULL
	    
	    }else{
	    fit0orig = glm(formula.new.withCov, data = data.new, family = binomial)
	    Xorig = model.matrix(fit0orig)
	    rm(fit0orig)
	    gc()
            fit0 = glm(formula.new, data = data.new, offset = covoffset, 
                family = binomial)
        }

    }else if(traitType == "quantitative"){
        cat(phenoCol, " is a quantitative trait\n")
        if (!isCovariateOffset) {
            fit0 = glm(formula.new, data = data.new, family = gaussian(link = "identity"))
            Xorig = NULL
        }else{
            fit0orig = glm(formula.new.withCov, data = data.new, family = gaussian(link = "identity"))
            Xorig = model.matrix(fit0orig)
            rm(fit0orig)
            gc()
            fit0 = glm(formula.new, data = data.new, offset = covoffset,
            family = gaussian(link = "identity"))
        }	
    }else if(traitType == "count"){
        cat(phenoCol, " is a count trait\n")
        miny = min(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)])
        if (miny < 0 ) {
            stop("ERROR! phenotype value needs to be non-negative \n")
        }

        if (!isCovariateOffset) {
            fit0 = glm(formula.new, data = data.new, family = "poisson")
            Xorig = NULL
        }else{
            fit0orig = glm(formula.new.withCov, data = data.new, family = "poisson")
            Xorig = model.matrix(fit0orig)
            rm(fit0orig)
            gc()
            fit0 = glm(formula.new, data = data.new, offset = covoffset,
            family = "poisson")
        }
    }
    

    cat("glm:\n")
    print(fit0)
    obj.noK = NULL


    #if(length(fit0$y) > 20000){
    isStoreSigma = FALSE
    #}else{
    #  isStoreSigma = TRUE
    #}
    print("isStoreSigma")
    print(isStoreSigma)
    set_store_sigma(isStoreSigma)    

    if (!skipModelFitting) {
        #setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
        cat("Start fitting the NULL GLMM\n")
        t_begin = proc.time()
        print(t_begin)


        tau = rep(0, k) 
	fixtau = rep(0, k)
	tauInit = tau

	set_isSparseGRM(useSparseGRMtoFitNULL)
        set_useGRMtoFitNULL(useGRMtoFitNULL)

        system.time(modglmm <- glmmkin.ai_PCG_Rcpp_multiV(bedFile, bimFile, famFile, Xorig, isCovariateOffset, 
                fit0, tau = tau, fixtau = fixtau, maxiter = maxiter, 
                tol = tol, verbose = TRUE, nrun = nrun, tolPCG = tolPCG, 
                maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno, 
                obj.noK = obj.noK, out.transform = out.transform, 
                tauInit = tauInit, memoryChunk = memoryChunk, 
                LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec, 
                chromosomeEndIndexVec = chromosomeEndIndexVec, 
                traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform, 
                isDiagofKinSetAsOne = isDiagofKinSetAsOne, 
		isLowMemLOCO = isLowMemLOCO, covarianceIdxMat = covarianceIdxMat, isStoreSigma = isStoreSigma, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL, useGRMtoFitNULL = useGRMtoFitNULL))
         
	modglmm$obj.glm.null$model <- data.frame(modglmm$obj.glm.null$model)
            
	for (x in names(modglmm$obj.glm.null)) {
            attr(modglmm$obj.glm.null[[x]], ".Environment") <- c()
        }
	    #modglmm$offset = covoffset
	if(isCovariateOffset){
	    modglmm$offset = covoffset
	 }else{
	    if(hasCovariate){
	       #data.new = data.new[,-c(1,2), drop=F]
	       #data.new = as.matrix(data.new[,-ncol(data.new), drop=F])
	       data.new.X = model.matrix(fit0)[,-1,drop=F]
	       print(head(data.new))
	       print(head(data.new.X))
	       print(head(modglmm$coefficients[-1]))
	       modglmm$offset = data.new.X%*%(as.vector(modglmm$coefficients[-1]))
	       if(LOCO){
	           for(j in 1:22){	
			if(modglmm$LOCOResult[[j]]$isLOCO){
			    modglmm$LOCOResult[[j]]$offset = data.new.X %*%(as.vector(modglmm$LOCOResult[[j]]$coefficients[-1]))
			}		       
		   }
	       }
	    }else{
		modglmm$offset = covoffset
	    }	
	}
       

        if(skipVarianceRatioEstimation & useSparseGRMtoFitNULL){
                family = fit0$family
                eta = modglmm$linear.predictors
                mu = modglmm$fitted.values
                mu.eta = family$mu.eta(eta)
                sqrtW = mu.eta/sqrt(family$variance(mu))
                W = sqrtW^2
                tauVecNew = modglmm$theta
                Sigma_iX =  getSigma_X_multiV(W, tauVecNew, modglmm$X, maxiterPCG, tolPCG, LOCO=FALSE)
                Sigma_iXXSigma_iX = Sigma_iX%*%(solve(t(modglmm$X)%*%Sigma_iX))
                modglmm$Sigma_iXXSigma_iX = Sigma_iXXSigma_iX
        }
    
        modglmm$useSparseGRMforVarRatio = useSparseGRMforVarRatio 	    
        #save(modglmm, file = modelOut)
        tau = modglmm$theta
        alpha0 = modglmm$coefficients

        if(!is.null(out.transform) & is.null(fit0$offset)){
                coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
                modglmm$coefficients = coef.alpha
        }


        if(LOCO & isLowMemLOCO){
                modglmm$LOCOResult = NULL
                modglmm$LOCO = FALSE
                chromosomeStartIndexVec = modglmm$chromosomeStartIndexVec
                chromosomeEndIndexVec = modglmm$chromosomeEndIndexVec
                modglmm$chromosomeStartIndexVec = NULL
                modglmm$chromosomeEndIndexVec = NULL
                modelOut = paste(c(outputPrefix,"_noLOCO.rda"), collapse="")
                save(modglmm, file = modelOut)
                modglmm$LOCO = TRUE
                modglmm$Y = NULL
                eta0 = modglmm$linear.predictors
                modglmm$linear.predictors = NULL
                modglmm$coefficients = NULL
                modglmm$cov = NULL
                modglmm$fitted.values = NULL
                modglmm$residuals = NULL
                modglmm$obj.noK = NULL
                offset0 = modglmm$offset
                modglmm$offset = NULL
                y = fit0$y
                gc()
                #save(modglmm, file = modelOut)
                set_Diagof_StdGeno_LOCO()
                modglmm$LOCOResult = list()
                for(j in 1:22){
                        startIndex = chromosomeStartIndexVec[j]
                        endIndex = chromosomeEndIndexVec[j]
                        if(!is.na(startIndex) && !is.na(endIndex)){
                                cat("leave chromosome ", j, " out\n")
                                setStartEndIndex(startIndex, endIndex, j-1)

                                re.coef_LOCO = Get_Coef_multiV(y, X=model.matrix(fit0), tau, family = fit0$family, alpha = alpha0, eta = eta0,  offset =  offset0, verbose=TRUE, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter, LOCO=TRUE)
                                cov = re.coef_LOCO$cov
                                alpha = re.coef_LOCO$alpha
                                eta = re.coef_LOCO$eta
                                Y = re.coef_LOCO$Y
                                mu = re.coef_LOCO$mu
        			if(family$family == "binomial"){
          				mu2 = mu * (1-mu)
        			}else if(family$family == "poisson"){
          				mu2 = mu
        			}else if(family$family == "gaussian"){
          				mu2 = rep((1/(tau[1])),length(res))
        			}
                                res = y - mu

                                if(!is.null(out.transform) & is.null(fit0$offset)){
                                        coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
                                }else{
                                        coef.alpha = alpha
                                }
                                if(!isCovariateOffset){
                                        obj.noK = ScoreTest_NULL_Model(mu, mu2, y, X)
                                }else{
                                        obj.noK = ScoreTest_NULL_Model(mu, mu2, y, Xorig)
                                }

                                modglmm$LOCOResult[[j]] = list(isLOCO = TRUE, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, obj.noK = obj.noK)
                                if(!isCovariateOffset & hasCovariate){
					 data.new.X = model.matrix(fit0)[,-1,drop=F]
                                         modglmm$LOCOResult[[j]]$offset = data.new.X %*%(as.vector(modglmm$LOCOResult[[j]]$coefficients[-1]))
                                }
                                modelOutbychr = paste(c(outputPrefix,"_chr",j,".rda"), collapse="")
                                if(j!=22){
                                        for(j1 in (j+1):22){
                                                modglmm$LOCOResult[[j1]] = list(NULL)
                                        }
                                }
                                save(modglmm, file = modelOutbychr)
                                modglmm$LOCOResult[[j]] = list(NULL)
                                gc()
                        }else{
                                modglmm$LOCOResult[[j]] = list(isLOCO = FALSE)
                        }

                }
                gc()
                #modelOut_nonauto = paste(c(outputPrefix,"_noLOCO.rda"), collapse="")
           }else{
            save(modglmm, file = modelOut)
           }


        t_end = proc.time()
        print(t_end)
        cat("t_end - t_begin, fitting the NULL model took\n")
        print(t_end - t_begin)
    }else{
	cat("Skip fitting the NULL GLMM\n")
	if(!file.exists(modelOut)){
            stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
        }
        load(modelOut)
        if(is.null(modglmm$LOCO)) {
           modglmm$LOCO = FALSE
        }
        setgeno(bedFile, bimFile, famFile, dataMerge_sort$IndexGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
        tau = modglmm$theta	    
        setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
    }

    if(!skipVarianceRatioEstimation) {
	if(LOCO){
    	    MsubIndVec = getQCdMarkerIndex()
            print(length(MsubIndVec))
    	    chrVec = data.table:::fread(bimFile, header = F)[,1]
	    print(length(chrVec))
   	    chrVec = chrVec[which(MsubIndVec == TRUE)]
    	    updatechrList = updateChrStartEndIndexVec(chrVec)
    	    LOCO = updatechrList$LOCO
    	    chromosomeStartIndexVec = updatechrList$chromosomeStartIndexVec
    	    chromosomeEndIndexVec = updatechrList$chromosomeEndIndexVec
            set_Diagof_StdGeno_LOCO()
  	}
        cat("Start estimating variance ratios\n")
	load(modelOut)
	extractVarianceRatio_multiV(obj.glmm.null = modglmm,
                obj.glm.null = fit0, maxiterPCG = maxiterPCG,
                tolPCG = tolPCG, numMarkers = numMarkersForVarRatio, varRatioOutFile = varRatioFile,
                ratioCVcutoff = ratioCVcutoff, testOut = SPAGMMATOut,
                bedFile=bedFile, bimFile=bimFile, famFile=famFile, chromosomeStartIndexVec = chromosomeStartIndexVec,
                chromosomeEndIndexVec = chromosomeEndIndexVec,
                isCateVarianceRatio = isCateVarianceRatio, cateVarRatioIndexVec = cateVarRatioIndexVec,
                useSparseGRMforVarRatio = useSparseGRMforVarRatio, sparseGRMFile = sparseGRMFile,
                sparseGRMSampleIDFile = sparseGRMSampleIDFile,
                numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
                relatednessCutoff = relatednessCutoff, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL,
                nThreads = nThreads, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
                cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
                minMAFforGRM = minMAFforGRM, isDiagofKinSetAsOne = isDiagofKinSetAsOne,
                includeNonautoMarkersforVarRatio = includeNonautoMarkersforVarRatio, isStoreSigma = isStoreSigma, useGRMtoFitNULL = useGRMtoFitNULL)
    }else{
        cat("Skip estimating variance ratios\n")
    }
    closeGenoFile_plink()

}



extractVarianceRatio_multiV = function(obj.glmm.null,
                                                    obj.glm.null,
                                                    maxiterPCG = 500,
                                                    tolPCG = 0.01,
                                                    numMarkers,
                                                    varRatioOutFile,
                                                    ratioCVcutoff,
                                                    testOut,
                                                    bedFile,
						    bimFile,
						    famFile,
                                                    chromosomeStartIndexVec,
                                                    chromosomeEndIndexVec,
                                                    isCateVarianceRatio,
                                                    cateVarRatioIndexVec,
                                                    useSparseGRMforVarRatio,
                                                    sparseGRMFile,
                                                    sparseGRMSampleIDFile,
                                                    numRandomMarkerforSparseKin,
                                                    relatednessCutoff,
                                                    useSparseGRMtoFitNULL,
                                                    nThreads,
                                                    cateVarRatioMinMACVecExclude,
                                                    cateVarRatioMaxMACVecInclude,
                                                    minMAFforGRM,
                                                    isDiagofKinSetAsOne,
                                                    includeNonautoMarkersforVarRatio,
						    isStoreSigma = FALSE, 
						    useGRMtoFitNULL = TRUE){

  obj.noK = obj.glmm.null$obj.noK
  if(file.exists(testOut)){file.remove(testOut)}
  bimPlink = data.frame(data.table:::fread(bimFile, header=F))
  if(sum(sapply(bimPlink[,1], is.numeric)) != nrow(bimPlink)){
    stop("ERROR: chromosome column in plink bim file is no numeric!\n")
  }

  family = obj.glm.null$family
  print(family)
  eta = obj.glmm.null$linear.predictors
  mu = obj.glmm.null$fitted.values
  mu.eta = family$mu.eta(eta)
  sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
  W = sqrtW^2   ##(mu*(1-mu) for binary)
  tauVecNew = obj.glmm.null$theta
  if(isStoreSigma){
         gen_sp_Sigma_multiV(W, tauVecNew)
  }
  X = obj.glmm.null$X


  set_isSparseGRM(useSparseGRMtoFitNULL)
  set_useGRMtoFitNULL(useGRMtoFitNULL)
  if(!useGRMtoFitNULL){
    useSparseGRMforVarRatio = FALSE
  }	  
  Sigma_iX_noLOCO = getSigma_X_multiV(W, tauVecNew, X, maxiterPCG, tolPCG, LOCO=FALSE)


  y = obj.glm.null$y

  if(any(duplicated(obj.glmm.null$sampleID))){
    dupSampleIndex = as.numeric(factor(obj.glmm.null$sampleID, levels =  unique(obj.glmm.null$sampleID)))
  }


  ##randomize the marker orders to be tested
 if(FALSE){
  if(useSparseGRMtoFitNULL | useSparseGRMforVarRatio){
    sparseSigma = getSparseSigma(bedFile = bedFile, bimFile = bimFile, famFile = famFile,
                outputPrefix=varRatioOutFile,
                sparseGRMFile=sparseGRMFile,
                sparseGRMSampleIDFile=sparseGRMSampleIDFile,
                numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
                relatednessCutoff = relatednessCutoff,
                minMAFforGRM = minMAFforGRM,
                nThreads = nThreads,
                isDiagofKinSetAsOne = isDiagofKinSetAsOne,
                obj.glmm.null = obj.glmm.null,
                W=W, tauVecNew=tauVecNew)
    if(length(tauVecNew)  > 2){
    	sparseSigma = sparseSigma + getProdTauKmat(tauVecNew[3:length(tauVecNew)])
    }

  }
 }


  mMarkers = gettotalMarker()
  listOfMarkersForVarRatio = list()
  MACvector = getMACVec()
  isVarianceRatioinGeno = getIsVarRatioGeno()

  if(isVarianceRatioinGeno){
         MACvector_forVarRatio =  getMACVec_forVarRatio()
         Indexvector_forVarRatio =  getIndexVec_forVarRatio()
         cat("length(MACvector): ", length(MACvector), "\n")
         cat("length(MACvector_forVarRatio): ", length(MACvector_forVarRatio), "\n")

         if(length(MACvector_forVarRatio) > 0){

                #MACdata = data.frame(MACvector = c(MACvector, MACvector_forVarRatio), geno_ind = c(rep(0, length(MACvector)), rep(1, length(MACvector_forVarRatio))), indexInGeno = c(seq(1,length(MACvector)), seq(1,length(MACvector_forVarRatio))))
		MACdata = data.frame(MACvector = MACvector_forVarRatio, geno_ind = rep(1, length(MACvector_forVarRatio)), indexInGeno = seq(1,length(MACvector_forVarRatio)))	
        }else{
                stop("No markers were found for variance ratio estimation. Please make sure there are at least 200 markers in each MAC category\n")
        }

  }else{
         MACdata = data.frame(MACvector = MACvector, geno_ind = rep(0, length(MACvector)), indexInGeno = seq(1,length(MACvector)))
  }

  if(!isCateVarianceRatio){
    cat("Only one variance ratio will be estimated using randomly selected markers with MAC >= 20\n")
    MACindex = 1:nrow(MACdata)
    listOfMarkersForVarRatio[[1]] = sample(MACindex, size = length(MACindex), replace = FALSE)
    cateVarRatioIndexVec=c(1)
  }else{
    cat("Categorical variance ratios will be estimated.\n")

    if(is.null(cateVarRatioIndexVec)){cateVarRatioIndexVec = rep(1, length(cateVarRatioMinMACVecExclude))}
    numCate = length(cateVarRatioIndexVec)
    for(i in 1:(numCate-1)){
      MACindex = which(MACdata$MACvector > cateVarRatioMinMACVecExclude[i] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[i])
      listOfMarkersForVarRatio[[i]] = sample(MACindex, size = length(MACindex), replace = FALSE)
    }

    if(length(cateVarRatioMaxMACVecInclude) == (numCate-1)){
      MACindex = which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate])
    }else{
      MACindex = which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[numCate])
    }

    listOfMarkersForVarRatio[[numCate]] = sample(MACindex, size = length(MACindex), replace = FALSE)

    for(k in 1:length(cateVarRatioIndexVec)){
      if(k <= length(cateVarRatioIndexVec)-1){
        if(cateVarRatioIndexVec[k] == 1){
          cat(cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k],"\n")
          if(length(listOfMarkersForVarRatio[[k]]) < numMarkers){
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], " is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }else{
        if(cateVarRatioIndexVec[k] == 1){
          cat(cateVarRatioMinMACVecExclude[k], "< MAC\n")
          if(length(listOfMarkersForVarRatio[[k]]) < numMarkers){
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC  is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }
    }


  }# if(!isCateVarianceRatio){

  freqVec = getAlleleFreqVec()

  Nnomissing = length(mu)
  varRatioTable = NULL


  for(k in 1:length(listOfMarkersForVarRatio)){
    if(cateVarRatioIndexVec[k] == 1){

      numMarkers0 = numMarkers
      varRatio_sparseGRM_vec = NULL
      varRatio_NULL_vec = NULL
      indexInMarkerList = 1
      numTestedMarker = 0
      ratioCV = ratioCVcutoff + 0.1

      while(ratioCV > ratioCVcutoff){
        while(numTestedMarker < numMarkers0){
          macdata_i = listOfMarkersForVarRatio[[k]][indexInMarkerList]
          i = (MACdata$indexInGeno)[macdata_i]
          genoInd = (MACdata$geno_ind)[macdata_i]
          cat(i, "th marker in geno ", genoInd, "\n")
          cat("MAC: ", (MACdata$MACvector)[macdata_i], "\n")
          if(genoInd == 0){
                G0 = Get_OneSNP_Geno(i-1)
          }else if(genoInd == 1){
                G0 = Get_OneSNP_Geno_forVarRatio(i-1)
          }


          cat("G0", G0[1:10], "\n")
          CHR = bimPlink[Indexvector_forVarRatio[i]+1,1]
	  cat("CHR ", CHR, "\n")
          #if(sum(G0)/(2*Nnomissing) > 0.5){
          if(sum(G0)/(2*length(G0)) > 0.5){
            G0 = 2-G0
          }


          #if(any(duplicated(obj.glmm.null$sampleID))){
	  #	G0 = G0[dupSampleIndex]		
          #}		   
          NAset = which(G0==0)
          AC = sum(G0)

         indexInMarkerList = indexInMarkerList + 1
         if((CHR >= 1 & CHR <= 22) | includeNonautoMarkersforVarRatio){
          AF = AC/(2*Nnomissing)
          if(CHR >= 1 & CHR <= 22){
               autoMarker=TRUE
          }else{
               autoMarker=FALSE
          }

          G = G0  -  obj.noK$XXVX_inv %*%  (obj.noK$XV %*% G0) # G1 is X adjusted
          g = G/sqrt(AC)
          q = innerProduct(g,y)
          eta = obj.glmm.null$linear.predictors
          mu = obj.glmm.null$fitted.values
          mu.eta = family$mu.eta(eta)
          sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
          W = sqrtW^2
          print("W[1:10]")
          print(W[1:10])
          print("mu[1:10]")
          print(mu[1:10])


	  set_isSparseGRM(useSparseGRMtoFitNULL)
  	  set_useGRMtoFitNULL(useGRMtoFitNULL)

          Sigma_iG = getSigma_G_multiV(W, tauVecNew, G, maxiterPCG, tolPCG, LOCO=FALSE)
          Sigma_iX = Sigma_iX_noLOCO

          var1a = t(G)%*%Sigma_iG - t(G)%*%Sigma_iX%*%(solve(t(X)%*%Sigma_iX))%*%t(X)%*%Sigma_iG
          var1 = var1a/AC
	  cat("AC ", AC, "\n")	
          m1 = innerProduct(mu,g)
	  S = q-m1

   if(useSparseGRMforVarRatio){
	set_isSparseGRM(useSparseGRMforVarRatio)
	Sigma_iG = getSigma_G_noV(W, tauVecNew, G, maxiterPCG, tolPCG, LOCO=FALSE)
	var2_a = t(g) %*% Sigma_iG/(sqrt(AC))
	var2sparseGRM = var2_a[1,1]
	cat("var2sparseGRM ", var2sparseGRM*AC, "\n")
	varRatio_sparseGRM_vec = c(varRatio_sparseGRM_vec, var1/var2sparseGRM)

  }else{
        varRatio_sparseGRM_vec = c(varRatio_sparseGRM_vec,1)
  }	

   p_exact = pchisq(S^2/var1, df=1, lower.tail=F)



    if(obj.glmm.null$traitType == "binary"){
         var2null = innerProduct(mu*(1-mu), g*g)
    varRatio_NULL_vec = c(varRatio_NULL_vec, var1/var2null)
    }else if(obj.glmm.null$traitType == "quantitative"){
         var2null = innerProduct(g, g)
    varRatio_NULL_vec = c(varRatio_NULL_vec, var1/var2null)
    }else if(obj.glmm.null$traitType == "count"){
         var2null = innerProduct(mu, g*g)
    varRatio_NULL_vec = c(varRatio_NULL_vec, var1/var2null)
    }
         p_approx = pchisq(S^2/(var2null), df=1, lower.tail=F)
         p_approx_true = pchisq(S^2/(var2null*0.0009604056), df=1, lower.tail=F)

    cat("AC ", AC, "\n")
    cat("S ", S, "\n")
    cat("var1 ", var1*AC, "\n")
    cat("var2null ", var2null*AC, "\n")
    cat("p_exact ", p_exact, "\n")
    cat("p_approx ", p_approx, "\n")
    cat("p_approx_true ", p_approx_true, "\n")
    varRatio_NULL_vec = c(varRatio_NULL_vec, var1/var2null)

    #indexInMarkerList = indexInMarkerList + 1
    numTestedMarker = numTestedMarker + 1

      }else{
        indexInMarkerList = indexInMarkerList + 1
      }

      if(indexInMarkerList-1 == length(listOfMarkersForVarRatio[[k]])){
        numTestedMarker = numMarkers0
      }
    }#end of while(numTestedMarker < numMarkers)


    print("varRatio_NULL_vec")
    print(varRatio_NULL_vec)

    ratioCV = calCV(varRatio_NULL_vec)

    if(ratioCV > ratioCVcutoff){
      cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " > ", ratioCVcutoff, "\n")
      numMarkers0 = numMarkers0 + 10
      cat("try ", numMarkers0, " markers\n")
    }else{
      cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " < ", ratioCVcutoff, "\n")
    }

    if(indexInMarkerList-1 == length(listOfMarkersForVarRatio[[k]])){
      ratioCV = ratioCVcutoff
      cat("no more markers are available in the MAC category ", k, "\n")
      print(indexInMarkerList-1)
    }

  }#end of while(ratioCV > ratioCVcutoff)

  if(length(varRatio_sparseGRM_vec) > 0){
    varRatio_sparse = mean(varRatio_sparseGRM_vec)
    cat("varRatio_sparse", varRatio_sparse, "\n")
    varRatioTable = rbind(varRatioTable, c(varRatio_sparse, "sparse", k))
  }
  varRatio_null = mean(varRatio_NULL_vec)
  cat("varRatio_null", varRatio_null, "\n")
  varRatioTable = rbind(varRatioTable, c(varRatio_null, "null", k))
  }else{# if(cateVarRatioVec[k] == 1)
    varRatioTable = rbind(varRatioTable, c(1, "null", k))
    if(length(varRatio_sparseGRM_vec) > 0){
      varRatioTable = rbind(varRatioTable, c(1, "sparse", k))
    }
  }

} #for(k in 1:length(listOfMarkersForVarRatio)){
  write.table(varRatioTable, varRatioOutFile, quote=F, col.names=F, row.names=F)
  data = read.table(varRatioOutFile, header=F)
  print(data)

}


#Fits the null glmm
glmmkin.ai_PCG_Rcpp_multiV = function(bedFile, bimFile, famFile, Xorig, isCovariateOffset, fit0, tau=c(0,0), fixtau = c(0,0), maxiter =20, tol = 0.02, verbose = TRUE, nrun=30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, indicatorGenoSamplesWithPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne, isLowMemLOCO, covarianceIdxMat = NULL, isStoreSigma = FALSE, useSparseGRMtoFitNULL = TRUE, useGRMtoFitNULL = TRUE) {
  #Fits the null generalized linear mixed model for a poisson, binomial, and gaussian
  #Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix
  #  fit0: glm model. Logistic model output (with no sample relatedness accounted for)
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the glmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  #Returns:
  #  model output for the null glmm

  t_begin = proc.time()
  print(t_begin)
  subSampleInGeno = subPheno$IndexGeno
  if(verbose){
    print("Start reading genotype plink file here")
  }

  #print("subSampleInGeno")
  #print(subSampleInGeno)

  print(head(subPheno))
  set_dup_sample_index(as.numeric(factor(subPheno$IID, levels =  unique(subPheno$IID))))

  #if((!useSparseGRMtoFitNULL & useGRMtoFitNULL) | (skipVarianceRatioEstimation)){
  if(bedFile != ""){
  re1 = system.time({setgeno(bedFile, bimFile, famFile, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)})
  }
  if(verbose){
    print("Genotype reading is done")
  }

  if (LOCO){
    MsubIndVec = getQCdMarkerIndex()
    chrVec = data.table:::fread(bimFile, header = F)[,1]
    chrVec = chrVec[which(MsubIndVec == TRUE)]
    updatechrList = updateChrStartEndIndexVec(chrVec)
    LOCO = updatechrList$LOCO
    chromosomeStartIndexVec = updatechrList$chromosomeStartIndexVec
    chromosomeEndIndexVec = updatechrList$chromosomeEndIndexVec
  }

  y = fit0$y
  n = length(y)
  X = model.matrix(fit0)
  offset = fit0$offset
  if(is.null(offset)){
    offset = rep(0, n)
  }

  family = fit0$family
  eta = fit0$linear.predictors
  mu = fit0$fitted.values
  mu.eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu)/mu.eta
  sqrtW = mu.eta/sqrt(fit0$family$variance(mu))
  W = sqrtW^2

  alpha0 = fit0$coef
  eta0 = eta


  cat("tauInit")
  print(tauInit)

  tau[1:length(tau)] = 0
  if(family$family %in% c("poisson", "binomial")) {
    tau[1] = 1
    fixtau[1] = 1
    tauInit[1] = 1
    idxtau <- which(fixtau == 0)
    cat("fixtau ", fixtau, "\n")
    cat("tauInit ", tauInit, "\n")
    if(sum(tauInit[idxtau]) == 0){
      tau[idxtau] = 0.1
    }else{
      tau[idxtau] = tauInit[idxtau]
    }
  }else{#  if(family$family %in% c("poisson", "binomial")) {
    idxtau <- which(fixtau == 0)	  
    if(sum(tauInit[idxtau]) == 0){
      tau[1] = 1
      #tauInit[1] = 1    
      tau[idxtau] = var(Y)/(length(tau))
      #tau[2] = 0
      #tau[2:length(tau)] = 0
      if (abs(var(Y)) < 0.1){
        stop("WARNING: variance of the phenotype is much smaller than 1. Please consider invNormalize=T\n")
      }
    }else{
      tau[fixtau == 0] = tauInit[fixtau == 0]
    }
  }

  cat("inital tau is ", tau,"\n")

  if(!is.null(covarianceIdxMat)){
	  idxtau2 <- intersect(covarianceIdxMat[, 1], idxtau)
	  print("covarianceIdxMat")
	  print(covarianceIdxMat)
	  print("idxtau2")
	  print(idxtau2)
	  if(length(idxtau2) > 0){
		tau[idxtau2] = 0
          }
	#i_kmat = get_numofV()
    	#if(i_kmat > 0){
        Kmatdiag = getMeanDiagofKmat(LOCO)
	print(Kmatdiag)
    	#}	  
        tau[2:length(tau)] = tau[2:length(tau)]/Kmatdiag
  }

    print("tau")
    print(tau)
	
    if(isStoreSigma){
      gen_sp_Sigma_multiV(W, tau)
    }
    
    re.coef = Get_Coef_multiV(y, X, tau, family, alpha0, eta0,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter, LOCO = FALSE)

    if(isStoreSigma){
      gen_sp_Sigma_multiV(re.coef$W, tau)
    }

    re = getAIScore_multiV(re.coef$Y, X, re.coef$W, tau, fixtau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG,tolPCG = tolPCG, traceCVcutoff = traceCVcutoff, LOCO = FALSE)
    tau0=tau
    tau0_q2 = tau[fixtau == 0]
    #tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace)/n)
    print("tau0_q2 a")
    print(tau0_q2)
    print("idxtau")
    print(idxtau)


    tau_q2 = pmax(0, tau0_q2 + tau0_q2^2 * (re$YPAPY - re$Trace)/n)	    
    tau[idxtau] = tau_q2

    if(!is.null(covarianceIdxMat)){
    	tau[idxtau[which(idxtau %in% idxtau2)]] = 0    
    } 
    print("re$YPAPY")
    print(re$YPAPY)
    print("re$Trace")
    print(re$Trace)

  if(verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }

  for (i in seq_len(maxiter)) {
    W = sqrtW^2

    if(verbose) cat("\nIteration ", i, tau, ":\n")
      alpha0 = re.coef$alpha
      tau0 = tau
      cat("tau0_v1: ", tau0, "\n")
      eta0 = eta
      # use Get_Coef before getAIScore
      t_begin_Get_Coef = proc.time()
      #if(isStoreSigma){
      #  gen_sp_Sigma_multiV(W, tau)
      #}
#cat("eta0 ", eta0, "\n")
      re.coef = Get_Coef_multiV(y, X, tau, family, alpha0, eta0,  offset, verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter, LOCO = FALSE)
      t_end_Get_Coef =  proc.time()
      cat("t_end_Get_Coef - t_begin_Get_Coef\n")
      print(t_end_Get_Coef - t_begin_Get_Coef)
      if(isStoreSigma){
        gen_sp_Sigma_multiV(re.coef$W, tau)
      }
      fit = fitglmmaiRPCG_multiV(re.coef$Y, X, re.coef$W, tau, fixtau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG, tol = tol, traceCVcutoff = traceCVcutoff, LOCO = FALSE)

	if(i == 1){
		save(fit, file="/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_newgit/SAIGE_1.1.3/SAIGE/extdata/SAIGE_fit.rda")
	}	

      t_end_fitglmmaiRPCG= proc.time()
      cat("t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n")
      print(t_end_fitglmmaiRPCG - t_end_Get_Coef)

      tau = as.numeric(fit$tau)
      cov = re.coef$cov
      alpha = re.coef$alpha
      eta = re.coef$eta
      Y = re.coef$Y
      mu = re.coef$mu

     mu.eta = family$mu.eta(eta)
     sqrtW = mu.eta/sqrt(fit0$family$variance(mu))
     W = sqrtW^2



     
      print(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol))
      cat("tau: ", tau, "\n")
      cat("tau0: ", tau0, "\n")

      if(family$family == "gaussian"){
        if(tau[1]<=0){
          stop("ERROR! The first variance component parameter estimate is 0\n")
        }
      }

      if(sum(tau[2:length(tau)]) == 0) break
      
      # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
      if(max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break

      if(max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
        i = maxiter
        break
      }
  }

  if(verbose) cat("\nFinal " ,tau, ":\n")


     if(isStoreSigma){
        gen_sp_Sigma_multiV(W, tau)
     }

    #added these steps after tau is estimated 04-14-2018
  re.coef = Get_Coef_multiV(y, X, tau, family, alpha, eta,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter, LOCO=FALSE)
  cov = re.coef$cov
  alpha = re.coef$alpha
  eta = re.coef$eta
  Y = re.coef$Y
  mu = re.coef$mu

  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu

  if(family$family == "binomial"){
    mu2 = mu * (1-mu)
    traitType = "binary"
  }else if(family$family == "poisson"){
    mu2 = mu
    traitType = "count"
  }else if(family$family == "gaussian"){
    mu2 = rep((1/(tau[1])),length(res))
    traitType = "quantitative"
  }

  #if(isCovariateTransform & hasCovariate){
  #if(!is.null(out.transform) & is.null(fit0$offset)){
  #  coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  #}else{
  #  coef.alpha = alpha
  #}

  #mu2 = mu * (1-mu)

  if(!isCovariateOffset){
    obj.noK = ScoreTest_NULL_Model(mu, mu2, y, X)
    glmmResult = list(theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = subPheno$IID, obj.noK=obj.noK, y = y, X = X, traitType=traitType, isCovariateOffset = isCovariateOffset)
  }else{
    obj.noK = ScoreTest_NULL_Model(mu, mu2, y, Xorig)
    glmmResult = list(theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = subPheno$IID, obj.noK=obj.noK, y = y, X = Xorig, traitType=traitType, isCovariateOffset = isCovariateOffset)
  }


  #LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme

  glmmResult$LOCO = LOCO
  t_end_null = proc.time()
  cat("t_end_null - t_begin, fitting the NULL model without LOCO took\n")
  print(t_end_null - t_begin)
  if(!isLowMemLOCO & LOCO){
     if(isStoreSigma){
        gen_sp_Sigma_multiV(re.coef$W, tau)
     }	  
    set_Diagof_StdGeno_LOCO()
    glmmResult$LOCOResult = list()
    for (j in 1:22){
      startIndex = chromosomeStartIndexVec[j]
      endIndex = chromosomeEndIndexVec[j]
      if(!is.na(startIndex) && !is.na(endIndex)){
        cat("leave chromosome ", j, " out\n")
        setStartEndIndex(startIndex, endIndex, j-1)
        t_begin_Get_Coef_LOCO = proc.time()
        re.coef_LOCO = Get_Coef_multiV(y, X, tau, family, alpha, eta,  offset,verbose=verbose, maxiterPCG=maxiterPCG, tolPCG = tolPCG, maxiter=maxiter, LOCO=TRUE)
        t_end_Get_Coef_LOCO = proc.time()
        cat("t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO\n")
        print(t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO)
        cov = re.coef_LOCO$cov
        alpha = re.coef_LOCO$alpha
        eta = re.coef_LOCO$eta
        Y = re.coef_LOCO$Y
        mu = re.coef_LOCO$mu
        #mu2 = mu * (1-mu)
        #mu2 = mu
        res = y - mu

        if(family$family == "binomial"){
          mu2 = mu * (1-mu)
        }else if(family$family == "poisson"){
          mu2 = mu
        }else if(family$family == "gaussian"){
          mu2 = rep((1/(tau[1])),length(res))
        }

        if(!is.null(out.transform) & is.null(fit0$offset)){
          coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
        }else{
          coef.alpha = alpha
        }
        if(!isCovariateOffset){
          obj.noK = ScoreTest_NULL_Model(mu, mu2, y, X)
        }else{
          obj.noK = ScoreTest_NULL_Model(mu, mu2, y, Xorig)
        }
        glmmResult$LOCOResult[[j]] = list(isLOCO = TRUE, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, obj.noK = obj.noK)
      }else{
        glmmResult$LOCOResult[[j]] = list(isLOCO = FALSE)
      }
    }
  }

  if(isLowMemLOCO & LOCO){
    glmmResult$chromosomeStartIndexVec = chromosomeStartIndexVec
    glmmResult$chromosomeEndIndexVec = chromosomeEndIndexVec
  }
  return(glmmResult)
}
