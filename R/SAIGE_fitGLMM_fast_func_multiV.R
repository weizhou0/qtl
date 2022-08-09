# Run iterations to get converged alpha and eta
Get_Coef_multiV = function(y, X, tau, family, alpha0, eta0,  offset, maxiterPCG, tolPCG,maxiter, verbose=FALSE, LOCO = FALSE){
  tol.coef = 0.1
  mu = family$linkinv(eta0)
  mu.eta = family$mu.eta(eta0)
  Y = eta0 - offset + (y - mu)/mu.eta

  sqrtW = mu.eta/sqrt(family$variance(mu))
  W = sqrtW^2

  for(i in 1:maxiter){
    re.coef = getCoefficients_multiV(Y, X, W, tau, maxiter=maxiterPCG, tol=tolPCG, LOCO)
    alpha = re.coef$alpha
    eta = re.coef$eta + offset

    if(verbose) {
      cat("Tau:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    mu = family$linkinv(eta)
    mu.eta = family$mu.eta(eta)

    Y = eta - offset + (y - mu)/mu.eta

    print("y[1]")
    print(y[1])
    print("y[2]")
    print(y[2])

    sqrtW = mu.eta/sqrt(family$variance(mu))
    W = sqrtW^2

    if( max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol.coef))< tol.coef){
        break
    }
      alpha0 = alpha
    }

    print("alpha")
    print(alpha)
    
    re = list(Y=Y, alpha=alpha, eta=eta, W=W, cov=re.coef$cov, sqrtW=sqrtW, Sigma_iY = re.coef$Sigma_iY, Sigma_iX = re.coef$Sigma_iX, mu=mu)
}


ScoreTest_NULL_Model = function(mu, mu2, y, X){
  V = as.vector(mu2)
  res = as.vector(y - mu)
  XV = t(X * V)
  XVX = t(X) %*% (t(XV))
  XVX_inv = solve(XVX)
  XXVX_inv = X %*% XVX_inv
  XVX_inv_XV = XXVX_inv * V
  S_a =  colSums(X * res)
  re = list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  class(re) = "SA_NULL"
  return(re)
}

##suggested by Shawn 01-19-2018
Covariate_Transform<-function(formula, data){
  X1<-model.matrix(formula,data=data)
#  X1=X1[,c(2:ncol(X1))] #remove intercept
  formula.frame<-model.frame(formula,data=data)
  Y = model.response(formula.frame, type = "any")
  X_name = colnames(X1)

  # First run linear regression to identify multi collinearity
  out.lm<-lm(Y ~ X1 - 1, data=data)
#  out.lm<-lm(Y ~ X1, data=data)
  idx.na<-which(is.na(out.lm$coef))
  if(length(idx.na)> 0){
        print(head(X1))
        X1<-X1[, -idx.na]
        print(head(X1))
        cat("Warning: multi collinearity is detected in covariates! ", X_name[idx.na], " will be excluded in the model\n")
        X_name = X_name[-idx.na]
        #cat("Warning: multi collinearity is detected in covariates! ", X_name[idx.na], " will be excluded in the model\n")
  }
  if(!(1 %in% idx.na)){
    X_name[1] = "minus1"
  }


 # QR decomposition
  Xqr = qr(X1)
  X1_Q = qr.Q(Xqr)
  qrr = qr.R(Xqr)

  N<-nrow(X1)

  # Make square summation=N (so mean=1)
  X1_new<-X1_Q * sqrt(N)
  Param.transform<-list(qrr=qrr, N=N, X_name = X_name, idx.na=idx.na)
  re<-list(Y =Y, X1 = X1_new, Param.transform=Param.transform)
}

# In case to recover original scale coefficients
# X \beta = Q R \beta = (Q \sqrt(N)) ( R \beta / \sqrt(N))
# So coefficient from fit.new is the same as R \beta / \sqrt(N)
Covariate_Transform_Back<-function(coef, Param.transform){
        #coef<-fit.new$coef; Param.transform=out.transform$Param.transform
        coef1<-coef * sqrt(Param.transform$N)
        coef.org<-solve(Param.transform$qrr, coef1)

        names(coef.org)<-Param.transform$X_name
        return(coef.org)
}


pcg<-function (A, b, M=NULL, maxiter = 1e+05, tol = 1e-06){

  # A<-a; b<-c1[,1]; M<-NULL;maxiter = 1e+05; tol = 1e-06
  if (is.null(M)) {
    dA <- diag(A)
    dA[which(dA == 0)] = 1e-04
    #print("dA")
    #print(dA)
    Minv = 1/dA
  } else Minv = solve(M)
  print("R")
  print(Minv)
  x = rep(0, length(b))
  r = b
  if(is.null(M)){
    z = Minv *r
  } else {
    z= Minv %*% r
  }
  p = z
  iter = 0
  sumr2 = sum(r^2)
  while (sumr2 > tol & iter < maxiter) {
    iter = iter + 1
#    cat("iter is ", iter, "\n")
    Ap = crossprod(p, A)[1,]
    a = as.numeric((t(r) %*% z)/(t(p) %*% Ap))
    x = x + a * p
    r1 = r - a * Ap

    if(is.null(M)){
      z1 = Minv * r1
    } else {
      z1 = Minv %*% r1
    }


    bet = as.numeric((t(z1) %*% r1)/(t(z) %*% r))

    p = z1 + bet * p

    z = z1
    r = r1
    sumr2 = sum(r^2)
  }
  if (iter >= maxiter)
    x = "pcg did not converge. You may increase maxiter number."
  return(x)
}



pcgSparse<-function (A, b, M=NULL, maxiter = 1e+05, tol = 1e-06){
  # A<-a; b<-c1[,1]; M<-NULL;maxiter = 1e+05; tol = 1e-06
  if (is.null(M)) {
    dA <- diag(A)
    dA[which(dA == 0)] = 1e-04
    #print("dA")
    #print(dA)
    Minv = 1/dA
  } else Minv = solve(M)
  x = rep(0, length(b))
  r = b
  if(is.null(M)){
    z = Minv *r
  } else {
    z= Minv %*% r
  }
  p = z
  iter = 0
  sumr2 = sum(r^2)
  print("psparse0")
  psparse =  Matrix:::sparseMatrix(i = rep(1, length(p)), j = c(1:length(p)), x = as.vector(p))
  cat("nrow(psparse) ", nrow(psparse), "\n")
  cat("ncol(psparse) ", ncol(psparse), "\n")
  print(class(psparse))

  while (sumr2 > tol & iter < maxiter) {
    iter = iter + 1
    #Ap = crossprod(p, A)[1,]
    print(class(psparse)[1])
    print(class(A)[1])
    print(nrow(psparse))
    print(ncol(psparse))
    print(nrow(A))
    print(ncol(A))

    #Ap = sparse_row_idx_mult(psparse, A)
    Ap = psparse%*%A[1,]
    print("psparse1")
    a = as.numeric((t(r) %*% z)/(t(p) %*% Ap))
    x = x + a * p
    r1 = r - a * Ap

    if(is.null(M)){
      z1 = Minv * r1
    } else {
      z1 = Minv %*% r1
    }


    bet = as.numeric((t(z1) %*% r1)/(t(z) %*% r))

    p = z1 + bet * p

    z = z1
    r = r1
    sumr2 = sum(r^2)
  }
  if (iter >= maxiter)
    x = "pcg did not converge. You may increase maxiter number."
  return(x)
}

getVmatSub <- function (sparseVFile = NULL, sparseVSampleIDFile = "",
    modelID = NULL)
{
  cat("Extract matrix in ", sparseVFile, "\n")
  sparseVLarge = Matrix:::readMM(sparseVFile)
  print(nnzero(sparseVLarge))
  sparseVLarge = sparseVLarge * 1
  sparseVSampleID = data.frame(data.table:::fread(sparseVSampleIDFile,
    header = F, stringsAsFactors = FALSE, colClasses = c("character")))
  colnames(sparseVSampleID) = c("sampleID")
  sparseVSampleID$IndexGRM = seq(1, nrow(sparseVSampleID),by = 1)
  if(nrow(sparseVSampleID) != dim(sparseVLarge)[1] |
     nrow(sparseVSampleID) != dim(sparseVLarge)[2]) {
    stop("ERROR! number of samples in ", sparseVFile, " is not the same to the number of sample IDs in the specified sparseVSampleIDFile ", sparseVSampleIDFile, "\n")
  }else{
    sampleInModel = NULL
    sampleInModel$IID = modelID
    sampleInModel = data.frame(sampleInModel)
    sampleInModel$IndexInModel = seq(1, length(sampleInModel$IID), by = 1)
    cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
    mergeID = merge(sampleInModel, sparseVSampleID, by.x = "IID", by.y = "sampleID")
    if(!any(duplicated(sampleInModel$IID))){ 	
      if (nrow(sampleInModel) > nrow(mergeID)) {
        stop("ERROR: ", nrow(sampleInModel) - nrow(mergeID), "samples are not in the matrix\n")
      }else{
        mergeID = mergeID[with(mergeID, order(IndexInModel)),]
        indexIDofV = mergeID$IndexGRM
        sparseV = sparseVLarge[indexIDofV, indexIDofV]
        rm(sparseVLarge)
        return(sparseV)
      }
    }else{
       if(length(unique(sampleInModel$IID)) > length(unique(mergeID$IID))){
                stop("ERROR: ", length(unique(mergeID$IID)) - length(unique(sampleInModel$IID)),
         "samples are not in the matrix\n")      
       }else{
	mergeID = mergeID[with(mergeID, order(IndexInModel)),]
        indexIDofV = mergeID$IndexGRM
        sparseV = sparseVLarge[indexIDofV, indexIDofV]
        rm(sparseVLarge)
        return(sparseV)
       }	       
    }	    
 }
}


set_Vmat_vec = function(VmatFilelist, VmatSampleFilelist, modelID, longlVarVec=NULL){	
   if(any(duplicated(modelID))){
        VmatsampleID = unique(modelID)
        #sparseVmat0 = Matrix::Diagonal(length(unique(modelID)))
        #sparseVmat0 = sparseVmat0 * 1
        #sparseVmat0 = diag(length(unique(modelID)))
        sparseVmat0 = Matrix:::sparseMatrix(i = c(1: length(VmatsampleID)), j = c(1:length(VmatsampleID)), x = as.vector(rep(1,  length(VmatsampleID))))
        cat("length(unique(modelID)) ", length(unique(modelID)), "\n")
        print(dim(sparseVmat0))
        sparseVmat0 = as(sparseVmat0, "sparseMatrix")
        sampleInModel = data.frame("IID" = modelID, IndexInModel = seq(1, length(modelID)))
        sparseVSampleID = data.frame("sampleID" = VmatsampleID, IndexInMat = seq(1, length(VmatsampleID)))
        mergeID = merge(sampleInModel, sparseVSampleID, by.x = "IID", by.y = "sampleID")
        mergeID = mergeID[with(mergeID, order(IndexInModel)),]
        indexIDofV = mergeID$IndexInMat
        print("okkk5")
        print(length(indexIDofV))
        sparseVmat = sparseVmat0[indexIDofV, indexIDofV]
        print("okkk6")
        addNewKat(sparseVmat)
	print(sparseVmat[1:10,1:10])
	print("okkkkkkk7")
	if(!is.null(longlVarVec)){
		sparseVmat_longl = sparseVmat * longlVarVec
		sparseVmat_longl_b = sparseVmat_longl + t(sparseVmat_longl)
		addNewKat(sparseVmat_longl_b)
		print(sparseVmat_longl_b[1:10,1:10])
		sparseVmat_longl = t(sparseVmat_longl)*longlVarVec
		addNewKat(sparseVmat_longl)
		print(sparseVmat_longl[1:10,1:10])
		rm(sparseVmat_longl)
		rm(sparseVmat_longl_b)
	}
    }
      
  if(VmatFilelist != ""){
    VmatFile_vec = unlist(strsplit(VmatFilelist, split=","))
    cat(length(VmatFile_vec), " additional variance covariance matrices are specified\n")
    VmatSampleFile_vec = unlist(strsplit(VmatSampleFilelist, split=","))
    if(length(VmatSampleFile_vec) != length(VmatFile_vec)){
      stop("Number of sample files in ", VmatSampleFilelist, " does not equal to number of matrix files in ", VmatFilelist, "\n")
    }else{
      for(i in 1:length(VmatFile_vec)){
        Vmatfile = VmatFile_vec[i]
        VmatSamplefile = VmatSampleFile_vec[i]
        if(!file.exists(Vmatfile)){
          stop(Vmatfile, " does not exist\n")
        }else if(!file.exists(VmatSamplefile)){
          stop(VmatSamplefile, " does not exist\n")
        }else{
          sparseVmat = getVmatSub(Vmatfile, VmatSamplefile, modelID)
          addNewKat(sparseVmat)
	  if(!is.null(longlVarVec)){
            sparseVmat_longl = sparseVmat * longlVarVec
	    sparseVmat_longl_b = sparseVmat_longl + t(sparseVmat_longl)
            addNewKat(sparseVmat_longl_b)
	    print(sparseVmat_longl_b[1:10,1:10])
            sparseVmat_longl = t(sparseVmat_longl) * longlVarVec
            addNewKat(sparseVmat_longl)
	    print(sparseVmat_longl[1:10,1:10])
	    rm(sparseVmat_longl_b)
	    rm(sparseVmat_longl)
          }
	  #Matrix::writeMM(sparseVmat, file="kin2_SAIGE.mtx")
        }
      }
    }
  }
}



set_Vmat_vec_orig = function(VmatFilelist, VmatSampleFilelist, modelID){

  if(VmatFilelist != ""){
    VmatFile_vec = unlist(strsplit(VmatFilelist, split=","))
    cat(length(VmatFile_vec), " additional variance covariance matrices are specified\n")
    VmatSampleFile_vec = unlist(strsplit(VmatSampleFilelist, split=","))
    if(length(VmatSampleFile_vec) != length(VmatFile_vec)){
      stop("Number of sample files in ", VmatSampleFilelist, " does not equal to number of matrix files in ", VmatFilelist, "\n")
    }else{
      modelID = unique(modelID)	    
      for(i in 1:length(VmatFile_vec)){
        Vmatfile = VmatFile_vec[i]
        VmatSamplefile = VmatSampleFile_vec[i]
        if(!file.exists(Vmatfile)){
          stop(Vmatfile, " does not exist\n")
        }else if(!file.exists(VmatSamplefile)){
          stop(VmatSamplefile, " does not exist\n")
        }else{          		
          sparseVmat = getVmatSub(Vmatfile, VmatSamplefile, modelID)
          addNewKat(sparseVmat)
          #Matrix::writeMM(sparseVmat, file="kin2_SAIGE.mtx")
        }
      }
    }
  }
}





getsubGRM_orig = function (sparseGRMFile = NULL, sparseGRMSampleIDFile = "", relatednessCutoff,
    modelID = NULL){

    cat("extract sparse GRM\n")
    sparseGRMLarge = Matrix:::readMM(sparseGRMFile)
    print(nnzero(sparseGRMLarge))
    cat("set elements in the sparse GRM <= ", relatednessCutoff,
        " to zero\n")
    sparseGRMLarge = Matrix:::drop0(sparseGRMLarge, tol = relatednessCutoff)
    sparseGRMLarge = sparseGRMLarge * 1
    if (!file.exists(sparseGRMSampleIDFile)) {
        stop("ERROR! sparseSigmaSampleIDFile ", sparseGRMSampleIDFile,
            " does not exist\n")
    }else{
        sparseGRMSampleID = data.frame(data.table:::fread(sparseGRMSampleIDFile,
            header = F, stringsAsFactors = FALSE, colClasses = c("character")))
        colnames(sparseGRMSampleID) = c("sampleID")
        sparseGRMSampleID$IndexGRM = seq(1, nrow(sparseGRMSampleID),
            by = 1)
        if (nrow(sparseGRMSampleID) != dim(sparseGRMLarge)[1] |
            nrow(sparseGRMSampleID) != dim(sparseGRMLarge)[2]) {
            stop("ERROR! number of samples in the sparse GRM is not the same to the number of sample IDs in the specified sparseGRMSampleIDFile ",
                sparseGRMSampleIDFile, "\n")
        }else {
            sampleInModel = NULL
	    modelID = unique(modelID) #model ID not duplicated
            sampleInModel$IID = modelID
            sampleInModel = data.frame(sampleInModel)
            sampleInModel$IndexInModel = seq(1, length(sampleInModel$IID),
                by = 1)
            cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
            mergeID = merge(sampleInModel, sparseGRMSampleID,
                by.x = "IID", by.y = "sampleID")
            if (nrow(sampleInModel) > nrow(mergeID)) {
                stop("ERROR: ", nrow(sampleInModel) - nrow(mergeID),
                  "samples used for model fitting are not in the specified GRM\n")
            }else {
                mergeID = mergeID[with(mergeID, order(IndexInModel)),
                  ]
                indexIDofGRM = mergeID$IndexGRM
                sparseGRM = sparseGRMLarge[indexIDofGRM, indexIDofGRM]
                setupSparseGRM_new(sparseGRM)
                rm(sparseGRMLarge)
                rm(sparseGRM)
                #return(sparseGRM)
            }
       }
   }
	
}



set_I_mat_inR = function(modelID){
	b = as.numeric(factor(modelID, levels =  unique(modelID)))
	I_mat = Matrix::sparseMatrix(i = 1:length(b), j = b, x = rep(1, length(b)))
	I_mat = 1.0 * I_mat
	set_I_longl_mat(I_mat, b-1)
}



set_T_mat_inR = function(modelID, longlVarVec=NULL){
    if(is.null(longlVarVec)){
	stop("longlVarVec is not specified\n")
    }else{	
        b = as.numeric(factor(modelID, levels =  unique(modelID)))
        I_mat = Matrix::sparseMatrix(i = 1:length(b), j = b, x = rep(1, length(b)))
        T_mat = I_mat * longlVarVec
        set_T_longl_mat(T_mat, longlVarVec)
    }	
}



getsubGRM <- function (sparseGRMFile = NULL, sparseGRMSampleIDFile = "", relatednessCutoff,
    modelID = NULL, longlVarVec=NULL)
{
    cat("extract sparse GRM\n")
    sparseGRMLarge = Matrix:::readMM(sparseGRMFile)
    print(nnzero(sparseGRMLarge))
    cat("set elements in the sparse GRM <= ", relatednessCutoff,
        " to zero\n")
    sparseGRMLarge = Matrix:::drop0(sparseGRMLarge, tol = relatednessCutoff)
    sparseGRMLarge = sparseGRMLarge * 1
    if (!file.exists(sparseGRMSampleIDFile)) {
        stop("ERROR! sparseSigmaSampleIDFile ", sparseGRMSampleIDFile,
            " does not exist\n")
    }else{
        sparseGRMSampleID = data.frame(data.table:::fread(sparseGRMSampleIDFile,
            header = F, stringsAsFactors = FALSE, colClasses = c("character")))
        colnames(sparseGRMSampleID) = c("sampleID")
        sparseGRMSampleID$IndexGRM = seq(1, nrow(sparseGRMSampleID),
            by = 1)
        if (nrow(sparseGRMSampleID) != dim(sparseGRMLarge)[1] |
            nrow(sparseGRMSampleID) != dim(sparseGRMLarge)[2]) {
            stop("ERROR! number of samples in the sparse GRM is not the same to the number of sample IDs in the specified sparseGRMSampleIDFile ",
                sparseGRMSampleIDFile, "\n")
        }else {
            sampleInModel = NULL
            sampleInModel$IID = modelID
            sampleInModel = data.frame(sampleInModel)
            sampleInModel$IndexInModel = seq(1, length(sampleInModel$IID),
                by = 1)
	   if(!any(duplicated(sampleInModel$IID))){
            cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
            mergeID = merge(sampleInModel, sparseGRMSampleID,
                by.x = "IID", by.y = "sampleID")
            if (nrow(sampleInModel) > nrow(mergeID)) {
                stop("ERROR: ", nrow(sampleInModel) - nrow(mergeID),
                  "samples used for model fitting are not in the specified GRM\n")
            }else {
                mergeID = mergeID[with(mergeID, order(IndexInModel)),
                  ]
                indexIDofGRM = mergeID$IndexGRM
                sparseGRM = sparseGRMLarge[indexIDofGRM, indexIDofGRM]
		setupSparseGRM_new(sparseGRM)
                rm(sparseGRMLarge)
		rm(sparseGRM)
                #return(sparseGRM)
            }	    
          }else{ #if(!any(duplicated(sampleInModel$IID))){
	     cat(nrow(sampleInModel), " observations have been used to fit the glmm null model\n")	
	      mergeID = merge(sampleInModel, sparseGRMSampleID,
                by.x = "IID", by.y = "sampleID")
  	    if(length(unique(sampleInModel$IID)) > length(unique(mergeID$IID))){    	  
		stop("ERROR: ", length(unique(mergeID$IID)) - length(unique(sampleInModel$IID)), 
			  "samples used for model fitting are not in the specified GRM\n")
			    
            }else{
		mergeID = mergeID[with(mergeID, order(IndexInModel)),
                  ]
                indexIDofGRM = mergeID$IndexGRM
                sparseGRM = sparseGRMLarge[indexIDofGRM, indexIDofGRM]
		setupSparseGRM_new(sparseGRM)
		print(sparseGRM[1:10,1:10])
                rm(sparseGRMLarge)
		if(!is.null(longlVarVec)){
                  sparseVmat_longl = sparseGRM * longlVarVec
		  sparseVmat_longl_b = sparseVmat_longl + t(sparseVmat_longl)
                  addNewKat(sparseVmat_longl_b)
		  print(sparseVmat_longl_b[1:10,1:10])
                  sparseVmat_longl = t(sparseVmat_longl) * longlVarVec
                  addNewKat(sparseVmat_longl)
		  print(sparseVmat_longl[1:10,1:10])
		  rm(sparseVmat_longl)
		  rm(sparseVmat_longl_b)
          	}
		rm(sparseGRM)
                #return(sparseGRM)
	    }			    
	  }	  
       }
   }
}	


checkPerfectSep<-function(formula, data, minCovariateCount){
  X1<-model.matrix(formula,data=data)
  X_name = colnames(X1)
  X1 = as.matrix(X1[,-1])
  X_name = X_name[-1]
  colnames(X1) = X_name
  formula.frame<-model.frame(formula,data=data)
  Y = model.response(formula.frame, type = "any")
  q = length(X_name)
  colnamesDelete = c()
  for(i in 1:q){
    if (length(unique(X1[,i])) == 2){
      sumTable = table(Y, X1[,i])
      if(sum(sumTable < minCovariateCount) > 0){
        colnamesDelete = c(colnamesDelete, X_name[i])
        cat("less than ", minCovariateCount, " samples in a covariate  detected! ", X_name[i], " will be excluded in the model\n")
      }

    #if(sum(sumTable == 0) > 0){
    #    colnamesDelete = c(colnamesDelete, X_name[i])
    #    cat("perfect seperation is detected! ", X_name[i], " will be excluded in the model\n")
    #  }
    }
  }

  return(colnamesDelete)
}
