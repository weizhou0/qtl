#Fits the null glmm
glmmkin.ai_PCG_Rcpp_multiV_NB = function(bedFile, bimFile, famFile, Xorig, isCovariateOffset, fit0, tau=c(0,0), fixtau = c(0,0), maxiter =20, tol = 0.02, verbose = TRUE, nrun=30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, indicatorGenoSamplesWithPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne, isLowMemLOCO, covarianceIdxMat = NULL, isStoreSigma = FALSE, useSparseGRMtoFitNULL = TRUE, useGRMtoFitNULL = TRUE) {
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
  if(is.null(subPheno$IndexGeno)){
	subSampleInGeno = subPheno$IndexPheno 
  }	  
  if(verbose){
    print("Start reading genotype plink file here")
  }

  #print("subSampleInGeno")
  #print(subSampleInGeno)

  set_dup_sample_index(as.numeric(factor(subPheno$IID, levels =  unique(subPheno$IID))))

  if(bedFile != "" & useGRMtoFitNULL){

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

  family = fit0$family
  y = fit0$y
  n = length(y)
  X = model.matrix(fit0)
  offset = fit0$offset
  
  print("offset")
  print(offset[1:10])

  if(is.null(offset)){
    offset = rep(0, n)
  }


  wts <- model.weights(fit0)
  if (is.null(wts)) wts <- rep(1, n)
  w = fit0$prior.weights
  eta = fit0$linear.predictors
  print("eta")
  print(eta[1:10])
  zz = eta + fit0$residuals - offset
  wz = fit0$weights
  invwt <- 1/(wz + 1e-04)
 
  print("invwt")
  print(invwt[1:10])

  y = fit0$y

  mu = fit0$fitted.values
  mu.eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu)/mu.eta
  print("Y")
  print(Y[1:10])
  zz <- eta + fit0$residuals - offset
  print("zz")
  print(zz[1:10])
  sqrtW = mu.eta/sqrt(fit0$family$variance(mu))
  W = sqrtW^2

  datanew = cbind(X, zz)
  datanew = as.data.frame(datanew)
  colnames(datanew)[ncol(datanew)] = "y"
  newformula = as.formula(paste0("y ~ ", strsplit(as.character(fit0$formula), split="~")[[3]]))

  print("head(datanew)")
  print(head(datanew))
  fit0_gaussian = glm(newformula, data = datanew,
                  family = gaussian(link = "identity"))

  family <- NegBin()
  fam <- family


  for (i in seq_len(maxiter)) {
    W = sqrtW^2
    if(verbose) cat("\nIteration ", i, tau, ":\n")
      datanew = cbind(X, zz)    
      datanew = as.data.frame(datanew)
      colnames(datanew)[ncol(datanew)] = "y"
      newformula = as.formula(paste0("y ~ ", strsplit(as.character(fit0$formula), split="~")[[3]]))
      fit0_gaussian = glm(newformula, data = datanew,
                  family = gaussian(link = "identity"), weights=invwt)



       system.time(modglmm_gaussian <- glmmkin.ai_PCG_Rcpp_multiV(bedFile, bimFile, famFile, Xorig, isCovariateOffset,
                fit0_gaussian, tau = tau, fixtau = fixtau, maxiter = maxiter,
                tol = tol, verbose = TRUE, nrun = nrun, tolPCG = tolPCG,
                maxiterPCG = maxiterPCG, subPheno = subPheno, indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno,
                obj.noK = obj.noK, out.transform = out.transform,
                tauInit = tauInit, memoryChunk = memoryChunk,
                LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec,
                chromosomeEndIndexVec = chromosomeEndIndexVec,
                traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform,
                isDiagofKinSetAsOne = isDiagofKinSetAsOne,
                isLowMemLOCO = isLowMemLOCO, covarianceIdxMat = covarianceIdxMat, isStoreSigma = isStoreSigma, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL, useGRMtoFitNULL = useGRMtoFitNULL)) 
      etaold = eta
      #eta = (modglmm_gaussian$fitted.values)*(modglmm_gaussian$varWeights) + offset
      eta = (modglmm_gaussian$fitted.values) + offset

      print("modglmm_gaussian$fitted.value")
      print(modglmm_gaussian$fitted.value[1:10])

      #cat("sum((eta - etaold)^2)/sum(eta^2) ", sum((eta - etaold)^2)/sum(eta^2), "\n")	

      #if (i > 1 & sum((eta - etaold)^2) < 1e-05 * sum(eta^2)) break
      mu <- fam$linkinv(eta)
      mu.eta.val <- fam$mu.eta(eta)
      mu.eta.val <- ifelse(mu.eta.val == 0, 1e-04, mu.eta.val) 
      #print("mu")
      #print(mu[1:20]) 
      varmu <- fam$variance(mu)
      #print("varmu")
      #print(varmu[1:20])
      varmu <- ifelse(varmu == 0, 1e-04, varmu)
      zz <- eta + (y - mu)/mu.eta.val - offset
      wz <- w * mu.eta.val^2/varmu
      wz <- ifelse(wz == 0, 1e-04, wz)
      invwt <- 1/wz
      #invwt <- wz

	print("zz")
       print(zz[1:10])
       print("mu")
       print(mu[1:10])
	print("invwt")
       print(invwt[1:10])

      if(i > 1){
	th0 = th
      }else{
	th0 = th = 0
      }	      
      th <- suppressWarnings( MASS::theta.ml(y=y, mu=mu, n=sum(wts), weights=wts, limit=10, trace=FALSE) )
      if (is.null(th)) th <- fam$theta

     if(i > 1){
      cat("th ", th, "\n")
      cat("th0 ", th0, "\n")
      cat("abs(th - th0)/(abs(th) + abs(th0) + tol) ", abs(th - th0)/(abs(th) + abs(th0) + tol), "\n")	
      if(abs(th - th0)/(abs(th) + abs(th0) + tol) < tol) break
     }
     fam <- NegBin(theta = th)

  }

  if(verbose) cat("\nFinal " ,modglmm_gaussian$theta, ":\n")

  print("fam$theta")
  print(fam$theta)

  print("eta[1:10]")
  print(eta[1:10])
  print("mu[1:10]")
  print(mu[1:10])
      mu <- fam$linkinv(eta)
      mu.eta.val <- fam$mu.eta(eta)
      mu.eta.val <- ifelse(mu.eta.val == 0, 1e-04, mu.eta.val)
      varmu <- fam$variance(mu)
      varmu <- ifelse(varmu == 0, 1e-04, varmu)
       


  converged = ifelse(i < maxiter, TRUE, FALSE)
  #res = y - mu
  #mu2 = fit0$family$variance(mu)
  #traitType = "count_nb"

  print(names(modglmm_gaussian))
  print(modglmm_gaussian$coefficients)
  print(modglmm_gaussian$alpha)

  glmmResult = modglmm_gaussian


  #print("names(modglmm_gaussian)")
  #print(names(modglmm_gaussian))
  #print("modglmm_gaussian$y[1:10]")
  #print(modglmm_gaussian$y[1:10])

  #mu2 = varmu
  #mu2 = rep(((1/(modglmm_gaussian$theta[1]))),length(modglmm_gaussian$y))
  #mu2 = mu2 * invwt
  #print("mu2 HEREHREHREHREHRERHE")
  #print(mu2[1:10])
  #print(length(mu2))
  #print("y[1:10]")
  #print(y[1:10])


  #mu2 = mu + (mu^2/fam$theta)
  #res = y - mu
  #y = modglmm_gaussian$y

  #if(!isCovariateOffset){
  #  obj.noK = ScoreTest_NULL_Model(mu, mu2, y, X)
  #  glmmResult$X = X
  #   glmmResult$obj.noK = obj.noK
    #glmmResult = list(theta=modglmm_gaussian$theta, coefficients=modglmm_gaussian$coefficients, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = modglmm_gaussian$sampleID, obj.noK=obj.noK, y = y, X = modglmm_gaussian$X, traitType=traitType, isCovariateOffset = modglmm_gaussian$isCovariateOffset)
  #}else{
  #  obj.noK = ScoreTest_NULL_Model(mu, mu2, y, Xorig)
  #  glmmResult$X = Xorig
  #  glmmResult$obj.noK = obj.noK
    #glmmResult = list(theta=modglmm_gaussian$theta, coefficients=modglmm_gaussian$coefficients, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,sampleID = modglmm_gaussian$sampleID, obj.noK=obj.noK, y = y, X = modglmm_gaussian$Xorig, traitType=traitType, isCovariateOffset = modglmm_gaussian$isCovariateOffset)
  #}
  #glmmResult$traitType = traitType
  #glmmResult$theta[1] = 1
  glmmResult$theta_overdispersion = fam$theta
  #glmmResult$varWeights = invwt 
  #glmmResult$varWeights = NULL
  #glmmResult$traitType = "quantitative"
  #if(isLowMemLOCO & LOCO){
  #  glmmResult$chromosomeStartIndexVec = modglmm_gaussian$chromosomeStartIndexVec
  #  glmmResult$chromosomeEndIndexVec = modglmm_gaussian$chromosomeEndIndexVec
  #}
  return(glmmResult)
}
