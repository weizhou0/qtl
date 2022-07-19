#' Construct a sparse GRM for a given data set 
#'
#' @param bedFile character. Path to plink file (bed) to be used for calculating the sparse GRM
#' @param bimFile character. Path to plink file (bim) to be used for calculating the sparse GRM
#' @param famFile character. Path to plink file (fam) to be used for calculating the sparse GRM
#' @param outputPrefix character. Path to the output files with prefix
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples for sparse GRM. By default, 1000
#' @param relatednessCutoff float. The threshold to treat two samples as unrelated if IsSparseKin is TRUE. By default, 0.125
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param isDiagofKinSetAsOne  logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param nThreads integer. Number of threads to be used. By default, 1 
#' @param minMAFforGRM numeric. Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01
#' @return a file ended with sampleIDs.txt that contains sample IDs for the sparse GRM and a file ended with .sparseGRM.mtx that contains the sparse GRM 
#' @export
createSparseGRM = function(bedFile = "", 
		bimFile = "",
		famFile = "",
		outputPrefix="",
                numRandomMarkerforSparseKin = 1000,
                relatednessCutoff = 0.125,
		memoryChunk = 2,
	        isDiagofKinSetAsOne = FALSE,
		nThreads = 1,
		minMAFforGRM = 0.01,
		maxMissingRateforGRM = 0.15,
		isSetGeno=TRUE,
		isWritetoFiles=TRUE
                ){

  if(nThreads > 1){
    RcppParallel:::setThreadOptions(numThreads = nThreads)
    cat(nThreads, " threads are set to be used ", "\n")
  }
  #cat("numRandomMarkerforSparseKin is ", numRandomMarkerforSparseKin, "\n")
  cat("sparse GRM will be created\n")
  setminMAFforGRM(minMAFforGRM)

  if (minMAFforGRM > 0) {
    cat("Markers in the Plink file with MAF < ", minMAFforGRM,
            " will be removed before constructing GRM\n")
  }
  setmaxMissingRateforGRM(maxMissingRateforGRM)
  if (maxMissingRateforGRM > 0){
    cat("Markers in the Plink file with missing rate > ", maxMissingRateforGRM, " will be removed before constructing GRM\n")
  }

  #  
  #famFile = paste0(plinkFile, ".fam")

  fam = data.frame(data.table:::fread(famFile, header=F, stringsAsFactors=FALSE, colClasses = c(rep("character",4), rep("numeric", 2))))

  sparseGRMSampleID = fam[,2]
  sparseGRMSampleIDFile = paste0(outputPrefix,"_relatednessCutoff_",relatednessCutoff,"_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")

  if(isWritetoFiles){
    cat("write sample IDs for the sparse GRM to ", sparseGRMSampleIDFile ,"\n")
    write.table(sparseGRMSampleID, sparseGRMSampleIDFile, quote=F, col.names=F, row.names=F)
  }

  genoSampleIndex = seq(1, nrow(fam))

  cat("isDiagofKinSetAsOne ", isDiagofKinSetAsOne, "\n")

  indicatorGenoSamplesWithPheno = rep(TRUE, nrow(fam))

  if(isSetGeno){
    setgeno(bedFile, bimFile, famFile, genoSampleIndex, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne)
  }
    freqVec = getAlleleFreqVec()
    if(minMAFforGRM > 0){
      MAFindex = which(freqVec >= minMAFforGRM & freqVec <= 1-minMAFforGRM)
      cat(length(MAFindex), " markers have MAF >= ", minMAFforGRM, "\n") 
    }else{
      MAFindex = which(freqVec > 0 & freqVec < 1)
      cat(length(MAFindex), " markers have MAF > ", minMAFforGRM, "\n") 
    } 

    cat(numRandomMarkerforSparseKin, " genetic markers are randomly selected to decide which samples are related\n")
    if(length(MAFindex) < numRandomMarkerforSparseKin){
     if(minMAFforGRM > 0){
      stop("ERROR! not enough genetic markers with MAF >= ", minMAFforGRM, " to detect which samples are related\n","Try include at least ", numRandomMarkerforSparseKin, " genetic markers with MAF >= ", minMAFforGRM, " in the plink file\n")
     }else{
      stop("ERROR! not enough genetic markers with MAF > ", minMAFforGRM, " to detect which samples are related\n","Try include at least ", numRandomMarkerforSparseKin, " genetic markers with MAF > ", minMAFforGRM, " in the plink file\n")
	

     }	
    }

    markerIndexforSparseM = sample(MAFindex, size = numRandomMarkerforSparseKin, replace=FALSE)

    cat("Start detecting related samples for the sparse GRM\n")
    ta = proc.time()
    setSubMarkerIndex(markerIndexforSparseM -1)
    tb = proc.time()
    cat("tb-ta\n")
    print(tb-ta)


    cat("Start creating sparse GRM\n")
    ta = proc.time()
    sparseMList = createSparseKinParallel(nblocks = nThreads, ncore = nThreads, relatednessCutoff)
    tb = proc.time()
    cat("tb-ta\n")
    print(tb-ta)



    cat("length(sparseMList$iIndex): ", length(sparseMList$iIndex), "\n")
    print(sparseMList$iIndex[1:102])
    cat("length(sparseMList$jIndex): ", length(sparseMList$jIndex), "\n")
    print(sparseMList$jIndex[1:102])
    cat("length(sparseMList$kinValue): ", length(sparseMList$kinValue), "\n")
    print(sparseMList$kinValue[1:102])
    sparseGRM = Matrix:::sparseMatrix(i = as.vector(sparseMList$iIndex), j = as.vector(sparseMList$jIndex), x = as.vector(sparseMList$kinValue), symmetric = TRUE)
    cat("nrow(sparseGRM): ", nrow(sparseGRM), "\n")
    cat("ncol(sparseGRM): ", ncol(sparseGRM), "\n")
    cat("ncol(sparseGRM): ", sum(sparseGRM != 0), "\n")

    tc = proc.time()
    cat("tc-tb\n")
    print(tc-tb)


  if(isWritetoFiles){
    sparseGRMFile = paste0(outputPrefix,"_relatednessCutoff_",relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx")
    cat("write sparse GRM to ", sparseGRMFile ,"\n")
    Matrix:::writeMM(sparseGRM, sparseGRMFile)
  }
  return(list(sparseGRMSampleID=sparseGRMSampleID, sparseGRM=sparseGRM))
}




#createSparseKinParallel = function(markerIndexVec, nblocks, ncore, relatednessCutoff, W, tauVecNew){
createSparseKinParallel = function(nblocks, ncore, relatednessCutoff){
  #get MAT
  #MAT = Get_MultiMarkersBySample_StdGeno_Mat(markerIndexVec)
  #MAT = Get_MultiMarkersBySample_StdGeno_Mat()
  setRelatednessCutoff(relatednessCutoff)
  Get_MultiMarkersBySample_StdGeno_Mat()
#  cat("dim(MAT) is ", dim(MAT), "\n")
  tp0 = proc.time()
#  indexVec = bigGRMPar(MAT, nblocks = nblocks, verbose = FALSE, ncore = nblocks, relatednessCutoff = relatednessCutoff)
#  indexVec = bigGRMPar_new(nblocks = nblocks, verbose = TRUE, ncore = nblocks, relatednessCutoff = relatednessCutoff)

  printComb(3)
  #indexVec = findIndiceRelatedSample()
  findIndiceRelatedSample()

  #print(indexVec)
  tp1 = proc.time()
  cat("tp1 - tp0: ", tp1-tp0, "\n")
#  cat(indexVec)
  #sparseKinList = refineKin(indexVec-1, relatednessCutoff, W, tauVecNew)
  sparseKinList = refineKin(relatednessCutoff)

#  sparseKinList$kinValue = sparseKinList$kinValue * tauVecNew[2]

  Nval = getNnomissingOut()

  sparseKinList$iIndex = c(sparseKinList$iIndex, seq(1:Nval))
  sparseKinList$jIndex = c(sparseKinList$jIndex, seq(1:Nval))
#  diagKin = getDiagOfSigma(W, tauVecNew)
#  diagKin = rep(1, Nval)
  diagKin = get_DiagofKin()
  sparseKinList$kinValue = c(sparseKinList$kinValue, diagKin)
  #sparseKinList = refineKin(indexVec, relatednessCutoff, W, tauVecNew)
  #GRMvec = refineKinPar(indexVec, relatednessCutoff = relatednessCutoff, W = W, tauVecNew = tauVecNew, nblocks = nblocks, verbose = TRUE, ncore= nblocks)
  #sparseKinList = shortenList(indexVec-1, GRMvec, relatednessCutoff, W, tauVecNew)
 tp2 = proc.time()
  cat("tp2 - tp1: ", tp2-tp1, "\n")
  return(sparseKinList)
}



getSparseSigma = function(bedFile,
                bimFile,
                famFile,
                outputPrefix="",
                sparseGRMFile=NULL,
                sparseGRMSampleIDFile="",
                numRandomMarkerforSparseKin = 1000,
                relatednessCutoff = 0.125,
                minMAFforGRM=0,
                nThreads = 1,
                isDiagofKinSetAsOne = FALSE,
                obj.glmm.null,
                W,
                tauVecNew){

  cat("sparse GRM will be used\n")
  #cat("sparseGRMFile is ", sparseGRMFile, "\n")
  #sparseGRMFile = paste0(outputPrefix, ".sparseGRM.mtx")
  if(is.null(sparseGRMFile)){
    cat("sparseGRMFile is not specified and the sparse GRM will be constructed\n")
    outputPrefix1 = paste0(outputPrefix, "_allPlinksamples")
    sparseGRMList = createSparseGRM(bedFile=bedFile, bimFile=bimFile, famFile=famFile,
    outputPrefix = outputPrefix1,
    numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
    relatednessCutoff = relatednessCutoff,
    isDiagofKinSetAsOne = isDiagofKinSetAsOne,
    nThreads = nThreads,
    minMAFforGRM = minMAFforGRM,
    isSetGeno = FALSE,
    isWritetoFiles = FALSE)
    sparseGRM = sparseGRMList$sparseGRM
    sparseGRMSampleID = data.frame(sampleID = sparseGRMList$sparseGRMSampleID)
    colnames(sparseGRMSampleID) = c("sampleID")
    rm(sparseGRMList)
    #sparseGRMFile = paste0(outputPrefix1,"_relatednessCutoff_",relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx")
    #cat("sparse GRM for all samples in the plink files is stored in ", sparseGRMFile, "\n")
    #sparseGRMSampleIDFile = paste0(outputPrefix1,"_relatednessCutoff_",relatednessCutoff,"_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
    #cat("sample IDs for the constructed sparse GRM are stored in ", sparseGRMSampleIDFile, "\n")
  }else{ # if(sparseGRMFile=="")
    cat("sparse GRM has been specified\n")
    cat("read in sparse GRM from ",sparseGRMFile,"\n")
    sparseGRMLarge = Matrix:::readMM(sparseGRMFile)
    #cat("sparseSigmaFile: ", sparseSigmaFile, "\n")
    if(sparseGRMSampleIDFile != ""){
      if(!file.exists(sparseGRMSampleIDFile)){
        stop("ERROR! sparseSigmaSampleIDFile ", sparseGRMSampleIDFile, " does not exsit\n")
      }else{
        sparseGRMSampleID = data.frame(data.table:::fread(sparseGRMSampleIDFile, header=F, stringsAsFactors=FALSE, colClasses=c("character")))
        colnames(sparseGRMSampleID) = c("sampleID")
      }
    }else{#end of if(sparseSigmaSampleIDFile != "")
      stop("ERROR! sparseSigmaSampleIDFile is not specified\n")
    }
  #}

      sparseGRMSampleID$IndexGRM = c(1:nrow(sparseGRMSampleID))
        cat("length(sparseGRMSampleID$IndexGRM): ", length(sparseGRMSampleID$IndexGRM), "\n")
        cat("nrow(sparseGRMSampleID): ", nrow(sparseGRMSampleID), "\n")
      sampleInModel = NULL
      sampleInModel$IID = obj.glmm.null$sampleID
      sampleInModel = data.frame(sampleInModel)
      sampleInModel$IndexInModel = seq(1,length(sampleInModel$IID), by=1)
      cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
      mergeID = merge(sampleInModel, sparseGRMSampleID, by.x="IID", by.y = "sampleID")
      mergeID = mergeID[with(mergeID, order(IndexInModel)), ]
      print(dim(mergeID))
      print(head(mergeID))
      indexIDofGRM=mergeID$IndexGRM
      #cat("indexIDofGRM = ", indexIDofGRM, "\n")
      #cat("Subset sparse GRM to be ", indexIDofSigma," by ", indexIDofSigma, "\n")
      sparseGRM = sparseGRMLarge[indexIDofGRM, indexIDofGRM]
      rm(sparseGRMLarge)
}
  #sparseGRMFile = paste0(outputPrefix,"_relatednessCutoff_",relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.subset.mtx")
  #cat("write sparse GRM to ", sparseGRMFile ,"\n")
  #Matrix:::writeMM(sparseGRM, sparseGRMFile)
  Nval = length(W)
  sparseSigma = sparseGRM * tauVecNew[2]
  #diag(sparseSigma) = getDiagOfSigma(W, tauVecNew)
  diag(sparseSigma) = diag(sparseSigma) + (1/W)*tauVecNew[1]
  #sparseSigmaFile = paste0(outputPrefix, "_relatednessCutoff_",relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseSigma.mtx")
  #cat("write sparse Sigma to ", sparseSigmaFile ,"\n")
  #Matrix:::writeMM(sparseSigma, sparseSigmaFile)
  return(sparseSigma)
}



