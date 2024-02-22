setSparseSigma = function(sparseSigmaFile){

  Check_File_Exist(sparseSigmaFile, "sparseSigmaFile")
  sparseSigma = Matrix:::readMM(sparseSigmaFile)
  locations = rbind(sparseSigma@i, sparseSigma@j)
  values = sparseSigma@x
  nSubj = dim(sparseSigma)[1]
  sigmaMatListR = list(locations = locations,
                     values = values,
                     nSubj = nSubj)
  return(sigmaMatListR)	
}	


setSparseSigma_new = function(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, sampleIDInModel, tauVec, W, traitType){

  Check_File_Exist(sparseGRMFile, "sparseGRMFile")
  Check_File_Exist(sparseGRMSampleIDFile, "sparseGRMSampleIDFile")

  sparseGRM = Matrix:::readMM(sparseGRMFile)
  if(class(sparseGRM) == "nsTMatrix"){
     sparseGRM = sparseGRM * 1
  }

  #print(sparseGRM[1:20,1:20])
  sparseGRMSampleID = data.frame(data.table:::fread(sparseGRMSampleIDFile, header=F, stringsAsFactors=FALSE, colClasses=c("character")))
  colnames(sparseGRMSampleID) = "sampleID"
  sparseGRMSampleID$IndexGRM = c(1:nrow(sparseGRMSampleID))
  cat("length(sparseGRMSampleID$IndexGRM): ", length(sparseGRMSampleID$IndexGRM), "\n")
  cat("nrow(sparseGRMSampleID): ", nrow(sparseGRMSampleID), "\n")
  sampleInModel = NULL
  sampleInModel$IID = sampleIDInModel
  #rm(sampleIDInModel)
  sampleInModel = data.frame(sampleInModel)
  sampleInModel$IndexInModel = seq(1,length(sampleInModel$IID), by=1)
  cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
  mergeID = merge(sampleInModel, sparseGRMSampleID, by.x="IID", by.y = "sampleID")
  mergeID = mergeID[with(mergeID, order(IndexInModel)), ]
  print(dim(mergeID))
  print(head(mergeID))
  indexIDofGRM=mergeID$IndexGRM
  sparseGRM = sparseGRM[indexIDofGRM, indexIDofGRM]
  if(length(indexIDofGRM) < nrow(sampleInModel)){
    stop(nrow(sampleInModel)-length(indexIDofGRM), " samples were not found in the sparse GRM\n")
  }else{
    print("Subsetting GRM")
  }
  sumSpGRM=summary(sparseGRM)

  print(head(sumSpGRM))
  
  removeIndex = which(sumSpGRM[,3] < relatednessCutoff)
  if(length(removeIndex) > 0){
	cat("Removing ", length(removeIndex), " elements in the sparse GRM < ", relatednessCutoff, ".\n")  
  	sumSpGRM = sumSpGRM[-removeIndex,,drop=F]
  }

  sumSpGRM[,3] = sumSpGRM[,3] * tauVec[2]
  #sparseSigma = sparseGRM * tauVec[2]
  if(traitType == "binary" | traitType == "count"){
	sumSpGRM[which(sumSpGRM[,1] == sumSpGRM[,2]), 3] = sumSpGRM[which(sumSpGRM[,1] == sumSpGRM[,2]), 3] + 1/W

   #diag(sparseSigma) = W + diag(sparseSigma)
  }else if(traitType == "quantitative"){
	sumSpGRM[which(sumSpGRM[,1] == sumSpGRM[,2]), 3] = sumSpGRM[which(sumSpGRM[,1] == sumSpGRM[,2]), 3] + tauVec[1]

  }
  #nSubj = dim(sparseGRM)[1]

  #locations = t(sumSpGRM[,c(1,2)])
  #print(dim(locations))
  #sigmaMatListR = list(locations = t(sumSpGRM[,c(1,2)]),
  #                   values = sumSpGRM[,3],
  #                   nSubj = nSubj)
  #return(sigmaMatListR)
  return(sumSpGRM)
}



                        #DosageCutoff_for_UltraRarePresence, 
			#method_to_CollapseUltraRare,
##working
SAIGE.Region = function(mu,
			OutputFile,
                        MACCutoff_to_CollapseUltraRare,
			groupFile, 
			annolist, 
			maxMAFlist,
			minMAFlist,
		        markers_per_chunk_in_groupTest,	
			genoType, 
			markerInfo,
			traitType,
			phenotype_name_vec,
			isImputation,
			isCondition,
			weight_cond,
			groups_per_chunk,
			r.corr,
			isOverWriteOutput,
			is_single_in_groupTest,
			BetaDist_weight_mat, 
			is_equal_weight_in_groupTest,
			is_output_markerList_in_groupTest,
			is_SKATO, 
			chrom, 
			is_fastTest,
			pval_cutoff_for_fastTest,
			is_output_moreDetails){
  cat("maxMAFlist ", maxMAFlist, "\n")
  cat("minMAFlist ", minMAFlist, "\n")

  OutputFileIndex = NULL	
  if(is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index")
   
  outList = checkOutputFile(OutputFile, OutputFileIndex, "Region", 1, isOverWriteOutput) # Check 'Util.R'

  cat("traitType ", traitType, "\n")

  indexChunk = outList$indexChunk


  Start = outList$Start
  End = outList$End	

  cat("Start ", Start, "\n")
  cat("End ", End, "\n")


  if(End)
  {
    message = paste0("The analysis has been completed in earlier analysis. Results are saved in '", OutputFile, "'. ",
                     "If you want to change parameters and restart the analysis, please use another 'OutputFile'.")
    return(message)
  }

  isappend=FALSE
 if(!Start){ 
  isappend=TRUE
 }

  n = length(mu) #sample size 


  if(r.corr==0){
    out.method = SKAT:::SKAT_Check_Method(method="optimal.adj", r.corr=0)
    method=out.method$method
    r.corr=out.method$r.corr
    cat("SKAT-O test will be performed. P-values for BURDEN and SKAT tests will also be output\n")
    regionTestType = "SKAT-O"
    is_single_in_groupTest = TRUE
  }else if(r.corr == 1){	  
    method = NULL
    cat("BURDEN test will be performed\n")
    regionTestType = "BURDEN"

    #output the result from Rcpp
    cat("isappend ", isappend, "\n")
    for(itt in 1:length(traitType)){
      if(length(traitType) == 1){
	OutputFile_itt = OutputFile
      }else{
	OutputFile_itt = paste0(OutputFile, "_", phenotype_name_vec[itt])
      }
        assign_g_outputFilePrefix(OutputFile_itt)
        isOpenOutFile = openOutfile(traitType[itt],isappend)
        if(!isOpenOutFile){
                stop("Output file ", OutputFile_itt, " can't be opened\n")
        }
    }

  }else{
    stop("r.corr needs to be either 1 (BURDEN test) or 0 (SKAT-O test)\n")
  }	  

  assign_g_outputFilePrefix0(OutputFile)

  if(is_single_in_groupTest){
  cat("is_single_in_groupTest = TRUE. Single-variant assoc tests results will be output\n")


    for(itt in 1:length(traitType)){
      if(length(traitType) == 1){
	OutputFile_itt = OutputFile
      }else{
        OutputFile_itt = paste0(OutputFile, "_", phenotype_name_vec[itt])
      } 
       assign_g_outputFilePrefix(OutputFile_itt)
       isOpenOutFile_singleinGroup = openOutfile_singleinGroup(traitType[itt], isImputation, isappend, is_output_moreDetails)
       if(!isOpenOutFile_singleinGroup){
           stop("Output file ", OutputFile_itt, ".singleAssoc.txt can't be opened\n")
       }
    }
  }else{
      cat("is_single_in_groupTest = FALSE. Single-variant assoc tests results will not be output\n")
  }

  cat("Number of phenotypes to test:\t", length(traitType), "\n")
  ##check group file
  region_list = checkGroupFile(groupFile)

  nRegions = region_list$nRegions
  is_weight_included = region_list$is_weight_included
  numberofWeightlists = region_list$numberofWeightlists
  nameofWeightlists = region_list$nameofWeightlists
  weightlist = nameofWeightlists
  #if(is_equal_weight_in_groupTest & is_weight_included){
  #  stop("is_equal_weight_in_groupTest = TRUE but weights are found in the group file.\n")
  #}

  if(sum(BetaDist_weight_mat) > 0){
	for(i in 1:nrow(BetaDist_weight_mat)){
		weightlist = c(weightlist, paste0("Beta_", BetaDist_weight_mat[i,1], "_", BetaDist_weight_mat[i,2]))
	}
  }

  if(is_equal_weight_in_groupTest){
     cat("Equal weights are used in the group test\n")
     weightlist = c(weightlist, "Equal_Weights")
  }

  



  if(is_weight_included){
    nline_per_gene = 2 + numberofWeightlists
  }else{
    nline_per_gene = 2
  }	  

  gf = file(groupFile, "r")

  skipline = indexChunk*nline_per_gene
  if(indexChunk > 0 & indexChunk < nRegions){
    for(k in 1:skipline){
        marker_group_line_temp = readLines(gf, n = 1) 
        rm(marker_group_line_temp)
    }
  }

  if(regionTestType != "BURDEN"){

  	cat("length(traitType) ", length(traitType), "\n")
	cat("markers_per_chunk_in_groupTest ", markers_per_chunk_in_groupTest, "\n")
  	P1Mat = matrix(0, markers_per_chunk_in_groupTest * length(traitType), n)
  	P2Mat = matrix(0, n, markers_per_chunk_in_groupTest * length(traitType))
  }else{
	P1Mat = matrix(0, 1, 1)
	P2Mat = matrix(0, 1, 1)  
  }	  

  chrom1 = "FakeCHR";

  gc()
  num_region = 0
  mth = 0

  numberRegionsInChunk = 0
  #cat("indexChunk ", indexChunk, "\n")
  #cat("nRegions ", nRegions, "\n")
  pval.Region.all = NULL
  pval.Region.all.list = list()
  
  OutList.all = NULL
  Output_MarkerList.all = NULL
  cth_chunk_to_output=1

  i = indexChunk+1


  cat("nRegions " , nRegions, "\n")


  while(i <= nRegions){
  #for(i in (indexChunk+1):nRegions){
   if(mth ==  numberRegionsInChunk){
      if(i + groups_per_chunk > nRegions){
  	      nregions_ro_read = nRegions - i + 1	      
      }else{
	      nregions_ro_read = groups_per_chunk
      }
      nlinetoread = nregions_ro_read * nline_per_gene
      marker_group_line = readLines(gf, n = nlinetoread)
      RegionList = SAIGE.getRegionList_new(marker_group_line, nline_per_gene, annolist, markerInfo, chrom, nameofWeightlists)
      #print(RegionList)
      cat("Read in ", nregions_ro_read, " region(s) from the group file.\n")
      mth = 0
      #numberRegionsInChunk = length(RegionList)
      numberRegionsInChunk = nregions_ro_read
    }

   mth = mth + 1
   if(!is.null(RegionList)){
    pval.Region = NULL
    region = RegionList[[mth]]

    #print("region")
    #print(names(region))


    annolistsub = region$annoVec 
    regionName = names(RegionList)[mth]
    i = i + 1

    if(!is.null(region$SNP) & length(annolistsub) > 0){

      SNP = region$SNP
      if(genoType == "vcf"){
        SNPlist = paste(c(regionName, SNP), collapse = "\t") 
        if(length(SNP) == 1){
		 fakem = strsplit(SNP, split=":")[[1]]
		 fakemb = paste(c(fakem[1], as.numeric(fakem[2
])+2, "N", "N"), collapse=":")
		 SNPlisttemp = paste(c(SNPlist, fakemb), collapse = "\t") 
		set_iterator_inVcf(SNPlisttemp, chrom1, 1, 250000000)
	}else{	
        	set_iterator_inVcf(SNPlist, chrom1, 1, 250000000)
	}

	isVcfEnd =  check_Vcf_end()
    	if(!isVcfEnd){
		region$genoIndex = rep("0", length(SNP))
		region$genoIndex_prev = rep("0", length(SNP))
    	}else{
        	warning("No markers in region ", regionName, " are found in the VCF file")
		next
    	}	
      }

      
      #if(is_equal_weight_in_groupTest){
	#WEIGHT = rep(1, length(SNP))
      #}else{	      
      #  WEIGHT = as.numeric(region$WEIGHT)
      #}
	if(is_weight_included){
		WEIGHT = matrix(as.numeric(unlist(region$WEIGHT)),    # Convert to numeric matrix
		                  ncol = length(region$WEIGHT))
		
	}else{
		WEIGHT = matrix(c(0,0), ncol=1)
	}

      annoIndicatorMat = region$annoIndicatorMat
      #cat("annoIndicatorMat ", annoIndicatorMat, "\n")


      #chrom = region$chrom
      print(paste0("Analyzing Region ", regionName, " (",i-1,"/",nRegions,")."))
      #tp1 = proc.time()
      #gc()

      if(!is_fastTest){
        set_flagSparseGRM_cur_SAIGE_org()
      }else{
        set_flagSparseGRM_cur_SAIGE(FALSE)
      }
      outList = mainRegionInCPP(genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat, weightlist, maxMAFlist, minMAFlist,  OutputFile, traitType, n, P1Mat, P2Mat, regionTestType, isImputation, BetaDist_weight_mat, WEIGHT, weight_cond, is_equal_weight_in_groupTest, is_single_in_groupTest, is_output_markerList_in_groupTest, annolistsub, regionName, is_fastTest, is_output_moreDetails)	

      if(regionTestType == "BURDEN" & is_fastTest){
	if(!is.null(outList$iswriteOutput)){
	  if(!(outList$iswriteOutput)){	
            set_flagSparseGRM_cur_SAIGE(TRUE)
	    outList = mainRegionInCPP(genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat, weightlist,maxMAFlist, minMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType, isImputation, BetaDist_weight_mat, WEIGHT, weight_cond, is_equal_weight_in_groupTest, is_single_in_groupTest, is_output_markerList_in_groupTest, annolistsub, regionName, is_fastTest, is_output_moreDetails)
          }
	}
      }	
#print("time_mainRegionInCPP")
#print(time_mainRegionInCPP)
     if(!is_fastTest){
      rm(region)
     }  
      #rm(genoIndex)
      #gc()
    #tb0 = proc.time()

#print("outList")
#print(outList)
#print(names(outList))
#print(outList)

q_multTrait = length(outList$pvalVec)/length(traitType)

for(tt in 1:length(traitType)){
	tt_start = q_multTrait*(tt-1)+1
	tt_end = q_multTrait*(tt)
	pval.Region = NULL
	pval.Region.all = NULL
      if(is_single_in_groupTest){
      #OutList = as.data.frame(outList$OUT_DF)

        #cat("length(outList$pvalVec) ", length(outList$pvalVec), "\n")
	#print(outList$pvalVec[tt_start:tt_end])
	
        noNAIndices = which(!is.na(outList$pvalVec[tt_start:tt_end]))
	#print("length(outList$pvalVec)")

	#print("WEIGHT")
	#print(WEIGHT)

        if(sum(WEIGHT) > 0){
          AnnoWeights = c(WEIGHT, rep(1, outList$numofUR))
        }
      }

      annoMAFIndicatorMat = outList$annoMAFIndicatorMat

	#print("tt_start")
	#print(tt_start)
	#print("tt_end")
	#print(tt_end)


      if((sum(outList$NumUltraRare_GroupVec) + sum(outList$NumRare_GroupVec)) > 0){

        if(regionTestType != "BURDEN"){
 	  #ta0 = proc.time()
          if(traitType[tt] == "binary" ){	     
            gyVec= outList$gyVec[tt_start:tt_end][noNAIndices]
          }

          if(isCondition){#check later
            VarMatAdjCond = outList$VarMatAdjCond[tt_start:tt_end,][noNAIndices,noNAIndices]
            TstatAdjCond = outList$TstatAdjCond[tt_start:tt_end][noNAIndices]
	    print(dim(outList$G1tilde_P_G2tilde_Weighted_Mat))
	    cat("tt_start ", tt_start, "\n")
	    cat("tt_end ", tt_end, "\n")
	
            G1tilde_P_G2tilde_Weighted_Mat = outList$G1tilde_P_G2tilde_Weighted_Mat[tt_start:tt_end,,drop=F][noNAIndices,,drop=F]

            weightMat_G2_G2 = outList$G2_Weight_cond %*% t(outList$G2_Weight_cond)
          }	  
        
       ### Get annotation maf indicators

       ### 3. Adjust for saddlepoint approximation
       #notNAindice = which(!is.na(outList$TstatVec_flip))
       #print("notNAindice")
       #print(notNAindice)
       #cat("length(outList$TstatVec_flip) ", length(outList$TstatVec_flip), "\n")
       StatVec = outList$TstatVec_flip[tt_start:tt_end][noNAIndices]
	#noNAIndices_start_end = c(tt_start:tt_end)[noNAIndices]

	#cat("length(StatVec) ", length(StatVec), "\n")
	#print(dim(outList$VarMat))

  	#ncolVarMat = ncol(outList$VarMat)/length(traitType)
 	#tt_start_sub = ncolVarMat*(tt-1)+1
 	#tt_end_sub = ncolVarMat*(tt)

  	length_varchunk = ncol(outList$VarMat)/length(traitType)
        startvarMat = length_varchunk*(tt-1)+1
        endvarMat = length_varchunk*tt


  	#cat("ncolVarMat ", ncolVarMat, "\n")
  	#cat("tt ", tt, "\n")
	#cat("tt_start_sub ", tt_start_sub, " tt_end_sub ", tt_end_sub, "\n")
	#print(dim(outList$VarMat[, tt_start_sub:tt_end_sub]))
       VarSVec = diag(outList$VarMat[,startvarMat:endvarMat])
	#cat("length(VarSVec) ", length(VarSVec), "\n")
       VarSVec = VarSVec[!is.na(VarSVec)]
	#cat("length(VarSVec) ", length(VarSVec), "\n")
       pvalVec0 = outList$pvalVec[tt_start:tt_end]
	#cat("length(pvalVec0) ", length(pvalVec0), "\n")
       adjPVec = pvalVec0[!is.na(pvalVec0)]
		

       #varTestedIndices = which(apply(annoMAFIndicatorMat, 1, isContainValue, val=1))
       #print("varTestedIndices")
       #print(varTestedIndices)
       #varTestedIndices = which(rowSums(annoMAFIndicatorMat) > 0)
       #annoMAFIndicatorMat = annoMAFIndicatorMat[varTestedIndices, , drop=F]
       #MAFVec = outList$MAFVec[varTestedIndices]
       annoMAFIndicatorMat = annoMAFIndicatorMat[noNAIndices, , drop=F]
       #MAFVec = outList$MAFVec[1:q_multTrait][noNAIndices]
       MAFVec = outList$MAFVec[noNAIndices]
       #AnnoWeights = dbeta(MAFVec,1,25)
       #if(sum(WEIGHT) > 0){
 #	 AnnoWeights = AnnoWeights[noNAIndices] 
 #      }else{
 #        AnnoWeights = dbeta(MAFVec,1,25)
 #      }

 #       weightMat = AnnoWeights %*% t(AnnoWeights)

	if(sum(WEIGHT) > 0){
		AnnoWeightsMat = WEIGHT[noNAIndices,]
	}else{
		AnnoWeightsMat = NULL
	}

  	if(sum(BetaDist_weight_mat) > 0){
        	for(beta in 1:nrow(BetaDist_weight_mat)){
			AnnoWeightsMat = cbind(AnnoWeightsMat, dbeta(MAFVec,BetaDist_weight_mat[beta,1], BetaDist_weight_mat[beta,2]))			
        	}
  	}
	if(is_equal_weight_in_groupTest){
		AnnoWeightsMat = cbind(AnnoWeightsMat, rep(1, nrow(AnnoWeightsMat)))

	}


	#print("AnnoWeightsMat")
	#print(AnnoWeightsMat)

	weightMat_G1_G2_Mat = list()
	wStatMat = NULL
	wadjVarSMat_list = list()
	wStatVec_cond_Mat = NULL
	wadjVarSMat_cond_list = list()


       if(isCondition){  #check  
    	 #weightMat_G1_G2 = AnnoWeights %*% t(outList$G2_Weight_cond)
	for(ia in 1:ncol(AnnoWeightsMat)){
    	 	weightMat_G1_G2_Mat[[ia]] = AnnoWeightsMat[,ia] %*% t(outList$G2_Weight_cond)
	 }
       }


	#length_varchunk = ncol(outList$VarMat)/length(traitType)
	#startvarMat = length_varchunk*(tt-1)+1
	#endvarMat = length_varchunk*tt

	#cat("startvarMat ", startvarMat, "\n")
	#cat("endvarMat ", endvarMat, "\n")


	for(ia in 1:ncol(AnnoWeightsMat)){	
                AnnoWeights = AnnoWeightsMat[,ia]
                wStatVec = StatVec * AnnoWeights
		weightMat = AnnoWeights %*% t(AnnoWeights)
                wadjVarSMat = (outList$VarMat)[,startvarMat:endvarMat] * weightMat

                wStatMat = cbind(wStatMat, wStatVec)
                wadjVarSMat_list[[ia]] = wadjVarSMat
		if(isCondition){
	  		wStatVec_cond = wStatVec - outList$TstatAdjCond[tt_start:tt_end][noNAIndices]
			print(dim(outList$VarMatAdjCond))
			print(dim(wadjVarSMat))
			print(length(wStatVec))
			print(length(outList$TstatAdjCond))
			print(tt_start)
			print(tt_end)
          		wadjVarSMat_cond = wadjVarSMat - outList$VarMatAdjCond[tt_start:tt_end,][noNAIndices,noNAIndices]
			wStatVec_cond_Mat = cbind(wStatVec_cond_Mat, wStatVec_cond)
			wadjVarSMat_cond_list[[ia]] = wadjVarSMat_cond	
     		}
	}	
      #print("outList$VarMat")
      #print(outList$VarMat)


    #gc()

    #print("q_maf_for_anno")
    #print(outList$q_maf_for_annoMat)

#	print("AnnoWeightsMat")
#	print(dim(AnnoWeightsMat))

#	print("annoMAFIndicatorMat")
#	print(dim(annoMAFIndicatorMat))
    annoMAFIndVec = c()
  for(r in 1:length(weightlist)){
    for(j in 1:length(annolistsub)){
	AnnoName = annolistsub[j]
	#maxMAF0 = outList$q_maf_for_annoVec[j]
	isPolyRegion = TRUE
	for(m in 1:length(maxMAFlist)){
	  maxMAF0 = outList$q_maf_for_annoMat[j, m]
	  jm = (j-1)*(length(maxMAFlist)) + m
    	  jmr = (r-1)*(length(annolistsub)) * (length(maxMAFlist)) + (j-1) * (length(maxMAFlist)) + m
	  maxMAFName = maxMAFlist[m]
	  minMAFName = minMAFlist[m]

	  if(maxMAF0 != 1){
          	pval.Region.temp = NULL
          	annoMAFIndVec.temp = NULL
          }

	#cat("jmr ", jmr, "\n")

	#print(weightlist)
	#print(dim(wadjVarSMat_list[[r]]))
	#print(dim(wStatMat))

       if(maxMAF0 != 1){
		
	    #if(m <= maxMAF0){
	       tempPos = which(annoMAFIndicatorMat[,jmr] == 1)
		#print("length(tempPos)")
		#print(length(tempPos))
	

	       if(length(tempPos) > 0){
	        isPolyRegion = TRUE
		annoMAFIndVec = c(annoMAFIndVec, jmr)
		annoMAFIndVec.temp = c(annoMAFIndVec.temp, jmr)
		weightName = weightlist[r]
		wadjVarSMat = wadjVarSMat_list[[r]] 
		wStatVec = wStatMat[,r] 
		AnnoWeights = AnnoWeightsMat[,r]

		
		#print("AnnoWeightsMat b")
		#print(AnnoWeightsMat)

		Phi = wadjVarSMat[tempPos, tempPos, drop=F]
		Score = wStatVec[tempPos]
		p.new = adjPVec[tempPos]

		#if(traitType == "binary" | traitType == "count"){
			if(traitType[tt] == "binary"){
				g.sum = outList$genoSumMat[,jmr]
				q.sum<-sum(gyVec[tempPos] * AnnoWeights[tempPos])
				mu.a = mu

				re_phi = get_newPhi_scaleFactor(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType)
		        	Phi = re_phi$val
                	}

			#print("AnnoWeights[tempPos]")
			#print(AnnoWeights[tempPos])
			Pvalue_ACATV = get_CCT_pvalue(p.new, AnnoWeights[tempPos])	
			if(is_SKATO){	
				groupOutList = get_SKAT_pvalue(Score, Phi, r.corr, regionTestType)
			#groupOutList = get_SKAT_pvalue_Burden_SKAT_ACATV(Score, Phi, p.new, AnnoWeights[tempPos])
			#print("weightName")	
			#print(weightName)
			#print("regionName")
			#print(regionName)
			#print("groupOutList")
			#print(groupOutList)
				Pvalue_ACATO = get_CCT_pvalue(c(Pvalue_ACATV, groupOutList$Pvalue_SKATO))
			}else{
				groupOutList = get_SKAT_pvalue_Burden_SKAT(Score, Phi)
				Pvalue_ACATO = get_CCT_pvalue(c(Pvalue_ACATV, groupOutList$Pvalue_SKAT, groupOutList$Pvalue_Burden))
				groupOutList$Pvalue_SKATO = NA
			}

			resultDF = data.frame(Region = regionName,
                                                    Group = AnnoName,
						    min_MAF = minMAFName,
                                                    max_MAF = maxMAFName,
						    Weight = weightName,
						    Pvalue = Pvalue_ACATO,
                                                    Pvalue_ACATV = Pvalue_ACATV,
						    Pvalue_SKATO = groupOutList$Pvalue_SKATO,
                                                    Pvalue_Burden = groupOutList$Pvalue_Burden,
                                                    Pvalue_SKAT = groupOutList$Pvalue_SKAT,
                                                    BETA_Burden = groupOutList$BETA_Burden,
                                                    SE_Burden = groupOutList$SE_Burden)
	      	     if(isCondition){


			if(traitType[tt] == "binary"){
				#print("length(outList$scalefactor_G2_cond)")
				#print(length(outList$scalefactor_G2_cond))
				G1tilde_P_G2tilde_Mat_scaled = t(t((outList$G1tilde_P_G2tilde_Weighted_Mat[tt_start:tt_end,][tempPos,,drop=F]) * sqrt(as.vector(re_phi$scaleFactor))) * sqrt(as.vector(outList$scalefactor_G2_cond)))
#t(t(b * sqrt(a1)) * sqrt(a2))
		        	adjCondTemp = G1tilde_P_G2tilde_Mat_scaled %*% outList$VarInvMat_G2_cond_scaled	
				VarMatAdjCond = adjCondTemp %*% t(G1tilde_P_G2tilde_Mat_scaled)
				TstatAdjCond = adjCondTemp %*% (outList$Tstat_G2_cond * outList$G2_Weight_cond)
				Phi_cond = re_phi$val - diag(VarMatAdjCond)
				Score_cond = Score - TstatAdjCond
			
			}else{
				wStatVec_cond = wStatVec_cond_Mat[,r]
				#print("length(wStatVec_cond)")
				#print(length(wStatVec_cond))
				#print(length(tempPos))
				wadjVarSMat_cond = wadjVarSMat_cond_list[[r]]
				#print(dim(wadjVarSMat_cond))
				Score_cond = wStatVec_cond[tempPos]
				Phi_cond = wadjVarSMat_cond[tempPos, tempPos]
			}
			P_cond = pchisq(Score_cond^2/diag(Phi_cond), df=1, lower.tail=F)
			#groupOutList_cond = get_SKAT_pvalue(Score_cond, Phi_cond, r.corr, regionTestType)
			#groupOutList_cond = get_SKAT_pvalue_Burden_SKAT_ACATV(Score_cond, Phi_cond, P_cond, AnnoWeights)

			Pvalue_ACATV_cond = get_CCT_pvalue(P_cond, AnnoWeights[tempPos])

			if(is_SKATO){
				groupOutList_cond = get_SKAT_pvalue(Score_cond, Phi_cond, r.corr, regionTestType)
				Pvalue_ACATO_cond = get_CCT_pvalue(c(Pvalue_ACATV_cond, groupOutList_cond$Pvalue_SKATO))	
			}else{
				groupOutList_cond = get_SKAT_pvalue_Burden_SKAT(Score_cond, Phi_cond)
				Pvalue_ACATO_cond = get_CCT_pvalue(c(Pvalue_ACATV_cond, groupOutList_cond$Pvalue_SKAT, groupOutList_cond$Pvalue_Burden))
				                                groupOutList_cond$Pvalue_SKATO = NA
			}

		#resultDF$Pvalue_cond = groupOutList_cond$Pvalue_SKATO
			resultDF$Pvalue_cond = Pvalue_ACATO_cond
			resultDF$Pvalue_ACATV_cond = Pvalue_ACATV_cond
			resultDF$Pvalue_SKATO_cond = groupOutList_cond$Pvalue_SKATO
			resultDF$Pvalue_Burden_cond = groupOutList_cond$Pvalue_Burden
			resultDF$Pvalue_SKAT_cond = groupOutList_cond$Pvalue_SKAT
			resultDF$BETA_Burden_cond = groupOutList_cond$BETA_Burden
			resultDF$SE_Burden_cond = groupOutList_cond$SE_Burden
	     	}#if(isCondition){
		pval.Region.temp = rbind(pval.Region.temp, resultDF)	
		pval.Region = rbind.data.frame(pval.Region, resultDF)

	   }else{#if(length(tempPos) > 0){
		isPolyRegion = FALSE
	   }	
		

	#}else{ #if(m <= maxMAF0){
	
		   
	}else{ #if(m <= maxMAF0){
		
	   if(isPolyRegion){
		annoMAFIndVec = c(annoMAFIndVec, annoMAFIndVec.temp)
		resultDF = pval.Region.temp
		print("pval.Region.temp")
		print(pval.Region.temp)
		#resultDF$Region = regionName
		#resultDF$Group = AnnoName
		resultDF$min_MAF = minMAFName	
		resultDF$max_MAF = maxMAFName	
		#resultDF$Weight = weightName	
		pval.Region = rbind.data.frame(pval.Region, resultDF)
			
		#pval.Region_temp = rbind.data.frame(pval.Region_temp, resultDF)
	   }

	}
    }#for(m in 1:length(maxMAFlist)){
  }#for(j in 1:length(annolist)){
}#for(r in 1:length(weightlist)){


gc()
#}
#print("pval.Region")
#print(pval.Region)
#print("outList$NumRare_GroupVec")
#print(outList$NumRare_GroupVec)
#print("annoMAFIndVec")
#print(annoMAFIndVec)


    if(length(annoMAFIndVec) > 0){
      pval.Region$MAC = outList$MAC_GroupVec[annoMAFIndVec]
      if(traitType[tt] == "binary"){
    	pval.Region$MAC_case = outList$MACCase_GroupVec[annoMAFIndVec]
    	pval.Region$MAC_control = outList$MACCtrl_GroupVec[annoMAFIndVec]
      }
      pval.Region$Number_rare = outList$NumRare_GroupVec[annoMAFIndVec]
      pval.Region$Number_ultra_rare = outList$NumUltraRare_GroupVec[annoMAFIndVec]
    }

if(length(annolistsub) > 1 | length(maxMAFlist) > 1 | length(weightlist) > 1){
   cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden)
   if(regionTestType != "BURDEN"){
     cctpval = get_CCT_pvalue(pval.Region$Pvalue)
     #cctpval = get_CCT_pvalue(c(pval.Region$Pvalue_Burden, pval.Region$Pvalue_SKAT, pval.Region$Pvalue_ACATV))
     cctpval_SKATO = get_CCT_pvalue(pval.Region$Pvalue_SKATO)
     cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT)
     cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden)
     cctpval_ACATV = get_CCT_pvalue(pval.Region$Pvalue_ACATV)
     #cctVec = c(regionName, "Cauchy", NA, cctpval, cctpval_Burden, cctpval_SKAT, NA, NA)
     cctVec = c(regionName, "Cauchy", NA, NA, NA, cctpval, cctpval_ACATV, cctpval_SKATO, cctpval_Burden, cctpval_SKAT,  NA, NA)
   }else{
	cctVec = c(regionName, "Cauchy", NA, NA, NA, cctpval_Burden, NA, NA)
   } 	   
   if(isCondition){
   	cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)

   	if(regionTestType != "BURDEN"){
     	  cctpval = get_CCT_pvalue(pval.Region$Pvalue_cond)
	  #cctpval = get_CCT_pvalue(c(pval.Region$Pvalue_Burden_cond, pval.Region$Pvalue_SKAT_cond, pval.Region$Pvalue_ACATV_cond))
	  cctpval_SKATO = get_CCT_pvalue(pval.Region$Pvalue_SKATO_cond)
	  cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT_cond)
	  cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)
	  cctpval_ACATV = get_CCT_pvalue(pval.Region$Pvalue_ACATV_cond)
	  cctVec = c(cctVec, cctpval, cctpval_ACATV, cctpval_SKATO, cctpval_Burden, cctpval_SKAT, NA, NA)
	}else{
	  cctVec = c(cctVec, cctpval_Burden, NA, NA)

	}	
   }
   cctVec = c(cctVec, NA)
   if(traitType[tt] == "binary"){
     cctVec = c(cctVec, NA, NA)
   }
   cctVec = c(cctVec, NA, NA)

   pval.Region = rbind(pval.Region, cctVec)
}else{
   cctpval = 1	
}
#ta1 = proc.time()
#print("ta0 - tb0")
#print(ta0 - tb0)
#print("ta1 - ta0")
#print(ta1 - ta0)




.f = function() {


if(is_fastTest){
  if(cctpval < pval_cutoff_for_fastTest){
#if(cctpval < 0.1 & is_fastTest){
      pval.Region = NULL
      cat("Non-fast test is performed\n")
      set_flagSparseGRM_cur_SAIGE(TRUE)
      outList = mainRegionInCPP(genoType, region$genoIndex_prev, region$genoIndex, annoIndicatorMat, weightlist,maxMAFlist, minMAFlist, OutputFile, traitType, n, P1Mat, P2Mat, regionTestType, isImputation, BetaDist_weight_mat, WEIGHT, weight_cond, is_equal_weight_in_groupTest, is_single_in_groupTest, is_output_markerList_in_groupTest, annolistsub, regionName, is_fastTest, is_output_moreDetails)
  if(is_single_in_groupTest){
      #OutList = as.data.frame(outList$OUT_DF)
      noNAIndices = which(!is.na(outList$pvalVec))
      #print(noNAIndices)
      annoMAFIndicatorMat = outList$annoMAFIndicatorMat
      if(sum(WEIGHT) > 0){
        AnnoWeights = c(WEIGHT, rep(1, outList$numofUR))
      }
  }

       if(traitType[tt] == "binary"){	     
         outList$gyVec = outList$gyVec[noNAIndices]
       }

       if(isCondition){
         outList$VarMatAdjCond = outList$VarMatAdjCond[noNAIndices,noNAIndices]
         outList$TstatAdjCond = outList$TstatAdjCond[noNAIndices]
         outList$G1tilde_P_G2tilde_Weighted_Mat = outList$G1tilde_P_G2tilde_Weighted_Mat[noNAIndices,,drop=F]
         weightMat_G2_G2 = outList$G2_Weight_cond %*% t(outList$G2_Weight_cond)
       }	  
        
       ### Get annotation maf indicators

       ### 3. Adjust for saddlepoint approximation
       #notNAindice = which(!is.na(outList$TstatVec_flip))
       #print("notNAindice")
       #print(notNAindice)
       StatVec = outList$TstatVec_flip[noNAIndices]
       VarSVec = diag(outList$VarMat)
       VarSVec = VarSVec[!is.na(VarSVec)]
       adjPVec = outList$pvalVec[!is.na(outList$pvalVec)]
		
       #varTestedIndices = which(apply(annoMAFIndicatorMat, 1, isContainValue, val=1))
       #print("varTestedIndices")
       #print(varTestedIndices)
       #varTestedIndices = which(rowSums(annoMAFIndicatorMat) > 0)
       #annoMAFIndicatorMat = annoMAFIndicatorMat[varTestedIndices, , drop=F]
       #MAFVec = outList$MAFVec[varTestedIndices]
       annoMAFIndicatorMat = annoMAFIndicatorMat[noNAIndices, , drop=F]
       MAFVec = outList$MAFVec[noNAIndices]
       #AnnoWeights = dbeta(MAFVec,1,25)

	print("AnnoWeights")
	print(AnnoWeights)

       #if(sum(WEIGHT) > 0){
 #	 #AnnoWeights = AnnoWeights[varTestedIndices] 
 #	 AnnoWeights = AnnoWeights[noNAIndices] 
 #      }else{
 #        AnnoWeights = dbeta(MAFVec,1,25)
 #      }
       #weightMat = AnnoWeights %*% t(AnnoWeights)


       #if(isCondition){    
    #	 weightMat_G1_G2 = AnnoWeights %*% t(outList$G2_Weight_cond)
    #   }	
 
 #     wStatVec = StatVec * AnnoWeights

#      wadjVarSMat = outList$VarMat * weightMat

#	if(isCondition){
#	  wStatVec_cond = wStatVec - outList$TstatAdjCond
#          wadjVarSMat_cond = wadjVarSMat - outList$VarMatAdjCond	
#      }

 	 if(isCondition){
        	for(ia in 1:ncol(AnnoWeightsMat)){
                	weightMat_G1_G2_Mat[[ia]] = AnnoWeightsMat[,ia] %*% t(outList$G2_Weight_cond)
         	}
       	}
        for(ia in 1:ncol(AnnoWeightsMat)){
                AnnoWeights = AnnoWeightsMat[,ia]
                wStatVec = StatVec * AnnoWeights
                weightMat = AnnoWeights %*% t(AnnoWeights)
                wadjVarSMat = outList$VarMat * weightMat

                wStatMat = cbind(wStatMat, wStatVec)
                wadjVarSMat_list[[ia]] = wadjVarSMat
                if(isCondition){
                        wStatVec_cond = wStatVec - outList$TstatAdjCond
                        wadjVarSMat_cond = wadjVarSMat - outList$VarMatAdjCond
                        wStatVec_cond_Mat = cbind(wStatVec_cond_Mat, wStatVec_cond)
                        wadjVarSMat_cond_list[[ia]] = wadjVarSMat_cond
                }
        }




    #gc()


    annoMAFIndVec = c()
    for(j in 1:length(annolistsub)){
	AnnoName = annolistsub[j]
	#maxMAF0 = outList$q_maf_for_annoVec[j]
	isPolyRegion = TRUE
	for(m in 1:length(maxMAFlist)){
		maxMAF0 = outList$q_maf_for_annoMat[j,m]
		jm = (j-1)*(length(maxMAFlist)) + m
		minMAFName = minMAFlist[m]
		maxMAFName = maxMAFlist[m]
		if(maxMAF0 != 1){
			pval.Region.temp = NULL	
		}


	for(r in 1:length(weightlist)){

         jmr = (j-1)*(length(maxMAFlist)) * (length(weightlist)) + (m-1) * (length(weightlist)) + r

	    #if(m <= maxMAF0){
	    if(maxMAF0 != 1){
	       tempPos = which(annoMAFIndicatorMat[,jmr] == 1)
	       if(length(tempPos) > 0){
	       isPolyRegion = TRUE
	        annoMAFIndVec = c(annoMAFIndVec, jmr)
		wadjVarSMat = wadjVarSMat_list[[r]]
		Phi = wadjVarSMat[tempPos, tempPos, drop=F]
		wStatVec = wStatMat[, r]
		Score = wStatVec[tempPos]


         	weightName = weightlist[r]
                AnnoWeights = AnnoWeightsMat[,r]

		if(traitType[tt] == "binary" | traitType[tt] == "count"){
			p.new = adjPVec[tempPos]
			g.sum = outList$genoSumMat[,jmr]
			q.sum<-sum(outList$gyVec[tempPos] * AnnoWeights[tempPos])
			mu.a = mu
			re_phi = get_newPhi_scaleFactor(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType)
		        Phi = re_phi$val
                }
                        Pvalue_ACATV = get_CCT_pvalue(p.new, AnnoWeights[tempPos])
                        Pvalue_ACATO = get_CCT_pvalue(c(Pvalue_ACATV, groupOutList$Pvalue_SKATO))
                        resultDF = data.frame(Region = regionName,
                                                    Group = AnnoName,
                                                    min_MAF = minMAFName,
                                                    max_MAF = maxMAFName,
                                                    Weight = weightName,
                                                    Pvalue = Pvalue_ACATO,
                                                    Pvalue_ACATV = Pvalue_ACATV,
                                                    Pvalue_SKATO = groupOutList$Pvalue_SKATO,
                                                    Pvalue_Burden = groupOutList$Pvalue_Burden,
                                                    Pvalue_SKAT = groupOutList$Pvalue_SKAT,
                                                    BETA_Burden = groupOutList$BETA_Burden,
                                                    SE_Burden = groupOutList$SE_Burden)


		#groupOutList = get_SKAT_pvalue(Score, Phi, r.corr, regionTestType)

		#resultDF = data.frame(Region = regionName,
                #                                    Group = AnnoName,
                #                                    max_MAF = maxMAFName,
                #                                    Pvalue = groupOutList$Pvalue_SKATO,
                #                                    Pvalue_Burden = groupOutList$Pvalue_Burden,
                #                                    Pvalue_SKAT = groupOutList$Pvalue_SKAT,
                #                                    BETA_Burden = groupOutList$BETA_Burden,
                #                                    SE_Burden = groupOutList$SE_Burden)
	      if(isCondition){
		if(traitType[tt] == "binary"){
			G1tilde_P_G2tilde_Mat_scaled = t(t((outList$G1tilde_P_G2tilde_Weighted_Mat[tempPos,,drop=F]) * sqrt(as.vector(re_phi$scaleFactor))) * sqrt(as.vector(outList$scalefactor_G2_cond)))
#t(t(b * sqrt(a1)) * sqrt(a2))
		        adjCondTemp = G1tilde_P_G2tilde_Mat_scaled %*% outList$VarInvMat_G2_cond_scaled	
			VarMatAdjCond = adjCondTemp %*% t(G1tilde_P_G2tilde_Mat_scaled)
			TstatAdjCond = adjCondTemp %*% (outList$Tstat_G2_cond * outList$G2_Weight_cond)
			Phi_cond = re_phi$val - diag(VarMatAdjCond)
			Score_cond = Score - TstatAdjCond
			
		}else{
			wStatVec_cond = wStatVec_cond_Mat[,r]
			Score_cond = wStatVec_cond[tempPos]
			wadjVarSMat_cond = wadjVarSMat_cond_list[[r]]
			Phi_cond = wadjVarSMat_cond[tempPos, tempPos]
		}
		 Pvalue_ACATV_cond = get_CCT_pvalue(P_cond, AnnoWeights[tempPos])

		if(is_SKATO){

			groupOutList_cond = get_SKAT_pvalue(Score_cond, Phi_cond, r.corr, regionTestType)
                	Pvalue_ACATO_cond = get_CCT_pvalue(c(Pvalue_ACATV_cond, groupOutList_cond$Pvalue_SKATO))
		}else{
			groupOutList_cond = get_SKAT_pvalue_Burden_SKAT(Score_cond, Phi_cond)
			Pvalue_ACATO_cond = get_CCT_pvalue(c(Pvalue_ACATV_cond, groupOutList_cond$Pvalue_SKAT, groupOutList_cond$Pvalue_Burden))
	                                                                groupOutList_cond$Pvalue_SKATO = NA
		}	


                #resultDF$Pvalue_cond = groupOutList_cond$Pvalue_SKATO
                        resultDF$Pvalue_cond = Pvalue_ACATO_cond
                        resultDF$Pvalue_ACATV_cond = Pvalue_ACATV_cond
                        resultDF$Pvalue_SKATO_cond = groupOutList_cond$Pvalue_SKATO
                        resultDF$Pvalue_Burden_cond = groupOutList_cond$Pvalue_Burden
                        resultDF$Pvalue_SKAT_cond = groupOutList_cond$Pvalue_SKAT
                        resultDF$BETA_Burden_cond = groupOutList_cond$BETA_Burden
                        resultDF$SE_Burden_cond = groupOutList_cond$SE_Burden




		#resultDF$Pvalue_cond = groupOutList_cond$Pvalue_SKATO
		#resultDF$Pvalue_Burden_cond = groupOutList_cond$Pvalue_Burden
		#resultDF$Pvalue_SKAT_cond = groupOutList_cond$Pvalue_SKAT
		#resultDF$BETA_Burden_cond = groupOutList_cond$BETA_Burden
		#resultDF$SE_Burden_cond = groupOutList_cond$SE_Burden
	     }#if(isCondition){
		pval.Region = rbind.data.frame(pval.Region, resultDF)
		pval.Region.temp = rbind.data.frame(pval.Region.temp, resultDF)
		print("pval.Region.temp a")
		print(pval.Region.temp)
	   }else{#if(length(tempPos) > 0){
		isPolyRegion = FALSE
	   }	

	}else{ #if(m <= maxMAF0){
	   if(isPolyRegion){
		annoMAFIndVec = c(annoMAFIndVec, jm)   
		#resultDF$Region = regionName
		#resultDF$Group = AnnoName
		#resultDF$min_MAF = minMAFName	
		#resultDF$max_MAF = maxMAFName
	  	pval.Region.temp$Region = regionName
                pval.Region.temp$Group = AnnoName
                pval.Region.temp$min_MAF = minMAFName
                pval.Region.temp$max_MAF = maxMAFName
		#print("pval.Region.temp b")
		#print(pval.Region.temp)
		#pval.Region = rbind.data.frame(pval.Region, resultDF)
		pval.Region = rbind.data.frame(pval.Region, pval.Region.temp)
	   }			   
	}#if(m <= maxMAF0){

        }#for(r in 1:length(weightlist)){


    }#for(m in 1:length(maxMAFlist)){
}#for(j in 1:length(annolist)){


gc()
#}



#if(regionTestType != "BURDEN"){

 if(length(annoMAFIndVec) > 0){
      pval.Region$MAC = outList$MAC_GroupVec[annoMAFIndVec]
      if(traitType == "binary"){
    	pval.Region$MAC_case = outList$MACCase_GroupVec[annoMAFIndVec]
    	pval.Region$MAC_control = outList$MACCtrl_GroupVec[annoMAFIndVec]
      }
      pval.Region$Number_rare = outList$NumRare_GroupVec[annoMAFIndVec]
      pval.Region$Number_ultra_rare = outList$NumUltraRare_GroupVec[annoMAFIndVec]
    }

 if(length(annolistsub) > 1 | length(maxMAFlist) > 1){

   cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden)
   if(regionTestType != "BURDEN"){
     cctpval = get_CCT_pvalue(pval.Region$Pvalue)
     cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT)
     cctpval_SKATO = get_CCT_pvalue(pval.Region$Pvalue_SKATO)
     cctpval_ACATV = get_CCT_pvalue(pval.Region$Pvalue_ACATV)
     cctVec = c(regionName, "Cauchy", NA, NA, cctpval, cctpval_ACATV, cctpval_SKATO, cctpval_Burden, cctpval_SKAT,  NA, NA)	
     #cctVec = c(regionName, "Cauchy", NA, cctpval, cctpval_Burden, cctpval_SKAT, NA, NA)
   }else{
	cctVec = c(regionName, "Cauchy", NA, cctpval_Burden, NA, NA)
   } 
   if(isCondition){
cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)

        if(regionTestType != "BURDEN"){
          cctpval = get_CCT_pvalue(pval.Region$Pvalue_cond)
          #cctpval = get_CCT_pvalue(c(pval.Region$Pvalue_Burden_cond, pval.Region$Pvalue_SKAT_cond, pval.Region$Pvalue_ACATV_cond))
          cctpval_SKATO = get_CCT_pvalue(pval.Region$Pvalue_SKATO_cond)
          cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT_cond)
          cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)
          cctpval_ACATV = get_CCT_pvalue(pval.Region$Pvalue_ACATV_cond)
          cctVec = c(cctVec, cctpval, cctpval_ACATV, cctpval_SKATO, cctpval_Burden, cctpval_SKAT, NA, NA)
        }else{
          cctVec = c(cctVec, cctpval_Burden, NA, NA)

        }
}

	   
  # if(isCondition){
  # 	cctpval_Burden = get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)

  # 	if(regionTestType != "BURDEN"){
#	  cctpval = get_CCT_pvalue(pval.Region$Pvalue_cond)
#	  cctpval_SKAT = get_CCT_pvalue(pval.Region$Pvalue_SKAT_cond)
#	  cctVec = c(cctVec, cctpval, cctpval_Burden, cctpval_SKAT, NA, NA)
#	}else{
#	  cctVec = c(cctVec, cctpval_Burden, NA, NA)

#	}	
#   }
   cctVec = c(cctVec, NA)
   if(traitType == "binary"){
     cctVec = c(cctVec, NA, NA)
   }
   cctVec = c(cctVec, NA, NA)

   pval.Region = rbind(pval.Region, cctVec)
}


}else{ #cctpval < 0.05)
   copy_singleInGroup()
}
}#if(is_fastTest){

}#.f = function() {

}#if(regionTestType != "BURDEN"){


  Output_MarkerList = NULL
  if(is_output_markerList_in_groupTest){
    for(j in 1:length(annolistsub)){
        AnnoName = annolistsub[j]
        for(m in 1:length(maxMAFlist)){
                jm = (j-1)*(length(maxMAFlist)) + m
                minMAFName = minMAFlist[m]
                maxMAFName = maxMAFlist[m]
                tempPos = which(outList$annoMAFIndicatorMat[,jm] > 0)
                marker_rare_pos = which(outList$markerIndcatorVec == 1)
                marker_ultrarare_pos = which(outList$markerIndcatorVec == 2)
                if(length(tempPos) > 0){
                        if(length(marker_rare_pos) > 0){
				markerind_b = which(marker_rare_pos %in% tempPos)
				if(length(markerind_b) > 0){
					markerind = marker_rare_pos[markerind_b]
					SNPlist_rare = paste(SNP[markerind], collapse=",")
				}else{
					SNPlist_rare = ""
				}	
                        }else{
                                SNPlist_rare = ""
                        }
                        if(length(marker_ultrarare_pos) > 0){
				markerind_UR_b = which(marker_ultrarare_pos %in% tempPos)
				if(length(markerind_UR_b) > 0){
					markerindUR = marker_ultrarare_pos[markerind_UR_b]
					SNPlist_Ultra_rare = paste(SNP[markerindUR], collapse=",")
				}else{
					SNPlist_Ultra_rare = ""
				}	
                        }else{
                                SNPlist_Ultra_rare = ""
                        }
                        Output_MarkerList = rbind(Output_MarkerList, c(regionName, AnnoName, maxMAFName, SNPlist_rare, SNPlist_Ultra_rare))
                }

        }
    }
  }

#ta2 = proc.time()	
#print("ta2 - ta1")
#print(ta2 - ta1)
   if(is_output_markerList_in_groupTest){
	colnames(Output_MarkerList) = c("Region", "Group", "max_MAF", "Rare_Variants", "Ultra_Rare_Variants")
	Output_MarkerList.all = rbind(Output_MarkerList.all, Output_MarkerList)   
   }else{
	Output_MarkerList.all = NULL
   }	   

  indexChunk = i

  #Start = (i==1)

#if(tt == length(traitType)){
  Start = (cth_chunk_to_output==1)
  End = (i==nRegions)
#}

#i = i + 1
  AnalysisType = "Region"
  nEachChunk = 1

if(regionTestType != "BURDEN"){  
    pval.Region.all = rbind(pval.Region.all, pval.Region)
}


if(is_output_markerList_in_groupTest){
     rm(Output_MarkerList)
 } 	   

if(tt == length(traitType)){
	rm(outList)
}
rm(pval.Region)
if(regionTestType != "BURDEN"){
     rm(resultDF)
} 
gc()
  
}#if(length(noNAIndices) > 0){ 
pval.Region.all.list[[tt]] = pval.Region.all
pval.Region.all = NULL
}##for(tt in 1:length(traitType)){


}else{#if(!is.null(region)){
    cat(regionName, " is empty.\n")
  }


# output
if(mth ==  numberRegionsInChunk){



  message1 = "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)

  #cat("indexChunk ", indexChunk, "\n")



  message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 = "Have completed the analyses of all chunks."
  #n1 = length(Output)
  #n2 = length(OutputFile)
  cat("write to output\n")
  #cat("n1 is ", n1, "\n")
  #cat("n2 is ", n2, "\n")

for(tt in 1:length(traitType)){
  pval.Region.all = pval.Region.all.list[[tt]]
  if(length(traitType) == 1){
	OutputFile_tt = OutputFile
  }else{
  	OutputFile_tt = paste0(OutputFile, "_", phenotype_name_vec[tt])
  }
  if(regionTestType != "BURDEN"){
      if(Start){
        if(!is.null(pval.Region.all)){
          fwrite(pval.Region.all, OutputFile_tt, quote = F, sep = "\t", append = F, col.names = T, row.names = F, na="NA")
        }
      }else{
        if(!is.null(pval.Region.all)){
          fwrite(pval.Region.all, OutputFile_tt, quote = F, sep = "\t", append = T, col.names = F, row.names = F, na="NA")
        }
        #write.table(Output, OutputFile, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
      }
  }
      if(is_output_markerList_in_groupTest){
        if(Start){
        if(!is.null(Output_MarkerList.all)){
          fwrite(Output_MarkerList.all, paste0(OutputFile, ".markerList.txt"), quote = F, sep = "\t", append = F, col.names = T, row.names = F, na="NA")
        }
      }else{
        if(!is.null(Output_MarkerList.all)){
          fwrite(Output_MarkerList.all, paste0(OutputFile, ".markerList.txt"), quote = F, sep = "\t", append = T, col.names = F, row.names = F, na="NA")
        }
        #write.table(Output, OutputFile, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
      }
    }


pval.Region.all = NULL
OutList.all = NULL
Output_MarkerList.all = NULL
cth_chunk_to_output = cth_chunk_to_output + 1
gc()
} #for(tt in 1:length(traitType)){

#if(FALSE){
  #print("write Output 2")
  if(Start){
    write.table(c(message1, message2, message3), OutputFileIndex,
                quote = F, sep = "\t", append = F, col.names = F, row.names = F)
  }
  if(!End){
  write.table(message4, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
  }

  if(End){
    write.table(message5, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
  }
#}#if(FALSE)



}#if(mth ==  numberRegionsInChunk){


}else{#if(!is.null(RegionList)){
     cat("The chunk is empty\n")	   
     #mth = 0
     mth = numberRegionsInChunk
     i = i + numberRegionsInChunk
     pval.Region = NULL
   }
	   
}

  message = paste0("Analysis done! The results have been saved to '", OutputFile,"' and '",
                   paste0(OutputFile, ".markerInfo"),"'.")

  message = "Analysis done!"

for(itt in 1:length(traitType)){

  if(length(traitType) == 1){
       OutputFile_itt = OutputFile
  }else{
       OutputFile_itt = paste0(OutputFile, "_", phenotype_name_vec[itt])
  } 


  message = paste0(message, " The set-based tests results have been saved to '", OutputFile_itt, "'.")
  if(is_output_markerList_in_groupTest){
	message = paste0(message, " The marker lists have been saved to '", OutputFile_itt, ".markerList.txt'.")
  }	  
  if(is_single_in_groupTest){
  	message = paste0(message, " The single-variant association tests results have been saved to '", OutputFile_itt, ".singleAssoc.txt'.")
  } 

 }#for(itt in 1:length(traitType)){

if(length(traitType) > 1){
       assign_g_outputFilePrefix(OutputFile)
       removeOutfile_inGroup()
}

 if(is_single_in_groupTest){
    for(itt in 1:length(traitType)){
     if(length(traitType) == 1){
	OutputFile_itt = OutputFile
     }else{
       OutputFile_itt = paste0(OutputFile, "_", phenotype_name_vec[itt])
     }	
       assign_g_outputFilePrefix(OutputFile_itt)
       removeOutfile_singleinGroup_temp()
    } 
 }


  return(message)

}



SAIGE.getRegionList_new = function(marker_group_line,
			 nline_per_gene,	   
                         annoVec, #c("lof","lof;missense"
                         markerInfo, 
			 chrom="", 
			 nameofWeightlists = NULL)
{
  chrom_nochr = gsub("CHR", "", chrom, ignore.case = T)
  # read group file
  ngroup<-length(marker_group_line)/nline_per_gene
  #cat("ngroup is ", ngroup, "\n")
  RegionData = NULL
  geneList = c() 
  for(i in 1:ngroup){
	  marker_group_line_list = strsplit(marker_group_line[1+(i-1)*nline_per_gene], split="[\ \t]+")[[1]]
          gene=marker_group_line_list[1]
          var=marker_group_line_list[3:length(marker_group_line_list)]
          marker_group_line_list_anno = strsplit(marker_group_line[2+(i-1)*nline_per_gene], split="[\ \t]+")[[1]]
	  anno=marker_group_line_list_anno[3:length(marker_group_line_list_anno)]
  	  geneData = cbind(rep(gene, length(var)), var, anno) 
	  if(nline_per_gene > 2){
	      for(j in 1:(nline_per_gene-2)){
              	marker_group_line_list_weight = strsplit(marker_group_line[2+j+(i-1)*nline_per_gene], split="[\ \t]+")[[1]]
              	weight = as.numeric(marker_group_line_list_weight[3:length(marker_group_line_list_weight)])
		geneData = cbind(geneData, weight)
	      }	
	      
	      
	      RegionData = rbind(RegionData, geneData)
	  }else if(nline_per_gene == 2){
	      RegionData = rbind(RegionData, cbind(rep(gene, length(var)), var, anno))
	  }

	  if(gene %in% geneList){
		stop(gene, " is duplicated in the group File\n")
          }else{		  
	  	geneList = c(geneList, gene)
	  }
  }
    if(nline_per_gene == 2){
    	colnames(RegionData) = c("REGION", "SNP", "ANNO")
    }else if(nline_per_gene == 3){
	colnames(RegionData) = c("REGION", "SNP", "ANNO", "WEIGHT")
    }else{
	if(is.null(nameofWeightlists)){
		stop("nameofWeightlists is empty\n")
	}else if(length(nameofWeightlists) != (nline_per_gene-2)){
		stop("The length of nameofWeightlists is not equal to nline_per_gene-2\n")
	}else{
		colnames(RegionData) = c("REGION", "SNP", "ANNO", nameofWeightlists)
	}
    }
    RegionData = as.data.frame(RegionData)
    setDT(RegionData)
    uRegion0 = unique(RegionData$REGION)    

    if(chrom != "" & is.null(markerInfo)){
      RegionData[, c("chr") := tstrsplit(RegionData$SNP, ":")[[1]] ]
      setkey(RegionData, "chr")
      RegionData = RegionData[chr == chrom | chr == chrom_nochr]
      RegionData[,chr:=NULL]
    }

if(nrow(RegionData) != 0){
if(!is.null(markerInfo)){
  setkey(RegionData, "SNP")
  RegionData = merge(RegionData, markerInfo, by.x = "SNP", by.y = "ID", all.x = T, sort = F) 

  if(!is.null(markerInfo$ID2)){
        RegionData = merge(RegionData, markerInfo, by.x = "SNP", by.y = "ID2", all.x = T, sort = F)
	#SNP REGION     ANNO CHROM.x POS.x genoIndex2.x  ID2 genoIndex_prev.x CHROM.y POS.y    ID genoIndex2.y genoIndex_prev.y
        setnames(RegionData, "genoIndex.x", "genoIndex")
        setnames(RegionData, "genoIndex.y", "genoIndex2")
        #setnames(RegionData, "genoIndex.y", "genoIndex")
        setnames(RegionData, "CHROM.x", "CHROM")
        setnames(RegionData, "CHROM.y", "CHROM2")
        setnames(RegionData, "POS.x", "POS")
        setnames(RegionData, "POS.y", "POS2")

	if(!is.null(RegionData$genoIndex_prev.y)){
		setnames(RegionData, "genoIndex_prev.x", "genoIndex_prev")
		setnames(RegionData, "genoIndex_prev.y", "genoIndex_prev2")
		#markerInfo[,genoIndex_prev.y:=NULL]

	}
	posNA = which(is.na(RegionData$genoIndex) & !is.na(RegionData$genoIndex2))

	if(length(posNA) != 0){
		RegionData$genoIndex[posNA] = RegionData$genoIndex2[posNA]
		RegionData$CHROM[posNA] = RegionData$CHROM2[posNA]
		RegionData$POS[posNA] = RegionData$POS2[posNA]
		if(!is.null(RegionData$genoIndex_prev)){
			RegionData$genoIndex_prev[posNA] = RegionData$genoIndex_prev2[posNA]			
		}
	}
	#RegionData$genoIndex2 = NULL	
  }
  posNA = which(is.na(RegionData$genoIndex))
  
    if(length(posNA) != 0){
      RegionData = RegionData[which(!is.na(RegionData$genoIndex))]
      cat(length(posNA)," markers in 'RegionFile' are not in 'GenoFile'.\n")
    }
   setorderv(RegionData, col=c("CHROM", "POS"))
}

}


if(nrow(RegionData) != 0){

  
  #HeaderInRegionData = colnames(RegionData)
  HeaderInRegionData = unique(RegionData$ANNO)
  RegionAnnoHeaderList = list()
  if(length(annoVec) == 0){
        stop("At least one annotation is required\n")
  }
  for(q in 1:length(annoVec)){
        RegionAnnoHeaderList[[q]] = strsplit(annoVec[q],";")[[1]]
  }

  RegionList = list()
  uRegion = unique(RegionData$REGION)
  #RegionData = as.data.frame(RegionData)


  for(r in uRegion0){
             #print(paste0("Analyzing region ",r,"...."))
    #print(RegionData$REGION)
    #print(r)
    #which(as.vector(RegionData$REGION) == r)
   
    posSNP = which(RegionData$REGION == r)
    if(length(posSNP) > 0){
       SNP = RegionData$SNP[posSNP]

       if(nline_per_gene == 3){	
         WEIGHT = RegionData[posSNP, 4, drop=F]
       }else if(nline_per_gene > 3){
	 WEIGHT = RegionData[posSNP, 4:(nline_per_gene+1), drop=F]
       }


      if(any(duplicated(SNP)))
        stop("Please check RegionFile: in region ", r,": duplicated SNPs exist.")

      if(!is.null(markerInfo)){
        genoIndex = as.numeric(RegionData$genoIndex[posSNP])
        if(!is.null(RegionData$genoIndex_prev)){
		genoIndex_prev = as.numeric(RegionData$genoIndex_prev[posSNP])	
        }
        chrom = RegionData$CHROM[posSNP]
      #uchrom = unique(chrom)
      #if(length(uchrom) != 1)
      #  stop("In region ",r,", markers are from multiple chromosomes.")
       }

    annoIndicatorMat = matrix(0, nrow=length(posSNP), ncol=length(annoVec))
    annoVecNew = c()
        for(q in 1:length(annoVec)){
        indiceVec = which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]])
        if(length(indiceVec) > 0){
                annoVecNew = c(annoVecNew, annoVec[q])
                annoIndicatorMat[indiceVec, q] = 1
                #annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
        }
    }

  RegionAnnoHeaderListNew = list()
  if(length(annoVecNew) == 0){
        warning("No markers are found for at least one annotation, so region ", r, " is skipped\n")
        #stop("At least one annotation is required\n")
  }else{
   if(length(annoVecNew) < length(annoVec)){
        annoIndicatorMat = matrix(0, nrow=length(posSNP), ncol=length(annoVecNew))
    for(q in 1:length(annoVecNew)){
        RegionAnnoHeaderListNew[[q]] = strsplit(annoVecNew[q],";")[[1]]
        indiceVec = which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderListNew[[q]])
               #if(length(indiceVec) > 0){
        annoIndicatorMat[indiceVec, q] = 1
                #annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
    }
   }else{
        annoVecNew = annoVec
   }


  }

  annoIndicatorMat_rmind = which(rowSums(annoIndicatorMat) == 0)
  if(length(annoIndicatorMat_rmind) > 0){
    SNP = SNP[-annoIndicatorMat_rmind]

    if(nline_per_gene >= 3){
      WEIGHT = WEIGHT[-annoIndicatorMat_rmind,,drop=F]
    }

    annoIndicatorMat = annoIndicatorMat[-annoIndicatorMat_rmind,,drop=F]
    if(!is.null(markerInfo)){
      genoIndex = genoIndex[-annoIndicatorMat_rmind]
     if(!is.null(RegionData$genoIndex_prev)){
                genoIndex_prev = genoIndex_prev[-annoIndicatorMat_rmind]
      }

    }
  }

  if(nline_per_gene < 3){
	WEIGHT = c(0)
  }	  

   if(!is.null(markerInfo)){
    RegionList[[r]] = list(SNP = SNP,
			   WEIGHT=WEIGHT,
                           annoIndicatorMat = annoIndicatorMat,
                           genoIndex =  as.character(format(genoIndex, scientific = FALSE)),
#                           chrom = uchrom,
                           annoVec = annoVecNew)
    if(!is.null(RegionData$genoIndex_prev)){
	RegionList[[r]]$genoIndex_prev = as.character(format(genoIndex_prev, scientific = FALSE)) 
     }else{
	RegionList[[r]]$genoIndex_prev = c("-1")	
     }

   }else{
    ##VCF
    if(length(SNP) > 0){
        orderposind = order(as.numeric(tstrsplit(SNP, ":")[[2]]))
        SNP = SNP[orderposind]
        annoIndicatorMat = annoIndicatorMat[orderposind,, drop=F]
        if(nrow(WEIGHT) == length(orderposind)){
            WEIGHT = WEIGHT[orderposind,,drop=F]
        }
    }    
    RegionList[[r]] = list(SNP = SNP,
			   WEIGHT=WEIGHT,
                           annoIndicatorMat = annoIndicatorMat,
 #                          chrom = uchrom,
                           annoVec = annoVecNew)
    }
  
   }else{ #if(length(posSNP) > 0){
     RegionList[[r]] = list(SNP = NULL)

   }
}

}else{#if(nrow(RegionData) == 0){
	RegionList = NULL
}


  return(RegionList)
}


mainRegionURV = function(NullModelClass = "SAIGE_NULL_Model",
                         genoType,
                         genoIndex,
                         n)
{
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.mainRegionURV = mainRegionURVInCPP("SAIGE", genoType, genoIndex, n)

  return(obj.mainRegionURV)
}

