#' Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.
#'
#' @param bgenFile character. Path to bgen file. Currently version 1.2 with 8 bit compression is supported
#' @param bgenFileIndex character. Path to the .bgi file (index of the bgen file)
#' @param sampleFile character. Path to the file that contains one column for IDs of samples in the bgen file. The file does not contain header lines.
#' @param vcfFile character. Path to vcf file
#' @param vcfFileIndex character. Path to vcf index file. Indexed by tabix. Path to index for vcf file by tabix, .csi file using 'tabix --csi -p vcf file.vcf.gz'
#' @param vcfField character. genotype field in vcf file to use. "DS" for dosages or "GT" for genotypes. By default, "DS".
#' @param savFile character. Path to sav file
#' @param savFileIndex character. Path to index for sav file .s1r
#' @param bedFile character. Path to bed file (PLINK)
#' @param bimFile character. Path to bim file (PLINK)
#' @param famFile character. Path to fam file (PLINK)
#' @param AlleleOrder character. alt-first or ref-first for bgen or PLINK files. By default, alt-first
#' @param idstoIncludeFile character. Path to a file containing variant ids to be included from the dosage file. The file does not have a header and each line is for a marker ID. Variant ids are in the format chr:pos_ref/alt
#' @param rangestoIncludeFile character. Path to a file containing genome regions to be included from the dosage file. The file contains three columns for chromosome, start, and end respectively with no header. Note for vcf and sav files, only the first line in the file will be used.
#' @param chrom character. If LOCO is specified, chrom is required. chrom is also required for VCF/BCF/SAV input. Note: the string needs to exactly match the chromosome string in the vcf/sav file. For example, 1 does not match chr1.
#' @param is_imputed_data logical. Whether the dosages/genotypes imputed are imputed. If TRUE, the program will output the imputed info score. By default, FALSE.
#' @param minMAC numeric. Minimum minor allele count of markers to test. By default, 0.5. The higher threshold between minMAC and minMAF will be used
#' @param minMAF numeric. Minimum minor allele frequency of markers to test. By default 0. The higher threshold between minMAC and minMAF will be used
#' @param minInfo numeric. Minimum imputation info of markers to test. By default, 0.
#' @param maxMissing numeric. Maximum missing rate for markers to be tested. By default, 0.15
#' @param impute_method character. Imputation method for missing dosages. best_guess, mean or minor. best_guess: missing dosages imputed as best guessed genotyes round(2*allele frequency). mean: missing dosages are imputed as mean (2*allele frequency). minor: missing dosages are imputed as minor allele homozygotes. By default, minor
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out option. If TRUE, --chrom is required. By default, TRUE
#' @param GMMATmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param varianceRatioFile character. Path to the input file containing the variance ratio, which is output from the previous step
#' @param GMMATmodel_varianceRatio_multiTraits_File character. Path to the input file containing 3 columns: phenotype name, model file, and variance ratio file. Each line is for one phenotype. This file is used when multiple phenotypes are analyzed simutaneously
#' @param SAIGEOutputFile character. Prefix of the output files containing assoc test results
#' @param markers_per_chunk character. Number of markers to be tested and output in each chunk in the single-variant assoc tests. By default, 10000
#' @param groups_per_chunk character. Number of groups/sets to be read in and tested in each chunk in the set-based assoc tests. By default, 100
#' @param is_output_moreDetails logical. Whether to output heterozygous and homozygous counts in cases and controls. By default, FALSE. If True, the columns homN_Allele2_cases, hetN_Allelelogical2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls will be output. By default, FALSE
#' @param is_overwrite_output logical. Whether to overwrite the output file if it exists. If FALSE, the program will continue the unfinished analysis instead of starting over from the beginining. By default, TRUE
#' @param maxMAF_in_groupTest. vector of numeric. Max MAF for markers tested in group test seperated by comma. e.g. c(0.0001,0.001,0.01). By default, c(0.01)
#' @param maxMAC_in_groupTest. vector of numeric. Max MAC for markers tested in group test seperated by comma. This vector will be combined with maxMAF_in_groupTest. e.g. c(1) to only test singletons. By default, c(0) and no Max MAC cutoffs are applied.
#' @param minGroupMAC_in_BurdenTest numeric. Only applied when only Burden tests are performed (r.corr=1). Minimum minor allele count in the Burden test for the psueodo marker. By default, 5
#' @param  annotation_in_groupTest. vector of character. annotations of markers to be tested in the set-based tests. using ; to combine multiple annotations in the same test. e.g. c("lof","missense;lof","missense;lof;synonymous")  will test lof variants only, missense+lof variants, and missense+lof+synonymous variants. By default:  c("lof","missense;lof","missense;lof;synonymous")
#' @param groupFile character. Path to the file containing the group information for gene-based tests. Each gene/set has 2 or 3 lines in the group file. The first element is the gene/set name. The second element in the first line is to indicate whether this line contains variant IDs (var), annotations (anno), or weights (weight). The line for weights is optional. If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use weights.beta to change the parameters for the Beta distribution. The variant ids must be in the format chr:pos_ref/alt. Elements are seperated by tab or space.
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file that was used in Step 1
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM
#' @param relatednessCutoff float. The threshold for coefficient of relatedness to treat two samples as unrelated in the sparse GRM. By default, 0
#' @param MACCutoff_to_CollapseUltraRare numeric. MAC cutoff to collpase the ultra rare variants (<= MACCutoff_to_CollapseUltraRare) in the set-based association tests. By default, 10.
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param weights.beta vector of numeric with two elements. parameters for the beta distribution to weight genetic markers in gene-based tests. By default, "c(1,25)".
#' @param r.corr numeric. bewteen 0 and 1. parameters for gene-based tests. If r.corr = 1, only Burden tests will be performed. If r.corr = 0, SKAT-O tests will be performed and results for Burden tests and SKAT tests will be output too.  By default, 0.
#' @param markers_per_chunk_in_groupTest numeric. Number of markers in each chunk when calculating the variance covariance matrix in the set/group-based tests. By default, 100.
#' @param condition character. For conditional analysis. Variant ids are in the format chr:pos_ref/alt and seperated by by comma. e.g."chr3:101651171:C:T,chr3:101651186:G:A".
#' @param weights_for_condition. matrix of numeric. weights for conditioning markers for gene- or region-based tests. The nrow equals to the number of conditioning markers and the ncol equals to the number of sets of weights specified in the group file, e.g. c(1,2,3). If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use weights.beta to change the parameters for the Beta distribution.
#' @param SPAcutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)
#' @param dosage_zerod_cutoff numeric. If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. By derault, 0.2
#' @param dosage_zerod_MAC_cutoff numeric. If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. By derault, 10
#' @param is_single_in_groupTest logical.  Whether to output single-variant assoc test results when perform group tests. Note, single-variant assoc test results will always be output when SKAT and SKAT-O tests are conducted with r.corr=0. This parameter should only be used when only Burden tests are condcuted with r.corr=1. By default, TRUE
#' @param is_equal_weight_in_groupTest logical. Whether equal weights are used in group Test. If TRUE, equal weights will be included in the group test. By default, FALSE
#' @param is_output_markerList_in_groupTest logical. Whether to output the marker lists included in the set-based tests for each mask. By default, FALSE
#' @param is_Firth_beta logical. Whether to estimate effect sizes using approx Firth, only for binary traits. By default, FALSE
#' @param pCutoffforFirth numeric. p-value cutoff to use approx Firth to estiamte the effect sizes. Only for binary traits. The effect sizes of markers with p-value <= pCutoffforFirth will be estimated using approx Firth. By default, 0.01.
#' @param IsOutputlogPforSingle logical. Whether to output log(Pvalue) for single-variant assoc tests. By default, FALSE. If TRUE, the log(Pvalue) instead of original P values will be output (Not activated)
#' @param X_PARregion character. ranges of (pseudoautosomal) PAR region on chromosome X, which are seperated by comma and in the format start:end. By default: '60001-2699520,154931044-155260560' in the UCSC build hg19. For males, there are two X alleles in the PAR region, so PAR regions are treated the same as autosomes. In the NON-PAR regions (outside the specified PAR regions on chromosome X), for males, there is only one X allele. If is_rewrite_XnonPAR_forMales=TRUE, genotypes/dosages of all variants in the NON-PAR regions on chromosome X will be multiplied by 2 (Not activated).
#' @param is_rewrite_XnonPAR_forMales logical. Whether to rewrite gentoypes or dosages of variants in the NON-PAR regions on chromosome X for males (multiply by 2). By default, FALSE. Note, only use is_rewrite_XnonPAR_forMales=TRUE when the specified VCF or Bgen file only has variants on chromosome X. When is_rewrite_XnonPAR_forMales=TRUE, the program does not check the chromosome value by assuming all variants are on chromosome X (Not activated)
#' @param sampleFile_male character. Path to the file containing one column for IDs of MALE samples in the bgen or vcf file with NO header. Order does not matter
#' @return SAIGEOutputFile
#' @export
SPAGMMATtest <- function(bgenFile = "",
                         bgenFileIndex = "",
                         sampleFile = "",
                         vcfFile = "",
                         vcfFileIndex = "",
                         vcfField = "DS",
                         savFile = "",
                         savFileIndex = "",
                         bedFile = "",
                         bimFile = "",
                         famFile = "",
                         AlleleOrder = "alt-first", # new
                         idstoIncludeFile = "",
                         rangestoIncludeFile = "",
                         chrom = "", # for vcf file
                         max_missing = 0.15, # new
                         impute_method = "best_guess", # "mean", "minor", "best_guess"     #new
                         min_MAC = 0.5,
                         min_MAF = 0,
                         min_Info = 0,
                         is_imputed_data = FALSE, # new
                         GMMATmodelFile = "",
                         LOCO = TRUE,
                         varianceRatioFile = "",
                         GMMATmodel_varianceRatio_multiTraits_File = "",
                         cateVarRatioMinMACVecExclude = c(10.5, 20.5),
                         cateVarRatioMaxMACVecInclude = c(20.5),
                         SPAcutoff = 2,
                         SAIGEOutputFile = "",
                         markers_per_chunk = 10000,
                         groups_per_chunk = 100,
                         markers_per_chunk_in_groupTest = 100, # new
                         condition = "",
                         sparseGRMFile = "",
                         sparseGRMSampleIDFile = "",
                         VmatFilelist = "",
                         VmatSampleFilelist = "",
                         relatednessCutoff = 0,
                         groupFile = "",
                         weights.beta = c("1,25"),
                         weights_for_condition = NULL,
                         r.corr = 0,
                         dosage_zerod_cutoff = 0.2,
                         dosage_zerod_MAC_cutoff = 10,
                         is_output_moreDetails = FALSE, # new
                         MACCutoff_to_CollapseUltraRare = 10,
                         annotation_in_groupTest = c("lof", "missense;lof", "missense;lof;synonymous"), # new
                         maxMAF_in_groupTest = c(0.01),
                         minMAF_in_groupTest_Exclude = NULL,
                         maxMAC_in_groupTest = c(0),
                         minMAC_in_groupTest_Exclude = NULL,
                         minGroupMAC_in_BurdenTest = 5,
                         is_Firth_beta = FALSE,
                         pCutoffforFirth = 0.01,
                         is_overwrite_output = TRUE,
                         is_single_in_groupTest = TRUE,
                         is_SKATO = FALSE,
                         is_equal_weight_in_groupTest = FALSE,
                         is_output_markerList_in_groupTest = FALSE,
                         pval_cutoff_for_fastTest = 0.05,
                         is_fastTest = FALSE,
                         pval_cutoff_for_gxe = 0.001,
                         is_noadjCov = TRUE,
                         is_sparseGRM = TRUE,
                         max_MAC_use_ER = 4,
                         is_EmpSPA = FALSE) {
  # cat("r.corr is ", r.corr, "\n")
  if (!(impute_method %in% c("best_guess", "mean", "minor"))) {
    stop("impute_method should be 'best_guess', 'mean' or 'minor'.")
  }

  checkArgsListBool(
    is_imputed_data = is_imputed_data,
    LOCO = LOCO,
    is_output_moreDetails = is_output_moreDetails,
    is_overwrite_output = is_overwrite_output
  )
  # is_rewrite_XnonPAR_forMales = is_rewrite_XnonPAR_forMales)
  cat("dosage_zerod_cutoff ", dosage_zerod_cutoff, "\n")
  checkArgsListNumeric(
    start = 1,
    end = 250000000,
    max_missing = max_missing,
    min_MAC = min_MAC,
    min_MAF = min_MAF,
    min_Info = min_Info,
    SPAcutoff = SPAcutoff,
    dosage_zerod_cutoff = dosage_zerod_cutoff,
    dosage_zerod_MAC_cutoff = dosage_zerod_MAC_cutoff,
    markers_per_chunk = markers_per_chunk,
    groups_per_chunk = groups_per_chunk,
    minGroupMAC_in_BurdenTest = minGroupMAC_in_BurdenTest,
    max_MAC_use_ER = max_MAC_use_ER
  )


  # if(file.exists(SAIGEOutputFile)) {print("ok -2 file exist")}
  # print("setSAIGEobjInCPP -3")
  # print_g_n_unique()


  ## check and create the output file
  # Check_OutputFile_Create(SAIGEOutputFile)
  OutputFile <- SAIGEOutputFile
  OutputFileIndex <- NULL
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }


  if (!is.null(minMAF_in_groupTest_Exclude)) {
    min_MAF <- min(min_MAF, min(minMAF_in_groupTest_Exclude))
  }
  if (!is.null(minMAC_in_groupTest_Exclude)) {
    min_MAC <- min(min_MAF, min(minMAC_in_groupTest_Exclude))
  }

  setAssocTest_GlobalVarsInCPP(
    impute_method,
    max_missing,
    min_MAF,
    min_MAC,
    min_Info,
    dosage_zerod_cutoff,
    dosage_zerod_MAC_cutoff,
    OutputFile,
    max_MAC_use_ER
  )

  if (groupFile == "") {
    isGroupTest <- FALSE
    cat("single-variant association test will be performed\n")

    setMarker_GlobalVarsInCPP(
      is_output_moreDetails,
      markers_per_chunk
    )
  } else {
    isGroupTest <- TRUE
    Check_File_Exist(groupFile, "groupFile")
    cat("group-based test will be performed\n")
    # checkArgsList_for_Region(method_to_CollapseUltraRare,
    # order the max MAF from lowest to highest
    # maxMAF_in_groupTest = maxMAF_in_groupTest[order(maxMAF_in_groupTest)]
    # maxMAC_in_groupTest = maxMAC_in_groupTest[order(maxMAC_in_groupTest)]

    cat("maxMAF_in_groupTest ", maxMAF_in_groupTest, "\n")
    cat("minMAF_in_groupTest_Exclude ", minMAF_in_groupTest_Exclude, "\n")


    checkArgsList_for_Region(
      MACCutoff_to_CollapseUltraRare,
      # DosageCutoff_for_UltraRarePresence,
      maxMAF_in_groupTest = maxMAF_in_groupTest,
      minMAF_in_groupTest_Exclude = minMAF_in_groupTest_Exclude,
      maxMAC_in_groupTest = maxMAC_in_groupTest,
      minMAC_in_groupTest_Exclude = minMAC_in_groupTest_Exclude,
      markers_per_chunk_in_groupTest = markers_per_chunk_in_groupTest
    )



    # if(file.exists(SAIGEOutputFile)) {print("ok -1 file exist")}

    IsOutputlogPforSingle <- FALSE # to check
    # OUT_Filename_Single<-sprintf("%s.single",SAIGEOutputFile)
    # Check_OutputFile_Create(OUT_Filename_Single)
    # if (sum(weights.beta.rare != weights.beta.common) > 0) {
    #  cat("WARNING:The option for weights.beta.common is not fully developed\n")
    #  cat("weights.beta.common is set to be equal to weights.beta.rare\n")
    #  weights.beta.common = weights.beta.rare
    # }

    # method_to_CollapseUltraRare,
    # DosageCutoff_for_UltraRarePresence,
    setRegion_GlobalVarsInCPP(
      maxMAF_in_groupTest,
      markers_per_chunk_in_groupTest,
      MACCutoff_to_CollapseUltraRare,
      minGroupMAC_in_BurdenTest
    )
    # cat("dosage_zerod_cutoff is ", dosage_zerod_cutoff, "\n")
    # cat("dosage_zerod_MAC_cutoff is ", dosage_zerod_MAC_cutoff, "\n")
  }
  # print("setSAIGEobjInCPP -2")
  # print_g_n_unique()
  if (GMMATmodel_varianceRatio_multiTraits_File != "") {
    if (GMMATmodelFile != "" | varianceRatioFile != "") {
      stop("GMMATmodel_varianceRatio_multiTraits_File is specified while varianceRatioFile and/or GMMATmodelFile are also specified. Please check\n")
    } else {
      modelvrfile_data <- data.table::fread(GMMATmodel_varianceRatio_multiTraits_File, header = F, data.table = F)
      if (ncol(modelvrfile_data) != 3) {
        stop("GMMATmodel_varianceRatio_multiTraits_File needs to have 3 columns: phenotype name, model file, and variance ratio file\n")
      } else {
        GMMATmodelFile_vec <- modelvrfile_data[, 2]
        GMMATmodelFile <- paste(GMMATmodelFile_vec, collapse = ",")
        varianceRatioFile_vec <- modelvrfile_data[, 3]
        varianceRatioFile <- paste(varianceRatioFile_vec, collapse = ",")
        phenotype_name_vec <- as.character(modelvrfile_data[, 1])
      }
    }
  } else {
    GMMATmodelFile_vec <- unlist(strsplit(GMMATmodelFile, split = ","))
    phenotype_name_vec <- as.character(seq(1, length(GMMATmodelFile_vec)))
  }

  obj.model.List <- ReadModel_multiTrait(GMMATmodelFile, chrom, LOCO, is_Firth_beta, is_EmpSPA, espa_nt = 9999, espa_range = c(-20, 20)) # readInGLMM.R8


  if (obj.model.List[[1]]$traitType == "binary") {
    if (max_MAC_use_ER > 0) {
      cat("P-values of genetic variants with MAC <= ", max_MAC_use_ER, " will be calculated via effecient resampling.\n")
    }
  } else {
    max_MAC_use_ER <- 0
  }

  if (!LOCO) {
    # 	LOCO = FALSE
    print("LOCO = FASLE and leave-one-chromosome-out is not applied")
  }

  isSparseGRM <- is_sparseGRM


  cat("isSparseGRM ", isSparseGRM, "\n")

  set_Vmat_vec_orig(VmatFilelist, VmatSampleFilelist, obj.model.List[[1]]$sampleID)

  ratioVecList <- Get_Variance_Ratio_multiTrait(varianceRatioFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest, isSparseGRM) # readInGLMM.R
  # print("ratioVecList")
  # print(ratioVecList)


  # pval_cutoff_for_fastTest = 0

  if (!is_fastTest) {
    pval_cutoff_for_fastTest <- 1
  }


  nsample <- length(unique(obj.model.List[[1]]$sampleID))
  cateVarRatioMaxMACVecInclude <- c(cateVarRatioMaxMACVecInclude, nsample)

  # print("setSAIGEobjInCPP -1b")
  # print_g_n_unique()

  # in Geno.R
  objGeno <- setGenoInput(
    bgenFile = bgenFile,
    bgenFileIndex = bgenFileIndex,
    vcfFile = vcfFile, # not activate yet
    vcfFileIndex = vcfFileIndex,
    vcfField = vcfField,
    savFile = savFile,
    savFileIndex = savFileIndex,
    sampleFile = sampleFile,
    bedFile = bedFile,
    bimFile = bimFile,
    famFile = famFile,
    idstoIncludeFile = idstoIncludeFile,
    rangestoIncludeFile = rangestoIncludeFile,
    chrom = chrom,
    AlleleOrder = AlleleOrder,
    sampleInModel = obj.model.List[[1]]$sampleID
  )

  genoType <- objGeno$genoType
  if (condition != "") {
    isCondition <- TRUE
  } else {
    isCondition <- FALSE
  }

  # print("setSAIGEobjInCPP -1a")
  # print_g_n_unique()
  condition_genoIndex_a <- c(-1)
  condition_genoIndex <- c(-1)
  if (isCondition) {
    cat("Conducting conditional analysis. Please specify the conditioning markers in the order as they are store in the genotype/dosage file.\n")
    condition_genoIndex <- extract_genoIndex_condition(condition, objGeno$markerInfo, genoType)
    # }else{
    condition_genoIndex_a <- condition_genoIndex$cond_genoIndex
    print("condition_genoIndex_a")
    print(condition_genoIndex_a)
  }
  # set up the SAIGE object based on the null model results
  # print("SigmaMat_sp")
  # print(SigmaMat_sp)


  # print(names(obj.model))
  # print(names(obj.model$obj.noK))
  # obj.model$varWeights = rep(1, length(obj.model$y))
  # print(obj.model$obj_cc$res.out)


  # print("SigmaMat_sp")
  # print(SigmaMat_sp)

  b <- as.numeric(factor(obj.model.List[[1]]$sampleID, levels = unique(obj.model.List[[1]]$sampleID)))
  I_mat <- Matrix::sparseMatrix(i = 1:length(b), j = b, x = rep(1, length(b)))
  I_mat <- 1.0 * I_mat
  #     set_I_longl_mat_SAIGEtest(I_mat, b-1)
  if (!is.null(obj.model.List[[1]]$T_longl_vec)) {
    T_longl_mat <- I_mat * (obj.model.List[[1]]$T_longl_vec)
    # set_T_longl_mat_SAIGEtest(T_longl_mat, obj.model$T_longl_vec)
  } else {
    obj.model.List[[1]]$T_longl_vec <- rep(1, length(b))
    T_longl_mat <- I_mat * (obj.model.List[[1]]$T_longl_vec)
  }

  # print("obj.model$spSigma ")
  # print(obj.model$spSigma)
  if (!is.null(obj.model.List[[1]]$spSigma)) {
    # SigmaMat_sp = getSparseSigma_new()
    isSparseGRM <- TRUE
    SigmaMat_sp <- NULL
    # SigmaMat_sp = chol2inv(chol(obj.model.List[[1]]$spSigma))
    for (gm in 1:length(obj.model.List)) {
      SigmaMat_sp <- cbind(SigmaMat_sp, obj.model.List[[gm]]$spSigma)
      print(dim(SigmaMat_sp))
      # SigmaMat_sp = SigmaMat_sp %*% I_mat
      # SigmaMat_sp = t(I_mat)%*%SigmaMat_sp
      # print(dim(SigmaMat_sp))
    }
    cat("isSparseGRM 2 ", isSparseGRM, "\n")
  } else {
    SigmaMat_sp <- Matrix:::sparseMatrix(i = c(1, 1, 2, 2), j = c(1, 2, 1, 2), x = as.vector(c(0, 0, 0, 0)))
  }




  # print("sum(!duplicated(obj.model$X))")
  # print(sum(!duplicated(obj.model$X)))
  # print("length(unique(obj.model$sampleID))")
  # print(length(unique(obj.model$sampleID)))
  eMat <- NULL
  isgxe_vec <- NULL
  for (oml in 1:length(obj.model.List)) {
    # print("dim(I_mat)")
    # print(dim(I_mat))
    # print("obj.model.List[[oml]]$eMat)")
    # print(dim(obj.model.List[[oml]]$eMat))
    # eMat = cbind(eMat, t(I_mat)%*%(obj.model.List[[oml]]$eMat))
    eMat <- cbind(eMat, (obj.model.List[[oml]]$eMat))
    isgxe_vec <- c(isgxe_vec, obj.model.List[[oml]]$isgxe)
  }

  if (sum(isgxe_vec) != 0 && sum(isgxe_vec) != length(isgxe_vec)) {
    stop("isgxe for all traits needs to be the same (perform or not perform dynamic qtl analyses)\n")
  }

  cat("pval_cutoff_for_gxe ", pval_cutoff_for_gxe, "\n")
  # print("dim(eMat)")
  # print(dim(eMat))
  # print(isgxe_vec)
  eMat <- as.matrix(eMat)
  # setAssocTest_GlobalVarsInCPP_GbyE(eMat, TRUE, 0.001)
  XV_gxe <- NULL
  X_gxe <- NULL
  XVX_inv_XV_gxe <- NULL
  XVX_gxe <- NULL
  S_a_gxe <- NULL
  XXVX_inv_gxe <- NULL
  y_gxe <- NULL
  res_gxe <- NULL
  mu2_gxe <- NULL
  mu_gxe <- NULL

  varWeights_gxe <- NULL

  ## Set up the sample level summary statistics. This is specifically for data sets with repeated measurements
  if (sum(duplicated(obj.model.List[[1]]$sampleID)) > 0) {
    Xsample <- list()
    Vsample <- list()
    XVsample <- list()
    XVXsample <- list()
    XXVXsample_inv <- list()
    XVX_inv_XVsample <- list()
    Sigma_iXXSigma_iX <- list()
    res_sample <- NULL
    mu_sample <- NULL
    mu2_sample <- NULL
    S_a_sample <- NULL
    theta <- NULL
    y <- NULL
    offset <- NULL
    obj_cc_res.out <- NULL
    tauVal_sp <- NULL
    traitType <- NULL
    varWeights_sample <- NULL

    for (oml in 1:length(obj.model.List)) {
      obj.model <- obj.model.List[[oml]]
      Xsample0 <- obj.model$sampleXMat ## from step 1
      Xsample[[oml]] <- data.table::as.data.table(Xsample0)
      # print(dim(I_mat))
      # print(length(obj.model$obj.noK$V))
      Vsample0 <- as.vector(t(obj.model$obj.noK$V) %*% I_mat)
      Vsample[[oml]] <- data.table::as.data.table(Vsample0)
      # print(dim(Xsample))
      # print(length(Vsample))
      XVsample0 <- t(Xsample0 * Vsample0)
      XVsample[[oml]] <- data.table::as.data.table(XVsample0)
      XVXsample0 <- t(Xsample0) %*% (t(XVsample0))
      XVXsample_inv0 <- solve(XVXsample0)
      XVXsample[[oml]] <- data.table::as.data.table(XVXsample0)
      XXVXsample_inv0 <- Xsample0 %*% XVXsample_inv0
      XVX_inv_XVsample0 <- XXVXsample_inv0 * Vsample0

      XXVXsample_inv[[oml]] <- data.table::as.data.table(XXVXsample_inv0)

      XVX_inv_XVsample[[oml]] <- data.table::as.data.table(XVX_inv_XVsample0)

      Sigma_iXXSigma_iX0 <- obj.model$Sigma_iXXSigma_iX
      Sigma_iXXSigma_iX[[oml]] <- data.table::as.data.table(Sigma_iXXSigma_iX0)

      # print(dim(XVXsample))
      # print("dim(XXVXsample_inv)")
      # print(dim(XXVXsample_inv))
      # print("dim(XVsample)")
      # print(dim(XVsample))
      # print("dim(XVX_inv_XVsample)")
      # print(dim(XVX_inv_XVsample))
      # print("dim(Sigma_iXXSigma_iX)")
      # print(dim(Sigma_iXXSigma_iX))
      # print("dim(Xsample)")
      # print(dim(Xsample))

      res_sample0 <- as.vector(t(I_mat) %*% (obj.model$residuals))
      mu_sample0 <- as.vector(t(I_mat) %*% (obj.model$mu))
      mu2_sample0 <- as.vector(t(I_mat) %*% (obj.model$mu2))
      S_a_sample0 <- rowSums(t(Xsample0) * res_sample0)

      res_sample <- cbind(res_sample, res_sample0)
      mu_sample <- cbind(mu_sample, mu_sample0)
      mu2_sample <- cbind(mu2_sample, mu2_sample0)
      S_a_sample <- cbind(S_a_sample, S_a_sample0)

      # print("dim(S_a_sample)")
      # print(dim(S_a_sample))
      # print("dim(res_sample)")
      # print(dim(res_sample))
      # print("dim(mu2_sample)")
      # print(dim(mu2_sample))
      # print("dim(mu_sample)")
      # print(dim(mu_sample))
      uniqsampleind <- which(!duplicated(obj.model$sampleID))
      varWeights_sample <- cbind(varWeights_sample, obj.model$varWeights[uniqsampleind])

      theta <- cbind(theta, obj.model$theta)
      y <- cbind(y, obj.model$y)
      offset <- cbind(offset, obj.model$offset)
      obj_cc_res.out <- cbind(obj_cc_res.out, obj.model$obj_cc$res.out)
      # tauVal_sp = cbind(tauVal_sp, obj.model$tauVal_sp)
      traitType <- c(traitType, obj.model$traitType)
      if (isgxe_vec[1]) {
        XV_gxe <- rbind(XV_gxe, obj.model$obj.noK$XV)
        XXVX_inv_gxe <- rbind(XXVX_inv_gxe, obj.model$obj.noK$XXVX_inv)
        # X_gxe = rbind(X_gxe, obj.model$X)
        # XVX_inv_XV_gxe = rbind(XVX_inv_XV_gxe, obj.model$obj.noK$XVX_inv_XV)
        # XVX_gxe = rbind(XVX_gxe,  obj.model$obj.noK$XVX)
        X_gxe <- matrix(1)
        XVX_inv_XV_gxe <- matrix(1)
        XVX_gxe <- matrix(1)
        y_gxe <- cbind(y_gxe, obj.model$y)
        res_gxe <- cbind(res_gxe, obj.model$residuals)
        mu2_gxe <- cbind(mu2_gxe, obj.model$mu2)
        mu_gxe <- cbind(mu_gxe, obj.model$mu)
        # S_a_gxe_sub = (as.matrix(obj.model$X)) * (as.vector(obj.model$residuals))
        # S_a_gxe = cbind(S_a_gxe, colSums(S_a_gxe_sub))
        S_a_gxe <- matrix(1)
        varWeights_gxe <- cbind(varWeights_gxe, obj.model$varWeights)
      } else {
        XV_gxe <- matrix(1)
        XXVX_inv_gxe <- matrix(1)
        X_gxe <- matrix(1)
        XVX_inv_XV_gxe <- matrix(1)
        XVX_gxe <- matrix(1)
        y_gxe <- matrix(1)
        res_gxe <- matrix(1)
        mu2_gxe <- matrix(1)
        mu_gxe <- matrix(1)
        S_a_gxe <- matrix(1)
        varWeights_gxe <- matrix(1)
      }
    }
    Xsample <- as.matrix(data.table::rbindlist(Xsample))
    Vsample <- as.matrix(data.table::rbindlist(Vsample))
    XVsample <- as.matrix(data.table::rbindlist(XVsample))
    XVXsample <- as.matrix(data.table::rbindlist(XVXsample))
    XXVXsample_inv <- as.matrix(data.table::rbindlist(XXVXsample_inv))
    XVX_inv_XVsample <- as.matrix(data.table::rbindlist(XVX_inv_XVsample))
    Sigma_iXXSigma_iX <- as.matrix(data.table::rbindlist(Sigma_iXXSigma_iX))
  }
  gc()


  # print("dim(ratioVecList$ratioVec_sparse)")
  # print(dim(ratioVecList$ratioVec_sparse))

  # print(ratioVecList)
  # print(SPAcutoff)
  # print(theta)
  # print(dim(varWeights_sample))
  # print(traitType)
  # print(dim(y))
  # print(is_noadjCov)
  # print(pval_cutoff_for_fastTest)
  # print(condition_genoIndex)
  # print(is_Firth_beta)
  # print(pCutoffforFirth)

  # print(dim(offset))
  # print(dim(obj_cc_res.out))
  # print(dim(SigmaMat_sp))
  # print(obj.model$tauVal_sp)
  # print(dim(I_mat))
  # print(b-1)
  # print(dim(T_longl_mat))
  # print(length(obj.model$T_longl_vec))
  # print(obj.model$T_longl_vec)
  # print("ok")
  # print(obj.model$cumul)

  # print("XXVXsample_inv")
  # print(XXVXsample_inv)
  if (sum(duplicated(obj.model.List[[1]]$sampleID)) > 0) {
    if (FALSE) {
      print("XXVXsample_inv")
      print(length(XVXsample))
      print(length(XXVXsample_inv))
      print(length(XVsample))
      print(XXVXsample_inv)
      print(XVsample)
      print(length(XVX_inv_XVsample))
      print(length(Sigma_iXXSigma_iX))
      print(length(Xsample))
      print(length(S_a_sample))
      print(length(res_sample))
      print(length(mu2_sample))
      print(length(mu_sample))
      print(ratioVecList)
      print(length(as.matrix(ratioVecList$ratioVec_sparse)))
      print(length(as.matrix(ratioVecList$ratioVec_null)))
      print(length(as.matrix(ratioVecList$ratioVec_null_noXadj)))
      print(length(as.matrix(ratioVecList$ratioVec_null_eg)))
      print(length(as.matrix(ratioVecList$ratioVec_sparse_eg)))
      print(length(cateVarRatioMinMACVecExclude))
      print(length(cateVarRatioMaxMACVecInclude))
      print(length(SPAcutoff))
      print(length(theta))
      print(length(varWeights_sample))
      print(length(traitType))
      print(length(y))
      print(length(impute_method))
      print(length(isSparseGRM))
      print(length(is_noadjCov))
      print(length(pval_cutoff_for_fastTest))
      print(length(isCondition))
      print(length(condition_genoIndex_a))
      print(length(is_Firth_beta))
      print(length(pCutoffforFirth))
      print(length(offset))
      print(length(obj_cc_res.out))
      print(length(SigmaMat_sp))
      print(length(obj.model$tauVal_sp))
      print(obj.model$tauVal_sp)
      print("ok")
      print(length(I_mat))
      print(length(b - 1))
      print(length(T_longl_mat))
      print(length(obj.model.List[[1]]$T_longl_vec))
      print(length(is_EmpSPA))
      print(obj.model$cumul)
      print(isgxe_vec[1])
      print(dim(XV_gxe))
      print(dim(X_gxe))
      print(dim(XVX_inv_XV_gxe))
      print(dim(XVX_gxe))
      print(dim(S_a_gxe))
      print(dim(XXVX_inv_gxe))
      print(dim(y_gxe))
      print(dim(res_gxe))
      print(dim(mu2_gxe))
      print(dim(mu_gxe))
      print(varWeights_gxe)
    }

    # print("OKKKK")
    # print(traitType)


    # if(FALSE){
    setSAIGEobjInCPP(
      t_XVX = XVXsample,
      t_XXVX_inv = XXVXsample_inv,
      t_XV = XVsample,
      t_XVX_inv_XV = XVX_inv_XVsample,
      t_Sigma_iXXSigma_iX = Sigma_iXXSigma_iX, ## specifically for sparse V, check later
      t_X = Xsample,
      t_S_a = S_a_sample,
      t_res = res_sample,
      t_mu2 = mu2_sample,
      t_mu = mu_sample,
      t_varRatio_sparse = as.matrix(ratioVecList$ratioVec_sparse),
      t_varRatio_null = as.matrix(ratioVecList$ratioVec_null),
      t_varRatio_null_sample = as.matrix(ratioVecList$ratioVec_null_sample),
      t_varRatio_null_noXadj = as.matrix(ratioVecList$ratioVec_null_noXadj),
      t_varRatio_null_eg = as.matrix(ratioVecList$ratioVec_null_eg),
      t_varRatio_sparse_eg = as.matrix(ratioVecList$ratioVec_sparse_eg),
      t_cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
      t_cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
      t_SPA_Cutoff = SPAcutoff,
      t_tauvec = theta,
      t_varWeightsvec = varWeights_sample,
      t_traitType = as.vector(traitType),
      t_y = as.matrix(y),
      t_impute_method = impute_method,
      t_flagSparseGRM = isSparseGRM,
      t_isnoadjCov = is_noadjCov,
      t_pval_cutoff_for_fastTest = pval_cutoff_for_fastTest,
      t_isCondition = isCondition,
      t_condition_genoIndex = condition_genoIndex_a,
      t_is_Firth_beta = is_Firth_beta,
      t_pCutoffforFirth = pCutoffforFirth,
      t_offset = offset, ## check later
      t_resout = obj_cc_res.out, ## for ER, check later
      t_SigmaMat_sp = SigmaMat_sp,
      t_tauVal_sp = obj.model$tauVal_sp,
      t_Ilongmat = I_mat,
      t_I_longl_vec = b - 1,
      t_Tlongmat = T_longl_mat,
      t_T_longl_vec = obj.model.List[[1]]$T_longl_vec,
      t_is_EmpSPA = is_EmpSPA,
      t_cumul = obj.model$cumul,
      t_is_gxe = isgxe_vec[1],
      t_XV_gxe = XV_gxe,
      t_X_gxe = X_gxe,
      t_XVX_inv_XV_gxe = XVX_inv_XV_gxe,
      t_XVX_gxe = XVX_gxe,
      t_S_a_gxe = S_a_gxe,
      t_XXVX_inv_gxe = XXVX_inv_gxe,
      t_y_gxe = y_gxe,
      t_res_gxe = res_gxe,
      t_mu2_gxe = mu2_gxe,
      t_mu_gxe = mu_gxe,
      t_varWeights_gxe = varWeights_gxe
    )
    # }



    # print("setSAIGEobjInCPP 0")
    # print_g_n_unique()

    # print("ratioVecList")
    # print(ratioVecList)
  } else {
    # }
    X <- NULL
    V <- NULL
    XV <- NULL
    XVX <- NULL
    XXVX_inv <- NULL
    XVX_inv_XV <- NULL
    Sigma_iXXSigma_iX <- NULL
    res <- NULL
    mu <- NULL
    mu2 <- NULL
    S_a <- NULL
    theta <- NULL
    y <- NULL
    offset <- NULL
    obj_cc_res.out <- NULL
    tauVal_sp <- NULL
    traitType <- NULL
    varWeights <- NULL
    for (oml in 1:length(obj.model.List)) {
      obj.model <- obj.model.List[[oml]]
      traitType <- c(traitType, obj.model$traitType)
      # XVX = rbind(XVX, obj.model$obj.noK$XVX)
      XXVX_inv <- rbind(XXVX_inv, obj.model$obj.noK$XXVX_inv)
      XV <- rbind(XV, obj.model$obj.noK$XV)
      # XVX_inv_XV = rbind(XVX_inv_XV, obj.model$obj.noK$XVX_inv_XV)
      Sigma_iXXSigma_iX <- rbind(Sigma_iXXSigma_iX, obj.model$Sigma_iXXSigma_iX)
      # X = rbind(X, obj.model$X)
      # S_a = cbind(S_a, obj.model$obj.noK$S_a)
      res <- cbind(res, obj.model$residuals)
      mu2 <- cbind(mu2, obj.model$mu2)
      mu <- cbind(mu, obj.model$mu)

      theta <- cbind(theta, obj.model$theta)
      y <- cbind(y, obj.model$y)
      offset <- cbind(offset, obj.model$offset)
      obj_cc_res.out <- cbind(obj_cc_res.out, obj.model$obj_cc$res.out)
      tauVal_sp <- cbind(tauVal_sp, obj.model$tauVal_sp)
      varWeights <- cbind(varWeights, obj.model$varWeights)
    }

    if (isgxe_vec[1]) {
      XV_gxe <- XV
      XXVX_inv_gxe <- XXVX_inv
      # X_gxe = X
      X_gxe <- matrix(1)
      # XVX_inv_XV_gxe = XVX_inv_XV
      XVX_inv_XV_gxe <- matrix(1)
      # XVX_gxe = XVX
      XVX_gxe <- matrix(1)
      y_gxe <- y
      res_gxe <- res
      mu2_gxe <- mu2
      mu_gxe <- mu
      varWeights_gxe <- varWeights
      # S_a_gxe = S_a
      S_a_gxe <- matrix(1)
    } else {
      XV_gxe <- matrix(1)
      XXVX_inv_gxe <- matrix(1)
      X_gxe <- matrix(1)
      XVX_inv_XV_gxe <- matrix(1)
      XVX_gxe <- matrix(1)
      y_gxe <- matrix(1)
      res_gxe <- matrix(1)
      mu2_gxe <- matrix(1)
      mu_gxe <- matrix(1)
      varWeights_gxe <- matrix(1)
      S_a_gxe <- matrix(1)
    }

    mu_sample <- mu

    setSAIGEobjInCPP(
      t_XVX = XVX,
      t_XXVX_inv = XXVX_inv,
      t_XV = XV,
      t_XVX_inv_XV = XVX_inv_XV,
      t_Sigma_iXXSigma_iX = Sigma_iXXSigma_iX,
      t_X = X,
      t_S_a = S_a,
      t_res = res,
      t_mu2 = mu2,
      t_mu = mu,
      t_varRatio_sparse = as.matrix(ratioVecList$ratioVec_sparse),
      t_varRatio_null = as.matrix(ratioVecList$ratioVec_null),
      t_varRatio_null_sample = as.matrix(ratioVecList$ratioVec_null_sample),
      t_varRatio_null_noXadj = as.matrix(ratioVecList$ratioVec_null_noXadj),
      t_varRatio_null_eg = as.matrix(ratioVecList$ratioVec_null_eg),
      t_varRatio_sparse_eg = as.matrix(ratioVecList$ratioVec_sparse_eg),
      t_cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
      t_cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
      t_SPA_Cutoff = SPAcutoff,
      t_tauvec = theta,
      t_varWeightsvec = varWeights,
      t_traitType = traitType,
      t_y = y,
      t_impute_method = impute_method,
      t_flagSparseGRM = isSparseGRM,
      t_isnoadjCov = is_noadjCov,
      t_pval_cutoff_for_fastTest = pval_cutoff_for_fastTest,
      t_isCondition = isCondition,
      t_condition_genoIndex = condition_genoIndex_a,
      t_is_Firth_beta = is_Firth_beta,
      t_pCutoffforFirth = pCutoffforFirth,
      t_offset = offset,
      t_resout = obj_cc_res.out,
      t_SigmaMat_sp = SigmaMat_sp,
      t_tauVal_sp = obj.model$tauVal_sp,
      t_Ilongmat = I_mat,
      t_I_longl_vec = b - 1,
      t_Tlongmat = T_longl_mat,
      t_T_longl_vec = obj.model.List[[1]]$T_longl_vec,
      t_is_EmpSPA = is_EmpSPA,
      t_cumul = obj.model$cumul,
      t_is_gxe = isgxe_vec[1],
      t_XV_gxe = XV_gxe,
      t_X_gxe = X_gxe,
      t_XVX_inv_XV_gxe = XVX_inv_XV_gxe,
      t_XVX_gxe = XVX_gxe,
      t_S_a_gxe = S_a_gxe,
      t_XXVX_inv_gxe = XXVX_inv_gxe,
      t_y_gxe = y_gxe,
      t_res_gxe = res_gxe,
      t_mu2_gxe = mu2_gxe,
      t_mu_gxe = mu_gxe,
      t_varWeights_gxe = varWeights_gxe
    )
  }

  # if(any(duplicated(obj.model$sampleID))){
  # 	b = as.numeric(factor(obj.model$sampleID, levels =  unique(obj.model$sampleID)))
  # 	I_mat = Matrix::sparseMatrix(i = 1:length(b), j = b, x = rep(1, length(b)))
  # 	I_mat = 1.0 * I_mat
  # 	set_I_longl_mat_SAIGEtest(I_mat, b-1)
  # 	if(!is.null(obj.model$T_longl_vec)){
  # 		T_longl_mat = I_mat * (obj.model$T_longl_vec)
  # 		set_T_longl_mat_SAIGEtest(T_longl_mat, obj.model$T_longl_vec)
  # 	}
  # }

  # rm(sparseSigmaRList)
  gc()


  setAssocTest_GlobalVarsInCPP_GbyE(eMat, isgxe_vec[1], as.numeric(pval_cutoff_for_gxe), XV_gxe, XXVX_inv_gxe, y_gxe, res_gxe, mu2_gxe, mu_gxe, varWeights_gxe)

  # process condition
  if (isCondition) {
    # n = length(obj.model$y) #sample size
    n <- ncol(I_mat)
    # n_uniq = length(unique(obj.model))
    ## re-order the conditioning markers
    ## condition_original = unlist(strsplit(condition, ","))
    # condition_genoIndex=extract_genoIndex_condition(condition, objGeno$markerInfo, genoType)
    print("condition_genoIndex")
    print(condition_genoIndex)

    if (isGroupTest) {
      if (!is.null(weights_for_condition)) {
        condition_weights <- as.matrix(weights_for_condition)
        print(condition_weights)
        # print(condition_genoIndex$cond_genoIndex)
        # condition_weights = as.numeric(unlist(strsplit(weights_for_condition, ",")))
        if (nrow(condition_weights) != length(condition_genoIndex$cond_genoIndex)) {
          stop("The length of the provided weights for conditioning markers is not equal to the number of conditioning markers\n")
        }
      } else {
        condition_weights <- matrix(rep(0, length(condition_genoIndex$cond_genoIndex)), ncol = 1)
      }


      if (!is.null(weights.beta)) {
        BetaDist_weight_mat <- NULL
        for (i in 1:length(weights.beta)) {
          weightsbeta_val_vec <- as.numeric(unlist(strsplit(weights.beta[i], split = ",")))
          if (length(weightsbeta_val_vec) == 2) {
            BetaDist_weight_mat <- rbind(BetaDist_weight_mat, weightsbeta_val_vec)
          } else {
            stop("The ", i, "th element in weights.beta does not have 2 elements\n")
          }
        }
        BetaDist_weight_mat <- as.matrix(BetaDist_weight_mat)
      } else {
        BetaDist_weight_mat <- matrix(c(0, 0), ncol = 2)
      }
    } else { # if(isGroupTest){
      BetaDist_weight_mat <- matrix(c(0, 0), ncol = 2)
      condition_weights <- matrix(rep(0, length(condition_genoIndex$cond_genoIndex)), ncol = 1)
    }


    condition_genoIndex_a <- as.character(format(condition_genoIndex$cond_genoIndex, scientific = FALSE))
    condition_genoIndex_prev_a <- as.character(format(condition_genoIndex$cond_genoIndex_prev, scientific = FALSE))
    # 	print("OKKK")

    print("condition_genoIndex_prev_a")
    print(condition_genoIndex_prev_a)
    print("condition_genoIndex_a")
    print(condition_genoIndex_a)
    print("condition_weights")
    print(condition_weights)
    print("BetaDist_weight_mat")
    print(BetaDist_weight_mat)
    BetaDist_weight_mat <- as.matrix(BetaDist_weight_mat)
    print("BetaDist_weight_mat")
    print(dim(BetaDist_weight_mat))



    assign_conditionMarkers_factors(genoType, condition_genoIndex_prev_a, condition_genoIndex_a, n, condition_weights, BetaDist_weight_mat, is_equal_weight_in_groupTest)

    # 	print("OKKK2")
    if (obj.model$traitType[1] == "binary" & isGroupTest) {
      outG2cond <- RegionSetUpConditional_binary_InCPP(condition_weights)
      G2condList_list <- NULL
      for (oml in 1:length(obj.model.List)) {
        startcond <- (oml - 1) * length(condition_genoIndex$cond_genoIndex) + 1
        endcond <- oml * length(condition_genoIndex$cond_genoIndex)


        G2condList <- get_newPhi_scaleFactor(q.sum = outG2cond$qsum_G2_cond[oml], mu.a = mu_sample[, oml], g.sum = outG2cond$gsum_G2_cond[, oml], p.new = outG2cond$pval_G2_cond[startcond:endcond], Score = outG2cond$Score_G2_cond[startcond:endcond], Phi = outG2cond$VarMat_G2_cond[, startcond:endcond], "SKAT-O")
        scaleFactorVec <- as.vector(G2condList$scaleFactor)
        G2condList$scaleFactorVec <- scaleFactorVec
        G2condList_list[[oml]] <- G2condList
        assign_conditionMarkers_factors_binary_region_multiTrait(scaleFactorVec, oml - 1)
      }
      # print(G2condList)
      # print(scaleFactorVec)
    }
  } else {
    condition_weights <- c(0)
  }

  # traitType = obj.model$traitType
  mu <- as.vector(t(I_mat) %*% (obj.model$mu))
  isgxe <- obj.model$isgxe
  rm(obj.model)
  gc()
  # print(gc(v=T))
  # if(file.exists(SAIGEOutputFile)) {print("ok 0 file exist")}


  # cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
  # cat("Number of markers in each chunk:\t", numLinesOutput, "\n")
  # cat("Number of chunks for all markers:\t", nChunks, "\n")
  # }
  # print("SAIGE.Marker -1")
  # print_g_n_unique()

  if (!isGroupTest) {
    OutputFile <- SAIGEOutputFile

    # if(file.exists(SAIGEOutputFile)) {print("ok 2 file exist")}
    if (!is.null(objGeno$markerInfo$CHROM)) {
      setorderv(objGeno$markerInfo, col = c("CHROM", "POS"))
    }

    SAIGE.Marker(
      traitType,
      phenotype_name_vec,
      genoType,
      objGeno$markerInfo$genoIndex_prev,
      objGeno$markerInfo$genoIndex,
      objGeno$markerInfo$CHROM,
      OutputFile,
      OutputFileIndex,
      markers_per_chunk,
      is_output_moreDetails,
      is_imputed_data,
      is_Firth_beta,
      LOCO,
      chrom,
      isCondition,
      is_overwrite_output,
      objGeno$anyInclude,
      isgxe
    )
  } else {
    MAFlimitMat <- NULL
    if (is.null(minMAF_in_groupTest_Exclude)) {
      minMAF_in_groupTest <- rep(0, length(maxMAF_in_groupTest))
    } else {
      minMAF_in_groupTest <- minMAF_in_groupTest_Exclude
    }
    MAFlimitMat <- cbind(minMAF_in_groupTest, maxMAF_in_groupTest)
    if (length(maxMAF_in_groupTest) == 1) {
      MAFlimitMat <- matrix(c(minMAF_in_groupTest, maxMAF_in_groupTest), ncol = 2, nrow = 1)
    }


    # cat("maxMAF_in_groupTest b ", maxMAF_in_groupTest, "\n")
    # cat("minMAF_in_groupTest b ", minMAF_in_groupTest, "\n")
    # cat("MAFlimitMat ", MAFlimitMat, "\n")

    maxMACbinind <- which(maxMAC_in_groupTest > 0)
    if (length(maxMACbinind) > 0) {
      maxMAC_in_groupTest_to_MAF <- (maxMAC_in_groupTest) / (2 * length(mu))
      cat("maxMAC_in_groupTest: ", maxMAC_in_groupTest, " is specified, corresponding to max MAF ", maxMAC_in_groupTest_to_MAF, "\n")
      for (i in 1:length(maxMACbinind)) {
        checkArgNumeric(maxMAC_in_groupTest_to_MAF[i], deparse(substitute(maxMAC_in_groupTest_to_MAF[i])), 0, 0.5, FALSE, TRUE)
      }
      # maxMAF_in_groupTest = unique(c(maxMAF_in_groupTest, maxMAC_in_groupTest_to_MAF))
      # maxMAF_in_groupTest = maxMAF_in_groupTest[order(maxMAF_in_groupTest)]
      # cat("max MAF cutoff ", maxMAF_in_groupTest, "will be applied\n")


      if (is.null(minMAC_in_groupTest_Exclude)) {
        minMAC_in_groupTest <- rep(0, length(maxMAC_in_groupTest))
      } else {
        minMAC_in_groupTest <- minMAC_in_groupTest_Exclude
        minMAC_in_groupTest_to_MAF <- minMAC_in_groupTest / (2 * length(mu))
        cat("minMAC_in_groupTest: ", minMAC_in_groupTest, " is specified, corresponding to min MAF ", minMAC_in_groupTest_to_MAF, "\n")
      }
      MAFlimitMat <- rbind(MAFlimitMat, cbind(minMAC_in_groupTest_to_MAF, maxMAC_in_groupTest_to_MAF))
    } # if(length(maxMACbinind) > 0){

    print(MAFlimitMat)
    MAFlimitMat <- MAFlimitMat[!duplicated(MAFlimitMat), , drop = F]
    print(MAFlimitMat)
    MAFlimitMat <- MAFlimitMat[order(MAFlimitMat[, 2], MAFlimitMat[, 1]), , drop = F]
    minMAF_in_groupTest <- MAFlimitMat[, 1]
    maxMAF_in_groupTest <- MAFlimitMat[, 2]
    cat("max MAF cutoff ", maxMAF_in_groupTest, "will be applied\n")
    cat("corresponding min MAF cutoff (exclude) ", minMAF_in_groupTest, "will be applied\n")


    # cat("MAFlimitMat b ", MAFlimitMat, "\n")


    # method_to_CollapseUltraRare,
    # DosageCutoff_for_UltraRarePresence,
    if (r.corr == 0 & is_Firth_beta) {
      print("WARNING: Note that the Firth correction has not been implemented for Burden test effect sizes when SKAT-O test is conducted. If the corrected Burden test effect sizes is needed, please use r.corr=1 to only conduct Burden test.")
      if (is_single_in_groupTest) {
        print("Firth correction will be used for effect sizes of single variant tests")
      }
    }

    BetaDist_weight_mat <- matrix(c(0, 0), ncol = 2)

    if (!is.null(weights.beta)) {
      BetaDist_weight_mat <- NULL
      for (i in 1:length(weights.beta)) {
        weightsbeta_val_vec <- as.numeric(unlist(strsplit(weights.beta[i], split = ",")))
        if (length(weightsbeta_val_vec) == 2) {
          BetaDist_weight_mat <- rbind(BetaDist_weight_mat, weightsbeta_val_vec)
        } else {
          stop("The ", i, "th element in weights.beta does not have 2 elements\n")
        }
      }
    }

    print("condition_weights")
    print(condition_weights)

    SAIGE.Region(
      mu,
      OutputFile,
      MACCutoff_to_CollapseUltraRare,
      groupFile,
      annotation_in_groupTest,
      maxMAF_in_groupTest,
      minMAF_in_groupTest,
      markers_per_chunk_in_groupTest,
      genoType,
      objGeno$markerInfo,
      traitType,
      phenotype_name_vec,
      is_imputed_data,
      isCondition,
      condition_weights,
      groups_per_chunk,
      r.corr,
      is_overwrite_output,
      is_single_in_groupTest,
      BetaDist_weight_mat,
      is_equal_weight_in_groupTest,
      is_output_markerList_in_groupTest,
      is_SKATO,
      chrom,
      is_fastTest,
      pval_cutoff_for_fastTest,
      is_output_moreDetails
    )
  }
}
