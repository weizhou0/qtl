#if there is no sample relatedness in the data set, GRM is not needed to fit the null model
##--useGRMtoFitNULL=FALSE
#if there is sample relatedness in the data set, use a sparse GRM
##--sparseGRMFile=nfam_5_nindep_0.mtx     \
##--sparseGRMSampleIDFile=nfam_5_nindep_0.mtx.sampleID    \
##--useSparseGRMtoFitNULL=TRUE    \


#we always use a plink file that contains 2000 independent markers (LD pruned) to estimate a variance ratio 
##--plinkFile=

Rscript step1_fitNULLGLMM_qtl.R	\
	--sparseGRMFile=./input/nfam_5_nindep_0.mtx	\
	--sparseGRMSampleIDFile=./input/nfam_5_nindep_0.mtx.sampleID	\
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./input/seed_1_nfam_5_nindep_0_ncell_100_lambda_50_withCov_Poisson.txt	\
        --phenoCol=y      \
        --covarColList=PC1	\
        --sampleIDColinphenoFile=IND_ID        \
        --traitType=count       \
        --outputPrefix=./output/poisson_test        \
        --skipVarianceRatioEstimation=FALSE      \
        --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE     \
        --skipModelFitting=FALSE        \
        --plinkFile=./input/nfam_500_nindep_5000_step1-nodupSNPs	\
        --IsOverwriteVarianceRatioFile=TRUE &> ./output/poisson_test.log

Rscript step1_fitNULLGLMM_qtl.R \
        --useGRMtoFitNULL=FALSE    \
        --phenoFile=./input/seed_1_nfam_5_nindep_0_ncell_100_lambda_50_withCov_Poisson.txt      \
        --phenoCol=y      \
        --covarColList=PC1      \
        --sampleIDColinphenoFile=IND_ID        \
        --traitType=count       \
        --outputPrefix=./output/poisson_test        \
        --skipVarianceRatioEstimation=FALSE      \
        --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE     \
        --skipModelFitting=FALSE        \
        --plinkFile=./input/nfam_500_nindep_5000_step1-nodupSNPs        \
        --IsOverwriteVarianceRatioFile=TRUE &> ./output/poisson_test.log





#Run a single variant association tests
Rscript step2_tests_qtl.R	\
        --bedFile=./input/nfam_500_nindep_5000_01.bed	\
        --bimFile=./input/nfam_500_nindep_5000_01.bim	\
	--famFile=./input/nfam_500_nindep_5000_01.fam	\
        --AlleleOrder=alt-first \
        --SAIGEOutputFile=./output/nfam_500_nindep_5000_01_poisson_single    \
        --chrom=1       \
        --minMAF=0 \
        --minMAC=5 \
        --LOCO=FALSE    \
        --GMMATmodelFile=./output/poisson_test.rda  \
	--varianceRatioFile=./output/poisson_test.varianceRatio.txt	\
        --is_noadjCov=TRUE      \
        --is_sparseGRM=FALSE    \
        --markers_per_chunk=10000       \
        --pval_cutoff_for_fastTest=0.05 \
        --SPAcutoff=10000 \
        --is_EmpSPA=FALSE &>  ./output/nfam_500_nindep_5000_01_poisson_single.log


#Run a group tests
Rscript step2_tests_qtl.R       \
        --bedFile=./input/nfam_500_nindep_5000_01.bed   \
        --bimFile=./input/nfam_500_nindep_5000_01.bim   \
        --famFile=./input/nfam_500_nindep_5000_01.fam   \
	--AlleleOrder=alt-first \
        --SAIGEOutputFile=./output/nfam_500_nindep_5000_01_poisson_set    \
        --chrom=1       \
        --minMAF=0 \
        --minMAC=5 \
        --LOCO=FALSE    \
        --GMMATmodelFile=./output/poisson_test.rda  \
        --varianceRatioFile=./output/poisson_test.varianceRatio.txt     \
        --is_noadjCov=TRUE      \
        --is_sparseGRM=FALSE    \
        --markers_per_chunk=10000       \
        --pval_cutoff_for_fastTest=0.05 \
        --SPAcutoff=10000 \
        --is_EmpSPA=FALSE \
	--annotation_in_groupTest=lof:missense  \
        --maxMAF_in_groupTest=0.5       \
	--groupFile=./input/nfam_500_nindep_5000_01_group.txt &> ./output/nfam_500_nindep_5000_01_poisson_set.log 
