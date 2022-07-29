    Rscript step1_fitNULLGLMM_1.1.4.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --qCovarColList=x2  \
        --sampleIDColinphenoFile=IID \
        --invNormalize=TRUE     \
        --traitType=quantitative        \
        --outputPrefix=./output/example_quantitative_1.1.4 \
        --nThreads=24	\
        --IsOverwriteVarianceRatioFile=TRUE


    Rscript step1_fitNULLGLMM_1.1.4.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --qCovarColList=x2  \
        --sampleIDColinphenoFile=IID \
        --invNormalize=TRUE     \
        --traitType=quantitative        \
        --outputPrefix=./output/example_quantitative_multiV_1.1.4 \
        --nThreads=24   \
        --IsOverwriteVarianceRatioFile=TRUE	\
	--VmatFilelist="bdiagMM.mtx"	\
	--VmatSampleFilelist="bdiagMM_sampleID.txt"


    Rscript step1_fitNULLGLMM_1.1.4.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --qCovarColList=x2  \
        --sampleIDColinphenoFile=IID \
        --invNormalize=TRUE     \
        --traitType=quantitative        \
        --outputPrefix=./output/example_quantitative_multiV_1.1.4 \
        --nThreads=24   \
        --IsOverwriteVarianceRatioFile=TRUE     \
        --VmatFilelist="bdiagMM.mtx"    \
        --VmatSampleFilelist="bdiagMM_sampleID.txt"




        Rscript step1_fitNULLGLMM_1.1.4.R     \
	        --sparseGRMFile=./GMMAT_GRM.mtx	\
		--sparseGRMSampleIDFile=./GMMAT_sampleID.txt	\
        --useSparseGRMtoFitNULL=TRUE    \
	--phenoFile=./GMMAT_pheno.txt	\
	--phenoCol=disease	\
	--covarColList=age,sex	\
	--sampleIDColinphenoFile=id	\
	--traitType=binary	\
	--outputPrefix=GMMAT_binary_multiV_1.1.4	\
	--VmatFilelist="GMMAT_M10.mtx"	\
	--VmatSampleFilelist="./GMMAT_sampleID.txt"	\
	--skipVarianceRatioEstimation=TRUE	\
	--isCovariateOffset=FALSE	\
	--isCovariateTransform=FALSE &> step1_multiV.log



	 Rscript step1_fitNULLGLMM_1.1.3.R     \
 --sparseGRMFile=./GMMAT_GRM.mtx \
                --sparseGRMSampleIDFile=./GMMAT_sampleID.txt    \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./GMMAT_pheno.txt   \
        --phenoCol=disease      \
        --covarColList=age,sex  \
        --sampleIDColinphenoFile=id     \
        --traitType=binary      \
        --outputPrefix=GMMAT_binary_multiV_1.1.4        \
        --skipVarianceRatioEstimation=TRUE      \
       --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE


        Rscript step1_fitNULLGLMM_1.1.4.R     \
                --sparseGRMFile=./GMMAT_GRM.mtx \
                --sparseGRMSampleIDFile=./GMMAT_sampleID.txt    \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./GMMAT_pheno.txt   \
        --phenoCol=disease      \
        --covarColList=age,sex  \
        --sampleIDColinphenoFile=id     \
        --traitType=binary      \
        --outputPrefix=GMMAT_binary_multiV_1.1.4        \
        --skipVarianceRatioEstimation=TRUE	\
       --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE	&> step1_multiV_singleV.log


        Rscript step1_fitNULLGLMM_1.1.4.R     \
                --sparseGRMFile=./GMMAT_GRM.mtx \
                --sparseGRMSampleIDFile=./GMMAT_sampleID.txt    \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./GMMAT_pheno_quantitative.txt   \
        --phenoCol=y.repeated      \
        --covarColList=sex  \
        --sampleIDColinphenoFile=id     \
        --traitType=quantitative      \
	--invNormalize=FALSE	\
        --outputPrefix=GMMAT_quantitative_1.1.4        \
        --skipVarianceRatioEstimation=TRUE      \
       --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE    &> step1_singleV_quantitative.log


	        Rscript step1_fitNULLGLMM_1.1.4.R     \
                --sparseGRMFile=./GMMAT_GRM.mtx \
                --sparseGRMSampleIDFile=./GMMAT_sampleID.txt    \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./GMMAT_pheno_quantitative_repeat.txt   \
        --phenoCol=y.trend      \
	--longlCol=time	\
        --covarColList=sex  \
        --sampleIDColinphenoFile=id     \
        --traitType=quantitative      \
        --invNormalize=FALSE    \
        --outputPrefix=GMMAT_quantitative_1.1.4        \
        --skipVarianceRatioEstimation=TRUE      \
       --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE




 Rscript step1_fitNULLGLMM_1.1.4.R                     --sparseGRMFile=./GMMAT_GRM.mtx                 --sparseGRMSampleIDFile=./GMMAT_sampleID.txt            --useSparseGRMtoFitNULL=TRUE            --phenoFile=./GMMAT_pheno_quantitative.txt           --phenoCol=y.repeated              --covarColList=sex          --sampleIDColinphenoFile=id             --traitType=quantitative              --invNormalize=FALSE            --outputPrefix=GMMAT_quantitative_1.1.4                --skipVarianceRatioEstimation=TRUE             --isCovariateOffset=FALSE               --isCovariateTransform=FALSE --nrun=100 &> step1_fitNULLGLMM_1.1.4.log





        Rscript step1_fitNULLGLMM_1.1.4.R     \
        --sparseGRMFile=./GMMAT_GRM.mtx \
        --sparseGRMSampleIDFile=./GMMAT_sampleID.txt    \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./GMMAT_pheno_quantitative.txt   \
        --phenoCol=y.repeated      \
        --covarColList=sex  \
        --sampleIDColinphenoFile=id     \
        --traitType=quantitative      \
        --invNormalize=FALSE    \
        --outputPrefix=GMMAT_quantitative_multiV_1.1.4        \
        --skipVarianceRatioEstimation=TRUE      \
        --VmatFilelist="GMMAT_M10.mtx"  \
        --VmatSampleFilelist="./GMMAT_sampleID.txt"     \
       --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE    &> step1_singleV_quantitative.log



 Rscript step1_fitNULLGLMM_1.1.4.R     \
                --sparseGRMFile=./GMMAT_GRM.mtx \
                --sparseGRMSampleIDFile=./GMMAT_sampleID.txt    \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./GMMAT_pheno_quantitative_repeat.txt   \
        --phenoCol=y.repeated      \
        --covarColList=sex  \
        --sampleIDColinphenoFile=id     \
        --traitType=quantitative      \
        --invNormalize=FALSE    \
        --outputPrefix=GMMAT_quantitative_multiV_1.1.4        \
        --skipVarianceRatioEstimation=TRUE      \
        --VmatFilelist="I.mtx"  \
        --VmatSampleFilelist="./GMMAT_sampleID.txt"     \
       --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE    &> step1_repeat_quantitative.log


  Rscript step1_fitNULLGLMM_1.1.4.R     \
                --sparseGRMFile=./GMMAT_GRM.mtx \
                --sparseGRMSampleIDFile=./GMMAT_sampleID.txt    \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./GMMAT_pheno_quantitative_repeat.txt   \
        --phenoCol=y.repeated      \
        --covarColList=sex  \
        --sampleIDColinphenoFile=id     \
        --traitType=quantitative      \
        --invNormalize=FALSE    \
        --outputPrefix=GMMAT_quantitative_multiV_1.1.4        \
        --skipVarianceRatioEstimation=TRUE      \
       --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE



 Rscript step1_fitNULLGLMM_1.1.3.R     \
                --sparseGRMFile=./GMMAT_GRM.mtx \
                --sparseGRMSampleIDFile=./GMMAT_sampleID.txt    \
        --useSparseGRMtoFitNULL=TRUE    \
        --phenoFile=./GMMAT_pheno_quantitative.txt   \
        --phenoCol=y.repeated      \
        --covarColList=sex  \
        --sampleIDColinphenoFile=id     \
        --traitType=quantitative      \
        --invNormalize=FALSE    \
        --outputPrefix=GMMAT_quantitative_1.1.3        \
        --skipVarianceRatioEstimation=TRUE      \
       --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE




	                Rscript step1_fitNULLGLMM_1.1.4.R     \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
        --useSparseGRMtoFitNULL=TRUE    \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr_random1000 \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --skipVarianceRatioEstimation=FALSE \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --qCovarColList=x2  \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_sparseGRM_vr_1.1.4 \
        --IsOverwriteVarianceRatioFile=TRUE
