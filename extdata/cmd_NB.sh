        Rscript step1_fitNULLGLMM_1.1.3_Poisson.R	\
	--useGRMtoFitNULL=FALSE	\	
        --phenoFile=/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_newgit/SAIGE_1.1.3/NB_test/NBZIMM_pheno.txt   \
        --phenoCol=y      \
        --covarColList=Days,Age,Race,preg  \
        --sampleIDColinphenoFile=subject     \
        --traitType=count_nb      \
        --outputPrefix=NB_1.1.4_new        \
        --skipVarianceRatioEstimation=TRUE      \
        --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE 	\
	--offsetCol=logN
	
	
	&> step1_multiV.log





	        Rscript step1_fitNULLGLMM_1.1.3_Poisson.R       \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/SAIGE_newgit/SAIGE_1.1.3/NB_test/NBZIMM_pheno.txt   \
        --phenoCol=y      \
        --covarColList=Days,Age,Race,preg  \
        --sampleIDColinphenoFile=subject     \
        --traitType=count_nb      \
        --outputPrefix=NB_1.1.4_new_092022       \
        --skipVarianceRatioEstimation=TRUE      \
        --isCovariateOffset=FALSE       \
        --isCovariateTransform=FALSE	\
	--offsetCol=logN
