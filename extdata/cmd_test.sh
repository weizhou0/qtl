   Rscript step1_fitNULLGLMM.R     \
   Rscript step1_fitNULLGLMM_1.1.3_Poisson.R	\
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --useSparseGRMtoFitNULL=FALSE    \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --qCovarColList=x2  \
        --sampleIDColinphenoFile=IID \
        --invNormalize=TRUE     \
        --traitType=quantitative        \
        --outputPrefix=./output/example_quantitative_test \
        --nThreads=24	\
        --IsOverwriteVarianceRatioFile=TRUE	\
	--useGRMtoFitNULL=TRUE


       Rscript step2_SPAtests.R        \
        --vcfFile=./input/genotype_100markers.vcf.gz    \
        --vcfFileIndex=./input/genotype_100markers.vcf.gz.csi     \
        --vcfField=GT   \
	--chrom=1	\
        --SAIGEOutputFile=./output/genotype_100markers_marker_vcf_test.txt \
        --minMAF=0 \
        --minMAC=20 \
	--LOCO=FALSE	\
        --GMMATmodelFile=./output/example_quantitative_test.rda \
        --is_output_moreDetails=TRUE	\
        --varianceRatioFile=./output/example_quantitative_test.varianceRatio.txt
