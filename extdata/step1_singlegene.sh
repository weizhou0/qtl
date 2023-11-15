
    i=$1
    echo $i
    step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_${i}
    Rscript step1_fitNULLGLMM_qtl.R         --useSparseGRMtoFitNULL=FALSE          --useGRMtoFitNULL=FALSE         --phenoFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt         --phenoCol=gene_${i}               --covarColList=X1,X2,pf1,pf2            --sampleCovarColList=X1,X2              --sampleIDColinphenoFile=IND_ID         --traitType=count         --outputPrefix=${step1prefix}           --skipVarianceRatioEstimation=FALSE          --isRemoveZerosinPheno=FALSE         --isCovariateOffset=FALSE          --isCovariateTransform=TRUE          --skipModelFitting=FALSE          --tol=0.00001           --plinkFile=./input/n.indep_100_n.cell_1_01.step1               --IsOverwriteVarianceRatioFile=TRUE 
