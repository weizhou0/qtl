##The example data contain 100 cells for each of the 100 individuals without sample relatedness


##cis-eQTL analysis for a single gene

step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1

Rscript step1_fitNULLGLMM_qtl.R	\
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt	\
        --phenoCol=gene_1	\
        --covarColList=X1,X2,pf1,pf2	\
        --sampleCovarColList=X1,X2	\
        --sampleIDColinphenoFile=IND_ID \
        --traitType=count \
	--outputPrefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1	\
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=./input/n.indep_100_n.cell_1_01.step1	\
        --IsOverwriteVarianceRatioFile=TRUE

step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1
step2prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis
regionFile=./input/gene_1_cis_region.txt
echo -e "2\t300001\t610001" > ${regionFile}

Rscript step2_tests_qtl.R	\
	--bedFile=./input/n.indep_100_n.cell_1.bed	\
	--bimFile=./input/n.indep_100_n.cell_1.bim	\
	--famFile=./input/n.indep_100_n.cell_1.fam	\
        --SAIGEOutputFile=${step2prefix}     \
        --chrom=2       \
        --minMAF=0 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --rangestoIncludeFile=${regionFile}     \
        --markers_per_chunk=10000

Rscript step3_gene_pvalue_qtl.R	\
	--assocFile=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis	\
	--geneName=gene_1	\
	--genePval_outputFile=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_genePval


##test rare variants eQTLs
##make group file for the gene
groupFile=${regionFile}.grp
Rscript makeGroupFile.R	\
	--bedFile=./input/n.indep_100_n.cell_1.bed	\
	--bimFile=./input/n.indep_100_n.cell_1.bim	\
	--famFile=./input/n.indep_100_n.cell_1.fam	\
	--regionFile=${regionFile}	\
	--outputPrefix=${regionFile}.grp

Rscript step2_tests_qtl.R       \
        --bedFile=./input/n.indep_100_n.cell_1.bed      \
        --bimFile=./input/n.indep_100_n.cell_1.bim      \
        --famFile=./input/n.indep_100_n.cell_1.fam      \
        --SAIGEOutputFile=${step2prefix}_rare     \
        --chrom=2       \
        --maxMAF_in_groupTest=0.1   \
        --minMAF_in_groupTest_Exclude=0       \
	--groupFile=${groupFile}	\
	--annotation_in_groupTest=null  \
	--MACCutoff_to_CollapseUltraRare=10	\
	--is_single_in_groupTest=TRUE  \
	--is_equal_weight_in_groupTest=TRUE    \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --markers_per_chunk=10000



##trans-eQTL analysis for multiple genes


### step 1 can be run independently for each gene, so each gene job can use one CPU. 

for i in {1..100}
do 
    echo $i
    step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_${i}
    /bin/time -o ${step1prefix}.runinfo.txt -v Rscript step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt \
        --phenoCol=gene_${i}       \
        --covarColList=X1,X2,pf1,pf2    \
        --sampleCovarColList=X1,X2      \
        --sampleIDColinphenoFile=IND_ID \
        --traitType=count \
        --outputPrefix=${step1prefix}	\
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=./input/n.indep_100_n.cell_1_01.step1       \
        --IsOverwriteVarianceRatioFile=TRUE &> ${step1prefix}.log
done

##to use muliple CPUs for parallel computing. 
##1. create a bash file 

echo -e "
    i=\$1
    echo \$i
    step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_\${i}
    Rscript step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/input/seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt \
        --phenoCol=gene_\${i}       \
        --covarColList=X1,X2,pf1,pf2    \
        --sampleCovarColList=X1,X2      \
        --sampleIDColinphenoFile=IND_ID \
        --traitType=count \
        --outputPrefix=\${step1prefix}   \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=./input/n.indep_100_n.cell_1_01.step1       \
        --IsOverwriteVarianceRatioFile=TRUE " | tee step1_singlegene.sh

##2. write a job file with each line for one gene
touch step1_singlegene.job
rm step1_singlegene.job 
for i in {1..100}
do
  echo -e "bash step1_singlegene.sh $i" >> step1_singlegene.job
done 

##3. if there are multiple CPUs available, use them for step 1 jobs, with each for a gene
##e.g. rum 30 jobs simutaneously

cat step1_singlegene.job | xargs -I CMD --max-procs=30 bash -c CMD


##step 2: test eQTLs on chromosome 2  for all 100 genes
#use --GMMATmodel_varianceRatio_multiTraits_File
#first create the file containing 3 columns: phenotype name, model file, and variance ratio file. Each line is for one phenotype. This file is used when multiple phenotypes are analyzed simutaneously 

touch ./input/step1_output_formultigenes.txt
rm ./input/step1_output_formultigenes.txt
for i in {1..100}
do
  step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_${i}
  echo -e "gene_${i} ${step1prefix}.rda ${step1prefix}.varianceRatio.txt" >> ./input/step1_output_formultigenes.txt
done


step2prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_chr2
Rscript step2_tests_qtl.R       \
        --bedFile=./input/n.indep_100_n.cell_1.bed      \
        --bimFile=./input/n.indep_100_n.cell_1.bim      \
        --famFile=./input/n.indep_100_n.cell_1.fam      \
        --SAIGEOutputFile=${step2prefix}     \
	--GMMATmodel_varianceRatio_multiTraits_File=./input/step1_output_formultigenes.txt	\
        --chrom=2       \
        --minMAF=0.05 \
        --LOCO=FALSE    \
        --SPAcutoff=2 \
        --markers_per_chunk=10000
