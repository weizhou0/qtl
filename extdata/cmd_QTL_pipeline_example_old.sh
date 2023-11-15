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
echo -e "2\t100001\t310001" > ${regionFile}

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

Rscript step3_gene_pvalue.R	\
	--assocFile=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis	\
	--geneName=gene_1	\
	--genePval_outputFile=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_genePval




##trans-eQTL analysis for multiple single genes


### step 1 can be run independently for each gene, so each gene job can use one CPU. 

for i in {1..100}
do 
    echo $i
    step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_${i}
    Rscript step1_fitNULLGLMM_qtl.R \
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
        --IsOverwriteVarianceRatioFile=TRUE
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




Rscript /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/makeGroupFile.R \
        --bedFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bed \
        --bimFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bim \
        --famFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.fam \
	--regionFile=${regionFile}	\
	--outputPrefix=${groupFile}
#'

prefix=${outpath}${genename}_${cellType}_${traitType}_uncond_${windowsize}

Rscript /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2_tests_qtl_new.R      \
        --bedFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bed \
        --bimFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bim \
        --famFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.fam \
        --SAIGEOutputFile=${prefix}	\
        --chrom=${i}       \
        --minMAF=0 \
        --minMAC=5 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda	\
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt	\
        --groupFile=${groupFile}	\
        --annotation_in_groupTest=null  \
        --maxMAF_in_groupTest=0.1,0.5   \
        --minMAF_in_groupTest_Exclude=0,0       \
        --is_fastTest=FALSE      \
        --markers_per_chunk=10000       \
        --pval_cutoff_for_fastTest=0.000000000005 \
        --SPAcutoff=2 \
        --is_noadjCov=TRUE      \
        --is_single_in_groupTest=FALSE  \
        --is_equal_weight_in_groupTest=FALSE	\
	--groups_per_chunk=1    \
        --weights.beta="1;25" 


#        --outputPrefix=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/trans_region_list/bychr_grpfile/dose.filtered.R2_0.8.GRCh37_hg19_region_chr${i}_grpfile
        
	
	
	#--vcfFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz  \
        #--vcfFileIndex=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz.csi  \
