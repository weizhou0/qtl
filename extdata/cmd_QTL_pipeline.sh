genename=$1
windowsize=$2 

#1000000 1Mb
cellType="B_IN"
traitType="count"
#genename=SH3BGRL3

outpath=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/gene_group/output/
phenofile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/B_IN_pheno_cov.newheader.nonzero_1perc.tsv
geneLocationFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/GeneLocations.renamegene.tsv



###step 1
#phenofile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/B_IN_pheno_cov_weirdgenes_doubletQC.tsv
#step1outpath=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/
#step1prefix=${step1outpath}${genename}_${cellType}_${traitType}_uncond_test2

step1prefix=${outpath}${genename}_${cellType}_${traitType}_step1

#: '


Rscript /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1_fitNULLGLMM_qtl_new.R    \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=${phenofile}        \
        --phenoCol=${genename}  \
        --covarColList=age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pf1,pf2  \
        --sampleCovarColList=age,sex,pc1,pc2,pc3,pc4,pc5,pc6    \
        --sampleIDColinphenoFile=individual \
        --traitType=count \
	--outputPrefix=${step1prefix}	\
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/raw_genotype/OneK1K_AllChr_pruned_10k_randomMarkers_forVR  \
        --IsOverwriteVarianceRatioFile=TRUE

#'

i=$(awk -v gene=$genename '$1 == gene {print $3}' $geneLocationFile)

echo "$i"
regionFile=${outpath}Gene_${genename}_${windowsize}.region

awk -v gene=$genename -v chr=$i -v wd=$windowsize '$1 == gene {print $1" "$3" "$4-wd" "$5+wd}' ${geneLocationFile} | awk '$3 < 0 {$3 = 1}1' > ${regionFile}

groupFile=${outpath}Gene_${genename}_${windowsize}.grp

#awk -v chr=$i -v wd=$windowsize '$3 == chr {print $1" "$3" "$4-wd" "$5+wd}' ${geneLocationFile} | awk '$3 < 0 {$3 = 1}1' > ${regionFile}

head "${regionFile}"

#: '
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
        --weights.beta="1;25,1;1"   \
        --is_SKAT=FALSE


#        --outputPrefix=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/trans_region_list/bychr_grpfile/dose.filtered.R2_0.8.GRCh37_hg19_region_chr${i}_grpfile
        
	
	
	#--vcfFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz  \
        #--vcfFileIndex=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz.csi  \
