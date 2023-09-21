i=1

step1path=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/
prefix10=${step1path}ZRANB1_B_IN_count_uncond_test2

#test the groups in the gorup file for one gene's RNA expression
Rscript /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/step2_tests_qtl_new.R      \
        --bedFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bed \
        --bimFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bim \
        --famFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.fam \
        --SAIGEOutputFile=./multiTrait.output_chr1_grpfile_1    \
        --chrom=${i}       \
        --minMAF=0 \
        --minMAC=5 \
        --LOCO=FALSE    \
        --annotation_in_groupTest=null  \
        --maxMAF_in_groupTest=0.1,0.5   \
        --minMAF_in_groupTest_Exclude=0,0       \
        --markers_per_chunk=10000       \
        --SPAcutoff=2 \
        --is_equal_weight_in_groupTest=TRUE    \
        --groups_per_chunk=1    \
        --weights.beta="1;25"   \
        --GMMATmodelFile=${prefix10}.rda        \
        --varianceRatioFile=${prefix10}.varianceRatio.txt       \
        --groupFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/trans_region_list/bychr_grpfile/dose.filtered.R2_0.8.GRCh37_hg19_region_chr1_grpfile_1_first2

##group file example
#1_10539_1376204 var 1:715265:C:T 1:715367:A:G 1:717485:C:A
#1_10539_1376204 anno null null null
#1_1376205_2215495 var 1:1376214:C:G 1:1376567:A:G 1:1377151:G:T
#1_1376205_2215495 anno null null null

##group file can contain user specified weights
#1_10539_1376204 var 1:715265:C:T 1:715367:A:G 1:717485:C:A
#1_10539_1376204 anno null null null
#1_10539_1376204 weight:a 1 2 3
#1_10539_1376204 weight:b 4 5 6

#test the groups in the gorup file for 50 gene's RNA expression
Rscript /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/step2_tests_qtl_new.R      \
        --bedFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bed \
        --bimFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bim \
        --famFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.fam \
        --SAIGEOutputFile=./multiTrait.output_chr1_grpfile_1    \
        --chrom=${i}       \
        --minMAF=0 \
        --minMAC=5 \
        --LOCO=FALSE    \
        --annotation_in_groupTest=null  \
        --maxMAF_in_groupTest=0.1,0.5   \
        --minMAF_in_groupTest_Exclude=0,0       \
        --markers_per_chunk=10000       \
        --SPAcutoff=2 \
        --is_equal_weight_in_groupTest=TRUE    \
        --groups_per_chunk=1    \
        --weights.beta="1;25"   \
	--GMMATmodel_varianceRatio_multiTraits_File=./traitModelVRfile.txt.50	\
        --groupFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/trans_region_list/bychr_grpfile/dose.filtered.R2_0.8.GRCh37_hg19_region_chr1_grpfile_1_first2
