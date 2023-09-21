i=1

step1path=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/
prefix11=${step1path}ISG15_B_IN_count_uncond_test2

##test all variants with MAF >= 5 on chromosome 1 for the RNA experssion of one gene ISG15

Rscript /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/step2_tests_qtl_new.R     \
        --bedFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bed \
        --bimFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bim \
        --famFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.fam \
        --SAIGEOutputFile=./B_IN_count_uncon_variant    \
        --chrom=${i}       \
        --minMAF=0 \
        --minMAC=5 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${prefix11}.rda        \
        --varianceRatioFile=${prefix11}.varianceRatio.txt       \
        --markers_per_chunk=10000       \
        --SPAcutoff=2


##test all variants with MAF >= 5 on chromosome 1 for the RNA experssion of 50 genes
Rscript /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/tool_dev/qtl/extdata/step2_tests_qtl_new.R     \
        --bedFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bed \
        --bimFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bim \
        --famFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.fam \
        --SAIGEOutputFile=./B_IN_count_uncon_variant    \
        --chrom=${i}       \
        --minMAF=0 \
        --minMAC=5 \
        --LOCO=FALSE    \
	--GMMATmodel_varianceRatio_multiTraits_File=./traitModelVRfile.txt.50	\
        --markers_per_chunk=10000       \
        --SPAcutoff=2


#--GMMATmodel_varianceRatio_multiTraits_File: file containing 3 columns: phenotype name, model file, and variance ratio file. Each line is for one phenotype


##head ./traitModelVRfile.txt.50
#AP006222.2 /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/AP006222.2_B_IN_count_uncond_test2.rda /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/AP006222.2_B_IN_count_uncond_test2.varianceRatio.txt
#LINC00115 /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/LINC00115_B_IN_count_uncond_test2.rda /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/LINC00115_B_IN_count_uncond_test2.varianceRatio.txt
#FAM41C /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/FAM41C_B_IN_count_uncond_test2.rda /humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step1/output/FAM41C_B_IN_count_uncond_test2.varianceRatio.txt

