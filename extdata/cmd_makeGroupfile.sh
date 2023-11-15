i=$1
Rscript makeGroupFile.R	\
	--vcfFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz  \
        --vcfFileIndex=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz.csi  \
	--regionFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/trans_region_list/bychr/blocks_s2500_m25_f1_w200.GRCh37_hg19_region_chr${i}	\
	--outputPrefix=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/trans_region_list/bychr_grpfile/dose.filtered.R2_0.8.GRCh37_hg19_region_chr${i}_grpfile
