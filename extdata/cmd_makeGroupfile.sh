i=1
Rscript makeGroupFile.R	\
	--vcfFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz  \
        --vcfFileIndex=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz.csi  \
	--regionFile=./region.txt	\
	--outputPrefix=./regionmake_grp
