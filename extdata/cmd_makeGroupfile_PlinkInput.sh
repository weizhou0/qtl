i=$1 ##chromosome
windowsize=$2 #1000000 1Mb

outpath=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/gene_group/

geneLocationFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/GeneLocations.renamegene.tsv
regionFile=${outpath}Gene_chr${i}_${windowsize}.region
groupFile=${outpath}Gene_chr${i}_${windowsize}.grp

awk -v chr=$i -v wd=$windowsize '$3 == chr {print $1" "$3" "$4-wd" "$5+wd}' ${geneLocationFile} | awk '$3 < 0 {$3 = 1}1' > ${regionFile}


Rscript makeGroupFile.R \
        --bedFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bed \
        --bimFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.bim \
        --famFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/computationCost/GMMAT/full/full_genome_chr${i}.fam \
	--regionFile=${regionFile}	\
	--outputPrefix=${groupFile}
#        --outputPrefix=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/input_files/step2/trans_region_list/bychr_grpfile/dose.filtered.R2_0.8.GRCh37_hg19_region_chr${i}_grpfile
        
	
	
	#--vcfFile=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz  \
        #--vcfFileIndex=/humgen/atgu1/fin/wzhou/projects/eQTL_method_dev/realdata/oneK1K/AnnaCuomo_Yavar/genotype/filter_vcf_r08/chr${i}.dose.filtered.R2_0.8.vcf.gz.csi  \
