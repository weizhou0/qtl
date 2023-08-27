#' Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.
#'
#' @param bgenFile character. Path to bgen file. Currently version 1.2 with 8 bit compression is supported
#' @param bgenFileIndex character. Path to the .bgi file (index of the bgen file)
#' @param sampleFile character. Path to the file that contains one column for IDs of samples in the bgen file. The file does not contain header lines.
#' @param vcfFile character. Path to vcf file
#' @param vcfFileIndex character. Path to vcf index file. Indexed by tabix. Path to index for vcf file by tabix, .csi file using 'tabix --csi -p vcf file.vcf.gz'
#' @param vcfField character. genotype field in vcf file to use. "DS" for dosages or "GT" for genotypes. By default, "DS".
#' @param savFile character. Path to sav file
#' @param savFileIndex character. Path to index for sav file .s1r
#' @param bedFile character. Path to bed file (PLINK)
#' @param bimFile character. Path to bim file (PLINK)
#' @param famFile character. Path to fam file (PLINK)
#' @param AlleleOrder character. alt-first or ref-first for bgen or PLINK files. By default, alt-first
#' @param regionFile character. Path to a file containing genome regions to extract genetic markers for making group files
#' @param outputPrefix character. Path to the output files with prefix
#' @return one or multiple group files, each for one chromosome
#' @export
makeGroupFileforRegions = function(bgenFile = "",
                 bgenFileIndex = "",
                 vcfFile = "",
                 vcfFileIndex = "",
                 vcfField = "DS",
                 savFile = "",
                 savFileIndex = "",
                 bedFile="",
                 bimFile="",
                 famFile="",
                 AlleleOrder = "alt-first", #new
		 regionFile="",
		 outputPrefix=""){

	Check_File_Exist(regionFile, "regionFile")
	regionData = extractRegions(regionFile)	
	dosageFileType = checkGenoInput(bgenFile = bgenFile,
                 bgenFileIndex = bgenFileIndex,
                 vcfFile = vcfFile,
                 vcfFileIndex = vcfFileIndex,
                 vcfField = vcfField,
                 savFile = savFile,
                 savFileIndex = savFileIndex,
                 bedFile = bedFile,
                 bimFile = bimFile,
                 famFile = famFile)

 if(dosageFileType == "plink"){
    if(is.null(AlleleOrder)) AlleleOrder = "alt-first"

    cat("allele order in the plink file is ", AlleleOrder, ".\n")

    if(bimFile == ""){
        bimFile = gsub("bed$", "bim", bedFile)
    }
    if(famFile == ""){
        famFile = gsub("bed$", "fam", bedFile)
    }
    markerInfo = data.table::fread(bimFile, header = F, select = c(1, 2, 4, 5, 6))
    #markerInfo = as.data.frame(markerInfo)
    #markerInfo is a data.table
    if(AlleleOrder == "alt-first")
      names(markerInfo) = c("CHROM", "ID", "POS", "ALT", "REF")
      #markerInfo = markerInfo[,c(1,2,3,5,4)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    if(AlleleOrder == "ref-first")
      names(markerInfo) = c("CHROM", "ID", "POS", "REF", "ALT")
      markerInfo$genoIndex = 1:nrow(markerInfo) - 1  # -1 is to convert 'R' to 'C++'
      markerInfo$genoIndex_prev = rep(0, length(markerInfo$genoIndex))
    #if(chrom != ""){
    #  markerInfo = markerInfo[which(markerInfo[,1] == chrom), ]
    #}
    markerInfo$ID2 = paste0(markerInfo$CHROM,":",markerInfo$POS, ":", markerInfo$REF,":", markerInfo$ALT)
    markerInfo[,REF:=NULL]
    markerInfo[,ALT:=NULL]


    makeGroupFilewithRegionData(markerInfo, regionData, outputPrefix)

  }	

  if(dosageFileType == "bgen"){


    if(is.null(AlleleOrder)) AlleleOrder = "ref-first"

    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgenFileIndex)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    markerInfo = dplyr::tbl(db_con, "Variant")
    markerInfo = as.data.table(markerInfo, keep.rownames = FALSE)
    rmcol = setdiff(names(markerInfo), c("chromosome", "position", "rsid", "allele1", "allele2", 'file_start_position', 'size_in_bytes'))
    markerInfo = markerInfo[, !..rmcol]
    if(AlleleOrder == "alt-first")
      names(markerInfo) = c("CHROM", "POS", "ID", "ALT", "REF", "genoIndex", "size_in_bytes")
      #markerInfo = markerInfo[,c(1,2,3,6,5,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    if(AlleleOrder == "ref-first")
      names(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT", "genoIndex", "size_in_bytes")
    markerInfo$ID2 = paste0(markerInfo$CHROM,":", markerInfo$POS ,":", markerInfo$REF, ":", markerInfo$ALT)
    markerInfo$genoIndex_prev = markerInfo$genoIndex + markerInfo$size_in_bytes
    markerInfo[,REF:=NULL]
    markerInfo[,ALT:=NULL]
    markerInfo[,size_in_bytes:=NULL]
  
        makeGroupFilewithRegionData(markerInfo, regionData, outputPrefix)

}

 if(dosageFileType == "vcf"){

    if(vcfFile != ""){
        vcfFileIndex = paste(vcfFile, ".csi", sep = "")
    }else if(savFile != ""){
        vcfFile = savFile
        vcfFileIndex = paste(vcfFile, ".s1r", sep = "")
    }

      sampleInModel = as.character()
      CHROM1 = regionData[,1]
      START = regionData[,2]
      END = regionData[,3]
      inSNPlist=""

      chromList = unique(CHROM1)
      for(chr in chromList){

	outfile=paste0(outputPrefix, "_", chr)
	if(file.exists(outfile)){
		file.remove(outfile)
	}
      }	
      setVCFobjInCPP(vcfFile, vcfFileIndex, vcfField, t_SampleInModel = sampleInModel)
      for(i in 1:nrow(regionData)){	
        in_chrom=as.character(CHROM1[i])
        in_beg_pd=START[i]
        in_end_pd=END[i]
	regionid = paste0(in_chrom, "_", in_beg_pd, "_", in_end_pd)
	cat("regionid ", regionid, "\n")
	isread = TRUE
        set_iterator_inVcf(inSNPlist, in_chrom, in_beg_pd, in_end_pd)
	REF="x";
	ALT="x"
	ID="x"
	POS=0
	CHROM="1"
	IDList=NULL
	while(isread){
		isreadList = getOneMarkerID_VCF(REF, ALT, ID, POS, CHROM)
		isread = isreadList$isReadVariant
		CHROM = isreadList$CHROM
		POS = isreadList$POS
		REF = isreadList$REF
		ALT = isreadList$ALT
		if(isread){
			IDList = c(IDList, paste(c(CHROM, POS, REF, ALT), collapse=":"))
		}
	}
	if(!is.null(IDList)){
		IDstring = c(regionid, "var", IDList)
		annostring = c(regionid, "anno", rep("null", length(IDList)))
        	groupData = rbind(IDstring, annostring)
        	write.table(groupData, paste0(outputPrefix, "_", in_chrom), sep=" ", quote=F, col.names=F, row.names=F, append=T)
	}else{
		warning("No markers are in the region ", regionid, "\n")
	}
      }


  }

	
}


extractRegions = function(regionFile){

	regionData = data.table::fread(regionFile, header=F, data.table=F)
	if(ncol(regionData) != 3){
		stop("regionFile does not contain three columns\n")
	}
	return(regionData)
}

makeGroupFilewithRegionData = function(markerInfo, regionData, outputPrefix){
     
      chromList = unique(regionData[,1])
      for(chr in chromList){

        outfile=paste0(outputPrefix, "_", chr)
        if(file.exists(outfile)){
                file.remove(outfile)
        }
      }

	for(i in 1:nrow(regionData)){
		 chrom = regionData[i,1]
		 start = regionData[i,2]
		 end = regionData[i,3]
		 regionid = paste0(chrom, "_", start, "_", end)
		         cat("regionid ", regionid, "\n")
		 IDList = markerInfo$ID2[which(markerInfo$CHROM == chrom & markerInfo$POS <= end & markerInfo$POS >= start)]
		 if(!is.null(IDList)){
		 IDstring = c(regionid, "var", IDList)
		 annostring = c(regionid, "anno", rep("null", length(IDList)))
		 groupData = rbind(IDstring, annostring)
		 write.table(groupData, paste0(outputPrefix, "_", chrom), sep=" ", quote=F, col.names=F, row.names=F, append=T)
		 }else{
			warning("No markers are in the region ", regionid, "\n")
		 }
	}
}
