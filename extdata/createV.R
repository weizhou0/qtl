library(Matrix)

#a=matrix(1, nrow=20, ncol=20)
a = diag(rep(1,20))
b=list()
for(i in 1:100){
	b[[i]] = a
}

d = bdiag(b)
writeMM(d, "bdiagMM.mtx")
d=paste0("1a", c(1:2000))
write.table(d, "bdiagMM_sampleID.txt", col.names=F, row.names=F, quote=F)



M10 <- matrix(0, 400, 400)
for(i in 1:40) M10[(i-1)*10+(1:10), (i-1)*10+(1:10)] <- 1
M10 = as(M10, "sparseMatrix")
writeMM(M10, "GMMAT_M10.mtx")
d=c(1:400)
write.table(d, "GMMAT_sampleID.txt", col.names=F, row.names=F, quote=F)


library(GMMAT)
pheno.file <- system.file("extdata", "pheno.txt", package = "GMMAT")

pheno <- read.table(pheno.file, header = TRUE)

write.table(pheno, "GMMAT_pheno.txt", col.names=T, row.names=F, quote=F)


GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
GRM = as(GRM, "sparseMatrix")
writeMM(GRM, "GMMAT_GRM.mtx")

IMatrix = diag(1, nrow=400, ncol=400)
IMatrix = as(IMatrix, "sparseMatrix")
writeMM(IMatrix, "I.mtx")


pheno2.file <- system.file("extdata", "pheno2.txt", package = "GMMAT")
pheno2 <- read.table(pheno2.file, header = TRUE)
pheno3 = pheno2[which(pheno2$time == 1),]
write.table(pheno3, "GMMAT_pheno_quantitative.txt", col.names=T, row.names=F, quote=F)
 write.table(pheno2, "GMMAT_pheno_quantitative_repeat.txt" , col.names=T, row.names=F, quote=F)


