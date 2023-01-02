library(GMMAT)
pheno.file <- system.file("extdata", "pheno.txt", package = "GMMAT")
pheno <- read.table(pheno.file, header = TRUE)
GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))
M10 <- matrix(0, 400, 400)
for(i in 1:40) M10[(i-1)*10+(1:10), (i-1)*10+(1:10)] <- 1
rownames(M10) <- colnames(M10) <- 1:400
Mats <- list(GRM, M10)
#Mats <- list(GRM)
#model3 <- glmmkin(fixed = disease ~ age + sex, data = pheno, id = "id",kins = Mats, family = binomial(link = "logit"), verbose=T)
#model3 <- glmmkin(fixed = disease ~ age + sex, data = pheno, id = "id",kins = GRM, family = binomial(link = "logit"), verbose=T)
#model3 <- glmmkin(fixed = disease ~ age + sex, data = pheno, id = "id",kins = GRM, family = binomial(link = "logit"), verbose=T)


print(M10[1:20,1:20])

pheno2.file <- system.file("extdata", "pheno2.txt", package = "GMMAT")
pheno2 <- read.table(pheno2.file, header = TRUE)
pheno3 = pheno2[which(pheno2$time == 1),]

setseed(1)
wts = abs(rnorm(length(pheno2$y.repeated)))
pheno2$wts = wts
write.table(pheno2, "./GMMAT_pheno_quantitative_repeat.txt_withwts", quote=F, col.names=T, row.names=F`)
model4 <- glmmkin(y.repeated ~ sex, data = pheno2, kins = GRM, id = "id", family = gaussian(link = "identity"), weights=wts)

geno.file <- strsplit(system.file("extdata", "geno.bed", package = "GMMAT"), ".bed", fixed = TRUE)[[1]]
glmm.score(model4, infile = geno.file, outfile = "glmm.score.bed.testoutfile.txt")



library(Matrix)
GRM = readMM("/net/csgspare3/snowwhite.archive/zczhao/SAIGE-GENE-UPDATE/eqtl/PLINK_SLC_EXCLUDE_HWE_MAF0.05_MIND0.1_GENO0.05_VCF_rename_random1000_relatednessCutoff_0.05_1000_randomMarkersUsed.sparseGRM.mtx") 
gnames = read.table("/net/csgspare3/snowwhite.archive/zczhao/SAIGE-GENE-UPDATE/eqtl/PLINK_SLC_EXCLUDE_HWE_MAF0.05_MIND0.1_GENO0.05_VCF_rename_random1000_relatednessCutoff_0.05_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt", header=F)
rownames(GRM) <- colnames(GRM) <- gnames[,1]
pheno_eQTL = read.table("/net/dumbo/home/zhowei/projects/eqtl/pheno_random_10_genes_rand2000.txt", header=T)
library(GMMAT)
model4_eQTL <- glmmkin(ENSG00000117226 ~ Cohort_chem+RNA_snn_res.0.5+Age+Sex+PC1+PC2+PC3+PC4+PC5, data = pheno_eQTL, kins = GRM, id = "geno_ID", family = poisson(link = "log"))




#model4 <- glmmkin(y.repeated ~ sex, data = pheno2, kins = GRM, id = "id", family = gaussian(link = "identity"), random.slope = "time", verbose=T)
#model4 <- glmmkin(y.repeated ~ sex, data = pheno3, kins = Mats, id = "id", family = gaussian(link = "identity"), verbose=T)
#model4 <- glmmkin(y.repeated ~ sex, data = pheno3, kins = GRM, id = "id", family = gaussian(link = "identity"), verbose=T)



#model5 <- glmmkin(y.trend ~ sex, data = pheno2, kins = GRM, id = "id", random.slope = "time", family = gaussian(link = "identity"))
#model5 <- glmmkin(y.trend ~ sex, data = pheno2, kins = GRM, id = "id", family = gaussian(link = "identity"))
