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
model3 <- glmmkin(fixed = disease ~ age + sex, data = pheno, id = "id",kins = Mats, family = binomial(link = "logit"), verbose=T)



pheno2.file <- system.file("extdata", "pheno2.txt", package = "GMMAT")
pheno2 <- read.table(pheno2.file, header = TRUE)
pheno3 = pheno2[which(pheno2$id == 1),]
write.table(pheno3, "GMMAT_pheno_quantitative.txt", col.names=T, row.names=F, quote=F)
model4 <- glmmkin(y.repeated ~ sex, data = pheno2, kins = GRM, id = "id", family = gaussian(link = "identity"))
