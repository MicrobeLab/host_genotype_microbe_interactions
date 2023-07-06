library(glmnet)

Args <- commandArgs(TRUE)
microbe_file <- Args[1]
genotype_file <- Args[2]
expr_file <- Args[3]
output_prefix <- Args[4]

# input files 
taxa <- read.csv(microbe_file, check.names = F, row.names = 1)
dosage <- read.csv(genotype_file, check.names = F, row.names = 1)
gene_expr <- read.csv(expr_file, check.names = F, row.names = 1)

# predictor tab
keep_samples <- intersect(rownames(taxa), rownames(dosage))
taxa <- taxa[keep_samples,]
dosage <- dosage[keep_samples,]
x_predictor <- as.matrix(cbind(taxa, dosage))

# response tab
y_response <- gene_expr[rownames(x_predictor),]

# train the model without interaction terms (the first pass)
cvfit <- cv.glmnet(x_predictor, y_response[1], alpha = 0.5, nfolds = 10)

# get non-zero model coefficients
df_coef <- as.data.frame(as.matrix(coef(cvfit, s="lambda.min")))
colnames(df_coef) <- 'coefficient_EN'
df_coef <- subset(df_coef, coefficient_EN != 0)
stopifnot(nrow(df_coef) > 0)
df_coef$feature <- rownames(df_coef)
df_coef <- df_coef[, c('feature', 'coefficient_EN')]
df_coef <- subset(df_coef, feature != '(Intercept)')

# check if the model has non-zero coefficients for SNPs and microbes
feature_taxa <- df_coef$feature[grep('microbe', df_coef$feature, fixed = T)]
feature_snp <- df_coef$feature[grep('SNP', df_coef$feature, fixed = T)]
stopifnot(length(feature_taxa) > 0 & length(feature_snp) > 0)

# create interaction terms
pick_taxa <- taxa[, feature_taxa, drop = FALSE]
pick_snp <- dosage[, feature_snp, drop = FALSE]
mtx_interact <- matrix(nrow = nrow(pick_snp), ncol = length(feature_snp)*length(feature_taxa))
rownames(mtx_interact) <- rownames(pick_snp)
feature_interaction <- c()
k <- 0
for(i in 1:length(feature_taxa)){
  for(j in 1:length(feature_snp)){
    feature_interaction <- c(feature_interaction, paste0(feature_taxa[i], '*', feature_snp[j]))
    k <- k + 1
    mtx_interact[,k] <- pick_taxa[,i] * pick_snp[,j]
  }
}
colnames(mtx_interact) <- feature_interaction

# predictor tab with interactions
x_predictor <- cbind(mtx_interact, x_predictor)

# train the model with interaction terms (the second pass)
cvfit <- cv.glmnet(x_predictor, y_response[1], alpha = 0.5, nfolds = 10)

# check the presence of interactions with non-zero coefficients
df_coef <- as.data.frame(as.matrix(coef(cvfit, s="lambda.min")))
colnames(df_coef) <- 'coefficient_EN'
df_coef <- subset(df_coef, coefficient_EN != 0)
stopifnot(nrow(df_coef) > 0)
df_coef$feature <- rownames(df_coef)
df_coef <- df_coef[, c('feature', 'coefficient_EN')]
df_coef <- subset(df_coef, feature != '(Intercept)')
feature_inter <-  df_coef$feature[grep('*', df_coef$feature, fixed = T)]
stopifnot(length(feature_inter) > 0)

# save the model
saveRDS(cvfit, file = paste0(output_prefix, '.rds')) 
write.csv(df_coef, file = paste0(output_prefix, '.csv'))
