library(glmnet)

Args <- commandArgs(TRUE)
microbe_file <- Args[1]
genotype_file <- Args[2]
model_coef <- Args[3]
model_rds <- Args[4]
gene_id <- Args[5]

# input files 
taxa <- read.csv(microbe_file, check.names = F, row.names = 1)
dosage <- read.csv(genotype_file, check.names = F, row.names = 1)

# predictor tab
keep_samples <- intersect(rownames(taxa), rownames(dosage))
taxa <- taxa[keep_samples,]
dosage <- dosage[keep_samples,]
x_new <- as.matrix(cbind(taxa, dosage))

# create interaction terms
df_coef <- read.csv(model_coef)
feature_inter <-  df_coef$feature[grep('*', df_coef$feature, fixed = T)]
feature_taxa <- sapply(strsplit(feature_inter, '*', fixed = T), function(x) x[1])
feature_snp <- sapply(strsplit(feature_inter, '*', fixed = T), function(x) x[2])
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
x_new <- cbind(mtx_interact, x_new)

# load model
fit <- readRDS(model_rds)
predictors <- rownames(coef(fit))
x_new <- x_new[, predictors[-1]]

# predict expression
pred_min <- predict(fit, newx = x_new, s="lambda.min")
pred_min <- as.data.frame(pred_min)
colnames(pred_min) <- c(gene_id)
write.csv(pred_min, file = paste0('predicted_expr_', gene_id, '.csv'))
