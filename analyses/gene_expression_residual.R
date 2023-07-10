Args <- commandArgs(TRUE)
gene_id <- Args[1]

df_tmm <- read.csv(paste0('tmm/', gene_id, '.csv'))
df_covar <- read.csv('covariate.csv')
df_data <- merge(df_tmm, df_covar, by = 'sample_id')

md <- lm(gene_expr ~ age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + 
           InferredCov1 + InferredCov2 + InferredCov3 + InferredCov4 + InferredCov5 +
           InferredCov6 + InferredCov7 + InferredCov8 + InferredCov9 + InferredCov10 +
           InferredCov11 + InferredCov12 + InferredCov13 + InferredCov14 + InferredCov15, data = df_data)
md_resid <- residuals(md)

df_resid <- data.frame(gene_resid = md_resid, row.names = df_data$sample_id)
write.csv(df_resid, file = paste0('residual/', gene_id, '.csv'), row.names = F)
