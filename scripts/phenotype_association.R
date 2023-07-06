Args <- commandArgs(TRUE)
predicted_expr <- Args[1]
phenotype <- Args[2]  # 0/1

predicted_expr <- read.csv(predicted_expr, row.names = 1, check.names = F)
phenotype <- read.csv(phenotype, row.names = 1, check.names = F)

predicted_expr$sample_id <- rownames(predicted_expr)
phenotype$sample_id <- rownames(phenotype)

df <- merge(predicted_expr, phenotype, by = 'sample_id')
colnames(df) <- c('sample_id', 'expr', 'phe')

mod <- summary(glm(phe ~ expr, data = df, family = binomial))
print(mod$coefficients)
