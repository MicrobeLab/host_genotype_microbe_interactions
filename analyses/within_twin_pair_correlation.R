library(edgeR)

get_fc <- function(gene_line){
  # input is one line of gene-cpm table (one gene across samples)
  num_samples <- length(gene_line)
  ind1_vec <- gene_line[seq(1, num_samples, 2)]
  ind2_vec <- gene_line[seq(2, num_samples, 2)]
  fc_vec <- log2((ind1_vec+0.01) / (ind2_vec+0.01))
  return(fc_vec)
}

gene_count <- read.csv('twins_gene_count.csv', row.names = 1)
gene_cpm <- cpm(gene_count)
fc_gene <- t(apply(gene_cpm, 1, get_fc))

fc_gene_abs <- abs(fc_gene)
avg_gene_fc <- rowMeans(fc_gene_abs)
med_gene_fc <- apply(fc_gene_abs, 1, median)

fc_gene <- fc_gene[avg_gene_fc >= 2 & med_gene_fc >= 2,]

microbe_cpm <- read.csv('twins_microbe_cpm.csv', row.names = 1)
fc_microbe <- t(apply(microbe_cpm, 1, get_fc))

stopifnot(colnames(fc_gene) == colnames(fc_microbe))

microbe <- c()
gene <- c()
spearman_pval <- c()
spearman_rho <- c()

for(i in 1:nrow(fc_microbe)){
  for(j in 1:nrow(fc_gene)){
    microbe <- c(microbe, rownames(fc_microbe)[i])
    gene <- c(gene, rownames(fc_gene)[j])
    vec_microbe <- as.numeric(fc_microbe[i,])
    vec_gene <- as.numeric(fc_gene[j,])
    idx_pres <- vec_microbe != 0 & vec_gene != 0
    vec_microbe <- vec_microbe[idx_pres]
    vec_gene <- vec_gene[idx_pres]
    spear_cor <- cor.test(vec_microbe, vec_gene, method = 'spearman')
    spearman_pval <- c(spearman_pval, spear_cor$p.value)
    spearman_rho <- c(spearman_rho, spear_cor$estimate)
  }
  print(i)
}

df_cor <- data.frame(gene, microbe, spearman_rho, spearman_pval)
write.csv(df_cor, file = 'gene_microbe_cor.csv', row.names = F)
