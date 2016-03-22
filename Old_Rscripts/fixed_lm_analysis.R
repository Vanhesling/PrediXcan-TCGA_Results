library(dplyr)
library(beepr)
source('read_Data.R')

lm_analysis <- function(sample_type = readline(prompt = "Enter 'tumor' or 'normal'"))

## Check if expression data exists in global environment, if not, read it in
  if (exists("prdx_tumor") == FALSE){
    if (exists("prdx_normal") == FALSE){
      if(exists("normal_samples") == FALSE){
        if (exists("tumor_samples") == FALSE){
          read_all_data()
        }
      }
    }
  }

## sort by barcode so predicted expression and actual expression are in same order
prdx_normal <- arrange(prdx_normal, barcode)
normal_samples <- arrange(normal_samples, barcode)

prdx_tumor <- arrange(prdx_tumor, barcode)
tumor_samples <- arrange(tumor_samples, barcode)

##comment out for tumor or normal
#prdx_expression <- prdx_tumor[,-(1:2)]
prdx_expression <- prdx_normal[,-(1:2)]
#tumor_expression <- tumor_samples[,-(1:2)]
tumor_expression <- normal_samples[,-(1:2)]

## Remove columns where std dev == 0
prdx_expression <- prdx_expression[, sapply(prdx_expression,function(x) { sd(x) != 0} )]
tumor_expression <- tumor_expression[, sapply(tumor_expression, function(x) { sd(x) != 0} )]

genes_in_prdx <- colnames(prdx_expression)
genes_in_RNAseq <- colnames(tumor_expression)

TCGA_pca <- prcomp(tumor_expression, center = TRUE, scale. = TRUE)

all_pc_list <- list()
for (i in c(1:8)){
  assign(paste0("PC",i), TCGA_pca$x[,i])
  all_pc_list <- c(all_pc_list, paste0("PC",i))
}

lm_fits <- data.frame()

for (gene in genes_in_prdx){
  if (gene %in% genes_in_RNAseq){
    fit <- lm(as.formula(paste0("tumor_expression[[gene]] ~ prdx_expression[[gene]] + ",all_pc_list)))
    estimate <- coef(summary(fit))[2,1]
    pvalue <- summary(fit)$coefficients[2,4] 
    Rsquared <- summary(fit)$r.squared
    lm_fits <- rbind(lm_fits, data.frame(gene, estimate, Rsquared, pvalue))
  }
}

adj_pvalue <- p.adjust(lm_fits$pvalue, method = "BH")
lm_fits <- cbind(lm_fits, adj_pvalue)

write.table(lm_fits, file = "BRCA_Normal_linearModel_8C_03142016.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
beep(3)

