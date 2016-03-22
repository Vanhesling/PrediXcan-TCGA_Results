library(dplyr)
library(beepr)

prdx_tumor <- read.table("BRCA_tumor_PrediXcan-expression.txt", header = TRUE, stringsAsFactors = FALSE)
tumor_samples <- read.table("BRCA_tumor_TCGA-RNAseq-expression.txt", header = TRUE, stringsAsFactors = FALSE)

prdx_normal <- read.table("BRCA_normal_PrediXcan-expression.txt", header = TRUE, stringsAsFactors = FALSE)
normal_samples <- read.table("BRCA_normal_TCGA-RNAseq-expression.txt", header = TRUE, stringsAsFactors = FALSE)

prdx_normal <- arrange(prdx_normal, barcode)
normal_samples <- arrange(normal_samples, barcode)

prdx_tumor <- arrange(prdx_tumor, barcode)
tumor_samples <- arrange(tumor_samples, barcode)


##comment out for tumor or normal
prdx_expression <- prdx_tumor[,-(1:2)]
#prdx_expression <- prdx_normal[,-(1:2)]
tumor_expression <- tumor_samples[,-(1:2)]
#tumor_expression <- normal_samples[,-(1:2)]

## Remove columns where std dev == 0
prdx_expression <- prdx_expression[, sapply(prdx_expression,function(x) { sd(x) != 0} )]
tumor_expression <- tumor_expression[, sapply(tumor_expression, function(x) { sd(x) != 0} )]

genes_in_prdx <- colnames(prdx_expression)
genes_in_RNAseq <- colnames(tumor_expression)

TCGA_pca <- prcomp(tumor_expression, center = TRUE, scale. = TRUE)

lm.beta <- function (MOD) {
  b <- summary(MOD)$coef[2, 1]
  sx <- sd(as.matrix(MOD$model[2]))
  sy <- sd(as.matrix(MOD$model[1]))
  beta <- b * sx/sy
  return(beta)
}


lm_fits <- data.frame()

## First time no PCs
for (gene in genes_in_prdx){
  if (gene %in% genes_in_RNAseq){
    fit <- lm(tumor_expression[[gene]] ~ prdx_expression[[gene]])
    beta <- lm.beta(fit)
    pvalue <- summary(fit)$coefficients[2,4] 
    Rsquared <- summary(fit)$r.squared
    lm_fits <- rbind(lm_fits, data.frame(gene, beta, Rsquared, pvalue))
  }
}

adj_pvalue <- p.adjust(lm_fits$pvalue, method = "BH")
lm_fits <- cbind(lm_fits, adj_pvalue)


all_pc_list <- list()
for (i in c(1:420)){
  assign(paste0("PC",i), TCGA_pca$x[,i])
  all_pc_list <- c(all_pc_list, paste0("PC",i))
}

for (ii in seq(50, 400, by = 50)){
  print(paste(ii,":",Sys.time()))
  pc_list <- all_pc_list[1:ii]
  new_fits <- data.frame()
  for (gene in genes_in_prdx){
    if (gene %in% genes_in_RNAseq){
      pc_parts <- paste(pc_list, collapse=" + ")
      fit <- lm(as.formula(paste0("tumor_expression[[gene]] ~ prdx_expression[[gene]] + ",pc_parts)))
      beta <- lm.beta(fit)[[1]][1]
      pvalue <- summary(fit)$coefficients[2,4] 					
      Rsquared <- summary(fit)$r.squared
      new_fits <- rbind(new_fits, data.frame(beta, Rsquared, pvalue))
    }
  }
  adj_pvalue <- p.adjust(new_fits$pvalue, method = "BH")
  new_fits <- cbind(new_fits, adj_pvalue)
  lm_fits <- cbind(lm_fits, new_fits)	
}


write.table(lm_fits, file = "BRCA_Tumor_linearModel_150PC_03142016.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
beep(3)

