library(dplyr)
library(beepr)
source('R_scripts/read_Data.R')

#sample_type = 'normal' or 'tumor'
#pc_len = number of PCs to include in analysis, starting with 0
# onlyMatched -- by default, set to FALSE. If set to true, include only the tumor samples that 
# have a matched, normal sample

lm_analysis_toFile <- function(sample_type, pc_len, onlyMatched = FALSE){
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
  
  ## Use to get ONLY tumor samples with matched normal
  if (onlyMatched == TRUE){
    prdx_tumor <- semi_join(prdx_tumor, prdx_normal, by = "barcode")
    tumor_samples <- semi_join(tumor_samples, normal_samples, by = "barcode")    
  }
  
 
  ##set expression based on tumor or normal, remove first two columns with barcode etc.
  if (sample_type == 'normal'){
    prdx_expression <- prdx_normal[,-(1:2)]
    tissue_expression <- normal_samples[,-(1:2)]    
  }
  if (sample_type == 'tumor'){
    prdx_expression <- prdx_tumor[,-(1:2)]
    tissue_expression <- tumor_samples[,-(1:2)]
  }

  ## Remove columns where std dev == 0
  prdx_expression <- prdx_expression[, sapply(prdx_expression,function(x) { sd(x) != 0} )]
  tissue_expression <- tissue_expression[, sapply(tissue_expression, function(x) { sd(x) != 0} )]

# Genes in predicted and actual expression data
genes_in_prdx <- colnames(prdx_expression)
genes_in_RNAseq <- colnames(tissue_expression)

## Calculate the principal compoments
TCGA_pca <- prcomp(tissue_expression, center = TRUE, scale. = TRUE)


## If no PCs to be included, follow this loop
if (pc_len == 0){
  lm_fits <- data.frame()
  
  ## Go through and fit lm for each gene in common between predicted and RNAseq data
  for (gene in genes_in_prdx){
    if (gene %in% genes_in_RNAseq){
      fit <- lm(tissue_expression[[gene]] ~ prdx_expression[[gene]])
      estimate <- coef(summary(fit))[2,1]
      pvalue <- summary(fit)$coefficients[2,4] 
      Rsquared <- summary(fit)$r.squared
      lm_fits <- rbind(lm_fits, data.frame(gene, estimate, Rsquared, pvalue))
    }
  }
  
  adj_pvalue <- p.adjust(lm_fits$pvalue, method = "BH")
  lm_fits <- cbind(lm_fits, adj_pvalue)
  
  beep(3) ## Include if you want to be notified when done running
  filename <- readline(prompt = "Enter file name: ")
  write.table(lm_fits, file = filename, col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(paste0("created ",filename,". No PCs included in linear model."))
}

all_pc_list <- list()
for (i in c(1:pc_len)){
  assign(paste0("PC",i), TCGA_pca$x[,i])
  all_pc_list <- c(all_pc_list, paste0("PC",i))
}

lm_fits <- data.frame()

for (gene in genes_in_prdx){
  if (gene %in% genes_in_RNAseq){
    fit <- lm(as.formula(paste0("tissue_expression[[gene]] ~ prdx_expression[[gene]] + ",all_pc_list)))
    estimate <- coef(summary(fit))[2,1]
    pvalue <- summary(fit)$coefficients[2,4] 
    Rsquared <- summary(fit)$r.squared
    lm_fits <- rbind(lm_fits, data.frame(gene, estimate, Rsquared, pvalue))
  }
}

#Adjust p-value
adj_pvalue <- p.adjust(lm_fits$pvalue, method = "BH")
lm_fits <- cbind(lm_fits, adj_pvalue)

beep(3) ## Include if you want to be notified when done running
filename <- readline(prompt = "Enter file name: ")
write.table(lm_fits, file = filename, col.names = TRUE, row.names = FALSE, quote = FALSE)

}
