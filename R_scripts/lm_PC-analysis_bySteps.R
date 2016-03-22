library(dplyr)
library(beepr)
source('R_scripts/plot_lm-PC-analysis_results.R')
source('R_scripts/read_Data.R')

# sample_type = 'normal' or 'tumor'
# start_pc, stop_pc, step_pc = start, end and step in sequence function 

# this function will generate linear model fits for the predicted gene expression based 
# on actual tissue expression, and add increasing numbers of principal components (determined
# by the stop_ and step_pc defined) to see how that improves model estimate and p-value

# onlyMatched -- by default, set to FALSE. If set to true, include only the tumor samples that 
# have a matched, normal sample

lm_analysis_byStep <- function(sample_type, stop_pc, step_pc, onlyMatched = FALSE){
  start_pc <- step_pc
  
  # Check to see if expression data frames exist in global environment, if not, read data files
  if (exists("prdx_tumor") == FALSE){
    if (exists("prdx_normal") == FALSE){
      if(exists("normal_samples") == FALSE){
        if (exists("tumor_samples") == FALSE){
          read_all_data()
        }
      }
    }
  }
  
  
  ## Uncomment these two lines to get ONLY tumor samples with matched normal
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
  
  ## Generate linear models without including any PCs
  lm_fits <- data.frame()
  
  for (gene in genes_in_prdx){
    if (gene %in% genes_in_RNAseq){
      fit <- lm(tissue_expression[[gene]] ~ prdx_expression[[gene]])
      estimate <- coef(summary(fit))[2,1]
      pvalue <- summary(fit)$coefficients[2,4] 
      Rsquared <- summary(fit)$r.squared
      lm_fits <- rbind(lm_fits, data.frame(gene, estimate, Rsquared, pvalue))
    }
  }
  
  #Adjust p-value
  adj_pvalue <- p.adjust(lm_fits$pvalue, method = "BH")
  lm_fits <- cbind(lm_fits, adj_pvalue)
  
  ## create list with all of the PCs
  all_pc_list <- list()
  for (i in c(1:stop_pc)){
    assign(paste0("PC",i), TCGA_pca$x[,i])
    all_pc_list <- c(all_pc_list, paste0("PC",i))
  }
  
## index PC list by steps, fit linear model based on new PCs
  for (ii in seq(start_pc, stop_pc, by = step_pc)){
    print(paste(ii,":",Sys.time()))
    pc_list <- all_pc_list[1:ii]
    new_fits <- data.frame()
    for (gene in genes_in_prdx){
      if (gene %in% genes_in_RNAseq){
        pc_parts <- paste(pc_list, collapse=" + ")
        fit <- lm(as.formula(paste0("tissue_expression[[gene]] ~ prdx_expression[[gene]] + ",pc_parts)))
        estimate <- coef(summary(fit))[2,1]
        pvalue <- summary(fit)$coefficients[2,4] 					
        Rsquared <- summary(fit)$r.squared
        new_fits <- rbind(new_fits, data.frame(estimate, Rsquared, pvalue))
      }
    }
    adj_pvalue <- p.adjust(new_fits$pvalue, method = "BH")
    new_fits <- cbind(new_fits, adj_pvalue)
    lm_fits <- cbind(lm_fits, new_fits)	
  }
  
  beep(3) ## Include if you want to be notified when user input needed
  filename <- readline(prompt = "Enter file name: ")
  write.table(lm_fits, file = filename, col.names = TRUE, row.names = FALSE, quote = FALSE)
  plot_PC_results(filename, start_pc, stop_pc, step_pc)
  
}
