read_all_data <- function(){
  #browser()
  ## Read in predicted and actual gene expression data
  prdx_tumor <<- read.table("expression_data/BRCA_tumor_PrediXcan-expression.txt", header = TRUE, stringsAsFactors = FALSE)
  tumor_samples <<- read.table("expression_data/BRCA_tumor_TCGA-RNAseq-expression.txt", header = TRUE, stringsAsFactors = FALSE)
  print("Tumor data read")
  prdx_normal <<- read.table("expression_data/BRCA_normal_PrediXcan-expression.txt", header = TRUE, stringsAsFactors = FALSE)
  normal_samples <<- read.table("expression_data/BRCA_normal_TCGA-RNAseq-expression.txt", header = TRUE, stringsAsFactors = FALSE)
  print("Normal data read")
}