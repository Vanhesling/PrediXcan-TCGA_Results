library(ggplot2)
library(grid)
source('R_scripts/read_Data.R')

plot_expression <- function(sample_type, gene){
  
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
  
  if (sample_type == "normal"){
    actual_expression <- normal_samples
    predicted_expression <- prdx_normal
  }
  if (sample_type == "tumor"){
    actual_expression <- tumor_samples
    predicted_expression <- prdx_tumor
  }
  if (sample_type != "normal" & sample_type != "tumor"){
    return ("Not a valid sample type, enter 'normal' or 'tumor'")
  }
  
  gene_df <- data.frame()
  gene_df <- cbind(actual_expression[[gene]], predicted_expression[[gene]])
  colnames(gene_df) <- c("Observed_Expression", "Predicted_Expression")
  gene_df <- data.frame(gene_df)
  
  gene_fit <- lm(gene_df$Observed_Expression ~ gene_df$Predicted_Expression)
  g_plot <-ggplot(data = gene_df, aes(x = Predicted_Expression, y = Observed_Expression)) + 
            geom_point() + geom_smooth(method = lm)
  r <- grobTree(textGrob(paste("R squared: ", round(summary(gene_fit)$r.squared,3)), x = 0.01, hjust = 0, y = .99, vjust = 1))
  g_plot <- g_plot + ggtitle(paste0(gene,", ",sample_type)) + annotation_custom(r)

  return(g_plot)
}

