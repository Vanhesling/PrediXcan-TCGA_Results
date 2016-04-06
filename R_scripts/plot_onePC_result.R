library(dplyr)
library(ggplot2)
library(grid)

## Analyze the results from the linear model analysis function
## Plots the model_parameter (Rsquared or estimate) for tumor vs. normal

plot_lm_results <- function(normal_file, tumor_file, model_parameter, statsSig = FALSE){

  normal_pc <- read.table(normal_file, header = TRUE, stringsAsFactors = FALSE)
  tumor_pc <- read.table(tumor_file, header = TRUE, stringsAsFactors = FALSE)

  #filter for genes with adjusted p-value less than 0.05
  if (statsSig == TRUE){
    norm_sigP <- filter(normal_pc, adj_pvalue < 0.05)
    tum_sigP <- filter(tumor_pc, adj_pvalue < 0.05)  
  }
  if (statsSig == FALSE){
    norm_sigP <- normal_pc
    tum_sigP <- tumor_pc  
  }
  
  
  #browser()

  # Make individual data frames
  norm_plot <- data.frame(group = "normal", model_results = norm_sigP[[model_parameter]])
  tum_plot <- data.frame(group = "tumor", model_results =  tum_sigP[[model_parameter]])
  
  # Combine into one long data frame
  plot.data <- rbind(norm_plot, tum_plot)
  
  wt <- wilcox.test(norm_plot[["model_results"]], tum_plot[["model_results"]], paired=FALSE) 
  wt_pvalue <- wt$p.value
  if (wt_pvalue > 0.001){
    wt_pvalue <- round(wt_pvalue, 4)
  }
  p <- grobTree(textGrob(paste("P-value:",wt_pvalue), x = .5, hjust = .5 , y = 1, vjust = 1), gp = gpar(fontsize = 12))
  
  bp <- ggplot(plot.data, aes(x=group, y=model_results, fill=group)) + geom_boxplot() + ylab(model_parameter) 

  plot_title <- readline(prompt = "Enter plot title: ") # Enter title for boxplot
  bp <- bp + ggtitle(plot_title) + annotation_custom(p) +theme_classic() 

  print(bp)
  
  change_ylim <- readline(prompt = "Change y axis? Enter 'YES' or 'NO' or 'LOG': ")
  if (change_ylim == "YES"){
    ymin <- readline(prompt = "Enter y-min: ")
    ymax <- readline(prompt = "Enter y-max: ")
    #browser()
    bp <- bp + ylim(as.numeric(ymin), as.numeric(ymax))
  }
  if (change_ylim == "LOG"){
    bp <- bp + scale_y_log10()
  }
  
  return(bp)
}