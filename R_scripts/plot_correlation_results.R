library(dplyr)
library(ggplot2)
library(grid)

plot_correlations <- function(normal_file_name, tumor_file_name, test_method){
      normal_is <- TRUE
      for (corr_file in c(normal_file_name, tumor_file_name)){
        correlations <- read.table(corr_file, header = TRUE)
        if (test_method != "Pearson" & test_method != "Spearman"){
          return("Not a correct method, enter either Pearson or Spearman") #exit functiod
        }
        if (test_method == "Pearson"){
          my_data <- data.frame(cbind(correlations[,1:3], correlations[,6]))
          }  
        if (test_method == "Spearman"){
          my_data <- data.frame(cbind(correlations[,1], correlations[,4:5], correlations[,7]))
          }
        colnames(my_data) <- c("gene","r","pvalue","adj_pvalue")  
        sig_data <- filter(my_data, adj_pvalue < 0.05)
        if (normal_is == TRUE){
          n_data <- data.frame(sample = "normal", r = sig_data$r)
          plot.data <- rbind(n_data)
          normal_mean <- round(mean(sig_data$r), 3)
          normal_is <- FALSE
        }
      t_data <- data.frame(sample = "tumor", r = sig_data$r)
      plot.data <- rbind(plot.data, t_data)
      tumor_mean <- round(mean(sig_data$r), 3)
      }
  wt <- wilcox.test(n_data$r, t_data$r, paired = FALSE)
  wt_pvalue <- wt$p.value
  p <- grobTree(textGrob(paste("P-value:",wt_pvalue), x = .5, hjust = .5, y = 1, vjust = 1), gp = gpar(fontsize = 12))
  t_mean <- grobTree(textGrob(paste("tumor, mean r:",tumor_mean), x = 1, hjust = 1, y = .1, vjust = 0), gp = gpar(fontsize = 12))
  n_mean <- grobTree(textGrob(paste("normal, mean r:",normal_mean), x = .45, hjust = .5, y = .1, vjust = 0), gp = gpar(fontsize = 12))
  boxplot <- ggplot(plot.data, aes(x = sample, y = r, fill = sample)) + geom_boxplot()
  boxplot <- boxplot + annotation_custom(p) + theme_classic() + ggtitle(paste("Comparison of r, method = ",test_method))
  boxplot + annotation_custom(t_mean) + annotation_custom(n_mean)
}



