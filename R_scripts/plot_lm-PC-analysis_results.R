library(ggplot2)

plot_PC_results <- function(my_file, start_pc, stop_pc, step_pc){
  step_list <- seq(start_pc, stop_pc, by = step_pc)  
  results <- read.table(file = my_file, header = TRUE, stringsAsFactors = FALSE)
  plot_title <- readline(prompt = 'Enter title for plot: ')
 
  pv0 <- results$pvalue

  for (i in 1:(length(step_list))){
    assign(paste0("pv",step_list[i]), results[[paste0("pvalue.",i)]])
  }

  sum_list <- data.frame()

  pval_list <- seq(0, stop_pc, step_pc)
  for (i in pval_list){
    sum_list <- rbind(sum_list, sum(get(paste0("pv",i)) < 0.05))
  }

  pval_list <- cbind(pval_list, sum_list)
  colnames(pval_list) <- c("Num_PCs", "Num_Sig_Genes")
  pval_plot <- ggplot(data = pval_list, aes( x = Num_PCs, y = Num_Sig_Genes)) + geom_point() + ggtitle(paste0(plot_title,", p-value"))
  print(pval_plot)
  
### Adjust P-vals ###

  pva0 <- results$adj_pvalue

  for (i in 1:(length(step_list))){
    assign(paste0("pva",step_list[i]), results[[paste0("adj_pvalue.",i)]])
  }

  sum_list <- data.frame()

  pvala_list <- seq(0, stop_pc, step_pc)
  for (i in pvala_list){
    sum_list <- rbind(sum_list, sum(get(paste0("pva",i)) < 0.05))
  }
  
  pvala_list <- cbind(pvala_list, sum_list)
  colnames(pvala_list) <- c("Num_PCs", "Num_Sig_Genes")
  adj_pval_plot <- ggplot(data = pvala_list, aes( x = Num_PCs, y = Num_Sig_Genes)) + geom_point() + ggtitle(paste0(plot_title,", Adj p-value"))
  print(adj_pval_plot)
}