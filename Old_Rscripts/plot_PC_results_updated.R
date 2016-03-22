library(ggplot2)

results <- read.table("BRCA_Tumor_linearModel_400PCby50s_03102016.txt", header = TRUE, stringsAsFactors = FALSE)

start <- 50
num_pc <- 400
step <- 50
step_list <- seq(start, num_pc, by = step)
 
pv0 <- results$pvalue

for (i in 1:(length(step_list))){
#print(paste0("pv",step_list[i]))
assign(paste0("pv",step_list[i]), results[[paste0("pvalue.",i)]])
}

sum_list <- data.frame()

pval_list <- seq(0, num_pc, step)
for (i in pval_list){
sum_list <- rbind(sum_list, sum(get(paste0("pv",i)) < 0.05))
}

pval_list <- cbind(pval_list, sum_list)
#plot(pval_list)
colnames(pval_list) <- c("Num_PCs", "Num_Sig_Genes")
ggplot(data = pval_list, aes( x = Num_PCs, y = Num_Sig_Genes)) + geom_point() + ggtitle("BRCA Tumor")

### Adjust P-vals ###


pva0 <- results$adj_pvalue

for (i in 1:(length(step_list))){
  #print(paste0("pv",step_list[i]))
  assign(paste0("pva",step_list[i]), results[[paste0("adj_pvalue.",i)]])
}

sum_list <- data.frame()

pvala_list <- seq(0, num_pc, step)
for (i in pvala_list){
  sum_list <- rbind(sum_list, sum(get(paste0("pva",i)) < 0.05))
}

pvala_list <- cbind(pvala_list, sum_list)
#plot(pval_list)
colnames(pvala_list) <- c("Num_PCs", "Num_Sig_Genes")
ggplot(data = pvala_list, aes( x = Num_PCs, y = Num_Sig_Genes)) + geom_point() + ggtitle("BRCA Tumor, Adj P-value")
