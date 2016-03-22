library(ggplot2)

start <- 10
num_pc <- 100
step <- 10

pv0 <- results$pvalue

for (i in seq(start, num_pc, by = step)){
assign(paste0("pv",i), results[[paste0("pvalue.",i)]])
}

sum_list <- data.frame()
#pval_list <- cbind(0:num_pc)
pval_list <- seq(0,num_pc, by = )
for (i in 0:num_pc){
sum_list <- rbind(sum_list, sum(get(paste0("pv",i)) < 0.05))
}

pval_list <- cbind(pval_list, sum_list)
#plot(pval_list)
colnames(pval_list) <- c("Num_PCs", "Num_Sig_Genes")
ggplot(data = pval_list, aes( x = Num_PCs, y = Num_Sig_Genes)) + geom_point() + ggtitle("BRCA Tumor")

### Adjust P-vals ###

pva0 <- results$adj_pvalue

for (i in 1:num_pc){
assign(paste0("pva",i), results[[paste0("adj_pvalue.",i)]])
}

sum_list <- data.frame()
pvala_list <- cbind(0:num_pc)
for (i in 0:num_pc){
sum_list <- rbind(sum_list, sum(get(paste0("pva",i)) < 0.05))
}

pvala_list <- cbind(pvala_list, sum_list)
#plot(pvala_list)
colnames(pvala_list) <- c("Num_PCs", "Num_Sig_Genes")
ggplot(data = pvala_list, aes( x = Num_PCs, y = Num_Sig_Genes)) + geom_point() + ggtitle("Pval_adj BRCA Tumor")
