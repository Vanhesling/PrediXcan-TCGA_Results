
#ggplot(data = ERAP2, aes(x = Predicted, y = Observed)) + geom_point() + geom_smooth(method = lm) + ggtitle("ERAP2") +geom_text(data = NULL, x = -.75, y = 6000, label = paste("R^2: ",summary(erap2_fit)$r.squared,sep=""), parse = TRUE


gene      beta  Rsquared       pvalue   adj_pvalue    beta.1 Rsquared.1     pvalue.1 adj_pvalue.1
1  XRRA1 0.8689420 0.7550602 6.370777e-12 2.747940e-08  75.21601  0.7724873 3.924498e-12 1.239829e-08
2 POMZP3 0.8680432 0.7534990 7.104268e-12 2.747940e-08  59.72629  0.7642464 1.122393e-11 1.761708e-08
3  C1QL3 0.8647627 0.7478145 1.050436e-11 2.747940e-08  74.47940  0.7552769 1.550748e-11 2.028379e-08
4   UTS2 0.8443111 0.7128612 9.757306e-11 1.914383e-07 384.78749  0.7148822 2.129756e-10 1.857147e-07
5  TDGF1 0.8388601 0.7036863 1.675521e-10 2.629898e-07 123.68889  0.7168510 1.725550e-10 1.692765e-07
6   LDHC 0.8249859 0.6806017 6


plot_correlation <- function(GENE){
GENE_list <- data.frame()
GENE_list <- cbind(normal_samples[[GENE]], prdx_normal[[GENE]])
print(head(GENE_list))
#GENE_list <- data.frame(GENE_list)
print(colnames(GENE_list)) #<- c("Observed_Expression", "Predicted_Expression")
#GENE_fit <- lm(GENE_list$Observed_Expression ~ GENE_list$Predicted_Expression)
#ggplot(data = GENE_list, aes(x = Predicted_Expression, y = Observed_Expression)) + geom_point() + geom_smooth(method = lm) + ggtitle(GENE) + geom_text(data = NULL, x = 0, y = 0, label = paste("R^2: ", summary(GENE_fit)$r.squared), parse = TRUE)
}


XRRA1 <- cbind(normal_samples$XRRA1, prdx_normal$XRRA1)
XRRA1 <- data.frame(XRRA1)
colnames(XRRA1) <- c("Observed_Expression", "Predicted_Expression")
XRRA1_fit <- lm(XRRA1$Observed_Expression ~ XRRA1$Predicted_Expression)
ggplot(data = XRRA1, aes(x = Predicted_Expression, y = Observed_Expression)) + geom_point() + geom_smooth(method = lm) + ggtitle("XRRA1")
+ geom_text(data = NULL, x = -.250, y = 700, label= paste("R^2: ",summary(XRRA1_fit)$r.squared), parse = TRUE)

           
POMZP3 <- cbind(normal_samples$POMZP3, prdx_normal$POMZP3)
POMZP3 <- data.frame(POMZP3)
colnames(POMZP3) <- c("Observed_Expression", "Predicted_Expression")
POMZP3_fit <- lm(POMZP3$Observed_Expression ~ POMZP3$Predicted_Expression)
popmzp3_plot <- ggplot(data = POMZP3, aes(x = Predicted_Expression, y = Observed_Expression)) + geom_point() + geom_smooth(method = lm) + ggtitle("POMZP3")
ggplot(data = POMZP3, aes(x = Predicted_Expression, y = Observed_Expression)) + geom_point() + geom_smooth(method = lm) + ggtitle("POMZP3") + geom_text(data = NULL, x = 0.1, y = 700, label = paste("R^2: ", summary(POMZP3_fit)$r.squared), parse = TRUE)

C1QL3 <- cbind(normal_samples$C1QL3, prdx_normal$C1QL3)
C1QL3 <- data.frame(C1QL3)
colnames(C1QL3) <- c("Observed_Expression", "Predicted_Expression")
C1QL3_fit <- lm(C1QL3$Observed_Expression ~ C1QL3$Predicted_Expression)
#C1QL3_plot <- ggplot(data = POMZP3, aes(x = Predicted_Expression, y = Observed_Expression)) + geom_point() + geom_smooth(method = lm) + ggtitle("POMZP3")
ggplot(data = C1QL3, aes(x = Predicted_Expression, y = Observed_Expression)) + geom_point() + geom_smooth(method = lm) + ggtitle("C1QL3") + geom_text(data = NULL, x = 0.25, y = 60, label = paste("R^2: ", summary(C1QL3_fit)$r.squared), parse = TRUE)

