library(dplyr)
library(ggplot2)

normal_pc <- read.table("BRCA_Normal_linearModel_30PC_03092016.txt", header = TRUE, stringsAsFactors = FALSE)
pc8_norm <- data.frame()
pc8_norm <- cbind(data.frame(normal_pc$gene, normal_pc$pvalue.8, normal_pc$adj_pvalue.8, normal_pc$Rsquared.8, normal_pc$beta.8, stringsAsFactors = FALSE))

norm8_sigP <- filter(pc8_norm, normal_pc.adj_pvalue.8 < 0.05)


tumor_pc <- read.table("BRCA_Tumor_linearModel_200PCby10s_03102016.txt", header = TRUE, stringsAsFactors = FALSE)
pc150_tum <- data.frame()
pc150_tum <- cbind(data.frame(tumor_pc$gene, tumor_pc$pvalue.15, tumor_pc$adj_pvalue.15, tumor_pc$Rsquared.15, tumor_pc$beta.15, stringsAsFactors = FALSE))

tum150_sigP <- filter(pc150_tum, tumor_pc.adj_pvalue.15 < 0.05)

# Make individual data frames
norm_beta <- data.frame(group = "normal", value = norm8_sigP$normal_pc.beta.8)
tum_beta <- data.frame(group = "tumor", value =  tum150_sigP$tumor_pc.beta.15)
# Combine into one long data frame
plot.data <- rbind(norm_beta, tum_beta)

bp <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot()
bp + ylim(-10,200)

wilcox.test(norm_beta$value, tum_beta$value, paired=FALSE) 

common_genes <- inner_join(tum150_sigP$tumor_pc.gene, norm8_sigP$normal_pc.gene)
