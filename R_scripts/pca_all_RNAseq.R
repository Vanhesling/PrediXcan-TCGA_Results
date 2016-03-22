library(dplyr)
library(ggplot2)
library(ggfortify)
source('R_scripts/multiplot.R')
source('R_scripts/read_Data.R')

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


tumor_expression <- tumor_samples[,-(1:2)]
normal_expression <- normal_samples[,-(1:2)]    

## Remove columns where std dev == 0
tumor_expression <- tumor_expression[, sapply(tumor_expression,function(x) { sd(x) != 0} )]
normal_expression <- normal_expression[, sapply(normal_expression, function(x) { sd(x) != 0} )]

t_tumor_expression <- t(tumor_expression)
t_tumor_expression <- data.frame(t_tumor_expression, stringsAsFactors = FALSE)
t_tumor_expression <- cbind(rownames(t_tumor_expression), t_tumor_expression)
colnames(t_tumor_expression)[1] <- "gene"
t_tumor_expression$gene <- as.character(t_tumor_expression$gene)

t_normal_expression <- t(normal_expression)
t_normal_expression <- data.frame(t_normal_expression, stringsAsFactors = FALSE)
t_normal_expression <- cbind(rownames(t_normal_expression), t_normal_expression)
colnames(t_normal_expression)[1] <- "gene"
t_normal_expression$gene <- as.character(t_normal_expression$gene)

matched_tumor <- semi_join(t_tumor_expression, t_normal_expression, by = "gene") 
matched_normal <- semi_join(t_normal_expression, t_tumor_expression, by = "gene")

matched_tumor <- arrange(matched_tumor, gene)
matched_normal <- arrange(matched_normal, gene)

gene_list <- matched_normal$gene

matched_tumor <- as.data.frame(t(matched_tumor[,(-1)]))
matched_normal <- as.data.frame(t(matched_normal[,(-1)]))

matched_tumor$sample <- "tumor" ## add sample category
matched_normal$sample <- "normal"

all_expression <- rbind(matched_normal, matched_tumor) 
#colnames(all_expression) <- gene_list ## replace colnames with genes
colnames(all_expression)[length(colnames(all_expression))] <- "sample"

all_expression.data <- all_expression[,1:(NCOL(all_expression)-1)]
all_expression.sample <- all_expression[,(NCOL(all_expression))]

#browser()

all_TCGA_pca <- prcomp(all_expression.data, center = TRUE, scale. = TRUE)
#all_TCGA_pca <- prcomp(all_expression[,1:length(gene_list)], center = TRUE, scale. = TRUE)


scores <- data.frame(all_expression.sample, all_TCGA_pca$x[,1:4])
colnames(scores)[1] <- "sample"

pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=factor(sample))+ theme(legend.position="none")
pc1.3 <- qplot(x=PC1, y=PC3, data=scores, colour=factor(sample))+ theme(legend.position="none")
pc1.4 <- qplot(x=PC1, y=PC4, data=scores, colour=factor(sample))+ theme(legend.position="none")
pc2.3 <- qplot(x=PC2, y=PC3, data=scores, colour=factor(sample))+ theme(legend.position="none")
pc2.4 <- qplot(x=PC2, y=PC4, data=scores, colour=factor(sample))+ theme(legend.position="none")
pc3.4 <- qplot(x=PC3, y=PC4, data=scores, colour=factor(sample))+ theme(legend.position="none")

multiplot(pc1.2, pc1.3, pc1.4, pc2.3, pc2.4, pc3.4, cols = 2)
