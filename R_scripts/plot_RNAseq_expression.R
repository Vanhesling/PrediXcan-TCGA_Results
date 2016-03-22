
## Plot the expression data for 4 TCGA participants from the tumor and matched normal
## tissue to make sure expression between two are similar

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


## Select only tumor data with matched normal
tumor_with_normal <- semi_join(tumor_samples, normal_samples, by="barcode")

## data frame with participant ID for plotting
part_data <- data.frame(stringsAsFactors = FALSE)

part1_tumor <- data.frame(group = "tumor", part = "participant1", value = t(tumor_with_normal[1,-(1:2)]))
part1_normal <- data.frame(group = "normal", part = "participant1", value = t(normal_samples[1,-(1:2)]))
colnames(part1_normal)[3] <- "value"
colnames(part1_tumor)[3] <- "value"
part_data <- rbind(part1_normal, part1_tumor)

part2_tumor <- data.frame(group = "tumor", part = "participant2", value = t(tumor_with_normal[2,-(1:2)]))
part2_normal <- data.frame(group = "normal", part = "participant2", value = t(normal_samples[2,-(1:2)]))
colnames(part2_normal)[3] <- "value"
colnames(part2_tumor)[3] <- "value"
part_data <- rbind(part_data, part2_normal, part2_tumor)

part3_tumor <- data.frame(group = "tumor", part = "participant3", value = t(tumor_with_normal[3,-(1:2)]))
part3_normal <- data.frame(group = "normal", part = "participant3", value = t(normal_samples[3,-(1:2)]))
colnames(part3_normal)[3] <- "value"
colnames(part3_tumor)[3] <- "value"
part_data <- rbind(part_data, part3_normal, part3_tumor)

part4_tumor <- data.frame(group = "tumor", part = "participant4", value = t(tumor_with_normal[4,-(1:2)]))
part4_normal <- data.frame(group = "normal", part = "participant4", value = t(normal_samples[4,-(1:2)]))
colnames(part4_normal)[3] <- "value"
colnames(part4_tumor)[3] <- "value"
part_data <- rbind(part_data, part4_normal, part4_tumor)

bp <- ggplot(part_data, aes(x = part, y= value, fill = group)) + geom_boxplot()
bp <- bp + ylim(0, 5000) + xlab(NULL)
bp
