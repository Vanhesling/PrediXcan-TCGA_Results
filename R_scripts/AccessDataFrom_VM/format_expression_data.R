## Run within /mnt/Analysis/ directory
#swd("/mnt/Analysis/")

library(dplyr)

## Info from TCGA, includes other information including UUID (not same as RNAseq data), collection site etc
barcode_disease <- read.table('../TCGA_data/Participant_Info/barcode_disease.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t') 

##PrediXcan results from model with Breast Mammary Tissue
prdx_breast <- read.table('../PrediXcan_Results/TW_Breast-MammaryTissue_ElasticNet.0.5.db_results_20160226', header = TRUE, stringsAsFactors = FALSE) 

## data frame with just the barcodes and disease
BD_df <- as.data.frame(cbind(barcode_disease$Disease, barcode_disease$Barcode), stringsAsFactors = FALSE) 
colnames(BD_df)[2] <- "barcode"
colnames(BD_df)[1] <- "disease"

## This transposed and transformed header from the dosage files has all the participant barcodes in order that the appear in the PrediXcan results file 
ordered_barcodes <- read.table('/mnt/zz/new_dosages/header_files/participants.txt', header = FALSE, stringsAsFactors = FALSE) 

prdx_breast <-cbind(ordered_barcodes, prdx_breast) ## Add the barcodes to the PrediXcan file results

colnames(prdx_breast)[1] <- "barcode"
merged <- inner_join(BD_df, prdx_breast, by="barcode")

BRCA_predixcan <- filter(merged, disease == "BRCA")

gene_map <- read.table('~/PrediXcan/Software/PredictDB/GeneList/Breast_genes.csv',header = TRUE, sep = ',', stringsAsFactors = FALSE)
gene_names <- gene_map$genename
ENSG_names <- gene_map$gene

names(BRCA_predixcan)[match(ENSG_names,names(BRCA_predixcan))] <- gene_names ## match ENSG_names to gene_names

BRCA_predixcan <- BRCA_predixcan[ , !(duplicated(colnames(BRCA_predixcan)) | duplicated(colnames(BRCA_predixcan), fromLast = TRUE))] ## Remove duplicate genes!!! 

### All good up until here! Now let's get the RNAseq data in and with ONLY THE TUMOR or NORMAL DATA to start with

TCGA_RNAseq <- read.table('/mnt/TCGA_data/RNAseq/BRCA/BRCA_RNAseqV2_all.txt', header = TRUE, stringsAsFactors = FALSE)

TCGA_RNAseq <- TCGA_RNAseq[ , !(duplicated(colnames(TCGA_RNAseq)) | duplicated(colnames(TCGA_RNAseq), fromLast = TRUE))] ## Remove any duplicated genes

TCGA_RNAseq$barcode <- sapply(TCGA_RNAseq$barcode,function(xx) paste(strsplit(xx, "-")[[1]][1:4],collapse="-")) ## remove last part to make for easier grep-ing

## find just tumor/normal samples -01A/-01B = tumor; -11A/-11B = normal
normal_samples <- filter(TCGA_RNAseq, grepl('11A$|-11B$', barcode)) 
tumor_samples <- filter(TCGA_RNAseq, grepl('01A$|-01B$', barcode)) 

normal_samples$barcode <- sapply(normal_samples$barcode,function(xx) paste(strsplit(xx, "-")[[1]][1:3],collapse="-")) ## remove last part of barcode so it matches with PrediXcan barcodes
tumor_samples$barcode <- sapply(tumor_samples$barcode,function(xx) paste(strsplit(xx, "-")[[1]][1:3],collapse="-"))

prdx_normal <- semi_join(BRCA_predixcan, normal_samples, by = "barcode") ## Want just the barcodes in RNAseq
prdx_tumor <- semi_join(BRCA_predixcan, tumor_samples, by = "barcode") ## Want just the barcodes in RNAseq

 ## Put the barcodes in the same order 
prdx_normal <- arrange(prdx_normal, barcode)
normal_samples <- arrange(normal_samples, barcode)

prdx_tumor <- arrange(prdx_tumor, barcode)
tumor_samples <- arrange(tumor_samples, barcode)

