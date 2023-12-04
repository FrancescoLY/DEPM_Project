library(TCGAbiolinks)
library(SummarizedExperiment)

proj <- "TCGA-HNSC"
dir.create(file.path(proj))

rna.query.All <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts")
GDCdownload(query = rna.query.All, directory = "GDCdata_All", method = "api")
rna.data.All <- GDCprepare(rna.query.All, directory = "GDCdata_All")

rna.expr.data.All <- assay(rna.data.All)
View(BiocGenerics::as.data.frame(rowRanges(rna.data.All)))
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.All))

rna.data.All.colData <- BiocGenerics::as.data.frame(rna.data.All@colData)
rna.data.All.colData.N = rna.data.All.colData[rna.data.All.colData$sample_type == "Solid Tissue Normal", ]
rna.data.All.colData.C = rna.data.All.colData[rna.data.All.colData$sample_type %in% c("Primary Tumor", "Metastatic"), ]
patients_C <- rna.data.All.colData.C$patient
patients_N <- rna.data.All.colData.N$patient
patients_C_and_N <- intersect(patients_C, patients_N)       # 43
rna.data.All.colData.C_2 <- rna.data.All.colData.C[rna.data.All.colData.C$patient %in% patients_C_and_N, ]
barcodes_C <- rna.data.All.colData.C_2$barcode              # 43
rna.data.All.colData.N_2 <- rna.data.All.colData.N[rna.data.All.colData.N$patient %in% patients_C_and_N, ]
barcodes_N <- rna.data.All.colData.N_2$barcode              # 43
rna.expr.data.All.C <- rna.expr.data.All[,barcodes_C]
dim(rna.expr.data.All.C)                                    # (60660, 43)
colnames(rna.expr.data.All.C) <- substr(colnames(rna.expr.data.All.C), 1,12)
rna.expr.data.All.N <- rna.expr.data.All[,barcodes_N]
dim(rna.expr.data.All.N)                                    # (60660, 43)
colnames(rna.expr.data.All.N) <- substr(colnames(rna.expr.data.All.N), 1,12)
all(rownames(rna.expr.data.All.C) == rownames(rna.expr.data.All.N))
full.data <- cbind(rna.expr.data.All.N, rna.expr.data.All.C)
dim(full.data)                                              # (60660, 86)
for (col in colnames(full.data)) {                          # Checking the column partition
  a = which(colnames(full.data)==col)                       
  print(paste(a[1]<=43, a[2] > 43))
}
full.data.df <- data.frame(full.data)
metad <- rep("cancer", 86)
metad[1:43] <- "normal"
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data.df)
colnames(metad)[1] <- "condition"
# metad[,1]
metad[,1] <- as.factor(metad[,1])
full.data.df <- cbind(rownames(full.data.df), full.data.df)

dds <- DESeqDataSetFromMatrix(countData=full.data.df, 
                              colData=metad, 
                              design= ~condition,
                              tidy=TRUE)

# View(counts(dds))
dim(counts(dds))                                              # (60660, 86)

# filtering: at least ten counts on 90% of patients
ninety_precent_of_patients = (86*90)/100
keep <- rowSums(counts(dds) >= 10) >= ninety_precent_of_patients
dds <- dds[keep,]
dim(counts(dds))                                              # (14778, 86)

# filtering: removing rows with at least one 0 value
keep <- rowSums(counts(dds) == 0) == 0
dds <- dds[keep,]
dim(counts(dds))                                              # (14570, 86)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == 43)                    # no null rows
filtr.expr.n <- as.data.frame(normalized_counts[, 1:43])
filtr.expr.c <- as.data.frame(normalized_counts[, 44:86])

# renaming cancer expression data columns
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1,12)  

# Identifying DEGs
fc <-  log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n) ) 

# names(fc)
# a <- rownames(filtr.expr.c)
# identical(names(fc),a)
# names(fc) <- rownames(filtr.expr.c)
# head(fc)

# t.test(filtr.expr.c[3,], filtr.expr.n[3,] )   # Welch Two Sample t-test

a=0                                             # Checking row ordering
for (i in 1:nrow(filtr.expr.n)) {                          
  if (rownames(filtr.expr.n)[i]!=rownames(filtr.expr.c)[i]) {
    a <- a+1
  }
}
a   # a = 0

# pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(filtr.expr.c[i,], filtr.expr.n[i,] ))$p.value)
# # With method="fdr", p.adjust() implements the Benjamini–Hochberg correction
# pval.fc.fdr <- p.adjust(pval.fc, method="fdr")
# expr.table <- data.frame(cbind(fc, pval.fc.fdr))
# # expr.table[,1] <- round(expr.table[,1],2)       # Round FC to second decimal digit
# 
# deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.2 & expr.table$pval.fc.fdr <=0.05,])
# length(deg.genes)   # 1590

# This time we try Bonferroni correction instead
pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(filtr.expr.c[i,], filtr.expr.n[i,] ))$p.value)
pval.fc.fdr <- p.adjust(pval.fc, method="bonferroni")
expr.table <- data.frame(cbind(fc, pval.fc.fdr))
# expr.table[,1] <- round(expr.table[,1],2)       # Round FC to second decimal digit

deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.2 & expr.table$pval.fc.fdr <=0.05,])
length(deg.genes)   # 650

# # This time we try Benjamini–Yekutieli correction instead
# pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(filtr.expr.c[i,], filtr.expr.n[i,] ))$p.value)
# pval.fc.fdr <- p.adjust(pval.fc, method="BY")
# expr.table <- data.frame(cbind(fc, pval.fc.fdr))
# # expr.table[,1] <- round(expr.table[,1],2)       # Round FC to second decimal digit
# 
# deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.2 & expr.table$pval.fc.fdr <=0.05,])
# length(deg.genes)   # 1238