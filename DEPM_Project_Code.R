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
rna.data.All.colData.C_and_N <- rna.data.All.colData[rna.data.All.colData$patient %in% patients_C_and_N, ]
barcodes_C_and_N <- rna.data.All.colData.C_and_N$barcode    # 86
rna.expr.data.All.C_and_N <- rna.expr.data.All[,barcodes_C_and_N]
dim(rna.expr.data.All.C_and_N)                              # (60660, 86)