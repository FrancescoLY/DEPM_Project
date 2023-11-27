library(BiocGenerics) 
library(DESeq2)
library(psych) 
library(NetworkToolbox)
library(ggplot2);
library(GGally);library(sna);library(network)

#Part 1

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)

#Create a directory for the project
proj <- "TCGA-HNSC"
dir.create(file.path(proj))

#Cancer sample
queryC <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts", sample.type = "Primary Tumor")


GDCdownload(query = queryC, directory = "GDCdata", method = "api")

dataC <- GDCprepare(queryC, directory = "GDCdata")
expr.dataC <- assay(dataC)

View(BiocGenerics::as.data.frame(rowRanges(dataC)))
genes.info1 <- BiocGenerics::as.data.frame(rowRanges(dataC))

#Normal sample
queryN <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "STAR - Counts", sample.type = "Solid Tissue Normal")

GDCdownload(query = queryN, directory = "GDCdata", method = "api")

dataN <- GDCprepare(queryN,directory = "GDCdata")
expr.dataN <- assay(dataN)

View(BiocGenerics::as.data.frame(rowRanges(dataN)))
genes.info2 <- BiocGenerics::as.data.frame(rowRanges(dataN))

all(na.omit(genes.info2) == na.omit(genes.info1))

#Data cleaning
dim(expr.dataC)
length(unique(substr(colnames(expr.dataC), 1,12))) #no duplicates

dim(expr.dataN)
length(unique(substr(colnames(expr.dataN),1,12))) #no duplicates

patientsCancer <- substr(colnames(expr.dataC), 1,12)
sort(table(patientsCancer)) 

#We want only patients with one sample
uniquepatientsCancer <- names(which(table(patientsCancer) == 1))
idxuniqueP <- match(uniquepatientsCancer, substr(colnames(expr.dataC), 1,12) )

exprC <- as.data.frame(expr.dataC[,idxuniqueP])
exprN <- as.data.frame(expr.dataN)

#Rename the patients
colnames(exprC) <- substr(colnames(exprC), 1,12)
unique(colnames(exprC))
colnames(exprN) <- substr(colnames(exprN), 1,12)
unique(colnames(exprN))

#There is a normal sample that do not have a cancerous sample. We have to remove it.
intersect(colnames(exprN), colnames(exprC)) 
setdiff(colnames(exprN), colnames(exprC))

#idx of the patient to remove
match(setdiff(colnames(exprN), colnames(exprC)), colnames(exprN))

colnames(exprN)[21]
