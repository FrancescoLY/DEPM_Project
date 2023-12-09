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
exprN <- exprN[,-c(25)]

#We check if there are null values
any(is.na(exprC)) 
any(is.nan(as.matrix(exprC)))

any(is.na(exprN)) 
any(is.nan(as.matrix(exprN)))

exprC <- exprC[, colnames(exprN)]

clinical.query<- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)

#Normalize data
all(rownames(exprC) == rownames(exprN))
full.data <- cbind(exprN, exprC)

dim(full.data)
full.data <- data.frame(full.data)


colData <- rep("cancer",86)
colData[1:43] <- "normal"
colData
colData <- data.frame(TCGA_id=colnames(full.data))
colnames(colData)[1] <- "condition"
colData[,1] <- as.factor(colData[,1])
full.data <- cbind(rownames(full.data), full.data)

dds <- DESeqDataSetFromMatrix(countData=full.data, 
                              colData=colData, 
                              design= ~condition,
                              tidy=TRUE)

View(counts(dds))
dim(counts(dds))

# filtering: at least ten counts on 90% of patients
(86*90)/100
keep <- rowSums(counts(dds) >= 10) >= 78
dds <- dds[keep,]
dim(counts(dds))

#delete the row with 0 values
keep <- rowSums(counts(dds) == 0) == 0
dds <- dds[keep,]
dim(counts(dds)) 

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == 43) #no null rows

filtr_expr_n <- as.data.frame(normalized_counts[, 1:43])
filtr_expr_c <- as.data.frame(normalized_counts[, 44:86])
colnames(filtr_expr_c) <- substr(colnames(filtr_expr_c), 1,12)


#DEGs

#fold change
fc <-  log2(rowMeans(filtr_expr_c) / rowMeans(filtr_expr_n) ) 
names(fc) <- rownames(filtr_expr_c)
head(fc)

#p-value
pval <- sapply(1:nrow(filtr_expr_c), function(i) (t.test(filtr_expr_c[i,], filtr_expr_n[i,] ))$p.value)
pval_adj <- p.adjust(pval, method="bonferroni")

expr.table <- data.frame(cbind(fc, pval_adj))
expr.table[,1] <- round(expr.table[,1],2)

deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 1.2 & expr.table$pval_adj <=0.05,]) 
deg.genes


#---Part 3---
#Co-expression networks

#Remove non-deg genes
idxc <- c(match(setdiff(rownames(filtr_expr_c), deg.genes), rownames(filtr_expr_c)))
filtr_expr_c <- filtr_expr_c[-c(idxc),]
View(filtr_expr_c)

idxn <- c(match(setdiff(rownames(filtr_expr_n),deg.genes),rownames(filtr_expr_n)))
filtr_expr_n <- filtr_expr_n[-c(idxn),]
View(filtr_expr_n)

#cancer adjacency matrix
corr.c <- corr.test(t(filtr_expr_c), use = "pairwise", 
                       method = "pearson",adjust="holm",ci=FALSE)

rho.c <- corr.c$r
diag(rho.c) <- 0
qval.c <- corr.c$p
#qvals are reported on the upper triangle only
qval.c[lower.tri(qval.c)] <- t(qval.c)[lower.tri(qval.c)]

adj.mat.c <- rho.c * (qval.c <= 0.7)

#normal adjacency matrix
corr.n <- corr.test(t(filtr_expr_n),use = "pairwise", method = "pearson", adjust = "holm",ci=FALSE)

rho.n <- corr.n$r
diag(rho.n) <- 0

qval.n <- corr.n$p
qval.n[lower.tri(qval.n)] <- t(qval.n)[lower.tri(qval.n)]

adj.mat.n <- rho.n*(qval.n <= 0.7)

#cancer network
net.c <- network(adj.mat.c, matrix.type="adjacency",ignore.eval = FALSE, 
                 names.eval = "weights", directed = F)

network.density(net.c)
network.size(net.c)
network.edgecount(net.c)
clustcoeff(adj.mat.c, weighted = FALSE)$CC

sum(adj.mat.c != 0)
#how many positive/negative correlations? 
sum(adj.mat.c > 0) 
sum(adj.mat.c < 0) 

degree.c <- rowSums(adj.mat.c != 0)
names(degree.c) <- rownames(adj.mat.c)
degree.c <- sort(degree.c,decreasing = TRUE)
head(degree.c,10)
sum(degree.c == 0) #unconnected nodes 

hist(degree.c)
x <- quantile(degree.c[degree.c>0],0.95) #how big is the degree of the most connected nodes?
x
hist(degree.c)
abline(v=x, col="red")

hubs.c <- degree.c[degree.c>=x]
names(hubs.c) #24 hubs

net.c %v% "type" = ifelse(network.vertex.names(net.c) %in% names(hubs.c),"hub", "non-hub")
net.c %v% "color" = ifelse(net.c %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.c, "edgecolor", ifelse(net.c %e% "weights" > 0, "red", "blue"))

ggnet2(net.c, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 


#normal network
net.n <- network(adj.mat.n, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

network.density(net.n)
network.size(net.n)
network.edgecount(net.n)
clustcoeff(adj.mat.n, weighted = FALSE)$CC


sum(adj.mat.n != 0)
#how many positive/negative correlations? 
sum(adj.mat.n > 0) 
sum(adj.mat.n < 0) 

degree.n <- rowSums(adj.mat.n != 0)
names(degree.n) <- rownames(adj.mat.n)
degree.n <- sort(degree.n, decreasing = TRUE)
head(degree.n,10)
sum(degree.n == 0) #unconnected nodes 

hist(degree.n)
sampledata <- data.frame(x = degree.n, y = 1:654)
plot_data <- data.frame(x= log(sampledata$x), y= log(sampledata$y))

ggplot(plot_data, aes(x=x,y=y)) + 
  geom_point()+ labs(title = 'log-log plot', x = 'degree', y = 'frequency')

y <- quantile(degree.n[degree.n>0],0.95) #how big is the degree of the most connected nodes?
y
hist(degree.n)
abline(v=y, col="red")

hubs.n <- degree.n[degree.n>=y]
names(hubs.n)#30 hubs

net.n %v% "type" = ifelse(network.vertex.names(net.n) %in% names(hubs.n),"hub", "non-hub")
net.n %v% "color" = ifelse(net.n %v% "type" == "hub", "tomato", "deepskyblue3")
set.edge.attribute(net.n, "edgecolor", ifelse(net.n %e% "weights" > 0, "red", "blue"))

ggnet2(net.n, color = "color", alpha = 0.7, size = 2,
       edge.color = "edgecolor", edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 


intersect(names(hubs.c),names(hubs.n))
