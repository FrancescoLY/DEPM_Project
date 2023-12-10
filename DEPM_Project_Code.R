library(TCGAbiolinks)
library(SummarizedExperiment)
library(psych)
library(DescTools)
library(network)
library(igraph)
library(GGally)

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


# 3. Co-expression networks -----------------------------------------------

# Using Pearson correlation coefficient, Bonferroni correction, 0.7 threshold

filtr.expr.c.deg <- filtr.expr.c[deg.genes,]
# dim(filtr.expr.c.deg) # (650x43)
cor.mat.c <- corr.test(t(filtr.expr.c.deg), use = "pairwise", 
                       method = "pearson",adjust="bonferroni",ci=FALSE)
rho.c <- cor.mat.c$r
# dim(rho.c)            # (650x650) 
diag(rho.c) <- 0

filtr.expr.n.deg <- filtr.expr.n[deg.genes,]
# dim(filtr.expr.n.deg) # (650x43)
cor.mat.n <- corr.test(t(filtr.expr.n.deg), use = "pairwise", 
                       method = "pearson",adjust="bonferroni",ci=FALSE)
rho.n <- cor.mat.n$r
# dim(rho.n)            # (650x650) 
diag(rho.n) <- 0


# 4. Differential Co-expressed Network  -----------------------------------

fisher.z.c <- FisherZ(rho.c)
# dim(fisher.z.c)

fisher.z.n <- FisherZ(rho.n)
# dim(fisher.z.n)

# Assuming both fractions in the formula denominator are under the same 
# square root
sample_size = 43
Z.scores <- (fisher.z.c-fisher.z.n)/sqrt(1/(sample_size-3)+1/(sample_size-3))
# dim(Z.scores)
# isSymmetric(Z.scores) 
diag(Z.scores) <- 0
adj.mat.diff <- Z.scores * (abs(Z.scores) >= 3)
# Check that all elements with |Z| < 3 are removed
j = 0
for (i in adj.mat.diff) {
  # j <- j+1
  if (i != 0 && abs(i) < 3) {
    j <- j+1
  }
}
j  # == 0
# adj.mat.diff <- matrix(0, ncol = 650, nrow = 650)
adj.mat.diff[adj.mat.diff != 0] <- 1
# sum(diag(adj.mat.diff))                                             # == 0
# length(which(adj.mat.diff!=0)) - length(which(filtr.Z.scores!=0))   # == 0
# Check all nonzero elements are 1
j = 0
k = 0
for (i in adj.mat.diff) {
  # j <- j+1
  if (i != 0) {
    j <- j+1
    k <- k+i
  }
}
j-k  # == 0
net.diff <- network(adj.mat.diff, matrix.type="adjacency", directed = F)
# network.density(net.diff)
digree.diff <- rowSums(adj.mat.diff != 0)
names(degree.diff) <- rownames(adj.mat.diff)
hist(digree.diff)
hist(digree.diff, breaks = 20)
hist(digree.diff, breaks = 40)
hist(digree.diff, breaks = 80)
hist(digree.diff, breaks = 160)
hist(digree.diff, breaks = network.size(net.diff))
# plot(table(degree.diff))
degree.diff.freq <- as.data.frame(table(degree.diff))
# # colnames(degree.diff.freq)
plot(as.numeric(degree.diff.freq$degree.diff),
     as.numeric(degree.diff.freq$Freq), log = "xy", type='p')

# # Fitting a power law over the histogram
# power.law.params <- power.law.fit(digree.diff)
# hist(digree.diff, prob=TRUE, breaks = 80)
# fun <- function(x) {
#   x**(-1*power.law.params$alpha)
# }
# x2 <- seq(min(digree.diff), max(digree.diff), length = 80)
# curve(fun, from = 157, to = 250, add= TRUE, lwd = 2)     

# Assuming the network is a scale free network
# sum(degree.diff == 0)
# q <- quantile(digree.diff[digree.diff>0],0.95)
q <- quantile(digree.diff,0.95)
# q
# hist(degree.diff)
# abline(v=q, col="red")

hubs.diff <- degree.diff[degree.diff>=q]
names(hubs.diff)

# hubs.diff <- sort(hubs.diff, decreasing = F)
# hubs.diff <- hubs.diff[1]

hubs.diff.ids <- vector("integer",length(hubs.diff))
for (i in 1:length(hubs.diff)){hubs.diff.ids[i] <- match(names(hubs.diff)[i],
                                                      rownames(adj.mat.diff))}
hubs.diff.ids
#identifying the neighborhood of the hubs
hubs.diff.neigh.ids <- c()
for (hub in hubs.diff.ids){
  hubs.diff.neigh.ids <- append(hubs.diff.neigh.ids, 
                                get.neighborhood(net.diff, hub))
}
hubs.diff.neigh.ids <- unique(hubs.diff.neigh.ids)
hubs.diff.neigh.ids
hubs.diff.neigh.names <- rownames(adj.mat.diff[hubs.diff.neigh.ids,])
subnet.diff <- unique(c(names(hubs.diff), hubs.diff.neigh.names))
hubs.diff.adj <- adj.mat.diff[subnet.diff, subnet.diff]
names.hubs.diff <-names(hubs.diff)

net.hubs.diff <- network(hubs.diff.adj, matrix.type="adjacency")
net.hubs.diff %v% "type" = ifelse(network.vertex.names(net.hubs.diff) %in% 
                                    names.hubs.diff,"hub", "non-hub")
net.hubs.diff %v% "color" = ifelse(net.hubs.diff %v% "type" == "hub", 
                                   "tomato", "deepskyblue3")
net.hubs.diff %v% "color" <- sort(net.hubs.diff %v% "color")

# ggnet2(net.hubs.diff,  color = "color", size = "degree", # max_size=250,
#        edge.alpha = 0.7,  edge.size = 0.1,
#        # node.label = names.hubs.diff, label.color = "black", # label.size = 4
#        ) + guides(size = "none")

ggnet2(net.hubs.diff, color = "color", alpha = 0.7, size = 2,
       # edge.color = "edgecolor", 
       edge.alpha = 1, edge.size = 0.15)+
  guides(size = "none") 
