# R version 3.6.0

## TODO:
# modify plots to show p-value and change above a line connecting the two relevant normalized count

####  load/install required packages  ####
packages = c("gplots","RColorBrewer","BiocManager","utils", "ggplot2", "ggrepel", "factoextra")
package.check <- function(x){
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
    }
}
for(pkg in packages) {
  package.check(pkg)
}
if(!require('DESeq2',character.only = TRUE)) {
  BiocManager::install("DESeq2")
  library("DESeq2")
}

####  data preprocessing ####
# set working directory to script directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load gene id and gene names
load("r.geneid.ws271")

# load raw count data generated by featureCount
ocyts <- read.delim("./oocyte", comment.char="#", stringsAsFactors=FALSE)

# remove Chr, Start, End, Length, Strand attributes
ocyts <- ocyts[,-c(2:6)]

# remove .cutadapt.sam from sample names
colnames(ocyts) <- c('geneid',gsub(".*(m..c.).*","\\1",colnames(ocyts)[2:ncol(ocyts)]))

# save geneid in temporary vector ids
ids <- ocyts$geneid

# delete column of geneid now, format matrix into integer matrix
ocyts <- ocyts[,-1]
ocyts <- apply(ocyts,2,as.integer)

# set geneid as rownames
row.names(ocyts) <- ids

# ocyts with row names
ocytsR <- ocyts

temp <- geneid[row.names(ocytsR),]
row.names(ocytsR) <- temp$gene.name
ocytsR <- cbind(ocytsR, temp$gene.name)

####  basic sanity check  ####
# check read depth for each cell
# conclusion: all cells got >8M reads
count.in.millions <- colSums(ocyts)/1e6
par(las=2)
par(mar=c(5,5,5,5))
barplot(count.in.millions, main="Counts per Oocyte", xlab="Oocyte", ylab="Counts (1e6)")
par(mfrow=c(1,1))
hist(count.in.millions)

# check if read depth is sufficient to detect all genes
# conclusion: read depth is sufficient, ~8000 genes detected at a threshold of 1 CPM
# calculate CPM
cpm <- t(t(ocyts)/count.in.millions)

# number of gene detected by different threshold
th.100 <- colSums(cpm>100)
th.10 <- colSums(cpm>10)
th.1 <- colSums(cpm>1)
th.0.1 <- colSums(cpm>0.1)
temp1 <- data.frame(gene.num=th.100,threshold=100)
rownames(temp1) <- c()
temp2 <- data.frame(gene.num=th.10,threshold=10)
rownames(temp2) <- c()
temp3 <- data.frame(gene.num=th.1,threshold=1)
rownames(temp3) <- c()
temp4 <- data.frame(gene.num=th.0.1,threshold=0.1)
rownames(temp4) <- c()
temp <- rbind(temp1,temp2,temp3,temp4)
par(mfrow=c(1,1))
boxplot(gene.num~threshold,data = temp, xlab='threshold (cpm)', ylab='gene detected', main='gene detection by threshold')

# check if read depth is sufficient
# when read depth is insufficient, the number of gene detected increase with count, especially for low detection threshold, e.g. 0.1 CPM
par(mfrow=c(1,4))
plot(colSums(ocyts)/1e6,th.100,xlab='sequencing depth (M)',ylab='gene detected', main='100 cpm')
plot(colSums(ocyts)/1e6,th.10,xlab='sequencing depth (M)',ylab='gene detected',main='10 cpm')
plot(colSums(ocyts)/1e6,th.1,xlab='sequencing depth (M)',ylab='gene detected',main='10 cpm')
plot(colSums(ocyts)/1e6,th.0.1,xlab='sequencing depth (M)',ylab='gene detected',main='0.1 cpm')

#### Unfiltered Unsupervised clustering #### KAT
# Remove any rows in which the mean would be 0
ocytsP <- ocyts[rowMeans(ocyts) > 0,]

# Transpose the counts matrix to enable PCA on a per oocyte basis
pca_res <- prcomp(ocytsP, scale. = TRUE)

biplot(pca_res)

# Transpose counts matrix for k-means analysis
ocytsT <- as.data.frame(t(ocytsP))

set.seed(123)
ocyts.km.res <- kmeans(ocytsT, 4, nstart=25)
ocytsT$ocyt <- c('m1','m1','m1','m1','m3','m3','m3','m3','m5','m5','m5','m5')
plot(ocyts.km.res, data=ocytsT, class=NULL)



####  differential expression analysis ####
## get normalized counts
# basic filtering, keep genes with >2 counts in at least 8 (out of 12) oocytes
cts <- ocyts[rowSums(ocyts>2)>8,]

# setup design matrix for DESeq2
coldata <- as.matrix(substr(colnames(cts),1,2))
colnames(coldata) <- c('condition')
rownames(coldata) <- colnames(cts)
summary(coldata)
# run DESeq2, should take <2 minutes
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)

# extract normalized counts
cd <- counts(dds,normalize=TRUE)
temp <- geneid[row.names(cd),]
row.names(cd) <- temp$gene.name

## do differential expression between -5 and -3
# basic filtering, keep genes with >2 counts in at least 8 (out of 12) oocytes
cts <- ocyts[rowSums(ocyts>2)>8,]
# remove columes containing -3 oocytes
# note: in theory, DESeq2 can take >2 groups for differential expression analysis, but in reality it never works well.
# therefore, we only compare 2 groups at a time. -5 vs -1 first, then -3 vs -1.
cts <- cts[,c(-1,-2,-3,-4)]

# setup design matrix for DESeq2
coldata <- as.matrix(substr(colnames(cts),1,2))
colnames(coldata) <- c('condition')
rownames(coldata) <- colnames(cts)
summary(coldata)
# run DESeq2, should take <2 minutes
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
# extract results
res <- results(dds)
# get DEGs
# order by adjusted p value
res <- res[order(res$padj),]
res <- as.data.frame(res)
res <- res[complete.cases(res),]
m3.m5.all <- res
rownames(m3.m5.all) <- geneid[rownames(m3.m5.all),]$gene.name
# only keep significantly upregulated genes in -1
keep <- (res$pvalue<0.05) & (res$log2FoldChange < -1)
# note: these DEGs are selected only based on -5 and -1
m3.m5.degs <- res[keep,]
m3.m5.degs <- m3.m5.degs[complete.cases(m3.m5.degs),]
rownames(m3.m5.degs) <- geneid[rownames(m3.m5.degs),]$gene.name
m3.m5.degs$test <- '-3 vs -5'


## do differential expression between -5 and -1
# basic filtering, keep genes with >2 counts in at least 8 (out of 12) oocytes
cts <- ocyts[rowSums(ocyts>2)>8,]
# remove columes containing -3 oocytes
# note: on theory, DESeq2 can take >2 groups for differential expression analysis, but in fact it never works well.
# therefore, we only compare 2 groups at a time. -5 vs -1 first, then -3 vs -1.
cts <- cts[,c(-5,-6,-7,-8)]
# setup design matrix for DESeq2
coldata <- as.matrix(substr(colnames(cts),1,2))
colnames(coldata) <- c('condition')
rownames(coldata) <- colnames(cts)
summary(coldata)
# run DESeq2, should take <2 minutes
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
# extract results
res <- results(dds)
# get DEGs
# order by adjusted p value
res <- res[order(res$padj),]
res <- as.data.frame(res)
res <- res[complete.cases(res),]
m1.m5.all <- res
rownames(m1.m5.all) <- geneid[rownames(m1.m5.all),]$gene.name
# only keep significantly upregulated genes in -1
keep <- (res$pvalue<0.05) & (res$log2FoldChange < -1)
# note: these DEGs are selected only based on -5 and -1
m1.m5.degs <- res[keep,]
m1.m5.degs <- m1.m5.degs[complete.cases(m1.m5.degs),]
rownames(m1.m5.degs) <- geneid[rownames(m1.m5.degs),]$gene.name
m1.m5.degs$test <- '-3 vs -5'

## m1.m5.degs contain genes upregulated in -1 compared to -5
## there could also be genes upregulated in -1 compared to -3, but have similar levels between -5 and -3
## then, do differential expression between -3 and -1
# basic filtering, keep genes with >2 counts in at least 8 (out of 12) oocytes
cts <- ocyts[rowSums(ocyts>2)>8,]
# remove column containing -5 oocytes
cts <- cts[,c(-9,-10,-11,-12)]
# setup design matrix for DESeq2
coldata <- as.matrix(substr(colnames(cts),1,2))
colnames(coldata) <- c('condition')
rownames(coldata) <- colnames(cts)
summary(coldata)
# run DESeq2, should take <2 minutes
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
# extract results
res <- results(dds)
# get DEGs
# order by adjusted p value
res <- res[order(res$padj),]
res <- as.data.frame(res)
res <- res[complete.cases(res),]
m1.m3.all <- res
rownames(m1.m3.all) <- geneid[rownames(m1.m3.all),]$gene.name
# only keep significantly upregulated genes in -1
keep <- (res$pvalue<0.05) & (res$log2FoldChange < -1)
# note: these DEGs are selected only based on -3 and -1
m1.m3.degs <- res[keep,]
m1.m3.degs <- m1.m3.degs[complete.cases(m1.m3.degs),]
rownames(m1.m3.degs) <- geneid[rownames(m1.m3.degs),]$gene.name
m1.m3.degs$test <- '-1 vs -3'

# combine these 2 together to get gene of interest
# to combine
toremove <- intersect(rownames(m1.m5.degs),rownames(m1.m3.degs))
m3.m5.intersect <- m1.m3.degs[rownames(m1.m3.degs)%in%toremove,]
m1.m3.degs <- m1.m3.degs[!rownames(m1.m3.degs)%in%toremove,]
GOI <- unique(rbind(m1.m3.degs,m1.m5.degs))
GOI <- GOI[order(GOI$pvalue),]




#### figure style visualizations
### volcano plot of each comparison (m3 or m5 to m1)
## plot for m1.m3.all. 
# Since a negative fold change reflects an increase when going from -3 -> -1,
# Negate all of log2FoldChange so that positive reflects an increase when going from -3 -> -1
m1.m3.all$log2FoldChange <- -(m1.m3.all$log2FoldChange)
#Add column for differential expression
m1.m3.all$diffexpressed <- 'NO'
#if log2FoldChange > 1, pvalue < 0.05 set as 'up'
m1.m3.all$diffexpressed[m1.m3.all$log2FoldChange > 1 & m1.m3.all$padj < 0.05] <- "UP"
#if log2FoldChange < -1, pvalue < 0.05 set as 'down'
m1.m3.all$diffexpressed[m1.m3.all$log2FoldChange < -1 & m1.m3.all$padj < 0.05] <- "DOWN"

# Set a label column
m1.m3.all$delabel <- NA
m1.m3.all$delabel[m1.m3.all$diffexpressed != "NO"] <- row.names(m1.m3.all)[m1.m3.all$diffexpressed != "NO"]
# make initial volcano plot
m3plot <- ggplot(data=m1.m3.all, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) + geom_point(size=0.5)
# print gene names as text
m3plot <- m3plot + theme_bw() #+ geom_label_repel()
# set lines for volcano plot
m3plot <- m3plot + geom_vline(xintercept=c(-1, 1),) + geom_hline(yintercept=-log10(0.05))
# set colors for volcano plot ("Down", "No", "UP")
m3plot <- m3plot + scale_color_manual(values=c("blue", "black", "red")) + labs(title='m1 vs m3')
# set axis scale numbers
m3plot <- m3plot + scale_y_continuous(breaks=seq(0,4,by=1)) + scale_x_continuous(breaks=seq(-6,6,by=1))
# print out the plot
m3plot

## Do it again for m1.m5
# Negate all of log2FoldChange so that positive reflects an increase when going from -5 -> -1
m1.m5.all$log2FoldChange <- -(m1.m5.all$log2FoldChange)
#Add column for differential expression
m1.m5.all$diffexpressed <- 'NO'
#if log2FoldChange > 1, pvalue < 0.05 set as 'up'
m1.m5.all$diffexpressed[m1.m5.all$log2FoldChange > 1 & m1.m5.all$padj < 0.05] <- "UP"
#if log2FoldChange < -1, pvalue < 0.05 set as 'down'
m1.m5.all$diffexpressed[m1.m5.all$log2FoldChange < -1 & m1.m5.all$padj < 0.05] <- "DOWN"
# Set a label column
m1.m5.all$delabel <- NA
m1.m5.all$delabel[m1.m5.all$diffexpressed != "NO"] <- row.names(m1.m5.all)[m1.m5.all$diffexpressed != "NO"]
# make initial volcano plot

m5plot <- ggplot(data=m1.m5.all, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) + geom_point(size=0.5)
# print gene names as text
m5plot <- m5plot + theme_bw()
# set lines for volcano plot
m5plot <- m5plot + geom_vline(xintercept=c(-1, 1),) + geom_hline(yintercept=-log10(0.05))
# set colors for volcano plot ("Down", "No", "UP")
m5plot <- m5plot + scale_color_manual(values=c("blue", "black", "red")) + labs(title="m5")
# set axis scale numbers
m5plot <- m5plot + scale_y_continuous(breaks=seq(0,20,by=1)) + scale_x_continuous(breaks=seq(-6,6,by=1))
# print out the plot
m5plot


### Visualizations for oocyte read data ###
## From Above ##
####  basic sanity check  ####
# check read depth for each cell
# conclusion: all cells got >8M reads
oocyte.read.counts.in.millions <- data.frame(colnames(ocyts),colSums(ocyts)/1e6)
names(oocyte.read.counts.in.millions) <- c('name','counts')
ocyts.read.counts.plot <- ggplot(data=oocyte.read.counts.in.millions, aes(x=name, y=counts)) + geom_bar(stat='identity')
ocyts.read.counts.plot <- ocyts.read.counts.plot + theme_bw()
ocyts.read.counts.plot

#barplot(count.in.millions, main="Counts per Oocyte", xlab="Oocyte", ylab="Counts (1e6)")

# check if read depth is sufficient to detect all genes
# conclusion: read depth is sufficient, ~8000 genes detected at a threshold of 1 CPM
# calculate CPM
count.in.millions <- colSums(ocyts)/1e6
cpm <- t(t(ocyts)/count.in.millions)
# number of gene detected by different threshold
th.100 <- colSums(cpm>100)
th.10 <- colSums(cpm>10)
th.1 <- colSums(cpm>1)
th.0.1 <- colSums(cpm>0.1)
temp1 <- data.frame(colnames(ocyts),gene.num=th.100,oocyte.read.counts.in.millions$counts,rep('> 100 cpm'))
rownames(temp1) <- c()
names(temp1)[4] <- 'threshold'
temp2 <- data.frame(colnames(ocyts),gene.num=th.10,oocyte.read.counts.in.millions$counts,rep('> 10 cpm'))
rownames(temp2) <- c()
names(temp2)[4] <- 'threshold'
temp3 <- data.frame(colnames(ocyts),gene.num=th.1,oocyte.read.counts.in.millions$counts,rep('> 1 cpm'))
rownames(temp3) <- c()
names(temp3)[4] <- 'threshold'
temp4 <- data.frame(colnames(ocyts),gene.num=th.0.1,oocyte.read.counts.in.millions$counts,rep('> 0.1 cpm'))
rownames(temp4) <- c()
names(temp4)[4] <- 'threshold'
temp <- rbind(temp1,temp2,temp3,temp4)
names(temp) <- c('oocyte', 'gene.count', 'read.count', 'threshold')


ocyts.threshold.gene.counts <- ggplot(data=temp, aes(x=threshold, y=gene.count, group=threshold)) + geom_boxplot(outlier.color = 'black')
ocyts.threshold.gene.counts <- ocyts.threshold.gene.counts + facet_grid(cols=vars(threshold), scales='free_x') + scale_x_discrete()
ocyts.threshold.gene.counts <- ocyts.threshold.gene.counts + theme_bw()
ocyts.threshold.gene.counts

#boxplot(gene.num~threshold,data = temp, xlab='threshold (cpm)', ylab='gene detected', main='gene detection by threshold')

# check if read depth is sufficient
# when read depth is insufficient, the number of gene detected increase with count, especially for low detection threshold, e.g. 0.1 CPM

sequencing.depth.plot <- ggplot(data=temp, aes(x=read.count, y=gene.count)) + geom_point() + facet_wrap(. ~ threshold, scales='free')
sequencing.depth.plot <- sequencing.depth.plot + theme_bw()
sequencing.depth.plot

#par(mfrow=c(1,4))
#plot(colSums(ocyts)/1e6,th.100,xlab='sequencing depth (M)',ylab='gene detected', main='100 cpm')
#plot(colSums(ocyts)/1e6,th.10,xlab='sequencing depth (M)',ylab='gene detected',main='10 cpm')
#plot(colSums(ocyts)/1e6,th.1,xlab='sequencing depth (M)',ylab='gene detected',main='10 cpm')
#plot(colSums(ocyts)/1e6,th.0.1,xlab='sequencing depth (M)',ylab='gene detected',main='0.1 cpm')

# Look at specific row values:
m1.m5.all["alg-5",]

m1.thresh.10 <- cd[,c(1,2,3,4)]
m1.thresh.10 <- m1.thresh.10[min(m1.thresh.10[1],m1.thresh.10[2],m1.thresh.10[3],m1.thresh.10[4]) > 10,]
m3.thresh.10 <- cd[,c(5,6,7,8)]
m3.thresh.10 <- m3.thresh.10[min(m3.thresh.10[1],m3.thresh.10[2],m3.thresh.10[3],m3.thresh.10[4]) > 10,]
m5.thresh.10 <- cd[,c(9,10,11,12)]
m5.thresh.10 <- m5.thresh.10[min(m5.thresh.10[1],m5.thresh.10[2],m5.thresh.10[3],m5.thresh.10[4]) > 10,]

m1.no.counts <- cd[min(m1.thresh.10[1],m1.thresh.10[2],m1.thresh.10[3],m1.thresh.10[4])==0,c(1,2,3,4)]

# Separate out UP and DOWN expression data for each group

m1.m5.up <- m1.m5.all[m1.m5.all$diffexpressed == 'UP',]
m1.m5.down <- m1.m5.all[m1.m5.all$diffexpressed == 'DOWN',]

# Extract data for each set and write to csv file
all.up.frame <- data.frame(
  gene_name = m1.m5.up,
  m1m5_fold_change = m1.m5.all[m1.m5.up,2],
  m1m5_pvalue = m1.m5.all[m1.m5.up,5],
  stringsAsFactors = FALSE
)
write.csv(all.up.frame, 'm1m5_up.csv', row.names=FALSE)

all.down.frame <- data.frame(
  gene_name = m1.m5.down,
  m1m5_fold_change = m1.m5.all[m1.m5.down,2],
  m1m5_pvalue = m1.m5.all[m1.m5.down,5],
  stringsAsFactors = FALSE
)
write.csv(all.down.frame, 'm1m5_down.csv', row.names=FALSE)

# write csv files for all of the DESEQ2 results for -1 vs -5
m1.m5.all.out <- data.frame(
  geneName = row.names(m1.m5.all),
  baseMean = m1.m5.all[,1],
  log2FoldChange = m1.m5.all[,2],
  lfcSE = m1.m5.all[,3],
  stat = m1.m5.all[,4],
  pvalue = m1.m5.all[,5],
  padj = m1.m5.all[,6],
  stringsAsFactors = FALSE
)
write.csv(m1.m5.all.out, 'deseq_m1m5.csv', row.names=FALSE)

# write csv files for all of the DESEQ2 results for -1 vs -3
m1.m3.all.out <- data.frame(
  geneName = row.names(m1.m3.all),
  baseMean = m1.m3.all[,1],
  log2FoldChange = m1.m3.all[,2],
  lfcSE = m1.m3.all[,3],
  stat = m1.m3.all[,4],
  pvalue = m1.m3.all[,5],
  padj = m1.m3.all[,6],
  stringsAsFactors = FALSE
)
write.csv(m1.m3.all.out, 'deseq_m1m3.csv', row.names=FALSE)

# write csv files for all of the DESEQ2 results for -3 vs -5
m3.m5.all.out <- data.frame(
  geneName = row.names(m3.m5.all),
  baseMean = m3.m5.all[,1],
  log2FoldChange = m3.m5.all[,2],
  lfcSE = m3.m5.all[,3],
  stat = m3.m5.all[,4],
  pvalue = m3.m5.all[,5],
  padj = m3.m5.all[,6],
  stringsAsFactors = FALSE
)
write.csv(m3.m5.all.out, 'deseq_m3m5.csv', row.names=FALSE)