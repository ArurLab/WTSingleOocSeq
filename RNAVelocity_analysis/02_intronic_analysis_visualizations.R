### Intronic gene count plots: spliced vs. unspliced

packages = c("BiocManager","ggplot2")
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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load files
filt.unspliced.m3.m1 <- read.csv("tsv_matrices/matrix_filtered_unspliced_m3-m1.tsv", sep='\t')
filt.unspliced.m5.m1 <- read.csv("tsv_matrices/matrix_filtered_unspliced_m5-m1.tsv", sep='\t')
filt.spliced.m3.m1 <- read.csv("tsv_matrices/matrix_filtered_spliced_m3-m1.tsv", sep='\t')
filt.spliced.m5.m1 <- read.csv("tsv_matrices/matrix_filtered_spliced_m5-m1.tsv", sep='\t')

# Use deseq2 to normalize counts
group1="m1"
group2="m3"

#this DESeq object is a combination of  
# 1) a count-matrix
#             treated1 treated2 treated3 untreated1 untreated2 untreated3
# FBgn0000003        0        0        1          0          0          0
# FBgn0000008      140       88       70         92        161         76


#condition
#type
#treated1	treated	single-read
#treated2	treated	paired-end
#treated3	treated	paired-end
#untreated1	untreated	single-read
#untreated2	untreated	single-read
#untreated3	untreated	paired-end
#untreated4	untreated

# 2) a matrix which assigns each of the headers in the count matrix to a group

m3m1 <- filt.spliced.m3.m1
rownames(m3m1) <- m3m1$X
m5m1 <- filt.spliced.m5.m1
rownames(m5m1) <- m5m1$X

all.counts <- data.frame(m3m1[,c(6,7,8,9)], m3m1[,c(2,3,4,5)], m5m1[,c(2,3,4,5)])

Oocyte_ID <- colnames(all.counts)
condition<-c('m1','m1','m1','m1','m3','m3','m3','m3','m5','m5','m5','m5')

coldata<-data.frame(condition)
rownames(coldata) <- Oocyte_ID

dds<-DESeqDataSetFromMatrix(countData=all.counts, colData=coldata, design= ~ condition)

dds<-DESeq(dds)

spliced.cnts<-counts(dds, normalize=TRUE)

m3m1 <- read.csv("tsv_matrices/matrix_filtered_unspliced_m3-m1.tsv", sep='\t')
rownames(m3m1) <- m3m1$X
m5m1 <- read.csv("tsv_matrices/matrix_filtered_unspliced_m5-m1.tsv", sep='\t')
rownames(m5m1) <- m5m1$X

all.counts <- data.frame(m3m1[,c(6,7,8,9)], m3m1[,c(2,3,4,5)], m5m1[,c(2,3,4,5)])

Oocyte_ID <- colnames(all.counts)
condition<-c('m1','m1','m1','m1','m3','m3','m3','m3','m5','m5','m5','m5')

coldata<-data.frame(condition)
rownames(coldata) <- Oocyte_ID

dds<-DESeqDataSetFromMatrix(countData=all.counts, colData=coldata, design= ~ condition)

dds<-DESeq(dds)

unspliced.cnts<-counts(dds, normalize=TRUE)

#### Create dataframes with the normalized counts for each position ####
filt.m1.counts.unspliced <- data.frame(rownames(unspliced.cnts),unspliced.cnts[,c(1,2,3,4)])
filt.m3.counts.unspliced <- data.frame(rownames(unspliced.cnts),unspliced.cnts[,c(5,6,7,8)])
filt.m5.counts.unspliced <- data.frame(rownames(unspliced.cnts),unspliced.cnts[,c(9,10,11,12)])
filt.m1.counts.spliced <- data.frame(rownames(spliced.cnts),spliced.cnts[,c(1,2,3,4)])
filt.m3.counts.spliced <- data.frame(rownames(spliced.cnts),spliced.cnts[,c(5,6,7,8)])
filt.m5.counts.spliced <- data.frame(rownames(spliced.cnts),spliced.cnts[,c(9,10,11,12)])
names(filt.m1.counts.unspliced)[1] = 'X'
names(filt.m3.counts.unspliced)[1] = 'X'
names(filt.m5.counts.unspliced)[1] = 'X'
names(filt.m1.counts.spliced)[1] = 'X'
names(filt.m3.counts.spliced)[1] = 'X'
names(filt.m5.counts.spliced)[1] = 'X'

intronic_cutoff <- 10
filt.m1.unspliced.list <- filt.m1.counts.unspliced$X
filt.m3.unspliced.list <- filt.m3.counts.unspliced$X
filt.m5.unspliced.list <- filt.m5.counts.unspliced$X
for (i in 2:5) {
  filt.m1.unspliced.list <- intersect(filt.m1.unspliced.list, filt.m1.counts.unspliced[filt.m1.counts.unspliced[i] > intronic_cutoff,1])
  filt.m3.unspliced.list <- intersect(filt.m3.unspliced.list, filt.m3.counts.unspliced[filt.m3.counts.unspliced[i] > intronic_cutoff,1])
  filt.m5.unspliced.list <- intersect(filt.m5.unspliced.list, filt.m5.counts.unspliced[filt.m5.counts.unspliced[i] > intronic_cutoff,1])
}

spliced_cutoff <- 10
filt.m1.spliced.list <- filt.m1.counts.spliced[rowMeans(filt.m1.counts.spliced[c(2,3,4,5)]) > spliced_cutoff,1]
filt.m3.spliced.list <- filt.m3.counts.spliced[rowMeans(filt.m3.counts.spliced[c(2,3,4,5)]) > spliced_cutoff,1]
filt.m5.spliced.list <- filt.m5.counts.spliced[rowMeans(filt.m5.counts.spliced[c(2,3,4,5)]) > spliced_cutoff,1]

filt.filtered.gene.list.m1 <- intersect(filt.m1.unspliced.list,filt.m1.spliced.list)
filt.filtered.gene.list.m3 <- intersect(filt.m3.unspliced.list,filt.m3.spliced.list)
filt.filtered.gene.list.m5 <- intersect(filt.m5.unspliced.list,filt.m5.spliced.list)

filt.filtered.counts.m1.unspliced <- filt.m1.counts.unspliced[filt.m1.counts.unspliced$X %in% filt.filtered.gene.list.m1,]
filt.filtered.counts.m3.unspliced <- filt.m3.counts.unspliced[filt.m3.counts.unspliced$X %in% filt.filtered.gene.list.m3,]
filt.filtered.counts.m5.unspliced <- filt.m5.counts.unspliced[filt.m5.counts.unspliced$X %in% filt.filtered.gene.list.m5,]
filt.filtered.counts.m1.spliced <- filt.m1.counts.spliced[filt.m1.counts.spliced$X %in% filt.filtered.gene.list.m1,]
filt.filtered.counts.m3.spliced <- filt.m3.counts.spliced[filt.m3.counts.spliced$X %in% filt.filtered.gene.list.m3,]
filt.filtered.counts.m5.spliced <- filt.m5.counts.spliced[filt.m5.counts.spliced$X %in% filt.filtered.gene.list.m5,]

filt.m1.counts.spliced.avg <- data.frame(filt.filtered.counts.m1.spliced[1],rep('spliced', nrow(filt.filtered.counts.m1.spliced)),rowMeans(filt.filtered.counts.m1.spliced[c(2,3,4,5)]))
filt.m3.counts.spliced.avg <- data.frame(filt.filtered.counts.m3.spliced[1],rep('spliced', nrow(filt.filtered.counts.m3.spliced)),rowMeans(filt.filtered.counts.m3.spliced[c(2,3,4,5)]))
filt.m5.counts.spliced.avg <- data.frame(filt.filtered.counts.m5.spliced[1],rep('spliced', nrow(filt.filtered.counts.m5.spliced)),rowMeans(filt.filtered.counts.m5.spliced[c(2,3,4,5)]))
filt.m1.counts.unspliced.avg <- data.frame(filt.filtered.counts.m1.unspliced[1],rep('unspliced', nrow(filt.filtered.counts.m1.spliced)),rowMeans(filt.filtered.counts.m1.unspliced[c(2,3,4,5)]))
filt.m3.counts.unspliced.avg <- data.frame(filt.filtered.counts.m3.unspliced[1],rep('unspliced', nrow(filt.filtered.counts.m3.unspliced)),rowMeans(filt.filtered.counts.m3.unspliced[c(2,3,4,5)]))
filt.m5.counts.unspliced.avg <- data.frame(filt.filtered.counts.m5.unspliced[1],rep('unspliced', nrow(filt.filtered.counts.m5.unspliced)),rowMeans(filt.filtered.counts.m5.unspliced[c(2,3,4,5)]))

names(filt.m1.counts.spliced.avg) <- c('gene','type','count')
names(filt.m3.counts.spliced.avg) <- c('gene','type','count')
names(filt.m5.counts.spliced.avg) <- c('gene','type','count')
names(filt.m1.counts.unspliced.avg) <- c('gene','type','count')
names(filt.m3.counts.unspliced.avg) <- c('gene','type','count')
names(filt.m5.counts.unspliced.avg) <- c('gene','type','count')
filt.m1.counts.merged.avg <- rbind(filt.m1.counts.spliced.avg, filt.m1.counts.unspliced.avg)
filt.m3.counts.merged.avg <- rbind(filt.m3.counts.spliced.avg, filt.m3.counts.unspliced.avg)
filt.m5.counts.merged.avg <- rbind(filt.m5.counts.spliced.avg, filt.m5.counts.unspliced.avg)

title <- '-1 Oocyte Intronic Gene Counts'
filt.m1.intronic.genes.counts.plot <- ggplot(data=filt.m1.counts.merged.avg, aes(x=gene, y=count, fill=type)) + geom_bar(position='dodge', stat='identity')
filt.m1.intronic.genes.counts.plot <- filt.m1.intronic.genes.counts.plot + guides(x=guide_axis(angle=60)) + theme_bw() + labs(title=title)
filt.m1.intronic.genes.counts.plot <- filt.m1.intronic.genes.counts.plot + scale_fill_grey(start=0.5, end=0.8)
filt.m1.intronic.genes.counts.plot

title <- '-3 Oocyte Intronic Gene Counts'
filt.m3.intronic.genes.counts.plot <- ggplot(data=filt.m3.counts.merged.avg, aes(x=gene, y=count, fill=type)) + geom_bar(position='dodge', stat='identity')
filt.m3.intronic.genes.counts.plot <- filt.m3.intronic.genes.counts.plot + guides(x=guide_axis(angle=60)) + theme_bw() + labs(title=title)
filt.m3.intronic.genes.counts.plot <- filt.m3.intronic.genes.counts.plot + scale_fill_grey(start=0.5, end=0.8)
filt.m3.intronic.genes.counts.plot

title <- '-5 Oocyte Intronic Gene Counts'
filt.m5.intronic.genes.counts.plot <- ggplot(data=filt.m5.counts.merged.avg, aes(x=gene, y=count, fill=type)) + geom_bar(position='dodge', stat='identity')
filt.m5.intronic.genes.counts.plot <- filt.m5.intronic.genes.counts.plot + guides(x=guide_axis(angle=60))
filt.m5.intronic.genes.counts.plot <- filt.m5.intronic.genes.counts.plot + scale_fill_grey(start=0.5, end=0.8) + theme_bw() + labs(title=title)
filt.m5.intronic.genes.counts.plot

## Create a plot showing exonic counts in M1 vs M3 for m3 intronic genes

# get intronic genes in common
m1.m3.common.genes <- intersect(filt.filtered.gene.list.m1,filt.filtered.gene.list.m3)
m1.m3.common.gene.counts <- data.frame(filt.m1.counts.spliced[filt.m1.counts.spliced$X %in% filt.filtered.gene.list.m3,],filt.m3.counts.spliced[filt.m3.counts.spliced$X %in% filt.filtered.gene.list.m3,c(2,3,4,5)])
names(m1.m3.common.gene.counts)[1] = 'gene'
m1.m3.common.gene.counts$m1Sum <- rowSums(m1.m3.common.gene.counts[c(2,3,4,5)])
m1.m3.common.gene.counts$m1Avg <- rowMeans(m1.m3.common.gene.counts[c(2,3,4,5)])
m1.m3.common.gene.counts$m3Sum <- rowSums(m1.m3.common.gene.counts[c(6,7,8,9)])
m1.m3.common.gene.counts$m3Avg <- rowMeans(m1.m3.common.gene.counts[c(6,7,8,9)])

#### m1 vs m3 plot
filt.m1.m3.int.genes.sum.plot <- ggplot(data=m1.m3.common.gene.counts, aes(x=m3Sum, y=m1Sum)) + geom_point()
filt.m1.m3.int.genes.sum.plot <- filt.m1.m3.int.genes.sum.plot + theme_bw() 
filt.m1.m3.int.genes.sum.plot <- filt.m1.m3.int.genes.sum.plot + coord_cartesian(xlim=c(110,2390), ylim=c(110,2390)) 

title <- 'm1 vs m3 for m3 intronic genes with .5 fold change'
filt.m1.m3.int.genes.sum.plot.1.5 <- filt.m1.m3.int.genes.sum.plot + labs(title=title) + geom_abline(slope=1, intercept=0) + geom_abline(slope=1.5, intercept=0, linetype='dashed') + geom_abline(slope=0.6666, intercept=0, linetype='dashed')
filt.m1.m3.int.genes.sum.plot.1.5

## Now for m1 vs m5
m1.m5.common.genes <- intersect(filt.filtered.gene.list.m1,filt.filtered.gene.list.m5)
m1.m5.common.gene.counts <- data.frame(filt.m1.counts.spliced[filt.m1.counts.spliced$X %in% filt.filtered.gene.list.m5,],filt.m5.counts.spliced[filt.m5.counts.spliced$X %in% filt.filtered.gene.list.m5,c(2,3,4,5)])
names(m1.m5.common.gene.counts)[1] = 'gene'
m1.m5.common.gene.counts$m1Sum <- rowSums(m1.m5.common.gene.counts[c(2,3,4,5)])
m1.m5.common.gene.counts$m1Avg <- rowMeans(m1.m5.common.gene.counts[c(2,3,4,5)])
m1.m5.common.gene.counts$m5Sum <- rowSums(m1.m5.common.gene.counts[c(6,7,8,9)])
m1.m5.common.gene.counts$m5Avg <- rowMeans(m1.m5.common.gene.counts[c(6,7,8,9)])

#### m1 vs m5 plots
filt.m1.m5.int.genes.sum.plot <- ggplot(data=m1.m5.common.gene.counts, aes(x=m5Sum, y=m1Sum)) + geom_point()
filt.m1.m5.int.genes.sum.plot <- filt.m1.m5.int.genes.sum.plot + theme_bw()
filt.m1.m5.int.genes.sum.plot <- filt.m1.m5.int.genes.sum.plot + coord_cartesian(xlim=c(110,2390), ylim=c(110,2390)) 

title <- 'm1 vs m5 for m5 intronic genes with .5 fold change'
filt.m1.m5.int.genes.sum.plot.1.5 <- filt.m1.m5.int.genes.sum.plot + labs(title=title) + geom_abline(slope=1, intercept=0) + geom_abline(slope=1.5, intercept=0, linetype='dashed') + geom_abline(slope=0.6666, intercept=0, linetype='dashed')
filt.m1.m5.int.genes.sum.plot.1.5

