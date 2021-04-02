#this to create box (violin) plots of counts (filtered, normalised for library size)
#compare X expression in kd v u6m2 (side be side)
#use counts_normal_cpm [<- cpm(counts_filtered, normalized.lib.sizes = TRUE)]
#maybe alternative in EdgeR to use raw counts? 
#But DGEListobject$counts just lists original raw counts, library sizes change to "effective" lib size with normalisation
#With "effective" library sizes, then norm. factor goes to 1
##compare autosome expression in kd v u6m2 (side be side)
setwd("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/")
library(ggplot2)
library(tidyverse)
library(dplyr)
#use counts_filtered (DGEList object - from previous- can't save DGEList object as a whole)
counts_cpm <- data.frame(counts_filtered$genes, counts_normal_cpm)
#counts_filtered$count

#OR read in data from counts normalisation, create dataframe 
#counts_normal_cpm <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/counts_normal_cpm.txt", header=TRUE, sep="\t")
#counts_cpm <- data.frame(counts_filtered[ ,1:3], counts_normal_cpm[ ,1:6])

write.table(counts_cpm, quote=FALSE, sep="\t", file = "counts_cpm.txt", row.names=TRUE)

#subset genes
xgenes <- counts_cpm %>% filter(chromosome=="X")
c1genes <- counts_cpm %>% filter(chromosome=="1")
c2genes <- counts_cpm %>% filter(chromosome=="2")
c3genes <- counts_cpm %>% filter(chromosome=="3")
c4genes <- counts_cpm %>% filter(chromosome=="4")
c5genes <- counts_cpm %>% filter(chromosome=="5")
c6genes <- counts_cpm %>% filter(chromosome=="6")
c7genes <- counts_cpm %>% filter(chromosome=="7")
c8genes <- counts_cpm %>% filter(chromosome=="8")
autogenes <- counts_cpm %>% filter(chromosome!="X")
#NB already filtered out all genes not attached to chromosomes

#create dataframe, with kd listed in column
u6x <- data.frame(group = "u6m2", id = "X", value = log2(xgenes$u6m2))
hkx <- data.frame(group = "hnrnpk", id = "X", value = log2(xgenes$hnrnpk))
ckx <- data.frame(group = "ckap4", id = "X", value = log2(xgenes$ckap4))
syx <- data.frame(group = "syncrip", id = "X", value = log2(xgenes$syncrip))
cpx <- data.frame(group = "caprin1", id = "X", value = log2(xgenes$caprin1))
nox <- data.frame(group = "nono", id = "X", value = log2(xgenes$nono))
#rowbind into long list
df_boxplot_x <- rbind(u6x,hkx,ckx,syx,cpx,nox)

#same for autosomes
u6a <- data.frame(group = "u6m2", id = "A", value = log2(autogenes$u6m2))
hka <- data.frame(group = "hnrnpk", id = "A", value = log2(autogenes$hnrnpk))
cka <- data.frame(group = "ckap4", id = "A", value = log2(autogenes$ckap4))
sya <- data.frame(group = "syncrip", id = "A", value = log2(autogenes$syncrip))
cpa <- data.frame(group = "caprin1", id = "A", value = log2(autogenes$caprin1))
noa <- data.frame(group = "nono", id = "A", value = log2(autogenes$nono))
df_boxplot_a <- rbind(u6a,hka,cka,sya,cpa,noa)

#combine both into one df
df_boxplot <-rbind(df_boxplot_x, df_boxplot_a)


# grouped boxplot: NB, this shows median line
plot_box_log2_XandA <- ggplot(df_boxplot, aes(x=group, y=value, fill=id)) + 
  geom_boxplot()+
  xlab("knockdown")+
  ylab("log2_cpm")+
  labs(fill = "")
plot_box_log2_XandA

#violinplot
plot_violin_log2_XandA <- ggplot(df_boxplot, aes(x=group, y=value, fill=id)) + 
  geom_violin()+
  xlab("knockdown")+
  ylab("log2_cpm")+
  labs(fill = "")
plot_violin_log2_XandA

#same but raw counts (cpm) - no good, median counts too low
#range(xgenes$hnrnpk)
#median(xgenes$hnrnpk)
#mean(xgenes$hnrnpk)
##needto ignore outliers - they are distorting the plots to baseline 


#create dataframe, with kd listed in column
u6x2 <- data.frame(group = "u6m2", id = "X", value = (xgenes$u6m2))
hkx2 <- data.frame(group = "hnrnpk", id = "X", value = (xgenes$hnrnpk))
ckx2 <- data.frame(group = "ckap4", id = "X", value = (xgenes$ckap4))
syx2 <- data.frame(group = "syncrip", id = "X", value = (xgenes$syncrip))
cpx2 <- data.frame(group = "caprin1", id = "X", value = (xgenes$caprin1))
nox2 <- data.frame(group = "nono", id = "X", value = (xgenes$nono))
#rowbind into long list
df_boxplot_x2 <- rbind(u6x2,hkx2,ckx2,syx2,cpx2,nox2)

#same for autosomes
u6a2 <- data.frame(group = "u6m2", id = "A", value = (autogenes$u6m2))
hka2 <- data.frame(group = "hnrnpk", id = "A", value = (autogenes$hnrnpk))
cka2 <- data.frame(group = "ckap4", id = "A", value = (autogenes$ckap4))
sya2 <- data.frame(group = "syncrip", id = "A", value = (autogenes$syncrip))
cpa2 <- data.frame(group = "caprin1", id = "A", value = (autogenes$caprin1))
noa2 <- data.frame(group = "nono", id = "A", value = (autogenes$nono))
df_boxplot_a2 <- rbind(u6a2,hka2,cka2,sya2,cpa2,noa2)

#combine both into one df
df_boxplot2 <-rbind(df_boxplot_x2, df_boxplot_a2)


# grouped boxplot: NB, this shows median line
plot_box_XandA <- ggplot(df_boxplot2, aes(x=group, y=value, fill=id)) + 
  geom_boxplot()+
  xlab("knockdown")+
  ylab("cpm")+
  labs(fill = "")
plot_box_XandA

#violinplot
plot_violin_XandA <- ggplot(df_boxplot2, aes(x=group, y=value, fill=id)) + 
  geom_violin()+
  xlab("knockdown")+
  ylab("cpm")+
  labs(fill = "")
plot_violin_XandA

