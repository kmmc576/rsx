setwd("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/")
#read in log2 ratios (kd/U6M2) for x and autosomes from log2 scatterplot script
counts_auto_sort <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/counts_log2_auto_sort.txt", header=TRUE, sep="\t")
counts_log2_x <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/counts_log2_x.txt", header=TRUE, sep="\t")
#create boxplot
library(ggplot2)
library(tidyverse)
library(dplyr)
#create dataframe, with kd listed in column
hkx <- data.frame(group = "hnrnpk", id = "X_chr", value = counts_log2_x$hnrnpk)
ckx <- data.frame(group = "ckap4", id = "X_chr", value = counts_log2_x$ckap4)
syx <- data.frame(group = "syncrip", id = "X_chr", value = counts_log2_x$syncrip)
cpx <- data.frame(group = "caprin1", id = "X_chr", value = counts_log2_x$caprin1)
nox <- data.frame(group = "nono", id = "X_chr", value = counts_log2_x$nono)

df_boxplot_x <- rbind(hkx,ckx,syx,cpx,nox)
#same for autosomes
hka <- data.frame(group = "hnrnpk", id = "Autosomes", value = counts_auto_sort$hnrnpk)
cka <- data.frame(group = "ckap4", id = "Autosomes", value = counts_auto_sort$ckap4)
sya <- data.frame(group = "syncrip", id = "Autosomes", value = counts_auto_sort$syncrip)
cpa <- data.frame(group = "caprin1", id = "Autosomes", value = counts_auto_sort$caprin1)
noa <- data.frame(group = "nono", id = "Autosomes", value = counts_auto_sort$nono)

df_boxplot_a <- rbind(hka,cka,sya,cpa,noa)
#combine both into one df
df_boxplot <-rbind(df_boxplot_x, df_boxplot_a)

# grouped boxplot: NB, this shows median line
plot_XandA_log2 <- ggplot(df_boxplot, aes(x=group, y=value, fill=id)) + 
  geom_boxplot()+
  xlab("knockdown")+
  ylab("log2[kd_cpm/u6m2_cpm]")+
  labs(fill = "")
plot_XandA_log2 
 
#violinplot
plot_XandA_log2_violin <- ggplot(df_boxplot, aes(x=group, y=value, fill=id)) + 
  geom_violin()+
  xlab("knockdown")+
  ylab("log2[kd_cpm/u6m2_cpm]")+
  labs(fill = "")
plot_XandA_log2_violin
