#1) scatter plot log2 ratio (cf u6m2), sorted on TSS for X - combine all 5 kd in one facet grid
#2) scatter plot log2 ratio (cf u6m2), sorted on TSS for each autosome - combine all 5 kd in one facet grid for each autosome
#3) boxplot log2 ratio (cf u6m2) X and A paired, for each kd (all on one graph)
#4) violin plots of (3)
#5) boxplot log2 ratio u6m2 x gene and autosomal gene cpm (counts, not ratio) for each kd
#6) identify 2x, 3x, 5x, 10x fold changes (outliers)

#for comparison of expression of gene in U6M2 cf kd sample, don't need to normalise for gene length (constant)
#plot expression of all X genes in U6M2 cf each kd 
#plot (log2) ratio of expression in U6M2 cf each kd for each X gene. check why log2 rather than natural log?
#do scatterplot of log2 of ratios, map along loci locations on X chromosome - get from BioMart??
library(ggplot2)
library(tidyverse)
library(plotly)

#use counts_filtered (DGEList object)
#OR read in data from counts normalisation, create dataframe 
#counts_filtered <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/counts_filtered_genes.txt", header=TRUE, sep="\t")
#counts_normal <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/counts_normal_cpm.txt", header=TRUE, sep="\t")

#counts_cpm <- data.frame(counts_filtered[ ,1:3], counts_normal[ ,1:6])
# or for DGEList object, use: 
counts_cpm <- data.frame(counts_filtered$genes, counts_normal_cpm)
names(counts_cpm)[13] <- "nono"

write.table(counts_cpm, quote=FALSE, sep="\t", file = "counts_cpm.txt", row.names=TRUE)

#subset genes, then order based on TSS
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

#Sort xgenes, c1, etc  based on TSS (lowest to highest)
#Sort autogenes  based on: 1) chromosome number, then 2) TSS
xgenes <- xgenes[order(xgenes$TSS),]
c1genes <- c1genes[order(c1genes$TSS),]
c2genes <- c2genes[order(c2genes$TSS),]
c3genes <- c3genes[order(c3genes$TSS),]
c4genes <- c4genes[order(c4genes$TSS),]
c5genes <- c5genes[order(c5genes$TSS),]
c6genes <- c6genes[order(c6genes$TSS),]
c7genes <- c7genes[order(c7genes$TSS),]
c8genes <- c8genes[order(c8genes$TSS),]
autogenes <- autogenes[order(autogenes$chromosome, autogenes$TSS),]


#DATAFRAMES
#create a dataframe of log2 expression ratio for each kd (kd/u6m)_x expression
hnrnpk <- c(log2(xgenes$hnrnpk/xgenes$u6m2))
ckap4 <- c(log2(xgenes$ckap4/xgenes$u6m2))
syncrip <- c(log2(xgenes$syncrip/xgenes$u6m2))
caprin1 <- c(log2(xgenes$caprin1/xgenes$u6m2))
nono <- c(log2(xgenes$nono/xgenes$u6m2))

log2_x <- data.frame(xgenes[ ,1:7], hnrnpk, ckap4,syncrip, caprin1, nono)
write.table(log2_x, quote=FALSE, sep="\t", file = "counts_log2_x.txt", row.names=TRUE)

#check range of log2values for plot y axis
rg_hnrnpk <- range(log2_x$hnrnpk,na.rm=FALSE, finite=TRUE)
rg_ckap4 <- range(log2_x$ckap4,na.rm=FALSE, finite=TRUE)
rg_syncrip <- range(log2_x$syncrip,na.rm=FALSE, finite=TRUE)
rg_caprin1 <- range(log2_x$caprin1,na.rm=FALSE, finite=TRUE)
rg_nono <- range(log2_x$nono,na.rm=FALSE, finite=TRUE)

#autosomes
#create a dataframe of log2 expression ratio for each kd (kd/u6m)_auto expression
hnrnpk <- c(log2(autogenes$hnrnpk/autogenes$u6m2))
ckap4 <- c(log2(autogenes$ckap4/autogenes$u6m2))
syncrip <- c(log2(autogenes$syncrip/autogenes$u6m2))
caprin1 <- c(log2(autogenes$caprin1/autogenes$u6m2))
nono <- c(log2(autogenes$nono/autogenes$u6m2))

log2_auto <- data.frame(autogenes[ ,1:7], hnrnpk, ckap4,syncrip, caprin1, nono)
write.table(log2_auto, quote=FALSE, sep="\t", file = "counts_log2_auto.txt", row.names=TRUE)

#check range of log2values for plot y axis
rg_hnrnpk <- range(log2_auto_sort$hnrnpk,na.rm=FALSE, finite=TRUE)
rg_ckap4 <- range(log2_auto_sort$ckap4,na.rm=FALSE, finite=TRUE)
rg_syncrip <- range(log2_auto_sort$syncrip,na.rm=FALSE, finite=TRUE)
rg_caprin1 <- range(log2_auto_sort$caprin1,na.rm=FALSE, finite=TRUE)
rg_nono <- range(log2_auto_sort$nono,na.rm=FALSE, finite=TRUE)


#PLOTS
#X plots
#plot log2 ratios as scatterplots, ordered by TSS

plot_log2_x_hnrnpk <- ggplot(log2_x, aes(TSS, hnrnpk))+
  geom_point(aes(color=hnrnpk), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("") 
plot_log2_x_hnrnpk
#plotly_log2_x_hnrnpk <- ggplotly(plot_log2_x_hnrnpk)
#plotly_log2_x_hnrnpk

plot_log2_x_ckap4 <- ggplot(log2_x, aes(TSS, ckap4))+
  geom_point(aes(color=ckap4), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
plot_log2_x_ckap4
#plotly_log2_x_ckap4 <- ggplotly(plot_log2_x_ckap4)
#plotly_log2_x_ckap4

plot_log2_x_syncrip <- ggplot(log2_x, aes(TSS, syncrip))+
  geom_point(aes(color=syncrip), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
plot_log2_x_syncrip
#plotly_log2_x_syncrip <- ggplotly(plot_log2_x_syncrip)
#plotly_log2_x_syncrip

plot_log2_x_caprin1 <- ggplot(log2_x, aes(TSS, caprin1))+
  geom_point(aes(color=caprin1), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
plot_log2_x_caprin1
#plotly_log2_x_caprin1 <- ggplotly(plot_log2_x_caprin1)
#plotly_log2_x_caprin1

plot_log2_x_nono <- ggplot(log2_x, aes(TSS, nono))+
  geom_point(aes(color=nono), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0 ,limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
plot_log2_x_nono
#plotly_log2_x_nono <- ggplotly(plot_log2_x_nono)
#plotly_log2_x_nono

#combine all 5 xgene plots using facet
#make into data into long dataframe with kd at a group
#create dataframe, with kd listed in column
hkx <- data.frame(group = "hnrnpk", id = "X", TSS = xgenes$TSS, value = log2(xgenes$hnrnpk/xgenes$u6m2))
ckx <- data.frame(group = "ckap4", id = "X", TSS = xgenes$TSS, value = log2(xgenes$ckap4/xgenes$u6m2))
syx <- data.frame(group = "syncrip", id = "X", TSS = xgenes$TSS, value = log2(xgenes$syncrip/xgenes$u6m2))
cpx <- data.frame(group = "caprin1", id = "X", TSS = xgenes$TSS, value = log2(xgenes$caprin1/xgenes$u6m2))
nox <- data.frame(group = "nono", id = "X", TSS = xgenes$TSS, value = log2(xgenes$nono/xgenes$u6m2))
#rowbind into long list
df_facet_x <- rbind(hkx,ckx,syx,cpx,nox)

#use factor, levels to reorder facet graphs (Default is alphabetical) 
df_facet_x$group <- factor(df_facet_x$group, levels = c("hnrnpk", "ckap4", "syncrip", "caprin1", "nono"))
plot_log2_x_all_facet <- ggplot(df_facet_x, aes(TSS, value))+
  geom_point(aes(color=value), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0 ,limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("TSS_X_chromosome")+
  ylab("log2[cpm_kd/cpm_U6M2")+
  facet_grid(group ~.)
plot_log2_x_all_facet


#AUTOSOME plots
#plot log2 ratios as scatterplots
#list programmed colours using grDevices::colors()
#plot_log2_auto_all <- ggplot()+
#  geom_point(aes(x=log2_auto_sort$ensemble_id, y=log2_auto_sort$hnrnpk, color="hnrnpk"))+
#  geom_point(aes(x=log2_auto_sort$ensemble_id, y=log2_auto_sort$ckap4, color="ckap4"))+
#  geom_point(aes(x=log2_auto_sort$ensemble_id, y=log2_auto_sort$syncrip, color="syncrip"))+
#  geom_point(aes(x=log2_auto_sort$ensemble_id, y=log2_auto_sort$caprin1, color="caprin1"))+
#  geom_point(aes(x=log2_auto_sort$ensemble_id, y=log2_auto_sort$nono, color="nono"))

#plot_log2_auto_all
#arrange by chromosome
#to add vertical line use  geom_vline(xintercept = 928262)
#add x label with number of genes
#this plots all point in vertical line - use violin plot instead??
plot_log2_auto_hnrnpk_vertical <- ggplot(log2_auto, aes(x=chromosome, y=hnrnpk))+
  geom_point(aes(color=hnrnpk), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
plot_log2_auto_hnrnpk_vertical
#ggplotly(plot_log2_auto_hnrnpk)

plot_log2_auto_ckap4_vertical <- ggplot(log2_auto, aes(x=chromosome, y=ckap4))+
  geom_point(aes(color=ckap4), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
  plot_log2_auto_ckap4_vertical

plot_log2_auto_syncrip_vertical <- ggplot(log2_auto, aes(x=chromosome, y=syncrip))+
  geom_point(aes(color=syncrip), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
  plot_log2_auto_syncrip_vertical

plot_log2_auto_caprin1_vertical <- ggplot(log2_auto, aes(x=chromosome, y=caprin1))+
  geom_point(aes(color=caprin1), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
  plot_log2_auto_caprin1_vertical

plot_log2_auto_nono_vertical <- ggplot(log2_auto, aes(x=chromosome, y=nono))+
  geom_point(aes(color=nono), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")
  plot_log2_auto_nono_vertical

#arrange horizontal, sorted on TSS, but grouped by chromosome
  #Hnrnpk - chromosomes 1 to 8
plot_log2_auto_hnrnpk_grouped1 <- ggplot(c1genes, aes(TSS, y=log2(hnrnpk/u6m2)))+
    geom_point(aes(color=log2(hnrnpk/u6m2)), shape=18)+
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
    ylim(-4.2,4.2)+
    xlab("Hnrnpk chr1")+
    ylab("")+
   theme(legend.position = "none")
  plot_log2_auto_hnrnpk_grouped1  
  
plot_log2_auto_hnrnpk_grouped2 <- ggplot(c2genes, aes(TSS, y=log2(hnrnpk/u6m2)))+
    geom_point(aes(color=log2(hnrnpk/u6m2)), shape=18)+
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
    ylim(-4.2,4.2)+
    xlab("chr2")+
    ylab("")+
  theme(legend.position = "none", axis.text.y = element_blank())
plot_log2_auto_hnrnpk_grouped2  


plot_log2_auto_hnrnpk_grouped3 <- ggplot(c3genes, aes(TSS, y=log2(hnrnpk/u6m2)))+
  geom_point(aes(color=log2(hnrnpk/u6m2)), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("chr3")+
  ylab("")+
  theme(legend.position = "none", axis.text.y = element_blank())
plot_log2_auto_hnrnpk_grouped3

plot_log2_auto_hnrnpk_grouped4 <- ggplot(c4genes, aes(TSS, y=log2(hnrnpk/u6m2)))+
  geom_point(aes(color=log2(hnrnpk/u6m2)), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("chr4")+
  ylab("")+
  theme(legend.position = "none", axis.text.y = element_blank())
plot_log2_auto_hnrnpk_grouped4

plot_log2_auto_hnrnpk_grouped5 <- ggplot(c5genes, aes(TSS, y=log2(hnrnpk/u6m2)))+
  geom_point(aes(color=log2(hnrnpk/u6m2)), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("chr5")+
  ylab("")+
  theme(legend.position = "none", axis.text.y = element_blank())
plot_log2_auto_hnrnpk_grouped5

plot_log2_auto_hnrnpk_grouped6 <- ggplot(c6genes, aes(TSS, y=log2(hnrnpk/u6m2)))+
  geom_point(aes(color=log2(hnrnpk/u6m2)), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("chr6")+
  ylab("")+
  theme(legend.position = "none", axis.text.y = element_blank())
plot_log2_auto_hnrnpk_grouped6

plot_log2_auto_hnrnpk_grouped7 <- ggplot(c7genes, aes(TSS, y=log2(hnrnpk/u6m2)))+
  geom_point(aes(color=log2(hnrnpk/u6m2)), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("chr7")+
  ylab("")+
  theme(legend.position = "none", axis.text.y = element_blank())
plot_log2_auto_hnrnpk_grouped7

plot_log2_auto_hnrnpk_grouped8 <- ggplot(c8genes, aes(TSS, y=log2(hnrnpk/u6m2)))+
  geom_point(aes(color=log2(hnrnpk/u6m2)), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("chr8")+
  ylab("")+
  theme(legend.position = "none", axis.text.y = element_blank())
plot_log2_auto_hnrnpk_grouped8 

#combine 8 plots in single
#install.packages("ggpubr")
library("ggpubr")
library(gridExtra)
plot_log2_auto_hnrnpk_grouped_all <- grid.arrange(plot_log2_auto_hnrnpk_grouped1, plot_log2_auto_hnrnpk_grouped2, plot_log2_auto_hnrnpk_grouped3, plot_log2_auto_hnrnpk_grouped4, plot_log2_auto_hnrnpk_grouped5, plot_log2_auto_hnrnpk_grouped6, plot_log2_auto_hnrnpk_grouped7,  plot_log2_auto_hnrnpk_grouped8, nrow = 1)
plot_log2_auto_hnrnpk_grouped_all
#need to add common legend, title
#need to scale all x axes the same

#Facet_grid, based on all autogenes
plot_log2_auto_hnrnpk <- ggplot(log2_auto, aes(x=TSS, y=hnrnpk, group=chromosome))+
  geom_point(aes(color=hnrnpk), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_log2_auto_hnrnpk

plot_log2_auto_ckap4 <- ggplot(log2_auto, aes(x=TSS, y=ckap4, group=chromosome))+
  geom_point(aes(color=ckap4), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_log2_auto_ckap4

plot_log2_auto_syncrip <- ggplot(log2_auto, aes(x=TSS, y=syncrip, group=chromosome))+
  geom_point(aes(color=syncrip), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_log2_auto_syncrip

plot_log2_auto_caprin1 <- ggplot(log2_auto, aes(x=TSS, y=caprin1, group=chromosome))+
  geom_point(aes(color=caprin1), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_log2_auto_caprin1

plot_log2_auto_nono <- ggplot(log2_auto, aes(x=TSS, y=nono, group=chromosome))+
  geom_point(aes(color=nono), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_log2_auto_nono

