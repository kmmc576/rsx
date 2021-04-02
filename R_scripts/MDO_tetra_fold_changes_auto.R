#identify 2x, 4x, 8x fold changes (outliers)
#read in log2 ratio spreadsheets (from log2 ratio scatterplots)
#Outliers: log2 fold change (ratio kd/u6m2)
#log2(1)=0, log2(2)=1, log2(4)=2, log2(8)=3, log2(16)=4
log2_auto <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/counts_log2_auto.txt",header=TRUE, sep="\t")
log2_hka <- log2_auto[,1:8]
fc_pos2_hka <- log2_hka %>% filter(hnrnpk>=1)
fc_pos4_hka <- log2_hka %>% filter(hnrnpk>=2)
fc_pos8_hka <- log2_hka %>% filter(hnrnpk>=3)
fc_neg2_hka <- log2_hka %>% filter(hnrnpk<=-1)   
fc_neg4_hka <- log2_hka %>% filter(hnrnpk<=-2) 
fc_neg8_hka <- log2_hka %>% filter(hnrnpk<=-3)
fc_hka_all <- rbind(fc_pos2_hka, fc_neg2_hka, fc_neg4_hka, fc_neg8_hka)
#remove duplicates
fc_hka_all <- fc_hka_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

log2_cka <- subset(log2_auto, select = -c(hnrnpk, syncrip, caprin1, nono))
fc_pos2_cka <- log2_cka %>% filter(ckap4>=1)
fc_pos4_cka <- log2_cka %>% filter(ckap4>=2)
fc_pos8_cka <- log2_cka %>% filter(ckap4>=3)
fc_neg2_cka <- log2_cka %>% filter(ckap4<=-1)   
fc_neg4_cka <- log2_cka %>% filter(ckap4<=-2) 
fc_neg8_cka <- log2_cka %>% filter(ckap4<=-3)
fc_cka_all <- rbind(fc_pos2_cka, fc_pos4_cka, fc_pos8_cka, fc_neg2_cka, fc_neg4_cka, fc_neg8_cka)
#remove duplicates
fc_cka_all <- fc_cka_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

log2_sya <- subset(log2_auto, select = -c(hnrnpk, ckap4, caprin1, nono))
fc_pos2_sya <- log2_sya %>% filter(syncrip>=1)
fc_pos4_sya <- log2_sya %>% filter(syncrip>=2)
fc_pos8_sya <- log2_sya %>% filter(syncrip>=3)
fc_neg2_sya <- log2_sya %>% filter(syncrip<=-1)   
fc_neg4_sya <- log2_sya %>% filter(syncrip<=-2) 
fc_neg8_sya <- log2_sya %>% filter(syncrip<=-3)
fc_sya_all <- rbind(fc_pos2_sya, fc_pos4_sya, fc_pos8_sya, fc_neg2_sya, fc_neg4_sya, fc_neg8_sya)
#remove duplicates
fc_sya_all <- fc_sya_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

log2_caa <- subset(log2_auto, select = -c(hnrnpk, ckap4, syncrip, nono))
fc_pos2_caa <- log2_caa %>% filter(caprin1>=1)
fc_pos4_caa <- log2_caa %>% filter(caprin1>=2)
fc_pos8_caa <- log2_caa %>% filter(caprin1>=3)
fc_neg2_caa <- log2_caa %>% filter(caprin1<=-1)   
fc_neg4_caa <- log2_caa %>% filter(caprin1<=-2) 
fc_neg8_caa <- log2_caa %>% filter(caprin1<=-3)
fc_caa_all <- rbind(fc_pos2_caa, fc_pos4_caa, fc_pos8_caa, fc_neg2_caa, fc_neg4_caa, fc_neg8_caa)
#remove duplicates
fc_caa_all <- fc_caa_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

log2_noa <- subset(log2_auto, select = -c(hnrnpk, ckap4, syncrip, caprin1))
fc_pos2_noa <- log2_noa %>% filter(nono>=1)
fc_pos4_noa <- log2_noa %>% filter(nono>=2)
fc_pos8_noa <- log2_noa %>% filter(nono>=3)
fc_neg2_noa <- log2_noa %>% filter(nono<=-1)   
fc_neg4_noa <- log2_noa %>% filter(nono<=-2) 
fc_neg8_noa <- log2_noa %>% filter(nono<=-3)
fc_noa_all <- rbind(fc_pos2_noa, fc_pos4_noa, fc_pos8_noa, fc_neg2_noa, fc_neg4_noa, fc_neg8_noa)
#remove duplicates
fc_noa_all <- fc_noa_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

#combine all fc (xgene) datatables into single, with id=kd
fc_hka <- data.frame(ensemble_id = fc_hka_all$ensemble_id, 
                     gene = fc_hka_all$gene, 
                     chromosome = fc_hka_all$chromosome,
                     length = fc_hka_all$length,
                     TSS = fc_hka_all$TSS,
                     kd_group = "hnrnpk",
                     log2_ratio = fc_hka_all$hnrnpk)

fc_cka <- data.frame(ensemble_id = fc_cka_all$ensemble_id, 
                     gene = fc_cka_all$gene, 
                     chromosome = fc_cka_all$chromosome,
                     length = fc_cka_all$length,
                     TSS = fc_cka_all$TSS,
                     kd_group = "ckap4",
                     log2_ratio = fc_cka_all$ckap4)

fc_sya <- data.frame(ensemble_id = fc_sya_all$ensemble_id, 
                     gene = fc_sya_all$gene, 
                     chromosome = fc_sya_all$chromosome,
                     length = fc_sya_all$length,
                     TSS = fc_sya_all$TSS,
                     kd_group = "syncrip",
                     log2_ratio = fc_sya_all$syncrip)

fc_caa <- data.frame(ensemble_id = fc_caa_all$ensemble_id, 
                     gene = fc_caa_all$gene, 
                     chromosome = fc_caa_all$chromosome,
                     length = fc_caa_all$length,
                     TSS = fc_caa_all$TSS,
                     kd_group = "caprin1",
                     log2_ratio = fc_caa_all$caprin1)

fc_noa <- data.frame(ensemble_id = fc_noa_all$ensemble_id, 
                     gene = fc_noa_all$gene, 
                     chromosome = fc_noa_all$chromosome,
                     length = fc_noa_all$length,
                     TSS = fc_noa_all$TSS,
                     kd_group = "nono",
                     log2_ratio = fc_noa_all$nono)

#rowbind into long list
fc_all_auto <- rbind(fc_hka, fc_cka, fc_sya, fc_caa, fc_noa)

#scatterplot all together
plot_fc_all_auto <- ggplot(fc_all_auto, aes(TSS, log2_ratio))+
  geom_point(aes(color=kd_group), shape=18)+
  ylim(-4.2,4.2)+
  xlab("") 
plot_fc_all_auto

#Facet_grid, based on all autogenes
plot_fc_auto_hnrnpk <- ggplot(fc_hka_all, aes(x=TSS, y=hnrnpk, group=chromosome))+
  geom_point(aes(color=hnrnpk), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_fc_auto_hnrnpk

plot_fc_auto_ckap4 <- ggplot(fc_cka_all, aes(x=TSS, y=ckap4, group=chromosome))+
  geom_point(aes(color=ckap4), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_fc_auto_ckap4

plot_fc_auto_syncrip <- ggplot(fc_sya_all, aes(x=TSS, y=syncrip, group=chromosome))+
  geom_point(aes(color=syncrip), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_fc_auto_syncrip

plot_fc_auto_caprin1 <- ggplot(fc_caa_all, aes(x=TSS, y=caprin1, group=chromosome))+
  geom_point(aes(color=caprin1), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_fc_auto_caprin1

plot_fc_auto_nono <- ggplot(fc_noa_all, aes(x=TSS, y=nono, group=chromosome))+
  geom_point(aes(color=nono), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("")+
  facet_grid(chromosome ~.)
plot_fc_auto_nono


#export fc_all_auto
setwd("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/")
write.table(fc_all_auto, quote=FALSE, sep="\t", file = "fold_change_all_auto.txt", row.names=FALSE, col.names = TRUE)

#export gene list for GO analysis
write.table(fc_all_auto$ensemble_id, quote=FALSE, sep="\t", file = "fold_change_genes_auto.txt", row.names=FALSE, col.names = FALSE)

