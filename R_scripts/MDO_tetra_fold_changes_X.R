#identify 2x, 4x, 8x fold changes (outliers)
#read in log2 ratio spreadsheets (from log2 ratio scatterplots)
#Outliers: log2 fold change (ratio kd/u6m2)
#log2(1)=0, log2(2)=1, log2(4)=2, log2(8)=3, log2(16)=4
log2_x <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/counts_log2_x.txt",header=TRUE, sep="\t")
log2_hkx <- log2_x[,1:8]
fc_pos2_hkx <- log2_hkx %>% filter(hnrnpk>=1)
#highest log2 ratio =0.97
#nb noval lncRNA ENSMODG00000050503 is upregulated: log2 ratio 0.97
fc_neg2_hkx <- log2_hkx %>% filter(hnrnpk<=-1)   
fc_neg4_hkx <- log2_hkx %>% filter(hnrnpk<=-2) 
fc_neg8_hkx <- log2_hkx %>% filter(hnrnpk<=-3)
fc_hkx_all <- rbind(fc_pos2_hkx, fc_neg2_hkx, fc_neg4_hkx, fc_neg8_hkx)
#remove duplicates
fc_hkx_all <- fc_hkx_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

log2_ckx <- subset(log2_x, select = -c(hnrnpk, syncrip, caprin1, nono))
fc_pos2_ckx <- log2_ckx %>% filter(ckap4>=1)
#highest log2 ratio =3.39
fc_pos4_ckx <- log2_ckx %>% filter(ckap4>=2)
fc_pos8_ckx <- log2_ckx %>% filter(ckap4>=3)
fc_neg2_ckx <- log2_ckx %>% filter(ckap4<=-1)   
fc_neg4_ckx <- log2_ckx %>% filter(ckap4<=-2) 
fc_neg8_ckx <- log2_ckx %>% filter(ckap4<=-3)
fc_ckx_all <- rbind(fc_pos2_ckx, fc_pos4_ckx, fc_pos8_ckx, fc_neg2_ckx, fc_neg4_ckx, fc_neg8_ckx)
#remove duplicates
fc_ckx_all <- fc_ckx_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

log2_syx <- subset(log2_x, select = -c(hnrnpk, ckap4, caprin1, nono))
fc_pos2_syx <- log2_syx %>% filter(syncrip>=1)
#highest log2 ratio =3.05
fc_pos4_syx <- log2_syx %>% filter(syncrip>=2)
fc_pos8_syx <- log2_syx %>% filter(syncrip>=3)
fc_neg2_syx <- log2_syx %>% filter(syncrip<=-1)   
fc_neg4_syx <- log2_syx %>% filter(syncrip<=-2) 
fc_neg8_syx <- log2_syx %>% filter(syncrip<=-3)
fc_syx_all <- rbind(fc_pos2_syx, fc_pos4_syx, fc_pos8_syx, fc_neg2_syx, fc_neg4_syx, fc_neg8_syx)
#remove duplicates
fc_syx_all <- fc_syx_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

log2_cax <- subset(log2_x, select = -c(hnrnpk, ckap4, syncrip, nono))
fc_pos2_cax <- log2_cax %>% filter(caprin1>=1)
#highest log2 ratio =2.76
fc_pos4_cax <- log2_cax %>% filter(caprin1>=2)
fc_pos8_cax <- log2_cax %>% filter(caprin1>=3)
fc_neg2_cax <- log2_cax %>% filter(caprin1<=-1)   
fc_neg4_cax <- log2_cax %>% filter(caprin1<=-2) 
fc_neg8_cax <- log2_cax %>% filter(caprin1<=-3)
fc_cax_all <- rbind(fc_pos2_cax, fc_pos4_cax, fc_pos8_cax, fc_neg2_cax, fc_neg4_cax, fc_neg8_cax)
#remove duplicates
fc_cax_all <- fc_cax_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

log2_nox <- subset(log2_x, select = -c(hnrnpk, ckap4, syncrip, caprin1))
fc_pos2_nox <- log2_nox %>% filter(nono>=1)
#highest log2 ratio =2.76
fc_pos4_nox <- log2_nox %>% filter(nono>=2)
fc_pos8_nox <- log2_nox %>% filter(nono>=3)
fc_neg2_nox <- log2_nox %>% filter(nono<=-1)   
fc_neg4_nox <- log2_nox %>% filter(nono<=-2) 
fc_neg8_nox <- log2_nox %>% filter(nono<=-3)
fc_nox_all <- rbind(fc_pos2_nox, fc_pos4_nox, fc_pos8_nox, fc_neg2_nox, fc_neg4_nox, fc_neg8_nox)
#remove duplicates
fc_nox_all <- fc_nox_all %>%
  distinct(ensemble_id, .keep_all = TRUE)

#combine all fc (xgene) datatables into single, with id=kd
fc_hkx <- data.frame(ensemble_id = fc_hkx_all$ensemble_id, 
                     gene = fc_hkx_all$gene, 
                     chromosome = fc_hkx_all$chromosome,
                     length = fc_hkx_all$length,
                     TSS = fc_hkx_all$TSS,
                     kd_group = "hnrnpk",
                     log2_ratio = fc_hkx_all$hnrnpk)

fc_ckx <- data.frame(ensemble_id = fc_ckx_all$ensemble_id, 
                     gene = fc_ckx_all$gene, 
                     chromosome = fc_ckx_all$chromosome,
                     length = fc_ckx_all$length,
                     TSS = fc_ckx_all$TSS,
                     kd_group = "ckap4",
                     log2_ratio = fc_ckx_all$ckap4)

fc_syx <- data.frame(ensemble_id = fc_syx_all$ensemble_id, 
                     gene = fc_syx_all$gene, 
                     chromosome = fc_syx_all$chromosome,
                     length = fc_syx_all$length,
                     TSS = fc_syx_all$TSS,
                     kd_group = "syncrip",
                     log2_ratio = fc_syx_all$syncrip)

fc_cax <- data.frame(ensemble_id = fc_cax_all$ensemble_id, 
                     gene = fc_cax_all$gene, 
                     chromosome = fc_cax_all$chromosome,
                     length = fc_cax_all$length,
                     TSS = fc_cax_all$TSS,
                     kd_group = "caprin1",
                     log2_ratio = fc_cax_all$caprin1)

fc_nox <- data.frame(ensemble_id = fc_nox_all$ensemble_id, 
                     gene = fc_nox_all$gene, 
                     chromosome = fc_nox_all$chromosome,
                     length = fc_nox_all$length,
                     TSS = fc_nox_all$TSS,
                     kd_group = "nono",
                     log2_ratio = fc_nox_all$nono)

#rowbind into long list
fc_all_x <- rbind(fc_hkx, fc_ckx, fc_syx, fc_cax, fc_nox)

#scatterplot all together
plot_fc_all_x <- ggplot(fc_all_x, aes(TSS, log2_ratio))+
  geom_point(aes(color=kd_group), shape=18)+
  ylim(-4.2,4.2)+
  xlab("") 
plot_fc_all_x

#facet scatterplots
plot_fc_all_x_facet <- ggplot(fc_all_x, aes(TSS, log2_ratio))+
  geom_point(aes(color=log2_ratio), shape=18)+
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0 ,limits=c(-4,4))+
  ylim(-4.2,4.2)+
  xlab("TSS_X_chromosome")+
  ylab("log2[cpm_kd/cpm_U6M2")+
  facet_grid(kd_group ~.)
plot_fc_all_x_facet

#export fc_all_x
setwd("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/")
write.table(fc_all_x, quote=FALSE, sep="\t", file = "fold_change_all_x.txt", row.names=FALSE, col.names = TRUE)

#export gene list for GO analysis
write.table(fc_all_x$ensemble_id, quote=FALSE, sep="\t", file = "fold_change_genes_x.txt", row.names=FALSE, col.names = FALSE)





