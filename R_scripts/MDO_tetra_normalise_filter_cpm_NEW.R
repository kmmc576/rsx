#This to filter and normalise htseq counts, then create DGEList object
#This normalises first, then filters on cpm (PDW says use this order)  
#install edgeR
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install("edgeR")

#load edgeR library
library("edgeR")
library(dplyr)

#import htseq txt files + also ensemble_chromosomes (ensembl id, gene name, chromosome location) + MDO_genes_transcript_length.txt


#NB for Illumina Stranded mRNA prep, use -s (stranded) "reverse" (not "yes") in htseq_count
getwd()
setwd("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/")
#read in id info, with ensemble_id as row name (column 1)
ensemble_chromosomes <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/ensemble_chromosomes.txt", header=FALSE, sep="\t")
hiseq_u6m2 <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/MDO_sequence_files/MDO_tetra_htseq/MDO_tetra_hiseq_U6M2_strdrev.txt",header=FALSE, sep="\t")
hiseq_hnrnpk <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/MDO_sequence_files/MDO_tetra_htseq/MDO_tetra_hiseq_Hnrnpk_strdrev.txt",header=FALSE, sep="\t")
hiseq_ckap4 <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/MDO_sequence_files/MDO_tetra_htseq/MDO_tetra_hiseq_Ckap4_strdrev.txt",header=FALSE, sep="\t")
hiseq_syncrip <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/MDO_sequence_files/MDO_tetra_htseq/MDO_tetra_hiseq_Syncrip_strdrev.txt",header=FALSE, sep="\t")
hiseq_caprin1 <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/MDO_tetra_htseq/MDO_tetra_hiseq_Caprin1_strdrev.txt",header=FALSE, sep="\t")
hiseq_nono <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/MDO_sequence_files/MDO_tetra_htseq/MDO_tetra_hiseq_Nono_strdrev.txt",header=FALSE, sep="\t")


#combine _counts files (dataframes),with ensembl_id and remove last 5 rows of summary data from each
#use second (count) column of kd files
all_counts <- cbind(ensemble_chromosomes, hiseq_u6m2[1:34985,2], hiseq_hnrnpk[1:34985,2], hiseq_ckap4[1:34985,2], hiseq_syncrip[1:34985,2], hiseq_caprin1[1:34985,2], hiseq_nono[1:34985,2]) 

#rename columns
names(all_counts)[1] <- "ensemble_id"
names(all_counts)[2] <- "gene"
names(all_counts)[3] <- "chromosome"
names(all_counts)[4] <- "u6m2"
names(all_counts)[5] <- "hnrnpk"
names(all_counts)[6] <- "ckap4"
names(all_counts)[7] <- "syncrip"
names(all_counts)[8] <- "caprin1"
names(all_counts)[9] <- "nono"

#calculate total mapped counts for each
tot <- colSums(all_counts[ ,4:9])
tot
mean(tot)
#    u6m2   hnrnpk    ckap4  syncrip  caprin1     nono 
#23727731 24519553 22994956 23342913 23360607 20619671 

##extract and sort htseq statistics
htseq_stats <- cbind(hiseq_u6m2[34986:34990, ], hiseq_hnrnpk[34986:34990,2], hiseq_ckap4[34986:34990,2], hiseq_syncrip[34986:34990,2], hiseq_caprin1[34986:34990,2], hiseq_nono[34986:34990,2]) 
names(htseq_stats)[1] <- "htseq_total"
names(htseq_stats)[2] <- "u6m2"
names(htseq_stats)[3] <- "hnrnpk"
names(htseq_stats)[4] <- "ckap4"
names(htseq_stats)[5] <- "syncrip"
names(htseq_stats)[6] <- "caprin1"
names(htseq_stats)[7] <- "nono"
#insert total mapped reads (caluclated using colSums --> tot)
htseq_stats [6, 2:7] <- tot
htseq_stats [6, 1] <- "mapped reads"
mean_mapped_reads <- mean(tot)
mean_mapped_reads
#[1] 23094238

#check all_counts
write.table(all_counts, file = "all_counts_strdrev.txt", quote = FALSE, sep = "\t")

#filter to retain only genes with attached chromosome annotation (ie, not scaffolds)
all_counts_att <- all_counts %>% filter(chromosome=="X" | chromosome=="1" | chromosome=="2" | chromosome=="3" | chromosome=="4" | chromosome=="5" | chromosome=="6" | chromosome=="7" | chromosome=="8")

#convert all_counts (dataframe) to matrix 
#
#convert to DGEList object
#normalise (TMM)
#then filter on cpm

#convert to DGEList object
group = c("u6m2", "hnrnpk", "ckap4", "syncrip", "caprin1", "nono")
all_counts_DGE <- DGEList(counts=all_counts_att[ ,4:9], group=group, genes=all_counts_att[ ,1:3])

#normalise TMM (normalises to calculate effective library sizes). Does not change count values 
counts_norm <- calcNormFactors(all_counts_DGE)

#list counts
head(counts_norm$counts)
#list group:NB NULL because no replicates 
head(counts_norm$group)
#list gene id and chromosome
head(counts_norm$gene)
#list samples (groups), with library size- same as tot (calculated above)
head(counts_norm$samples)

el<- effectiveLibSizes(counts_norm)
max(el)
#Filter (cpm) for low counts after normalization
#aim for 5-10 reads per library. 

#Largest effective library size 23848036 mapped reads. In cpm, 5 reads is this is 0.210; 10 reads is 0.419
#use 0.4 - this is 9.54 reads in largest library
#TMM normalization in calcNormFactors removes composition biases that remain after library size normalization. 
#BUT low counts reduce the accuracy of TMM normalization,so best to do it after filtering for low counts?? PDW says no
#filter for cpm>0.4, in at least 1 sample (maybe too lax) 
#>=1 sample - F(19639) T(13752)
#>=2 samples- F(20068) T(13323)
#>=5 samples - F(20789) T(12602)
#>=6 samples - F(21081) T(12310)

keep <- rowSums(cpm(counts_norm)>1) >= 1
summary(keep)
#  for cpm>0.4
# FALSE    TRUE 
#  19639   13752 
#for cpm>1 - lose 1511
# FALSE    TRUE 
# 21150   12241
#for cpm>10
#FALSE    TRUE 
# 24656    8735

counts_filtered <- counts_norm[keep, , keep.lib.sizes=FALSE]

#Add in TSS:
#ensembl BioMart annotates TSS and also "transcript start" and "transcript end".
#what is diff betwen TSS and transcription start??- use TSS
head(counts_filtered)
write.table(counts_filtered$genes$ensemble_id, quote=FALSE, sep="\t", file = "counts_filtered_gene_list.txt", row.names=FALSE)
#extract TSS data from ensembl BioMart andimport
genes_TSS <- read.table("/Users/kimmcintyre/Desktop/mart_export-2.txt", sep = "\t", header=TRUE)

names(genes_TSS)[1] <- "ensemble_id"
names(genes_TSS)[2] <- "length"
names(genes_TSS)[3] <- "transcript_start"
names(genes_TSS)[4] <- "transcript_end"
names(genes_TSS)[5] <- "TSS"

dim(genes_TSS)
#filter for longest transcript for each id
filt_genes_TSS <- genes_TSS %>% 
  group_by(ensemble_id) %>%
  filter(length == max(length))
#filter for duplicates
filt_genes_TSS <- filt_genes_TSS %>% 
  distinct(ensemble_id, .keep_all = TRUE)
#add length. Dplyr: inner_join() combines the columns from both data frames, but only keeps rows where the value in the key column matches in both data frames.
counts_filtered$genes <- inner_join(counts_filtered$genes, filt_genes_TSS, by="ensemble_id")
head(counts_filtered)


class(counts_filtered)
dim(counts_filtered)

#counts_filtered$counts
#counts_filtered$genes
counts_filtered$samples

#NB library size now adjusted (cf tot) by TMM
#DGEListobject$counts just lists original raw counts, library sizes change to "effective" lib size with normalisation

head(counts_filtered)

counts_normal_cpm <- cpm(counts_filtered, normalized.lib.sizes = TRUE)
#counts_filtered remains a DGEList object
#counts_normal_cpm is matrix of CPM
head(counts_normal_cpm)
class(counts_normal_cpm)
write.table(counts_filtered$genes, quote=FALSE, sep="\t", file = "counts_filtered_genes_cpm1.txt", row.names=TRUE)
write.table(counts_normal_cpm, quote=FALSE, sep="\t", file = "counts_normal_cpm_cpm1.txt", row.names=TRUE)
write.table(counts_filtered$counts, quote=FALSE, sep="\t", file = "counts_filtered_counts_cpm1.txt", row.names=TRUE)

