#This to filter and normalise htseq counts, then create DGEList object
#This filters on cpm before EdgeR normalisation
#PDW says to normalise first, then filter  
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

#aim for 5-10 reads per library. Average library size 23,094,239 mapped reads. 
#In cpm, this is 0.217 - 0.433 -> use 0.45

#check all_counts
write.table(all_counts, file = "all_counts_strdrev.txt", quote = FALSE, sep = "\t")

#import to excel to view/check using sort

#convert all_counts (dataframe) to matrix 
#convert to DGEList object
#filter on cpm
#then normalise (TMM)

#convert to DGEList object
group = c("u6m2", "hnrnpk", "ckap4", "syncrip", "caprin1", "nono")
all_counts_DGE <- DGEList(counts=all_counts[ ,4:9], group=group, genes=all_counts[ ,1:3])

#list counts
all_counts_DGE$counts
#list group:NB NULL because no replicates 
all_counts_DGE$group
#list gene id and chromosome
all_counts_DGE$gene
#list samples (groups), with library size- same as tot (calculated above)
all_counts_DGE$samples


#Filter (cpm) for low counts before normalization
#TMM normalization in calcNormFactors removes composition biases that remain after library size normalization. 
#BUT low counts reduce the accuracy of TMM normalization,so best to do it after filtering for low counts
#filter for cpm>0.45, in at least 1 sample (maybe too lax) - use cpm > 0.5
#>=1 sample - F(21007) T(13978)
#>=2 samples- F(21433) T(13552)
#>=5 samples - F(22090) T(12895)
#>=6 samples - F(22396) T(12589)

keep <- rowSums(cpm(all_counts_DGE)>0.5) >= 1
summary(keep)
#  Mode   FALSE    TRUE 
# logical   21185   13800 

counts_filtered <- all_counts_DGE[keep, , keep.lib.sizes=FALSE]
class(counts_filtered)
dim(counts_filtered)
#normalise TMM (normalises to calculate effective library sizes). Does not change count values 
counts_filtered <- calcNormFactors(counts_filtered)
counts_filtered$counts
counts_filtered$genes
counts_filtered$samples

#NB library size now adjusted (cf tot) by TMM

head(counts_filtered)

counts_normal_cpm <- cpm(counts_filtered, normalized.lib.sizes = TRUE)

head(counts_normal_cpm)

write.table(counts_filtered$genes, quote=FALSE, sep="\t", file = "counts_filtered_genes.txt", row.names=TRUE)
write.table(counts_normal_cpm, quote=FALSE, sep="\t", file = "counts_normal_cpm.txt", row.names=TRUE)
write.table(counts_filtered$counts, quote=FALSE, sep="\t", file = "counts_filtered_counts.txt", row.names=TRUE)
