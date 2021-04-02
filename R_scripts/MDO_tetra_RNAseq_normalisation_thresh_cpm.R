#install edgeR
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install("edgeR")

#load edgeR library
library("edgeR")

#import htseq txt files + also ensemble_chromosomes (ensembl id, gene name, chromosome location)
#NB for Illumina Stranded mRNA prep, use -s (stranded) "reverse" (not "yes") in htseq_count

#_counts files are dataframes
class(MDO_tetra_hiseq_U6M2)

#combine _counts files (dataframes),with ensembl_id and remove last 5 rows of summary data from each
#use both columns of U6M2,and only second (count) column of kd files
all_counts <- cbind(ensemble_chromosomes, MDO_tetra_hiseq_U6M2[1:34985,2], MDO_tetra_hiseq_Hnrnpk[1:34985,2], MDO_tetra_hiseq_Ckap4[1:34985,2], MDO_tetra_hiseq_Syncrip[1:34985,2], MDO_tetra_hiseq_Caprin1[1:34985,2], MDO_tetra_hiseq_Nono[1:34985,2]) 

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

#check
write.table(all_counts, file = "all_counts.txt", quote = FALSE, sep = "\t")

#import to excel to view/check using sort

#convert all_counts (dataframe) to matrix -> DGEList
#convert to DGElist
counts = as.matrix(all_counts[ ,4:9])
group = c("u6m2", "hnrnpk", "ckap4", "syncrip", "caprin1", "nono")
genes = all_counts[ ,1:3]
all_counts_DGE <- DGEList(counts=counts, group=group, genes=genes)
#list counts
all_counts_DGE$counts
#list group:NB NULL because no replicates 
all_counts_DGE$group
#list gene id and chromosome
all_counts_DGE$gene
#list samples (groups), library size, norm.factors
all_counts_DGE$samples

write.table(genes, file = "genes.txt", quote=FALSE,sep="\t")



#aim for 5-10 reads per library. Average library size ~320,000 mapped reads. Check why so low??
#In cpm, this is 15.6 - 31.25 - use 30 for row sum (across all samples) +  - maybe a bit too low??

#NB: edgeR filterByExpr filters based in cpm above k (min.count) in n samples (n determined by group sample size in design matrix - one here)
#result TRUE returns only where all samples have cpm >=5
#need to input design: design <- model.matrix(~0+group)
#BUT I want to keep where ANY sample is cpm >= 5

#calculate cpm using normalised library size (??default TMM - trimmed mean of M values (weighted according to inverse variances)??)
#then apply filtering

all_counts_cpm <- cpm(all_counts_DGE, normalized.lib.sizes = TRUE)

#all_counts_cpm is a matrix
class(all_counts_cpm)
#check dimensions of matrix
dim(all_counts_cpm)

#export as txt to open in excel
write.table(all_counts_cpm, file = "all_counts_cpm.txt", quote=FALSE,sep="\t")

#the following returns logical matrix TRUE or FALSE for each 
thresh <- all_counts_cpm >= 15

#table of how many rows (genes) have TRUE
table(rowSums(thresh))

#keep rows with at least 2 TRUE (cpm >=15)
keep <- rowSums(thresh) >=2
counts_cpm_keep <- all_counts_cpm[keep,]
summary(keep)
dim(counts_cpm_keep)

#merge filtered count rows with ensemble id and chromosome number and gene names
#use merge (for data frames), merging on row number
#convert filtered counts (matrix) to dataframe
counts_filt_df <- as.data.frame(counts_cpm_keep)

write.table(counts_filt_df, file = "counts_filt_df.txt", quote=FALSE,sep="\t")

filtered <- merge(genes, counts_filt_df)
#ran merge on R console- works differently can rename column of original row names as new column, then sort by
#vector memory exhausted - try on Katana??

