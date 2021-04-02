#calculate X:A ratio by dividing the median expression (TPM) of X-linked genes with the median expression (TPM) of the autosomal genes.
#ref: D Chandel, C H Naik,  View ORCID ProfileS Gayen
#read in gene lengths (transcript length including UTR and CDS)from ensemble BioMart
#TPM = transcripts per million
#1)Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#2)Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#3)Divide the RPK values by the “per million” scaling factor. This gives you TPM.
#4)similar to RPKM, but sum of all TPMs in each sample are the same. This makes it easier to compare the proportion of reads that mapped to a gene in each sample. 
#In contrast, with RPKM and FPKM, the sum of the normalized reads in each sample may be different
#when calculating TPM, the only difference is that you normalize for gene length first, and then normalize for sequencing depth second.

#First, get transcript lengths from BioMart
#create list of genes
filt_genes <- data.frame(counts_filtered$genes$ensemble_id)
head(filt_genes)
write.table(filt_genes, file="/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/filtered_genes.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#read in transcript lengths
transcript_length <- read.table("/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/transcript_length.txt", header=TRUE, sep="\t")
#rename columns
names(transcript_length)[1] <- "ensemble_id"
names(transcript_length)[2] <- "length"

#ensemble_info transcript lengths: 26139 obs cf 13800 filtered genes (more than one for some ids) --> select longest for each id
#confirm 13800 unique elements
unique_id <- unique(transcript_length[ ,1])

#filter gene length to retain only longest transcript for each ensemble id
#returns 13865 elements, 13800 unique
#inspect/edit manuallu
filt_transcript_length <- transcript_length %>% 
  group_by(ensemble_id) %>%
  filter(length == max(length))
unique_id2 <- unique(filt_transcript_length[ ,1])


#inspection in excel shows duplicates -> filter for duplicates - yields 13800 unique ids
filt_transcript_length <- filt_transcript_length %>% 
  distinct(ensemble_id, .keep_all = TRUE)
unique_id2 <- unique(filt_transcript_length[ ,1])

write.table(filt_transcript_length, file="/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seq/filtered_transcript_length.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#convert gene length to kb
filt_transcript_length$length <- filt_transcript_length$length/1000
head(filt_transcript_length) 

#add length. Dplyr: inner_join() combines the columns from both data frames, but only keeps rows where the value in the key column matches in both data frames.
counts_filtered$genes <- inner_join(counts_filtered$genes, filt_transcript_length, by="ensemble_id")
head(counts_filtered$genes)

# Start from feature summarization data (as an edgeR DGEList object)

#counts_filtered (DGEList of filtered counts), now with lengths added to $genes element

# Filter-out low expression features (already done) and perform sample normalization
counts_filtered <- edgeR::calcNormFactors(counts_filtered)

# Calculate TPM from normalized counts

calc_tpm <- function(x, gene.length) {
  x <- as.matrix(x)
  len.norm.lib.size <- colSums(x / gene.length)
  return((t(t(x) / len.norm.lib.size) * 1e06) / gene.length)
}

tpm_man <- calc_tpm(counts_filtered, gene.length = counts_filtered$genes$length)

head(tpm_man)
head(counts_filtered$counts)
write.table(tpm_man, file = "tpm.txt", quote = FALSE, sep = "\t")
write.table(counts_filtered$counts, file = "counts_filtered$counts.txt", quote = FALSE, sep = "\t")
write.table(counts_filtered$genes, file = "counts_filtered$genes.txt", quote = FALSE, sep = "\t")
#checked TPM calculation manually (in excel)- identical: "/Users/kimmcintyre/Desktop/PhD/RNA_seq/MDO_tetra_RNA-seqMDO_tetra_TPM_calcs"

