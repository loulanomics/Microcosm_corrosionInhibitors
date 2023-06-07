########################################
### V4 rRNA sequencing processing for 
### for Kimbell, et al., 2023
### Lou LaMartina, finalized Jun 7, 2023
########################################


setwd("~/Desktop/16S/Kimbell")
library(dada2)
library(ggplot2)


###################################
### prepare data for processing ###
###################################

# set working directory
setwd("~/Desktop/Lab/Projects/Others/Kimbell/Microcosm")


# set file paths
path <- "./cutadapt"
pathF <- "./cutadapt/fastqF"
pathR <- "./cutadapt/fastqR"


# set paths for filtered reads
filtered_pathF <- "./cutadapt/fastqF/Filtered"
filtered_pathR <- "./cutadapt/fastqR/Filtered"


# sort forward and reverse reads into file paths
fastqFs <- sort(list.files(pathF, pattern = "_R1_", full.names = TRUE))
fastqRs <- sort(list.files(pathR, pattern = "_R2_", full.names = TRUE))


# extract file names
samples <- sapply(strsplit(basename(fastqFs), "_"), '[', 1)




################################
### Inspect sequence quality ###
################################

# visualize quality of reads
qualityF <- plotQualityProfile(fastqFs[1:4])#; qualityF
qualityR <- plotQualityProfile(fastqRs[1:4])#; qualityR


# save quality profiles
ggsave("./Plots/qualityF.pdf", plot = qualityF, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/qualityR.pdf", plot = qualityR, device = "pdf", width = 12, height = 8, units = "in")


# check: if there are not the same number of F and R files, stop.
if(length(fastqFs) != length(fastqRs)) 
  stop("Forward and reverse files do not match.")


# inspect sequence distribution - this takes a few mins
fastq_lengths1.ls <- list()
for(i in 1:length(fastqFs)){
  fastq_lengths1.ls[[i]] <- data.frame(table(nchar(getSequences(fastqFs[i]))))
}

fastq_lengths1 <- do.call(rbind, fastq_lengths1.ls)
colnames(fastq_lengths1) <- c("SeqLength", "Frequency")
fastq_lengths1$SeqLength <- as.numeric(as.character(fastq_lengths1$SeqLength))


# visualize
dist1 <-
  ggplot(fastq_lengths1, aes(x = SeqLength, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black")
dist1

ggsave("./Plots/seq_distribution_beforeFilt.pdf", plot = dist1, device = "pdf",
       width = 10, height = 4, units = "in")




#######################
### Quality control ###
#######################

# give filtered files new names and paths
filtFs <- file.path(filtered_pathF, paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_pathR, paste0(samples, "_R_filt.fastq.gz"))


# filter based on quality and read length
start = Sys.time()
filtered_out <- filterAndTrim(fastqFs, filtFs, fastqRs, filtRs, matchIDs = TRUE,
                              maxEE = 2, maxN = 0, truncQ = 10,
                              rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)
Sys.time() - start
# 17.6167 mins


# inspect how many reads were filtered out of each sample
mean((1 - (filtered_out[,2] / filtered_out[,1])) * 100)
# [1] 24.03875 % removed


# plot quality profiles of filtered reads
filtF.plot <- plotQualityProfile(filtFs[1:4]); filtF.plot
filtR.plot <- plotQualityProfile(filtRs[1:4]); filtR.plot


# save quality profiles
ggsave("./Plots/filt_qualityF.pdf", plot = filtF.plot, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/filt_qualityR.pdf", plot = filtR.plot, device = "pdf", width = 12, height = 8, units = "in")


# set sample names to the ID only
names(filtFs) <- samples
names(filtRs) <- samples


# inspect sequence distribution
fastq_lengths2.ls <- list()
for(i in 1:length(filtFs)){
  fastq_lengths2.ls[[i]] <- data.frame(table(nchar(getSequences(filtFs[i]))))
}

fastq_lengths2 <- do.call(rbind, fastq_lengths2.ls)
colnames(fastq_lengths2) <- c("SeqLength", "Frequency")
fastq_lengths2$SeqLength <- as.numeric(as.character(fastq_lengths2$SeqLength))


# visualize
dist2 <-
  ggplot(fastq_lengths2, aes(x = SeqLength, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black")
dist2

ggsave("./Plots/seq_distribution_afterFilt.pdf", plot = dist2, device = "pdf", 
       width = 10, height = 4, units = "in")


# save progress
save.image("./RData/Microcosm_dada2_env.RData")




############################
### Learning error rates ###
############################

# learn and visualize error rates of F reads
start = Sys.time()
errorF <- learnErrors(filtFs, multithread = TRUE)
errF.plot <- plotErrors(errorF, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errF.plot


# learn and visualize error rates of R reads
errorR <- learnErrors(filtRs, multithread = TRUE)
errR.plot <- plotErrors(errorR, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errR.plot
Sys.time() - start
# 12.9828 mins


# save error plots
ggsave("./Plots/errorF.pdf", plot = errF.plot, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/errorR.pdf", plot = errR.plot, device = "pdf", width = 12, height = 8, units = "in")


# save progress
save.image("./RData/Microcosm_dada2_env.RData")




################################
### Merging paired-end reads ###
################################

# create list of merged reads
mergers <- vector("list", length(samples))
names(mergers) <- samples


# sample inference and merging paired-end reads
start = Sys.time()
for(i in samples) {
  cat("\nProcessing:", i, "(", match(i, samples), ") :", format(Sys.time(), "%H:%M %p"), "\n")
  derepF <- derepFastq(filtFs[[i]])
  dadaF <- dada(derepF, err = errorF, multithread = TRUE)
  derepR <- derepFastq(filtRs[[i]])
  dadaR <- dada(derepR, err = errorR, multithread = TRUE)
  Merger <- mergePairs(dadaF, derepF, dadaR, derepR)
  mergers[[i]] <- Merger
}
Sys.time() - start
# 1.006945 hours


# removing dereps to save memory
rm(derepF, derepR)


# contruct a sequence table
counts <- makeSequenceTable(mergers)


# dimensions: num rows x num columns
dim(counts)
# [1]    81 12888


# save progress
save.image("./RData/Microcosm_dada2_env.RData")




##################################
### Quality control: processed ###
##################################


########
### trim

# inspect sequence distribution
seq_distribution <- data.frame(table(nchar(getSequences(counts)))) 
colnames(seq_distribution) <- c("SeqLength", "Frequency")
seq_distribution$SeqLength <- as.numeric(as.character(seq_distribution$SeqLength))


# visualize and save
dist3 <- 
  ggplot(seq_distribution, aes(x = SeqLength, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black")
dist3

ggsave("./Plots/seq_distribution_afterMerge.pdf", plot = dist3, device = "pdf", width = 10, height = 4, units = "in")


# remove reads of non target length, 5% above and below the median 
median(nchar(getSequences(counts)))
# [1] 253
min_len <- floor(median(nchar(getSequences(counts))) * 0.95); min_len
# [1] 240
max_len <- ceiling(median(nchar(getSequences(counts))) * 1.05); max_len
# [1] 266


# modify sequence table with new guidelines
counts_trimmed <- counts[ ,nchar(colnames(counts)) 
                                     %in% seq(min_len, max_len)]


dim(counts_trimmed)
# [1]    81 11383



###################
### Remove chimeras

# removing chimeras with denovo screening
start = Sys.time()
counts_nochim <- removeBimeraDenovo(counts_trimmed, method = "consensus",
                                            multithread = TRUE, verbose = TRUE)
Sys.time() - start
#  53.06822 secs


# how many unique sequences were moved?
ncol(counts_trimmed)               # [1] 11383
ncol(counts_nochim)        # [1] 11246


# what percentage of reads were identified as chimeras?
(1 - sum(counts_nochim) / sum(counts_trimmed)) * 100
# [1] 0.2946269 %


# save progress
save.image("./RData/Microcosm_dada2_env.RData")




#######################
### Assign taxonomy ###
#######################

taxa <- assignTaxonomy(counts_nochim,
                                      "~/Desktop/Lab/Projects/Misc/silva_nr_v132_train_set.fa.gz",
                                      multithread = TRUE)


# save progress
save.image("./RData/Microcosm_dada2_env.RData")


# identify non-bacterial/archael ASVs
taxa.df <- data.frame(taxa)
taxa.df <- data.frame(FASTA = rownames(taxa.df), taxa.df)
nonspecifics <- rbind(subset(taxa.df, Kingdom == "Eukaryota"),
                      subset(taxa.df, Order == "Chloroplast"),
                      subset(taxa.df, Family == "Mitochondria"))
nonspecifics$Source <- "not bacteria"


### no negative (no template) control was made ###

### no mock community was made ###




#########################
### organize and save ###
#########################

# remove contaminants
counts.df <- data.frame(counts_nochim)
identical(rownames(taxa.df), colnames(counts.df))


# simplify ASV names
taxa.df <- data.frame(ASV = paste0("ASV", sprintf("%05d", 1:ncol(counts.df))), taxa.df)


# create FASTA
uniquesToFasta(as.matrix(counts.df),
               ids = paste0(taxa.df$ASV, "__",
                            taxa.df$Kingdom, "__",
                            taxa.df$Phylum, "__",
                            taxa.df$Class, "__",
                            taxa.df$Order, "__",
                            taxa.df$Family, "__",
                            taxa.df$Genus),
               fout = "./RData/Microcosm_v4.fasta")


# change ASVs in counts
identical(rownames(taxa.df), colnames(counts.df))
colnames(counts.df) <- taxa.df$ASV
counts.df <- data.frame(Sample_name = rownames(counts.df), counts.df)


# save
write.csv(counts.df, "./RData/Microcosm_ASV_counts.csv", row.names = F, na = "")
write.csv(taxa.df, "./RData/Microcosm_ASV_taxa.csv", row.names = F, na = "")


# save progress
save.image("./RData/Microcosm_dada2_env.RData")
