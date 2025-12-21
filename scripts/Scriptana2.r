
library(dada2)
library(phyloseq)
library(ggplot2)
library(icesTAF)
library(ips)
library(tidyverse)
library(metagMisc)
library(decontam)
library(DECIPHER)
library(phangorn)
# set dir for results
path = choose.dir(default = "", caption = "Choose local project folder")
setwd(path)
# Define the path to your demultiplexed samples
list.files(path)

## Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern = "_L001_R1_001.fastq.gz"))
fnRs <- sort(list.files(path, pattern = "_L001_R2_001.fastq.gz"))
## Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
## Specify the full path to the fnFs
fnFs <- file.path(path, fnFs)
fnFs
fnRs <- file.path(path, fnRs)
fnRs

setwd("*/Results")
f1 <- plotQualityProfile(fnFs[1:12])
r1 <- plotQualityProfile(fnRs[1:12])
f1
r1
ggsave(paste("./f1.pdf"), f1)
ggsave(paste("./r1.pdf"), r1)
# First we define the filenames for the filtered fastq.gz files:
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
  trimLeft = c(21, 20), truncLen = c(285, 220),
  maxN = 0, maxEE = c(2, 2),
  truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = TRUE
)
head(out)
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
setwd("$path/Results/EE1")
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
# It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:
errFplot <- plotErrors(errF, nominalQ = TRUE) ### inspect the fit between the observed error rates (black points) and the fitted error rates (black lines)
errRplot <- plotErrors(errR, nominalQ = TRUE)
errFplot
errRplot
ggsave(paste("./errorplot.pdf"), errFplot)
ggsave(paste("./errorplot.pdf"), errRplot)
saveRDS(errF, "./errF.rds")
saveRDS(errR, "./errR.rds")
setwd("/Results")
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
saveRDS(dadaFs, "./dadaFs.rds")
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)
saveRDS(dadaRs, "./dadaRs.rds")
setwd("./Results")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE) # Non-overlapping reads are supported, but not recommended, with mergePairs(..., justConcatenate=TRUE)
head(mergers[[1]]) # Inspect the merger data.frame from the first sample
saveRDS(mergers, "./mergers.rds")
setwd("./Results")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
write.table(seqtab, "./seqtab.txt", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(seqtab, "./seqtab.rds")
table(nchar(getSequences(seqtab)))
setwd("./Results")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)
setwd("./Results")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "./track_reads.txt", sep = "\t", row.names = FALSE, quote = FALSE)
setwd("./Results")
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load(config$taxonomydatabase_IDTAXA_bac_dir_local)
ids <- IdTaxa(dna, trainingSet, strand = "top", processors = 20, verbose = FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
saveRDS(taxid, "./taxid.rds")
write.table(taxid, "./taxid_silva.txt", sep = "\t", row.names = FALSE, quote = FALSE)
setwd("./Results")
taxa_si <- assignTaxonomy(seqtab.nochim, config$taxonomydatabase_NavBaysian_dir_local, multithread = FALSE, outputBootstraps = TRUE)
write.table(taxa_si, "./taxa_si.txt", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(taxa_si, "./taxa_si.rds")
### Add species level assignment
taxa_si_sp <- addSpecies(taxa_si[[1]], config$taxonomydatabase_NavBaysiansp_dir_local)
write.table(taxa_si_sp, "./taxa_si_sp.txt", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(taxa_si_sp, "./taxa_si_sp.rds")
if (config$mock == 1) {
  setwd(config$input_dir_local)
  unqs.mock <- seqtab.nochim["MockB",]
  unqs.mock <- sort(unqs.mock[unqs.mock > 0], decreasing = TRUE) # Drop ASVs absent in the Mock
  cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
  mock.ref <- getSequences(file.path(config$mockcommunity_dir_local, "mock.fasta"))
  match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
  cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
}
setwd("./Results")
uniquesToFasta(seqtab.nochim, "./rep_set.fna", ids = paste0("SV_", seq(length(getUniques(seqtab.nochim)))))
setwd(config$input_dir_local)
samdf <- read.table("mapping.txt", header = TRUE, sep = "\t", row.names = 1)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(samdf), tax_table(taxa_si_sp))
# substitute the actual sequence by ASV. It does also change it from the taxonomy table. https://github.com/joey711/phyloseq/issues/833
# taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
t <- t(otu_table(ps))
d <- cbind(rownames(t), data.frame(t, row.names = NULL))
colnames(d)[1] <- "#OTU ID" # Add '#OTUID' to the header (required by later biom)
setwd("./Results")
write.table(d, "./ASV_table_ps.txt", sep = "\t")
ps
saveRDS(ps, "./ps.rds")
final_table_taxa_si_sp <- phyloseq_to_df(ps) # you can't have - in sample names
write.table(final_table_taxa_si_sp, "./final_table_taxa_si_sp.txt", sep = "\t", quote = F, row.names = FALSE)
setwd(config$input_dir_local)
samdf <- read.table("mapping.txt", header = TRUE, sep = "\t", row.names = 1)
ps_taxid <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(samdf), tax_table(taxid))
# substitute the actual sequence by ASV. It does also change it from the taxonomy table. https_taxid://github.com/joey711/phyloseq/issues/833
# taxa_names(ps_taxid) <- paste0("ASV_", seq(ntaxa(ps_taxid)))
t <- t(otu_table(ps_taxid))
d <- cbind(rownames(t), data.frame(t, row.names = NULL))
colnames(d)[1] <- "#OTU ID" # Add '#OTUID' to the header (required by later phyloseq)
setwd("./Results")
write.table(d, "./ASV_table_ps_taxid.txt", sep = "\t")
ps_taxid
saveRDS(ps_taxid, "./ps_taxid.rds")
final_table_taxid <- phyloseq_to_df(ps_taxid) # you can't have - in sample names
write.table(final_table_taxid, "./final_table_taxid.txt", sep = "\t", quote = F, row.names = FALSE)
if (config$decontamblank == 1) {
  setwd("./Results")
  # Inspect library sizes
  df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  # Contaminant identification using prevalence
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
  contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg")
  table(contamdf.prev$contaminant)
  head(which(contamdf.prev$contaminant))
  contamdf.prev05 <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5) # Note that as before, the default threshold for a contaminant is that it reaches a probability of 0.1 in the statistical test being performed. In the prevalence test there is a special value worth knowing, threshold=0.5, that will identify as contaminants all sequences thare are more prevalent in negative controls than in positive samples.
  table(contamdf.prev05$contaminant)
  # Make phyloseq object of presence-absence in negative controls and true samples
  ps.pa <- transform_sample_counts(ps, function(abund) 1 * (abund > 0))
  ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
  ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
  # Make data.frame of prevalence in positive and negative samples
  df.pa <- data.frame(
    pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg),
    contaminant = contamdf.prev$contaminant
  )
  plot <- ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
    geom_point() +
    xlab("Prevalence (Negative Controls)") +
    ylab("Prevalence (True Samples)")
  ggsave(paste("./Prevalence_neg_controls_truesamples.pdf", sep = ""), plot)
  # now let's remove them
  ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
  ps.noncontam
  saveRDS(ps.noncontam, "./ps.noncontam.rds")
  final_table_taxa_si_sp_noncontam <- phyloseq_to_df(ps.noncontam)
  write.table(final_table_taxa_si_sp_noncontam, "./final_table_taxa_si_sp_noncontam.txt", sep = "\t", quote = F, row.names = FALSE)
}
if (config$decontamblank == 1) {
  setwd("./Results")
  # Inspect library sizes
  df <- as.data.frame(sample_data(ps_taxid)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps_taxid)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  # Contaminant identification using prevalence
  sample_data(ps_taxid)$is.neg <- sample_data(ps_taxid)$Sample_or_Control == "Control Sample"
  contamdf.prev <- isContaminant(ps_taxid, method = "prevalence", neg = "is.neg")
  table(contamdf.prev$contaminant)
  head(which(contamdf.prev$contaminant))
  contamdf.prev05 <- isContaminant(ps_taxid, method = "prevalence", neg = "is.neg", threshold = 0.5) # Note that as before, the default threshold for a contaminant is that it reaches a probability of 0.1 in the statistical test being performed. In the prevalence test there is a special value worth knowing, threshold=0.5, that will identify as contaminants all sequences thare are more prevalent in negative controls than in positive samples.
  table(contamdf.prev05$contaminant)
  # Make phyloseq object of presence-absence in negative controls and true samples
  ps_taxid.pa <- transform_sample_counts(ps_taxid, function(abund) 1 * (abund > 0))
  ps_taxid.pa.neg <- prune_samples(sample_data(ps_taxid.pa)$Sample_or_Control == "Control Sample", ps_taxid.pa)
  ps_taxid.pa.pos <- prune_samples(sample_data(ps_taxid.pa)$Sample_or_Control == "True Sample", ps_taxid.pa)
  # Make data.frame of prevalence in positive and negative samples
  df.pa <- data.frame(
    pa.pos = taxa_sums(ps_taxid.pa.pos), pa.neg = taxa_sums(ps_taxid.pa.neg),
    contaminant = contamdf.prev$contaminant
  )
  plot <- ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
    geom_point() +
    xlab("Prevalence (Negative Controls)") +
    ylab("Prevalence (True Samples)")
  ggsave(paste("./Prevalence_neg_controls_truesamples_taxid.pdf", sep = ""), plot)
  ps_taxid.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps_taxid)
  ps_taxid.noncontam
  saveRDS(ps_taxid.noncontam, "./ps_taxid.noncontam.rds")
  final_table_taxid_noncontam <- phyloseq_to_df(ps_taxid.noncontam)
  write.table(final_table_taxid_noncontam, "./final_table_taxid_noncontam.txt", sep = "\t", quote = F, row.names = FALSE)
}
if (config$phylo == 1) {
  setwd("./Results")
  seqs <- getSequences(seqtab.nochim)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
  saveRDS(alignment, "./alignment.rds")
  # Phangorn
  phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit <- pml(treeNJ, data = phang.align)
  # negative edges length changed to 0!
  fitGTR <- update(fit, k = 4, inv = 0.2)
  saveRDS(fitGTR, "./fitGTR.rds")
  fitGTR2 <- optim.pml(fitGTR,
    model = "GTR", optInv = TRUE, optGamma = TRUE,
    rearrangement = "stochastic", control = pml.control(trace = 0)
  )
  saveRDS(fitGTR2, "./fitGTR2.rds")
}