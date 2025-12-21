renv::load()
library(dada2);
library(phyloseqGraphTest);
library(phyloseq)
library(ggplot2)
library(icesTAF)
library(ips)
library(tidyverse)
library(metagMisc)
library(decontam)
library(DECIPHER)
library(phangorn)
#set paths to preference or unassign (delete lines) to prompt selection
wdpath <- "C:/Users/m/Google Drive/GAPMananas/Trials/First_Trial"
dbpath <- "C:/Users/m/Documents/Databases"
if (!exists("wdpath")) {
  wdpath <- choose.dir(default = "", caption = "Choose work directory")
}
if (!exists("dbpath")) {
  dbpath <- choose.dir(default = "", caption = "Choose database path")
}
inpath <- file.path(wdpath, "Input")
respath <- file.path(wdpath, "Results")
rrpath <- file.path(wdpath, "Raw_reads")
#tspath <- system.file(dbpath, "silva_nr_v132_train_set.fa.gz", package = "dada2")
tspath <- file.path(dbpath, "silva_nr_v132_train_set.fa.gz")
sppath <- file.path(dbpath, "silva_species_assignment_v132.fa.gz")
shpath <- file.path(dbpath, "microgreen_algaebase_biocomPipe.fasta")
unpath <- file.path(dbpath, "UNITE_public.gz")
chopath <- shpath
if (!exists("chopath")) {
  x <- as.numeric(readline('for sh: "1" \nfor unite: "2" \n'))
  if (x == 1) { chopath <- shpath }
  if (x == 2) { chopath <- unpath }
}
filtpath <- file.path(rrpath, "filtered")

setwd(respath)
list.files(rrpath)
fnFs <- sort(list.files(rrpath, pattern = "_L001_R1_001.fastq.gz"));
fnRs <- sort(list.files(rrpath, pattern = "_L001_R2_001.fastq.gz"));
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1);
fnFs <- file.path(rrpath, fnFs);
fnFs;
fnRs <- file.path(rrpath, fnRs);
fnRs;
if (!file.exists("./dadaRs.rds")) {
  if (!file.exists("./r1.pdf")) {
    f1 = plotQualityProfile(fnFs[1:12]);
    r1 = plotQualityProfile(fnRs[1:12]);
    f1;
    r1;
    ggsave(paste("./f1.pdf"), f1);
    ggsave(paste("./r1.pdf"), r1);
  }
  if (!file.exists("./filtered")) {
    filtpath <- file.path(respath, "filtered")
    filtFs <- file.path(filtpath, paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(filtpath, paste0(sample.names, "_R_filt.fastq.gz"))
    derepFs <- derepFastq(filtFs, verbose = TRUE)
    derepRs <- derepFastq(filtRs, verbose = TRUE)
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(21, 20), truncLen = c(285, 220),
              maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
              compress = TRUE, multithread = FALSE)
    head(out)
  }
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(220, 200), maxN = 0, maxEE = c(3, 4), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
  head(out)
  errF <- learnErrors(filtFs, multithread = FALSE)
  errR <- learnErrors(filtRs, multithread = FALSE)
  errFplot <- plotErrors(errF, nominalQ = TRUE)
  errFplot
  errRplot <- plotErrors(errR, nominalQ = TRUE)
  errRplot
  ggsave(paste("./errorplot.pdf"), errFplot)
  ggsave(paste("./errorplot.pdf"), errRplot)
  saveRDS(errF, "./errF.rds")
  saveRDS(errR, "./errR.rds")
  dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
  saveRDS(dadaFs, "./dadaFs.rds")
  dadaFs[[1]]
  dadaRs <- dada(derepRs, err = errR, multithread = FALSE)
  saveRDS(dadaRs, "./dadaRs.rds")
}
filtpath <- file.path(path, "filtered")
filtFs <- sort(list.files(filtpath, pattern = "_F_filt.fastq.gz"))
filtRs <- sort(list.files(filtpath, pattern = "_R_filt.fastq.gz"))
Mock <- sort(list.files(filtpath, pattern = "Zymo-mock"))
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(21, 20), truncLen = c(285, 220),
              maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
              compress = TRUE, multithread = FALSE)
head(out)
dadaFs <- system.file(respath, "dadaFs.rds")
dadaRs <- system.file(respath, "dadaRs.rds")
dadaRs[[1]]

if (!file.exists("final_table_taxa_si_sp.txt")) {
  dadaFs <- readRDS("./dadaFs.rds")
  dadaRs <- readRDS("./dadaRs.rds")
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
  head(mergers[[1]])
  saveRDS(mergers, "./mergers.rds")
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)
  write.table(seqtab, "./seqtab.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(seqtab, "./seqtab.rds")
  table(nchar(getSequences(seqtab)))
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
  dim(seqtab.nochim)
  sum(seqtab.nochim) / sum(seqtab)
  write.table(seqtab.nochim, "./seqtab.nochim", sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(seqtab.nochim, "./seqtab.nochim.rds")
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "tabled", "nonchim")
  rownames(track) <- sample.names
  head(track)
  write.table(track, "./track_reads.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  taxa_si <- assignTaxonomy(seqtab.nochim, chopath, multithread = FALSE, outputBootstraps = TRUE)
  write.table(taxa_si, "./taxa_si.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(taxa_si, "./taxa_si.rds")
  taxa_si_sp <- addSpecies(taxa_si[[1]], sppath)
  write.table(taxa_si_sp, "./taxa_si_sp.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(taxa_si_sp, "./taxa_si_sp.rds")
  mock <- 1
  if (mock == 1) {
    setwd(rrpath)
    unqs.mock <- seqtab.nochim[Mock]
    unqs.mock <- sort(unqs.mock[unqs.mock > 0], decreasing = TRUE)
    cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
    mock.ref <- getSequences(Mock)
    match.ref <- as.integer(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
    cat("Of those,", match.ref, "were exact matches to the expected reference sequences.\n")
  }
  setwd(respath)
  uniquesToFasta(seqtab.nochim, "./rep_set.fna", ids = paste0("SV_", seq(length(getUniques(seqtab.nochim)))))
  setwd(inpath)
  samdf <- read.table("mapping.txt", header = TRUE, sep = "\t", row.names = 1)
  ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(samdf), tax_table(taxa_si_sp))
  t = t(otu_table(ps))
  d <- cbind(rownames(t), data.frame(t, row.names = NULL))
  colnames(d)[1] <- "#OTU ID"
  setwd(respath)
  write.table(d, "./ASV_table_ps.txt", sep = "\t")
  saveRDS(ps, "./ps.rds")
  final_table_taxa_si_sp <- phyloseq_to_df(ps) # you
  write.table(final_table_taxa_si_sp, "./final_table_taxa_si_sp.txt", sep = "\t", quote = F, row.names = FALSE)
  blank <- 1
  if (blank == 1) {
    setwd(respath)
    df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
    df$LibrarySize <- sample_sums(ps)
    df <- df[order(df$LibrarySize),]
    df$Index <- seq(nrow(df))
    sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
    contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg")
    table(contamdf.prev$contaminant)
    head(which(contamdf.prev$contaminant))
    contamdf.prev05 <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5) #Note that as before, the default threshold for a contaminant is that it reaches a probability of 0.1 in the statistical test being performed. In the prevalence test there is a special value worth knowing, threshold=0.5, that will identify as contaminants all sequences thare are more prevalent in negative controls than in positive samples.
    table(contamdf.prev05$contaminant)
    ps.pa <- transform_sample_counts(ps, function(abund) 1 * (abund > 0))
    ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
    ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
    df.pa <- data.frame(pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg),
                      contaminant = contamdf.prev$contaminant)
    plot <- ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
    ggsave(paste("./Prevalence_neg_controls_truesamples.pdf", sep = ""), plot)
    ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
    ps.noncontam
    saveRDS(ps.noncontam, "./ps.noncontam.rds")
    final_table_taxa_si_sp_noncontam <- phyloseq_to_df(ps.noncontam)
    write.table(final_table_taxa_si_sp_noncontam, "./final_table_taxa_si_sp_noncontam.txt", sep = "\t", quote = F, row.names = FALSE)
  }
  phylo <- 1
  if (phylo == 1) {
    setwd(respath)
    seqs <- getSequences(seqtab.nochim)
    names(seqs) <- seqs
    alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
    saveRDS(alignment, "./alignment.rds")
    phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm)
    fit = pml(treeNJ, data = phang.align)
    fitGTR <- update(fit, k = 4, inv = 0.2)
    saveRDS(fitGTR, "./fitGTR.rds")
    fitGTR2 <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
    saveRDS(fitGTR2, "./fitGTR2.rds")
  }
}
sessionInfo();
