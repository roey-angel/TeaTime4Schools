################################################################################
# DADA2 pipeline - ITS                                                         #
# This implementation of the DADA2 pipeline merges 2 read files per sample:    #
# *_R1_* and *_R2_*                                          .                 #
# Optimised for both MiSeq and MiniSeq runs                                    #
# Usage: preferebly run on a server using run_DADA2.sh                         #
#                                                                              #
# V8.2                                                                         #
# Roey Angel - July 2019                                                   #
# roey.angel@bc.cas.cz                                                         #
# Based on http://benjjneb.github.io/dada2/tutorial.html                       #
################################################################################


############################################################
# TO DO                                                    #
# 1. Generate meaningful names for samples                 #
# 2. Remove samples that dont classify to the class level  #
# 3. Output code lines to log                              #
############################################################

#!diagnostics suppress=biocLite
  
Sys.setenv(R_LIBS_USER = "~/R/library") #
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())) # Uncomment if you have no write access to R path

repo <- "http://cran.wu.ac.at"
userLocation <- Sys.getenv("R_LIBS_USER")

update.packages(
  userLocation, 
  repos = repo,
  ask = FALSE
)

.cran_libs <- c(
  "tidyverse", # for dplyr forcats ggplot2 readr tibble
  "broom", # Convert Statistical Analysis Objects into Tidy Data Frames
  "magrittr", # pipes
  "scales", # Generic plot scaling methods
  "svglite", # for svg files
  "seqinr" # Biological Sequences Retrieval and Analysis 
) 

.inst <- .cran_libs %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_libs[!.inst], 
                   lib = userLocation, 
                   repos = repo)
}

.bioc_libs <- c(
  "Biostrings", # Efficient manipulation of biological strings
  "dada2", # sample inference from amplicon sequencing data
  "ShortRead", # FASTQ input and manipulation 
  "microseq", # Basic Biological Sequence Analysis
  "microcontax", # consensus taxonomy for prokaryotes
  "DECIPHER" # Tools for curating, analyzing, and manipulating biological sequences 
)

# devtools::install_github("benjjneb/dada2", lib = userLocation, force = TRUE)

.bioc_inst <- .bioc_libs %in% installed.packages()
if (any(!.bioc_inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(ask = F, lib = lib.loc)  # upgrade bioC packages
  BiocManager::install(.bioc_libs[!.bioc_inst], ask = F, lib = lib.loc)
}

# Load packages into session, and print package version
loaded_libs <-
  sapply(c(.cran_libs, .bioc_libs), require, character.only = TRUE)
if (!all(loaded_libs)) {
  print(paste0(
    "Error: package(s) ",
    paste(names(loaded_libs[loaded_libs == FALSE]), collapse = ", "),
    " failed to load"
  ))
  break()
}
sapply(c(.cran_libs, .bioc_libs), packageVersion)

# Set general parameters ----------------------------
# grab supplied argument
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  paste("No args supplied; using default values")
  (pooling <- FALSE) # defualt - single sample analysis
  (data_path <- "../") # directory containing the zipped fastq files after untaring
  (db_path <- "~/DADA2/Resources/")
} else if (length(args) == 1) {
  paste("Pooling is supplied")
  (pooling <- args[1])
  (data_path <- "../") # directory containing the zipped fastq files after untaring
  (db_path <- "../../../../../Resources/")
} else if (length(args) == 2) {
  paste("Pooling and data path args supplied")
  (pooling <- args[1])
  (data_path <- args[2])
  (db_path <- "../../../../../Resources/")
} else if (length(args) == 3) {
  paste("Pooling and data path args supplied")
  (pooling <- args[1])
  (data_path <- args[2])
  (db_path <- args[3])
} else {
  stop("Only up to 3 arguments are allowed.\n", call. = FALSE)
}

# list.files(db_path)
unite.ref <- "sh_general_release_dynamic_02.02.2019.fasta"
# mock_ref_seqs <- "ATCC_MSA-2003.fa"
# theme_set(theme_bw(base_size = 14, base_family = "DejaVu Sans"))
#if (!file_test("-d", "Results"))
#  dir.create("Results") # generate a folder for DADA2 results
amp_size <- 350
#truncate <- c(200, 160) # set triming points for read 1 and read 2 NOT NEEDED FOR MINISEQ

# Filtering and Trimming ----------------------------------------------------------------------
fastqs <- list.files(data_path, pattern = "*.fastq.gz$", full.names = TRUE) # list gz fastq files only
if (isEmpty(fastqs)) { # make sure files are there
  print(paste0(
    "Error: could not find any fastq.gz files in ",
    data_path
  ))
  break()
}

fastqs <-
  sort(fastqs) # Sort ensures forward/reverse reads are in same order
fastq1s <-
  fastqs[grepl(".*_R1_001_noPrimers.fastq.gz|\\.1\\.fastq", fastqs)] # Just the read 1 files
fastq2s <-
  fastqs[grepl(".*_R2_001_noPrimers.fastq.gz|\\.2\\.fastq", fastqs)]# Just the read 2 files

# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
fastq1s %>% 
  str_replace(., file.path(data_path, "/*(.*)\\.[0-9]\\.fastq\\.gz$"), "\\1") %>% 
  str_replace(., file.path(data_path, "/*(.*)_L001_R1_001_noPrimers.fastq.gz$"), "\\1") ->
  sample_names
print(sample_names)
if (isEmpty(sample_names)) {
  print(paste0(
    "Error: could not match sample name pattern in fastq.gz file names in ",
    data_path
  ))
  break()
}

pdf("QualProf.pdf")
plotQualityProfile(fastq1s[1:2])
plotQualityProfile(fastq2s[1:2])
dev.off()

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(data_path, "filtered")
if (!file_test("-d", filt_path))
  dir.create(filt_path)
if (length(fastq1s) != length(fastq2s))
  stop("Forward and reverse files do not match.")
filt1s <-
  file.path(filt_path, paste0(sample_names, "_R1_filt.fastq.gz"))
filt2s <-
  file.path(filt_path, paste0(sample_names, "_R2_filt.fastq.gz"))
# Filter
filter_output <- filterAndTrim(
  fastq1s,
  filt1s,
  fastq2s,
  filt2s,
  # truncLen = truncate, # 
  minLen = 50, # redundatnt with cutadapt, but just in case.
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
filter_output <- bind_cols(sample_name = sample_names, file_name = rownames(filter_output), as.data.frame(filter_output)) 
print(filter_output)

# Scan filter path for files that passed filtering
final_filt1s <- list.files(file.path(data_path, "filtered"),  pattern = "*_R1_filt.fastq.gz", full.names = TRUE)
final_filt2s <- list.files(file.path(data_path, "filtered"), pattern = "*_R2_filt.fastq.gz", full.names = TRUE)
final_filt1s %>% 
  sub("_R1_filt.fastq.gz$", "", .) %>% 
  sub((file.path(data_path, "filtered/")), "", .) ->
  filt_sample_names 

filtered_out <- sample_names[!(sample_names %in% filt_sample_names)]
if (length(filtered_out) == 0) {
  print("No sample was filtered out")
} else if (length(filtered_out) == 1) {
  print(paste0("Sample: ", filtered_out, " was filtered out"))
  filter_output %<>% filter(!sample_name %in% filtered_out)
} else if (length(filtered_out) > 1) {
  print(paste0("Samples: ", paste(filtered_out, collapse = ", "), " were filtered out"))
  filter_output %<>% filter(!sample_name %in% filtered_out)
}

# Learn the Error Rates -----------------------------------------------------------------------
err1 <- learnErrors(final_filt1s, randomize = TRUE, multithread = TRUE)
err2 <- learnErrors(final_filt2s, randomize = TRUE, multithread = TRUE)

ggsave("Read1Errors.pdf", plotErrors(err1, nominalQ = TRUE))
ggsave("Read2Errors.pdf", plotErrors(err2, nominalQ = TRUE))

# Dereplication -------------------------------------------------------------------------------
derep1s <- derepFastq(final_filt1s, verbose = TRUE)
derep2s <- derepFastq(final_filt2s, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derep1s) <- filt_sample_names
names(derep2s) <- filt_sample_names

# Sample Inference ----------------------------------------------------------------------------
dada1s <- dada(derep1s, err = err1, pool = pooling, multithread = TRUE)
dada2s <- dada(derep2s, err = err2, pool = pooling, multithread = TRUE)
dada1s[[1]]

# Merge paired reads --------------------------------------------------------------------------
mergers <-
  mergePairs(dada1s, derep1s, dada2s, derep2s, minOverlap = 6, verbose = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# reverseComplement rev-primer reads (if analysing a 4 files per sample MiSeq run (i.e. DOME method))
#for (rev_sample in grep("\\.R", names(mergers))) {
#  mergers[[rev_sample]]$sequence <-
#    microseq::reverseComplement(mergers[[rev_sample]]$sequence)
#}

# Construct sequence table and remove chimeras ------------------------------------------------
seqtab <- makeSequenceTable(mergers)
#seqtab <- makeSequenceTable(mergers[grep("^Zymo", names(mergers), orderBy = TRUE)])

# Table dimensions
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Now merge fwd-primer and rev-primer samples (if analysing a 4 files per sample MiSeq run (i.e. DOME method))
# samples2merge <- gsub('\\.[FR]$', '', rownames(seqtab)) # remove FR to make identical names
# seqtab %>%
#   as.data.frame() %>%
#   mutate(., Sample = samples2merge) %>%
#   group_by(Sample) %>%
#   summarise_all(sum) %>% 
#   select(-matches("Sample")) %>% 
#   as.matrix ->
#   seqtab_merged
# rownames(seqtab_merged) <- unique(samples2merge)
# 
# duplicated_samples <-
#   duplicated(samples2merge) |
#   duplicated(samples2merge, fromLast = TRUE)
# 
# if (!all(duplicated_samples)) {
#   cat(
#     "The following sample(s) were missing either fwd or reverse primer reads:",
#     rownames(seqtab)[!duplicated_samples]
#   )
# }

# then change the downstream var from 'seqtab' to 'seqtab_merged'

# Remove sequences differing from excected lenght by more than 3 bp (not recommended for ITS)
# seqtab_clean <- seqtab[, nchar(colnames(seqtab)) %in% seq(amp_size - 3, amp_size + 3)]
seqtab_clean <- seqtab

# remove chimeras
seqtab_nochim <-
  removeBimeraDenovo(
    seqtab_clean,
    method = "consensus", 
    allowOneOff = TRUE,
    verbose = TRUE,
    multithread = TRUE
  )
dim(seqtab_nochim)
# What proportions passed through?
sum(seqtab_nochim) / sum(seqtab_clean)


# Save sequences and seq_table --------------------------------------------
# save ESVs to a fasta file
Seq_names <- paste0("Seq_", seq(1, ncol(seqtab_nochim)))
Seqs_nochim <- data.frame(
  Seq.name = Seq_names,
  Sequence = colnames(seqtab_nochim),
  stringsAsFactors = FALSE
)
write.fasta(
  as.list(Seqs_nochim$Sequence),
  Seqs_nochim$Seq.name,
  "DADA2.Seqs.fa",
  as.string = TRUE,
  nbchar = 1000
)

saveRDS(seqtab_nochim, "DADA2.seqtab_nochim.RDS") # save also an RDS object for easy merging

# colnames(seqtab_merged_nochim) %<>% 
# gsub(".*([mM]ock)", "\\1", .) # this is needed for uncross

# save seq table without sequences
seqtab_nochim %>% 
  t() %>% 
  as_tibble() %>% 
  #  bind_cols("Seq#" = Seq_names, .) 
  bind_cols("#OTU" = Seq_names, .) %>% # #OTU is needed for uncross
  write_delim(.,
              path = "DADA2.seqtab_nochim.tsv",
              delim = "\t",
              col_names = TRUE
  )

# Track reads through the pipeline ------------------------------------------------------------
# 
getN <- function(x) sum(getUniques(x))
track <- cbind(filter_output, 
                map_dbl(dada1s, getN), 
                map_dbl(mergers, getN), 
                rowSums(seqtab_clean), 
                rowSums(seqtab_nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dada1s, getN) with getN(dada1s)
colnames(track) <- c("sample_name", "file_name", "input", "filtered", "denoised", "merged", "tabled", "nonchim")
# track <- cbind(filt_sample_names, track)
head(track)

write_delim(track, 
            path = "DADA2.track_each.tsv",
            delim = "\t",
            col_names = TRUE)

# Assign taxonomy -----------------------------------------------------------------------------
# genus.species <- assignSpecies(seqs, "TrainingSets/silva_species_assignment_v132.fa.gz")
# cat(unname(head(genus.species)))

# Using assignTaxonomy - UNITE
taxa_unite <- assignTaxonomy(
  seqs = seqtab_nochim,
  refFasta = file.path(db_path, unite.ref),
  tryRC = TRUE,
  outputBootstraps = TRUE,
  multithread = TRUE
)

# save taxa table without sequences
taxa_unite %<>% 
  do.call(cbind, .) %>% 
  as.tibble() %>% 
  bind_cols("Seq#" = Seq_names, .) 

write_delim(taxa_unite, 
            path = "DADA2.taxa_unite.tsv",
            delim = "\t",
            col_names = TRUE)
