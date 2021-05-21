################################################################################
# Cutadapt for DADA2 pipeline - ITS                                            #
#                                                                              #
#                                                                              #                          
# This script removes the FWD or REV primer from each read and also the REVCOMP#
# version of its corresponding primer in case the read extends long enough     #
#                                                                              #
# V1.0                                                                         #
# Roey Angel - May 2019                                                        #
# roey.angel@bc.cas.cz                                                         #
# Based on https://benjjneb.github.io/dada2/ITS_workflow.html                  #
################################################################################


############################################################
# TO DO                                                    #
############################################################

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
  "stringr" # Simple, Consistent Wrappers for Common String Operations
  # "tidyverse", # for dplyr forcats ggplot2 readr tibble
  # "broom", # Convert Statistical Analysis Objects into Tidy Data Frames
  # "magrittr", # pipes
  # "scales", # Generic plot scaling methods
  # "svglite", # for svg files
  # "seqinr" # Biological Sequences Retrieval and Analysis 
) 

.inst <- .cran_libs %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_libs[!.inst], 
                   lib = userLocation, 
                   repos = repo)
}

.bioc_libs <- c(
  "dada2", # sample inference from amplicon sequencing data
  "ShortRead", # FASTQ input and manipulation 
  "Biostrings" # Efficient manipulation of biological strings
)

# devtools::install_github("benjjneb/dada2", lib = userLocation, force = TRUE)

.inst <- .bioc_libs %in% installed.packages()
if (any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(ask = F, lib = userLocation) # upgrade bioC packages
  biocLite(.bioc_libs[!.inst], ask = F, lib = userLocation, lib.loc = userLocation)
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

# Functions ---------------------------------------------------------------

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Set general parameters ----------------------------
# grab supplied argument
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  message("No args supplied; using default values")
  (path <- "./") # directory containing the folder with the fastq files
  (min_read_len <- 50)
  (FWD <- "CTTGGTCATTTAGAGGAAGTAA") # forward primer ITS1f
  (REV <- "GCTGCGTTCTTCATCGATGC") # reverse primer ITS2
} else if (length(args) == 1) {
  message("Data path is supplied")
  (path <- args[1]) # directory containing the folder with the fastq files
  (min_read_len <- 50)
  (FWD <- "CTTGGTCATTTAGAGGAAGTAA") # forward primer ITS1f
  (REV <- "GCTGCGTTCTTCATCGATGC") # reverse primer ITS2
} else if (length(args) == 2) {
  message("Data path is supplied")
  (path <- args[1]) # directory containing the folder with the fastq files
  (min_read_len <- args[2])
  (FWD <- "CTTGGTCATTTAGAGGAAGTAA") # forward primer ITS1f
  (REV <- "GCTGCGTTCTTCATCGATGC") # reverse primer ITS2
} else if (length(args) == 3) {
  stop("Only 1, 2 or 4 arguments are allowed.\n", call. = FALSE)
} else if (length(args) == 4) {
  message("Data path and primers are supplied")
  (path <-  args[1]) # directory containing the folder with the fastq files
  (min_read_len <- args[2])
  (FWD <- args[3]) # forward primer
  (REV <- args[4]) # reverse primer
} else {
  stop("Only up to 3 arguments are allowed.\n", call. = FALSE)
}

data_path <- paste0(path, "data_files")
fastqs <- list.files(data_path, pattern = "*.fastq.gz$", full.names = TRUE) # list gz fastq files only
if (isEmpty(fastqs)) { # make sure files are there
  message(
    "Error: could not find any fastq.gz files in ",
    data_path
  )
  break()
}

fastqs <-
  sort(fastqs) # Sort ensures forward/reverse reads are in same order
fastq1s <-
  fastqs[grepl(".*_R1_001.fastq.gz|\\.1\\.fastq", fastqs)] # Just the read 1 files
fastq2s <-
  fastqs[grepl(".*_R2_001.fastq.gz|\\.2\\.fastq", fastqs)]# Just the read 2 files

# Generate all possible primer orientations
(FWD.orients <- allOrients(FWD))
(REV.orients <- allOrients(REV))

# Filter out N-containing seqs, place rest in the filtN folder
fastq1s.filtN <- file.path(data_path, "filtN", basename(fastq1s)) # Put N-filterd files in filtN/ subdirectory
fastq2s.filtN <- file.path(data_path, "filtN", basename(fastq2s))
filterAndTrim(fastq1s, fastq1s.filtN, fastq2s, fastq2s.filtN, maxN = 0, multithread = TRUE)

# Count the number of times the primers appear in the forward and reverse read, 

rbind(FWD.1st.Reads = sapply(FWD.orients, primerHits, fn = fastq1s.filtN[[2]]), 
      FWD.2nd.Reads = sapply(FWD.orients, primerHits, fn = fastq2s.filtN[[2]]), 
      REV.1st.Reads = sapply(REV.orients, primerHits, fn = fastq1s.filtN[[2]]), 
      REV.2ndReads = sapply(REV.orients, primerHits, fn = fastq2s.filtN[[2]]))

# Remove Primers
cutadapt <- "~/.local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "noPrimers")
if (!dir.exists(path.cut)) dir.create(path.cut)
fastq1s.cut <- file.path(path.cut, 
                         str_replace(basename(fastq1s), ".fastq.gz", "_noPrimers.fastq.gz"))
fastq2s.cut <- file.path(path.cut, 
                         str_replace(basename(fastq2s), ".fastq.gz", "_noPrimers.fastq.gz"))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for (i in seq_along(fastq1s)) {
  system2(
    cutadapt,
    args = c(
        "-j", 0,
        R1.flags,
        R2.flags,
        "--minimum-length", min_read_len,
        # "--untrimmed-output untrimmed1.fa",
        # "--untrimmed-paired-output untrimmed2.fa",
        "--discard-untrimmed",
        "--pair-filter=any",
        "-n", 2, # -n 2 required to remove FWD and REV from reads
        "-o", fastq1s.cut[i], # output files
        "-p", fastq2s.cut[i], 
        fastq1s.filtN[i], # input files
        fastq2s.filtN[i]
    )
  ) # input files
}

# As a sanity check, we will count the presence of primers in the first cutadapt-ed sample:
  
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fastq1s.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fastq2s.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fastq1s.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fastq2s.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cut1s <- sort(list.files(path.cut, pattern = ".*_R1_001.fastq.gz", full.names = TRUE))
cut2s <- sort(list.files(path.cut, pattern = ".*_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cut1s, get.sample.name))
print(sample.names)
