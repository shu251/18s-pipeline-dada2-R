# library("devtools")
# devtools::install_github("benjjneb/dada2", ref="v1.16", force = TRUE)

library(tidyverse); 
library(dada2)

# Trim raw fastq reads ----------------------------------------------------

path_filtered <- file.path(getwd(), "filtered-seqs")
r1_fastq_files <- list.files(pattern = "R1.fastq.gz", path = "raw-seqs", full.names = TRUE)
r2_fastq_files <- list.files(pattern = "R2.fastq.gz", path = "raw-seqs", full.names = TRUE)

# fastq_files
# file.path(getwd(),fastq_files)
# ?filterAndTrim()
# filterAndTrim(file.path(getwd(),fastq_files), path_filtered, 
              # truncLen=240, maxEE=1, truncQ=11, rm.phix=TRUE,
              # compress=TRUE, verbose=TRUE, multithread=TRUE)

# File parsing
# filtpathF <- "raw-seqs" # CHANGE ME to the directory containing your filtered forward fastqs
# filtpathR <- "raw-seqs" # CHANGE ME ...
# filtFs <- list.files(filtpathF, pattern="R1.fastq.gz", full.names = TRUE)
# filtRs <- list.files(filtpathR, pattern="R2.fastq.gz", full.names = TRUE)
# filtFs
# filtRs

r1_fastq_files <- list.files(pattern = "R1.fastq.gz", path = "raw-seqs", full.names = TRUE)
r2_fastq_files <- list.files(pattern = "R2.fastq.gz", path = "raw-seqs", full.names = TRUE)

# devtools::install_github("tidyverse/tidyverse")
# library(tidyverse)
# str_split(basename(filtFs), "(?=_R1.fastq.gz)")
# str_split(basename(filtRs), "(?=_R2.fastq.gz)")

sample.names <- sapply(str_split(basename(r1_fastq_files), "(?=_R1.fastq.gz)"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

sample.namesR <- sapply(str_split(basename(r2_fastq_files), "(?=_R2.fastq.gz)"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(r1_fastq_files) <- sample.names
names(r2_fastq_files) <- sample.names
sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(r1_fastq_files, nbases=100, multithread=TRUE)
# ?learnErrors
# Learn reverse error rates
errR <- learnErrors(r2_fastq_files, nbases=100, multithread=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
head(mergers)

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "path/to/run1/output/seqtab.rds") # CHANGE ME to where you want sequence table saved