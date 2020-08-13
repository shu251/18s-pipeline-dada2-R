# 18s-pipeline-dada2-R

## Pipeline for running DADA2 in R specific for microbial eukaryotic work.

### (1) Perform initial QC & removal primers

For high-throughput sequence quality control and primer removal, I use snakemake.

```
# build conda environment to run snakemake
conda env create --name snake-qc --file envs/snake.yaml

# Enter this conda environment
conda activate snake-qc
```

**Set up working directory**
1. Modify your ```config.yaml``` file to tell snakemake where to look for raw sequences, the format of those sequences, and where you want the output trimmed reads. Use a text editor like _nano_ to modify this file. 
2. **Raw fastq sequences** should be


### (2) DADA2 pipeline in R



## Assign taxonomy to reference sequence

Extract reference sequences from qiime2 run:
```
# Enter qiime2 environment
qiime tools export --input-path ref-sequences.qza  
	--output-path ref-seqs/

```
Output will be a new directory called ```ref-seqs/```, inside is a .fna file with all reference sequences. Fasta headers are **Feature.ID**s from output ASV or OTU table.   


Now import those reads into R to assign taxonomy from dada2.
```
## R v3.6.1

# Read in reference sequence file, create dataframe with Feature.ID and reference sequence as columns

library(Biostrings)
fna_in <- readDNAStringSet("xxxx.fasta") #Import reference sequences
Feature.ID <- names(fna_in)
SEQUENCE <- paste(fna_in)
fna_df <- data.frame(Feature.ID, SEQUENCE)


# Assign taxonomy
library(dada2)

seqs <- as.character(fna_df$SEQUENCE) #extract sequences

# Assign taxonomy. Note that for the PR2 database in dada2, you need to use the taxLevels argument as is here (do not use the default)
taxa_pr2 <- assignTaxonomy(seqs, "/vortexfs1/omics/huber/shu/db/pr2-db/pr2_version_4.12.0_18S_dada2.fasta.gz", taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"), multithread = TRUE)

df_assigned$SEQUENCE <- row.names(df_assigned) 



