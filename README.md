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

```
library(dada2)

## pending ##

```

***

## Assign taxonomy to reference sequences  

Use ```assignTaxonomy()``` from the **dada2** package to assign taxonomy. Below code can also be used to assign taxonomy to other sequences.   

### (1) Bring sequences for assignment into R
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
# Dataframe with 1 column of sequence header and a 2nd column with sequence.
```

### (2) Assign taxonomy

```
library(dada2)

seqs <- as.character(fna_df$SEQUENCE) #extract sequences

# Assign taxonomy. Note that for the PR2 database in dada2, you need to use the taxLevels argument as is here (do not use the default)
taxa_pr2 <- assignTaxonomy(seqs, "/vortexfs1/omics/huber/shu/db/pr2-db/pr2_version_4.12.0_18S_dada2.fasta.gz", 
	taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"), 
	multithread = TRUE, minBoot = 0, outputBootstraps = TRUE)

# See "Compile-taxonomy-assignments.ipynb" to compile with exisiting count table
```

### (3) Compile with original fasta file or count table
See R notebook ```Compile-taxonomy-assignments.ipynb``` to compile taxonomy assignments with previous count table. This specific example compared output from qiime2 taxonomy and dada2 approach.

```
library(Biostrings); library(tidyverse)
## ASVs
ref_asv <- readDNAStringSet("/vortexfs1/omics/huber/shu/slo-pier-weekly/qiime2/asv/slo-pier-ref-seqs-asv.fna")
Feature.ID <- names(ref_asv)
ReferenceSequence <- paste(ref_asv)
fna_df <- data.frame(Feature.ID, ReferenceSequence)

load("Pier-assigned-refseqs.RData", verbose = T)
asv_table <- read.delim("/vortexfs1/omics/huber/shu/slo-pier-weekly/qiime2/asv/CountTable-wtax-2020-04-22.txt")

# Comile ASV results with reference sequence
asv_wtax <- data.frame(taxa_pr2) %>% 
    rownames_to_column(var = "ReferenceSequence") %>% 
    right_join(fna_df) %>% 
    unite(Taxon_dada2_boot0, starts_with("tax."), sep = ";") %>% 
    unite(Confidence_dada2, starts_with("boot."), sep = ";") %>% 
    left_join(asv_table) %>% 
    select(Feature.ID, Taxon_qiime2 = Taxon, Conf_qiime2 = Confidence, 
            Taxon_dada2_boot0, Conf_dada2_boot0 = Confidence_dada2, everything()) %>% 
    data.frame
# head(asv_wtax)
```


### (4) Explore & reassign taxonomy from PR2

Using the list of taxa and bootstrap values for all taxonomic assignment levels, we need to determine what the most appropriate threshold is.  

In R (v3.6.1)
```
library(tidyverse)

# Import R object that has FeatureID information, reference sequence, and PR2 assignment information from above
head(asv_wtax)


# Create new column of updated taxonomic list, based on a minBoot threshold of 70
asv_updated_tax <- asv_wtax %>% 
  type.convert(as.is = TRUE) %>%
  separate(Confidence_dada2, c("Kingdom_boot","Supergroup_boot","Division_boot","Class_boot","Order_boot","Family_boot","Genus_boot","Species_boot"), sep = ";", convert = TRUE) %>%
  separate(Taxon_dada2_boot0, c("Kingdom_lev","Supergroup_lev","Division_lev","Class_lev","Order_lev","Family_lev","Genus_lev","Species_lev"), sep = ";", convert = TRUE) %>%
  mutate(Taxon_updated = 
     case_when(
      Species_boot >= 70 ~
         select(., ends_with("_lev")) %>% 
            reduce(str_c, sep=";"),
      Species_boot < 70 & Genus_boot >= 70 ~ 
          select(., Kingdom_lev:Genus_lev) %>%
            reduce(str_c, sep=";"),
      Genus_boot < 70 & Family_boot >=70 ~
          select(., Kingdom_lev:Family_lev) %>%
            reduce(str_c, sep=";"),
      Family_boot < 70 & Order_boot >=70 ~
          select(., Kingdom_lev:Order_lev) %>%
            reduce(str_c, sep=";"),
      Order_boot < 70 & Class_boot >=70 ~
          select(., Kingdom_lev:Class_lev) %>%
            reduce(str_c, sep=";"),
      Class_boot < 70 & Division_boot >=70 ~
          select(., Kingdom_lev:Division_lev) %>%
            reduce(str_c, sep=";"),
      Division_boot < 70 & Supergroup_boot >=70 ~
          select(., Kingdom_lev:Supergroup_lev) %>%
            reduce(str_c, sep=";"),
      TRUE ~ Kingdom_lev)) %>% 
  data.frame

```
Output generates a new column with an updated taxonomic assignment that is concatenated (;) by the appropriate minimum bootstrap value. Option to set different levels for various taxonomic levels.

_Last updated - 19-08-2020 - SKH_
