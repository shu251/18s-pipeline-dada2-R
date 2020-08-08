# 18s-pipeline-dada2-R

Pipeline for running DADA2 in R specific for microbial eukaryotic work.



## Perform initial QC & removal primers

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
