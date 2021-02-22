# Age-dependent transcriptomic and proteomic responses to SARS-CoV-2 infection
RNA-seq and proteomic analysis of SARS-CoV-2 infected primary differentiated airway epithelial cultures in an age-dependent manner.

## Genome alignment

Reads were aligned to the human reference genome (GRCh38, GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna) combined with the SARS-CoV-2 assembly (NC_045512.2, GCF_009858895.2_ASM985889v3) with [subread-align](http://bioinf.wehi.edu.au/subread/) within the `Subread` package. Y-chromsosome-derived sequences were hard masked in the reference genome or alignments of female samples to ensure accuracy for RNA quantitation.

```
subread-align -t 0 --multiMapping -B 1 -i hg38.sarscov2_index -a hg38.sarscov2.refseq.gtf -r reads.fastq.gz -o subread_aligned.bam
```

## Read counting

Aligned reads were counted at genomic regions using featureCounts within the `Subread` package.

```
featureCounts -O -M -s 2 -a hg38.sarscov2.refseq.gtf -o counts.txt subread_aligned.bam
```

## Differential gene expression analysis

Differential gene expression analysis was perfomed in `R 4.0` using [EdgeR](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html) (v3.30.3) using the generalized log-linear model (glm) method. Plots were created using [ggplot2](https://ggplot2.tidyverse.org/) and [plotly](https://plotly.com/r/).

## Reference data

1. Human reference genome devoid of Alt-contigs was downloaded from [NCBI GenBank](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz)
1. Human RefSeq annotations were downloaded from [UCSC database](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)
1. SARS-CoV-2 assembly (GCF_009858895.2_ASM985889v3) downloaded from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz)
1. SARS-CoV-2 annotations (gtf file) was from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3)

