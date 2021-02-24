# Age-dependent transcriptomic and proteomic responses to SARS-CoV-2 infection
RNA-seq and proteomic analysis of SARS-CoV-2 infected primary differentiated airway epithelial cultures in an age-dependent manner.

## Transcriptomic analysis

Scripts used to analyse the RNA-seq data is available in the **Transcriptomics** folder.

### Genome alignment

Reads were aligned to the human reference genome (GRCh38, GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna) combined with the SARS-CoV-2 assembly (NC_045512.2, GCF_009858895.2_ASM985889v3) with `subread-align` within the [Subread](http://bioinf.wehi.edu.au/subread/) (v2.0.1) package. Y-chromsosome-derived sequences were hard masked in the reference genome or alignments of female samples to ensure accuracy for RNA quantitation.

```
subread-align -t 0 --multiMapping -B 1 -i hg38.sarscov2_index -a hg38.sarscov2.refseq.gtf -r reads.fastq.gz -o subread_aligned.bam
```

For intron retention analysis, reads were aligned to the human reference genome with the splice-aware alignment tool [STAR](https://github.com/alexdobin/STAR).

```
STAR --runMode alignReads --genomeLoad  LoadAndKeep --readFilesCommand zcat --outSAMtype BAM Unsorted --genomeDir hg38_index --readFilesIn reads.fastq.gz
```

### Read counting

Aligned reads were counted at genomic regions using featureCounts within the `Subread` package.

```
featureCounts -O -M -s 2 -a hg38.sarscov2.refseq.gtf -o counts.txt subread_aligned.bam
```

### Differential gene expression analysis

Differential gene expression analysis was perfomed in `R 4.0` using [EdgeR](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html) (v3.30.3) using the generalized log-linear model (glm) method. Plots were created using [ggplot2](https://ggplot2.tidyverse.org/).

### Reference data

1. Human reference genome devoid of Alt-contigs was downloaded from [NCBI GenBank](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz)
1. Human RefSeq annotations were downloaded from [UCSC database](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)
1. SARS-CoV-2 assembly (GCF_009858895.2_ASM985889v3) downloaded from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz)
1. SARS-CoV-2 annotations (gtf file) was from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3)

## Proteomics analysis

Scripts used to analyse the RNA-seq data is available in the **Proteomics** folder.

### Mass spectrometry protein identification and quantification

Raw mass spectrometry peak lists was analysed using [MaxQuant](https://www.maxquant.org/) (v.1.6.2.10) against the human Swiss-Prot database [UP000005640](https://www.uniprot.org/proteomes/UP000005640) and SARS-CoV-2 Swiss-Prot database [UP000464024](https://www.uniprot.org/proteomes/UP000464024). 

Differential protein expression analysis was performed with the `R` package [`DEP`](https://github.com/arnesmits/DEP) (v.1.10.0). 

Functional analysis of differentially abundant proteins was performed with IPA. Graphs were plotted in `R` using `ggplot2`. 
