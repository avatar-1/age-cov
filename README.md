# cov-rnaseq
RNAseq of SARS-CoV-2 infected organoids

## Trimming

Raw reads were trimmed using `Trimmomatic` (v0.38).

```
java -jar trimmomatic-0.38.jar SE -phred33 reads.fastq.gz reads_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
```

## Genome alignment

Trimmed reads were aligned to the human reference genome (hg38) combined with the SARS-CoV-2 assembly (GCF_009858895.2_ASM985889v3) with [subread-align](http://bioinf.wehi.edu.au/subread/) within the `Subread` package. 

```
subread-buildindex -o hg38_CoV_index hg38_CoV.fna.gz
subread-align -t 0 -T 16 -i hg38_CoV_index -a hg38_CoV.gtf -r reads_trimmed.fastq.gz -o subread_aligned.bam
```

## Read counting

Aligned reads were counted at genomic regions using featureCounts within the `Subread` package.

```
featureCounts -T 16 -t exon -g gene_id -a hg38_CoV.gtf -o counts.txt subread_aligned.bam
```

## Differential gene expression analysis

Differential gene expression analysis was perfomed in `R 4.0` using [EdgeR](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html) (v3.30.3). Plots were created using [ggplot2](https://ggplot2.tidyverse.org/) and [plotly](https://plotly.com/r/).

## Functional enrichment analysis

Functional enrichment analysis of differentially expressed genes was perfomed using [gProfiler](https://biit.cs.ut.ee/gprofiler/) and [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp).

## Reference data

1. Human reference genome devoid of Alt-contigs was downloaded from [NCBI GenBank](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz)
1. Human RefSeq annotations were downloaded from [UCSC database](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)
1. SARS-CoV-2 assembly (GCF_009858895.2_ASM985889v3) downloaded from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz)
1. SARS-CoV-2 annotations (gtf file) was from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3)

