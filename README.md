# cov-rnaseq
RNAseq of SARS-CoV-2 infected organoids

## Trimming

Raw reads were trimmed using Trimmomatic (v0.38)

```
java -jar trimmomatic-0.38.jar SE -phred33 reads.fastq.gz reads_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
```

## Genome alignment

Trimmed reads were aligned to the human reference genome (hg38) combined with the SARS-CoV-2 assembly (GCF_009858895.2_ASM985889v3) using [subread-align](http://subread.sourceforge.net/).

```
subread-buildindex -o hg38_CoV_index hg38_CoV.fna.gz
subread-align -t 0 -T 16 -i hg38_CoV_index -r reads_trimmed.fastq.gz -o subread results.bam
```

## Reference data
1. Human reference genome devoid of Alt-contigs was downloaded from [NCBI GenBank](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz)
1. SARS-CoV-2 assembly (GCF_009858895.2_ASM985889v3) downloaded from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz)
1. SARS-CoV-2 annotations (gff and gtf file) was from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3)
