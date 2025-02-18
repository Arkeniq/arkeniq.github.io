---
title: "🧬 Sarcoma RNA-Seq DEG Run"
date: 2025-02-12 14:00:00 +0200
last_modified_at: 2025-02-12 23:26:00 +0200
description: "RNA-Seq differential gene expression analysis from sarcoma patients"
tags: [Drylab, Wetlab, RNA-Seq]
---

## Experimental Overview 🧪

RNA-Seq analysis was performed to identify differentially expressed genes (DEGs) between and 

## Wetlab Phase 💧

<dl>
<dt>1. RNA Extraction</dt>
<dd>Total RNA was extracted using QIAGEN's RNeasy Fibrous Tissue Kit from tumoral tissue preserved in RNAlater. After extraction, Invitrogen Qubit RNA High Sensitivity Kit was used to assess the RNA concentration. Yield: <b>X ng/µL</b>, RIN: <b>X.X</b>.</dd>

<dt>2. Library Preparation</dt>
<dd>mRNA enrichment via poly(A) selection, fragmentation, cDNA synthesis, and adapter ligation using [kit name]. Libraries were QC'd with Bioanalyzer.</dd>

<dt>3. Sequencing</dt>
<dd>Paired-end (2 × 150 bp) sequencing on an Illumina [platform], generating ~X million reads per sample.</dd>
</dl>

## Drylab Time 🧬

### 1. Quality Control & Trimming 🚀

QC and adapter trimming were done using FastQC and Trimmomatic.

```bash
fastqc *.fastq.gz
trimmomatic PE -threads 8 input_R1.fastq.gz input_R2.fastq.gz \
  output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
  output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
```

### 2. Alignment 📌

Reads were aligned to the reference genome ([Genome version]) using HISAT2.

```bash
hisat2 -p 8 -x genome_index -1 output_R1_paired.fastq.gz -2 output_R2_paired.fastq.gz \
  -S aligned_reads.sam
```

SAM files were converted to sorted BAM files using SAMtools:

```bash
samtools view -bS aligned_reads.sam | samtools sort -o aligned_reads_sorted.bam
samtools index aligned_reads_sorted.bam
```

### 3. Gene Quantification 📊

FeatureCounts was used to assign reads to genes:

```bash
featureCounts -T 8 -a annotation.gtf -o counts.txt aligned_reads_sorted.bam
```

### 4. Differential Expression Analysis 🔬

DESeq2 was used for identifying DEGs.

```r
library(DESeq2)
counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(condition=c("A", "A", "B", "B"))
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, file="DEGs.csv")
```

### 5. Functional Enrichment Analysis 🔎

Gene Ontology (GO) and KEGG pathway enrichment were performed using g:Profiler.

```r
library(clusterProfiler)
library(org.Hs.eg.db)
enriched <- enrichGO(gene=rownames(res[res$padj < 0.05, ]), OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH")
```

## Results & Insights 📌

- **X** genes were differentially expressed (<i>padj</i> < 0.05).
- Upregulated pathways: [Pathway X, Pathway Y]
- Downregulated pathways: [Pathway Z]
- [Biological interpretation and next steps]

---

Next up: Validation and integration with other omics data! 🚀
