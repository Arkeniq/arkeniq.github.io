---
title: "🧬 Sarcoma RNA-Seq DEG Run"
date: 2024-07-12 14:00:00 +0200
last_modified_at: 2024-07-15 23:26:00 +0200
description: "RNA-Seq differential gene expression analysis from sarcoma patients"
tags: [Drylab, Wetlab, RNA-Seq]
---

## Experimental Overview 🧪

RNA-Seq analysis was performed to identify differentially expressed genes (DEGs) between 12 sarcoma tissue samples and 2 controls. The amount of physical samples was limited, especially the amount of control samples we had. To add to that problem, the two control samples we had available were different in tissue type too - one was mostly adipose and the other muscular tissue. Decided to give it a go and try my best to make some sort of DGE analysis out of the available data, just to learn some of that RNA-seq goodness that everybody is bragging about. Here goes...

## Wetlab Phase 💧

<dl>
<dt>1. RNA Extraction</dt>
<dd>Total RNA was extracted using QIAGEN's RNeasy Fibrous Tissue Kit from 40-60 mg of tumoral tissue preserved in RNAlater. After extraction, Invitrogen Qubit RNA High Sensitivity Kit was used to assess the RNA concentration. Yield approx. 20-60 ng/µL</b>.</dd>

<dt>2. Library Preparation</dt>
<dd>Library preparation was performed using TruSight RNA Pan-Cancer Set, guided by the manufacturer's protocol, and loading up the cartridge into our Illumina miniSeq Sequencing System. Final goal was to be able to use the data analysis offered by Illumina's BaseSpace cloud-based environment, but I wanted to also try my hand at some manual DGE.</dd>
</dl>

## Drylab Time 🧬

### 1. Quality control & Trimming 

QC and adapter trimming were done using FastQC and Trimmomatic. Quality score was pretty good, the last basepair of each FASTQ was pretty low so I decided to crop it out, followed by a MAXINFO which is "an adaptive quality trimmer which balances read length and error rate to maximise the value of each read".  

```bash
fastqc *.fastq.gz

for s in {1..8}; do
for r in {1..2}; do
java -jar trimmomatic/trimmomatic-0.39.jar PE -threads 48 \
-basein P${s}_S${s}_L001_R${r}_001.fastq.gz \
-baseout P${s}_${r}.fastq.gz \
CROP:75; done; done

for s in {1..8}; do
for r in {1..2}; do
java -jar trimmomatic/trimmomatic-0.39.jar PE -threads 48 \
-basein P${s}_S${s}_L001_R${r}_001.fastq.gz \
-baseout P${s}_${r}.fastq.gz \
MAXINFO:75:0.5; done; done

```

### 2. Alignment

Reads were aligned to the reference genome GRCh38/hg38 using HISAT2. I tried following the differential expression analysis instructions found in the [StringTie manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) but got stopped in my tracks after reading that Ballgown was not really designed for *gene*-level differential expression analysis -- it was written specifically to do *isoform*-level DE. Using DESeq2 with FeatureCounts is a much better-supported operation if your main interests are in gene-level DE. Count-based models like those in DESeq2 are totally appropriate for gene-level DE (whereas Stringtie and Ballgown are tools for when the count-based models *don't* work). So Hisat2 -> FeatureCounts -> DESeq2 would be a great workflow if you don't need isoform-level results. (Thanks to @Alyssa Frazee on Bioconductor Support for this info)

```bash
for s in {1..8}; do 
hisat2 -q -p 26 --rna-strandness RF \
-x grch38_snp_tran/genome_snp_tran -1 Fastq/P${s}_1.fastq.gz \
-2 Fastq/P${s}_2.fastq.gz --dta | samtools sort -o P${s}.bam
done
```

### 3. Gene Quantification 📊

FeatureCounts was used to assign reads to genes:

```bash
featureCounts -S 2 -p -T 48 -a Homo_sapiens.GRCh38.111.gtf -o FeatureCountsFinal.txt *.bam
```

```stl
solid MyCoolModel
  facet normal 0.0 -1.0 0.0
    outer loop
      vertex 0.0 0.0 0.0
      vertex 1.0 0.0 0.0
      vertex 0.0 0.0 1.0
    endloop
  endfacet
endsolid MyCoolModel
```

### 4. Differential Expression Analysis

For identifying DEGs, I used the popular R package DESeq2. This tool is excellent for normalizing count data and performing statistical analysis, especially with small sample sizes, though the heterogeneity of my control samples was a known challenge.

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
