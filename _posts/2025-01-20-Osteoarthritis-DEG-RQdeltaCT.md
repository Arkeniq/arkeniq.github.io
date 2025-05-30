---
title: "🔬 Differential gene expression in osteoarthritis using qPCR & RQdeltaCT"
date: 2025-01-20 11:00:00 +0200
last_modified_at: 2025-05-15 21:41:00 +0200
description: "A walkthrough of relative quantification of gene expression in an osteoarthritis dataset using the RQdeltaCT R package"
tags: [Bioinformatics, R, qPCR, Osteoarthritis, RQdeltaCT]
---

A colleague was working on an interesting osteoarthritis project and had a set of qPCR (Quantitative Polymerase Chain Reaction) data measuring Ct values for several genes in osteoarthritis and control samples. They were looking to perform relative quantification of gene expression, and I jumped in to help with the data analysis using R\! This post walks through the main steps we took, primarily leveraging a neat R package called [RQdeltaCT](https://github.com/Donadelnal/RQdeltaCT) which is designed to streamline this kind of analysis.

## 🛠️ Setting up the environment 

First things first, I needed to get R ready with the necessary packages and load the data. I started by ensuring all the required packages were installed. The star of the show was `RQdeltaCT`, which I installed directly from GitHub. We also used `tidyverse` for general data manipulation. A handful of genes were amplified multiple times using the Bio-Rad CFX Real-Time PCR system, from 28 osteoarthritis and 3 control samples, and an average Ct was calculated. My colleague, Mădălina, merged all the raw Ct data into a CSV file named `data_Ct_long.csv` which was loaded using the `read_Ct_long()` function after being molded into the required format.

```r
# Set the working directory 
setwd("S:/BioMol/RQdeltaCT")
install.packages("remotes") 
remotes::install_github("Donadelnal/RQdeltaCT")
library(RQdeltaCT)
library(tidyverse)
path <- "data_Ct_long.csv"

# Load the Ct data using the RQdeltaCT function
data.Ct <- read_Ct_long(path = path,
                        sep = ",",           # CSV separator
                        dec = ".",           # Decimal character
                        skip = 0,            # Lines to skip at the beginning
                        add.column.Flag = TRUE, # Add a flag column if not present
                        column.Sample = 2,   # Column index for Sample names
                        column.Gene = 3,     # Column index for Gene names
                        column.Ct = 4,       # Column index for Ct values
                        column.Group = 1,    # Column index for Group information
                        column.Flag = 5)     # Column index for Flag information

# Check the structure of the loaded data
str(data.Ct)
```

## ⚙️ Quality control & Data preparation

Before diving into quantification, it's crucial to perform some quality control on the Ct values. The RQdeltaCT package provides functions to quickly visualize the raw Ct values, which can help spot outliers or issues. I generated barplots for Ct values per sample and per gene. These also save as TIFF files if `save.to.tiff = TRUE`.

```r
# Barplot of Ct values per sample
sample.Ct.control <- control_Ct_barplot_sample(data = data.Ct,
                                             flag.Ct = "Undetermined", # How undetermined values are flagged
                                             maxCt = 38,               # Max Ct value considered reliable
                                             save.to.tiff = TRUE)

# Barplot of Ct values per gene
gene.Ct.control <- control_Ct_barplot_gene(data = data.Ct,
                                           flag.Ct = "Undetermined",
                                           maxCt = 38,
                                           save.to.tiff = TRUE)
```

The plots were showing reliable enough results to keep going forward with all the samples/genes.

---

## Heatmap attempt and filtering

I then tried to generate a heatmap of the raw Ct values using the `control_heatmap()` function. This initially failed with a "breaks are not unique" error. The RQdeltaCT manual notes: "The control_heatmap() works only if various numbers of replicates are in the data, otherwise error will appear because the inherited pheatmap() function can not deal with the situation where all values are equal." This suggested the raw data might not have enough variability for this specific heatmap function as initially configured.

```r
library(pheatmap)
colors <- c("#4575B4","#FFFFBF","#C32B23")
control_heatmap(data.Ct,
                sel.Gene = "all",
                colors = colors,
                show.colnames = TRUE,
                show.rownames = TRUE, 
                fontsize = 9,
                fontsize.row = 9,
                angle.col = 45,
save.to.tiff = TRUE)
```

Next, I looked for samples or genes with a high fraction of unreliable Ct values (e.g., "Undetermined" or >38).

```r
# Finding samples with >50% unreliable Ct values
low.quality.samples <- filter(sample.Ct.control[[2]], Not.reliable.fraction > 0.5)$Sample
low.quality.samples <- as.vector(low.quality.samples)
print(paste("Low quality samples to remove:", low.quality.samples))

# Finding genes with >50% unreliable Ct values in any given group
low.quality.genes <- filter(gene.Ct.control[[2]], Not.reliable.fraction > 0.5)$Gene
low.quality.genes <- unique(as.vector(low.quality.genes))
print(paste("Low quality genes to remove:", low.quality.genes))
```

In our case, both `low.quality.samples` and `low.quality.genes` lists were empty, meaning all samples and genes passed this particular QC check.

We then explored filtering the data using `filter_Ct` primarily to handle values flagged as "Undetermined" or those exceeding maxCt.

```r
data.CtF <- filter_Ct(data = data.Ct,
                      flag.Ct = "Undetermined", # String indicating undetermined Ct
                      maxCt = 38,               # Ct values above this are considered unreliable
                      flag = c("Undetermined"), # Other flags to consider unreliable
                      remove.Gene = low.quality.genes, # From previous step (was empty)
                      remove.Sample = low.quality.samples) # From previous step (was empty)

dim(data.Ct)
dim(data.CtF)
```

The dimensions were both 203, so no entries were removed. Great! Onwards! Next step was preparing data for ddCt analysis, using the `make_Ct_ready` function.

```r
data.Ct.ready <- make_Ct_ready(data = data.CtF,
                                imput.by.mean.within.groups = FALSE)
```

## Finding a reliable reference gene

For the ΔΔCt method, a stable reference (housekeeping) gene is essential for normalization. We used the ctrlGene package to help identify the best one from our candidates ("ACTB1", "GUSB2").

```r
library(ctrlGene) # For finding reference genes
ref_gene_analysis <- find_ref_gene(data = data.CtF.ready, 
                                   groups = c("Osteoarthritis","Control"),
                                   candidates = c("ACTB1", "GUSB2"),
                                   norm.finder.score = TRUE,
                                   genorm.score = TRUE,
                                   save.to.tiff = TRUE)
print(ref_gene_analysis[[2]]) # Show stability scores table
```

The analysis (including NormFinder and geNorm scores) suggested ACTB1 had good characteristics (like low variance) to be a reliable reference gene for this dataset. So, ACTB1 became our reference, and I had to use it to normalise the rest.

```r
data.dCt.exp <- delta_Ct(data = data.CtF.ready,
                         normalise = TRUE,
                         ref = "ACTB1",
                         transform = TRUE)

data.dCt <- delta_Ct(data = data.CtF.ready,
                     normalise = TRUE,
                     ref = "ACTB1",
                     transform = FALSE)
```

