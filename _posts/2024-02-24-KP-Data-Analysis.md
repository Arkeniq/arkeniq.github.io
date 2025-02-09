---
title: "🧫 Sequencing a peculiar <i>K. pneumoniae</i> isolate"
date: 2024-02-24 14:14:00 +0200
last_modified_at: 2024-02-25 23:26:00 +0200
description: "ARGs found in a weird looking KP culture"
tags: [Drylab, Wetlab]
---


I've started a sequencing run from an atypical *Klebsiella pneumoniae* isolate that displayed unusual growth characteristics and high antibiotic resistance. By sequencing its genome, I wanted to identify the resistance genes responsible for its phenotype and understand their genomic context.

## Wetlab Phase 💧

<dl>
<dt>1. DNA Extraction</dt>
<dd>Bacterial DNA was extracted directly from the bacterial culture using QIAGEN’s QiAMP Microbiome Kit (manual extraction, bead beating as cell lysis method), yielding approximately <b>80 ng/µL</b>.</dd>

<dt>2. Library prep</dt>
<dd>SQK-LSK110 ligation sequencing kit, NEBNext Ultra II End repair module, NEBNext Quick Ligation Module and AMPure XP Beads were used for library preparation.</dd>

<dt>3. Sequencing</dt>
<dd>Loading the library on an R9.4.1 flowcell on the MinION Mk1B was straight forward. <b>13.6 fmoles</b> of library used. The flowcell used was being stored at 4°C for months after being previously used, but it still had a decent amount of active pores, at least on the initial check - 920 out of the maximum advertised amount of 2048.</dd>

<dt>4. Basecalling</dt>
<dd>The next day, after 19 hours of run time, we were left with approx. 580 pores and enough reads (612k passed) to stop the run. The real-time basecalling (HAC model) generated <b>4Gb</b>, with an approx. <b>N50 of 16Kb</b>. Estimated coverage of <b>~742x</b>.</dd>

</dl>

## Drylab Time 🧬

I wanted to test RRWick's new tool - [Trycycler](https://github.com/rrwick/Trycycler). It's designed to generate high-quality bacterial consensus assemblies from long-read sequencing data. Unlike traditional assemblers that rely on a single method, Trycycler takes a hybrid approach by combining multiple assemblies of the same dataset, leveraging their strengths to produce a more accurate and contiguous final assembly. It seemed like the perfect moment to test this method so I had to give it a go.

![Trycycler workflow](assets/images/trycyclerworkflow.png)
_Trycycler workflow. Credits: Trycycler Github_

1. Running **EPI2ME's WIMP**:
<br>
Just to make sure our reads are *K. pneumoniae*.
<br>
![WIMP](assets/images/WIMPKP.png){: w="300" .center }

2. **QC and filtering** using Filtlong: <br>
QC scores were good, we kept fragments longer than 1Kb to help with the assembly.

    ```bash
   zcat *.fastq.gz > KP9_sup_seq.fastq
   fastqc KP9_sup_seq.fastq
   filtlong –min_length 1000 –keep_percent 95 KP9_sup_seq.fastq > KP9_sup_filtrat.fastq
    ```

3. **Subsampling reads**: <br>
This step creates the assemblies.We used Flye, Miniasm/Minipolish and Raven.

   ```bash
   trycycler subsample --reads KP9_sup_filtrat.fastq --out_dir subsets --genome_size 5.5m --threads 48

   #Flye:
   for s in 01 04 07 10; do
   flye --nano-hq subsets/sample_”$s”.fastq --threads 52 --out_dir assemblies/assembly”$s”
   done

   #Miniasm/Minipolish:
   for s in 02 05 08 11; do
   minimap2 -x ava-ont -t 52 subsets/sample_”$s”.fastq subsets/sample_”$s”.fastq > overlaps.paf
   miniasm -f subsets/sample_”$s”.fastq overlaps.paf > assembly.gfa
   minipolish -t 52 subsets/sample_”$s”.fastq assembly.gfa > polished.gfa
   any2fasta polished.gfa > assemblies/assembly”$s”/assembly”$s”.fasta
   done

   #Raven:
   for s in 03 06 09 12; do
   raven --threads 52 --disable-checkpoints subsets/sample_”$s”.fastq > assemblies/assembly”$s”/assembly”$s”.fasta
   done
   ```


4. **Manual curation and clustering**:
<br>
I used [Bandage](https://rrwick.github.io/Bandage/) to check the contig circularity. It looked fine, with a clear circular chromosome and what appeared to be 2 plasmids. 
<br>
![Bandage](assets/images/bandage.png){: w="400" }
<br>
Clustering the assemblies was straight-forward using Trycycler.

   ```bash
   trycycler cluster --threads 48 --assemblies assemblies/*.fasta --reads KP9_sup_filtrat.fastq --out_dir cluster
   ```
<br>
I then looked at [FigTree](https://evomics.org/resources/software/molecular-evolution-software/figtree/) to identify out of place contigs. The clustering was succesful - 4 clearly defined clusters - a chromosome of ~5.3Mb and 3 smaller contigs of ~153Kb, ~56Kb and ~4Kb respectively, which were probably plasmids. I ended up removing the D contig in the 3rd cluster since its size was very different to the other ones.
<br>
![FigTree](assets/images/figtree.png){: w="300" .center }


5. **Alignment and consensus generation**:
<br>
Before generating the consensus sequence, the contigs have to be reconciled, aligned and partitioned. **Reconciling contigs** refers to the process of merging or resolving discrepancies between overlapping contigs, whereas the **multiple sequence alignment (MSA)** aligns the fragments against each other to identify regions of similarity/difference. **Partinioning** organizes the sequences into distinct regions which can be processed independently, improving accuracy. Contig reconciliation was the only step that needed manual processing, having me remove a few contigs that had low pairwise identities to one another. I only kept those with higher than 99% pairwise identity. The latter steps were as straightforward as they can get.
<br>
   ```bash
   #Ran these separately for all the clusters
   trycycler reconcile --threads 48 --reads KP9_sup_filtrat.fastq --cluster_dir clusters/cluster_00*
   trycycler msa –threads 48 –cluster_dir clusters/cluster_00*
   trycycler partition –threads 48 –reads KP9_sup_filtrat.fastq –cluster_dir clusters/cluster_00*
   trycycler consensus –threads 48 –cluster_dir clusters/cluster_00*
   ```

6. **Long and short read polishing**:
<br>
Using some short Illumina reads from our friends over at the Cantacuzino Institute I was able to take the clusters not only through **Medaka** for long-read polishing, but also **Polypolish**.
<br>
   ```bash
   #Medaka for all the clusters, renaming the files and cleaning the outputs
   for c in clusters/cluster_*; do
   medaka_consensus -i “$c”/4_reads.fastq -d “$c”/7_final_consensus.fasta -o “$c”/medaka -m r941_min_sup_g507 -t 48
   mv “$c”/medaka/consensus.fasta “$c”/8_medaka.fasta
   rm -r “$c”/medaka “$c”/*.fai “$c”/*.mmi
   done
   #Combining the clusters  into one single .FASTA
   cat clusters/cluster_*/8_medaka.fasta > KP9S20_consensus.fasta
   #Polypolish using shortreads
   bwa index KP9S20_consensus.fasta
   bwa mem -t 48 -a KP9S20_consensus.fasta shortreads/kp9_sr1.fastq.gz > shortreads/pp_align_1.sam
   bwa mem -t 48 -a KP9S20_consensus.fasta shortreads/kp9_sr2.fastq.gz > shortreads/pp_align_2.sam
   polypolish KP9S20_consensus.fasta shortreads/pp_align_1.sam shortreads/pp_align_2.sam > shortreads/polypolishkp9.fasta
   ```



## ARGs Detection 🦠

To identify ARGs I've used another awesome tool - [StarAMR](https://github.com/phac-nml/staramr) which is able to scan the .FASTA against ResFinder/PointFinder/PlasmidFinder DBs, making a snazzy report of known antimicrobial resistance genes it detects.

```bash
staramr search polypolishkp9.fasta --output-dir ./AMR
```

In the chromosome, StarAMR spotted *blaTEM-1b*, which is linked to resistance against piperacillin. The more notable ARGs were found in the plasmid - *qnrB4*, *blaOXA-1*, *aac(6')-Ib-cr* and *aph(3')-Ia*, which gives resistance to multiple antibiotics. The coolest (and a bit worrying) part? It found ***blaDHA-1***, which hasn’t been documented in my region before.
<br>

![ARGChrom](assets/images/ARGchrom.png){: w="800" }_ARGs in the chromosome_
![ARGPlasmid](assets/images/ARGplasmid.png){: w="800" }_ARGs in the smaller, IncR plasmid_

I also found this super cool and easy to use web interface that lets you visualise and annotate bacterial genomes, with interactive circule and linear genome maps - [Proksee](https://proksee.ca/). I initially used pyCirclize to make some circos but this is just super fun to use. I'll make a separate journal entry about it since the whole Python deal is interesting enough to have its own post. 

Anyway, here are the circos generated with Proksee, as well as a zoomed view of the ***blaDHA-1*** area of the plasmid, with some mobile genomic elements (MGEs) I had Proksee annotate for me. 

<div style="display: flex; justify-content: space-between;">
   <figure style="text-align: center;">
      <img src="assets/images/prokseechrom.png" height="300" alt="K. pneumoniae chromosome">
      <figcaption style="font-size: 0.8em; opacity: 0.6;">K. pneumoniae chromosome</figcaption>
   </figure>
   <figure style="text-align: center;">
      <img src="assets/images/prokseeplasmid.png" height="300" alt="K. pneumoniae smaller, IncR plasmid">
      <figcaption style="font-size: 0.8em; opacity: 0.6;">K. pneumoniae plasmid</figcaption>
   </figure>
</div>


## Conclusions 🤔

Trycyler worked great and was rather easy to use. The idea behind it is very solid too, improving accurancy over other tools like Unicycler. Learned a lot from this and found some cool ARGs. Time to dig deeper and try to understand where this *Klebsiella pneumoniae* isolate got its weird plasmid from. 

---

Next up: phylogenetic analysis and comparison with other sequenced strains. Stay tuned for more genomic adventures! 🔬✨
