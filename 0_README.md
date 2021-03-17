# BIO722
This is the BJE portion of the BIO722 graduate course in Advanced Bioinformatics.  This section will be divided into two 2-hour interactive lectures. The first will be designed around a group project on macaque genomics that uses reduced representation genome sequence data and a reference genome. The second will focus on applications of kmers in bioinformatics.

## Background
Many interesting organisms have large genomes, making complete genome sequencing expensive and assembly challenging, especially for multiple individuals.  A relatively cost-efficient solution has been developed called "reduced representation genome sequencing" (RRGS).  This approach enables deep sequencing of multiple genetic samples from the same (homologous) genomic regions.  It takes advantage of next generation (Illumina) sequencing technology and requires relatively simple laboratory preparation (DNA extraction), which can be accomplished with a centrifuge and a heat block.  Other laboratory steps (library construction) can be outsourced or done in house, depending on the equipment, funds, and expertise that are available. This approach has many applications, including phylogenomics, population genomics, linkage mapping, and analysis of gene flow.

## Goals: Section 1
The goal of this section of the course is to introduce students to some basic aspects of read mapping and genotyping. Ideally this will provide a sufficient level of exposure so that students will be able to figure out how to learn more, and then be able to generate and analyze their own datasets. Many of the approaches to analyze RRGS data are identical to analysis of whole genome sequencing (WGS) data (e.g., demultiplexing, trimming, mapping to a reference, read filtering, genotyping), so even if you work on things with small genomes where you can generate whole genome sequences, the approaches will still be applicable.  All of the approaches here involve use of a reference genome.  De novo assembly and analysis of RRGS and WGS data is also possible, but not covered here.

Together we will 
* examine an example Illumina dataset, including an explanation of fastq format, de-multiplexing of multiplexed samples, and quality trimming.  
* discuss reference genomes, what they are and why we use them if possible
* align our example dataset to a reference genomes
* discuss the vcf file format and some applications of these data to population genomics.  
* we also will discuss indel re-alignment, base recalibration and multisample genotype calling with GATK, and (due to time constraints) only perform genotyping with bcftools

## Group project
We will do a group project that is based on a RRGS dataset from the Tonkean macaque monkey (*Macaca tonkeana*).  Tonkean macaques inhabit the central Indonesian island of Sulawesi and, like other papionin monkeys, have a social system characterized by strong female philopatry and obligate male migration.  Reproductive success is thought to be more variable among males than females.  If this is true, we  expect  molecular polymorphism on the X chromosome to be elevated relative to an expectation with equal variance in reproductive success among the sexes (sounds complicated, but Ben will explain this). To test this hypothesis, each student will map RRGS data from the Tonkean macaque to one of the chromosomes of the reference genome of a closely related species – the rhesus macaque (*Macaca mulatta*), whose genome has been completely sequenced.  We will then evaluate molecular diversity of the X chromosome and the autosomes of the Tonkean macaque, and explore what this can tell us about the social system of these super-cool monkeys.  
 
## OK, lets begin by clicking [here](https://github.com/evansbenj/BIO722.md/blob/main/1_Evaluating_your_data.md).

