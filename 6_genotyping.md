(Or you can go back to automating the readmap alignment for multiple samples [here](https://github.com/evansbenj/BIO722.md/blob/main/4_automating_readmapping_for_multiple_samples.md).


# Genotyping with bcftools

There are many programs that will call genotypes from `.bam` files; the ones I usually use are [GATK](https://gatk.broadinstitute.org/hc/en-us) and [bcftools](http://samtools.github.io/bcftools/bcftools.html). `GATK` has a procedure called base quality score recalibration (`BQSR`) which allows you to alter quality scores to more accurately reflect the error rate, which can vary in sample-specific ways depending on genomic context and technical variation. This procedure is perfoermed prior to genotyping (it generates new `.bam` files that potentially have new quality scores for every base pair of read. We're not going to go over this today, but please feel free to ask me about this later if you want.

For the sake of time and because it is somewhat simplier, today we will use bcftools to genotype our bam files. We have everything we need to make a genotype (variant call format - vcf) file with all of the samples, including:
* sorted bam files for each sample
* an index for each bam file (.bai file)
* a genome (actually a chromosome in our case for one chromosome - chrZZZ.fa) and also supporting files (chrZZZ.amb  chrZZZ.fa.ann  chrZZZ.fa.bwt  chrZZZ.fa.fai  chrZZZ.fa.pac  chrZZZ.fa.sa).

Please enter this command (after changing the ZZZ to match your chromosomes):

```
bcftools mpileup -Ou -f ../my_monkey_chromosome/chrZZZ.fa PF515_chrZZZ_sorted_realigned_dedup.bam PM561_chrZZZ_sorted_realigned_dedup.bam PM565_chrZZZ_sorted_realigned_dedup.bam PM566_chrZZZ_sorted_realigned_dedup.bam PM567_chrZZZ_sorted_realigned_dedup.bam PM582_chrZZZ_sorted_realigned_dedup.bam PM584_chrZZZ_sorted_realigned_dedup.bam PM592_chrZZZ_sorted_realigned_dedup.bam PM602_chrZZZ_sorted_realigned_dedup.bam | bcftools call -mv -Oz -o allsamples_chrZZZ_merged_sorted_realigned_dedup.bam.vcf.gz

```

This should finish quickly because we are working with small subsets of the data. 

The command above first uses the bcftools mpileup command. You can check out what the options are by typing `bcftools` and `bcftools mpileup`. You will see that the `mpileup` command generates genotype likelihoods from multiple `.bam` files. This is then piped to the `bcftools call` command, which does SNP and indel genotype calling.  The `-Ou` and `-Oz` command asks `bcftools` generate uncompressed BCF output and a compressed vcf output respectively.

Because we are working with such a small (subsetted) dataset, your bam files have only very few reads that mapped to your chromosome and they are at very low depth. At this point, let's work with files I made earlier using the complete dataset and pretend that you made them.  GREAT JOB!!!  

Please make a symbolic link to thes files like this:

```
ln -s /2/scratch/Bio722_BJE/chr1_only.recode.vcf.gz .
ln -s /2/scratch/Bio722_BJE/chr2_only.recode.vcf.gz .
ln -s /2/scratch/Bio722_BJE/chrX_only.recode.vcf.gz .
```

And now check out this file like this:

```
zmore chr1_only.recode.vcf.gz

```

Here you need to use `zmore` instead of `more` because the file are compressed (gzipped).

You can also check the header and the first few lines like this:
```
zcat chr1_only.recode.vcf.gz | grep '#CHR' -A 3 
```
This will give you the main comment line (which always begins with #CHR in vcf files) plus 3 lines of data that follow (-A 3).

# Vcf format

A common way that genotype information is conveyed is the `variant call format` – vcf, which is is described [here](https://en.wikipedia.org/wiki/Variant_Call_Format) and [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf). This is a widely used genotyping format that can be uses as input for many softwares (e.g. the R package [PopGenome](https://cran.r-project.org/web/packages/PopGenome/PopGenome.pdf)) or exported to other formats (e.g. using [PLINK](http://zzz.bwh.harvard.edu/plink/)). Another new format introduced by the [`Genome Analysis Toolkit`](https://software.broadinstitute.org/gatk/) of the Broad Institute is called the genomic variant call format – gvcf; this is described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format).

OK, now let's calculate some summary statistics from these files.  Please click [here](https://github.com/evansbenj/BIO722.md/blob/main/7_summary_statistics.md)



