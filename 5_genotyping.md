(Or you can go back to automating the readmap alignment for multiple samples [here](https://github.com/evansbenj/BIO722.md/blob/main/4_automating_readmapping_for_multiple_samples.md).


# Genotyping with bcftools

A common way that genotype information is conveyed is the `variant call format` – vcf, which is is described [here](https://en.wikipedia.org/wiki/Variant_Call_Format). Another new format introduced by the [`Genome Analysis Toolkit`](https://software.broadinstitute.org/gatk/) of the Broad Institute is called the genomic variant call format – gvcf; this is described [here](http://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf).

Now lets make a vcf file from our bam files. We have everything we need to make a genotype (vcf) file with all of the samples, including:
* sorted bam files for each sample
* an index for each bam file (.bai file)
* a genome (actually a chromosome in our case for one chromosome - chrZZZ.fa) and also supporting files (chrZZZ.amb  chrZZZ.fa.ann  chrZZZ.fa.bwt  chrZZZ.fa.fai  chrZZZ.fa.pac  chrZZZ.fa.sa).

Please start a screen and type this command:

```
bcftools mpileup -Ou -f ../my_monkey_chromosome/chr2.fa PF515_chr2_sorted.bam PM561_chr2_sorted.bam PM565_chr2_sorted.bam PM566_chr2_sorted.bam PM567_chr2_sorted.bam PM582_chr2_sorted.bam PM584_chr2_sorted.bam PM592_chr2_sorted.bam PM602_chr2_sorted.bam | bcftools call -mv -Oz -o allsamples_chr9_merged_sorted.bam.vcf.gz

```
This will take a while so please exit the screen by typing `ctrl-a` and `ctrl-d`

The control above first uses the bcftools mpileup command. You can check out what the options are by typing `bcftoos` and `bcftools mpileup`. You will see that the `mpileup` command generates genotype likelihoods from multiple `.bam` files. The `bcftools call` does SNP and indel genotype calling.  The `-Ou` and `-Oz` command asks `bcftools` generate uncompressed BCF output and a compressed vcf output respectively.

This should finish quickly because we are working with small subsets of the data. 

In the meantime, we can check out a file I made earlier.  Please make a symbolic link to this file like this:

```
ln -s /1/scratch/ben/allsamples_chr9_merged_sorted.bam.vcf.gz
```

And now check out this file like this:

```
zmore allsamples_chr9_merged_sorted.bam.vcf.gz

```

Here you need to use `zmore` instead of `more` because the file is compressed (gzipped).
