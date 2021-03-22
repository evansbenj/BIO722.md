# Goals

The pipeline we ran through last time took you from raw next gen data to a summary statistic of whole genome genotypes.  But we skipped some important steps which I'd like to run though today.  This includes:
* Deduping reads in bam files and
* Realigning indels

I'll then discuss BQSR using GATK but we won't actually do this.  After that, together we will
* Filter genotypes

# Dedupling a bam file

Many approaches to library construction have a PCR step, which means that some data you generate may have the same start and stop points for both pairs of a read.  This causes an overestimation of confidence and can lead to miscalled genotypes if one allele amplified more than the other.  Typically we therefore remove duplicates ("dedup") before genotyping.  We didn't do this last time (for the sake of time) but we are gonna do it now.

I've prepared an example `.bam` file that I'd like us to work with.  Please make a symbolic link to this file:
```
ln -s /2/scratch/Bio722_BJE/PM561_chr9_sorted.bam
```


