# Goals

The pipeline we ran through last time took you from raw next gen data to a summary statistic of whole genome genotypes.  But we skipped some important steps which I'd like to run though today.  This includes:
* Deduping reads in bam files and
* Realigning indels

I'll then discuss BQSR using GATK but we won't actually do this.  After that, together we will
* Filter genotypes

Ok let's get started. I've prepared an example `.bam` file that I'd like us to work with.
