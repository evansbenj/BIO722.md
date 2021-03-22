(Or you can go back to automating the readmap alignment for multiple samples [here](https://github.com/evansbenj/BIO722.md/blob/main/6_genotyping.md).


# Summary statistics

Now that we have our genotypes, there are lots of things we can do!  We don't have time to much, but as a simple example, let's calculate pairwise nucleotide diversity for our samples using `vcftools`.

This command outputs pairwise nucleotide diversity by site
```
vcftools --gzvcf chr1_only.recode.vcf.gz --chr chr1 --window-pi 10000 --out chr1_analysis
vcftools --gzvcf chr2_only.recode.vcf.gz --chr chr2 --window-pi 10000 --out chr2_analysis
vcftools --gzvcf chrX_only.recode.vcf.gz --chr chrX --window-pi 10000 --out chrX_analysis
```

We can get an average of the diversity across all variable sites like this:
```
cat chr1_analysis.windowed.pi | awk '{sum+=$5} END { print "Average = ",sum/NR}'
cat chr2_analysis.windowed.pi | awk '{sum+=$5} END { print "Average = ",sum/NR}'
cat chrX_analysis.windowed.pi | awk '{sum+=$5} END { print "Average = ",sum/NR}'
```

# Uhhh, why are we doing this again ?!?!

We can calculate the ratio of nucleotide diversity on the X versus the average on the autosomes (with lots of caveats) like this:
```
perl -E "say 3.58237/((4.67155+4.25618)/2)"
```

What does this suggest about the variance in reproductive success in male versus female macaque monkeys?
