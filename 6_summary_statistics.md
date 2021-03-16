(Or you can go back to automating the readmap alignment for multiple samples [here](https://github.com/evansbenj/BIO722.md/blob/main/5_genotyping.md).


# Summary statistics

Now that we have our genotypes, there are lots of things we can do!  We don't have time to much, but as a simple example, let's calculate pairwise nucleotide diversity for our samples using `vcftools`.

This command outputs pairwise nucleotide diversity by site
```
vcftools --gzvcf Example_chr19_macaque.vcf.gz --chr chr19 --site-pi --out chr19_analysis
```

We can get an average of the diversity across all variable sites like this:
```
cat chr19_analysis.sites.pi | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```
