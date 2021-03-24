Or go back to summary statistics [here](https://github.com/evansbenj/BIO722.md/blob/main/7_summary_statistics.md).

# Genotype filtering

It is often the case, despite our efforts to generate high quality genotype calls, that some genotypes just don't make sense.  For example, we might observe a heterozygous genotype on the male specific portion of the Y chromosome or we might see some genotypes from a female on the Y chromosome. Genotype errors arise from low or patchy coverage in some individuals, duplicated regions in the sample or reference, sequencing errors, mapping errors, and more. 

Generally before analysis, one filters vcf files to get rid of low quality genotypes.  There are many programs that can do this, including [vcftools](http://vcftools.sourceforge.net/man_latest.html), [bcftools](http://samtools.github.io/bcftools/bcftools.html), and [GATK](https://gatk.broadinstitute.org/hc/en-us).  Probably the most sophisticated approach is to use GATK's variant quality score recalibration (`VQSR`) which uses a model-based approach to assess quality scores and filter genotypes based on multi-dimensional criteria. This generally requires a large (>30 individuals) sample side and some fairly strong assumptions (e.g. a set of known variants). Instead, today, we will implement a hard-cutoff filter, also using GATK.

In order to work with vcf files, `GATK` needs us to index them first.  Please do this with `tabix`:
```
tabix -p vcf chr1_only.recode.vcf.gz
tabix -p vcf chr2_only.recode.vcf.gz
tabix -p vcf chrX_only.recode.vcf.gz
```
This creates a `.tbi` file for each of our vcf files.

Here's an example of a script that can do these steps.  Please copy this and make an executeable file and execute it on the subset data we have been working with. You could name it `Step_4_flag_and_filter.pl`.


``` perl
#!/usr/bin/perl
                                                                                                                         
use warnings;
use strict;

# This script will use GATK's SelectVariants to make a vcf file with only indels in it.                                  
# Then it will use VariantFiltration to mark indels and other potentially low quality sites near indels                  
# Then it will use SelectVariants to output an filtered vcf file from which filtered positions have been removed.        

# It will take two variables as input: (1) the prefix of the reference genome and                                        
# (2) the prefix of a gzipped vcf file                                                                                   
# execute like this:                                                                                                     
# Step_4_flag_and_filter.pl chr1 chr1_only.recode                                                                        

my $path_to_reference_genome="/2/scratch/USERNAME/monkey_directory/my_monkey_chromosome/";
my $reference_genome=$ARGV[0].".fa";
my $status;
my $vcffile=$ARGV[1].".vcf.gz";

# make a file with only indels using SelectVariants                                                                      
my $commandline = "java -Xmx2G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R ".$path_to_reference_genome.$reference_genome;
$commandline = $commandline." --variant ".$vcffile." -selectType INDEL -o ".$vcffile."_indels_only.vcf";
print $commandline,"\n";
$status = system($commandline); 


# flag the vcf file using the indel file and also flag some other low quality variants                                   
$commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R ".$path_to_reference_genome.$reference_genome;
$commandline = $commandline." -o ".$vcffile."_flagged.vcf --variant ".$vcffile;
$commandline = $commandline." --mask ".$vcffile."_indels_only.vcf --maskName INDEL --maskExtension 10";
$commandline = $commandline." --filterExpression \"QD < 2.0\" --filterName \"QD2\"";
$commandline = $commandline." --filterExpression \"QUAL < 30.0\" --filterName \"QUAL30\"";
$commandline = $commandline." --filterExpression \"SOR > 3.0\" --filterName \"SOR3\"";
$commandline = $commandline." --filterExpression \"FS > 60.0\" --filterName \"FS60\"";
$commandline = $commandline." --filterExpression \"MQ < 40.0\" --filterName \"MQ40\"";
$commandline = $commandline." --filterExpression \"MQRankSum < -12.5\" --filterName \"MQRankSum-12.5\"";
$commandline = $commandline." --filterExpression \"ReadPosRankSum < -8.0\" --filterName \"ReadPosRankSum-8\"";
print $commandline,"\n";
$status = system($commandline);                                                                                        

# output a new filtered genotype file using SelectVariants                                                                                                                                              
$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R ".$path_to_reference_genome.$reference_genome;
$commandline = $commandline." --variant ".$vcffile."_flagged.vcf -o ".$vcffile."_filtered.vcf.gz -select \'vc.isNotFiltered()\'";
print $commandline,"\n";
$status = system($commandline); 

```

Note that this part of the `VariantFiltration` commandline:

` --mask indels_only.vcf --maskName INDEL --maskExtension 10`

will flag indels and also sites +/- 10 bp from indels.

Other examples one could add to the VariantFilteration command line for Ben to discuss include:

`--filterExpression "DP < 5" --filterName "LowCoverage" `

which flags sites with an average depth of coverage less than 5,

`--filterExpression "MQ0 >= 4 && ((MQ0 / DP) > 0.1)" --filterName "HARD_TO_VALIDATE" `

which flags sites with at least 4 reads that map well to another part of the genome and where these reads comprise more than 10% of the reads at that position.


