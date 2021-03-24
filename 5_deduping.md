
Or you can go back to where we left off last time [here](https://github.com/evansbenj/BIO722.md/blob/main/4_automating_readmapping_for_multiple_samples.md)

# Goals

The pipeline we ran through last time took you from raw next gen data to a sorted bam file for each sample. Before we genotype these data, we need to do two more steps:

* Realigning indels
* Deduping reads and


After that, we will:
* genotype our data 
* Filter genotypes
* Calculate summary statistics

# First let's make sure we are all on the same page...

If you didn't manage to make it to the point we left off, please make a folder in your `/2/scratch/USERNAME` directory called `monkey_directory` and enter that directory like this:
```
cd /2/scratch/USERNAME; mkdir monkey_directory; cd monkey_directory
```

Now please make some more directories that will hold our bam files and reference genome (chromosome):
```
mkdir my_bam_files; mkdir my_monkey_chromosome
```

Now please enter your my_bam_files directory (`cd my_bam_files`) and make some symbolic links to these `.bam` and the corresponding index files (`.bai`) files that I made for you from within your `my_bam_files` folder:
```
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PF515_chr1_sorted.ba* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PM561_chr1_sorted.ba* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PM565_chr1_sorted.ba* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PM566_chr1_sorted.ba* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PM567_chr1_sorted.ba* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PM582_chr1_sorted.ba* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PM584_chr1_sorted.ba* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PM592_chr1_sorted.ba* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_bam_files/PM602_chr1_sorted.ba* .
```
Now please go over to your `my_monkey_chromosome` directory (`cd ../my_monkey_chromosome`) and make some links to these reference files and their supporting files with the same prefixes:
```
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_monkey_chromosome/chr1* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_monkey_chromosome/chr2* .
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_monkey_chromosome/chrX* .
```
For those of you that already have a `.bam` file for each sample for your reference chromosome (chrZZZ.fa) please enter your `my_monkey_chromosome` and, while you are in this directory, please make a `.fai` file for your reference chromosome like this:
```
samtools faidx chrZZZ.fa
```
or if you want to work with the ones I made, type this:
```
samtools faidx chr1.fa
samtools faidx chr2.fa
samtools faidx chrX.fa
```
Also, please make a `.dict` file like this:
```
java -jar /usr/local/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=chrZZZ.fa OUTPUT=chrZZZ.dict
```
or if you want to work with the ones I made, type this:
```
java -jar /usr/local/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=chr1.fa OUTPUT=chr1.dict
java -jar /usr/local/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=chr2.fa OUTPUT=chr2.dict
java -jar /usr/local/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=chrX.fa OUTPUT=chrX.dict
```
We will need both of these files for steps involving `GATK`..

OK now go back to the `my_bam_files` directory please (`cd ../my_bam_files`). 


# Realigning indels part 1: Identifying the indels

OK hopefully all of us are all set now. We will use `GATK` to identify indels that may be associated with inappropriate mapping differences among the individuals in our study, and then realign them across all individuals. This is done in two steps.  The first uses the `GATK` function called `RealignerTarget` to identify indels in our data. This function produces a text file that has information about the locations of all indels in any individual relative to the reference genome. For class we are working with a subsetted dataset, so this will be pretty quick (seconds).

Here is a perl script that will execute the `RealignerTarget` function in `GATK` on our data:

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute the GATK command "RealignerTargetCreator" on these files. 
# It takes as input a variable, which is the chromosome name (e.g. chr1)

my $path_to_reference_genome="/2/scratch/USERNAME/monkey_directory/my_monkey_chromosome/";
my $reference_genome=$ARGV[0].".fa";
my @files;
my $status;
   
@files = glob("*".$ARGV[0]."_sorted.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-R ".$path_to_reference_genome.$reference_genome." -o forIndelRealigner.$ARGV[0].intervals";

$status = system($commandline);
```

Please copy and paste this script into a new file, edit the path so that the `USERNAME` is your username, change the permissions to allow it to be executable (`chmod 755 your_file.pl`), and execute it on your samples from within the directory that contains your sorted bam files. I suggest naming your scripts using a sensible system, such as with descriptions and numbers in ascending order.  For example, you could name this script `Step_1_execute_GATK_RealignerTargetCreator.pl`. This script takes as input the name of your reference chromosome (e.g., `chr9`) so you can execute the script like this:
```
./Step_1_execute_GATK_RealignerTargetCreator.pl chrZZZ
```

In this script the `@files = glob("*_sorted.bam");` command uses the `glob` function to look for all files in your directory with that end with `chrZZZ_sorted.bam` and then it adds them to an array called `@files`. The commandline considers all of these `.bam` files together, so there is a loop to add each one to the commandline. The `$status = system($commandline);` line executes the test stored in the `$commandline` variable. The commandline executes the java `.jar` file and allocates a maximum of 1 Gb of memory for the Java virtual machine (`-Xmx1G`).

When it is done, please check out the file it made like this:

`more forIndelRealigner.chrZZZ.intervals`

You should see a list of intervals coordinates following the chromosome of your reference chromosome.

With the indel text file, we can then use a function called `IndelRealigner`, which takes as input this `vcf` file to realign bases when possible an minimize mis-called SNPs.

# Realigning indels part 2: Now realinging the indels across all samples

Now that we know where to move bits of the alignment around for all samples, let's do it!

Here is a perl script that executes the `IndelRealigner` function:

```perl
#!/usr/bin/perl                                                                                                                      
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and                                                           
# make and execute a GATK command IndelRealigner on these files.  
# It takes as input a variable, which is the chromosome name (e.g. chr1)
# and requires that you have run RealignerTargetCreator beforehand

my $path_to_reference_genome="/2/scratch/USERNAME/monkey_directory/my_monkey_chromosome/";
my $reference_genome=$ARGV[0].".fa";
my $status;
my @files;

@files = glob("*".$ARGV[0]."_sorted.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T IndelRealigner ";

foreach(@files){
    $commandline = $commandline." -I ".$_;
}

$commandline = $commandline." -R ".$path_to_reference_genome.$reference_genome." --targetIntervals forIndelRealigner.".$ARGV[0].".intervals --nWayOut _realigned.bam";

$status = system($commandline);

```

As above, please copy and paste this script, make it executable, and execute it. You could name this script `Step_2_execute_GATK_IndelRealigner.pl`. Please don't forget to modify the path of your reference chromosome as appropriate.  After making the script executable (`chmod 755 Step_2_execute_GATK_IndelRealigner.pl`), you can execute it like this:
```
./Step_2_execute_GATK_IndelRealigner.pl chrZZZ
```

If your run completed successfully, you should see new bam files that have the ending `*_sorted_realigned.bam`. These files should be about the same size as the `*sorted.bam` files.  Please check this by typing `ls -l`.

Some points worth noting:
* It may not be obvious, but we have been very careful to not overwrite our input files by having output files with the same name  
* Also very important is that when you are doing chromosome-specific steps, your supporting files should have these chromosomes integrated in their names.  This will prevent you accidentally using an input file from one chromosome for another chromosome.
* We also have named our `.bam` files in a way that provides information about where we are in the pipeline.  When you are happy with where you are, you can delete the upstream files later to save space.
* We have named our scripts in a way that indicates what order they were run in (e.g., Step_1_XXX, Step_2_XXX). We should have done this with the earlier alignment file too (we could go back and rename that one "Step_0_xxx).  This is better than naming your files `newnew_final_for_sure.pl`.
* We are using scripts that accept variables as input  so we can reuse them for multiple chromsomes, and that automatically detect files in a directory to make sure we don't accidentally miss one due to a typo.
* We have used comments to remind us what each script does


# Dedupling a bam file

Many approaches to library construction have a PCR step, which means that some data you generate may have the same start and stop points for both pairs of a read.  This causes an overestimation of confidence and can lead to miscalled genotypes if one allele amplified more than the other.  Typically we therefore remove duplicates ("dedup") before genotyping. 

```
#!/usr/bin/perl 
                                                                                                                 
use warnings;
use strict;

# This script will read in the *_sorted_realigned.bam file names in a directory, and                          
# dedups them using picard 
# It takes as input a variable, which is the chromosome name (e.g. chr1)


my $path_to_reference_genome="/2/scratch/USERNAME/monkey_directory/my_monkey_chromosome/";
my $reference_genome=$ARGV[0].".fa";
my $status;
my @files;

@files = glob("*".$ARGV[0]."_sorted_realigned.bam");

my $commandline = "java -Xmx1G -jar /usr/local/picard-tools/picard.jar MarkDuplicates ";

foreach(@files){
    my $commandline = "java -Xmx1G -jar /usr/local/picard-tools/picard.jar MarkDuplicates ";
    $commandline = $commandline." I= ".$_." O=".substr($_,0,-4)."_dedup.bam" ;
    $commandline = $commandline." M=".substr($_,0,-20)."marked_dup_metrics.txt";

    $status = system($commandline);
}

```

Please copy and paste this into a new file called `Step_3_dedup.pl`.  After modifying the reference genome path to fix the `USERNAME` and then making the script executable (`chmod 755 Step_3_dedup.pl`), you can execute it like this:
```
./Step_3_dedup.pl chrZZ
```

## OK now we are ready to combine our bam files and do some genotyping. Please click [here](https://github.com/evansbenj/BIO722.md/blob/main/6_genotyping.md)
