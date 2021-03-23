
Or you can go back to where we left off last time [here](https://github.com/evansbenj/BIO722.md/blob/main/4_automating_readmapping_for_multiple_samples.md)

# Goals

The pipeline we ran through last time took you from raw next gen data to a sorted bam file for each sample. Before we genotype these data, we need to do two more steps:

* Realigning indels
* Deduping reads and


After that, we will:
* genotype our data 
* Filter genotypes
* Calculate summary statistics

# Realigning indels

If you didn't manage to make it to the point we left off, please make symbolic links to these `.bam` and the corresponding index files (`.bai`) files from within your `my_bam_files` folder:
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
You will also need to go up a directory and into your `my_monkey_chromosome` directory and make some links to these reference files
```
ln -s /2/scratch/Bio722_BJE/monkey_directory/my_monkey_chromosome/chr1* .
```

While you are in this directory, please make a `.fai` file for your reference chromosome like this:
```
samtools faidx chrZZZ.fa
```

Also, please make a `.dict` file like this:
```
java -jar /usr/local/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=chrZZZ.fa OUTPUT=chrZZZ.dict
```


OK now go back to the `my_bam_files` directory please. 

We will use `GATK` to identify indels that may be associated with inappropriate mapping differences among the individuals in our study, and then realign them across all individuals. This is done in two steps.  The first uses the `GATK` function called `RealignerTarget` to identify indels in our data. This function produces a text file that has information about the locations of all indels in any individual relative to the reference genome. Because the second step takes a few dozen minutes to run with the full dataset, for class we will work with the subset datasets that you previously mapped to your chromosome.

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
   
@files = glob("*".$reference_genome."_sorted.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-R ".$path_to_reference_genome.$ARGV[0]." -o forIndelRealigner.intervals";

$status = system($commandline);
```

Please copy and paste this script, change the permissions to allow it to be executable, and execute it on your samples from within the directory that contains your sorted bam files. I suggest naming your scripts using a sensible system, such as with descriptions and numbers in ascending order.  For example, you could name this script `Step_1_execute_GATK_RealignerTargetCreator.pl`. You will need to adjust the `$reference_genome` variable to match your chromosome.

In this script the `@files = glob("*_sorted.bam");` command uses the `glob` function to look for all files in your directory with that end with "_sorted.bam" and then it adds them to an array called `@files`. The `$status = system($commandline);` line executes the test stored in the `$commandline` variable. The commandline executes the java `.jar` file and allocates a maximum of 1 Gb of memory for the Java virtual machine (`-Xmx1G`).

When it is done, please check out the file it made like this:

`more forIndelRealigner.chrZZZ.intervals`

You should see a list of intervals coordinates following the chromosome of your reference chromosome.

With the indel text file, we can then use a function called `IndelRealigner`, which takes as input this `vcf` file to realign bases when possible an minimize mis-called SNPs.

Here is a perl script that executes the `IndelRealigner` function:

```perl
#!/usr/bin/perl                                                                                                                      
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and                                                           
# make and execute a GATK command IndelRealigner on these files.  
# It takes as input a variable, which is the chromosome name (e.g. chr1)

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

As above, please copy and paste this script, make it executable, and execute it. You could name this script `Step_2_execute_GATK_IndelRealigner.pl`. Please don't forget to modify the name of your reference chromosome as appropriate.

If your run completed successfully, you should see new bam files that have the ending `*_sorted_realigned.bam`. These files should be about the same size as the `*sorted.bam` files.  Please check this by typing `ls -l`.


# Dedupling a bam file

Many approaches to library construction have a PCR step, which means that some data you generate may have the same start and stop points for both pairs of a read.  This causes an overestimation of confidence and can lead to miscalled genotypes if one allele amplified more than the other.  Typically we therefore remove duplicates ("dedup") before genotyping. 

```
#!/usr/bin/perl                                                                                                 \
                                                                                                                 
use warnings;
use strict;

# This script will read in the *_sorted_realigned.bam file names in a directory, and                          
# dedup them using picard                                                                                         

my $path_to_reference_genome="/2/scratch/evanslab/monkey_directory/my_monkey_chromosome/";
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



## OK now we are ready to combine our bam files and do some genotyping. Please click [here](https://github.com/evansbenj/BIO722.md/blob/main/6_genotyping.md)
