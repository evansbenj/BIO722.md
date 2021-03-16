(Or you can go back to the page about preparing your reference genome [here](https://github.com/evansbenj/BIO722.md/blob/main/2_preparing_a_reference_genome.md).

# Read mapping
Please make a new directory within your `my_monkey_chromosome` directory called `my_bam_files` and enter that directory like this:
```
cd ..; mkdir my_bam_files; cd my_bam_files
```

The `mem` algorithm of bwa is a great way to map your reads.  You can do this for the reads from one individual (PF515) as follows:

`bwa mem -M -t 16 -r "@RG\tID:FLOWCELL1.LANE6\tSM:PF515\tPL:illumina" ../my_monkey_chromosome/chrZZZ.fa ../demultiplexed_subsetted_fq/PF515.fastq.gz | samtools view -bSh - > PF515_chrZZZ.bam`

Where `_chrZZZ` should be replaced with whatever chr you are working with. If this executed without error, you should have a new bam file in your directory.

`.bam` files are a compressed binary version of a SAM file, whose format is explained [here](https://samtools.github.io/hts-specs/SAMv1.pdf).  The bit after the `-r` flag is called a read group.  This is required by GATK and is used for base recalibration, which we will discuss if time permits.

You can check out the contents of a bam file like this:
```
samtools view -h PF515_chrZZZ.bam | more
```

We still need to sort and index the bam files. We can do this as follows:
```
samtools sort PF515_chrZZZ.bam -o PF515_chrZZZ_sorted.bam
```
This generates a file called `PF515_chrZZZ_sorted.bam`.  You can index this file like this:
```
samtools index PF515_chrZZZ_sorted.bam
```

## Assessing coverage

Samtools can provide information on the number of reads for each position of the reference sequence for which there are data.  You can see this information by typing this:

`samtools depth XXX_sorted.bam`

Where `XXX` is the sample ID number.  If you want to know the average depth across all sites, you could type this:

`samtools depth XXX_chrZZZ_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'`

Here, as previously, the vertical bar `|` is a "pipe" that sends the information from the command before it to the command after it.  So the data you generated from `samtools` will be parsed with the unix `awk` command.  This will add the values of the third column `$3` to a variable called `sum` and then at the end (`END`) print out the word `Average` followed by the quotient `sum/NR` where `NR` is a built in variable that keeps track of the number of records.  A good description of `awk` is [here](http://www.folkstalk.com/2011/12/good-examples-of-awk-command-in-unix.html).

## Practice Problem 3 (for home)

Using the idxstats options of [samtools](http://www.htslib.org/doc/samtools-0.1.19.html) please check how many reads mapped to your chromosome for your mapped data and how many failed to map.  

According to the [help page](http://www.htslib.org/doc/samtools-idxstats.html): "The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments. It is written to stdout." 

Do you know why so many failed to map?

## OK, if this all went smoothly we are now ready to automate the alignments with a bash script.  Please click [here](https://github.com/evansbenj/BIO722.md/blob/main/4_automating_readmapping_for_multiple_samples.md) to go to the next page.
