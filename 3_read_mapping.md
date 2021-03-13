(Or you can go back to the page about preparing your reference genome [here](https://github.com/evansbenj/BIO722.md/blob/main/2_preparing_a_reference_genome.md).

# Read mapping
Please make a new directory within your `my_monkey_chromosome` directory called `my_bam_files` and enter that directory like this:
```
cd ..; mkdir my_bam_files; cd my_bam_files
```

The `mem` algorithm of bwa is a great way to map your reads.  You can do this as follows:

`bwa mem -M -t 16 -r "@RG\tID:FLOWCELL1.LANE6\tSM:PF515\tPL:illumina" ../my_monkey_chromosome/chr2.fa demultiplexed_subsetted_fq/PF515.fastq.gz | samtools view -bSh - > PF515.bam`

If this executed without error, you should have a new bam file in your directory.

`.bam` files are a compressed binary version of a SAM file, whose format is explained [here](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm#:~:text=A%20BAM%20file%20(*.,sequences%20up%20to%20128%20Mb.&text=Alignments%E2%80%94Contains%20read%20name%2C%20read,alignment%20information%2C%20and%20custom%20tags.).

You can check out the contents of a bam file like this:
```
samtools view PF515.bam | more
```


## Practice Problem 4: Assessing coverage

Samtools can provide information on the number of reads for each position of the reference sequence for which there are data.  You can see this information by typing this:

`samtools depth XXX_sorted.bam`

Where `XXX` is the sample ID number.  If you want to know the average depth across all sites, you could type this:

`samtools depth XXX_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'`

Here, as previously, the vertical bar `|` is a "pipe" that sends the information from the command before it to the command after it.  So the data you generated from `samtools` will be parsed with the unix `awk` command.  This will add the values of the third column `$3` to a variable called `sum` and then at the end (`END`) print out the word `Average` followed by the quotient `sum/NR` where `NR` is a built in variable that keeps track of the number of records.  A good description of `awk` is [here](http://www.folkstalk.com/2011/12/good-examples-of-awk-command-in-unix.html).

## Practice Problem 5 (for home): De-multiplexing the complete dataset and mapping the data to your reference chromosome for one individual

Using the same pipeline we have just gone through, please do the following:
* demultiplex the complete dataset
* rename the resulting fastq files to match the sample names instead of the barcode sequences
* map the complete data from one individual (e.g. PF515) to your reference chromosome.

Now, using the idxstats options of [samtools](http://www.htslib.org/doc/samtools-0.1.19.html) please check how many reads mapped to your chromosome for your mapped data and how many failed to map.  Do you know why so many failed to map?

## OK, if this all went smoothly we are now ready to automate the alignments with a bash script.  Please click [here](https://github.com/evansbenj/BIO720/blob/master/3_Lecture_3_Automating_alignment_with_bash.md) to go to the next page.
