# BJE Lecture 1

(Or you can go back to the Readme page [here](https://github.com/evansbenj/BIO722.md/blob/main/0_README.md)).

## A quick note about genetic samples

Your statistical power, precision, and accuracy will depend on the quantity and quality of your data.  These factors, of course, depend on the starting material you use for sequencing.  Excitingly, a bunch of recent studies suggest that RRGS can be used on sub-optimal samples, for example from museum specimens, fecal extractions, and the like.  Nonetheless, it is in your interest to use as high quality DNA as possible.  
- For animal tissue, I recommend using ethanol or RNAlater to preserve your tissues and I suggest chilling your samples (at -20 or -80 degrees) as soon as possible after collection.
- For DNA extraction, I've had success using Qiagen DNEasy extraction kits.  These kits can be used in any lab that has a heat block and a centrifuge.  When a small amount of starting material is being used I recommend reducing the volume of elution buffer you use in order to concentrate the DNA.  You can always dilute it later, and it is harder and less efficient to make a sample more concentrated.
- If you want to outsource the library preparation (I always have), I recommend using an agarose gel to normalize the concentration of your samples and ensure that the quality is as good as possible.  Different methods have different requirements in terms of the volume and concentration of gDNA, but most expect these parameters to be uniform across the samples you submit.
- If you have many samples, I suggest extracting and comparing more than you plan to run so you can choose the best ones. 

## How much does this cost and where can I do it?

I've done RADseq multiple times at [Floragenex](http://www.floragenex.com/), which is located in Oregon, USA.  I've also done Genotype by Sequencing at [Cornell University] in New York, USA and most recently at [Genome Quebec](https://www.genomequebec.com/). Plenty of places do GBS now. RADseq and GBS are both RRGS methods that differ in the way the libraries are constructed, but both generate similar data. I had a great experience with Genome Quebec and plan to use that facility for GBS in the future. In 2017, the prices for a 95 RADseq sample run, including library preparation but no bioinformatics was  ~US$3300. In 2020, the cost was less (~CAD$3000) for 150 bp single end reads, but for runs on a newer machine that yeilded with much higher coverage. I anticipate that the cost of these services will continue to go down even more over the next few years.

## A quick note about Markdown and Github

This website is written in a markup language called [Markdown](https://en.wikipedia.org/wiki/Markdown) and hosted by [Github](www.github.com).  I've found both of these tools to be easy to learn and very useful.

# Major goals of this exercise
* We will prepare a reference genome using [`bwa`](http://bio-bwa.sourceforge.net/bwa.shtml) using the `index` flag
* For multiple individuals, we will map a subsetted de-multiplexed dataset to individual chromosomes in this reference genome with `bwa` using the `mem`
* We will combine these alignments (`.bam` files) and call genotypes simultaneously for all samples using [`samtools`](http://www.htslib.org/) and [`bcftools`](http://www.htslib.org/)
* We will calculate some diversity statistics using [`vcftools`](https://vcftools.github.io/man_latest.html).


## Example data
The data we will be working with are single end 100 bp reads from one Illumina lane that was collected in 2014. The data are from 9 individuals that were barcoded and multiplexed on this lane (see below for more explanation). The path to the complete dataset is:
```
/home/ben/2021_BIO722/complete_data/forward.fastq
```

Please use the `ls -lh` command to find out how large the complete dataset is:
```
ls -lh /home/ben/2021_BIO722/complete_data/forward.fastq.gz
```

Here the `l` flag asks the list (`ls`) command to provide the *l*ong format that has the file size. The `h` flag asks it to list the file size in a *h*uman readable way (using acronyms for bits). As you (hopefully) can see, the compressed file is fairly large (~11 Gb). The uncompressed file is ~3 times larger. Because the tasks we will perform take a while with this much data, I made a smaller dataset (31Mb) to work with here:

`ls -lh /home/ben/2021_BIO722/complete_data/forward.fastq.gz`

In case you are interested, I made this using the unix `zcat` and `awk` commands as follows:

`zcat forward.fastq | awk 'NR >= 0  && NR <= 500000 { print }' > forward_subset.fastq`

Here the `zcat` command sends each line of 'forward.fastq.gz` to the `awk command. Then the `awk` command searches the number of records `NR` (i.e. the line numbers) from 0-500,000 and prints them to a file called `forward_subset.fastq`.  In this command it is important that the number of lines that you select be divisible by 4, otherwise you will end up with an incomplete fastq entry at the end of the file, which will cause issues.

**FYI, as with most things, I did not figure this out myself, I found it on the internet somewhere.**

## De-Multiplexing
Most RRGS methods rely on the "short" (101-251 bp) read sequencing platform.  These machines generate data using something called a "flowcell" that is divided up into eight "lanes". Small scale projects typically would run multiple samples (from different species or different individuals within a species) on one lane. Because the sequence methodology requires the ligation (attachment) of a linker (a bit of DNA) to each side of bits of DNA that will be sequenced, it is straightforward to combine multiple samples (multiplex) from different individuals in a single lane. This is done by adding a unique identifier sequence (a barcode) to the linker that is used on each sample.  

Note that this barcode is different from "DNA barcoding", the latter of which generally refers to the use of a small variable genomic region (such as the COI gene for animals) for species and population identification.

A first step in our analysis pipeline is to organize data from each of our samples that were run together on an Illumina lane (demultiplexing our data) and also to filter our data and trim off bits that have lots of errors or that have sequences from the laboratory procedures that were used to generate the data (Trimming/Quality control).  

When samples are run on an Illumina machine, DNA is broken up into many small fragments and a small bit of DNA called an adaptor is then added on each of the fragments. This adaptor allows the sequencing process to occur, essentially by making possible high-throughput primer extension with fluorescent terminator nucleotides (ask Ben about this if you are unfamiliar). To make possible the multiplexing of samples on one Illumina lane, each sample is linked to a unique adaptor that contains a "barcode" sequence that allows us to sort out which samples each sequence came from.  For our dataset, we have nine individuals from one species (the Tonkean macaque). Each of the samples received the following barcodes:

```
CCTCTTATCA
TATCGTTAGT
TAGTGCGGTC
GGCCGGTAAC
AGGAACCTCG
TTATCCGTAG
CGCTATACGG
CACGCAACGA
ATCCGTCTAC
```

These barcode correspond with the following sample names:
```
CCTCTTATCA	PF515
TATCGTTAGT	PM561
TAGTGCGGTC	PM565
GGCCGGTAAC	PM566
AGGAACCTCG	PM567
TTATCCGTAG	PM582
CGCTATACGG	PM584
CACGCAACGA	PM592
ATCCGTCTAC	PM602
```
To save time, and also because these days many sequencing centres will demultiplex your data for you, have already used this information to de-multiplex our data. I put these data here:
```
/2/scratch/Bio722_BJE/demultiplexed_subsetted_fq/
```

I used software called `Stacks` to de-multiplex the data.  This is actually a suite of programs and I used the application called `process_radtags` within `Stacks`.  `Stacks` has a very nice online manual [here](http://catchenlab.life.illinois.edu/stacks/manual). Other software that does demultiplexing of RADseq data includes [radpools](https://github.com/johnomics/RADtools) and [cutadapt](https://cutadapt.readthedocs.io/en/stable/).

For those of you that are curious, I have copied below [More Details On Demultiplexing](#more-details-on-demultiplexing).  


## Trimming and Quality Control of NextGen Data

We have already discussed [fastq](http://en.wikipedia.org/wiki/FASTQ_format) format, quality assessment with [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and trimming with [TROMMOMATIC](http://www.usadellab.org/cms/?page=trimmomatic). So I won't go over this again.


## Set up a directory on scratch and make symbolic links

Please login to info, connect to info115 (rsh info115) and navigate to the scratch directory as follows:

`cd /2/scratch/ZZZZ`, where `ZZZ` is your username (e.g. `gradstd13`).

Make a directory to work in, and enter that directory:
```
mkdir monkey_directory; cd monkey_directory
```

Next, please make a symbolic link to a subsetted dataset (`ln -s /2/scratch/Bio722_BJE/demultiplexed_subsetted_fq/`).  This is a useful way to manage information and avoid duplicating files and using up Brian's disk space. You can check out the data like this:
```
cd demultiplexed_subsetted_fq
```
and
```
zcat PM582.fastq.gz | more
```

OK, so now we should all have some data ready for us to work with.  Later, for funsies, you can try the [Practice Problems](#practice-problems) below.  For now, let's move on to [prepare a reference genome](https://github.com/evansbenj/BIO722.md/blob/main/2_preparing_a_reference_genome.md).


## More Details On Demultiplexing

A first step in analysis of Illumina data is to identify adaptor and barcode sequences in our data, sort sequences by the barcode, and remove adaptor and barcode sequences from the data.  We can also get rid of sequences that have ambiguous barcodes due to sequencing errors. We additionally can get rid of sequences with low quality scores and trim them all so they have the same length (this last step would not normally be done for RNAseq data but it is a reasonable thing to do for RADseq data).

Illumina generates sequences that have errors in base calls.  Errors typically become more common towards the end of the sequence read, and sometimes (but not always) an "N" is inserted in positions where the base pair is difficult to call.  But sometimes it makes an incorrect call as well. 

We will use software package called `Stacks` to de-multiplex and trim our data.  This is actually a suite of programs and we will be using the application called `process_radtags` within `Stacks`.  `Stacks` has a very nice online manual [here](http://catchenlab.life.illinois.edu/stacks/manual). FYI, other software that does trimming of RADseq data is available [here](https://github.com/johnomics/RADtools/blob/master/RADpools).

Brian has installed most of the software we need in a directory called `/usr/local/bin`. Before we de-multiplex our data subset, we need to make a directory for the de-multiplexed data to be stored in. From your current directory  (`/2/scratch/your_usrname`), please type this:

`mkdir samples`

The command to execute this program on our data is:

`/usr/local/bin/process_radtags -f <inputfile> -b <barcode_file> -o ./samples/ -e sbfI -t 75 -r -c -q --adapter_1 GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --adapter_mm 2 --filter_illumina`

As detailed in the online manual, the first part (`/usr/local/bin/process_radtags`) directs the computer to run the program process_radtags, which is in the director called `/usr/local/bin/`.  The `-f` flag specifies where the data are and <inputfile> provides the path and filename of the data. The `-b` flag specifies where the barcode file is that we made eariler and <barcode_file> provides the path and name for this file.  The `-o` flag tells `process_radtags` to put the demultiplexed data in a folder called `./samples`, which we should make in advance using the unix `mkdir` command.  The `-e` flag tells `process_radtags` that the restriction enzyme called sbfI was used to generate the data. The `-t` flag tells `process_radtags` to trim all sequences to be 75 base pairs long. The `-r`, `-c`, and `-q` flags directs `process_radtags` to respectively
- rescue barcodes and RADtags when possible allowing upto a default value of 2 mismatches (this number can be changed too if you want)
- clean the data and remove reads with any uncalled bases, and 
- discard reads with low quality scores.  

The other flags tell `process_radtags` to remove adapter sequences that are specified and to remove bad reads recognized by the Illumina sequencing software.  Other details, such as the type of quality scores are set at default values. All of this information is available, of course, in the manual that comes with the program. 

Here is an example of the commandline I used to de-multiplex the subset of the data:

`/usr/local/bin/process_radtags -f forward_subset.fastq -b monkey.barcodes -o ./samples/ -e sbfI -t 75 -r -c -q --adapter_1 GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --adapter_mm 2 --filter_illumina`

If you enter the `samples` directory, you should see your de-multiplexed files, each named by the barcode. Let's check the log file:

`more samples/process_radtags.log`

You should see statistics on the number of reads retained or rejected, why they were rejected, and how many reads per individual were retained. Rejection of reads happens if the software cannot find a SbfI restriction enzyme site (ambiguous RADtag), or cannot place the barcode (ambiguous barcode), if it is low quality (lots of Ns), or if it contains Illumina adaptor sequences.

Please rename each of these nine files to have the sample name instead of the barcode using the `mv` command.  For example, within the `samples` directory:

```
mv sample_CCTCTTATCA.fq PF515.fq
mv sample_TATCGTTAGT.fq PM561.fq
mv sample_TAGTGCGGTC.fq PM565.fq
mv sample_GGCCGGTAAC.fq PM566.fq
mv sample_AGGAACCTCG.fq PM567.fq
mv sample_TTATCCGTAG.fq PM582.fq
mv sample_CGCTATACGG.fq PM584.fq
mv sample_CACGCAACGA.fq PM592.fq
mv sample_ATCCGTCTAC.fq PM602.fq
```

# Practice Problems
## Practice problem 1: How many reads do we have for each individual?

As an exercise, please use the [`grep`](https://man7.org/linux/man-pages/man1/grep.1.html) command to count how many reads we have for each individual.  A hint is that using `grep`, you can count the number of times an identifier character for each sequence appears in each file for each individual.  Another hint is that you can get the manual for any `Unix` command by typing `man command`.  Which individual has the most reads?  Which has the least reads?  Can you think of a reason that some samples have lots of reads while others have less? 

## Practice problem 2: 
How would your `grep` command differ for a `fasta` file compared to a `fastq` file?


## OK, now we are ready to prepare a reference genome.  Please click [here](https://github.com/evansbenj/BIO722.md/blob/main/2_preparing_a_reference_genome.md).



