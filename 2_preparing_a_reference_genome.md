# Preparing a reference genome
(or you can go back to the evaluating and demultipexing data page [here](https://github.com/evansbenj/BIO722.md/blob/main/1_Evaluating_your_data.md).

Depending on your organism of study, there may or may not be a relatively closely related genome sequence to work with.  Depending on your research question, this may or may not be useful.  In our example study on Tonkean macaques, we are  interested in quantifying molecular polymorphism on the X chromosome and comparing it to polymorphism on the autosomes.  For this reason, the genomic location of the data is crucial and we can benefit from the complete genome sequence of a closely related species of macaque monkey, the rhesus macaque (*Macaca mulatta*).  We will use software called [bwa](http://sourceforge.net/projects/bio-bwa/files) and [samtools/bcftools](http://samtools.sourceforge.net/), to map our data to individual chromosomes of the rhesus macaque.  Normally one would map reads to an entire genome because the data were generated from a complete genome, but in our case we are doing only an example analysis and we will each work on an individual chromosome. Ben will assign each of you a chromosome to work on.

## A note about "completely" sequenced genomes

FYI, essentially all completely sequenced genomes are not in fact completely sequenced.  
- Regions such as centromeric and telomeric regions and some portions of sex-specific sex chromosomes contain many repetitive elements that pose challenges to sequencing and assembly.  
- Sometimes the individual sequenced is female, so no Y chromosome is available.  Sometimes (usually?) when a genome is said to be "complete" it actually is a bunch of "contigs", or contiguous sequence, that may or may not be assembled into "scaffolds" that include contigs plus Ns to represent intervenining regions that are not yet sequenced.  
- And even then, we can expect sequence and assembly errors in our reference genome that make it different from the real genome sequence.  
- On top of that, there is population level variation to contend with, including SNPs and insertion deletion events.  This makes our samples different from any reference genome as well.

As an example, let's look at some information on the "completely" sequenced genomes of [some frogs](http://www.xenbase.org/other/static-xenbase/ftpDatafiles.jsp).  Of interest is the N50 statistic of a genome assembly, which is defined [here](https://en.wikipedia.org/wiki/N50_statistic).

## Preparing your reference genome

Reference genomes for many sequences are available at multiple publicly available databases.  We can download the complete genome sequence for the rhesus macaque from the [USC genome browser](http://hgdownload.cse.ucsc.edu/downloads.html#rhesus).  I did this earlier because it takes a while.  The whole genome comes as a fasta-formatted file with individual entries corresponding to individual chromosomes. 
```
/2/scratch/Bio722_BJE/rhesus_genome/macaque_masked_chromosomes_ym.fasta
```
You can check out this genome like this:
```
head /2/scratch/Bio722_BJE/rhesus_genome/macaque_masked_chromosomes_ym.fasta
```
I split it up into individual fasta files corresponding with each of the chromosomes. This is easy to do with samtools if you know the names in the header for each chromosome, e.g., for chr2:
```
samtools faidx /2/scratch/Bio722_BJE/rhesus_genome/macaque_masked_chromosomes_ym.fasta chr2 > chr2.fast
```

I put the individual chromosomes in this directory:

```
/2/scratch/Bio722_BJE/rhesus_chromosomes
```

Now check out what is in this directory by typing this:

`ls /2/scratch/Bio722_BJE/rhesus_chromosomes`

Ben will assign you a chromosome to work with.  From the `/2/scratch/YOUR_USERNAME/monkey_directory/` directory, please make a symbolic link to this chromosome reference sequence in a new directory that you make like this:

`mkdir my_monkey_chromosome; cd my_monkey_chromosome`

`ln -s /2/scratch/Bio722_BJE/rhesus_chromosomes/chr`ZZZ`.fa .` 

Here and henceforth, you will need to change the `chrZZZ.fa` part to match whatever chromosome Ben assigned to you.  For example, if you are working on chromosome 9, you should type this:

`ln -s /1/scratch/ben/rhesus_chromosomes/chr9.fa .` 

Before we map our data to this reference genome, we need to generate some files that will be used in the mapping process.  This can be done in three steps:

1. Make an index file.   

    The `bwa` command tells the computer to execute the bwa program.  The `index` command tells `bwa` to generate index files from the rhesus genome file that is indicated by the `my_monkey_chromosome/chr`ZZZ`.fa`. The `-a bwtsw` flag specifies the indexing algorithm for `bwa` to use.  
  
  This step will take a few minutes. If you are feeling adventurous (and I hope you are!), you can do this in a screen like this:
  
  `screen -S make_an_index_file`
  
  then type this:
  
  `bwa index chr`ZZZ`.fa`
  
  Then exit the screen by typing `ctrl-a` then `ctrl-d`
  
  You can list the screens you have like this:
  
  `screen -ls`

  and you can return to the screen you started like this:
  
  `screen -r make_an_index_file`
  
  when it is done, you can exit and then kill the screen like this:
  
  `ctrl-a` then `ctrl-d` and then
  
  screen -X -S make_an_index_file kill

Or you can exist a screen while attached by typing `ctrl-a` and then typing `:quit`.



2. Some software need an `fai` file that can be generated using `samtools`.  Please type this:

  `samtools faidx chr`ZZZ`.fa`

  Here, the `samtools` command tells the computer to execute the `samtools` program.  The `faidx` option tells samtools to generate a file called `chr`ZZZ`.fai` in which each line has information for one the contigs within the reference genome, including the contig name, size, location and other information.  Our reference genome has a contig for each chromosome.

3.  If you are going to use GATK for genotyping (which we are not), you will also need a `.dict` file, which can be generated with a program called [`picard`](http://broadinstitute.github.io/picard/).  To do this, we need to use the latest version of java.  Please type this command to make a symbolic link to this version:

`ln -s /home/ben/jre1.8.0_111/bin/java`

OK, now run picard like this:

  `java -jar /usr/local/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=chr`ZZZ`.fa OUTPUT=chr`ZZZ`.dict`

  As before, you will need to change the `chr`ZZZ in this command to match the chromosome you are working with.  This should generate a file called `chr`ZZZ`.dict`
  
  
  ## OK, now we are ready to map reads to our reference genome (chromosome).  Please click [here](https://github.com/evansbenj/BIO722.md/blob/main/3_read_mapping.md)
