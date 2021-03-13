# 4 Automating readmapping for multiple samples

(Or you can go back to the page about preparing your reference genome [here](https://github.com/evansbenj/BIO722.md/blob/main/3_read_mapping.md).

Now that you have seen how to align data from one individual to a reference genome, we can automate the alignment of all individuals to the reference genome using a bash script. This is easier than going through all that stuff independently for each individual. We can accomplish this with a `bash` script by defining an `array` that contains the names of all of the individuals in the analysis, and then loop through this array and execute each of the commands for each individual.

Below is an example `bash` script that should run all of our analyses for each individual.  Please use a text editor to make this program.  In the beginning of the script 5 variables are defined that specify, respectively, the path for the bwa and samtools programs, the path to the data, the path to the reference chromosome, and the name of the chromosome you are working on. You will need to modify the variables somewhat to match the chromosome you are working on and the directory. For example you should use the arrow keys to scroll up to the line that says `chromosome="chrXXX.fa"` and change the part that says `chrXXX.fa` to correspond with whatever chromosome you are working on.  For example, if youa re working on chromosome 10, please change this to instead read `chr10.fa`. Also, in the `path_to_chromosome` variable, you will need to change the part that says `YYY` to match your home directory name.

```
#!/bin/bash                                                                                            

path_to_data="/2/scratch/evanslab/samples"
path_to_chromosome="/2/scratch/evanslab/my_monkey_chromosome/"
chromosome="chrZZZ"

individuals="PF515                                                                                     
PM561                                                                                                  
PM565                                                                                                  
PM566                                                                                                  
PM567                                                                                                  
PM582                                                                                                  
PM584                                                                                                  
PM592                                                                                                  
PM602"

for each_individual in $individuals
do
    echo ${each_individual}
    bwa mem -M -t 16 -r "@RG\tID:FLOWCELL1.LANE6\tSM:${each_individual}\tPL:illumina" ../my_monkey_chromosome/${chromosome}.fa demultiplexed_subsetted_fq/${each_individual}.fastq.gz | samtools view -bSh - > ${each_individual}_chrZZZ.bam
    samtools sort ${each_individual}_$chromosome.bam ${each_individual}_${chromosome}_sorted
    samtools index ${each_individual}_${chromosome}_sorted.bam
done

```



Now we need to make the file executable, so type this:

`chmod +x alignment_commando`

And now we should be able to execute the file.  Type this:

`./alignment_commando`

Here we needed to preceed the name of our `bash` script by `./` to tell the computer where to find our script (i.e. in the current working directory).

Now please check whether this worked by checking out the file size of the files in your `~/samples` directory like this:

`ls -l samples`

You should see lots of file including, for each sample a `_sorted.bam` file that is not of file size zero.

If this worked, you can now delete the intermediate files like this:

```
rm -f samples/*_chrZZZ.bam
rm -f samples/*.sai
```


##  Practice Problem 3 (for home): Please try out combining some of the tricks that Brian told you about last class to modify the bash script to run the commands without listing the file names in the script.

Imagine you were using the script above to work with a complete genome alignment.  This could take some time and you want to have some idea of how much progress the script has made.  Can you please use the Unix `echo` command to the bash script to keep you informed about which command is executing.  You could (for example) ask the script to tell you when command 1, 2...5 is done and for which individual it has been completed.



## Please click [here](https://github.com/evansbenj/XXX) to go to the next page.

