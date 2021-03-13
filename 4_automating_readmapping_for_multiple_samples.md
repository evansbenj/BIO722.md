
(Or you can go back to the page about preparing your reference genome [here](https://github.com/evansbenj/BIO722.md/blob/main/3_read_mapping.md).

# Automating readmapping for multiple samples with a bash script

Now that you have seen how to align data from one individual to a reference genome, we can automate the alignment of all individuals to the reference genome using a bash script. This is easier than going through all that stuff independently for each individual. We can accomplish this with a `bash` script by defining a vector that contains the names of all of the individuals in the analysis, and then loop through this vector and execute each of the commands for each individual.

Below is an example `bash` script that should run all of our analyses for each individual.  Please use a text editor to make this program.  In the beginning of the script 5 variables are defined that specify, respectively, the path for the bwa and samtools programs, the path to the data, the path to the reference chromosome, and the name of the chromosome you are working on. You will need to modify the variables to have your USERNAME in the paths. For example, change the `USERNAME` in the `path_to_data="/2/scratch/USERNAME/monkey_directory/demultiplexed_subsetted_fq"` to whatever your username is  (e.g. `gradstd13`) and the same for the `path_to_chromosome` variable.

I have set the chromosome to take a argument that we will pass to the script when it is executed.  It is a good idea to do this for computecanada SLURM scripts too...

```
#!/bin/bash                                                                                            

path_to_data="/2/scratch/evanslab/monkey_directory/demultiplexed_subsetted_fq/"
path_to_chromosome="/2/scratch/evanslab/monkey_directory/my_monkey_chromosome/"
chromosome=${1}

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
    # echo which individual we are on
    echo ${each_individual}

    # align the data to the reference gemome
    bwa mem -t 16 -R "@RG\tID:FLOWCELL1.LANE6\tSM:${each_individual}\tPL:illumina" ../my_monkey_chromosome/${
chromosome}.fa ../demultiplexed_subsetted_fq/${each_individual}.fastq.gz | samtools view -bSh - > ${each_indi
vidual}_${chromosome}.bam

    # sort the bam file
    samtools sort ${each_individual}_$chromosome.bam -o ${each_individual}_${chromosome}_sorted.bam
    
    # index the sorted bam file
    samtools index ${each_individual}_${chromosome}_sorted.bam
    
    # clean up unsorted file
    rm ${each_individual}_${chromosome}.bam
done

```



Now we need to make the file executable, so type this:

`chmod +x alignment_commando`

And now we should be able to execute the file.  Type this:

`./alignment_commando chrZZZ`

Here we needed to preceed the name of our `bash` script by `./` to tell the computer where to find our script (i.e. in the current working directory). Also, instead of `chrZZZ` you should add whatever your chromosome is (e.g. `chr2` or `chr9`).

Now please check whether this worked by checking out the file size of the `.bam` files in your directory like this:

`ls -l`

You should see lots of file including, for each sample a `_sorted.bam` file that is not of file size zero.  Each file should have an index (`.bai`) as well.




##  Practice Problem 3 (for home): Please try out combining some of the tricks that Brian told you about last class to modify the bash script to run the commands without listing the file names in the script.

Imagine you were using the script above to work with a complete genome alignment.  This could take some time and you want to have some idea of how much progress the script has made.  Can you please use the Unix `echo` command to the bash script to keep you informed about which command is executing.  You could (for example) ask the script to tell you when command 1, 2...5 is done and for which individual it has been completed.



## Please click [here](https://github.com/evansbenj/XXX) to go to the next page.

