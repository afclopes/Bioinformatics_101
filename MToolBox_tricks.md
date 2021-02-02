# What to know when running MToolBox:

You will need 3 files:

## Create a file for running sbatch. Eg. Mtoolbox_sbatch_v2.sh

#!/bin/bash   <-- needs to be at the start of every script to run in bash

#SBATCH -J oocytes_SingleCell_run1   <--name of the job that will be visible to everyone and you when doing *squeue* 

#SBATCH -t 48:00:00     <-- time allowed for the job to run

#SBATCH --output=oocytes_singleCell_run1.out      <-- output file for the sbatch run. Can give you information if the sbatch didnt work

#SBATCH -c 1    <-- number of cores used in the CPU

#SBATCH -N 1    <-- number of nodes used. Max 3-4 per person to be used at at time

#SBATCH --partition=ServerPartitionName      <-- partition where the job will run. First check with *sinfo* which partition is most free

export PATH=/path/MToolBox/:$PATH       <-- necessary to run MToolBox. Depends on where the tool is saved

INPUT_FILE=oocytes_SingleCell_config_v2.sh        <-- name of configuration file which includes information needed for the MToolBox to run

MToolBox.sh -i $INPUT_FILE &> ${INPUT_FILE}.log     <-- how to run MToolBox




## Create a configuration file needed for MToolBox to run (need to reference this in the sbatch file!). Eg. oocytes_SingleCell_config_v2.sh

#!/bin/bash   <-- needs to be at the start of every script to run in bash

mtdb_fasta=chrM.fa      <-- file name of the mitochondrial reference fasta sequence

hg19_fasta=hg19RCRS.fa    <-- file name of the nuclear reference sequence

mtdb=chrM   <-- name of the mitochondrial gsnap database

humandb=hg19RCRS      <-- name of the nuclear (+mitochondrial) gsnap database

input_path=/InputDirectoryPath/   <-- specify the FULL PATH of the input directory. Default is the current working directory

output_name=oocytes_SingleCell_SRR1248455_output    <-- specify the FULL PATH of the output directory if you want it to be somewhere very different. Otherwise just name a folder where the output data will be added to which will be MADE by MToolBox and place inside the directory defined in the input_path.

list=oocytes_SingleCell_SRR1248455_v2.lst   <-- Specify the file name which is a list of files of the correct format to be used 

vcf_name=oocytes_SingleCell_SRR1248455    <-- Specifiy the name to assign to the VCF file

input_type=fastq    <-- Specify the input file format extension. [fasta | bam | sam | fastq ]. fastq can be also fastq.gz or fq.gz

ref=RCRS    <-- Specify the mitochondrial reference to be used for the mapping step with mapExome. [RCRS | RSRS]

UseMarkDuplicates=false   <-- Used to get rid of PCR and sequencing duplication. Used more for DNA rather then RNA data. But not useful for sparse data

UseIndelRealigner=false   <-- Removes false positives in complex regions which could have misalignments, eg. realigns insertions and deletions. Can set to true. But this parameter will disappear from newer version.

MitoExtraction=false    <-- Only need to use when using bam files

hf_min=0.2    <-- specificy the minumum level of heteroplasmy for the alternative allele(s)

hf_max=0.8    <-- specify the maximum level of heteroplasmy for the alternative allele(s)

minrd=5   <-- specify the minimum read depth per position, so that the base read can be considered real [DEFAULT is 5]

minqual=25    <-- specify the minimum quality score per allele [DEFAULT is 25]. Measured during Illumina sequencing with flourophore for each base. Quality score is measured during base calling.



## Create a list file with the sample data eg. paired or unpaired data. Eg. oocytes_SingleCell_SRR1248455_v2.lst
  
SRR1248455.R1.fastq.gz    <-- no file path

SRR1248455.R2.fastq.gz    <-- no file path

   - no need to write the path to these samples so long as they are in the same directory as the Mtoolbox file and the config file.
   - usually paired end sequencing will produce 2 data samples. But, if you run trimming on the samples first, and then use a threshold of which size of 
   sequencing fragments should be kept, this might produce "orfan" sequences, which then will be transferred into a new file called "unpaired.fastq". 
   MToolBox accepts all 3 of these sequences at a time (2 paired and 1 unpaired sequencing file from the same data sample).
   - if files have adaptors, you need to remove these before you run them. Remove by trimming with trim-galore
   - files after trimming must have the same prefix when using MToolBox, so prefix.R1.fastq.gz Many of the files from trim-galore come out as 
   prefix_1_val_1.fastq.gz and prefix_2_val_2.fastq.gz. You must change the first number to be the same and hence become part of the prefix, and must 
   change the second number to be R1 or R2.
   
## Make sure all files .sh are executable
   - chmod +x filename
   - confirm with ls -l
   - .lst files do not need to be executable
    
# Extras to know:
  - if you dont get a *.vcf* file, then it means that mtDNA heteroplasmies were not computed. But, it could be that the alignment and base calling for 
  every base in every fragment has happened. You can check for these in other file outputs.
  - you cannot run MToolBox on multiple nodes. But you can run a job faster with other tricks. See the next section.
  - short reads can lead to unspecific mapping
  - each node has around max 32 cores. Never exceed using 3-4 nodes
  - there is no limit of cores that you can use, but there is limit of nodes.
  - check how many CPUs per node there are: *sinfo -Nl*
  - how to look at all the directories in suffolk: >> ls /suffolk/
  - it is possible to run trimmed files, and if your sample has adaptors, it is recommended to remove them, but if the trimming causes your sample 
  to end up with a file extension of .fq.gz, then you must change this to .fastq.gz
  
  
# Output files:
  1. outmt.sam: mitochondrial DNA alignment
  2. OUT.sam: mapping of reads to the mtDNA but not to the nuclear DNA
  3. OUT2.pileup: file that shows every base position for which there was a sequencing read. 
      - Every . (forward strand) or , (reverse strand) indicates a sequencing read that matches the reference sequence. When there is no match to 
      reference sequence, then you can see the base which it was changed to.
      - “$”: If this is the last position covered by the read
      - If this is the first position covered by the read, a “^” character followed by the alignment's mapping quality encoded as an ASCII character. 
      For more information, check http://www.htslib.org/doc/samtools-mpileup.html and http://www.asciitable.com/
  4. OUT2-sorted.bam and .bam.bai (which is the index of the bam file) can be used in the IGV software to look at the coverage of the mtDNA. Both files 
  need to be put in the same folder, and only .bam needs to be opened in IGV.
 
  
# Increasing speed of job with MToolBox:
  1. increase the Phreds using with gsnap. How?
      - inside Mtoolbox_sbatch_v2.sh, change line from *MToolBox.sh -i $INPUT_FILE &> ${INPUT_FILE}.log* to *MToolBox.sh -i $INPUT_FILE -m "t10 &> ${INPUT_FILE}.log* It can be changed to t10 or t12.
      - inside Mtoolbox_sbatch_v2.sh, change line from *#SBATCH -c 1* to *#SBATCH -c 5*. this number should be roughly the number of phreds chosen (10) divided by 2 (number of cores).
  2. use a sbatch job array:
      - inside Mtoolbox_sbatch_v2.sh, change line from *INPUT_FILE=oocytes_SingleCell_config_v2.sh* to *INPUT_FILE=\`awk "NR==$SLURM_ARRAY_TASK_ID" file.txt\`*
      - this file.txt will be a list of all the config.sh files needed to run all of your samples. (Remember that each config.sh file will tell which 2 or 3 files to use for a single sample).
      - when you run sbatch on the commandline, instead of just doing *sbatch Mtoolbox_sbatch_v2.sh*, you must do *sbatch --array=1-10%2 Mtoolbox_sbatch_v2.sh*, where 1-10 refers to the 
      number of lines of config.sh files in file.txt (in this example there are 10 config.sh names in each line, so 10 lines), and where %2 refers to running 2 MToolBox jobs at a time 
      (so, not splitting a job on a single sample, but running 2 sample jobs in parallel).
      - trick to make the file.txt file on the command line: >> ls \*config.sh > file.txt
    
# Data storage:

All data produced from MToolBox were kept on the cluster only.
    
