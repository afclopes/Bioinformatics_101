# How to run bismark:

   ## 1. Export, to inside your directory, the path to bismark and bowtie2
   
> export PATH=/path/bowtie2-2.3.2/:$PATH

> export PATH=/path/Bismark_v0.19.0/:$PATH

> export PATH=/path/samtools-1.3/bin/:$PATH
   
Important note: please check that the versions of the packages that you will use are the msot up-to-date. If not, then update it. When way to do this is to create an 
anaconda environment and install your packages in there (supposedly these packages will be of the latest version).

   ## 2. Download the reference genome and unpack it:
   
> gunzip file.gz

   ## 3. Run bismark preparation:
  
https://www.bioinformatics.babraham.ac.uk/training/Methylation_Course/Basic%20BS-Seq%20processing%20Exercises.pdf
https://github.com/FelixKrueger/Bismark/tree/master/Docs
  
> bismark_genome_preparation --verbose /path/Genomes/ 
 
'/path/Genomes/' --> this is just the name of the folder containing the genome of reference
 
This step takes roughly 2 hours.
  
   ## 4.  Run bismark alignment: Make a sbatch script:
  
> #!/bin/bash

> #SBATCH -J Bismark_singleCell
> #SBATCH -t 48:00:00
> #SBATCH --output=Bismark_oocytes_SingleCell_Deeper_R1_run.out
> #SBATCH -c 1
> #SBATCH -N 1
> #SBATCH --partition=ServerPartitionName

> export PATH=/path/bowtie2-2.3.2/:$PATH

> export PATH=/path/Bismark_v0.19.0/:$PATH

> bismark --genome --non-directional /path/Genome_Downloads  /path/SRR1411188_val.R1.fastq.gz

This step takes about 15 hours to complete.

Note: if .fastq files are not used (this is the default), then specify using -f for fasta files.

Note: bismark can use samtools, in that case, it will produce data in .bam format (which you will not be able to access)

   ### 4.5. Run bismark alignment for bulk cells:
   
   In theory it is the same as for single cells. But it will take much longer. The 48h limit that was previously set will not be long enough. You could try 
   using the setting *--multicore* for example, *--multicore 4*, but it didnt fully run on the cluster. Next attempt could be to run with this *--multicore 4* 
   and change  *#SBATCH -N 4* including 4 instead of 1 for the -N parameter. This time, maybe it will run...but then again..this is trial and error...

   I ended up running this by increasing the time that the job is allowed to run for. I increased it to the maximum (324 h). In this case however, only 
   between 1h 22min and 2h 17min was necessary to complete the entire job.


  ## 5. Run de-duplication to remove PCR artifacts, and produce .bam files
  
     - First, unpack the .sam.gz which was the result of the alignment.
  
     > gunzip file.gz
  
     - This takes around 2mins.
  
     - Second, make sure samtools is installed. Maybe export the path?
     - Third, run de-duplication:
  
     > deduplicate_bismark --bam -s SRR1411188_val.R1_bismark_bt2.sam
      
     --bam: is used to make output files in .bam format
     -s: is to tell you are using single end files
     -p: is to tell you are using paired end files
  
     Run this step for each of the two files (if your bismark was run for single-ended files). Or once for data coming from a paired-end alignment pipeline.
      
     For further options for the deduplicate_bismark step, check out the general bismark settings, eg for bismark_methylation_extractor:                                  https://github.com/FelixKrueger/Bismark/tree/master/Docs
  
     This step takes around 6 mins.
  
  ## 6. Merge bam files into one using samtools:
  
 > samtools merge SRR1411188_deduplicated_merged.bam SRR1411188_val.R1_bismark_bt2.deduplicated.bam SRR1411188_val.R2_bismark_bt2.deduplicated.bam

 This step takes around 18 mins.
   
   ## 7. Sort the bam file using samtools:
   
 > samtools sort -o SRR1411188_deduplicated_merged_sorted.bam SRR1411188_deduplicated_merged.bam
   
 This takes around 4mins to complete.
   
   ## 8. Remove BAM flags for secondary alignments + make sure mapping quality is above 20:
   
   First, stop. Have a look at your file. What type of flags are there in your .sam file? Extract the unique flags. Try: 
   *cut -f2 SRR1411188_val.R1_bismark_bt2.sam | sort -u* Many different flags? Which ones? Worth trying to see where these flags are? Are there any 
   flags that you need to remove?
   
  > samtools view -b -F 0x100 -q 20 SRR1411188_deduplicated_merged_sorted.bam > SRR1411188_deduplicated_merged_sorted_noFlag.bam
  
  There seems to be no difference between using -F 256 or -F256 or -F 0x100
  
  Doing a quality control here using 'mapping quality' is a controversial topic, with some papers doing it, and often different from each other, and 
  others not doing it at all.
  
   ### 8.5 Removing multiple BAM flags:
  
  Interesting ones to consider:
  - 0x100 (or 256): secondary alignments
  - 0x4 (or 4): segment unmapped
  - 0x10 (or 16): segment being reverse complemented
  - 0x800 (or 2048): supplementary alignments
  
 To use these BAM flags with samtools, you can either use:
  - -F: dont give output
  - -f: give output
  
  For example: -f16 will give reverse complementary strands and -F16 will give NOT reverse complementary strands, therefore forward strand.
  
  You can combine -F and -f in the same command, for example using -F4 -f16 to find only mapped segments that are reverse complementary. You can also 
  combine two -F BAM flags like -F260 which would mean both -F4 and -F256 so would find only segments that are mapped and without secondary alignments.
  
  > samtools view -b -F 0x100 
  
  Interesting mixing and matching that I have tried before:
  - -F260 = 4 + 256
  - -F272 = 16 + 256
  - -F256 -f 16
  
  Remember: if you remove the secondary alignments, you need to also remove the primary alignments of those secondary alignments because this read is not 
  reliable. Try this by creating a list of all the IDs with secondary alignments, and use these unique IDs to remove primary and secondary alignments from 
  the data.
  
  Also remember: if you have supplementary alignments, you need to deal with them too. One way is to just get rid of them. Another way is to go through them 
  and select which ones start at the end of the sequence and end at the beginning of the sequence of mtDNA, and vice versa. But remove the reads that have half 
  of the sequence in one direction, and the other half in the other direction.
   
   ## 9. Index the bam file using samtools:
   
 > samtools index SRR1411188_deduplicated_merged_sorted.bam SRR1411188_deduplicated_merged_sorted.bai
 
 It will create a .bai, not a .bam.bai (like the MToolBox produces), but you can still use these .bam and the .bai files in igv.
   
Takes only a few seconds!
  
  ## 10. Transfer it to the computer:
  
Make sure to connect to the intranet first via VPN.
  
> rsync -av afcl2@suffolk:/work/WorkGenomicsB/afcl2/BSseq-Bismark/SRR1411188_deduplicated_merged_sorted* .
  
This takes around 10mins.
  
  ## 11. Visualise it in IGV: by loading the .bam file but having the .bai file in the same folder.
