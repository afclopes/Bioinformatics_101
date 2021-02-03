# Pipeline design

You can design a pipeline that involves one step after the other as seperate steps. In this case if you can find a pattern in the names of the files that will be 
produced, you can add that in your pipeline. See examples later on marked as ex1.

You can also design a pipeline that takes the product of the previous task, which has not been written into a file (this is called a standard output), and runs 
that product in the next task. See examples later on marked as ex2.

Script:

#!/bin/bash

dir="thisRAWInputDir"
file=$1

awk 'BEGIN {print "!!!Next step starts here: renaming files!!!"}'

mv "${dir}/${file}"_1_val_1.fq.gz "${dir}/${file}"_val.R1.fastq.gz

mv "${dir}/${file}"_2_val_2.fq.gz "${dir}/${file}"_val.R2.fastq.gz

awk 'BEGIN {print "!!!Next step starts here: aligning sequences!!!"}'

#\ex1

bismark --non-directional --genome /path/Genome_C57BL6NJ_full -1 "${dir}/${file}"_val.R1.fastq.gz -2 "${dir}/${file}"_val.R2.fastq.gz

awk 'BEGIN {print "!!!Next step starts here: deduplication!!!"}'

#\ex1

deduplicate_bismark --bam -s "${file}"_val.R1_bismark_bt2_pe.bam 

awk 'BEGIN {print "!!!Next step starts here: sorting!!!"}'

#\ex1

samtools sort -o "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted.bam "${file}"_val.R1_bismark_bt2_pe.deduplicated.bam

awk 'BEGIN {print "!!!Next step starts here: mapping quality control!!!"}'

#\ex1

samtools view -b -q 20 "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted.bam > "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20.bam

awk 'BEGIN {print "!!!Next step starts here: pileup, variant calling!!!"}'

#\ex2

bcftools mpileup -B -Ou --max-depth 1000 -Q 30 -f /path/Genome_C57BL6NJ_full/GCA_001632555.1_C57BL_6NJ_v1_genomic.fa "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20.bam | bcftools call -m -Ov -o "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools.vcf

awk 'BEGIN {print "!!!Next step starts here: selecting mito reads from vcf file!!!"}'

#\ex2
cat "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools.vcf | grep -w chrM > "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_MitoReads.vcf

awk 'BEGIN {print "!!!Next step starts here: producing table with only desirable columns!!!"}'

#\ex2

awk 'BEGIN{print "CHROM POS ID REF ALT QUAL FILTER DP4-Ref-Fwd DP4-Ref-Rev DP4-Alt-Fwd DP4-Alt-Rev"};/DP4=/ {dap=index($0,"DP4=");posdap=index(substr($0,dap+4,256),";");print $1,$2,$3,$4,$5,$6,$7,substr($0,dap+4,posdap-1)}' "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_MitoReads.vcf | awk -F',' '{print $1,$2,$3,$4}' - > "${file}"_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_MitoReads_final.vcf


# How to check if the pipeline is working?

  1. check with *squeue* if the pipeline has managed to run for more than 5mins
  2. go into the folder that you are supposed to be creating new files. Are files being created? Compare these files to the script to guesstimate which step your 
  pipeline is at
  3. tail -f x.out: shows you the last sentence of the output. To exit Ctrl+C
  4. less x.out: See file. If press **F**, then *less* will go to the last line of the output and it will keep updating. If you press **Crtl+C**, then you will 
  go back just *less* where you can scroll up and down the output. If you press **F** again, then you go back to following the lasted output line of that file.

# Things to remember about pipelines:

  1. the server you are in can affect the pipeline submissions. For example, it could be that in one server, all the packages/toold that you need are available. 
  But in the other server that is not the case.
  2. it is best practice to have StepIDs throughout the pipeline to be able to spot the errors and which step the pipeline is currently running.
