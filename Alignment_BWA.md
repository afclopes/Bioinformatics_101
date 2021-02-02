# How to do alignment of sequencing data to reference genome using BWA?

  1. First, you need to prepare an index for bwa:

> bwa index ref.fa

   2. Then you need to run bwa:

> bwa mem ref.fa read1.fq read2.fq > Sample-name.sam

  3. To be able to later visualise the aligned sequences in for example, igv, you need to convert the produced .sam file into a .bam file:

> samtools view -Sbu ${SAMPLE_NAME}.sam | samtools sort > ${SAMPLE_NAME}.bam

Where -Sbu stands for:
S = automatically finds the sam files
B = output as bam file
u = uncompressed bam file produced

  4. Then you need to index this bam file
  
> samtools index test_sorted.bam test_sorted.bai

  5. Then you can view the .bam using the .bai file in igv

# Transfering files from the cluster to the computer:

There are two ways: either using the software which you install in your computer called Cyberduck. Or directly from the cluster using rsync.

Here is an example of how to use rsync:

> rsync -av username@gateway:/path/SRR1411188_deduplicated_merged_sorted* .

av = archive

For this, you need to be connected to the intranet. Use VPN to connect.
