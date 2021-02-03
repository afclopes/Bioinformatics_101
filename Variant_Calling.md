# How to compute variant calling using samtools and bcftools:

Websites for more information:

http://www.htslib.org/doc/1.3/samtools.html

https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html

https://github.com/samtools/bcftools

https://github.com/samtools/bcftools/issues/945

Websites for information on other studies of sequencing with mitochondria:

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224847

https://github.com/GRC-CompGen/mitochondrial_seq_pipeline

### Things to keep in mind:
  - bcftools mpileup is the new samtools mpileup.
  - bcftools call is the new bcftools view

  ## 1. Create a sbatch script for variant calling - Using bcftool call -m:

> #!/bin/bash

> #SBATCH -J bcftools_samtools

> #SBATCH -t 48:00:00

> #SBATCH --output=Bismark_oocytes_SingleCell_Deeper_Samtools_bcftools_run.out

> #SBATCH -c 1

> #SBATCH -N 1

> #SBATCH --partition=ServerPartitionName

> export PATH=/path/bowtie2-2.3.2/:$PATH

> export PATH=/path/Bismark_v0.19.0/:$PATH

> export PATH=/path/samtools-1.3/bin/:$PATH

> export PATH=/path/bcftools-conda/bin/:$PATH

> bcftools mpileup -Ou --max-depth 800 -f /path/Genome_Downloads/GRCm38.p6.genome.fa SRR1411188_deduplicated_merged_sorted.bam | bcftools call -m -Ob -o SRR1411188_deduplicated_merged_sorted_bcftools_raw.bcf

Different settings that can be used for mpileup:
  - -B: Disable probabilistic realignment for the computation of base alignment quality (BAQ). BAQ is the Phred-scaled probability of a 
  read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments. 
  - -q: Minimum mapping quality for an alignment to be used 
  - -Q: Minimum base quality for a base to be considered 
  - -O b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between 
  bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion. 
  - Check the bcftools mpileup --max-depth option, most likely it should be increased. Note that by default only 250 reads per-file are considered at a position!
  - -o: necessary if you want the output to be put into a file.

Paper looking into mitochondrial sequencing from data from NGS sequencing: 
https://github.com/GRC-CompGen/mitochondrial_seq_pipeline and https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224847

Different settings that can be used for call:
  - -m: multiallelic
  - -ploidy 1: mitochondria is haplotype
  - 

 ## 2. Create a sbatch script for *viewing* the variant calling:
  
> #!/bin/bash

> #SBATCH -J bcftools_view
> #SBATCH -t 48:00:00
> #SBATCH --output=Bismark_oocytes_SingleCell_Deeper_Bcftools_view_run.out
> #SBATCH -c 1
> #SBATCH -N 1
> #SBATCH --partition=ServerPartitionName

> export PATH=/path/bowtie2-2.3.2/:$PATH

> export PATH=/path/Bismark_v0.19.0/:$PATH

> export PATH=/path/samtools-1.3/bin/:$PATH

>export PATH=/path/bcftools-conda/bin/:$PATH

> bcftools view -Ov SRR1411188_deduplicated_merged_sorted_bcftools_raw.bcf

## 3. Select only the mitochondria data from this big file:

Output came out as the .out file! A little strange...If I do 'cat' of the file, it comes out as a massive file. Maybe I can select just the mitochondrial 
chromosome to look at: 

> cat Bismark_oocytes_SingleCell_Deeper_Bcftools_view_run.out | grep -w chrM > Bismark_oocytes_SingleCell_Deeper_Bcftools_view_run_MitoReads.out

To be able to read this file, we need to know the following information:

##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">

##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">

##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">

##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">

##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">

##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">

##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">

##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">

##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">

##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">

##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">

##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">

##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">

##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">

##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">

##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">

##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">

##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">

##bcftools_callVersion=1.9+htslib-1.9

##bcftools_callCommand=call -m -Ob -o SRR1411188_deduplicated_merged_sorted_bcftools_raw.bcf; Date=Tue Sep 29 16:20:01 2020

##bcftools_viewVersion=1.9+htslib-1.9

##bcftools_viewCommand=view -Ov SRR1411188_deduplicated_merged_sorted_bcftools_raw.bcf; Date=Thu Oct  1 12:43:19 2020

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SRR1411188_deduplicated_merged_sorted.bam


*Notice that the DP4 term is what we can use for heteroplasmy calling. By creating a ratio of ref-forward + ref-reverse : alt-forward + alt-reverse.*

 ## 4. Create a sbatch script for filtering:
  
  Read more on which filtering to do here: https://samtools.github.io/bcftools/howtos/variant-calling.html


## 5. Using bcftools call -c:

This is useful to see what different alleles for ALT there are. The metadata for a file generated this way will be as such:

##ALT=<ID=*,Description="Represents allele(s) other than observed.">

##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">

##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">

##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">

##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">

##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">

##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">

##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">

##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">

##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">

##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">

##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">

##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">

##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">

##INFO=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">

##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">

##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">

##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">

##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">

##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">

##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">

##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">



### Understanding AF1 and AC1:

AF1=Max-likelihood estimate of the first ALT allele frequency (assuming HWE).

HWE (Hardy–Weinberg equilibrium), see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/ for more details, is a principle stating that the genetic 
variation in a population will remain constant from one generation to the next in the absence of disturbing factors. For instance, mutations disrupt the 
equilibrium of allele frequencies by introducing new alleles into a population. The Hardy-Weinberg equilibrium can be disturbed by a number of forces, 
including mutations, natural selection, nonrandom mating, genetic drift, and gene flow. Because all of these disruptive forces commonly occur in nature, 
the Hardy-Weinberg equilibrium rarely applies in reality. Therefore, the Hardy-Weinberg equilibrium describes an idealized state, and genetic variations 
in nature can be measured as changes from this equilibrium state.

Therefore, I assume that AF1 will define how strongly the first ALT allele appears, as compared to the reference base and the other ALT base. For example, 
if there is no ALT base, then AF1=0. If there is 1 ALT base, then AF1 is around 0.5. If there is 2 ALT bases, and AF1 is 1, it means that the likelihood is that 
almost all bases will be the first ALT base.

AC1=Max-likelihood estimate of the first ALT allele count (no HWE assumption). Therefore, I assume that AC1 will give the number of reads for the first ALT allele.


# INDEL:

In a few instances, I notice the strange case of having in the Ref base something along the lines of CAAAA or TTAT, and in the ALT position there was nothing. 
These rows were recorded as INDEL. According to this source (https://samtools.github.io/hts-specs/VCFv4.2.pdf and 
https://genome.sph.umich.edu/wiki/Variant_classification and https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele- and https://www.biostars.org/p/387479/),
it is likely that those positions are deletions.


### Picking between bcftools call -c -m, -C:

All of these options are slight variations of the next.

-C calls genotypes when given alleles. You need to give a vcf file which will be used to make a list of alleles, which then is loaded together with your bam 
file to be able to call the genotypes. --> we dont have a vcf file like that.

-c calls variants that are biallelic. Mitochondrial DNA is not biallelic. Having said that, many other tools for calling variants uses the biallelic 
assumption (eg. MToolBox). (https://sourceforge.net/p/samtools/mailman/samtools-help/thread/5AF6BC62-7EF4-43BB-B90D-149C21FA0ECF%40cmmt.ubc.ca/#msg32931405)

-m calls variants that are multiallelic. This is the case for mitochondria.

If I use -c, it will give me a number of positions with different ALT alleles. If I use -m, no positions with different ALT alleles will be given. I find it 
strange that the multiple ALT alleles shows up with the biallelic setting. But in the end, only around 4 positions with -c showed this feature. So I have decided
to go with -m and not dwell on only these 4 positions.

# Ploidy

Mitochondria is not diploid, is it haploid. Hence when using these tools, sometimes it is possible to make a specific setting for haploid. Still, many tools 
do not make this specification and let the measures be carried out assuming diploidy (eg. MToolBox). Having that said, it depends on the type of analysis that 
you will make from the results. If for example, like in our case, we want to measure heteroplasmy, then the values that really matter to us are the reads. The 
reads will not be affected by the ploidy, therefore setting this particular parameter should not matter. Moreover, there were some occasions where adding the 
setting for haploidy, the C>T mutations disappeared which is quite concerning since bisulfite sequencing should cause them.

# Normalisation

It is common practice to left-realign the sequences that have been aligned to the reference genome. Further explanation for this is found here 
https://genome.sph.umich.edu/wiki/Variant_Normalization and here https://databricks.com/blog/2019/12/05/streamlining-variant-normalization-on-large-genomic-datasets-with-glow.html. 
In fact, every time that I ran the normalisation, the INDEL disappeared from my data, which I see as a good sign of being thorough with my data. Moreover, 
other tools also do not use normalisation, once again, one example for this is MToolBox.

Interestingly, in the case of the positions of the INDELS, there often appeared another same position without an INDEL. Aftern normalisation, in the positions 
where there used to be an INDEL in one of the cases, now there appears only a mutation, and hence there are two same positions with mutations (often the same 
mutation) in the data. Luckily, only 2 INDELS appeared in my data, suggestive that this should not be a case for concern in the future.

Normalisation is useful particularly when there are complex homopolymer areas (see link for explanation https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2291785/). 
I believe that the loss of the INDELS that were previously seen could be due to one of these areas, since for example one of the INDELS was CAAAA.
