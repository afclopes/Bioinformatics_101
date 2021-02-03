# SAM files

Opening a .sam file might seem daunting, considering the number of tabs and columns and numbers you will find in there. But if you know the structure, 
you understand it without issues.

When you open a file, you will find a few rows of of @SQ with SN and LN fields:
- @SQ: Reference sequence dictionary.  The order of @SQ lines defines the alignment sorting order.
  - SN: Reference  sequence  name.
  - LN: Reference sequence length.Range: [1,231−1]

You will be able to spot every chromosome, however you might also find extras like this: '@SQ     SN:KZ289080.1   LN:394982' These are the novel patches 
that were made to the current sequence that is available for the reference genome. Until a new reference genome is created and made available, then we will 
use the versions of the reference 
genome plus its available patches. 

  - @PG: Program.
  - CL: Command line.  UTF-8 encoding may be used.
  - ID: Program record identifier.  

After @PG, the different reads that were aligned will appear. This will be divided into a table that is tab delimited. The column, in order are:

    1. QNAME: Read Name
    2. FLAG: SAM flag --> might need decoding
    3. RNAME: contig name or * for unmapped
    4. POS: mapped position of base 1 of a read on the reference sequence
    5. MAPQ: mapping quality
    6. CIGAR: string describing insertions and deletions
    7. RNEXT: Name of mate
    8. PNEXT: Position of mate
    9. TLEN: Template length
    10. SEQ: Read Sequence
    11. QUAL: Read Quality
    12. TAGS: Additional information in TAG:TYPE:VALUE format
    
    FLAG options can be a combination of:
    Bit     Hex code    Description    
    1       0x1         template having multiple segments in sequencing
    2       0x2         each segment properly aligned according to the aligner
    4       0x4         segment unmapped
    8       0x8         next segment in the template unmapped
    16      0x10        SEQbeing reverse complemented
    32      0x20        SEQof the next segment in the template being reverse complemented
    64      0x40        the first segment in the template
    128     0x80        the last segment in the template
    256     0x100       secondary alignment
    512     0x200       not passing filters, such as platform/vendor quality controls
    1024    0x400       PCR or optical duplicate
    2048    0x800       supplementary alignment
    
    For each read/contig in a SAM file, it is required that one and only one line associated with theread satisfies ‘FLAG& 0x900 == 0’.  This line is 
    called the *primary line* of the read.


For more information, please check: https://samtools.github.io/hts-specs/SAMv1.pdf

## MAPQ: Mapping quality

Confusing and tricky topic. Read more on http://www.acgt.me/blog/tag/sam 3rd article.

Bottom line is that each different alignment tool creates a different mapping quality score, which is not well documented. So, first, you need to find out 
which alignment tool you used. Either you know this already, or, and actually regardless, you should confirm this. How can you confirm? Different tools have 
different ranges of mapping quality scores:
  - Bowtie2: 0-42
  - BWA: 0-37
  
Check your own data with *cut -f5 SRR1411188_val.R1_bismark_bt2.sam | sort -u* Are there scores higher than 37? If yes, then you used Bowtie2.

From the SAM specification: 

> MAPQ: MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality 
is not available.

So if you happened to know that the probability of correctly mapping some random read was 0.99, then the MAPQ score should be 20 (i.e. -10\*log10(0.01), 
remember 0.99 is the probability of mapping correctly, we want the probability of mapping incorrectly for this maths). If the probability of a correct match 
increased to 0.999, the MAPQ score would increase to 30. If the probability of a correct match decreased to 0.90, the MAPQ score would decrease to 10.

Check how many reads fall into each mapping quality category, using *cut -f5 SRR1411188_val.R1_bismark_bt2.sam | sort | uniq -c*

MAPQ is the abbreviation for mapping quality in the sam file. -q is the value referring to the mapping quality in bcftools.

## QUAL: Base quality

### Definition of bcftools:

QUAL is the abbreviation for base read quality in the same file. -Q is the value referring to the base read quality in bcftools.

Definition of QUAL: QUAL: ASCII of base QUALity plus 33 (same as the quality string in the Sanger FASTQ format).   
> A base quality is the phred-scaled base error probability which equals −10 log10Pr{base is wrong}.  This field can be a ‘\*’ when quality is not stored.

Some papers suggest to have 10 or 20 as a threshold. Eg. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005306#sec002 Other papers do 
not even mention base quality. In some data in sam files, it is not even possible to find QUAL values.

Another reference on base quality:
https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf

https://www.researchgate.net/post/How_exactly_are_QUAL_values_calculated_in_a_VCF_file1

### Definition of VCF files:

QUAL: Phred-scaled probability of all samples being homozygous reference.

I found an explanation of QUAL when produced by the GATK tool for variant calling: https://www.biostars.org/p/174075/

"The Phred-scaled probability that a REF/ALT polymorphism exists at this site given sequencing data. Because the Phred scale is -10 * log(1-p), a value of 10 
indicates a 1 in 10 chance of error, while a 100 indicates a 1 in 10^10 chance (see the FAQ article for a detailed explanation). These values can grow very 
large when a large amount of data is used for variant calling, so QUAL is not often a very useful property for evaluating the quality of a variant call."

## SAM Flags

SAM flags can be a combination of the flags described in the table above. It can be difficult to decipher the numbers when they dont fit perfectly into 
only one category. In these cases, it is best to find out which flags are there for example using *cut -f2 SRR1411188_val.R1_bismark_bt2_pe.sam | sort -u*, 
and then using the website https://broadinstitute.github.io/picard/explain-flags.html to input these flags and get the breakdown of different combined flags. 
For example:

SAM flag 147:
   read paired (0x1)
   read mapped in proper pair (0x2)
   read reverse strand (0x10)
   second in pair (0x80)
   
SAM flag 163:
    read paired (0x1)
    read mapped in proper pair (0x2)
    mate reverse strand (0x20)
    second in pair (0x80)
    
SAM flag 83:
    read paired (0x1)
    read mapped in proper pair (0x2)
    read reverse strand (0x10)
    first in pair (0x40)
    
SAM flag 99:
    read paired (0x1)
    read mapped in proper pair (0x2)
    mate reverse strand (0x20)
    first in pair (0x40)
