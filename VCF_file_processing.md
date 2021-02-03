# How to process a vcf file.

It took many tries to manage to get the vcf file in the format that I wanted. I will present different attempts in the hope that we can learn different commands,
**then I will present the final version**.

Understanding better a vcf file: https://samtools.github.io/hts-specs/VCFv4.2.pdf

## Setting a header to the table:

> awk 'BEGIN{print "CHROM  POS     ID      REF     ALT     QUAL    FILTER   DP4"}{print $1,$2,$3,$4,$5,$6,$7}' SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split.vcf > SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new.vcf

## Spliting columns made of strings:

> awk -F';' '{print $5}' SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split.vcf | awk -F'=' '{print $2}' | paste SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new.vcf - > SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new2.vcf

> awk  -F',' '{print $1,$2,$3,$4}' SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new.vcf | paste SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new2.vcf - > SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new3.vcf

--> This did not work. It printed the first section of the table and a repetition of the first section of the table plus the second table with the new adaptation. 
Which means that the -F makes the splitting and prints what you asked for, but also will automatically print the other things that were in the table too.

- : means to add the results of the previous pipeline



## Remove line of the vcf file:

Use Vim to remove any line from the file

## Use a variable:

> awk   '/DP4=/ {dap=index($0,"DP4=");print $1,$2,$3,$4,$5,$6,$7,substr($0,dap+4,7)}' SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split.vcf > SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new.vcf

What this is doing:
1. '/DP4=/  find each line contains DP4=
2. dap=index($0,"DP4=") is telling to save inside the variable dap: the position where DP4= is found
3. substr($0,dap+4,7) refers to the part of the whole line ($0 is whole line) starting at position dap and continuing seven characters afterwards

## Use two variables:

> awk   '/DP4=/ {dap=index($0,"DP4=");posdap=index(substr($0,dap+4,256),";");print $1,$2,$3,$4,$5,$6,$7,substr($0,dap+4,posdap-1)}' table0.dat

What this is doing:
  1.'/DP4=/  find each line contains DP4=
  2. ;  do whats before ; and what is after ;
  3. posdap=index(substr($0,dap+4,256),";") is telling to save inside the variable posdap: the position of the substraction of variable dap (a position) plus 4 
  characters (which is DP4=, giving us a new position), minus 256 characters (better understood as up to 256 characters), ending in ;. So posdap will hold the 
  position of the first ';' after DP4=
  4. substr($0,dap+4,posdap-1)  is telling to substract the position of dap+4 minus the position of posdap-1. why -1? because you dont want to count the 
  character ';'. 

# Final version:

> awk 'BEGIN{print "CHROM  POS     ID      REF     ALT     QUAL    FILTER   DP4-Ref-Fwd  DP4-Ref-Rev  DP4-Alt-Fwd DP4-Alt-Rev"};/DP4=/ {dap=index($0,"DP4=");posdap=index(substr($0,dap+4,256),";");print $1,$2,$3,$4,$5,$6,$7,substr($0,dap+4,posdap-1)}' SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split.vcf > SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new.vcf

Here we have managed to achieve a table with DP4, but we still need to split DP4 into 4 individual columns.

> awk  -F',' '{print $1,$2,$3,$4}' SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new.vcf > SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new4.vcf

# Final version after bcftools call -m + bcftools norm -f:

With this new setting, sometimes more than one alternative allele appears in a given position. When this happens, the awk steps that i was doing can get confused. 
So here are the next steps to do instead:

>awk '/DP4=/ {dap=index($0,"DP4=");posdap=index(substr($0,dap+4,256),";");print $1,$2,$3,$4,$5,$6,$7,substr($0,dap+4,posdap-1)}' SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads.vcf > SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_part1_1.vcf

>awk '{print$1,$2,$3,$4,$5,$6,$7}' SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_part1_1.vcf > SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_part1_2.vcf


>awk '{print $8}' SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_part1_1.vcf | awk -F',' '{print $1,$2,$3,$4}' - > SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_part2_3.vcf


>paste SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_part1_2.vcf SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_part2_3.vcf | awk 'BEGIN{print "CHROM POS ID REF ALT QUAL FILTER DP4-Ref-Fwd DP4-Ref-Rev DP4-Alt-Fwd DP4-Alt-Rev"};{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' -  > SRR1411188_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_final.vcf

You can put this in one file to run on bash and run it as a script pipeline.

# General rules about awk:

  1. it will make changes without worrying about what comes next. So you can set a header for 8 columns but actually have 50 columns. Awk doesnt care, that is 
  your problem.
  2. you can separate different arguments for it to process by ';' as if you would use an 'and' in your sentence.
  3. 'index' is used to set a position, where you search a string for a value, and stop the search at the first instance of this value. index(search string, 
  stop value)
  4. 'substr' is used to make a substraction, you need to add information of where you will make the substraction, what is your start position, and for what 
  length of characters

