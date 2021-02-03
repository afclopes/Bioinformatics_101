When measuring DNA methylation, the alignment steps including deduplication are the same for any data analysis looking at data produced with bisulfite treatment. It is often 
the later steps in the analysis which defines which results you are interested in getting.

For calling DNA methylation, use bismark_methylation_extraction.

Parameters to consider using for bismark extraction:

--no_overlap: specifying this option will extract the methylation calls of overlapping parts in the middle of paired-end reads only once (using the calls from the 
first read which is presumably the one with a lowest error rate). Deafult: ON.

--include_overlap: For paired-end data all methylation calls will be extracted irrespective of whether they overlap or not. Default: OFF.

--CX: extended to cytosines in any sequence context. This takes a long time to run.

--comprehensive: Specifying this option will merge all four possible strand-specific methylation info into context-dependent output files. The default contexts are:
(i) CpG context (ii) CHG context (iii) CHH context (Depending on the C content of the Bismark result file, the output file size might reach 10-30GB!).

--merge_non_CpG: This will produce two output files (in --comprehensive mode) or eight strand- specific output files (default) for Cs in (i) CpG context and 
(ii) any non-CpG context. The file sizes can also very very big for these.

--samtools_path: The path to your Samtools installation, e.g. /home/user/samtools/. Does not need to be specified explicitly if Samtools is in the PATH already.

--multicore: Sets the number of cores to be used for the methylation extraction process. If system resources are plentiful this is a viable option to speed up 
the extraction process (we observed a near linear speed increase for up to 10 cores used). Please note that a typical process of extracting a BAM file and writing 
out '.gz' output streams will in fact use ~3 cores per value of --multicore <int> specified (1 for the methylation extractor
itself, 1 for a Samtools stream, 1 for GZIP stream), so --multicore 10 is likely to use around 30 cores of system resources. This option has no bearing on the speed 
of the bismark2bedGraph or genome-wide cytosine report processes. Perhaps try multicore 4?

--report: Prints out a short methylation summary and the parameters used to run this script. Default: ON.

Parameters to consider when using bedGraph:

--bedGraph --counts

--cutoff: The minimum number of times a methylation state has to be seen for that nucleotide before its methylation percentage is reported. Default: 1 
(i.e. all covered cytosines). Try 5 to match the minimum of 5 reads for heteroplasmy.

--remove_spaces: Replaces whitespaces in the sequence ID field with underscores to allow sorting.

--buffer_size <string>: This allows you to specify the main memory sort buffer when sorting the methylation information. Either specify a percentage of physical 
memory by appending % (e.g. --buffer_size 50%) or a multiple of 1024 bytes, e.g. 'K' multiplies by 1024, 'M' by 1048576 and so on for 'T' etc. 
(e.g. --buffer_size 20G). For more information on sort, type 'info sort' on a command line. Defaults to 2G.



Other packages to try: 

  1. bismark2bedGraph
  2. bedGraph2cytosine:
  
--cytosine_report: produces a genome-wide methylation report for all cytosines in the genome. The output considers all Cs on both forward and reverse strands 
and reports their position, strand, trinucleotide content and methylation state.

--CX: The output file contains information on every single cytosine in the genome irrespective of its context. This applies to both forward and reverse strands. 
Please be aware that this will generate output files with > 1.1 billion lines for a mammalian genome such as human or mouse. Default: OFF (i.e. Default = CpG 
context only).

--split_by_chromosome: Writes the output into individual files for each chromosome instead of a single output file. Files will be named to include the input 
filename and the chromosome number.

  3. coverage2cytosine: every cytosine on both the top and bottom strands will be considered irrespective of whether they were actually covered by any reads 
  in the experiment or not.

Typical example:

> bismark_methylation_extractor -p --no_overlap --samtools_path <path_to_samtools> --multicore 4 --bedGraph --counts --buffer_size 10G --report s_1_sequence_bismark_pe.sam(or .bam)

Example including the genome-wide cytosine methylation report:

> bismark_methylation_extractor -p --no_overlap --samtools_path <path_to_samtools> --multicore 4 --bedGraph --counts --buffer_size 10G --report --cytosine_report --genome_folder <path_to_genome_folder> s_1_sequence_bismark_pe.sam(or .bam)


# Types of data produced:

  1. file.bedGraph
  This file can be read in IGV.
  
  2. file.cov.gz
  This file can be read with SeqMonk.
