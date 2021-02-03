## Write a script

1. Use *vim* to write a python script.
2. Save the script as .py
3. Run the script on the commandline with 
  > python filename.py

## Have a look at your data

Questions you need to ask yourself:
  1. What type of data clean-ups will you have to do?
  2. What columns/rows are you interested in?
  
Open your data file with *less -N filename* and find out which lines you will have to ignore for your python script.


## Options to tackle vcf file via python: using pyvcf:
### Install pyvcf:

  1. install pyvcf
  
> conda install -c conda-forge pyvcf

  2. Open python

> ~/miniconda2/bin/python

  3. Working with pyvcf

https://pyvcf.readthedocs.io/en/latest/_modules/vcf/parser.html

> import vcf

>Paired_MAPQ20_data = vcf.Reader(open('/path/SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools.vcf','r'))

> for record in Paired_MAPQ20_data:

>     print record

--> this way we can check that the data was indeed read by pyvcf

> Paired_MAPQ20_data.formats

--> see what kind of columns are saved under 'formats'

> Paired_MAPQ20_data.infos

--> see what kind of columns are saved under 'infos'. Found the column you were looking for? Confirm that it is indeed what you think:

> Paired_MAPQ20_data.infos['DP4']

--> select only data on the mitochondrial DNA: first need to install pysam

> Paired_MAPQ21_data_ChrM = Paired_MAPQ20_data.fetch("chrM")

....I started to struggle, so I found another tool to use: VCF.py:


## Option to tackle vcf file via python: using VCF.py

### Working with VCF.py

https://gist.github.com/slowkow/6215557

Firstly, write a python file called "VCF.py" and save it to the directory where you will open python. Once in python, do:

df_INFO=VCF.dataframe('/path/SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools.vcf',large=False)

--> save as a data frame on python? possible?


## Another interesting method to tackle vcf file in python:

https://www.biostars.org/p/82053/

Try perl script? or run as described in .split() in python.


### Always make sure your tools are installed before you start:

> conda install pandas

# Loading table in python:

I ended up using awk to fix my table before loading it onto python.

> import pandas as pd

> PE_MAPQ20_table=pd.read_table("/path/SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads_Split_new4.vcf",sep="\s+",header='infer')



# Important commands to remember for python:

  1. check data type:
  > table.dtypes
  
  2. check data dimensions:
  > table.shape
  
  3. if you send a command, and python replies with ... it means it is waiting for you to complete the command (you likely forgot to close a bracket for example)
  
  
# Visualising plots made on python:

> bash ~/imgcat.sh Coverage_Plot_CT.png

## Visualising results

To visualise results, for example including the plots seen with imgcat.sh in an .out file, it is best to use 'cat' instead of 'less' as less does not show the 
images.

## Processing fasta files with python

Use BioPython package. Once installed, access it in python with > from Bio import SeqIO


# To work on python and create images that you can visualise, you must:

  1. Log onto the server via ssh -X @gateway, then ssh -X @lambda
  2. You must use the 'terminal' which is associated to the XQuartz
  3. Trying running step by step on my python first (call python with typing: *python3.6*) eg. use script for PythonScript_HetTable_AND_Plots_Sbatch_v2.py
  4. Trick: if there is an error about display, but you are only trying to produce a plot to be saved into a file and you are inside python, try exiting the cluster and entering again. 
  Then try making the plot again from inside python.

