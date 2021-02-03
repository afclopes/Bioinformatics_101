# Accessing Reference Genomes

To access reference genomes often requires heavy downloading from pages like:
  - Sanger
  - USCS
  - Ensembl
  - GRC
  - EBI ENA
  
  Alternatively, if you are working on a cluster where other people have used any of the reference genomes before, it is likely that they have it already downloaded.
  This will save you a lot of time and effort. Just ask them for their permission to use their files and ask for the precise path to get to them.


# How to get reference genomes?

Different databases have different formats in how they present reference genomes. Overall websites that talk about this are:

https://hbctraining.github.io/Accessing_public_genomic_data/lessons/accessing_genome_reference_data_odyssey.html

https://genome.ucsc.edu/FAQ/FAQgenes.html#mito


It is possible to find the genomes, for example of mouse with the reference split up into the different chromosomes, like:

ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/

ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/

But often the types of tools that you will need to use require all the chromosomes to be together, and this can be found in places like:

ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/

https://www.gencodegenes.org/mouse/


# How to download the reference genome sequence to the cluster?

>> wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz

On the cluster, this download is done super fast.


## What is the haplogroup of the reference genome?

It appears that the mouse genome is based on the C57BL/6J strain.

https://www.ncbi.nlm.nih.gov/grc/mouse


One of the papers that I am currently working on works with the strain C57BL/6Babr, this strain name means that the mouse C57BL/6 that was used has been 
grown at Barbraham Institute, but it is likely that this strain has been at this institute before the strain subdivided.

## Where can we find reference genomes of other mouse strains?

Here you can find different strains, and you can download just one chromosome that you are interested or the whole genome.

https://www.ncbi.nlm.nih.gov/genome/browse#!/eukaryotes/52/

http://hgdownload.soe.ucsc.edu/hubs/mouseStrains/hubIndex.html

Note: file extensions .fa and .fna are both fasta files.

Here we can see what are the lastest versions:

https://www.ncbi.nlm.nih.gov/assembly/organism/10090/latest/

You can download the full genome from here:

https://www.ncbi.nlm.nih.gov/assembly/GCA_001632555.1/ 
OR https://www.ncbi.nlm.nih.gov/genome/52?genome_assembly_id=275024 
OR https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/632/555/GCA_001632555.1_C57BL_6NJ_v1/ 
OR https://www.ebi.ac.uk/ena/browser/view/GCA_001624185.1

Another option is to find the genome with UCSC, click on the mouse you want, then in Ensembl click on Whole Genome, and then click on Assembly, and download 
file from EBI ENA: http://www.ensembl.org/Mus_musculus_129S1_SvImJ/Location/Genome?db=core;r=4:136366473-136547301


# How to compare reference genomes?

Use Blast to be able to compare different reference genomes. You can add the accession number so you don't have to copy and paste the sequence in there (and risk 
human errors).

https://blast.ncbi.nlm.nih.gov/Blast.cgi?BLAST_SPEC=blast2seq&LINK_LOC=align2seq&PAGE_TYPE=BlastSearch

The accenssion numbers you can get from here: https://www.ncbi.nlm.nih.gov/genome/browse#!/eukaryotes/52/

Once the alignments have been done, you can go to Alignments, and select "Pairwise with dots for identities". In there you need to scroll through the entire 
sequence to spot the red font bases which are the differences between the two sequences.
