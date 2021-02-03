# Changing colours on shell terminal:

>> sh -c "$(wget https://raw.githubusercontent.com/ohmybash/oh-my-bash/master/tools/install.sh -O -)"

You can go in this website https://github.com/ohmybash/oh-my-bash/wiki/Themes and select which themes you like the best

Some times this new type of shell can affect the use of some of your packages inside bash, for this do:

>> source .bashrc.pre-oh-my-bash

# Running program on the commandline:

The program will likely start with '#!/bin/bash' for a bash program or '#!/usr/bin/env python3.6' for a python3.6 program. In order to be able to run these 
programs on the commandline, without submitting it to the cluster you must call the program file together with the name of the language it was 
written in. Eg.

> bash filename.sh
> python filename.py

# Commands worth remembering:

  - ps : see which processes are running in the background
  - hostname : useful to see in which computer/cluster/relay that you are currently in
  - pwd : see which directory you are curently in
  - locate : find a file or package that you know is somewhere in the computer/cluster
  - module :
  - wget : allows you to get files from web pages
  - file .bashrc : 
  - ssh : allows you to connect to another computer or cluster. You need to use >> ssh username@address
  - * : use it as a wild card. So \*.fastq means that every file before \* should be considered. val_1_* means that every file starting with val_1_ should be 
  considered. \*val\* means that every file containing val in its name should be considered
  - rm: remove
  - rm -R or rm -r or rm --recursive: remove the directory and all its files inside
  - rmdir: removes empty directories
  - cat: shows whats writen inside the file
  -  ls -l: result is -rw-r--r--. : first - would be d if its a directory, then we have r (read) w (write) x (execute). Then these rwx repeat. First trio 
  is for the user, then group, then others.
  - chmod +x filename : changes the permissions of the file to executable
  - #!/bin/bash : always at the top of the file, when you want it to be run in bash
  - less: shows your file:
      - /word: allows you to search for a work or term
      - g: go back to start of the file. every time you do a new search, you must go back to the start of the file, otherwise the search will start 
      wherever you are in the file at that moment
  - du: reports on disk space usage
      - du -h: reports on disk space in human readable form
      - du -sh: reports on disk space in human readable form of that directory
      - du -g or -m: gives space in gigabytes or megabytes
      - du -d 1: gives size in the first level of the directory
      - *Remember: 1 terabyte=1000 gigabytes; 1 gigabyte = 1000 megabytes; 1 megabyte= 1000 kilobytes*
  - cp and rm can be used to copy and delete more than one file at a time, wg. cp *sh *lst /path/to/directory
  - cp -R : copies directory and all its files inside
  - mv: useful to move files but also to rename files
  - ls \*endfoldername: allows us to access whats inside the folder by referring just to the end of the name of the folder
  - head -30 filename: see the top 30 lines of your file 
  - > : assigns what you did before to a new file
  - ln -s: makes a link to another folder and its files (its like a virtual passage to a far away folder where you can see and use files in that other folder)
    
# Writing on commandline:
   - Ctrl + a: go to start of the line
   - Ctrl + e: go to end of the line
   - writing file name: ls, click on the file, Cmd + v
   - Cmd + arrow: move between terminal tabs
   - Ctrl + arrow: move between desktops (not really to do with terminal, but super useful)
   - Ctrl + k: delete from cursor to end of line
   - Ctrl + u: delete from cursor to beginning of line
   
# Working in 'less':   
   - Shift +G: end of file
   - g: top of file
   - w: backwards one window
   - z: forwards one window
   - /: search term
   - n: repeat search forwards
   - N: repeat search backwards
   

# Working with multiple files:
   - making files of lists: *ls SRR1248482\* > ESC_SingleCell_SRR1248482.lst*
   - copying and renaming: *cp ESC_SingleCell_Multiple_SRR1248478_config.sh ESC_SingleCell_Multiple_SRR1248479_config.sh*
   

### Renaming files:
   - *rename new-file-name-section old-file-name-section old-file-name-complete* -> allows you to change just part of your file name
   - *rename new-file-name-section old-file-name-section old-\** -> allows you to change all the files which start with 'old-'
   - *rename new-file-name-section old-file-name-section \*-old* -> allows you to change all the files which end with '-old'
   
### Specific commands: awk, wc, cut, sort, uniq, grep

  Sometimes it is useful to use commands on bash that often come installed as a common basis.
  
  #### 1. awk
  
  awk is useful to find specific regions/terms in your file:
  > awk '$2 >= 6520 && $2 <= 11000'
  
  --> will look for the terms '6520' and '11000' in column 2 in both cases.
  
  > awk '$4 ~ /G/' filename.out > filename_Gs.tsv
  
  --> selects just rows where column 4 has a G.
  
  Note, there are special differences when searching for '.' .
  > $5~/./ searches in column 5 for just one character, any character 
  
  > $5~/\\./ searches in column 5 for exactly a dot
  
   *An example of awk: awk '$2 >= 6520 && $2 <= 11000' filename > newfilename*.
   
   You can also do this example without making a new file with the results, eg *awk '$4 ~/A/ && $5 ~/G/' filename.vcf*
    
   > awk '{print $4,$5}' SRR1411189_val.R1_bismark_bt2_pe.deduplicated_sorted_MAPQ20_bcftools_Normalised_MitoReads_final.vcf | sort | uniq -c
   
   --> does very similar to command seen later using 'cut'. Some times 'cut' gets confused about what are columns depending on how you created your data/table.
   
   
  #### 2. wc
  
  Allows you to count words, and by default it will actually count lines, words, and bytes. But you can tell it to just give you the count of lines: > wc -l
  
  
  #### 3. cut
  
  Allows you to select just a part of the data that you are interested in.
  
  *Example using cut in more than one column: cut -f4,5 SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_MitoReads.vcf | sort | uniq -c*
  
  Note: sometimes your data file will have line numbers (which you didnt include), and this will confuse the 'cut' function. So, do a few tests to check what is 
  in different columns and if it matches what you are expecting to be able to confirm where 'cut' starts counting your columns from.
  
  #### 4. sort
  
  Allows you to sort your data, for example to look for unique values
  
  *Example using cut and sort: cut -f2 SRR1411188_val.R1_bismark_bt2.sam | sort -u*
  
  #### 5. uniq
  
  Allows you to look only at the unique terms, and then do extra things with that, for example count them.
  
  *Example using cut, sort and uniq: cut -f5 SRR1411188_val.R1_bismark_bt2.sam | sort | uniq -c*
  
  *Another example. Sometimes cut gets confused with spaces and tabs in tables... in this case you can use awk instead of cut: 
  awk '{print int($6)}' SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_Q30_MitoReads_Split_new.vcf | sort -n | uniq -c*
  
  #### 6. grep
  
  Allows you to select things from a file.
  
  *Example: grep "^>" file*
  
  In this example, ^ is used to search only at the beginning of sentences. You can replace this by $ to search only at the end of sentences. Y
  ou can also not include ^ or $ and get all instances of your search in the entire file regardless of where in the sentence it can be found.
  
  #### 7. head and tail
  
  Allows you to see only the top or the end of your file.
  
  *Example: head -n 100 filename*  Where -n defines how many lines of the file you want to see.
  
  *Example: tail -n 100 filename*
  
  #### 8. sed
  
  Useful to edit a file without opening it (similar to awk, grep...). It is particularly useful when trying to do things like: replacing, deleting, 
  adding new string or patterns.
  
  # Stopping a job and resuming it on the commandline
  
  On the commandline:
  - jobs : to find out which job is running on the terminal
  - fg %[job number]: to resume job. Likely you will have paused only one or two jobs, so this job number will be very small and found in square 
  brackets, eg. [1]

  # Checking version of software
  
  software_name --version
  
  Note: if you are using a software that is in another path, remember to export to that path first, then to investigate the version of the software.
