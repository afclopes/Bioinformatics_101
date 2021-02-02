# To download a file to the cluster:

Make sure you are logged into the correct server. Some will not allow you to do downloads, while other will allow it. It could be that your directories are 
mirrored between servers, so you might have all your files in both servers, so it is just a matter of logging into the server that allows you to have access 
to the outside net and running the commandline from there. Often you cannot run an sbatch job for downloading files, you need to do it yourself from the 
commandline, and make sure that the terminal doesnt log you out mid download due to inactivity from your side.

Use wget to download a single file:
  1. Go to the website you want to download the file from.
  2. Right click on the file > Copy link location
  3. On command-line: >> wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/005/SRR1248455/SRR1248455_2.fastq.gz
  
Use wget to download all files from a page:
  1. Go to the website you want to download the file from.
  2. Right click on the file > Copy link location
  3. wget -r --level=2 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/006/SRR1248456
      - '-r' tells the terminal to download all of the files and and directories associated to that link
      - '--level=2' tells the terminal to download just the files immediately in the directory of the link and the folders that they are in
      - if we had used '--level=1' then the terminal would set up all the folders as found in the link until the place where to files are
      - if we had used '--level=3' then the terminal would download all of the files in the immediate current directory but also in the directories above
  4. To move the files to another folder, try >> mv oldfolder/file newfolder/file
  5. To delete old folder >> rmdir oldfolder OR if the old folder has a number of different folders inside it but you want to get rid of all of these, then 
  do >> rm -R oldfolder
  
Use wget to download a number of different files:
  1. Create a file on Vim with the location for your files. Eg. a file that says:
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR287/005/SRR2876085/SRR2876085_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR287/005/SRR2876085/SRR2876085_2.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/007/SRR1248497/SRR1248497_1.fastq.gz
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/007/SRR1248497/SRR1248497_2.fastq.gz
      - to create the file:
          - *vim* ->  Will open Vim. Alternatively, when you open vim, you can do so with the file name, eg. *vim filename.txt*
          - *i* -> Will allow you to start typing
          - *esc + :w filename.txt* -> Will save the file with the name you suggested
          - *:q* -> Will close Vim
          - Alternatively, when you want to close and save the file (without naming the file), do *esc + :wq*
  2. Check that the file has what you want -> *less filename.txt*  Exit less with -> *q*
  3. Go to the directory where you want to download the new files to, then 
  
  >*cat filePaths.txt | xargs -i wget {}*
  
      - 'cat' will print out what is writen inside the file
      - 'xargs -i' takes each line of the file as an input for 'wget', everytime replacing '{}' with a new input line.
      
  
