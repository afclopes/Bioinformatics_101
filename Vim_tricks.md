# Working with vim

## Opening a file on commandline:
   - > vim
   - > vim filename
    
## Saving a file:
   - :w saves a file
   - :w filename : saves a file with a new name
   - :wq saves a file and quits vim
    
## Entering and exiting "insert" mode:
   - i: enters "insert" mode, allowing you to make modifications to the text
   - esc: exits "insert" mode
    
## Manuvering in text:
   - arrows + i: to move and enter "insert" mode
   - 0: to move to beginning of line
   - $: move to end of line
   - A: move to end of line and get there already in "insert" mode
   - b: move to beginning of word
   - e: move to end of word
   - gg: first line
   - G: last line

## Extras:
   - u: undoing one step at a time
   - dd: deletes line that cursor is on
   - shift + c: deletes from the cursor until the end of the line, and puts you in insert mode
   - shift + d: deletes from the cursor until the end of the line
   - yy: copy line where cursor is on
   - p: paste line after cursor
   - P: paste line before cursor
   - <esc> :noh<return><esc>: stops search with highlighted word
   - :%s/foo/bah/gc: search in entire text for the string 'foo', and replace with 'bah' only after asking me everytime to confirm the change
   - :%s/phrase//gc: search in the entire text for the string 'phrase', and delete it but only after asking me everytime to confirm the change (which within the 
   options I can say y or a for all)
   
## Two files side by side:

   1. Open both files at the same time:
   
> vim -O +'set nu' 'windo set scrollbind' SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_Q30_MitoReads_Split_new_QC_Sbatch.vcf SRR1411188_pe_deduplicated_merged_sorted_MAPQ20_bcftools_Q30_MitoReads_Split_new_QC.vcf

Where *+'set nu'* sets the line numbers in the files on vim. You can also do this once inside vim, by going to the command mode and doing *:set number*. And 
*-O* opens the files side-by-side. The command *'windo'* makes sure everything is done on both windows. 

   2. Open the files one at a time in vim:
   
> vim filename1.vcf
> :e /path/to/filename2.vcf --> do this once inside vim
> : windo set scrollbind

   3. Comparing two files in the entire files:
   
> vimdiff filename1 filename2

   4. Comparing two files at specific columns:
   
> vimdiff <(awk '{print $8","$9","$10","$11}' filename1) <(awk '{print $8","$9","$10","$11}' path/filename2)
