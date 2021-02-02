# How to install:

1. Go to website http://software.broadinstitute.org/software/igv/ and download the one that matches your operational system.
2. Drag and drop the installed icon onto your desktop

# How to upload sequencing file:
1. Transfer your sequencing files .bam and .bam.bai into a same folder.
2. File > Load File > .bam
3. Select reference genome
4. Select chromosome of interest

# How to move aorund IGV:
  1. double click to zoom in, Alt + click to zoom out.
  
  
# What do the colours mean?
  1. grey refers to the read aligned to the reference genome
  2. coloured letters refer to the bases that did not match the reference. You can opt to colour all of the bases instead of just the ones that don't match.
  3. bolder letters refer to higher quality, and fainter letters refer to lower quality suggesting that they are less likely to be a real result.
      - if you want to remove this, click on "shade base by quality"
  4. show 'reads shading' by:
      - hidding the mismatch bases
      -  View > Preferences > Unclick 'Show center line'
      - Notice different colours of reads:
          + hallow reads: mapping quality of zero
          + different colour of reads: different mapping quality
          + purple I: represent insertions. If you zoom in on it, it will state how many bases may have been inserted, eg. 2. Clicking on the number will show which bases were inserted.
          + black dash: represent deletions. when multiple bases are deleted, the number of bases deleted is shown as a number in the middle of the dash.
 5. colour reads according to strand: Colour alignment by colour > Read strand:
      - red is forward
      - blue is reverse
 6. look for false positive: looking for strand bias
      - Sort alignments by colour > Read strand:
      - a true heterozygous genotype would show even distribution of the mismatch base across forward and reverse strands
 7. Read strand in pastels, red for positive rightward (5' to 3') DNA strand, blue for negative leftward (reverse-complement) DNA strand, and grey for unpaired mate, mate not mapped, or otherwise unknown status.
 8. The sequence is represented by colored bars or colored letters, depending on zoom level, with adenine in green, cytosine in blue, guanine in yellow, and thymine in red
 9. At the very bottom: Amino acids are displayed as blocks colored in alternating shades of gray. Methionines are colored green, and all stop codons are colored red. When you zoom all the way in, the amino acid symbols will appear. 
 10. Note that alignments that are displayed with light gray borders and transparent or white fill, have a mapping quality equal to zero
