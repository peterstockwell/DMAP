# illum2fasta.awk: read illumina text files and convert header & seq lines to fasta format.
#  shorten headers.
#
# to run: awk -f illum2fasta.awk s1_trim75.txt > s1_trim75.fa
#
BEGIN { if (hdrstr=="") hdrstr = "s1";}

FNR%4==1 {ns=split($1,s1,":");
  printf(">%s_%s_%s_%s\n",hdrstr,s1[ns-2],s1[ns-1],s1[ns]);
  }
FNR%4==2 {
  printf("%s\n",$0);
  }
