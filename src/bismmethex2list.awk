# bismmethex2list.awk - script to convert output from
# BisMark methylation_extractor to simple list in form
# Chr position +/-
# where +=> methylated C & -=> unmeth.
#
# Peter Stockwell: 3-Feb-2011
#
# revised Nov-2011 for CASAVA 1.8 output
#
# zcol can be defined to fit the output
# awk -f <path>bismmethex2list.awk zcol=6 <reads>_bismark.txt
#
# zcol defaults to 5
BEGIN {zcol=5;}

$zcol == "Z" { for (i=NF-2;i<NF;i++)printf("%s\t",$i); printf("+\n");}
$zcol == "z" { for (i=NF-2;i<NF;i++)printf("%s\t",$i); printf("-\n");}
