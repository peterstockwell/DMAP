# designed to be run in association with diffmeth
# in order to pull apart the detailed CpG info for
# a single sample and express the counts as
# %methylation for any CpG with non-zero counts.
#
# diffmeth should be run with option -e 40,220
# on one sample and the output piped through this
# script.
#
# Peter Stockwell: Jul-2013
BEGIN{chr="";}
$1!="Sample"{chr=$1;}
$1=="Sample" && $2=="1" {for (i=3;i<=NF;i++){split($i,si,":");split(si[2],si2,"+"); \
  meth=si2[1];unmeth=substr(si2[2],1,length(si2[2])-1); \
  if (meth+unmeth > 0){pc=100*meth/(meth+unmeth); printf("%s\t%s\t%.2f\n",chr,si[1],pc);} \
  else printf("%s\t%s\t-\n",chr,si[1]);}}
