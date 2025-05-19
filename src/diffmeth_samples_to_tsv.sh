#!/bin/bash

# diffmeth_samples_to_tsv.sh: to take a file of diffmeth output generated
# with the -B (anova more details) option and decompose final sample
# column into tab separated columns.
#
# Peter Stockwell May-2025
#
# run with no parameters for help info
#

if [[ -z $1 ]]; then

echo "diffmeth_samples_to_tsv.sh: to take a file of diffmeth output generated";
echo "with the -B (anova more details) option and decompose final sample";
echo "column into tab separated columns.";
echo "";
echo "Input parameter is a diffmeth output file for -B anova option";
echo "";
echo "Optional second parameter is the output file, defaults to <stdout>";
echo "";
echo "This script will generate two awk scripts to complete the work:";
echo "anovasmpls2stdout.awk & dmeth_to_cols.awk.  These can be deleted";
echo "after completion.";

exit 0;

fi

# Check we can read this file

if [[ ! -f $1 ]]; then

echo "Can't read file '$1'"

exit 1;

fi

# OK, so generate awk scripts

cat << 'SCRIPT1' > anovasmplstostdout.awk
# find_anova_samples.awk: script to scan diffmeth output for -B option
# (ANOVA with more detail) and find the number of samples for R & S
# groups
#
# write to a named file (defaults to ANOVA_group_list.txt
#
# Peter Stockwell May-2025

NR>1{nflds = split($NF,snf,";");
for (i = 2; i<= nflds; i++)
  decomposestr(snf[i]);
}

END{for (sn in smpllist)
  printf("%s\n",sn) | "sort";
}

function decomposestr(smplset)
# smpset is a series of comma delimited items like 'Ra=1.00'
# pull them apart and add to list
{
nsst = split(smplset,spltsmplset,",");
for (j = 1; j<= nsst; j++)
  addtocounts(spltsmplset[j]);
}

function addtocounts(smplstring)
# smplestring is a single item like 'Ra=1.00'
# add this to the smpllist array
{
split(smplstring,ssst,"=");
smpllist[ssst[1]]++;
}

SCRIPT1

cat << 'SCRIPT2' > dmeth_to_cols.awk
# dmeth_to_cols.awk: convert diffmeth output to tab delimited text,
# especially the last column of methylation proportions
#
# Version 1.1: to manage multi-group ANOVA runs
# Peter Stockwell May-2025

BEGIN{FS = "\t";
smplno = 0;
}

NR==1{printf("%s",$0);
cmd = sprintf("awk -f anovasmplstostdout.awk %s",FILENAME);
while (cmd | getline smplid > 0)
  smpllist[++smplno] = smplid;
close(cmd);

# look for number of groups: by scanning for R & S or just R in first char position,
# and further for 3 letter sample IDs (e.g. R1a,R2a, etc.) which indicate > 2 groups.
maxidlen = 0;
for (i = 1; i<= smplno; i++)
  {
  chararray[nthchar(smpllist[i],1)]++;
  if ((idlen = length(smpllist[i])) > maxidlen)
    maxidlen = idlen;
  }
nogroups = arraylength(chararray);
if ((nogroups < 2) || (maxidlen > 2))
  {
  for (aelmt in chararray)
    delete chararray[aelmt];
  for (i = 1; i <= smplno; i++)
    {
    chararray[nthchar(smpllist[i],2)]++;
    }
  nogroups = arraylength(chararray);
  }
# generate maxidlen-1 letter array of IDs
gno = 0;
prv2id = "";
for (i = 1; i <= smplno; i++)
  if ((new2id = substr(smpllist[i],1,maxidlen-1)) != prv2id)
    {
    id2letter[++gno] = new2id;
    prv2id = new2id;
    grpprop[gno] = "";
    }
for (i = 2; i <= nogroups; i++)
  printf("\t+");
for (i = 1; i <= smplno; i++)
  printf("\t%s",smpllist[i]);
printf("\n");
}

NR>1{for (i=1; i<NF; i++)
       printf("%s%s",$i,(i<(NF-1)?"\t":""));
nsnf = split($NF,snf,";")
# manage first part of overall group totals, 1st clear smplstr
for (i = 1; i <= nogroups; i++)
  grpprop[i] = "";
commatogrplist(snf[1]);
for (i = 1; i<= nogroups; i++)
  {
  if (substr(grpprop[i],1,maxidlen-1) == id2letter[i])
    printf("\t%s",grpprop[i]);
  else
    printf("\t-");
  }
for (i = 1; i <= smplno; i++)
  smplstr[smpllist[i]] = "";
for (k = 2; k <= nsnf; k++)
  commatosmpllist(snf[k]);
for (i = 1; i <= smplno; i++)
  if (smplstr[smpllist[i]] != "")
    printf("\t%s",smplstr[smpllist[i]]);
  else
    printf("\t-");
printf("\n");
}

function commatogrplist(cmmasepstr)
# break cmmasepstr is a set of comma delimited items like 'Ra=1.00'
# pull apart and store values in smplstr array
{
ncs = split(cmmasepstr,css,",");
for (j = 1; j<= ncs; j++)
  {
  split(css[j],scss,"=");
  for (gc = 1; gc <= nogroups; gc++)
    if (id2letter[gc] == scss[1])
      grpprop[gc] = css[j];
  }
}

function commatosmpllist(cmmasepstr)
# break cmmasepstr is a set of comma delimited items like 'Ra=1.00'
# pull apart and store values in smplstr array
{
ncs = split(cmmasepstr,css,",");
for (j = 1; j<= ncs; j++)
  {
  split(css[j],scss,"=");
  for (i = 1; i <= smplno; i++)
    if (smpllist[i] == scss[1])
      smplstr[smpllist[i]] = css[j];
  }
}

function comma2tabout(cstring)
# write out comma delimited fields of cstring tab delimited
# ensure complete list is output, dashes for absent 2 letter IDs
{
ocnt = 0;
ns = split(cstring,css,",");
for (j = 1; j<= ns; j++)
  if (ocnt++ == pos2charid(css[j]))
    printf("%s%s",css[j],(j<ns?"\t":""));
  else
    printf("-%s",(j<ns?"\t":""));
}

function arraylength(anarray)
# count number of elements in anarray
{
cnt = 0;
for (elmnt in anarray)
  cnt++;
return cnt;
}

function nthchar(astring,nchar)
# return an array of single chars in astring
{
if (nchar <= length(astring))
  return(substr(astring,nchar,1));
else
  return("\0");
}

function pos2charid(qid)
# return the position of the leading 2 chars of qid in id2letter array
# 0 if not
{
for (gno = 1; gno <= nogroups; gno++)
  if (id2letter[gno] == substr(qid,1,2))
    return(gno);
return 0;
}
SCRIPT2

# now call dmeth_to_cols.awk on input file

if [[ -z $2 ]]; then

awk -f dmeth_to_cols.awk $1

else

awk -f dmeth_to_cols.awk $1 > $2

fi

# done

exit 0;
