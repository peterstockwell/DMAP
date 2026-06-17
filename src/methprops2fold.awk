# methprops2fold.awk: to take diffmeth output with the final column
# headed by R=rprop,S=sprop (e.g. R=0.1250,S=0.2093;R#1=0.1000,R#5=0.1429;S#3=0.0833,S#5=0.2581)
# and return fold difference.
#
# Note: 0.0 meth values for R & S samples will return '-'
#
# Usage:
#
# awk -f methprops2fold.awk <diffmeth_output_file>
#
# or pipe diffmeth output directly to 'awk -f methprops2fold.awk'
#
# Peter Stockwell: 18-Jun-2026 (revised from 2022 version)
#

$1~/#/{printf("%s\tFoldDiff\n",$0);}

$1!~/#/{split($NF,snf,";");
  fdiff = split_fold_diff(snf[1]);
  if (fdiff < 0.0)
    printf("%s\t-\n",$0);
  else
    printf("%s\t%.4f\n",$0,fdiff);
  }

function fold_diff(x,y)
# return fold difference between x & y (real values), -99.0 for
# invalid.
{
if ((x == 0.0) || (y == 0.0))
  return(-99.0);
else
  if (x > y)
    return(x/y);
  else
    return(y/x);
}

function split_fold_diff(propstring)
# split propstring into R & S real values and
# return fold difference
{
ns = split(propstring,pstr,/[=,]/);
return(fold_diff(pstr[2]+0.0,pstr[ns]+0.0));
}
