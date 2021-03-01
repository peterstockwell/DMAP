/* rbs_fns.c: set of routines common to a number
of the bisulphite sequence programs: made into a 
separate file to simplify maintenance

Peter Stockwell 26-Sep-2010 */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "mrg_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"

/* global storage: to avoid making routines 'static' */
char glbl_chrstring[33];

char *any_chrno2str(int maxchrno,
                    int hasxy,
                    int chrno,
                    int terse)
/* depending on rdp, return a string for chrno */
{
if (chrno > maxchrno)
  snprintf(&glbl_chrstring[0],32,"%s?",(terse?"":"Chr"));
else
  if ((chrno <= maxchrno-2) || !hasxy)
    snprintf(&glbl_chrstring[0],32,"%s%d",(terse?"":"Chr"),chrno);
  else
    snprintf(&glbl_chrstring[0],32,"%s%s",(terse?"":"Chr"),
               (chrno==maxchrno?"Y":"X"));
return(&glbl_chrstring[0]);
}

char *rbc_chrno2str(RBC_CHRNO cno,
                    int terse)
{
return(any_chrno2str(ChrY,1,cno,terse));
}

int any_str2chrno(int maxchrno,
                  int hasxy,
                  char *ustr)
/* see if ustr matches any thing on chromosome
list with maxchrno/hasxy constraints.
Return that value, else 0 =(CHR_unk) */
{
int cp;

cp = 1;
while (cp <= maxchrno)
  if (strcmp(any_chrno2str(maxchrno,hasxy,cp,1),ustr) == 0)
    return(cp);
  else
    cp++;
return(0);
}

RBC_CHRNO rbc_str2chrno(char *cstr)
  /* return the most appropriate chromosome No for cstr */
{
return((RBC_CHRNO) any_str2chrno(ChrY,1,cstr));
}

int rbc_remaptileno(RBC_FLOWCELLVERSN fcv,
                    int fastqtileno)
/* remap flowcell V3 tiles to 1..24.  V2 just returns
the value */
{
int surface;  /* 1,2 */
int swath;  /* 1,2,3 */
int tdiv;

switch (fcv)
  {
  case RBC_fcv_3_2:
    surface = (int) fastqtileno/1000;
    fastqtileno %= 1000;
    swath = (int) fastqtileno/100;
    tdiv = fastqtileno%100;
    return(16*(surface-1)+8*(swath-1)+tdiv);
    break;
  case RBC_fcv_3_3:
    surface = (int) fastqtileno/1000;
    fastqtileno %= 1000;
    swath = (int) fastqtileno/100;
    tdiv = fastqtileno%100;
    return(24*(surface-1)+8*(swath-1)+tdiv);
    break;
  case RBC_fcv_2:
  default:
    return(fastqtileno);
    break;
  }
}

int rbc_intinrng(int b1,
                 int i,
                 int b2)
/* 1 if b1 <= i <= b2 or b2 <= i <= b1 */
{
if (b1 <= b2)
  return((b1 <= i) && (i <= b2));
else
  return(rbc_intinrng(b2,i,b1));
}

int rbc_modnwithn(int base,
                  int i)
/* return 1,2,...n, on a base of n */
{
int rtv;

if (base <= 0)
  return(0);
else
  if ((rtv = (i % base)) == 0)
    return(base);
  else
    return(rtv);
}

int rbc_invrtremaptileno(RBC_FLOWCELLVERSN fcv,
                         int tno)
/* tno is a tile number in range 1..120.  reverse the
remapping back to the original number.  Check that
tno is in appropriate range and return 0 if not */
{
int surface;
int swath;
int wosurf;

switch (fcv)
  {
  case RBC_fcv_3_2:
    if (rbc_intinrng(1,tno,32))
      {
      surface = (int) (tno - 1)/(8*2);
      wosurf = tno % 16;
      swath = (int)(wosurf + 7)/8;
      return(surface*1000 + swath*100 + rbc_modnwithn(8,tno));
      }
    else
      return(0);
    break;
  case RBC_fcv_3_3:
    if (rbc_intinrng(1,tno,48))
      {
      surface = (int) (tno - 1)/(8*3);
      wosurf = tno % 24;
      swath = (int)(wosurf + 7)/8;
      return(surface*1000 + swath*100 + rbc_modnwithn(8,tno));
      }
    else
      return(0);
    break;
  case RBC_fcv_2:
  default:
    if (rbc_intinrng(1,tno,MAXTILENO))
      return(tno);
    else
      return(0);
    break;
  }
}

int rbc_cntflds(char *buf,
                char fldsep)
/* scan buff for fldsep, return No of distinct fields
delimited by fldsep.  Note that leading or terminal 
occurrences are not considered */
{
char *bp;
int fcnt;
char prv;

bp = buf;
prv = '\0';
if (*bp == fldsep)
  fcnt = 0;
else
  fcnt = 1;
while ((prv = *bp) != '\0')
  {
  if (*bp == fldsep)
    fcnt++;
  bp++;
  }
if (prv == fldsep)
  fcnt--;
return(fcnt);
}

char *rbc_skiptofld(char *buf,
                    char fldsep,
                    int fldno)
/* jump to field no fldno in buf, based on
occurrences of fldsep.  return the point immediately
following the preceding delimiter, so if two successive
delimiters exist, the next will be found.
return NULL if the field doesn't exist */
{
int fcnt;
char *bp;

bp = buf;
if (*bp == fldsep)
  bp++;
fcnt = 1;
while (*bp != '\0')
  if (fcnt == fldno)
    return(bp);
  else
    {
    if (*bp == fldsep)
      fcnt++;
    bp++;
    }
return(NULL);
}

int rbc_readlntobuf(FILE *sfl,
                    char *buf,
                    int buflen)
/* read chars from sfl up to end of line, if not '\0' or
EOF, put into buf if it is not null up to buflen,
return no chars read */
{
int nc;
char *dp;
int rlen;

rlen = 0;
if ((dp = buf) != NULL)
  *dp = '\0';
while (((nc = fgetc(sfl)) != (int) '\n') && (nc != EOF))
  {
  if ((buf != NULL) && (rlen < buflen))
    {
    *dp = (char) nc;
    dp++;
    *dp = '\0';
    }
  rlen++;
  }
return(rlen);
}

