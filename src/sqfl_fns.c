/* sqfl_fns.c: c functions for sequence file accessing */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <libgen.h>
#include <sys/param.h>
#include <sys/stat.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"

SFMT_TYPE sqfl_chr2fmttp(char fchr)
  /* convert the letters s,q,f,n,u,r to corresponding file format type.
  currently case independent */
{
switch (tolower(fchr))
  {
  case 's':
    return(SFMT_staden);
    break;
  case 'q':
    return(SFMT_molgen);
    break;
  case 'f':
    return(SFMT_fasta);
    break;
  case 'n':
    return(SFMT_nbrf);
    break;
  case 'u':
    return(SFMT_gcg);
    break;
  case 'r':
    return(SFMT_raw);
    break;
  default:
    return(SFMT_undefined);
    break;
  }
}

SFMT_TYPE sqfl_getdeffmt()
  /* use env variable SQDEFFILEFMT to establish default file type */
{
char *enm;

if ((enm = getenv("SQDEFFILEFMT")) != NULL)
  return(sqfl_chr2fmttp(*enm));
else
  return(SFMT_undefined);
}

char *sqfl_fmttp2strng(SFMT_TYPE sfmt)
  /* return the pointer to a string defining sfmt */
{
switch (sfmt)
  {
  case SFMT_staden:
    return("Staden");
    break;
  case SFMT_molgen:
    return("Molgen/SEQ");
    break;
  case SFMT_fasta:
    return("FASTA");
    break;
  case SFMT_nbrf:
    return("NBRF");
    break;
  case SFMT_gcg:
    return("GCG");
    break;
  case SFMT_raw:
    return("RAW");
    break;
  case SFMT_undefined:
  default:
    return("Undefined");
    break;
  }
}

char sqfl_fmttp2chr(SFMT_TYPE sfmt)
  /* return a character sfmt */
{
switch (sfmt)
  {
  case SFMT_staden:
    return('S');
    break;
  case SFMT_molgen:
    return('Q');
    break;
  case SFMT_fasta:
    return('F');
    break;
  case SFMT_nbrf:
    return('N');
    break;
  case SFMT_gcg:
    return('U');
    break;
  case SFMT_raw:
    return('R');
    break;
  case SFMT_undefined:
  default:
    return('?');
    break;
  }
}

char sqfl_restype2chr(SQ_RESTYPE rt)
  /* return the character corresponding to rt */
{
switch (rt)
  {
  case RES_blk:
    return(' ');
    break;
  case RES_a:
    return('A');
    break;
  case RES_b:
    return('B');
    break;
  case RES_c:
    return('C');
    break;
  case RES_g:
    return('G');
    break;
  case RES_t:
    return('T');
    break;
  case RES_d:     /* Asp */
    return('D');
    break;
  case RES_f:
    return('F');
    break;
  case RES_h:
    return('H');
    break;
  case RES_i:
    return('I');
    break;
  case RES_k:
    return('K');
    break;
  case RES_l:
    return('L');
    break;
  case RES_m:
    return('M');
    break;
  case RES_n:
    return('N');
    break;
  case RES_p:
    return('P');
    break;
  case RES_q:
    return('Q');
    break;
  case RES_e:
    return('E');
    break;
  case RES_r:
    return('R');
    break;
  case RES_s:
    return('S');
    break;
  case RES_v:
    return('V');
    break;
  case RES_w:
    return('W');
    break;
  case RES_y:
    return('Y');
    break;
  case RES_z:
    return('Z');
    break;
  case RES_x:
    return('-');
    break;
  default:
    return('?');
    break;
  }
}

char *sqfl_restype2str(SQ_RESTYPE rt)
  /* return the address of a string corresponding to rt */
{
switch (rt)
  {
  case RES_blk:
    return("   ");
    break;
  case RES_a:
    return("Ala");
    break;
  case RES_b:
    return("Asx");
    break;
  case RES_c:
    return("Cys");
    break;
  case RES_g:
    return("Gly");
    break;
  case RES_t:
    return("Thr");
    break;
  case RES_d:     /* Asp */
    return("Asp");
    break;
  case RES_f:
    return("Phe");
    break;
  case RES_h:
    return("His");
    break;
  case RES_i:
    return("Ile");
    break;
  case RES_k:
    return("Lys");
    break;
  case RES_l:
    return("Leu");
    break;
  case RES_m:
    return("Met");
    break;
  case RES_n:
    return("Asn");
    break;
  case RES_p:
    return("Pro");
    break;
  case RES_q:
    return("Gln");
    break;
  case RES_e:
    return("Glu");
    break;
  case RES_r:
    return("Arg");
    break;
  case RES_s:
    return("Ser");
    break;
  case RES_v:
    return("Val");
    break;
  case RES_w:
    return("Trp");
    break;
  case RES_y:
    return("Tyr");
    break;
  case RES_z:
    return("Glx");
    break;
  case RES_x:
    return("---");
    break;
  default:
    return("???");
    break;
  }
}

SQ_RESTYPE sqfl_chr2restype(char res)
  /* return an internal restype for this character */
{
switch (toupper(res))
  {
  case 'A':
    return(RES_a);
    break;
  case 'B':
    return(RES_b);
    break;
  case 'C':
    return(RES_c);
    break;
  case 'G':
    return(RES_g);
    break;
  case 'U':
  case 'T':
    return(RES_t);
    break;
  case 'D':
    return(RES_d);
    break;
  case 'F':
    return(RES_f);
    break;
  case 'H':
    return(RES_h);
    break;
  case 'I':
    return(RES_i);
    break;
  case 'K':
    return(RES_k);
    break;
  case 'L':
    return(RES_l);
    break;
  case 'M':
    return(RES_m);
    break;
  case 'N':
    return(RES_n);
    break;
  case 'P':
    return(RES_p);
    break;
  case 'Q':
    return(RES_q);
    break;
  case 'E':
    return(RES_e);
    break;
  case 'R':
    return(RES_r);
    break;
  case 'S':
    return(RES_s);
    break;
  case 'V':
    return(RES_v);
    break;
  case 'W':
    return(RES_w);
    break;
  case 'Y':
    return(RES_y);
    break;
  case 'Z':
    return(RES_z);
    break;
  case ' ':
    return(RES_blk);
    break;
  default:
    return(RES_x);
    break;
  }
}

SQ_RESTYPE sqfl_chr2narestype(char rs)
  /* return RES_x..RES_t only */
{
SQ_RESTYPE rt;

if ((rt = sqfl_chr2restype(rs)) <= RES_t)
  return(rt);
else
  return(RES_x);
}

SQ_RESTYPE sqfl_chr2aarestype(char r)
  /* return the corresponding restype for r, assuming it is a valid AA code,
else return RES_x */
{
if (toupper(r) != 'U')
  return(sqfl_chr2restype(r));
else
  return(RES_x);
}

int sqfl_linelength(SFMT_TYPE sfmt)
  /* return normal linelength used by this format */
{
switch (sfmt)
  {
  case SFMT_molgen:
    return(70);
    break;
  case SFMT_gcg:
    return(50);
    break;
  case SFMT_raw:
    return(0);
    break;
  case SFMT_staden:
  case SFMT_fasta:
  case SFMT_nbrf:
  case SFMT_undefined:
  default:
    return(60);
    break;
  }
}

int sqfl_gcgchkinc(int sp,
                   char res)
/* return the incremental addition for the GCG check algorithm for res at
  position sp, (1..n) */
{
return((1 + --sp % SQFL_GCGCHK)*((int) toupper(res)));
}

int sqfl_gcgsqchksum(char *sq,
                     int sqlen)
/* return gcg checksum for (1..sqlen) chars of sq.  Stop at '\0'
if encountered ealier */
{
int csum;
char *sp;
int pt;

csum = 0;
sp = sq;
pt = 1;
while ((pt <= sqlen) && (*sp != '\0'))
  csum += sqfl_gcgchkinc(pt++,*sp++);
return(csum % 10000);
}

SQFL_SQTOS sqfl_scantos(char *sqbuf,
                        int sqln)
/* scan sqbuf, looking for evidence of NA or peptide.
  Algorithm looks for the letters E,F,I,O,X,Z: none of which should occur in
  DNA sequences, even using IUB or LDNA redundant base codes.  A proportion
  1% of strangers is allowed */
{
int ccnt[27];
char *sp;
int pt;
int dtot;
char uc;
int ptot;

sp = sqbuf;
for (pt = 0; pt < 27; pt++)
  ccnt[pt] = 0;
while ((*sp != '\0') && (sqln-- > 0))
  if (((uc = toupper(*sp++)) >= 'A') && (uc <= 'Z'))
    ccnt[((int) uc - 'A' + 1)]++;
  else
    ccnt[0]++;
dtot = ptot = 0;
for (pt = 1; pt < 27; pt++)  /* ignore unknowns */
  if (index("EFIOXZ",((char) pt + 'A' - 1)) != NULL) /* peptide residue */
    ptot += ccnt[pt];
  else
    dtot += ccnt[pt];
if (ptot > dtot/100)
  return(SQFL_peptide);
else
  if (ccnt[(int) ('T' - 'A' + 1)] > ccnt[(int) ('U' - 'A' + 1)])
    return(SQFL_dna);
  else
    return(SQFL_rna);
}

SQFL_SQTOS sqfl_scanfltos(SQFL_STRCT *src)
  /* scan contents of previously opened src, looking for evidence of 
NA or peptide.  Rewind source file on completion.
  Algorithm looks for the letters E,F,I,O,X,Z: none of which should occur in
  DNA sequences, even using IUB or LDNA redundant base codes.  A proportion
  1% of strangers is allowed */
{
int ccnt[27];
char nr;
int pt;
int dtot;
char uc;
int ptot;

if (sqfl_skipsqflhdr(src))
  {
  for (pt = 0; pt < 27; pt++)
    ccnt[pt] = 0;
  while ((nr = sqfl_getnxtres(src)) != '\0')
    if (((uc = toupper(nr)) >= 'A') && (uc <= 'Z'))
      ccnt[((int) uc - 'A' + 1)]++;
    else
      ccnt[0]++;
  sqfl_rewind(src);
  dtot = ptot = 0;
  for (pt = 1; pt < 27; pt++)  /* ignore unknowns */
    if (index("EFIOXZ",((char) pt + 'A' - 1)) != NULL) /* peptide residue */
      ptot += ccnt[pt];
    else
      dtot += ccnt[pt];
  if (ptot > dtot/100)    /* allow 1% "foreigners" */
    return(SQFL_peptide);
  else
    if (ccnt[(int) ('T' - 'A' + 1)] > ccnt[(int) ('U' - 'A' + 1)])
      return(SQFL_dna);
    else
      return(SQFL_rna);
  }
else
  return(SQFL_tosunknown);
}

char sqfl_tos2chr(SFMT_TYPE sfmt,
                  SQFL_SQTOS st)
  /* return the character for sequence st - default to D */
{
switch (sfmt)
  {
  case SFMT_nbrf:
  case SFMT_fasta:
    switch (st)
      {
      case SQFL_peptide:
        return('P');
        break;
      case SQFL_rna:
        return('R');
        break;
      case SQFL_dna:
      case SQFL_tosunknown:
      default:
        return('D');
        break;
      }
    break;
  case SFMT_gcg:
    switch (st)
      {
      case SQFL_peptide:
        return('P');
        break;
      case SQFL_rna:
      case SQFL_dna:
      case SQFL_tosunknown:
      default:
        return('N');
        break;
      }
    break;
  case SFMT_staden:
  case SFMT_molgen:
  case SFMT_undefined:
  case SFMT_raw:
  default:
    return('\0');
    break;
  }
return('D');
}

void sqfl_settos(SQFL_STRCT *st,
                 SQFL_SQTOS tos)
/* set tos field of st */
{
st->sqtos = tos;
}

void sqfl_settopol(SQFL_STRCT *st,
                   SQFL_SQTOPOL tpl)
/* set topology field of st to tpl */
{
st->stopol = tpl;
}

void sqfl_setsqdetails(SQFL_STRCT *st,
                       char *sqbuf,
                       int sqlen)
/* set as many as possible of st details by examination of sqbuf */
{
if ((sqbuf != NULL) && (sqlen > 0))
  {
  if (st->sqtos == SQFL_tosunknown)
    st->sqtos = sqfl_scantos(sqbuf,sqlen);
  if (st->stopol == SQTP_unknown)
    switch (st->sqtos)
      {
      case SQFL_peptide:
        st->stopol = SQTP_linear;
        break;
      case SQFL_dna:
      case SQFL_rna:
      case SQFL_tosunknown:
      default:
        break;
      }
  }
else
  st->sqtos = SQFL_tosunknown;
st->slength = sqlen;
/* if (st->flfmt == SFMT_gcg)
  st->gcgchksum = sqfl_gcgsqchksum(sqbuf,sqlen); */
}

void sqfl_headsfstr(SQFL_STRCT *sst,
                    char *sqname,
                    char *ann)
/* write sequence header info to previously opened sst->sfl, using string ann
 */
{
char *qp;

switch (sst->flfmt)
  {
  case SFMT_fasta:
    if (sqname != NULL)
      fprintf(sst->sfl,">%s - ",sqname);
    fprintf(sst->sfl,"%s\n",((ann != NULL)?ann:""));
    break;
  case SFMT_nbrf:
    if (sqname != NULL)
      fprintf(sst->sfl,">%c%c;%s\n",sqfl_tos2chr(sst->flfmt,sst->sqtos),
                sqfl_topol2chr(sst->flfmt,sst->stopol),sqname);
    else
      fprintf(sst->sfl,">DL;%s\n",((sst->filnam != NULL)?sst->filnam:"?SEQ"));
    if (ann != NULL)
      fprintf(sst->sfl,"%s\n",ann);
    else
      fputc('\n',sst->sfl);
    break;
  case SFMT_molgen:
    if (sqname != NULL)
      fprintf(sst->sfl,"; %s: ",sqname);
    qp = ann;
    while ((qp != NULL) && (*qp != '\0'))
      {
      if (*qp == '\n')
        fprintf(sst->sfl,"\n; ");
      else
        fputc(*qp,sst->sfl);
      qp++;
      }
    fputc('\n',sst->sfl);
    if (sqname != NULL)
      fprintf(sst->sfl,"%s\n",sqname);
    break;
  case SFMT_gcg:
  case SFMT_staden:
  case SFMT_undefined:
  case SFMT_raw:
  default:
    break;
  }
}

void sqfl_headsfstrct(SQFL_STRCT *sst,
                      char *sqname,
                      char *maker,
                      char *origin)
/* write sequence header info to previously opened sst->sfl */
{
char hbuf[HDRBUF_MAX+1];
char *bp;

if (sqname != NULL)
  sst->seqnam = bas_strdup(basename(sqname));
else
  sst->seqnam = NULL;
switch (sst->flfmt)
  {
  case SFMT_nbrf:
  case SFMT_fasta:
  case SFMT_molgen:
    bp = &hbuf[0];
    if (maker != NULL)
      {
      bas_appstr(&hbuf[0],&bp," Created by ",HDRBUF_MAX);
      bas_appstr(&hbuf[0],&bp,maker,HDRBUF_MAX);
      if (origin != NULL)
        {
        bas_appstr(&hbuf[0],&bp," from ",HDRBUF_MAX);
        bas_appstr(&hbuf[0],&bp,origin,HDRBUF_MAX);
        }
      bas_appstr(&hbuf[0],&bp," on ",HDRBUF_MAX);
      bas_appdate2str(&hbuf[0],&bp,HDRBUF_MAX);
      }
    if ((origin != NULL) && (maker == NULL))
      bas_appstr(&hbuf[0],&bp,origin,HDRBUF_MAX);
    if (sst->annot != NULL)
      {
      bas_appchr(&hbuf[0],&bp,' ',HDRBUF_MAX);
      bas_appstr(&hbuf[0],&bp,sst->annot,HDRBUF_MAX);
      }
    sqfl_headsfstr(sst,sqname,&hbuf[0]);
    break;
  case SFMT_gcg:
  case SFMT_staden:
  case SFMT_undefined:
  case SFMT_raw:
  default:
    break;
  }
}

void sqfl_headsfstrctann(SQFL_STRCT *sst,
                         char *sqname)
/* write sequence header info to previously opened sst->sfl. use annotation
  if set */
{
char hbuf[HDRBUF_MAX+1];
char *bp;

if (sst->annot != NULL)
  sqfl_headsfstr(sst,sqname,sst->annot);
else
  {
  bp = &hbuf[0];
  bas_appstr(&hbuf[0],&bp,"Written by sqfl_fns on ",HDRBUF_MAX);
  bas_appdate2str(&hbuf[0],&bp,HDRBUF_MAX);
  sqfl_headsfstr(sst,sqname,&hbuf[0]);
  }
}

void sqfl_newlnstart(FILE *sfl,
                     SFMT_TYPE sfmt,
                     int sp)
/* if necessary, write this number as a line start */
{
switch (sfmt)
  {
  case SFMT_gcg:
      if ((sfmt == SFMT_gcg) && (sp > 0))
    fprintf(sfl,"%8d  ",sp);
    break;
  case SFMT_undefined:
  case SFMT_staden:
  case SFMT_molgen:
  case SFMT_nbrf:
  case SFMT_fasta:
  case SFMT_raw:
  default:
    break;
  }
}

void sqfl_sqchr2fl(FILE *fl,
                   SFMT_TYPE sfmt,
                   char chr,
                   int llen,
                   int rp,
                   int *lcnt)
/* guts of writing sequence chars to a stream file */
{
if ((*lcnt >= llen) && (llen > 0))
  {
  fputc('\n',fl);
  *lcnt = 0;
  sqfl_newlnstart(fl,sfmt,rp);
  }
if (sfmt == SFMT_gcg)
  {
  if (rp == 1)
    sqfl_newlnstart(fl,SFMT_gcg,rp);
  if ((*lcnt % 10) == 0)
    fputc(' ',fl);
  }
fputc(chr,fl);
(*lcnt)++;
}

void sqfl_putchr(SQFL_STRCT *sqs,
                 char chr,
                 int rp,
                 int *lcnt)
/* put chr out to file, ignoring validity. observe line count etc */
{
sqs->scnt++;
if (sqs->flfmt == SFMT_gcg)
  {
  sqs->gcgchksum += sqfl_gcgchkinc(sqs->scnt,chr);
  fputc(chr,sqs->gfl);
  }
else
  sqfl_sqchr2fl(sqs->sfl,sqs->flfmt,chr,sqs->lnlen,rp,lcnt);
}

void sqfl_putres(SQFL_STRCT *sqs,
                 char res,
                 int rp,
                 int *lcnt)
/* put res out to file, keeping track of valid bases, line count etc */
{
if (sqs->validlut[(int) res])
  sqfl_putchr(sqs,res,rp,lcnt);
else
  if ((sqs->flfmt == SFMT_gcg) && (res == '-'))
    sqfl_putchr(sqs,'X',rp,lcnt);
}

int sqfl_fgetc(FILE *fl)
  /* return character (as int) from fl.  ignore '\r' chars totally */
{
int nc;

while ((nc = fgetc(fl)) == '\r');
return(nc);
}

void sqfl_termsqfl(SQFL_STRCT *sqs,
                   char *pname,
                   int *lcnt)
/* finish off sequence if necessary */
{
int sp;
int sc;

switch (sqs->flfmt)
  {
  case SFMT_staden:
    sqfl_putchr(sqs,'@',0,lcnt);
    break;
  case SFMT_molgen:
    sqfl_putchr(sqs,sqfl_topol2chr(sqs->flfmt,sqs->stopol),0,lcnt);
    break;
  case SFMT_nbrf:
    sqfl_putchr(sqs,'*',0,lcnt);
    break;
  case SFMT_gcg:
    if (sqs->annot != NULL)
      fprintf(sqs->sfl,"%s: %s",sqs->filnam,sqs->annot);
    else
      {
      if (pname != NULL)
        fprintf(sqs->sfl,"%s: Created by %s ",sqs->filnam,pname);
      else
        fprintf(sqs->sfl,"%s: Created ",sqs->filnam);
      fprintf (sqs->sfl,"on ");
      sqfl_date2file(sqs->sfl);
      }
    fprintf(sqs->sfl,"\n  Length: %d Type: %c Check: %d ..\n",
              sqs->scnt,
              sqfl_tos2chr(sqs->flfmt,sqs->sqtos),
              (sqs->gcgchksum % 10000));
    fflush(sqs->gfl);
    rewind(sqs->gfl);
    sp = 1;
    while ((sc = sqfl_fgetc(sqs->gfl)) != EOF)
      sqfl_sqchr2fl(sqs->sfl,SFMT_gcg,(char) sc,sqs->lnlen,sp++,lcnt);
    break;
  case SFMT_fasta:
  case SFMT_undefined:
  case SFMT_raw:
  default:
    break;
  }
if (*lcnt > 0)
  fputc('\n',sqs->sfl);
}

int sqfl_lookforchr(FILE *sfl,
                    char ec)
/* scan sfl, looking for ec.  stop and return 1.  Return 0 if EOF encountered
  before hand */
{
int nc;

while ((nc = sqfl_fgetc(sfl)) != EOF)
  if ((char) nc == ec)
    return(1);
return(0);
}

int sqfl_look4chrs(FILE *sfl,
                   char *ecs)
/* scan sfl, looking for anything in ecs.  stop and return that char.
  Return EOF if not encountered by EOF */
{
int nc;

while ((nc = sqfl_fgetc(sfl)) != EOF)
  if (index(ecs,(char) nc) != NULL)
    return(nc);
return(EOF);
}

SQFL_SQTOPOL sqfl_chr2topol(SFMT_TYPE sfmt,
                            char tc)
/* return the topology for char tc */
{
switch (sfmt)
  {
  case SFMT_fasta:
  case SFMT_nbrf:
    switch (toupper(tc))
      {
      case 'L':
        return(SQTP_linear);
        break;
      case 'C':
        return(SQTP_circular);
        break;
      default:
        return(SQTP_unknown);
        break;
      }
    break;
  case SFMT_molgen:
    switch (tc)
      {
      case '1':
        return(SQTP_linear);
        break;
      case '2':
        return(SQTP_circular);
        break;
      default:
        return(SQTP_unknown);
        break;
      }
    break;
  case SFMT_staden:
  case SFMT_gcg:
  case SFMT_undefined:
  case SFMT_raw:
  default:
    return(SQTP_unknown);
    break;
  }
return(SQTP_unknown);
}

char sqfl_topol2chr(SFMT_TYPE sfmt,
                    SQFL_SQTOPOL st)
/* return a character corresponding to st for sfmt */
{
switch (sfmt)
  {
  case SFMT_fasta:
  case SFMT_nbrf:
    switch (st)
      {
      case SQTP_circular:
        return('C');
        break;
      case SQTP_linear:
      case SQTP_unknown:
      default:
        return('L');
        break;
      }
    break;
  case SFMT_molgen:
    switch (st)
      {
      case SQTP_circular:
        return('2');
        break;
      case SQTP_linear:
      case SQTP_unknown:
      default:
        return('1');
        break;
      }
    break;
  case SFMT_gcg:
  case SFMT_staden:
  case SFMT_undefined:
  case SFMT_raw:
  default:
    return('\0');
    break;
  }
return('1');
}

SQFL_SQTOS sqfl_chr2tos(SFMT_TYPE sfmt,
                        char tc)
/* return the tos for tc - case independent */
{
switch (sfmt)
  {
  case SFMT_nbrf:
  case SFMT_fasta:
    switch (toupper(tc))
      {
      case 'P':
        return(SQFL_peptide);
        break;
      case 'D':
        return(SQFL_dna);
        break;
      case 'R':
        return(SQFL_rna);
        break;
      default:
        return(SQFL_tosunknown);
        break;
      }
    break;
  case SFMT_gcg:
    switch (toupper(tc))
      {
      case 'P':
        return(SQFL_peptide);
        break;
      case 'N':
        return(SQFL_dna);
        break;
      default:
        return(SQFL_tosunknown);
        break;
      }
    break;
  case SFMT_molgen:
  case SFMT_staden:
  case SFMT_raw:
  default:
    return(SQFL_tosunknown);
    break;
  }
}

void sqfl_skipeol(SQFL_STRCT *sqs)
  /* try to skip to end of current line - check if just seen eol already */
{
if ((sqs->lstinchr != '\n') && (sqs->lstinchr != '\r'))
  /* allow for PC-sourced files */
  bas_skipeol(sqs->sfl,NULL);
sqs->lstinchr = '\n';
}

void sqfl_skipeolnstor(SQFL_STRCT *sqs)
  /* try to skip to end of current line - check if just seen eol already.
    store skipped chars to annotation buffer if defined */
{
int nc;
char *bp;

bp = sqs->annot;
if ((sqs->lstinchr != '\n') && (sqs->lstinchr != '\r'))
  while (((nc = fgetc(sqs->sfl)) != EOF) && ((char) nc != '\n') &&
            ((char) nc != '\r'))   /* allow for PC-sourced files */
    if (sqs->annot != NULL)
      bas_appchr(sqs->annot,&bp,(char) nc,sqs->abuflen);
sqs->lstinchr = '\n';
}

void sqfl_forceeol(SQFL_STRCT *sqs)
  /* force skip to eol irrespective of what has been seen todate */
{
sqs->lstinchr = '\0';
sqfl_skipeol(sqs);
}

int sqfl_cpy2eol(SQFL_STRCT *sqs,
                 char *buf,
                 int bufmax,
                 int (*ufgetc)(FILE *fl))
/* copy to next line end to buf, up to bufmax, using ufgetc,
  return no of chars written */
{
int ccnt;
char *bp;
int nc;

ccnt = 0;
bp = buf;
while (((nc = (*ufgetc)(sqs->sfl)) != (char) EOF) && (nc != '\n'))
  if ((ccnt++ < bufmax) && (buf != NULL))
    *bp++ = (char) nc;
if (bp != NULL)
  *bp = '\0';
sqs->lstinchr = (char) nc;
return(ccnt);
}

int sqfl_skipsqflhdr(SQFL_STRCT *sqs)
  /* skip file header, getting to start of sequence.  Eventually this routine
should be able to extract data from the header.  Return 1 if nothing untoward
occured */
{
int nc;
char *np;
WRD_LUSTRCT gcglu;  /* for parsing gcg header words */
char gcgwrd[MAXPATHLEN + 1];    /* use as general atom */
WLU_CHRLUTBL toklu;  /* separators for gcg tokens */
int nl;
char *sp;

switch (sqs->flfmt)
  {
  case SFMT_fasta:    /* may or may not have nbrf-style headers.. */
    if ((sqs->lstinchr == '>') || (sqfl_lookforchr(sqs->sfl,'>')))
      {
      wlu_maktoklu(&toklu[0]," \t\n\r");
      if (wlu_gettokensep(sqs->sfl,gcgwrd,MAXPATHLEN,&toklu[0],&sqs->lstinchr) <= 0)
        return(0);
      else
        {
        if (index(gcgwrd,';') != NULL)     /* NBRF-style header */
          {                /* next char is sequence type */
          sp = &gcgwrd[0];
          sqs->sqtos = sqfl_chr2tos(sqs->flfmt,*sp++);
          if (*sp == '\0')
            return(0);
          else
            {
            sqs->stopol = sqfl_chr2topol(sqs->flfmt,*sp++);
            sp++;
            }
          }         
        else        /* not nbrf fmt, use entire string as name */
          {
          sqs->sqtos = SQFL_tosunknown;
          sqs->stopol = SQTP_unknown;
          sp = &gcgwrd[0];
          }
        }
/*      (void) strncpy(sqs->seqnam,basename(&gcgwrd[0]),(size_t) SQNAMEMAXLEN); */
      (void) strncpy(sqs->seqnam,&gcgwrd[0],(size_t) SQNAMEMAXLEN);
      if ((sqs->lstinchr != '\n') && (sqs->lstinchr != '\r'))
        (void) sqfl_cpy2eol(sqs,sqs->annot,sqs->abuflen,sqfl_fgetc);
      else
        sqs->annot = bas_strdup("");  /* must have malloc()ed this */
      np = sqs->annot;
      sqs->abuflen = strlen(np);
      sqs->annot = bas_strdup(np);
      memfree(np);      
/*      sqfl_skipeolnstor(sqs); */
      return(strlen(sqs->seqnam) > 0);
      }
    break;
  case SFMT_nbrf:
    if ((sqs->lstinchr == '>') || (sqfl_lookforchr(sqs->sfl,'>')))
      if ((nc = sqfl_fgetc(sqs->sfl)) == EOF)
        return(0);
      else
        {                /* next char is sequence type */
        sqs->lstinchr = (char) nc;
        sqs->sqtos = sqfl_chr2tos(sqs->flfmt,(char) nc);
        if ((nc = sqfl_fgetc(sqs->sfl)) == EOF)
          return(0);
        else
          {
          sqs->stopol = sqfl_chr2topol(sqs->flfmt,(char) nc);
          nl = bas_fgetatm_ufn(sqs->sfl,sqs->seqnam,SQNAMEMAXLEN,";\n\r",
                                 sqfl_fgetc);
          (void) sqfl_cpy2eol(sqs,sqs->annot,sqs->abuflen,sqfl_fgetc);
          np = sqs->annot;
          sqs->abuflen = strlen(np);
          sqs->annot = bas_strdup(np);
          memfree(np);      
          return(nl);
          }
        }
    break;
  case SFMT_molgen:
    while (((nc = sqfl_fgetc(sqs->sfl)) == ';') && (nc != EOF))
      if (sqfl_look4chrs(sqs->sfl,"\n\r") == EOF)
        return(0);
    if (nc == EOF)
      return(0);
    np = sqs->seqnam;
    if (np == NULL)
      return(1);
    while ((((char) nc != '\n') && ((char) nc != '\r')) &&
            ((np - sqs->seqnam) < SQNAMEMAXLEN))
      {
      *np++ = (char) nc;
      if ((nc = sqfl_fgetc(sqs->sfl)) == EOF)
        return(0);
      }
    if (((char) nc == '\n') || ((char) nc == '\r'))
 /* timed out through line terminator */
      {
      *np = '\0';
      return(1);
      }
    else
      if (sqfl_look4chrs(sqs->sfl,"\n\r") != EOF)
        return(1);
      else
        return(0);
    break;
  case SFMT_gcg:              /* look for .. */
    (void) strncpy(sqs->seqnam,sqs->filnam,SQNAMEMAXLEN);
    wlu_initlustrct(&gcglu,WLU_CASEDEP,SQFL_gcgunknown);
    wlu_addwrd(&gcglu,"Length:",SQFL_gcglength,NULL);
    wlu_addwrd(&gcglu,"Type:",SQFL_gcgtype,NULL);
    wlu_addwrd(&gcglu,"..",SQFL_gcgdots,NULL);
    wlu_maktoklu(&toklu[0]," \t\n");
    while (wlu_gettoken(sqs->sfl,gcgwrd,33,&toklu[0]))
      switch (wlu_chkwrd(&gcglu,gcgwrd))
        {
        case SQFL_gcgdots:      /* header scan completed ... */
          wlu_clrlustrct(&gcglu);
          return(1);
          break;
        case SQFL_gcglength:
          fscanf(sqs->sfl,"%d",&sqs->slength);
          break;
        case SQFL_gcgtype:
          if (wlu_gettoken(sqs->sfl,gcgwrd,33,&toklu[0]))
            sqs->sqtos = sqfl_chr2tos(SFMT_gcg,gcgwrd[0]);
          break;
        case SQFL_gcgunknown:
        default:
          break;
        }
    wlu_clrlustrct(&gcglu);
    return(0);              /* no .. by EOF, return false */
    break;
  case SFMT_staden:           /* do nothing - expect no header */
  case SFMT_undefined:
  case SFMT_raw:
  default:
    if (sqs->seqnam != NULL)
      memfree(sqs->seqnam);
    sqs->seqnam = bas_strdup(basename(sqs->filnam));
    return(1);
    break;
  }
return(0);
}

char *sqfl_valid4fmt(SFMT_TYPE fmt)
  /* return pointer to string of valid residues for fmt type of file.
  This is written as case-independent - uppercase only being shown */
{
switch (fmt)
  {
  case SFMT_fasta:
  case SFMT_nbrf:
  case SFMT_molgen:
    return("ABCDEFGHIJKLMNOPQRSTUVWXYZ-");
    break;
  case SFMT_gcg:
    return("ABCDEFGHIJKLMNOPQRSTUVWXYZ*");
    break;
  case SFMT_staden:
  case SFMT_undefined:
  case SFMT_raw:
  default:
    return("ABCDEFGHIJKLMNOPQRSTUVWXYZ-1234567890*");
    break;
  }
}  

char sqfl_normfillc(SFMT_TYPE sf)
  /* return the normal fill char used for this format */
{
switch (sf)
  {
  case SFMT_gcg:
    return('N');
    break;
  case SFMT_fasta:
  case SFMT_nbrf:
  case SFMT_molgen:
  case SFMT_staden:
  case SFMT_undefined:
  case SFMT_raw:
  default:
    return('-');
    break;
  }
}

char sqfl_getnxtres(SQFL_STRCT *sqs)
  /* return the next valid residue from sqs.  EOF or non-alphanumeric
   return NULL.  invalid chars are skipped */
{
int nc;
SQFL_SQTOPOL tc;

if ((nc = sqfl_fgetc(sqs->sfl)) == EOF)
  return('\0');
else
  if (sqs->validlut[nc])
    return(sqs->lstinchr = (char) nc);
  else
    {  /* may be topology char for molgen file, or end of nbrf/fasta, etc... */
    switch (sqs->flfmt)
      {
      case SFMT_molgen:
        if ((tc = sqfl_chr2topol(sqs->flfmt,(char) nc)) != SQTP_unknown)
          {
          sqs->stopol = tc;
          return('\0');
          }
        break;
      case SFMT_nbrf:
      case SFMT_fasta:
        switch ((char) nc)
          {
          case '*':       /* EOS char - may not be present tho */
            sqs->lstinchr = '\0';
            return('\0');
            break;
          case '>':       /* next entry - must note this */
            sqs->lstinchr = '>';
            return('\0');
            break;
          default:
            sqs->lstinchr = nc;
            break;
            }
        break;
      case SFMT_staden:
        if ((char) nc == '@')
          return('\0');
        break;
      case SFMT_raw:
        if (nc == EOF)
          return('\0');
        break;
      case SFMT_gcg:
      case SFMT_undefined:
      default:
        break;
      }
    return(sqfl_getnxtres(sqs));
    }
}

int loadsrcrng(SQFL_STRCT *sqs,
               char *seqbuf,
               int fstrt,
               int fstop)
/* read successive chars from sfl, write valid bases in fstrt..fstop into
  seqbuf if it is not NULL.  fstrt&fstop are sequence positions (1..n)
  Don't put in terminal null char.
  return no of chars that would be inserted (whether stored or not) */
{
char *sp;
char nb;
int sc;
int icnt;

icnt = sc = 0;
if (sqfl_skipsqflhdr(sqs))
  {
  sp = seqbuf;
  while ((nb = sqfl_getnxtres(sqs)) != '\0')
    if (++sc > fstop)
      return(icnt);
    else
      if (sc >= fstrt)
        {
        if (seqbuf)
          *sp++ = nb;
        icnt++;
        }
  }
return(icnt);
}

int loadsrcsq(SQFL_STRCT *sqs,
              char *seqbuf)
/* read successive chars from sfl, write valid bases into seqbuf if it is 
  not NULL.  Don't put in terminal null char.
  return highest position read (whether stored or not) */
{
char *sp;
char nb;
int sc;

sc = 0;
if (sqfl_skipsqflhdr(sqs))
  {
  sp = seqbuf;
  while ((nb = sqfl_getnxtres(sqs)) != '\0')
    {
    sc++;
    if (seqbuf)
      *sp++ = nb;
    }
  }
return(sc);
}

void sqfl_clssqstrct(SQFL_STRCT *sqs)
  /* close this sequence data structure, freeing the allocated memory */
{
if (sqs != NULL)
  {
  if (sqs->sfl != NULL)
    fclose(sqs->sfl);
  if (sqs->gfl != NULL)
    {      /* must close aux file and delete it */
    fclose(sqs->gfl);
    remove(sqs->gnam);
    memfree(sqs->gnam);
    }
  if (sqs->thisuse != NULL)
    memfree(sqs->thisuse);
  if (sqs->filnam != NULL)
    memfree(sqs->filnam);
  if (sqs->seqnam != NULL)
    memfree(sqs->seqnam);
  if (sqs->annot != NULL)
    memfree(sqs->annot);
  memfree(sqs);
  }
}

void sqfl_clsstrctfil(SQFL_STRCT *sqs)
  /* close the file component of *sqs, retaining other information */
{
if (sqs != NULL)
  if (sqs->sfl != NULL)
    {
    fclose(sqs->sfl);
    sqs->sfl = NULL;
    }
}

SQFL_STRCT *sqfl_creatsqstrctann(FILE *ufl,
                                 char *fnam,
                                 SFMT_TYPE sfmt,
                                 char *myuse,
                                 char *annotat)
/* use an open sequence file fnam, returning pointer to a new structure
  if successful.  NULL if not.
  annotat is pointer to annotation string (most useful for output files) to
  be used if non-NULL.  myuse is as for fopen */
{
SQFL_STRCT *sqs;
int cp;
char *vld;

if (ufl == NULL)
  return(NULL);
else
  {
  sqs = (SQFL_STRCT *) getmemory(sizeof(SQFL_STRCT),
                                   "Sequence data structure");
  sqs->sfl = ufl;
  sqs->flfmt = sfmt;
  if (fnam == NULL)
    sqs->filnam = bas_strdup("Null-name");
  else
    sqs->filnam = bas_strdup(fnam);
  sqs->seqnam = (char *) getmemory((SQNAMEMAXLEN+1),"Seq name buffer");
  sqs->stopol = SQTP_unknown;
  sqs->sqtos = SQFL_tosunknown;
  sqs->lstinchr = '\0';
  sqs->lnlen = sqfl_linelength(sfmt);
  vld = sqfl_valid4fmt(sfmt);
  sqs->validlut[0] = 0;
  sqs->gcgchksum = sqs->scnt = 0;
  if (annotat != NULL)
    {
    sqs->annot = bas_strdup(annotat);
    sqs->abuflen = strlen(sqs->annot) + 1;;
    }
  else
    {
    sqs->annot = (char *) getmemory((MAXPATHLEN + 1),"Annotation");
    sqs->abuflen = MAXPATHLEN + 1;
    *sqs->annot = '\0';
    }
  for (cp = 1; cp < MAXSQCHRVAL; cp++)
    sqs->validlut[cp] = (index(vld,(char) cp) != NULL) ||
                          (index(vld,toupper((char) cp)) != NULL);
  sqs->gfl = NULL;
  sqs->gnam = NULL;
  sqs->thisuse = bas_strdup(myuse);
  switch (sfmt)
    {
    case SFMT_gcg:
      if (bas_sharchrs(sqs->thisuse,"w") != NULL)
        {
        sqs->gnam = (char *) getmemory((strlen(sqs->filnam)+2),
                                         "Aux seq. file name");
        *(sqs->gnam) = '\0';
        (void) strcat(sqs->gnam,sqs->filnam);
        (void) strcat(sqs->gnam,"_");
        if ((sqs->gfl = fopen(sqs->gnam,"w+")) == NULL)
          {
          sqfl_clssqstrct(sqs);
          return(NULL);
          }
        }
      break;
    case SFMT_nbrf:
    case SFMT_fasta:
    case SFMT_molgen:
    case SFMT_raw:
    default:
      break;
    }
  return(sqs);
  }
}  

SQFL_STRCT *sqfl_opnsqstrctann(char *fnam,
                               SFMT_TYPE sfmt,
                               char *myuse,
                               char *annotat)
/* attempt to open a sequence file fnam, returning pointer to a new structure
  if successful.  NULL if not.  myuse is a normal parameter to fopen, though
  it is likely that "r" and "w" are the most useful of them,
  annotat is pointer to annotation string (most useful for output files) to
  be used if non-NULL */
{
FILE *nf;

if (fnam == NULL)
  return(NULL);
else
  if ((nf = fopen(fnam,myuse)) == NULL)
    return(NULL);
  else
    return(sqfl_creatsqstrctann(nf,fnam,sfmt,myuse,annotat));
}

SQFL_STRCT *sqfl_opnsqstrct(char *fnam,
                            SFMT_TYPE sfmt,
                            char *myuse)
/* attempt to open a sequence file fnam, returning pointer to a new structure
  if successful.  NULL if not.  myuse is a normal parameter to fopen, though
  it is likely that "r" and "w" are the most useful of them */
{
return(sqfl_opnsqstrctann(fnam,sfmt,myuse,NULL));
}  

SQFL_STRCT *sqfl_opntoread(char *fnam,
                           SFMT_TYPE sfmt)
/* attempt to open a sequence file fnam for reading, returning pointer to a 
  new structure if successful, NULL if not.  Perform header skipping, so
  seq can be read directly */
{
SQFL_STRCT *sqs;

if ((sqs = sqfl_opnsqstrct(fnam,sfmt,"r")) != NULL)
  if (sqfl_skipsqflhdr(sqs))
    return(sqs);
  else   /* bombed in skipping header, clear out sqs and return NULL */
    {
    sqfl_clssqstrct(sqs);
    return(NULL);
    }
else
  return(NULL);
}  

void sqfl_rewind(SQFL_STRCT *sqs)
{
if ((sqs != NULL) && (sqs->sfl != NULL))
  rewind(sqs->sfl);
}

off_t sqfl_filelength(FILE *fl)
  /* use fstat to return the length of file fl in bytes */
{
struct stat st;

if (fstat(fileno(fl),&st) < 0)
  {
  fprintf(stderr,"Can't stat file, length unknown\n");
  return(0);
  }
else
  return(st.st_size);
}

int readsrcsq(SQFL_STRCT *sqs,
              char *seqbuf)
/* read successive chars from sqs, write valid bases into seqbuf if it is 
  not NULL.
  return highest position read (whether stored or not) */
{
int sc;

sc = loadsrcsq(sqs,seqbuf);
if (seqbuf)
  *(seqbuf+sc) = '\0';
return(sc);
}

char *sqfl_tos2strng(SQFL_SQTOS st)
  /* return pointer to character string defining st */
{
switch (st)
  {
  case SQFL_dna:
    return("Nucleic acid: DNA");
    break;
  case SQFL_rna:
    return("Nucleic acid: RNA");
    break;
  case SQFL_peptide:
    return("Peptide");
    break;
  case SQFL_tosunknown:
  default:
    return("unknown");
    break;
  }
}

SFMT_TYPE sqfl_prmpt4sfmt(SFMT_TYPE def,
                          SFMT_TYPE rngex,
                          SFMT_TYPE rngto)
/* prompt for a sequence file format, offer default def, if defined */
{
SFMT_TYPE ff;
char ubuf[9];

for (ff = rngex; ff <= rngto; ff++)
  {
  fprintf(stdout,"%c=%s",sqfl_fmttp2chr(ff),sqfl_fmttp2strng(ff));
  if (ff < rngto)
    fputs(", ",stdout);
  }
fputs("\nFile type ",stdout);
if (def != SFMT_undefined)
  fprintf(stdout,"[%c] ",sqfl_fmttp2chr(def));
if (bas_ugetstr(">",&ubuf[0],8) <= 0)
  return(def);
else
  if (((ff = sqfl_chr2fmttp(ubuf[0])) >= rngex) && (ff <= rngto))
    return(ff);
  else
    {
    fprintf(stdout,"No such sequence type: '%c'\n",ubuf[0]);
    return(sqfl_prmpt4sfmt(def,rngex,rngto));
    }
}

char *sqfl_defextnsns(SFMT_TYPE sf)
  /* return a suggested extension for type sf */
{
switch (sf)
  {
  case SFMT_molgen:
    return(".seq");
    break;
  case SFMT_nbrf:
    return(".nbr");
    break;
  case SFMT_fasta:
    return(".fst");
    break;
  case SFMT_gcg:
    return(".gcg");
    break;
  case SFMT_raw:
  case SFMT_undefined:
  case SFMT_staden:
  default:
    return(".dat");
    break;
  }
}

void sqfl_sfmts2bf(char *ubuf,
                   int blim)
/* write sfmt details to ubuf */
{
SFMT_TYPE sfp;
char *bp;

bp = ubuf;
for (sfp = SFMT_staden; (sfp <= SFMT_gcg) && ((bp - ubuf) < blim); sfp++)
  {
  if (sfp > SFMT_staden)
    {
    *bp++ = ',';
    *bp++ = ' ';
    }
  *bp++ = sqfl_fmttp2chr(sfp);
  *bp++ = '=';
  *bp++ = '\0';
  strcat(ubuf,sqfl_fmttp2strng(sfp));
  while ((*bp != '\0') & ((bp - ubuf) < blim))
    bp++;
  }
}

int sqfl_wrtsq2fl(SQFL_STRCT *sfl,
                  char *ubuf,
                  int blen,
                  int *lcnt)
/* put blen chars of ubuf as sequence to sfl file.  Return res count */
{
char *bp;
int rc;

bp = ubuf;
rc = *lcnt = 0;
while (rc <= blen)
  {
  sqfl_putres(sfl,*bp,rc,lcnt);
  bp++;
  rc++;
  }
return(rc);
}

int sqfl_wrtsqbuf2fl(SQFL_STRCT *sfl,
                     char *sqnam,
                     char *maker,
                     char *origin,
                     char *ubuf,
                     int blen)
/* write ubuf up to blen out to sfl, (previously opened).  Do headers & 
terminators.  Don't close sfl, leave to calling routine.  Return no res
written */
{
int rc;
int lcnt;

sqfl_headsfstrct(sfl,sqnam,maker,origin);
rc = sqfl_wrtsq2fl(sfl,ubuf,blen,&lcnt);
sqfl_termsqfl(sfl,maker,&lcnt);
return(rc);
}
