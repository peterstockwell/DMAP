/* bin_cnts: to take generate binned lists of counts from 
appropriately-formatted input files (chrno, position, strand/status
where strand/status = '+' or '-')

Process requires scanning human genome, chromosome by chromosome
to get lengths.  Linked lists of bins are used to manage
variable Chr lengths & variable counts */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "mrg_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"
#include "bin_cnts.h"
#include "fsm_ops.h"

/* global debug & out style variables, for simplicity of access */
RBC_DBG debuglevel;
BC_OUTMODE omode;
RBC_SRC_STYLE srcstyle;
int glblreadlen;          /* global value for readlength */
/* BC_CHR_BIN *chrbinlsts[ChrY]; */ /* bin lists for each chromosome */
/* BC_CHR_BIN *chrbinsalt[ChrY]; */  /* alternate set of character bins */

int err_msg(char *fmt,
            ...)
/* write user error message.  Return 0 for err return status */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
return(0);
}

void err_msg_die(char *fmt,
                 ...)
/* write user error message then exit with error status */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
exit(1);
}

void say_usage(FILE *fl,
               char *pnam)
{
fprintf(fl,"%s: create binned counts for chromosomal positions\n",pnam);
fputs("Options:\n",fl);
fputs("     -r <posfile> read <posfile> as set of chr posit strand/meth\n",fl);
fputs("     -R <posfile2> as -r but for 2nd position file\n",fl);
fprintf(fl,"     -b <binlength>: set bin length (def=%d)\n",BC_DEF_BINLEN);
fputs("     -g <genomehead> dir and file string to locate genomic seq files by adding n.fa\n",fl);
fputs("     -l/-L list bins (-L=>only nonzero bins)\n",fl);
fputs(
"     -m scan for diff meth regions restricted rep (Li, et al. (2010) PLOSBiology,11,e1000533)\n"
,fl);
fputs("     -M <j,k> scan for restricted rep fragment sizes between j & k residues. Make bins\n",fl);
fputs("     -N as for -M but list bins to stdout\n",fl);
fputs("     -k as for -M, but note reads which don't map into restricted rep bins\n",fl);
fputs("     -K as for -k, but only print totals for meth & unmeth counts\n",fl);
fputs("     -S <dirheader> write .dat files to <dirheader> for RR genome for SeqMonk\n",fl);
fputs("     -c <n> restrict effort to Chromosome <n> (def = all chromosomes)\n",fl);
fputs("     -C <n> restrict bins to those with <n> or more CpGs\n",fl);
fputs("     -A attempt to amalgamate restr rep regions that might otherwise fail CpG criteria\n",fl);
fputs("     -x <excludefile> exclude regions (fmt: Chrno regionstart regionend)\n",fl);
}

BC_CHR_BIN *rbc_appndbin(BC_CHR_BIN **lstrt,
                          int bstart,
                          int len)
/* create and append a new element to *lstrt,
init count to 0.
 Return address of new element */
{
BC_CHR_BIN *prev, *end_ptr;

if (lstrt != NULL)
  {
/*  if (debuglevel > RBC_dbg_none)
    fprintf(stdout,"Bin %d..%d, len=%d\n",bstart,bstart+len-1,len); */
  prev = end_ptr = *lstrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtbin;
    }
  end_ptr = (BC_CHR_BIN *) getmemory(sizeof(BC_CHR_BIN),"bin elt");
  end_ptr->nxtbin = NULL;
  end_ptr->spos = bstart;
  end_ptr->binlen = len;
  end_ptr->bcntrev = end_ptr->bcntfwd = end_ptr->cpgcnt = 0;
  if (*lstrt == NULL)
    {
    *lstrt = end_ptr;
    end_ptr->prvbin = NULL;
    }
  else
    {
    prev->nxtbin = end_ptr;
    end_ptr->prvbin = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void rbc_delcntbin(BC_CHR_BIN *ep,
                   BC_CHR_BIN **lstrt)
/* delete ep from list *lstrt */
{
BC_CHR_BIN *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvbin) == NULL)
    *lstrt = ep->nxtbin;
  else
    pt->nxtbin = ep->nxtbin;
  if ((pt = ep->nxtbin) != NULL)
    pt->prvbin = ep->prvbin;
  memfree(ep);
  }
}

void rbc_clrallcntbins(BC_CHR_BIN **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  rbc_delcntbin(*lstrt,lstrt);
}

int bc_cntcntbins(BC_CHR_BIN *clst)
  /* recursively read list elements */
{
if (clst == NULL)
  return(0);
else
  return(bc_cntcntbins(clst->nxtbin) + 1);
}

int bc_sumcntfwdbins(BC_CHR_BIN *blst)
  /* recursively sum the counts in blst */
{
if (blst == NULL)
  return(0);
else
  {
  if (debuglevel > RBC_dbg_on)
    fprintf(stdout,"%d %d %lx\n",blst->spos,blst->bcntfwd,(long int) blst->nxtbin);
  return(blst->bcntfwd + bc_sumcntfwdbins(blst->nxtbin));
  }
}

int bc_itsumcntfwdbins(BC_CHR_BIN *blst)
  /* iteratively sum counts in blst... */
{
BC_CHR_BIN *bp;
int cnt;

bp = blst;
cnt = 0;
while (bp != NULL)
  {
  cnt += bp->bcntfwd;
  bp = bp->nxtbin;
  }
return(cnt);
}

int bc_itsumcntrevbins(BC_CHR_BIN *blst)
  /* iteratively sum counts in blst... */
{
BC_CHR_BIN *bp;
int cnt;

bp = blst;
cnt = 0;
while (bp != NULL)
  {
  cnt += bp->bcntrev;
  bp = bp->nxtbin;
  }
return(cnt);
}

BC_CHR_BIN *cpg_lastcntbin(BC_CHR_BIN *clst)
  /* iterate thru clst, returning
last element, NULL if none */
{
BC_CHR_BIN *ep;

if ((ep = clst) == NULL)
  return(NULL);
else
  {
  while (ep->nxtbin != NULL)
    ep = ep->nxtbin;
  return(ep);
  }
}

int bc_sumbinlens(BC_CHR_BIN *blst)
  /* recursively sum blst */
{
if (blst == NULL)
  return(0);
else
  return(bc_sumbinlens(blst->nxtbin) + blst->binlen);
}

int int_in_rng(int b1,
               int x,
               int b2)
/* return 1 if b1 <= x <= b2 or
b2 <= x <= b1 */
{
if (b1 <= b2)
  return((b1 <= x) && (x <= b2));
else
  return(int_in_rng(b2,x,b1));
}

BC_REGN_ELT *bc_appndrgnelt(BC_REGN_ELT **lstrt,
                            int rgstart,
                            int rgstop)
/* create and append a new element to *lstrt,
 Return address of new element */
{
BC_REGN_ELT *prev, *end_ptr;

if (lstrt != NULL)
  {
  prev = end_ptr = *lstrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtregn;
    }
  end_ptr = (BC_REGN_ELT *) getmemory(sizeof(BC_REGN_ELT),"region elt");
  end_ptr->nxtregn = NULL;
  end_ptr->rstart = rgstart;
  end_ptr->rstop = rgstop;
  if (*lstrt == NULL)
    {
    *lstrt = end_ptr;
    end_ptr->prvregn = NULL;
    }
  else
    {
    prev->nxtregn = end_ptr;
    end_ptr->prvregn = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void bc_delregnelt(BC_REGN_ELT *ep,
                    BC_REGN_ELT **lstrt)
/* delete ep from list *lstrt */
{
BC_REGN_ELT *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvregn) == NULL)
    *lstrt = ep->nxtregn;
  else
    pt->nxtregn = ep->nxtregn;
  if ((pt = ep->nxtregn) != NULL)
    pt->prvregn = ep->prvregn;
  memfree(ep);
  }
}

void bc_clrallregnelts(BC_REGN_ELT **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  bc_delregnelt(*lstrt,lstrt);
}

int bc_cntregnelts(BC_REGN_ELT *clst)
  /* recursively count list elements */
{
if (clst == NULL)
  return(0);
else
  return(bc_cntregnelts(clst->nxtregn) + 1);
}

BC_REGN_ELT *bc_regnelt4pos(BC_REGN_ELT *rlst,
                            int pos)
/* return the region element corresponding to pos, if any */
{
BC_REGN_ELT *rp;

rp = rlst;
while (rp != NULL)
  if (int_in_rng(rp->rstart,pos,rp->rstop))
    return(rp);
  else
    rp = rp->nxtregn;
/* fell off end, reurn NULL */
return(NULL);
}

BC_CHR_BIN *bc_cntbin4pos(BC_CHR_BIN *blst,
                          int posn)
/* return a pointer to the bin that contains posn.
NULL if none */
{
BC_CHR_BIN *bp;

bp = blst;
while (bp != NULL)
  if (int_in_rng(bp->spos,posn,(bp->spos + bp->binlen -1)))
    return(bp);
  else
    bp = bp->nxtbin;
/* fell off end, return NULL */
return(NULL);
}

int bc_readsrcfl(BC_CHR_BIN *bnlsts[],
                 FILE *sfl,
                 RBC_CHRNO uchrno,
                 BC_REGN_ELT *exclud[],
                 RBC_CHRNO mxchr)
/* use fscanf to read successive lines from sfl.
If uchrno is nonzero, then only do that chromosome */
{
char chrstr[5];
int sqpos;
char sensestr[5];
int matcnt;
int chrno;
int scnt;
BC_CHR_BIN *bp;
int prvchrno;
BC_REGN_ELT *xcld;

matcnt = 0;

#ifdef LEGACY_MODE

while ((scnt = fscanf(sfl,"%s %d %s",&chrstr[0],&sqpos,&sensestr[0])) != EOF)
  if ((scnt == 3) && ((chrno = rbc_str2chrno(&chrstr[0])) > Chr_unk) &&
        (chrno <= mxchr) && ((uchrno == 0) || (chrno == uchrno)))
    {
    if (((xcld = bc_regnelt4pos(exclud[chrno-1],sqpos)) == NULL) &&
         ((bp = bc_cntbin4pos(bnlsts[chrno-1],sqpos)) != NULL))
      {   
      if (sensestr[0] == '+')
        bp->bcntfwd++;
      else
        bp->bcntrev++;
      matcnt++;
      }
    }

#else    /* LEGACY_MODE */

prvchrno = 0;
bp = NULL;
while ((scnt = fscanf(sfl,"%s %d %s",&chrstr[0],&sqpos,&sensestr[0])) != EOF)
  if ((scnt == 3) && ((chrno = rbc_str2chrno(&chrstr[0])) > Chr_unk) &&
        (chrno <= mxchr) && ((uchrno == 0) || (chrno == uchrno)))
      {
      if (prvchrno != chrno)
        {
        bp = NULL;
        prvchrno = 0;
        }
      if ((bp != NULL) && !int_in_rng(bp->spos,sqpos,(bp->spos + bp->binlen -1)))
        {
        bp = NULL;
        prvchrno = 0;
        }
      xcld = bc_regnelt4pos(exclud[chrno-1],sqpos);
      if ((bp == NULL) && (xcld == NULL))
        bp = bc_cntbin4pos(bnlsts[chrno-1],sqpos);
      if ((bp != NULL) && (xcld == NULL))
        {
        if (sensestr[0] == '+')
          bp->bcntfwd++;
        else
          bp->bcntrev++;
        matcnt++;
        prvchrno = chrno;
        }
      }

#endif    /* LEGACY_MODE (else) */

return(matcnt);
}

int bc_chknreadsrcfl(BC_CHR_BIN *bnlsts[],
                     FILE *sfl,
                     RBC_CHRNO uchrno,
                     BC_REGN_ELT *exclud[],
                     RBC_CHRNO mxchr)
/* check sfl for openness, then call bc_readsrcfl, returning number
elements read */
{
if (sfl == NULL) /* can't do anything */
  return(0);
else
  return(bc_readsrcfl(bnlsts,sfl,uchrno,exclud,mxchr));
}

int rbc_getgenomesizs(char *hdrstr,
                      int seqlens[],
                      int chrmax,
                      RBC_CHRNO uchrno)
/* use hdrstr to create a series of file names, one
for each chromosome, open each as a Fasta sequence
file and scan for the sequence length.  Return the
number of chromosomes processed. if uchrno is nonzero
then restrict activity to that alone */
{
char *sqfnam;
int chno;
int nblen;
SQFL_STRCT *chsqfl;
int rcnt;

rcnt = 0;
sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq file name buf");
for (chno = 1; chno <= chrmax; chno++)
  {
  seqlens[chno-1] = 0;
  if ((uchrno == 0) || (chno == uchrno))
    {
    snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,rbc_chrno2str((RBC_CHRNO) chno,1));
    if ((chsqfl = sqfl_opnsqstrct(sqfnam,SFMT_fasta,"r")) != NULL)
      {
      seqlens[chno-1] = readsrcsq(chsqfl,NULL);
      if (debuglevel > RBC_dbg_on)
        {
        fprintf(stdout,"%s: %d res\n",sqfnam,seqlens[chno-1]);
        fflush(stdout);
        }
      sqfl_clssqstrct(chsqfl);
      rcnt++;
      }
    else
      err_msg("Can't open chromosome file %s\n",sqfnam);
    }
  }
return(rcnt);
}

int rbc_getgenomesqs(char *hdrstr,
                     char *seqsarr[],
                     int seqlens[],
                     int chrmax,
                     RBC_CHRNO uchrno)
/* use hdrstr to create a series of file names, one
for each chromosome, open each as a Fasta sequence
file and create a buffer for each.  Return the
number of chromosomes processed */
{
char *sqfnam;
int chno;
int nblen;
SQFL_STRCT *chsqfl;
int rcnt;

rcnt = 0;
sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq file name buf");
for (chno = 1; chno <= chrmax; chno++)
  if ((uchrno == 0) || (uchrno == chno))
    {
    snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,rbc_chrno2str((RBC_CHRNO) chno,1));
    if ((chsqfl = sqfl_opnsqstrct(sqfnam,SFMT_fasta,"r")) != NULL)
      {
      seqlens[chno-1] = readsrcsq(chsqfl,NULL);
      if (debuglevel > RBC_dbg_on)
        {
        fprintf(stdout,"%s: %d res,",sqfnam,seqlens[chno-1]);
        fflush(stdout);
        }
      seqsarr[chno-1] = (char *) getmemory(seqlens[chno-1]+1,"Chr Buff");
      sqfl_rewind(chsqfl);
      (void) readsrcsq(chsqfl,seqsarr[chno-1]);
      sqfl_clssqstrct(chsqfl);
      if (debuglevel > RBC_dbg_on)
        {
        fprintf(stdout,"'%.10s...'\n",seqsarr[chno-1]);
        fflush(stdout);
        }
      rcnt++;
      }
    else
      {
      err_msg("Can't open chromosome file %s\n",sqfnam);
      seqlens[chno-1] = 0;
      }
    }
  else
    {
    seqlens[chno-1] = 0;
    seqsarr[chno-1] = NULL;
    }
return(rcnt);
}

char tr_int2nares(int iv)
  /* return a nucleic acid residue for iv */
{
switch (iv)
  {
  case 0:
    return('a');
    break;
  case 1:
    return('c');
    break;
  case 2:
    return('g');
    break;
  case 3:
    return('t');
    break;
  default:
    return('?');
    break;
  }
}

int tr_nares2int(char res)
  /* return an int value 0..3 for res, -1 for unknown */
{
switch (toupper(res))
  {
  case 'A':
    return(0);
    break;
  case 'C':
    return(1);
    break;
  case 'G':
    return(2);
    break;
  case 'T':
  case 'U':
    return(3);
    break;
  default:
    return(-1);
    break;
  }
}

int rbc_scangenomesqs(char *hdrstr,
                      int chrmax,
                      RBC_CHRNO uchrno,
                      FS_FSMSTRCT *cpgfsmp)
/* use hdrstr to create a series of file names, one
for each chromosome, open each as a Fasta sequence
file and create a buffer for each.  Return the
number of chromosomes processed */
{
char *sqfnam;
int chno;
int nblen;
SQFL_STRCT *chsqfl;
int rcnt;
char nxtres;
int cpos;
FS_RESELT *frp;
FS_DATPRELT *dpep;

rcnt = 0;
sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq file name buf");
for (chno = 1; chno <= chrmax; chno++)
  if ((uchrno == 0) || (uchrno == chno))
    {
    snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,rbc_chrno2str((RBC_CHRNO) chno,1));
    if (((chsqfl = sqfl_opnsqstrct(sqfnam,SFMT_fasta,"r")) != NULL) &&
            sqfl_skipsqflhdr(chsqfl))
      {
      fs_initrun(cpgfsmp);
      cpos = 0;
      while ((nxtres = sqfl_getnxtres(chsqfl)) != '\0')
        {
        cpos++;
        if ((frp = fs_procchr(cpgfsmp,nxtres,tr_nares2int)) != NULL)
          {
          dpep = (FS_DATPRELT *) frp->action;
          fprintf(stdout,"C%s: %d\n",rbc_chrno2str(chno,1),
                    (cpos-dpep->ldstr+1));
          }
        }
      rcnt++;
      }
    else
      err_msg("Can't open chromosome file %s\n",sqfnam);
    }
memfree(sqfnam);
return(rcnt);
}

FILE *bc_opnsqmonkdatfl(char *dirhdr,
                        RBC_CHRNO cno,
                        char *ext)
/* attempt to create a destination file for cno
with dirhdr + ext and open for writing */
{
FILE *ofl;
char *dfname;

if (asprintf(&dfname,"%s%s%s",((dirhdr!=NULL)?dirhdr:""),
               rbc_chrno2str(cno,1),((ext!=NULL)?ext:".dat")) > 0)
  {
  ofl = fopen(dfname,"w");
  free(dfname);
  return(ofl);
  }
else
  return(NULL);
}

int bc_redrepgenomsqs(char *hdrstr,
                      int chrmax,
                      RBC_CHRNO uchrno,
                      FS_FSMSTRCT *cpgfsmp,
                      BC_CHR_BIN *chrbins[],
                      char *seqsarr[],
                      int seqlens[])
/* use hdrstr to create a series of file names, one
for each chromosome, open each as a Fasta sequence
file.  Scan with cpgfsmp (CCGG) and create bins, giving
a count of CpG found with CG.  Return the
number of chromosomes processed.
if seqsarr[] & seqlens != NULL, then
store seq therein and set final length */
{
char *sqfnam;
int chno;
int nblen;
SQFL_STRCT *chsqfl;
int rcnt;
char nxtres;
int cpos;
FS_RESELT *frp;
FS_DATPRELT *dpep;
int thismat;
int prvmat;
int thislen;
BC_CHR_BIN *binsend;
int cpgs;

rcnt = 0;
sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq file name buf");
for (chno = 1; chno <= chrmax; chno++)
  if ((uchrno == 0) || (uchrno == chno))
    {
    snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,rbc_chrno2str((RBC_CHRNO) chno,1));
    if ((chsqfl = sqfl_opnsqstrct(sqfnam,SFMT_fasta,"r")) != NULL)
      {
      thislen = readsrcsq(chsqfl,NULL);
      if (seqsarr != NULL)
        seqsarr[chno-1] = (char *) getmemory(thislen+1,"SeqBuf");
      if (seqlens != NULL)
        seqlens[chno-1] = thislen;
      sqfl_rewind(chsqfl);
      (void) sqfl_skipsqflhdr(chsqfl);
      fs_initrun(cpgfsmp);
      cpos = 0;
      prvmat = 1;
      binsend = NULL;
      cpgs = 0;
      while ((nxtres = sqfl_getnxtres(chsqfl)) != '\0')
        {
        if (seqsarr != NULL)
          {
          *(seqsarr[chno-1] + cpos) = nxtres;
          *(seqsarr[chno-1] + cpos + 1) = '\0';
          }
        cpos++;
        if ((frp = fs_procchr(cpgfsmp,nxtres,tr_nares2int)) != NULL)
          {
          dpep = (FS_DATPRELT *) frp->action;
/*          if (debuglevel > RBC_dbg_none)
             fprintf(stdout,"Match: %s @ %d\n",dpep->dstrng,cpos); */
          if (dpep->ldstr == 4)      /* MspI site */
            {
            thismat = cpos-dpep->ldstr+2;  /* MspI cuts C/CGG */
            thislen = thismat - prvmat;
/*            if (debuglevel > RBC_dbg_none)
              fprintf(stdout,"saving: %d %d\n",prvmat,thislen); */
            binsend = rbc_appndbin(&binsend,prvmat,thislen);
            binsend->cpgcnt = cpgs;
            cpgs = 0;
            if (chrbins[chno-1] == NULL)
              chrbins[chno-1] = binsend;
            prvmat = thismat;
            }
          else
            if (dpep->ldstr == 2)  /* CpG */
              cpgs++;
          }
        }
      rcnt++;
      }
    else
      err_msg("Can't open chromosome file %s\n",sqfnam);
    }
memfree(sqfnam);
return(rcnt);
}

void rbc_rptchrout(FILE *ofl,
                   char c,
                   int ccnt)
/* put c out ccnt times to ofl */
{
while (ccnt-- > 0)
  fputc(c,ofl);
}

double bin_coeff(int n,
                 int k)
/* return as double the binomial
coefficient of n,k.  Multiplicative
formula */
{
double prodct;
int i;

prodct = 1.0;
for (i = 1; i <= k; i++)
  prodct *= (double) (n - k + i)/((double) i);
return(prodct);
}

double fishers_exact_t1(int a,
                        int b,
                        int c,
                        int d)
/* return a value for Fisher's exact
statistic - one tailed */
{
double rval;
int tot;

tot = a + b + c + d;
rval = (bin_coeff(a+b,a) / bin_coeff(tot,a+c)) * bin_coeff(c+d,c);
return(rval);
}

int bc_getexcludregns(FILE *efl,
                      BC_REGN_ELT *eregns[],
                      RBC_CHRNO maxcno)
/* read lines from efl (fmt: chrno regstart regstop\n)
and append to appropriate chromosome no list */
{
RBC_CHRNO cno;
int rgstart;
int rgstop;
int rdval;
char cstr[32];
int rcnt;

rcnt = 0;
while ((rdval = fscanf(efl,"%s %d %d",&cstr[0],&rgstart,&rgstop)) != EOF)
  if ((rdval == 3) && ((cno = rbc_str2chrno(&cstr[0])) >= Chr1) && (cno <= ChrY))
    {
    (void) bc_appndrgnelt(&eregns[cno-1],rgstart,rgstop);
    rcnt++;
    }
return(rcnt);
}

void bc_commalst2ints(char *str,
                      int *v1,
                      int *v2)
/* expect str to contain two base10 integers separated by a comma.
attempt to read these and return as v1 & v2.  zero return on
failures */
{
char *ep;
char *fp;

*v1 = *v2 = 0;
*v1 = (int) strtol(str,&ep,10);
if (ep != str)
  *v2 = strtol(ep+1,&fp,10);
}

int bc_fragsizeok(int fragmin,
                  int fragmax,
                  int fragsize)
/* return if fragsize fragment is acceptable,
noting that zero limits will always return true. */
{
if ((fragmin != 0) && (fragmax != 0))
  return(int_in_rng(fragmin,fragsize,fragmax));
else
  return(1);
}

int bc_cpglmttest(int cpgmin,
                  int cpgcnt)
/* return 1 if cpgcnt passes cpgmin, noting
that 0 cpgmin will always return true */
{
if (cpgmin <= 0)
  return(1);
else
  return(cpgcnt >= cpgmin);
}

int bc_binundercpglmt(BC_CHR_BIN *bp,
                      int cpgmin)
/* apply minimum CpG count test to bin, noting
that cgpmin <= 0 will always fail */
{
if (bp != NULL)
  return(!bc_cpglmttest(cpgmin,bp->cpgcnt));
else
  return(1);
}

int bc_prunebins4len(BC_CHR_BIN *binlsts[],
                     RBC_CHRNO maxchr,
                     int minfrag,
                     int maxfrag)
/* scan all or selected chromosome bin lists, removing all
bins which don't conform to minfrag..maxfrag (inclusive).
Return count of deleted bins */
{
RBC_CHRNO chrno;
BC_CHR_BIN *bp;
int delcnt;
BC_CHR_BIN *nxt;

delcnt = 0;
for (chrno = Chr1; chrno <= maxchr; chrno++)
  {
  bp = binlsts[chrno-1];
  while (bp != NULL)
    {
    nxt = bp->nxtbin;
    if (!bc_fragsizeok(minfrag,maxfrag,bp->binlen))
      {
/*      if (debuglevel > RBC_dbg_none)
        fprintf(stdout,"Del: %d..%d (%d), CpGs=%d\n",bp->spos,bp->spos + bp->binlen -1,bp->binlen,bp->cpgcnt); */
      rbc_delcntbin(bp,&binlsts[chrno-1]);
      delcnt++;
      }
    bp = nxt;
    }
  }
return(delcnt);
}

int bc_prunebinslenok(BC_CHR_BIN *binlsts[],
                      RBC_CHRNO maxchr,
                      int minfrag,
                      int maxfrag)
/* scan all or selected chromosome bin lists, removing all
bins which do conform to minfrag..maxfrag (inclusive).
Return count of deleted bins */
{
RBC_CHRNO chrno;
BC_CHR_BIN *bp;
int delcnt;
BC_CHR_BIN *nxt;

delcnt = 0;
for (chrno = Chr1; chrno <= maxchr; chrno++)
  {
  bp = binlsts[chrno-1];
  while (bp != NULL)
    {
    nxt = bp->nxtbin;
    if (bc_fragsizeok(minfrag,maxfrag,bp->binlen))
      {
/*      if (debuglevel > RBC_dbg_none)
        fprintf(stdout,"Del: %d..%d (%d), CpGs=%d\n",bp->spos,bp->spos + bp->binlen -1,bp->binlen,bp->cpgcnt); */
      rbc_delcntbin(bp,&binlsts[chrno-1]);
      delcnt++;
      }
    bp = nxt;
    }
  }
return(delcnt);
}

int bc_prunebins4cpgmin(BC_CHR_BIN *binlsts[],
                        RBC_CHRNO maxchr,
                        int mincpg)
/* scan all or selected chromosome bin lists, removing all
bins which don't meet mincpg criterion.
Return count of deleted bins */
{
RBC_CHRNO chrno;
BC_CHR_BIN *bp;
int delcnt;
BC_CHR_BIN *nxt;

delcnt = 0;
for (chrno = Chr1; chrno <= maxchr; chrno++)
  {
  bp = binlsts[chrno-1];
  while (bp != NULL)
    {
    nxt = bp->nxtbin;

    if (bc_binundercpglmt(bp,mincpg))
      {
/*      if (debuglevel > RBC_dbg_none)
        fprintf(stdout,"Del: %d..%d (%d), CpGs=%d\n",bp->spos,bp->spos + bp->binlen -1,bp->binlen,bp->cpgcnt); */
      rbc_delcntbin(bp,&binlsts[chrno-1]);
      delcnt++;
      }
    bp = nxt;
    }
  }
return(delcnt);
}

int bc_losefailedbins(BC_CHR_BIN *binlsts[],
                      RBC_CHRNO maxchr,
                      int minfrag,
                      int maxfrag,
                      int cpgmin)
/* scan all or selected chromosome bin lists, removing all
bins which do conform to minfrag..maxfrag (inclusive) and
fail any cpgmin criterion.
Return count of deleted bins */
{
return(bc_prunebins4len(binlsts,maxchr,minfrag,maxfrag) +
         bc_prunebins4cpgmin(binlsts,maxchr,cpgmin));
}

int bc_adjacentbins(BC_CHR_BIN *b1,
                    BC_CHR_BIN *b2)
/* return true if b1 & b2 are adjacent with
no intervening gap */
{
if ((b1 == NULL) || (b2 == NULL))
  return(0);
else
  if (b1->spos > b2->spos)
    return(bc_adjacentbins(b2,b1));
  else
    return((b1->spos + b1->binlen) == b2->spos);
}

int bc_amlgmatelocpgbins(BC_CHR_BIN *cbinlst[],
                         RBC_CHRNO maxchr,
                         int minfrag,
                         int maxfrag,
                         int cpgmin)
/* scan cbinlst for bins which may fail a CpG count
but which could be amalgamated with adjacent bins to
pass the CpG count criterion implied by cpgmin.
Amalgamate valid bins.  Note that this may leave bins
that exceed length criteriam but if these are
made by joining valid adjacent bins, then the
size selection during restricted representation library
preparation should have left valid reads for those
regions.  Return number of bins pruned or amalgamated */
{
int abincnt;
RBC_CHRNO chrno;
BC_CHR_BIN *bp;
BC_CHR_BIN *nxt;

abincnt = bc_prunebins4len(cbinlst,maxchr,minfrag,maxfrag);
for (chrno = Chr1; chrno <= maxchr; chrno++)
  {
  bp = cbinlst[chrno-1];
  while (bp != NULL)
    {
    nxt = bp->nxtbin;
    if (bc_adjacentbins(bp,nxt) &&
          (bc_binundercpglmt(bp,cpgmin) || bc_binundercpglmt(nxt,cpgmin))
             && bc_cpglmttest(cpgmin,(bp->cpgcnt + nxt->cpgcnt)))
      {  /* checks that bins are adjacent, one or other fails cpg criterion & together they pass it */
/*          if (debuglevel > RBC_dbg_none)
            fprintf(stdout,"C%s amalgamating %d..%d (%d) & %d..%d (%d) with %d+%dCpGs\n",
                      rbc_chrno2str(chrno,1),bp->spos,(bp->spos+bp->binlen-1),bp->binlen,
                      nxt->spos,(nxt->spos+nxt->binlen-1),nxt->binlen,bp->cpgcnt,nxt->cpgcnt); */
      bp->binlen += nxt->binlen;
      bp->cpgcnt += nxt->cpgcnt;
      bp->bcntfwd += nxt->bcntfwd;
      bp->bcntrev += nxt->bcntrev;
      rbc_delcntbin(nxt,&cbinlst[chrno-1]);
      abincnt++;
      }
    bp = bp->nxtbin;
    }
/* now prune failed cpg bins */
  abincnt += bc_prunebins4cpgmin(cbinlst,maxchr,cpgmin);
  }
return(abincnt);
}

BC_CHR_BIN *bc_dupbinlst(BC_CHR_BIN *blstp)
  /* for every bin in blstp generate a new bin, init the
contents of each bin to those in blstp list */
{
BC_CHR_BIN *bp;
BC_CHR_BIN *ep;
BC_CHR_BIN *newlst;

bp = blstp;
ep = newlst = NULL;
while (bp != NULL)
  {
  ep = rbc_appndbin(&ep,bp->spos,bp->binlen);
  if (newlst == NULL)
    newlst = ep;
  ep->cpgcnt = bp->cpgcnt;
  ep->bcntfwd = bp->bcntfwd;
  ep->bcntrev = bp->bcntrev;
  bp = bp->nxtbin;
  }
return(newlst);
}

BC_CHR_BIN *bc_genlstforgaps(BC_CHR_BIN *blstp)
  /* generate a bin for each gap in blstp..., return the
start of the list */
{
BC_CHR_BIN *bp;
BC_CHR_BIN *ep;
BC_CHR_BIN *newlst;
BC_CHR_BIN *nxt;
int newstrt;
int newlen;

bp = blstp;
newlst = NULL;
ep = NULL;
while (bp != NULL)
  {
  if ((nxt = bp->nxtbin) != NULL)
    {
    if (!bc_adjacentbins(bp,nxt))
      {
      newstrt = bp->spos + bp->binlen;
      newlen = nxt->spos - newstrt;
      ep = rbc_appndbin(&ep,newstrt,newlen);
      if (newlst == NULL)
        newlst = ep;
      }  
    bp = nxt;
    }
  else
    bp = NULL;
  }
return(newlst);
}

int bc_hypometh(BC_CHR_BIN *bp)
  /* return 1 if this is not hypomethylated:
The requirement is stated by Li, et al. (2010) PLOSBiology,11,e1000533,
but is not actually defined by them.  Let's try some simple scheme
for now */
{
return(bp->bcntfwd <= 0);
}

int bc_2foldmethdiff(BC_CHR_BIN *b1p,
                     BC_CHR_BIN *b2p)
/* test if the methylation difference between
b1p & b2p is >= 2 fold */
{
float b1prop;
float b2prop;
float flddiff;

if ((b1p->bcntfwd > 0) && (b2p->bcntfwd > 0)
     && (b1p->bcntrev > 0) && (b2p->bcntrev > 0)) 
  {
  b1prop = (float) b1p->bcntfwd/b1p->bcntrev;
  b2prop = (float) b2p->bcntfwd/b2p->bcntrev;
  flddiff = b1prop/b2prop;
  return((flddiff >= 2.0) || (flddiff <= 0.5));
  }
else
  return(0);
}

void bc_scandiffmethbins(FILE *ofl,
                         BC_CHR_BIN *bins1p[],
                         BC_CHR_BIN *bins2p[],
                         RBC_CHRNO maxchr)
/* scan binsNp (N=1 & 2) which are expected to have
identical partitions.  Perform tests of
Li, et al. (2010) PLOSBiology,11,e1000533 for
differential methylation - print results to ofl */
{
BC_CHR_BIN *b1p;
BC_CHR_BIN *b2p;
RBC_CHRNO chrno;

for (chrno = Chr1; chrno <= maxchr; chrno++)
  {
  b1p = bins1p[chrno-1];
  b2p = bins2p[chrno-1];
  while ((b1p != NULL) && (b2p != NULL))
    {
    if ((b1p->spos != b2p->spos) || (b1p->binlen != b2p->binlen))
      fprintf(ofl,"Inconsistent Dmeth bins: Chr %s (%d..%d) vs (%d..%d)\n",
                rbc_chrno2str(chrno,1),b1p->spos,(b1p->spos+b1p->binlen-1),
                b2p->spos,(b2p->spos+b2p->binlen-1));
    else
      if ((b1p->cpgcnt >= 5) && !bc_hypometh(b1p) && !bc_hypometh(b2p) && bc_2foldmethdiff(b1p,b2p))
        {
        fprintf(ofl,"C%s\t%d..%d\t%d\t%d CpGs\t1: %d+ %d-\t2: %d+ %d-",
                  rbc_chrno2str(chrno,1),b1p->spos,(b1p->spos+b1p->binlen-1),
                  b1p->binlen,b1p->cpgcnt,b1p->bcntfwd,b1p->bcntrev,
                  b2p->bcntfwd,b2p->bcntrev);
        fprintf(ofl,"\tPr= %g",
                  fishers_exact_t1(b1p->bcntfwd,b1p->bcntrev,b2p->bcntfwd,b2p->bcntrev));
        fputc('\n',ofl);
        }
    b1p = b1p->nxtbin;
    b2p = b2p->nxtbin;
    }
  }
}

void bc_mksqmnkhdr(FILE *fl,
                   char *gnamstr,
                   RBC_CHRNO cno,
                   int sqlen,
                   int minfrag,
                   int maxfrag)
/* put a header section to fl in style of SeqMonk
.dat files */
{
fprintf(fl,"ID   C%sRR; RR 1; linear; genomic DNA; STD; PRO: %d BP.\n",
          rbc_chrno2str(cno,1),sqlen);
fprintf(fl,"AC   chromosome:RedRep%d_%d:Chr%s:1:%d:1\n",minfrag,maxfrag,
          rbc_chrno2str(cno,1),sqlen);
fprintf(fl,"DE   Reduced Representation (%d-%dbp) of %s Chr%s\n",
          minfrag,maxfrag,((gnamstr==NULL)?"GRCh37":gnamstr),
          rbc_chrno2str(cno,1));
fprintf(fl,"FH   Key             Location/Qualifiers\nFH\n");
}

void bc_bascnt4regn(char *seq,
                    int sqlen,
                    int rstart,
                    int rstop,
                    int bascnt[])
/* scan region rstart..rstop of seq, adding residues to
bascnt[] */
{
int sp;
SQ_RESTYPE rt;

sp = rstart;
rstop = imin(rstop,sqlen);
while (sp <= rstop)
  {
  rt = sqfl_chr2narestype(*(seq+sp));
  bascnt[rt]++;
  sp++;
  }
}

int main(int argc,
         char *argv[])
{
int ap;
char op;
int ecnts;
FILE *srcfl;
BC_CHR_BIN *chrbinlsts[ChrY];  /* bin lists for each chromosome */
int chrlens[ChrY];       /* length each chromosome */
char *genomsqhdrstr;         /* header string for genomic sequences */
int chrcnt;
BC_CHR_BIN *bp;
RBC_CHRNO chrno;
BC_CHR_BIN *chrbinends[ChrY];  /* bin ends for faster processing */
int sqpos;
int ubinlen;
int plustot;
int minustot;
int uchrno;    /* user has specified a chromosome, 0=>all */
char *fendp;
char *chromseq[ChrY];  /*seqs of each chromosome */
char *cp;
int cpgcnt;
int bstart;
BC_REGN_ELT *exclud[ChrY];
FILE *xcldfl;
FS_FSMSTRCT *cpgfsmp;    /* ptr to CpG searching fsm */
int minfrag;
int maxfrag;
int totlen;
int binend;
int prvbinend;
FILE *srcfl2;
int bincnt;
int cpgmin;
int amlgmate;
BC_CHR_BIN *chrbinsalt[ChrY];  /* alternate set of character bins */
int sumplustot;
int summinustot;
char *sqmonkhdrstr;
FILE *dfile;
int fcnt;
int frgpos;
int bascnt[RES_t+1];
SQ_RESTYPE resp;

debuglevel = RBC_dbg_none;
omode = BC_out_none;
srcstyle = RBC_src_rmapbs;
ecnts = 0;
srcfl = srcfl2 = NULL;
sqmonkhdrstr = genomsqhdrstr = NULL;
uchrno = 0;
ubinlen = BC_DEF_BINLEN;
cpgfsmp = fs_initnewfsm(4,1,FS_inv_ignor);
minfrag = maxfrag = cpgcnt = cpgmin = 0;
amlgmate = 0;
dfile = NULL;
for (chrno = Chr1; chrno <= ChrY; chrno++)
  exclud[chrno] = NULL;
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 'b':    /* bin length */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          ubinlen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (ubinlen <= 0)
              err_msg_die("Invalid bin length '%s'\n",argv[ap]);
          }
        break;
      case 'r':
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          if ((srcfl = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open .BED file %s\n",argv[ap]);
        break;
      case 'R':
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          if ((srcfl2 = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open .BED file %s\n",argv[ap]);
        break;
      case 'g':
        if (++ap > argc)
          err_msg_die("-%c needs header string\n",op);
        else
          genomsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'd':   /* debug on */
        debuglevel = RBC_dbg_on;
        break;
      case 'D':   /* debug^2 on */
        debuglevel = RBC_dbg_serious;
        break;
      case 'c':   /* a chromosome */
        if (++ap > argc)
          err_msg_die("-%c needs a chromosome identifier (1..20,X,Y)\n",op);
        else
          if ((uchrno = rbc_str2chrno(argv[ap])) == Chr_unk)
            err_msg_die("Can't determine Chromosome '%s'\n",argv[ap]);
        break;
      case 'l':         /* give listing */
        omode = BC_out_list;
        break;
      case 'L':         /* list non zero bins */
        omode = BC_out_listnz;
        break;
      case 'S':         /* generate SeqMonk .dat files */
        if (++ap > argc)
          err_msg_die("-%c needs dir header string for .dat files\n",op);
        else
          {
          sqmonkhdrstr = bas_strdup(argv[ap]);
          omode = BC_out_sqmonkdat;
          }
        break;
      case 'm':         /* diff methylated regions */
      case 'M':         /* make restricted representation bins */
      case 'N':         /* ditto, but list them */
      case 'n':         /* ditto, but include gaps */
      case 'k':         /* ditto, but show reads which don't map to RR bins */
      case 'K':         /* like -k but only show totals */
        switch (op)
          {
          case 'm':
            omode = BC_out_dmthrstrep;
            break;
          case 'N':
            omode = BC_out_rstrreplst;
            break;
          case 'k':
            omode = BC_out_rstrrepmiss;
            break;
          case 'K':
            omode = BC_out_rrmisstots;
            break;
          case 'M':
          default:
            omode = BC_out_rstrrepbins;
            break;
          }
        if (++ap < argc)
          if (*argv[ap] != '-')
            bc_commalst2ints(argv[ap],&minfrag,&maxfrag);
          else
            ap--;
        fs_adddatprs(cpgfsmp,"CCGG","CCGG");
        fs_adddatprs(cpgfsmp,"CG","CG");
        (void) fs_bldfsm(cpgfsmp,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
        break;
      case 'C':       /* -C set minimum CpG no for bins */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          {
          cpgmin = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (cpgmin < 0)
              err_msg_die("Invalid CpG min No.'%s'\n",argv[ap]);
          }
        break;
      case 'x':
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          if ((xcldfl = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open exlude file %s\n",argv[ap]);
          else
            (void) bc_getexcludregns(xcldfl,&exclud[0],ChrY);
        break;
      case 'A':
        amlgmate = 1;
        break;
      case 'h':
        say_usage(stdout,argv[0]);
        exit(0);
        break;
      default:
        err_msg("Unknown Option: '%s\n",argv[ap]);
        say_usage(stderr,argv[0]);
        exit(1);
        break;
      }
/* should be ready to go now */
for (chrno = Chr1; chrno <= ChrY; chrno++)
  {
  chrbinends[chrno-1] = chrbinlsts[chrno-1] = chrbinsalt[chrno-1] = NULL;
  chrlens[chrno-1] = 0;
  chromseq[chrno-1] = NULL;
  }
if ((genomsqhdrstr != NULL) && (ubinlen > 0))
  {
  switch (omode)
    {
    case BC_out_dmthrstrep:
      if ((srcfl == NULL) || (srcfl2 == NULL))
        err_msg_die("-m needs 2 read files (-r & -R)\n");
      else
        {
        chrcnt = bc_redrepgenomsqs(genomsqhdrstr,(int) ChrY,uchrno,cpgfsmp,
                                        &chrbinlsts[0],NULL,NULL);
        if (debuglevel > RBC_dbg_none)
          fprintf(stdout,"%d chromosomes read\n",chrcnt);
        if (amlgmate)
          (void) bc_amlgmatelocpgbins(&chrbinlsts[0],ChrY,minfrag,maxfrag,cpgmin);
        else
          (void) bc_losefailedbins(&chrbinlsts[0],ChrY,minfrag,maxfrag,cpgmin);
        for (chrno = Chr1; chrno <= ChrY; chrno++)
          if (chrbinlsts[chrno-1] != NULL)
            chrbinsalt[chrno-1] = bc_dupbinlst(chrbinlsts[chrno-1]);
        if ((ecnts = bc_chknreadsrcfl(&chrbinlsts[0],srcfl,uchrno,&exclud[0],ChrY)) <= 0)
          err_msg_die("Zero reads for -r source file\n");
        else
          {
          fclose(srcfl);
          if ((ecnts = bc_chknreadsrcfl(&chrbinsalt[0],srcfl2,uchrno,&exclud[0],ChrY)) <= 0)
            err_msg_die("Zero reads for -R source file\n");
          else
            {
            fclose(srcfl2);
            bc_scandiffmethbins(stdout,&chrbinlsts[0],&chrbinsalt[0],ChrY);
            }
          }
        }
      break;
    case BC_out_rstrreplst:
      chrcnt = bc_redrepgenomsqs(genomsqhdrstr,(int) ChrY,uchrno,cpgfsmp,
                                      &chrbinlsts[0],NULL,NULL);
      if (amlgmate)
        (void) bc_amlgmatelocpgbins(&chrbinlsts[0],ChrY,minfrag,maxfrag,cpgmin);
      else
        (void) bc_losefailedbins(&chrbinlsts[0],ChrY,minfrag,maxfrag,cpgmin);
      for (chrno = Chr1; chrno <= ChrY; chrno++)
        if ((bp = chrbinlsts[chrno-1]) != NULL)
          {
          prvbinend = 0;
          totlen = 0;
          while (bp != NULL)
            {
            binend = bp->spos + bp->binlen - 1;
            fprintf(stdout,"C%s: %d..%d, %d res %d cpG\n",rbc_chrno2str(chrno,1),
                      bp->spos,binend,bp->binlen,bp->cpgcnt);
            totlen += bp->binlen;
            bp = bp->nxtbin;
            }
          fprintf(stdout,"#C%s bin total=%d res in %d bins\n",rbc_chrno2str(chrno,1),totlen,
                    bc_cntcntbins(chrbinlsts[chrno-1]));
          }
      break;
    case BC_out_rstrrepbins:
    case BC_out_rstrrepmiss:
    case BC_out_rrmisstots:
    case BC_out_sqmonkdat:
      if (omode != BC_out_sqmonkdat)
        chrcnt = bc_redrepgenomsqs(genomsqhdrstr,(int) ChrY,uchrno,cpgfsmp,
                                     &chrbinlsts[0],NULL,NULL);
      else
        chrcnt = bc_redrepgenomsqs(genomsqhdrstr,(int) ChrY,uchrno,cpgfsmp,
                                     &chrbinlsts[0],&chromseq[0],&chrlens[0]);
      if (amlgmate)
        (void) bc_amlgmatelocpgbins(&chrbinlsts[0],ChrY,minfrag,maxfrag,cpgmin);
      else
        (void) bc_losefailedbins(&chrbinlsts[0],ChrY,minfrag,maxfrag,cpgmin);
      if ((omode == BC_out_rstrrepmiss) || (omode == BC_out_rrmisstots))
        for (chrno = Chr1; chrno <= ChrY; chrno++)
          {
          chrbinsalt[chrno-1] = bc_genlstforgaps(chrbinlsts[chrno-1]);
          rbc_clrallcntbins(&chrbinlsts[chrno-1]);
          chrbinlsts[chrno-1] = chrbinsalt[chrno-1];
          }
      switch (omode)
        {
        case BC_out_sqmonkdat:
          fcnt = 0;
          for (chrno = Chr1; chrno <= ChrY; chrno++)
            if ((bp = chrbinlsts[chrno-1]) != NULL)
              {
              if ((dfile = bc_opnsqmonkdatfl(sqmonkhdrstr,chrno,".dat")) == NULL)
                err_msg_die("Can't open SeqMonk .dat file for Chr%s\n",rbc_chrno2str(chrno,1));
              bc_mksqmnkhdr(dfile,"GRCh37",chrno,bc_sumbinlens(bp),minfrag,maxfrag);
              frgpos = 1;
              for (resp = RES_x; resp <= RES_t; resp++)
                bascnt[resp] = 0;
              while (bp != NULL)
                {
                fcnt++;
                fprintf(dfile,"FT   RR_FRAG         %d..%d\n",frgpos,
                          (frgpos + bp->binlen-1));
                fprintf(dfile,"FT                   /name=\"C%s_RRfrag%d\"\n",
                          rbc_chrno2str(chrno,1),fcnt);
                fprintf(dfile,
"FT                   /description=\"C%sFragment %d for %d..%d Red. Rep. genomic: %d..%d\"\n",
                          rbc_chrno2str(chrno,1),fcnt,minfrag,maxfrag,bp->spos,
                          (bp->spos+bp->binlen-1));
                frgpos += bp->binlen;
                bc_bascnt4regn(chromseq[chrno-1],chrlens[chrno-1],bp->spos,
                                 (bp->spos+bp->binlen-1),&bascnt[RES_x]);
                bp = bp->nxtbin;
                }
              fprintf(dfile,"SQ   Sequence %d BP;",frgpos-1);
              for (resp = RES_a; resp <= RES_t; resp++)
                fprintf(dfile,"  %d %c;",bascnt[resp],sqfl_restype2chr(resp));
              fprintf(dfile,"  %d other;\n//\n",bascnt[RES_x]);
              fclose(dfile);
              dfile = NULL;
              }
          break;
        case BC_out_rstrrepbins:
        case BC_out_rstrrepmiss:
        case BC_out_rrmisstots:
        default:
          if ((ecnts = bc_chknreadsrcfl(chrbinlsts,srcfl,uchrno,&exclud[0],ChrY)) > 0)
            {
            sumplustot = summinustot = 0;
            for (chrno = Chr1; chrno <= ChrY; chrno++)
             if ((bp = chrbinlsts[chrno-1]) != NULL)
                {
                totlen = bincnt = plustot = minustot = 0;
                while (bp != NULL)
                  {
                  switch (omode)
                    {
                    case BC_out_rstrrepmiss:
                      fprintf(stdout,"C%s: %d..%d (%d) %d meth, %d unmeth\n",rbc_chrno2str(chrno,1),
                                bp->spos,(bp->spos+bp->binlen-1),bp->binlen,bp->bcntfwd,bp->bcntrev);
                      break;
                    case BC_out_rrmisstots:
                      break;
                    default:
                      fprintf(stdout,"C%s: %d..%d (%d) %d meth, %d unmeth %d CpGs\n",rbc_chrno2str(chrno,1),
                                bp->spos,(bp->spos+bp->binlen-1),bp->binlen,bp->bcntfwd,bp->bcntrev,
                                bp->cpgcnt);
                      break;
                    }
                  plustot += bp->bcntfwd;
                  minustot += bp->bcntrev;
                  totlen += bp->binlen;
                  bincnt++;
                  bp = bp->nxtbin;
                  }
                fprintf(stdout,"#C%s totals: %d bins covering %d res +=%d -=%d\n",rbc_chrno2str(chrno,1),
                          bincnt,totlen,plustot,minustot);
                sumplustot += plustot;
                summinustot += minustot;
                }
            }
          if (omode == BC_out_rrmisstots)
            fprintf(stdout,"Totals missing reduced rep bins: +=%d -=%d both=%d\n",
                             sumplustot,summinustot,(sumplustot+summinustot));
          break;
        }
      break;
    case BC_out_list:
    case BC_out_listnz:
    case BC_out_none:
    default:
      chrcnt = rbc_getgenomesizs(genomsqhdrstr,&chrlens[0],(int) ChrY,uchrno);
      if (debuglevel > RBC_dbg_none)
        fprintf(stdout,"%d chromosomes read\n",chrcnt);
      for (chrno = Chr1; chrno <= ChrY; chrno++)
        if ((uchrno == 0) || (chrno == uchrno))
          {
          sqpos = 1;
          while (sqpos <= chrlens[chrno-1])
            {
            chrbinends[chrno-1] = rbc_appndbin(&chrbinends[chrno-1],sqpos,ubinlen);
            if (chrbinlsts[chrno-1] == NULL)
              chrbinlsts[chrno-1] = chrbinends[chrno-1];
            sqpos += ubinlen;
            }
          }
      switch (omode)
        { 
        case BC_out_list:
        case BC_out_listnz:
          if ((ecnts = bc_chknreadsrcfl(chrbinlsts,srcfl,uchrno,&exclud[0],ChrY)) > 0)
            {
            if (uchrno == 0)
              chrno = Chr1;
            else
              chrno = uchrno;
            while (chrno <= ChrY)
              {
              bp = chrbinlsts[chrno-1];
              while (bp != NULL)
                {
                switch (omode)
                  {
                  case BC_out_listnz:
                    if ((bp->bcntfwd > 0) || (bp->bcntrev > 0))
                      fprintf(stdout,"%s: %d..%d: +=%d -=%d\n",rbc_chrno2str(chrno,0),
                                bp->spos,(bp->spos+bp->binlen-1),bp->bcntfwd,bp->bcntrev);
                    break;
                  case BC_out_list:
                  default:
                    fprintf(stdout,"%s: %d..%d: +=%d -=%d\n",rbc_chrno2str(chrno,0),
                              bp->spos,(bp->spos+bp->binlen-1),bp->bcntfwd,bp->bcntrev);
                    break;
                  }
                bp = bp->nxtbin;
                }
              if (uchrno != 0)
                chrno = ChrY;
              chrno++;
              }
            }
          break;
        case BC_out_none:
        default:
          if (srcfl != NULL)
            {
            ecnts = bc_readsrcfl(chrbinlsts,srcfl,uchrno,&exclud[0],ChrY);
            if (debuglevel > RBC_dbg_none)
              {
              fprintf(stdout,"stored %d counts\n",ecnts);
              for (chrno = Chr1; chrno <= ChrY; chrno++)
                fprintf(stdout,"Chr%s: %dbp %d +hits %d -hits in %d bins\n",
                          rbc_chrno2str(chrno,1),chrlens[chrno-1],
                          bc_itsumcntfwdbins(chrbinlsts[chrno-1]),
                          bc_itsumcntrevbins(chrbinlsts[chrno-1]),
                          bc_cntcntbins(chrbinlsts[chrno-1]));
              if (ecnts > 0)
                {
                plustot = minustot = 0;
                for (chrno = Chr1; chrno <= ChrY; chrno++)
                  {
                  plustot += bc_itsumcntfwdbins(chrbinlsts[chrno-1]);
                  minustot += bc_itsumcntrevbins(chrbinlsts[chrno-1]);
                  }
                fprintf(stdout,"Totals: %d+ %d-\n",plustot,minustot);
                }
              }
            }
          break;
        }
      break;
    }
  }
exit(0);
}
