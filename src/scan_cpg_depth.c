/* scan_cpg_depth: take appropriately-formatted input files
(chrno, position, strand/status where strand/status = '+' or '-')
and count the number of +/- hits on each CpG for all or selected
chromosomes */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>
#include <sys/param.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "mrg_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"

#ifndef NO_ZLIB

#include <zlib.h>
/* #include <sys/types.h> */
#include "bam_fns.h"

#endif

#include "fsm_ops.h"

/* local defines */
/* #define PROG_VERSION 0.00 */
/* original */
/* #define PROG_VERSION 0.01 */
/* alternative restriction sites: Oct 2015 */
#define PROG_VERSION 0.02
/* work around '*' chrIDs in .bam files: Oct-2017 */

/* the length of sequence for each linked list cluster */
#define CHR_CLUSTER_DEF 1000

typedef struct SCD_cpgelt  /* element in linked list of CpGs */
  {
  int cpgpos;              /* position this CpG */
  int metcnt;              /* count methylated CpGs */
  int unmetcnt;            /* count unmethylated */
  struct SCD_cpgelt *nxtcpgelt;    /* forward link */
  struct SCD_cpgelt *prvcpgelt;    /* rev link */
  }
SCD_CPGELT;

typedef enum SCD_datum     /* different things we can return for SCD_CPGEELT list */
  {
  SCD_datm_count,
  SCD_datm_maxpos,
  SCD_datm_metcnt,
  SCD_datm_unmetcnt,
  SCD_datm_totcnt,
  SCD_datm_maxmetcnt,
  SCD_datm_maxunmetcnt,
  SCD_datm_maxtotcnt,
  SCD_datm_metunmetcpgpos,
  SCD_datm_metcpgpos,
  SCD_datm_unmetcpgpos
  }
SCD_DATUM;

typedef enum DM_read_sens
  {
  DM_rdsens_5p = 0,
  DM_rdsens_3p
  }
DM_READ_SENS;

typedef struct DM_chr_info     /* information about & bins for a chromosome */
  {
  char *chromseq;              /* chromosome seq pointer */
  int chrlen;
  char *chrflname;             /* file name for fa sequencs */
  char *chrid;                 /* id for this chr */
  }
DM_CHR_INFO;

typedef struct DM_sites
  {
  char *sitestr;
  int offset;
  struct DM_sites *nxtsite;
  struct DM_sites *prvsite;
  }
DM_SITES;

typedef enum SCD_srcmode
  {
  SCD_src_chrpos = 0,
  SCD_src_sam,
  SCD_src_bam
  }
SCD_SRCMODE;

typedef enum SCD_outstyle
  {
  SCD_out_stats = 0,
  SCD_out_list,
  SCD_out_histogram
  }
SCD_OUTSTYLE;

typedef struct SCD_runpars
  {
  int allowoutby1;
  FILE *srcfl;
  char *srcflnam;
  SCD_SRCMODE srcmode;
  int compress;
#ifndef NO_ZLIB
  BF_RUNDATA *bamrdp;   /* pointer to bam fns rundata */
#endif
  int outstyle;
  char *genomesqhdrstr;
  char *uchrno;       /* NULL => all */
  int clstrsiz;
  int omitzerocnt;
  FILE *hitfl;
  FILE *missfl;
  FS_FSMSTRCT *cpgfsmp;    /* ptr to CpG searching fsm */
  int maxchrno;
  int havexy;              /* 1 if mammalian type */
  char *chrinfoflnam;      /* name of chr info file */
  FILE *chrinfofl;         /* the chr info file handle */
  char *genomsqhdrstr;     /* hdr for numbered (+XY) chromosome files */
  char *destsqhdrstr;      /* hdr for output files */
  DM_CHR_INFO *chrinfo;    /* array of information on chromsomes */
  WRD_LUSTRCT *chridlu;    /* chromosome ID lookup */
  SCD_CPGELT **chrcpglsts;  /* cpg element lists */
  SCD_CPGELT ***clulists;   /* lookup lists for each chromosome */
  int sambuflen;
  }
SCD_RUNPARS;

#define DEF_SAMBUFLEN 150

typedef enum SAM_bsmk_fld
  {
  SAM_fld_header = 0,
  SAM_fld_sense,
  SAM_fld_chromosome,
  SAM_fld_positn,
  SAM_fld_flag,
  SAM_fld_cigar,
  SAM_fld_rnext,
  SAM_fld_pnext,
  SAM_fld_tlen,
  SAM_fld_seq,
  SAM_fld_qual,
  SAM_fld_nmtag,
  SAM_fld_xxtag,
  SAM_fld_methcall,
  SAM_fld_xrtag,
  SAM_fld_xgtag
  }
SAM_BSMK_FLD;

/* global debug & out style variables, for simplicity of access */
RBC_DBG debuglevel;
/* int allowoutby1;  */ /* set in order to accept out-by-one positions */

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
fprintf(fl,"%s v%.2f: scan for CpG read depth for all or selected chromosomes\n",pnam,
          PROG_VERSION);
fputs("Options:\n",fl);
fputs("     -r <posfile> read <posfile> as set of chr posit strand/meth\n",fl);
#ifdef NO_ZLIB
fputs(" or  -R <samfile> input from sam file\n",fl);
#else
fputs(" or  -R <sambamfile> input from sam (-Z) or bam (-z) file\n",fl);
#endif
fputs("     -g <genomehead> dir and file string to locate genomic seq files by adding n.fa\n",fl);
fputs(" or  -G <chr_info_file> file of chromosome IDs and filenames\n",fl);
fprintf(fl,"     -C <m> use cluster size of <m> for chromosome positions, Def=%d\n",CHR_CLUSTER_DEF);
#ifndef NO_ZLIB
fputs("     -z/-Z switch between sam (-Z) & bam (-z) for -R input - positional (def=-Z)\n",fl);
#endif
fputs("     -k <maxchrno> allow upto maxchrno different chromosomes (def = ChrY, 24)\n",fl);
fputs("     -Y don't expect XY chromosomes, just numbers (def=do)\n",fl);
fputs("     -l list each CpG to stdout with counts\n",fl);
fputs("     -c <n> restrict to Chromosome <n> 1..22,X,Y. Def=all\n",fl);
fprintf(fl,"     -b <sambuflen> max read length for .sam files (def=%d)\n",DEF_SAMBUFLEN);
fputs("     -p Permit out-by-one positions (e.g. Bismark complementary strand CpGs) def=don't\n",fl);
fputs("     -S generate statistics (range, mean, std deviation, etc. for counts\n",fl);
fputs("     -H generate histogram of counts\n",fl);
fputs("     -m list missed CpG lines\n",fl);
fputs("     -N list CpG hits (Nonmisses)\n",fl);
fputs("     -n omit zero count from histogram\n",fl);
}

SCD_CPGELT *scd_appndcpgelt(SCD_CPGELT **lstrt,
                            int cposn)
/* create and append a new element to *lstrt,
init counts to 0 & set position
Return address of new element */
{
SCD_CPGELT *prev, *end_ptr;

if (lstrt != NULL)
  {
  prev = end_ptr = *lstrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtcpgelt;
    }
  end_ptr = (SCD_CPGELT *) getmemory(sizeof(SCD_CPGELT),"CpG elt");
  end_ptr->nxtcpgelt = NULL;
  end_ptr->cpgpos = cposn;
  end_ptr->metcnt = end_ptr->unmetcnt = 0;
  if (*lstrt == NULL)
    {
    *lstrt = end_ptr;
    end_ptr->prvcpgelt = NULL;
    }
  else
    {
    prev->nxtcpgelt = end_ptr;
    end_ptr->prvcpgelt = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void scd_delcpgelt(SCD_CPGELT *ep,
                   SCD_CPGELT **lstrt)
/* delete ep from list *lstrt */
{
SCD_CPGELT *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvcpgelt) == NULL)
    *lstrt = ep->nxtcpgelt;
  else
    pt->nxtcpgelt = ep->nxtcpgelt;
  if ((pt = ep->nxtcpgelt) != NULL)
    pt->prvcpgelt = ep->prvcpgelt;
  memfree(ep);
  }
}

void scd_clrallcntbins(SCD_CPGELT **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  scd_delcpgelt(*lstrt,lstrt);
}

int scd_iterlstscan(SCD_CPGELT *clst,
                    SCD_DATUM datm)
/* iteratively traverse list clst, returning the appropriate
count/total */
{
int retval;
SCD_CPGELT *lp;

retval = 0;
lp = clst;
while (lp != NULL)
  {
  switch (datm)
    {
    case SCD_datm_count:
      retval++;
      break;
    case SCD_datm_maxpos:
      retval = imax(retval,lp->cpgpos);
      break;
    case SCD_datm_metcnt:
      retval += lp->metcnt;
      break;
    case SCD_datm_unmetcnt:
      retval += lp->unmetcnt;
      break;
    case SCD_datm_totcnt:
      retval += lp->metcnt + lp->unmetcnt;
      break;
    case SCD_datm_maxmetcnt:
      retval = imax(retval,lp->metcnt);
      break;
    case SCD_datm_maxunmetcnt:
      retval = imax(retval,lp->unmetcnt);
      break;
    case SCD_datm_maxtotcnt:
      retval = imax(retval,(lp->metcnt+lp->unmetcnt));
      break;
    case SCD_datm_metunmetcpgpos:
      retval++;
      break;
    case SCD_datm_metcpgpos:
      if (lp->metcnt > 0)
        retval++;
      break;
    case SCD_datm_unmetcpgpos:
      if (lp->unmetcnt > 0)
        retval++;
      break;
    default:
      break;
    }
  lp = lp->nxtcpgelt;
  }
return(retval);
}

int scd_datm4allchr(SCD_CPGELT *chrcpglsts[],
                    RBC_CHRNO mxchr,
                    SCD_DATUM datm)
/* sum datum test to all chromosomes to mxchr */
{
RBC_CHRNO chno;
int retval;

retval = 0;
for (chno = Chr1; chno <= mxchr; chno++)
  retval += scd_iterlstscan(chrcpglsts[chno-1],datm);
return(retval);
}

int scd_cntcpgelts(SCD_CPGELT *clst)
  /* recursively read list elements */
{
return(scd_iterlstscan(clst,SCD_datm_count));
}

SCD_CPGELT *scd_lastcpgelt(SCD_CPGELT *clst)
  /* iterate thru clst, returning
last element, NULL if none */
{
SCD_CPGELT *ep;

if ((ep = clst) == NULL)
  return(NULL);
else
  {
  while (ep->nxtcpgelt != NULL)
    ep = ep->nxtcpgelt;
  return(ep);
  }
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

int scd_pos2luindex(int clstrsiz,
                    int spos,
                    int widenby1)
/* return a lookup list index for spos */
{
if (clstrsiz != 0)
  {
  if (widenby1)
    spos = imax(1,spos-1);
  if ((spos % clstrsiz) == 0)
    return(imax(1,((int)( spos/clstrsiz) - 1)));
  else
    return((int) spos/clstrsiz);
  }
else
  return(0);
}

SCD_CPGELT **scd_mklookuplst(SCD_CPGELT *cpgelst,
                            int clstrsiz)
/* examine max position in cpgeltlst and
return a lookup list based on clstrsize to
find entry positions fast.  It is assumed here that
cpgelst is in ascending order of position */
{
int maxpos;
int ecnt;
SCD_CPGELT **lulst;
SCD_CPGELT *clp;
int luarrp;

if (clstrsiz == 0)
  return(NULL);
else
  {
  maxpos = scd_iterlstscan(cpgelst,SCD_datm_maxpos);
  ecnt = 1 + scd_pos2luindex(clstrsiz,maxpos,0);
  lulst = (SCD_CPGELT **) getmemory(ecnt*sizeof(SCD_CPGELT*),"List Lookup Array");
  for (luarrp = 0; luarrp < ecnt; luarrp++)
    *(lulst + luarrp) = NULL;
  clp = cpgelst;
  while (clp != NULL)
    {
    luarrp = scd_pos2luindex(clstrsiz,clp->cpgpos,0);
    if (*(lulst + luarrp) == NULL)    /* haven't seen this one yet */
      *(lulst + luarrp) = clp;
    clp = clp->nxtcpgelt;
    }
  return(lulst);
  }
}

SCD_CPGELT *scd_pos2cpgelt(SCD_RUNPARS *rpp,
                           int chrlen,
                           SCD_CPGELT **clulist,
                           int clstrsiz,
                           int posn)
/* cpgeltlst is an array of cpgelt pointers, one
for each clstrsize positions in chrlen.  look
up the appropriate cpg position and return a pointer to
the element for that position if it exists. NULL otherwise */
{
SCD_CPGELT *strtp;
int clstrno;

if ((posn <= chrlen) && (clstrsiz > 0))
  {
  strtp = *(clulist + (clstrno = scd_pos2luindex(clstrsiz,posn,rpp->allowoutby1)));
  while (strtp != NULL)
    if ((strtp->cpgpos == posn) ||
          (rpp->allowoutby1 && (abs(strtp->cpgpos-posn) == 1)))
      return(strtp);
    else
      if (strtp->cpgpos > posn+1) /* missed it, stop */
        return(NULL);
      else
        strtp = strtp->nxtcpgelt;
  return(NULL);
  }
else
  return(NULL);
}

int scd_readsrcfl(SCD_RUNPARS *rpp,
                  SCD_CPGELT **celsts,
                  DM_CHR_INFO *chrinfo,
                  SCD_CPGELT ***clulists,
                  int mxchr,
                  int clstrsiz,
                  FILE *sfl,
                  int *miscnt,
                  FILE *misfl,
                  FILE *hitfl)
/* use fscanf to read successive lines from sfl.
If rpp->uchrno is nonNULL, then only do that chromosome.

if miscnt is non-NULL then return to it the number of
missed counts */
{
char chrstr[5];
int sqpos;
char sensestr[5];
int matcnt;
int chrno;
int scnt;
SCD_CPGELT *bp;

matcnt = 0;
if (miscnt != NULL)
  *miscnt = 0;
while ((scnt = fscanf(sfl,"%s %d %s",&chrstr[0],&sqpos,&sensestr[0])) != EOF)
  if ((scnt == 3) && ((chrno = wlu_chkwrd(rpp->chridlu,&chrstr[0])) != Chr_unk) &&
        (chrno < mxchr) && ((rpp->uchrno == NULL) || (chrno == wlu_chkwrd(rpp->chridlu,rpp->uchrno))))
    {
    if ((bp = scd_pos2cpgelt(rpp,(chrinfo+chrno)->chrlen,*(clulists+chrno),clstrsiz,sqpos)) != NULL)
      {
      if (sensestr[0] == '+')
        bp->metcnt++;
      else
        bp->unmetcnt++;
      matcnt++;
      if (hitfl != NULL)  /* tell it about this hit */
        fprintf(hitfl,"Hit: C%s\t%d\t%c\n",rbc_chrno2str(chrno,1),sqpos,sensestr[0]);
      }
    else
      {
      if (miscnt != NULL)
        *miscnt += 1;
      if (misfl != NULL)
        fprintf(misfl,"Missed: C%s\t%d\t%c\n",rbc_chrno2str(chrno,1),sqpos,sensestr[0]);
      }
    }
return(matcnt);
}

int scd_chknreadsrcfl(SCD_RUNPARS *rpp,
                      SCD_CPGELT **celsts,
                      DM_CHR_INFO *chrinfo,
                      SCD_CPGELT ***clulists,
                      int mxchr,
                      int clstrsiz,
                      FILE *sfl,
                      int *miscnt,
                      FILE *misfl,
                      FILE *hitfl)
/* check sfl for openness, then call bc_readsrcfl, returning number
elements read */
{
if (sfl == NULL) /* can't do anything */
  return(0);
else
  return(scd_readsrcfl(rpp,celsts,chrinfo,clulists,mxchr,clstrsiz,sfl,miscnt,misfl,hitfl));
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

int scd_scangenomesqs(SCD_RUNPARS *rpp,
                      int chrmax,
                      FS_FSMSTRCT *cpgfsmp,
                      int clstrsiz,
                      SCD_CPGELT **chrcpglsts,
                      DM_CHR_INFO *chrinf)
/* using previously created chromosomal file names, 
open each as a Fasta sequence
file and create a buffer for each.  Return the
number of chromosomes processed.  Generate cpglsts for each
chromosome scanned, adding observed CpG positions thereto */
{
int chno;
int nblen;
SQFL_STRCT *chsqfl;
int rcnt;
char nxtres;
int cpos;
FS_RESELT *frp;
FS_DATPRELT *dpep;
int clstlen;
int clstp;
SCD_CPGELT *clstep;
int hitpos;
int clstrno;
int slen;
int hcnt;
char *sqp;

rcnt = 0;
for (chno = 0; chno < chrmax; chno++)
  if ((rpp->uchrno == NULL) || (chno == wlu_chkwrd(rpp->chridlu,rpp->uchrno)))
    {
    if ((chsqfl = sqfl_opnsqstrct((chrinf+chno)->chrflname,SFMT_fasta,"r")) != NULL)
      {
      (chrinf+chno)->chrlen = slen = readsrcsq(chsqfl,NULL);
      clstlen = (int) slen/clstrsiz + 1;
      if (rpp->srcmode > SCD_src_chrpos)
        sqp = (chrinf+chno)->chromseq =
               (char *) getmemory(sizeof(char) *((chrinf+chno)->chrlen+1),"chromosomebuff");
      sqfl_rewind(chsqfl);
      (void) sqfl_skipsqflhdr(chsqfl);
      fs_initrun(cpgfsmp);
      cpos = 0;
      hcnt = 0;
      clstep = NULL;
      while ((nxtres = sqfl_getnxtres(chsqfl)) != '\0')
        {
        cpos++;
        if ((frp = fs_procchr(cpgfsmp,nxtres,tr_nares2int)) != NULL)
          {
          dpep = (FS_DATPRELT *) frp->action;
          hitpos = cpos - dpep->ldstr + 1;
          hcnt++;
          clstep = scd_appndcpgelt(&clstep,hitpos);
          if (*(chrcpglsts + chno) == NULL)
            *(chrcpglsts + chno) = clstep;
/*           if ((rpp->chrcpglsts + chno) == NULL)
                 *(rpp->chrcpglsts + chno) = clstep; */
          }
        if (rpp->srcmode > SCD_src_chrpos)
          {
          *sqp = nxtres;
          sqp++;
          *sqp = '\0';
          }
        }
      rcnt++;
      }
    else
      err_msg("Can't open chromosome file %s\n",(chrinf+chno)->chrflname);
    }
return(rcnt);
}

double scd_medianxcntsarr(int *cnts,
                          int maxcnt)
/* cnts is an array of integers 0..maxcnt.  total
them and work out the median. */
{
int ttl;
int cp;
int medcnt;
int medcp1;
int tcnt;
int nxtnz;

cp = 0;
ttl = 0;
while (cp <= maxcnt)
  {
  ttl += *(cnts+cp);
  cp++;
  }
medcp1 = 0;
medcnt = (int) ttl/2.0;
if ((ttl % 2) == 0)
  medcp1 = medcnt + 1;
else
  medcnt++;
tcnt = 0;
cp = 0;
while (cp <= maxcnt)
  {
  if (tcnt >= medcnt)
    {
    cp--;
    if ((medcp1 == 0) /* odd */ ||
          (tcnt != medcnt)) /* even but don't need to look for next cnt bucket */
      return((double) cp);
    else             /* even, but need to look for next filled count bucket */
      {
      nxtnz = cp + 1;
      while ((nxtnz <= maxcnt) && (*(cnts+nxtnz) == 0))
        nxtnz++;
      return((cp + imin(nxtnz,maxcnt))/2.0);
      }
    }
  tcnt += *(cnts+cp);
  cp++;
  }
return((double) maxcnt);
}

int scd_modexcntsarr(int *cnts,
                     int maxcnt)
/* cnts is an array of integers 0..maxcnt. scan
for largest no and return position as the mode */
{
int biggest;
int cp;
int bigp;

biggest = 0;
cp = 0;
bigp = 0;
while (cp <= maxcnt)
  {
  if (*(cnts+cp) > biggest)
    {
    biggest = *(cnts+cp);
    bigp = cp;
    }
  cp++;
  }
return(bigp);
}

int dm_scanchrinfofl(FILE *cifile,
                     DM_CHR_INFO *cinfo,
                     WRD_LUSTRCT *chridlu)
/* scan lines in cifile as chromosomeID Filenam.
return number.  If cinfo is non-NULL, then put
values into fields. '#' causes remainder of line
to be treated as comment (unless quoted). */
{
int lno;
char *cid;
char *chrfile;
char nc;
char *bp;
int onquote;
char *target;
char prvchr;
int incmt;

cid = (char *) getmemory(MAXPATHLEN + 1,"chrid");
chrfile = (char *) getmemory(MAXPATHLEN + 1,"chrflnam");
target = bp = cid;
*chrfile = *bp = '\0';
lno = 0;
incmt = onquote = 0;
prvchr = '\0';
while ((nc = fgetc(cifile)) != EOF)
  switch (nc)
    {
    case '\n':
      if ((strlen(cid) > 0) && (strlen(chrfile) > 0))
        {
        if (cinfo != NULL)
          {
          (cinfo+lno)->chrid = bas_strdup(cid);
          (cinfo+lno)->chrflname = bas_strdup(chrfile);
          }
        if (chridlu != NULL)
          wlu_addwrd(chridlu,cid,(int) lno, NULL);
        *chrfile = *cid = '\0';
        lno++;
        }
/* reset to start */
      target = bp = cid;
      incmt = onquote = 0;
      prvchr = '\0';
      break;
    case '"':    /* start or end of quoted string */
      if (!incmt)
        onquote = !onquote;
      break;
    case ' ':     /* white space */
    case '\t':
      if (onquote)
        {
        if ((target == cid) && (!incmt))
          bas_appchr(cid,&bp,nc,MAXPATHLEN);
        else
          if ((target == chrfile) && (!incmt))
            bas_appchr(chrfile,&bp,nc,MAXPATHLEN);
        }
      else 
        {
        if (!incmt)
          {
          if (!isblank(prvchr))
            if (target == cid)
              target = bp = chrfile;
            else
              target = NULL;
          prvchr = nc;
          }
        }
      break;
    case '#':    /* start of comment field */
      if (!onquote)
        {
        incmt = 1;
        target = NULL;
        }
    default:
      if ((target == cid) && !incmt)
        bas_appchr(cid,&bp,nc,MAXPATHLEN);
      else
        if ((target == chrfile) && !incmt)
          bas_appchr(chrfile,&bp,nc,MAXPATHLEN);
      prvchr = nc;
      break;
    }
memfree(chrfile);
memfree(cid);
return(lno);
}      

void dm_bldchrfilenames(char *hdrstr,
                        SCD_RUNPARS *rpp)
/* create set of file names based on hdrstr.
put into rpp->chrinfo array (pre-existing) */
{
int cno;
char *sqfnam;
int nblen;

sqfnam = (char *) getmemory((nblen = (int) strlen(hdrstr) + 16),"Sq filename buff");
for (cno = 0; cno < rpp->maxchrno; cno++)
  {
  (rpp->chrinfo+cno)->chrid = bas_strdup(any_chrno2str(rpp->maxchrno,rpp->havexy,(cno+1),1));
  wlu_addwrd(rpp->chridlu,(rpp->chrinfo+cno)->chrid,(int) cno,NULL);
  snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,(rpp->chrinfo+cno)->chrid);
  (rpp->chrinfo+cno)->chrflname = bas_strdup(sqfnam);
  }
memfree(sqfnam);
}

int dm_notemethstat4chrnpos(SCD_RUNPARS *rpp,
                            int chrno,
                            int sqpos,
                            int methstat,
                            int fst3pcpg)
/* the guts of noting methylation status.  Note that
chrno is expected from 0 */
{
SCD_CPGELT *cep;

if ((cep = scd_pos2cpgelt(rpp,(rpp->chrinfo+chrno)->chrlen,
                            *(rpp->clulists+chrno),rpp->clstrsiz,sqpos)) != NULL)
  {
  if (methstat)
    cep->metcnt++;
  else
    cep->unmetcnt++;
  return(1);
  }
else
  return(0);
}

int dm_cmpsamread(SCD_RUNPARS *rpp,
                  char *sqread,
                  int readlen,
                  int rdpos,
                  int chrno,
                  DM_READ_SENS rdsens,
                  int (*proccpgpos)(SCD_RUNPARS *xrpp,
                                    int xchrno,
                                    int xsqpos,
                                    int xmetstat,
                                    int fst3pcpg))
/* compare this sequence read with rdsens, return
cpg count.  chrno now 0 based */
{
char *cspt;
char *rdpt;
int cpgcnt;
int rdlen;
int cpgpos;
int climit;

if (((rpp->uchrno == NULL) || (chrno == wlu_chkwrd(rpp->chridlu,rpp->uchrno))) && (chrno >= 0))
  {  /* reject MT */
  rdlen = readlen;
  if (rdsens == DM_rdsens_5p)
    {
    cpgpos = rdpos;
    rdpt = sqread;
    cspt = ((rpp->chrinfo+chrno)->chromseq+rdpos-1);
    climit = (rpp->chrinfo+chrno)->chrlen;
    while ((rdlen > 0) && rbc_intinrng(1,cpgpos,climit))
      {
      if ((toupper(*cspt) == 'C') && (toupper(*(cspt+1)) =='G'))
        switch (toupper(*rdpt))
          {
          case 'T':      /* unmeth */
            (void)(* proccpgpos)(rpp,chrno,cpgpos,0,0);
            break;
          case 'C':      /* meth */
            (void)(* proccpgpos)(rpp,chrno,cpgpos,1,0);
            break;
          default:
            break;
          }
      cpgpos++;
      rdpt++;
      cspt++;
      rdlen--;
      }
    }
  else
    {
    cpgcnt = 0;
    rdpt = sqread + rdlen - 1;
    cspt = ((rpp->chrinfo+chrno)->chromseq+rdpos+rdlen-2);
    climit = (rpp->chrinfo+chrno)->chrlen;
    cpgpos = rdpos + rdlen - 2;
    while ((rdpt >= sqread) && rbc_intinrng(1,cpgpos,climit))
      {  /* working in base complement space this time */
      if ((toupper(*cspt) == 'G') && (toupper(*(cspt-1)) =='C'))
        {
        switch (toupper(*rdpt))
          {
          case 'A':   /* unmeth */
            (void)(* proccpgpos)(rpp,chrno,cpgpos,0,cpgcnt==0);
            break;
          case 'G':
            (void)(* proccpgpos)(rpp,chrno,cpgpos,1,cpgcnt==0);
            break;
          default:
            break;
          }
        cpgcnt++;
        }
      cpgpos--;
      rdpt--;
      cspt--;
      }
    }
  return(1);
  }
else
  return(0);
}

int dm_readsamfl(FILE *sfl,
                 SCD_RUNPARS *rpp)
/* sfl is a bismark .sam file: use fscanf to read 
successive lines from it. */
{
char *readbuf;
SAM_BSMK_FLD atmno;
char nc;
char *bp;
DM_READ_SENS rdsens;
int sqpos;
int chrno;
int matcnt;
int lno;

matcnt = 0;
readbuf = (char *) getmemory(rpp->sambuflen+1,"Sam Read buf");
atmno = SAM_fld_header;
bp = readbuf;
/* rpp->srclno = 0; */
lno = 0;
while ((nc = fgetc(sfl)) != EOF)
  switch (nc)
    {
    case '\n':    /* end of line, just reset things */
      lno++;
/*      rpp->srclno++; */
      matcnt += dm_cmpsamread(rpp,readbuf,strlen(readbuf),sqpos,chrno,rdsens,
                                dm_notemethstat4chrnpos);
      atmno = SAM_fld_header;
      bp = readbuf;
      break;
    case '\t':    /* end of atom: decide what we should do next */
      switch (atmno)
        {
        case SAM_fld_sense:
          switch ((int)strtol(readbuf,NULL,10))
            {
            case 16:
              rdsens = DM_rdsens_3p;
              break;
            case 0:
            default:
              rdsens = DM_rdsens_5p;
              break;
            }
          bp = readbuf;
          atmno++;
          break;
        case SAM_fld_chromosome:
          chrno = wlu_chkwrd(rpp->chridlu,readbuf);
/*          chrno = any_str2chrno(rpp->maxchrno,rpp->havexy,readbuf); */
          bp = readbuf;
          atmno++;
          break;
        case SAM_fld_positn:
          sqpos = (int) strtol(readbuf,NULL,10);
          bp = readbuf;
          atmno++;
          break;
        case SAM_fld_header:
          bp = readbuf;
          atmno++;
          break;
        case SAM_fld_seq:     /* got what we want for now, process & bail out this line */
          lno++;
/*          rpp->srclno++; */
          matcnt += dm_cmpsamread(rpp,readbuf,strlen(readbuf),sqpos,chrno,rdsens,
                                    dm_notemethstat4chrnpos);
          bas_skipeol(sfl,NULL);
          atmno = SAM_fld_header;
          bp = readbuf;
          break;
        default:
          bp = readbuf;
          atmno++;
          break;
        }
      break;
    case '@':
      if ((atmno == SAM_fld_header) && (bp == readbuf)) /* still on start of line, bail out */
        {
        lno++;
/*        rpp->srclno++; */
        bas_skipeol(sfl,NULL);
        bp = readbuf;
        break;
        }
    default:
      bas_appchr(readbuf,&bp,nc,rpp->sambuflen);
      break;
    }
memfree(readbuf);
return(matcnt);
}

#ifndef NO_ZLIB

int dm_readbamfl(SCD_RUNPARS *rpp)
/* read from previously opened .bam file at rpp->bamrpd.
read successive successive records, return number of matches
*/
{
DM_READ_SENS rdsens;
int chrno;
int matcnt;
BF_EXPNDALIGNREC *earp;

matcnt = 0;
while ((earp = bf_nxtalignrec(rpp->bamrdp,1,err_msg_die)) != NULL)
  {
  if (earp->arec->refID >=0) /* avoid '*' chrIDs */
    {
    chrno = wlu_chkwrd(rpp->chridlu,(rpp->bamrdp->refinfo+earp->arec->refID)->refname);
    switch (earp->arec->FLAG)
      {
      case 16:
        rdsens = DM_rdsens_3p;
        break;
      case 0:
      default:
        rdsens = DM_rdsens_5p;
        break;
      }
    matcnt += dm_cmpsamread(rpp,earp->seq,earp->arec->l_seq,(earp->arec->pos+1),
                              chrno,rdsens,dm_notemethstat4chrnpos);
    }
  bf_disposeexpalignrec(earp);
  }
return(matcnt);
}

#endif

void scd_tellchrinfo(FILE *ofl,
                     SCD_RUNPARS *rpp)
/* tell ofl about the chromosome information
for this run */
{
if (rpp->genomesqhdrstr != NULL)
  if (rpp->uchrno != NULL)
    fprintf(ofl,"%d Chromosome:%s%s.fa\n",rpp->maxchrno,rpp->genomesqhdrstr,
              rpp->uchrno);
  else
    fprintf(ofl,"%d Chromosome%s:%s<n>.fa\n",rpp->maxchrno,(rpp->maxchrno>1?"s":""),
              rpp->genomesqhdrstr);
else
  if (rpp->chrinfoflnam != NULL)
    fprintf(ofl,"%d Chromosome%s from %s\n",rpp->maxchrno,(rpp->maxchrno>1?"s":""),
              rpp->chrinfoflnam);
  else
    fputs("Chromosome info undefined\n",ofl);
}

void scd_tellruninfo(FILE *ofl,
                     SCD_RUNPARS *rpp,
                     int matcnts,
                     int missedincnts)
/* same data for all outputs... */
{
fprintf(ofl,"Source:%s\n",rpp->srcflnam);
scd_tellchrinfo(ofl,rpp);
fprintf(ofl,"hit counts =%d, missed counts=%d, OutByOne=%s\n",
          matcnts,missedincnts,(rpp->allowoutby1?"yes":"no"));
}

int main(int argc,
         char *argv[])
{
int ap;
char op;
int ecnts;
FILE *srcfl;
/* SCD_CPGELT *chrcpglsts[ChrY]; */ /* bin lists for each chromosome */
/*int chrlens[ChrY];   */    /* length each chromosome */
/* char *genomsqhdrstr; */        /* header string for genomic sequences */
int chrcnt;
RBC_CHRNO chrno;
/* int uchrno; */   /* user has specified a chromosome, 0=>all */
int cpgcnt;
/* FS_FSMSTRCT *cpgfsmp; */  /* ptr to CpG searching fsm */
int listoutput;
int statsoutput;
SCD_CPGELT *cp;
int clstrno;
int clstp;
/* SCD_CPGELT **clulists[ChrY]; */
int maxcount;
double sumx;
double sumx2;
int *distn;
int dptr;
int curcnt;
int mincnt;
int histgrmoutput;
int biggstcnt;
int nast;
int ndigs;
int zcnt;
int hdigs;
int missedincnts;
int matcnts;
SCD_RUNPARS rpars;
int chrp;
char *fendp;

debuglevel = RBC_dbg_none;
rpars.allowoutby1 = 0;      /* don't by default */
ecnts = 0;
rpars.srcfl = NULL;
rpars.srcmode = SCD_src_chrpos;
rpars.compress = 0;
#ifndef NO_ZLIB
rpars.bamrdp = NULL;
#endif
rpars.genomsqhdrstr = NULL;
rpars.uchrno = NULL;
rpars.clstrsiz = CHR_CLUSTER_DEF;
rpars.cpgfsmp = fs_initnewfsm(4,1,FS_inv_ignor);
rpars.srcflnam = "";
rpars.omitzerocnt = 0;
rpars.hitfl;
rpars.missfl = NULL;
rpars.maxchrno = ChrY;
rpars.havexy = 1;
rpars.outstyle = 0;
rpars.sambuflen = DEF_SAMBUFLEN;
rpars.chrcpglsts = NULL;
rpars.clulists = NULL;
/* listoutput = statsoutput = histgrmoutput = 0; */
missedincnts = 0;
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 'R':
        rpars.allowoutby1 = 0;
        if (rpars.compress)
          rpars.srcmode = SCD_src_bam;
        else
          rpars.srcmode = SCD_src_sam;
      case 'r':
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          rpars.srcflnam = bas_strdup(argv[ap]);
        break;
      case 'g':
        if (++ap > argc)
          err_msg_die("-%c needs header string\n",op);
        else
          rpars.genomsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'G':   /* file name for chromosome spec file */
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          {
          rpars.chrinfoflnam = bas_strdup(argv[ap]);
          if ((rpars.chrinfofl = fopen(rpars.chrinfoflnam,"r")) == NULL)
            err_msg_die("Can't open chromosome information file '%s' for -%c\n",
                          rpars.chrinfoflnam,op);
          }
        break;
      case 'd':   /* debug on */
        debuglevel = RBC_dbg_on;
        break;
      case 'D':   /* debug^2 on */
        debuglevel = RBC_dbg_serious;
        break;
      case 'c':   /* a chromosome */
        if (++ap > argc)
          err_msg_die("-%c needs a chromosome identifier (%s..%s,%s,%s)\n",op,
                        any_chrno2str(rpars.maxchrno,rpars.havexy,1,1),
                        any_chrno2str(rpars.maxchrno,rpars.havexy,rpars.maxchrno-2,1),
                        any_chrno2str(rpars.maxchrno,rpars.havexy,rpars.maxchrno-1,1),
                        any_chrno2str(rpars.maxchrno,rpars.havexy,rpars.maxchrno,1));
        else
          rpars.uchrno = bas_strdup(argv[ap]);
        break;
      case 'l':         /* give listing */
        rpars.outstyle = rpars.outstyle | 1 << SCD_out_list;
        break;
      case 'S':         /* statistics output */
        rpars.outstyle = rpars.outstyle | 1 << SCD_out_stats;
        break;
      case 'H':         /* histogram */
        rpars.outstyle = rpars.outstyle | 1 << SCD_out_histogram;
        break;
      case 'm':         /* write misses to stdout */
        rpars.missfl = stdout;
        break;
      case 'N':        /* write nonmisses to stdout */
        rpars.hitfl = stdout;
        break;
      case 'p':        /* permit out-by-one position errors */
        rpars.allowoutby1 = 1;
        break;
      case 'n':        /* omit zero count bin from histogram */
        rpars.omitzerocnt = 1;
        break;
      case 'Y':
        rpars.havexy = 0;
        break;
      case 'k':       /* different chromosome Number */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          {
          rpars.maxchrno = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (rpars.maxchrno < 0)
              err_msg_die("Invalid chromosome number '%s'\n",argv[ap]);
          }
        break;
      case 'z':     /* bam format */
#ifdef NO_ZLIB
        err_msg_die("Compiled with NO_ZLIB set: can't read bam files\n");
#else
        rpars.compress = 1;
#endif
        break;
      case 'Z':     /* sam format */
        rpars.compress = 0;
        break;
      case 'b':      /* sam buffer length */
        if (++ap > argc)
          err_msg_die("-%c needs bufferlength\n",op);
        else
          {
          rpars.sambuflen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert int '%s' for -%c\n",argv[ap],op);
          }
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
if (rpars.chrinfofl != NULL)
  rpars.maxchrno = dm_scanchrinfofl(rpars.chrinfofl,NULL,NULL);
/* create chromosome information array */
rpars.chrinfo = (DM_CHR_INFO *) getmemory(sizeof(DM_CHR_INFO)*rpars.maxchrno,
                                         "Chr bin info");
/* chr ID lookup stuff */
rpars.chridlu = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"ChrIDlookup");
wlu_initlustrct(rpars.chridlu,WLU_CASEIND,-1);
rpars.chrcpglsts = (SCD_CPGELT **) getmemory(sizeof(SCD_CPGELT *)*rpars.maxchrno,
                                             "CpGlists");
rpars.clulists = (SCD_CPGELT ***) getmemory(sizeof(SCD_CPGELT **)*rpars.maxchrno,
                                             "CpGlistslookup");
/* init info array */
for (chrp = 0; chrp < rpars.maxchrno; chrp++)
  {
  (rpars.chrinfo+chrp)->chromseq = NULL;
  (rpars.chrinfo+chrp)->chrlen = 0;
  (rpars.chrinfo+chrp)->chrflname = NULL;
  (rpars.chrinfo+chrp)->chrid = NULL;
  *(rpars.chrcpglsts+chrp) = NULL;
  *(rpars.clulists+chrp) = NULL;
  }
/* establish file names, etc. either from params or the info file */
if ((rpars.chrinfofl == NULL) && (rpars.genomsqhdrstr != NULL))
  dm_bldchrfilenames(rpars.genomsqhdrstr,&rpars);
else
  {
  rewind(rpars.chrinfofl);
  (void) dm_scanchrinfofl(rpars.chrinfofl,rpars.chrinfo,rpars.chridlu);
  fclose(rpars.chrinfofl);
  }
if (debuglevel > RBC_dbg_none)
  fprintf(stdout,"info for %d chromosomes found\n",rpars.maxchrno);
if (((rpars.genomsqhdrstr != NULL) || (rpars.chrinfofl != NULL)) &&
       (rpars.outstyle) && (rpars.clstrsiz > 0))
  {
  fs_adddatprs(rpars.cpgfsmp,"CG","CG");
  (void) fs_bldfsm(rpars.cpgfsmp,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
  (void) scd_scangenomesqs(&rpars,rpars.maxchrno,rpars.cpgfsmp,
                             rpars.clstrsiz,rpars.chrcpglsts,rpars.chrinfo);
  if (debuglevel > RBC_dbg_none)
    for (chrno = 0; chrno < rpars.maxchrno; chrno++)
      if ((rpars.chrcpglsts + chrno) != NULL)
        fprintf(stdout,"C%s %d elts\n",*(rpars.chrinfo + chrno)->chrid,
                  scd_iterlstscan(*(rpars.chrcpglsts + chrno),SCD_datm_count));
  for (chrno = 0; chrno < rpars.maxchrno; chrno++)
    if ((rpars.chrcpglsts + chrno) != NULL)
      *(rpars.clulists + chrno) = scd_mklookuplst(*(rpars.chrcpglsts + chrno),rpars.clstrsiz);
  }
/* check source file */
if (rpars.srcflnam== NULL)
  err_msg_die("no source file name\n");
else
  {
  switch (rpars.srcmode)
    {
    case SCD_src_bam:
#ifdef NO_ZLIB
      err_msg_die("BAM source unavailable (zlib problem)\n");
#else
      if ((rpars.bamrdp = bf_opnchkbamflmsg(rpars.srcflnam,err_msg_die)) == NULL)
        exit(1);
      else
        {
        matcnts = dm_readbamfl(&rpars);
        bf_disposerundata(rpars.bamrdp);
        rpars.bamrdp = NULL;
        }
#endif
      break;
    case SCD_src_sam:
    case SCD_src_chrpos:
    default:
      switch (rpars.srcmode)
        {
#ifndef NO_ZLIB
        case SCD_src_bam:
          if ((rpars.bamrdp = bf_opnchkbamflmsg(rpars.srcflnam,err_msg_die)) == NULL)
            err_msg_die("Can't open source file '%s'\n",rpars.srcflnam);
          else
            matcnts = dm_readbamfl(&rpars);
          break;
#endif
        case SCD_src_sam:
          if ((rpars.srcfl = fopen(rpars.srcflnam,"r")) == NULL)
            err_msg_die("Can't open source file '%s'\n",rpars.srcflnam);
          else
            matcnts = dm_readsamfl(rpars.srcfl,&rpars);
          break;
        case SCD_src_chrpos:
        default:
          if ((rpars.srcfl = fopen(rpars.srcflnam,"r")) == NULL)
            err_msg_die("Can't open source file '%s'\n",rpars.srcflnam);
          else
            matcnts = scd_chknreadsrcfl(&rpars,rpars.chrcpglsts,rpars.chrinfo,rpars.clulists,rpars.maxchrno,
                                        rpars.clstrsiz,rpars.srcfl,&missedincnts,
                                        rpars.missfl,rpars.hitfl);
          break;
        }
    }
  if (matcnts <= 0)
    err_msg_die("Run failed at sourcefile ('%s') read\n",rpars.srcflnam);
  if (rpars.outstyle & 1 << SCD_out_list)
    {
    scd_tellruninfo(stdout,&rpars,matcnts,missedincnts);
    for (chrno = 0; chrno < rpars.maxchrno; chrno++)
      {
      cp = *(rpars.chrcpglsts + chrno);
      while (cp != NULL)
        {
        fprintf(stdout,"C%s\t%d\t+=%d\t-=%d\tsum=%d\n",(rpars.chrinfo + chrno)->chrid,
                  cp->cpgpos,cp->metcnt,cp->unmetcnt,(cp->metcnt+cp->unmetcnt));
        cp = cp->nxtcpgelt;
        }
      }
    }
  if (rpars.outstyle & 1 << SCD_out_stats)
    {
    maxcount = 0;
    scd_tellruninfo(stdout,&rpars,matcnts,missedincnts);
    fprintf(stdout,"TotMethCpG=%d at %dMCpGs, TotUnMethCpG=%d at %dUnMCpGs, (CpGs Meth+UnMeth=%d)\n",
              scd_datm4allchr(rpars.chrcpglsts,rpars.maxchrno,SCD_datm_metcnt),
              scd_datm4allchr(rpars.chrcpglsts,rpars.maxchrno,SCD_datm_metcpgpos),
              scd_datm4allchr(rpars.chrcpglsts,rpars.maxchrno,SCD_datm_unmetcnt),
              scd_datm4allchr(rpars.chrcpglsts,rpars.maxchrno,SCD_datm_unmetcpgpos),
              scd_datm4allchr(rpars.chrcpglsts,rpars.maxchrno,SCD_datm_metunmetcpgpos));
    for (chrno = 0; chrno < rpars.maxchrno; chrno++)
      maxcount = imax(maxcount,scd_iterlstscan(*(rpars.chrcpglsts + chrno),SCD_datm_maxtotcnt));
    distn = (int *) getmemory((maxcount + 1) * sizeof(int),"Distribution");
    sumx = sumx2 = 0.0;
    cpgcnt = 0;
    mincnt = maxcount;
    for (dptr = 0; dptr <= maxcount; dptr++)
      *(distn + dptr) = 0;
    for (chrno = 0; chrno < rpars.maxchrno; chrno++)
      {
      cp = *(rpars.chrcpglsts + chrno);
      while (cp != NULL)
        {
        cpgcnt++;
        curcnt = cp->metcnt + cp->unmetcnt;
        sumx += (double) curcnt;
        sumx2 += (double) curcnt * curcnt;
        if (curcnt <= maxcount)
          (*(distn+curcnt))++;
        mincnt = imin(mincnt,curcnt);
        cp = cp->nxtcpgelt;
        }
      }
    fprintf(stdout,"Counts: min=%d max=%d mean=%.2f sdev=%.2f\n",
              mincnt,maxcount,(double)sumx/cpgcnt,sqrt((sumx2 - sumx*sumx/cpgcnt)/(cpgcnt-1)));
    fprintf(stdout,"Median = %.1f mode=%d\n",scd_medianxcntsarr(distn,maxcount),
              scd_modexcntsarr(distn,maxcount));
    }
  if (rpars.outstyle & 1 << SCD_out_histogram)
    {
    scd_tellruninfo(stdout,&rpars,matcnts,missedincnts);
    maxcount = 0;
    for (chrno = 0; chrno < rpars.maxchrno; chrno++)
      maxcount = imax(maxcount,scd_iterlstscan(*(rpars.chrcpglsts + chrno),SCD_datm_maxtotcnt));
    distn = (int *) getmemory((maxcount + 1) * sizeof(int),"Distribution");
    cpgcnt = 0;
    mincnt = maxcount;
    for (dptr = 0; dptr <= maxcount; dptr++)
      *(distn + dptr) = 0;
    for (chrno = 0; chrno < rpars.maxchrno; chrno++)
      {
      cp = *(rpars.chrcpglsts + chrno);
      while (cp != NULL)
        {
        cpgcnt++;
        curcnt = cp->metcnt + cp->unmetcnt;
        if (curcnt <= maxcount)
          (*(distn+curcnt))++;
        mincnt = imin(mincnt,curcnt);
        cp = cp->nxtcpgelt;
        }
      }
    biggstcnt = 0;
    if (rpars.omitzerocnt)
      dptr = 1;
    else
      dptr = 0;
    while (dptr <= maxcount)
      {
      biggstcnt = imax(*(distn+dptr),biggstcnt);
      dptr++;
      }
    if (rpars.omitzerocnt)
      dptr = 1;
    else
      dptr = 0;
    ndigs = bas_digitsin(biggstcnt);
    hdigs = bas_digitsin(maxcount);
    zcnt = 0;
    while (dptr <= maxcount)
      {
      if (*(distn+dptr) > 0)
        zcnt = 0;
      else
        zcnt++;
      if (zcnt <= 3)
        {
        fprintf(stdout,"%*d %*d |",hdigs,dptr,ndigs,*(distn+dptr));
        if (rpars.omitzerocnt && (dptr == 0))
          nast = 50 - bas_digitsin(*(distn)) - 1;
        else
          nast = (int)((*(distn+dptr)*50)/ biggstcnt);
        while (nast > 0)
          {
          fputc('*',stdout);
          nast--;
          }
        if (rpars.omitzerocnt & (dptr == 0))
          fprintf(stdout,">%d",*(distn));
        fputc('\n',stdout);
        }
      else
        if (zcnt == 4)
          fputs("[...]\n",stdout);
      dptr++;
      }
    }
  }
exit(0);
}
