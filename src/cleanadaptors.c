/* cleanadaptors.c: attempt at writing something to scan for adaptor sequences
in Illumina reads.  PAS: May-2011 */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>

#ifndef NO_ZLIB
#include <zlib.h>
#endif

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "mrg_fns.h"
#include "rmapbsbed2cpg.h"
#include "bin_cnts.h"
#include "fsm_ops.h"

/* local defines */
/* #define PROG_VERSN 1.00 */
/* first release: June 2011 */
/* #define PROG_VERSN 1.01 */
/* don't write 0 length seqs: Jul-2011 */
/* #define PROG_VERSN 1.02 */
/* further refine length exclusion: Jul-2011 */
/* #define PROG_VERSN 1.03 */
/* allow N-filling of trims: Jul-2011 */
/* #define PROG_VERSN 1.04 */
/* allow quality lines starting with '>' */
/* #define PROG_VERSN 1.04 */
/* allow trimming back before adaptor: Jan-2012 */
/* #define PROG_VERSN 1.10 */
/* trimming paired-end reads from two files, maintaining pair matching:
 Jul-2012,Nov-2012,Feb-2014 */
/* #define PROG_VERSN 1.11 */
/* option to 5p truncate adaptors: Mar-2014 */
/* #define PROG_VERSN 1.20 */
/* incorporate quality trimming as well: Mar-2014 */
/* #define PROG_VERSN 1.21 */
/* incorporate length trim: Apr-2014 */
/* #define PROG_VERSN 1.22 */
/* incorporate zlib compression: Apr-2014 */
/* #define PROG_VERSN 1.23 */
/* modify to work around odd zlib problem: May-2014 */
/* #define PROG_VERSN 1.24 */
/* to skip gzbuffer() call if we don't have new enough zlib (1.2.5 or later): Jun-2014 */
/* #define PROG_VERSN 1.25 */
/* attempt to simplify gzbuffer issue with NO_GZBUFFER def: Sep-2014 */
/* #define PROG_VERSN 1.26 */
/* replace pointer value format specs with %p: Dec-2014 */
/* #define PROG_VERSN 1.27 */
/* correct final fasta format adaptor read: Apr-2015 */
/* #define PROG_VERSN 1.28 */
/* produce help message if no parameters: Sep-2015 */
/* #define PROG_VERSN 1.30 */
/* repeated scans on each read: Mar-2016 */
/* #define PROG_VERSN 1.31 */
/* correct null char fault: Sep-2016 */
/* #define PROG_VERSN 1.32 */
/* correct readlength issue for no adaptor match reads: Oct-2016 */
/* #define PROG_VERSN 1.33 */
/* track down uninitialised read start position problem: Oct-2016 */
/* #define PROG_VERSN 1.34 */
/* correct no_zlib compilation: Apr-2017 */
#define PROG_VERSN 1.35
/* correct read2 hitlist init problem: Jun-2017 */

#define MAX_ADAPTOR_LEN 256
#define MIN_LEAD_TRAIL 6
#define DEF_READBUFLEN 1024
#define DEF_3PMARGIN 0
#define DEF_MATCHCRITERION 85.0
#define HASH_LENGTH 6

#ifndef NO_ZLIB
/* compression lib buffer size: = 128k */
#define ZLIB_BUFFERLEN 131072
#endif

typedef struct CA_adhit    /* information about an adaptor hit */
  {
  int pos;                 /* where in seq line */
  FS_RESELT *resp;         /* result ptr */
  struct CA_adhit *nxthit;
  struct CA_adhit *prvhit;
  }
CA_ADHIT;   

typedef struct CA_adaptor   /* linked list element of original adaptor seqs */
  {
  char *adaptstr;          /* the sequence */
  int asqno;               /* its number */
  void *adnam;             /* its name if any */
  int alen;                /* its length */
  unsigned int hval;       /* hash value of start */
  void *rlist;             /* list of results pertaining to this adaptor */
  struct CA_adaptor *nxtadaptor;
  struct CA_adaptor *prvadaptor;
  }
CA_ADAPTOR;

typedef enum CA_adfilefmt   /* format of adaptor source file */
  {
  CA_adflfmt_list = 0,      /* simple list of seqs */
  CA_adflfmt_fasta          /* fasta format file */
  }
CA_ADFILEFMT;

typedef struct ca_adptrelt  /* element for a list of adaptor pointers */
  {
  CA_ADAPTOR *aptr;
  struct ca_adptrelt *nxtapelt;
  struct ca_adptrelt *prvapelt;
  }
CA_ADPTRELT;

typedef enum CA_fastqlntype
  {
  CA_fqlt_hdrline = 0,
  CA_fqlt_sqline,
  CA_fqlt_qulhdr,
  CA_fqlt_qualln,
  CA_fqlt_unk
  }
CA_FASTQLNTYPE;

typedef enum CA_outstyle
  {
  CA_out_fsmlist = 0,
  CA_out_listreads,
  CA_out_listhitreads,
  CA_out_trimhits,
  CA_out_nfill,           /* put 'N's in trimmed region */
  }
CA_OUTSTYLE;

typedef enum CA_srctype
  {
  CA_src_fastq = 0,
  CA_src_fasta
  }
CA_SRCTYPE;

typedef struct CA_srcdst_fl  /* datatype for io object */
  {
  FILE *txtfile;
#ifndef NO_ZLIB
  gzFile gzfile;
#endif
  int comprssd;       /* 1 if compressed: ie use gzfile */
  int srclincnt;      /* troubleshooting: count source lines */
  }
CA_SRCDST_FL;

typedef struct CA_runpars
  {
  int trimcnt;
  CA_SRCTYPE srcstyle;
  CA_OUTSTYLE ostyle;
  int minleadtrail;
  int readflbuflen;
  int margn3p;
  int maxadlen;
  CA_ADAPTOR *adaptlst;
  CA_ADAPTOR *rawalist;  /* adaptors first loaded here */
  int mismatches;
  int skipres;
  float matcriterion;
  int mintrimlen;     /* minimum length for output of trimmed seqs, 0=>no limit */
  int trimback;
  int pairedends;     /* 1 if we are doing a paired run */
  CA_SRCDST_FL *read1fl;
  CA_SRCDST_FL *read2fl;
  char *r1bufs[CA_fqlt_unk];    /* buffers for lines, read 1 */
  char *r2bufs[CA_fqlt_unk];    /* & read 2 */
  CA_SRCDST_FL *dstfl1;                 /* for paired run, need to define output files */
  CA_SRCDST_FL *dstfl2;
  int phredqlimit;              /* phred score quality limit: 0 -> no test */
  int qualbase;                    /* base for score quality - def to 33 */
  int p5trunc;                  /* length of each adaptor to 5' trim */
  int maxoutlen;                /* -l limit of output read length. 0=>no limit */
  CA_ADFILEFMT adflfmt;
  int compress;                 /* -z,-Z changes this on,off */
  int repeatscan;               /* -1 -> repeat scans till length stable */
  }
CA_RUNPARS;

/* local globals */

int debuglevel;
/* int trimcnt; */

/* code ... */

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
fprintf(fl,"%s v%.2f: scan Illumina reads for adaptor seqs: Trim FASTA/FASTQ files\n",
          pnam,PROG_VERSN);
fputs("Options:\n",fl);
fputs("     -i <adaptorfile> file of adaptor seqs (1/line)\n",fl);
fputs("     -I <adaptorfile> fasta file of adaptor seqs\n",fl);
fprintf(fl,"     -a max adaptor length (def=%d)\n",MAX_ADAPTOR_LEN);
fprintf(fl,"     -m min leading match for adaptor (def=%d)\n",MIN_LEAD_TRAIL);
fprintf(fl,"     -M margin at read end for trimming (def=%d)\n",DEF_3PMARGIN);
fprintf(fl,"     -R readfile buffer length (def=%d)\n",DEF_READBUFLEN);
fputs("     -s <skipres> skip <skipres> on each read before checking matches\n",fl);
fprintf(fl,"     -p <%%> %% match threshold with adaptor sequence for hit (def=%.1f%%)\n",
          DEF_MATCHCRITERION);
fputs("     -S enable single base mismatches (def=disallow)\n",fl);
fputs("     -L print fsm to stdout\n",fl);
fputs("     -F <readfile>: run scan on <readfile> FASTQ/FASTA fmt, trim matching ends\n",fl);
fputs("     -x <lengthlimit>: only save trimmed reads exceeding length limit (def=1)\n",fl);
fputs("     -t <3'trimlength>: take further 3'trimlength bases before adaptor match (def=0)\n",fl);
fputs("     -T <5'adaptortrim> delete 5' bp from the start of each adaptor (def=don't)\n",fl);
fputs("     -q <phredscore> quality trim reads to longest section > Q<phredscore> (def=don't)\n",fl);
fputs("     -Q <qualbase> change base for quality scores (def=33)\n",fl);
fputs("     -r <repeatno> restrict scans to <repeatno> for each read (def=2)\n",fl);
fputs("     -N <readfile>: fill lines with 'N's rather than trimming\n",fl);
fputs("     -f <readfile>: run scan on <readfile> FASTQ/FASTA fmt, show all reads, indicate matches\n",fl);
fputs("     -H <readfile>: run scan on <readfile> FASTQ/FASTA fmt, indicate matches on hit reads only\n",fl);
fputs("     -G <read2file>: paired end input: trim both ends, retain valid paired reads\n",fl);
fputs("     -l <trimlength> max length for output reads (def=nolimit)\n",fl);
#ifndef NO_ZLIB
fputs("     -z enable gzip compression input & output. Positional (def=nocompression)\n",fl);
fputs("     -Z disable gzip compression. Positional: effects following -F,-G,-o,-O files (def=-Z)\n",fl);
#endif
fputs("     -o <dstfile>: write output to dstfile (read1 data for paired ends) def=stdout\n",fl);
fputs("     -O <r2dstfile>: output file for read2 - needed for paired run\n",fl);
fputs("       if <readfile> is '-', then use stdin, for uncompressed input\n",fl);
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

int ca_exponent(int base,
                int p)
/* recursively determine base**p for p>=0
& integer base */
{
if (p <= 0)
  return(1);
else
  return(base*ca_exponent(base,p-1));
}

unsigned int ca_hashnastring(char *sq,
                             int slen,
                             int hlen)
/* return a hash value for sq up to hlen
residues */
{
unsigned int hval;
int i;
int sp;

i = hlen - 1;
hval = 0;
sp = 0;
/*
while (i >= 0)
  {
  if (sp < slen)
    hval += (tr_nares2int(*(sq+sp)+1) * ca_exponent(4,i);
  i--;
  sp++;
  }
*/
while (i >= 0)
  {
  if (sp < slen)
    hval = hval*4 + tr_nares2int(*(sq+sp)) + 1;
  i--;
  sp++;
  }
return(hval);
}

CA_ADAPTOR *ca_appndadaptor(CA_ADAPTOR **alst,
                            char *astr,
                            void *dptr,
                            int asno)
/* append a new string element to *alst.
strdup astr in order that a conserved 
value is created.  Return address of new element
in case it is useful */
{
CA_ADAPTOR *prev, *end_ptr;

if (alst != NULL)
  {
  prev = end_ptr = *alst;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtadaptor;
    }
  end_ptr = (CA_ADAPTOR *) getmemory(sizeof(CA_ADAPTOR),"Adaptor elt");
  end_ptr->nxtadaptor = NULL;
  end_ptr->adaptstr = bas_strdup(astr);
  end_ptr->adnam = dptr;
  end_ptr->alen = strlen(astr);
  end_ptr->hval = ca_hashnastring(astr,end_ptr->alen,HASH_LENGTH);
  end_ptr->asqno = asno;
  end_ptr->rlist = NULL;
  if (*alst == NULL)
    {
    *alst = end_ptr;
    end_ptr->prvadaptor = NULL;
    }
  else
    {
    prev->nxtadaptor = end_ptr;
    end_ptr->prvadaptor = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void ca_deladaptelt(CA_ADAPTOR *ep,
                    CA_ADAPTOR **lstrt)
/* delete ep from list *lstrt */
{
CA_ADAPTOR *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvadaptor) == NULL)
    *lstrt = ep->nxtadaptor;
  else
    pt->nxtadaptor = ep->nxtadaptor;
  if ((pt = ep->nxtadaptor) != NULL)
    pt->prvadaptor = ep->prvadaptor;
  if (ep->adaptstr != NULL)
    memfree(ep->adaptstr);
  memfree(ep);
  }
}

void ca_clralladaptors(CA_ADAPTOR **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  ca_deladaptelt(*lstrt,lstrt);
}

CA_ADAPTOR *ca_adaptor4sq(CA_ADAPTOR *alst,
                          char *sq)
/* iteratively scan alst for sq, returning the element
containing the match */
{
CA_ADAPTOR *lp;
unsigned int sqhash;
int sqlen;

lp = alst;
sqlen = strlen(sq);
sqhash = ca_hashnastring(sq,sqlen,HASH_LENGTH);
while (lp != NULL)
  if ((sqlen == lp->alen) && (sqhash == lp->hval) &&
        (strcasecmp(lp->adaptstr,sq) == 0))    /* have a match */
    return(lp);
  else
    lp = lp->nxtadaptor;
/* bombed, no match */
return(NULL);
}

CA_ADPTRELT *ca_appndadptr(CA_ADPTRELT **aplst,
                           CA_ADAPTOR *ap)
/* append a new adaptor element to *aplst.
Return address of new element
in case it is useful */
{
CA_ADPTRELT *prev, *end_ptr;

if (aplst != NULL)
  {
  prev = end_ptr = *aplst;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtapelt;
    }
  end_ptr = (CA_ADPTRELT *) getmemory(sizeof(CA_ADPTRELT),"Adaptor ptrN elt");
  end_ptr->nxtapelt = NULL;
  end_ptr->aptr = ap;
  if (*aplst == NULL)
    {
    *aplst = end_ptr;
    end_ptr->prvapelt = NULL;
    }
  else
    {
    prev->nxtapelt = end_ptr;
    end_ptr->prvapelt = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

CA_ADPTRELT *ca_adptrelt4aptr(CA_ADPTRELT *alist,
                              CA_ADAPTOR *ap)
/* if an element pertaining to ap is on alist
return the address, else NULL */
{
CA_ADPTRELT *pp;

pp = alist;
while (pp != NULL)
  if (pp->aptr == ap)
    return(pp);
  else
    pp = pp->nxtapelt;
return(NULL);
}

CA_ADPTRELT *ca_chkappndadptr(CA_ADPTRELT **alist,
                              CA_ADAPTOR *ap)
/* check if this ap is already in *alist, if not
append it, else return the element */
{
CA_ADPTRELT *pp;

if ((pp = ca_adptrelt4aptr(*alist,ap)) != NULL)
  return(pp);
else
  return(ca_appndadptr(alist,ap));
}

int ca_buildadaptorfsm(int sqno,
                       CA_ADAPTOR *adpt,
                       CA_RUNPARS *rpars,
                       FS_FSMSTRCT *adfsm)
/* scan adpt->adaptstr as an adaptor sequence, making
FSM adfsm in accordance with parameters.
return No. of fsm entries */
{
int sqlen;
char *leadfragbuf;
int sqp;
SQ_RESTYPE sres;
int frgp;
char cacheres;
char subres;
int fsmcnt;
CA_ADAPTOR *ap;

fsmcnt = 0;
if ((sqlen = strlen(adpt->adaptstr)) > 0)
  {
  leadfragbuf = (char *) getmemory(sqlen+1,"lead seq buf");
  for (sqp = rpars->minleadtrail; sqp <= sqlen; sqp++)
    {
    (void) memcpy(leadfragbuf,adpt->adaptstr,sqp);
    *(leadfragbuf + sqp) = '\0'; 
    for (frgp = 0; frgp < sqp; frgp++)
      if (rpars->mismatches)
        {
        cacheres = *(leadfragbuf+frgp);
        for (sres = RES_a; sres <= RES_t; sres++)
          {
          *(leadfragbuf+frgp) = sqfl_restype2chr(sres);
          if ((ap = ca_adaptor4sq(rpars->adaptlst,leadfragbuf)) == NULL)
            {
            ap = ca_appndadaptor(&rpars->adaptlst,leadfragbuf,adpt,fsmcnt+sqno);
            (void) ca_chkappndadptr((CA_ADPTRELT**) &ap->rlist,adpt);
            fs_adddatprs(adfsm,leadfragbuf,ap->rlist);
            fsmcnt++;
            }
          else
            (void) ca_chkappndadptr((CA_ADPTRELT**) &ap->rlist,adpt);
          }
        *(leadfragbuf+frgp) = cacheres;
        }
      else
        {
        if ((ap = ca_adaptor4sq(rpars->adaptlst,leadfragbuf)) == NULL)
          {
          ap = ca_appndadaptor(&rpars->adaptlst,leadfragbuf,adpt,fsmcnt+sqno);
          (void) ca_chkappndadptr((CA_ADPTRELT **) &ap->rlist,adpt);
          fs_adddatprs(adfsm,leadfragbuf,ap->rlist);
          fsmcnt++;
          }
        else
          (void) ca_chkappndadptr((CA_ADPTRELT **) &ap->rlist,adpt);
        }
     }
  }
memfree(leadfragbuf);
return(fsmcnt);
}

CA_ADAPTOR *ca_chkappndadaptor(CA_ADAPTOR **alst,
                               char *astr,
                               char *aname,
                               int *asno)
/* check if this astr is new: if so append a new 
string element to *alst.  If not new, return the
list pointer for it */
{
CA_ADAPTOR *ap;

if ((ap = ca_adaptor4sq(*alst,astr)) == NULL)
  {
  (*asno)++;
  return(ca_appndadaptor(alst,astr,aname,*asno));
  }
else
  return(ap);
}

int ca_readadaptorfile(FILE *sfl,
                       CA_RUNPARS *rpars,
                       FS_FSMSTRCT *adfsm)
/* scan sfl for adaptor sequences one seq/line.  Load
each original seq and all partitions starting from
position minleadtail up to minleadtail from end into
adfsm. return no of seqs read */
{
int sqno;
char *sqbuf;
CA_ADAPTOR *alstptr;
char *sbp;
char *hdr;
char *rdbuf;
char *sqp;

sqno = 0;
sqbuf = (char *) getmemory(rpars->maxadlen+1,"adaptor seq buf");
*sqbuf = '\0';
switch (rpars->adflfmt)
  {
  case CA_adflfmt_fasta:
    hdr = (char *) getmemory(rpars->maxadlen+1,"adaptor hdr buf");
    rdbuf = (char *) getmemory(rpars->maxadlen+1,"adaptor input buf");
    while ((sbp = fgets(rdbuf,rpars->maxadlen,sfl)) != NULL)
      {
      while (*sbp != '\0')
        {
        if (*sbp == '\n')
          *sbp = '\0';
        sbp++;
        }
      switch (*rdbuf)
        {
        case '>':      /* hdr line, note name */
          if (strlen(sqbuf) > rpars->p5trunc)
            {
            alstptr = ca_chkappndadaptor(&rpars->rawalist,(sqbuf+rpars->p5trunc),
                                           bas_strdup(hdr),&sqno);
            *sqbuf = '\0';
            }
          strncpy(hdr,rdbuf+1,rpars->maxadlen);
          break;
        default:      /* everything else */
          strncat(sqbuf,rdbuf,rpars->maxadlen);
          break;
        }
      if ((sqno > 0) && (strlen(sqbuf) > rpars->p5trunc))
        alstptr = ca_chkappndadaptor(&rpars->rawalist,(sqbuf+rpars->p5trunc),
                                       bas_strdup(hdr),&sqno);
      }
/* do we have a remaining adaptor? */
    if (strlen(sqbuf) > rpars->p5trunc)
      alstptr = ca_chkappndadaptor(&rpars->rawalist,(sqbuf+rpars->p5trunc),
                                     bas_strdup(hdr),&sqno);
    memfree(hdr);
    memfree(rdbuf);
    break;
  case CA_adflfmt_list:
  default:
    while ((sbp = fgets(sqbuf,rpars->maxadlen,sfl)) != NULL)
      {
      while (*sbp != '\0')
        {
        if (*sbp == '\n')
          *sbp = '\0';
        sbp++;
        }
      if (strlen(sqbuf) > rpars->p5trunc)
        alstptr = ca_chkappndadaptor(&rpars->rawalist,(sqbuf+rpars->p5trunc),"",&sqno);
      }
    break;
  }
memfree(sqbuf);
rpars->adaptlst = NULL;
alstptr = rpars->rawalist;
sqno = 0;
while (alstptr != NULL)
  {
  sqno += ca_buildadaptorfsm(sqno,alstptr,rpars,adfsm);
  alstptr = alstptr->nxtadaptor;
  }
return(sqno);
}

void ca_nl(FILE *fl)
{
fputc('\n',fl);
}

void ca_prtctxt(CA_SRCDST_FL *fl,
                void *p)
{
if (!fl->comprssd)
  (void) fs_lstcstr(fl->txtfile,p);
}

void ca_prtlval(CA_SRCDST_FL *fl,
                void *p)
/* print *p as a hexadec long int */
{
#ifndef NO_ZLIB
if (fl->comprssd)
  gzprintf(fl->gzfile,"%p", p);
else
#endif
  fprintf(fl->txtfile,"%p", p);
}

int ca_inrankordr(CA_ADHIT *e1p,
                  CA_ADHIT *e2p)
/* if e1 preceeds e2 in rank order, then return 1. */
{
if (e1p == NULL)
  return(0);
else
  if (e2p == NULL)
    return(1);
  else
    return((e1p->resp == e2p->resp) && (e1p->pos < e2p->pos));
}

CA_ADHIT *ca_insertnewhit(CA_ADHIT **alst,
                          int lpos,
                          FS_RESELT *rep,
                          int (* preceeds_fn)(CA_ADHIT *e1p,
                                              CA_ADHIT *e2p))
/* scan through *alst, if exists, finding the first element
which is after lpos, insert a new element before it, returning
the address of the new element.  If we fall off the end of the list,
then append the new element there. */
{
CA_ADHIT *prvelt;
CA_ADHIT *ep;
CA_ADHIT *newp;

prvelt = NULL;
ep = *alst;
newp = (CA_ADHIT *) getmemory(sizeof(CA_ADHIT),"ADhitelt");
newp->pos = lpos;
newp->resp = rep;
while (ep != NULL)
  if ((*preceeds_fn)(newp,ep))
    {
    newp->nxthit = ep;
    if (ep->prvhit == NULL)
      {
      *alst = newp;
      newp->prvhit = NULL;
      }
    else
      {
      ep->prvhit->nxthit = newp;
      newp->prvhit = ep->prvhit;
      }
    ep->prvhit = newp;
    return(newp);
    }
  else
    {
    prvelt = ep;
    ep = ep->nxthit;
    }
/* haven't found it, stick on end */
newp->nxthit = NULL;
newp->prvhit = prvelt;
if (prvelt == NULL)
  *alst = newp;
else
  prvelt->nxthit = newp;
return(newp);
}

CA_ADHIT *ca_hit4pos5p(CA_ADHIT *hlst,
                       int hpos)
/* return the first hit element on hlst which
matches 5' hpos */
{
CA_ADHIT *hp;

hp = hlst;
while (hp != NULL)
  if (hp->pos == hpos)
    return(hp);
  else
    hp = hp->nxthit;
return(NULL);
}

int ca_len4fs_reselt(FS_RESELT *rep)
  /* return the data length for rep */
{
FS_DATPRELT *dpep;

if ((rep != NULL) && ((dpep = (FS_DATPRELT *) rep->action) != NULL))
  return(dpep->ldstr);
else
  return(0);
}

void ca_killhitelt(CA_ADHIT **alst,
                   CA_ADHIT *hep)
/* remove *hep from alst */
{
if (hep != NULL)
  {
  if (hep->prvhit == NULL)
    *alst = hep->nxthit;
  else
    hep->prvhit->nxthit = hep->nxthit;
  if (hep->nxthit != NULL)
    hep->nxthit->prvhit = hep->prvhit;
  memfree(hep);
  }
}

void ca_killhitlst(CA_ADHIT **alst)
  /* iteratively kill all of alst */
{
while (*alst != NULL)
  ca_killhitelt(alst,*alst);
}

CA_ADHIT *ca_hit4pos3p(CA_ADHIT *hlst,
                       int pos3p)
/* return the first hit element on hlst which
matches 3' pos3p */
{
CA_ADHIT *hp;
int p3p;

hp = hlst;
while (hp != NULL)
  {
  p3p = hp->pos + ca_len4fs_reselt(hp->resp) - 1;
  if (p3p == pos3p)
    return(hp);
  else
    hp = hp->nxthit;
  }
return(NULL);
}

CA_ADHIT *ca_hitat3pend(CA_ADHIT *hlst,
                        int sqlen,
                        int margin)
/* if hlst contains any hits which lie within
margin of end of seq line (sqlen long), then
return point to that hit element, else NULL */
{
CA_ADHIT *hp;

hp = hlst;
while (hp != NULL)
  if ((hp->pos + ca_len4fs_reselt(hp->resp)) >= (sqlen-margin))
    return(hp);
  else
    hp = hp->nxthit;
return(NULL);
}

float ca_hitpercntmatch(CA_RUNPARS *rpp,
                        char *sqread,
                        CA_ADHIT *hitp,
                        int *olapcnt,
                        int *matcnt)
/* scan this hit vs sqread, counting
number of matched residues.  Return
% of matches by end of read, slen
residues or original
adaptor seq.  If olapcnt & matcnt are non-NULL
then return appropriate counts to them.
return 0.0 for issues */
{
int ocnt;
int mcnt;
char *ap;
char *sp;
FS_DATPRELT *dpep;
CA_ADPTRELT *paep;

if (hitp != NULL)
  {
  dpep = (FS_DATPRELT *) hitp->resp->action;
  paep = (CA_ADPTRELT *) dpep->stxt;
  ap = paep->aptr->adaptstr; /* + rpp->p5trunc; */
  sp = sqread + hitp->pos;
  ocnt = mcnt = 0;
  while ((*ap != '\0') && (*sp != '\0'))
    {
    if (toupper(*ap++) == toupper(*sp++))
      mcnt++;
    ocnt++;
    }
  if (olapcnt != NULL)
    *olapcnt = ocnt;
  if (matcnt != NULL)
    *matcnt = mcnt;
  if (ocnt > 0)
    return((float) mcnt*100.0/ocnt);
  else
    return(0.0);
  }
else
  return(0.0);
}

int ca_hitmeetscriteria(char *sqread,
                        CA_ADHIT *hitp,
                        CA_RUNPARS *rpars)
/* return 1 if this match meets criteria for a match */
{
float mprcnt;
int p3gap;
FS_DATPRELT *dpep;
CA_ADPTRELT *paep;
CA_ADAPTOR *ap;

mprcnt = ca_hitpercntmatch(rpars,sqread,hitp,NULL,NULL);
/* if (hitp != NULL)
  fprintf(stdout,"%.2f at %d\n",mprcnt,hitp->pos); */
if ((hitp!= NULL) && (mprcnt >= rpars->matcriterion))
  if (rpars->margn3p <= 0)
    return(1);
  else
    {
    dpep = (FS_DATPRELT *) hitp->resp->action;
    paep = (CA_ADPTRELT *) dpep->stxt;
    p3gap = hitp->pos + strlen(paep->aptr->adaptstr) - strlen(sqread);
    return(p3gap <= rpars->margn3p);
    }
return(0);
}

int ca_ahitmeetscriteria(char *sqread,
                         CA_ADHIT *hlst,
                         CA_RUNPARS *rpars)
/* check all elements of hitlst and return if any
matches criteria for a match */
{
if (hlst == NULL)
  return(0);
else
  if (ca_hitmeetscriteria(sqread,hlst,rpars))
    return(1);
  else
    return(ca_ahitmeetscriteria(sqread,hlst->nxthit,rpars));
}

CA_ADHIT *ca_5pmostvalidhit(CA_ADHIT *hlst,
                            char *sqread,
                            CA_RUNPARS *rpars)
/* scan any hits in hlst, checking
for validity.  Return the 5'-most such
hit, NULL if none */
{
CA_ADHIT *h5p;
CA_ADHIT *hp;

h5p = NULL;
hp = hlst;
while (hp != NULL)
  {
  if (ca_hitmeetscriteria(sqread,hp,rpars))
    if ((h5p == NULL) || (h5p->pos > hp->pos))
      h5p = hp;
  hp = hp->nxthit;
  }
return(h5p);
}

char ca_phredscore2chr(CA_RUNPARS *rpp,
                       int pscore)
/* return the character corresponding to phred
score pscore. */
{
return((char) (pscore+rpp->qualbase));
}

int ca_phredscorechr2int(CA_RUNPARS *rpp,
                         char pchr)
/* return the phredscore corresponding to
pchr. */
{
return((int) pchr - rpp->qualbase);
}

int ca_readisvalid(CA_RUNPARS *rpp,
                   char *sqread,
                   char *qline,
                   CA_ADHIT *hlst,
                   int *vstart,
                   int *vend)
/* do final length checks for most 5' valid
hit.  If vlen is non-null, return the effective
valid length to it.  Also do quality scan
if required and adjust vstart/vend as
appropriate */
{
CA_ADHIT *h5p;
int fnllend;
int hiqstart;
int hiqend;
int hiqlen;
int sp;
int beststart;
int bestend;

if ((h5p = ca_5pmostvalidhit(hlst,sqread,rpp)) != NULL)
  fnllend = h5p->pos - rpp->trimback;
else
  fnllend = strlen(sqread);
hiqstart = 0;
if (rpp->phredqlimit > 0)
  {
  beststart = bestend = hiqend = hiqlen = 0;
  for (sp=0; sp < fnllend; sp++)
    {
    if (ca_phredscorechr2int(rpp,*(qline + sp)) < rpp->phredqlimit)
      {
      if (hiqlen > (bestend - beststart + 1))
        {
        beststart = hiqstart;
        bestend = hiqend;
        }
      hiqend = hiqstart = sp + 1;
      hiqlen = 0;
      }
    else
      {
      hiqend = sp + 1;
      hiqlen++;
      }
    }
  if (hiqlen > (bestend - beststart + 1))
    {
    bestend = hiqend;
    beststart = hiqstart;
    }
  fnllend = imin(bestend,fnllend);
  }
if (vstart != NULL)
  {
  if (rpp->phredqlimit <=0)
    *vstart = 0;
  else
    *vstart = beststart;
  }
if (vend != NULL)
  *vend = fnllend;
return(((rpp->phredqlimit <= 0) && (fnllend >= rpp->mintrimlen)) ||
         ((rpp->phredqlimit > 0) && (fnllend - beststart) >= rpp->mintrimlen));
}

int ca_opensrcdststrct(CA_SRCDST_FL *srcdst,
                       int comprss,
                       char *fname,
                       char *fmode)
/* depending on comprss, either open fname in
fmode as txt file or gzfile - return 1 if OK.
if compressed, set buffer to ZLIB_BUFFERLEN */
{
srcdst->srclincnt = 0;
if (comprss)
#ifndef NO_ZLIB
  if ((srcdst->gzfile = gzopen(fname,fmode)) != NULL)
    {
    srcdst->comprssd = 1;
#ifndef NO_GZBUFFER
/* hack to avoid zlib version problems */
    (void) gzbuffer(srcdst->gzfile,ZLIB_BUFFERLEN);
#endif
    return(1);
    }
  else
    return(0);
#else
  err_msg_die("Compression not enabled for this compile: nozlib set\n");
#endif
else
  if ((srcdst->txtfile = fopen(fname,fmode)) != NULL)
    {
    srcdst->comprssd = 0;
    return(1);
    }
  else
    return(0);
}

void ca_closesrcdststrct(CA_SRCDST_FL *srcdst)
  /* depending on settings, close the stream */
{
#ifndef NO_ZLIB
if (srcdst->comprssd)
  {
  if (srcdst->gzfile != NULL)
    gzclose(srcdst->gzfile);
  srcdst->gzfile = NULL;
  }
else
#endif
  {
  if (srcdst->txtfile != NULL)
    fclose(srcdst->txtfile);
  srcdst->txtfile = NULL;
  }
}

CA_SRCDST_FL *ca_getsrcdststrct(char *what)
  /* allocate and zero a source/dest
data structure */
{
CA_SRCDST_FL *casdfl;

casdfl = (CA_SRCDST_FL *) getmemory(sizeof(CA_SRCDST_FL),what);
casdfl->txtfile = NULL;
#ifndef NO_ZLIB
casdfl->gzfile = NULL;
#endif
casdfl->comprssd = 0;
return(casdfl);
}

int ca_getc(CA_SRCDST_FL *srcdst)
  /* fgetc equivalent for txt/gz */
{
#ifndef NO_ZLIB
if (srcdst->comprssd)
  return(gzgetc(srcdst->gzfile));
else
#endif
  return(fgetc(srcdst->txtfile));
}

int ca_gets(CA_SRCDST_FL *srcdst,
            char *buf,
            int blen)
/* call appropriate gets function.
return 1 for success, else 0.
remove <nl> from string */
{
char *bp;

#ifndef NO_ZLIB
if (srcdst->comprssd)
  bp = gzgets(srcdst->gzfile,buf,blen);
else
#endif
  bp = fgets(buf,blen,srcdst->txtfile);
if (bp == NULL)        /* EOF or error */
  return(0);
else
  {     /* contains <nl> so replace it */
  while ((*bp != '\n') && (*bp != '\0') &&
           ((bp - buf) < blen))
    bp++;
  if ((bp - buf) < blen)
    *bp = '\0';
  srcdst->srclincnt++;  
  return(1);
  }
}

int ca_puts(char *str,
            CA_SRCDST_FL *srcdst)
/* fputs equivalent for txt/gz */
{
#ifndef NO_ZLIB
if (srcdst->comprssd)
  return(gzputs(srcdst->gzfile,str));
else
#endif
  return(fputs(str,srcdst->txtfile));
}

int ca_putc(char c,
            CA_SRCDST_FL *srcdst)
/* fputc equivalent for txt/gz */
{
#ifndef NO_ZLIB
if (srcdst->comprssd)
  return(gzputc(srcdst->gzfile,c));
else
#endif
  return(fputc(c,srcdst->txtfile));
}

void ca_rptchrout(CA_SRCDST_FL *dfl,
                  char outc,
                  int ccnt)
/* write ccnt outcX to dfl */
{
int cc;

cc = ccnt;
while (cc > 0)
  {
  ca_putc(outc,dfl);
  cc--;
  }
}

void ca_writelstreadout(CA_SRCDST_FL *dfile,
                        char *linbuf[],
                        CA_ADHIT *hitlst,
                        CA_RUNPARS *rpars)
/* logic to write out the data from
linbuf[] according to requirements */
{
CA_ADHIT *hlp;
int bcnt;
FS_DATPRELT *dpep;
CA_FASTQLNTYPE fqlt;
int slen;
char *sp;
char *ap;
int olapcnt;
int matcnt;
float pcntmat;
int origlen;
CA_FASTQLNTYPE mxlt;
CA_ADPTRELT *paep;

fqlt = CA_fqlt_hdrline;
switch (rpars->srcstyle)
  {
  case CA_src_fasta:
    mxlt = CA_fqlt_sqline;
    break;
  case CA_src_fastq:
  default:
    mxlt = CA_fqlt_qualln;
    break;
  }
origlen = slen = strlen(linbuf[CA_fqlt_sqline]);
if (slen >= 0)
  {
  while (fqlt <= mxlt)
    {
    switch (fqlt)
      {
      case CA_fqlt_hdrline:
      case CA_fqlt_qulhdr:
        ca_puts(linbuf[fqlt],dfile);
        ca_putc('\n',dfile);
        break;
      case CA_fqlt_sqline:
        bcnt = slen;
        sp = linbuf[fqlt];
        while (bcnt > 0)
          {
          ca_putc(*sp++,dfile);
          bcnt--;
          }
        ca_putc('\n',dfile);
        hlp = hitlst;
        while (hlp != NULL)
          {
          if (ca_hitmeetscriteria(linbuf[CA_fqlt_sqline],hlp,rpars))
            {
            if (hlp->pos > 1)
              ca_putc('|',dfile);
            ca_rptchrout(dfile,' ',(hlp->pos - 1));
            ca_rptchrout(dfile,'^',ca_len4fs_reselt(hlp->resp));
            ca_putc('\n',dfile);
            dpep = (FS_DATPRELT *) hlp->resp->action;
            paep = (CA_ADPTRELT *) dpep->stxt;
            ap = paep->aptr->adaptstr;
            if (hlp->pos > 1)
              ca_putc('|',dfile);
            ca_rptchrout(dfile,' ',(hlp->pos - 1));
            sp = linbuf[fqlt] + hlp->pos;
            while ((*ap != '\0') && (*sp != '\0'))
              {
              if (toupper(*ap) == toupper(*sp++))
                ca_putc(toupper(*ap),dfile);
              else
                ca_putc(tolower(*ap),dfile);
              ap++;
              }
            pcntmat = ca_hitpercntmatch(rpars,linbuf[fqlt],hlp,&olapcnt,&matcnt);
#ifndef NO_ZLIB
            if (dfile->comprssd)
              gzprintf(dfile->gzfile," %d/%d match (%.1f%%) ",matcnt,olapcnt,pcntmat);
            else
#endif
              fprintf(dfile->txtfile," %d/%d match (%.1f%%) ",matcnt,olapcnt,pcntmat);
            ca_prtctxt(dfile,(char *)paep->aptr->adnam);
            ca_putc('\n',dfile);
            }
          hlp = hlp->nxthit;
          }
        break;
      case CA_fqlt_qualln:
        bcnt = slen;
        sp = linbuf[fqlt];
        while (bcnt > 0)
          {
          ca_putc(*sp++,dfile);
          bcnt--;
          }
        if (rpars->ostyle == CA_out_nfill)
          ca_rptchrout(dfile,';',(origlen - slen));
        ca_putc('\n',dfile);
        break;
      default:
        break;
      }
    fqlt++;
    }
  }
}

void ca_writetrmreadout(CA_SRCDST_FL *dfile,
                        char *linbuf[],
                        CA_RUNPARS *rpars,
                        int sqstart,
                        int sqend)
/* logic to write out the data from
linbuf[] according to requirements */
{
int bcnt;
CA_FASTQLNTYPE fqlt;
char *sp;
int origlen;
CA_FASTQLNTYPE mxlt;
int qlen;

fqlt = CA_fqlt_hdrline;
switch (rpars->srcstyle)
  {
  case CA_src_fasta:
    mxlt = CA_fqlt_sqline;
    break;
  case CA_src_fastq:
  default:
    mxlt = CA_fqlt_qualln;
    break;
  }
origlen = strlen(linbuf[CA_fqlt_sqline]);
if ((sqend-sqstart) < origlen)
  rpars->trimcnt++;
while (fqlt <= mxlt)
  {
  switch (fqlt)
    {
    case CA_fqlt_hdrline:
    case CA_fqlt_qulhdr:
      ca_puts(linbuf[fqlt],dfile);
      ca_putc('\n',dfile);
      break;
    case CA_fqlt_sqline:
      if (rpars->ostyle == CA_out_nfill)
        ca_rptchrout(dfile,'N',sqstart);
      bcnt = sqend - sqstart;
      if (rpars->maxoutlen > 0)
        bcnt = imin(bcnt,rpars->maxoutlen);
      sp = linbuf[fqlt] + sqstart;
      qlen = bcnt = imin(bcnt,strlen(sp));
      while (bcnt > 0)
        {
        if (*sp != '\0')
          ca_putc(*sp++,dfile);
        bcnt--;
        }
      if (rpars->ostyle == CA_out_nfill)
        ca_rptchrout(dfile,'N',(origlen - sqend));
      ca_putc('\n',dfile);
      break;
    case CA_fqlt_qualln:
      if (rpars->ostyle == CA_out_nfill)
        bcnt = origlen;
      else
        bcnt = imin((sqend - sqstart),qlen);
      if (rpars->maxoutlen > 0)
        bcnt = imin(bcnt,rpars->maxoutlen);
      sp = linbuf[fqlt] + sqstart;
      while (bcnt > 0)
        {
        if (*sp != '\0')
          ca_putc(*sp++,dfile);
        bcnt--;
        }
      ca_putc('\n',dfile);
      break;
    default:
      break;
    }
  fqlt++;
  }
/* if (!dfile->comprssd)
  fflush(dfile->txtfile); */
}

int ca_readflread(CA_SRCDST_FL *rdfl,
                  FS_FSMSTRCT *adfsm,
                  CA_RUNPARS *rpars,
                  int *lcnt,
                  char *linbufs[])
/* read a complete read from rdfl as either fasta or fastq fmt
each read from skipres.  Return 1 if anything read  */
{
char nc;
int lpos;
CA_FASTQLNTYPE fqlt;
char *dstp;
char *sqp;

fqlt = CA_fqlt_hdrline;
if (adfsm != NULL)
  fs_initrun(adfsm);
dstp = linbufs[fqlt];
while (ca_gets(rdfl,dstp,rpars->readflbuflen) > 0)
  {
  switch (fqlt)
    {
    case CA_fqlt_sqline:
      if (rpars->srcstyle == CA_src_fasta) /* bail out now */
        return(1);
      break;
    case CA_fqlt_hdrline:
      if (*lcnt == 0)
        switch (*linbufs[CA_fqlt_hdrline])
          {
          case '>':
            rpars->srcstyle = CA_src_fasta;
            break;
          case '@':
          default:
            rpars->srcstyle = CA_src_fastq;
            break;
          }
      break;
    case CA_fqlt_qualln:
      return(1);
      break;
    case CA_fqlt_qulhdr:
    default:
      break;
    }
  (*lcnt)++;
  fqlt++;
  dstp = linbufs[fqlt];
  }
return(0);
}

CA_ADHIT *ca_scanstring(CA_RUNPARS *rpp,
                        char *str,
                        int slen,
                        FS_FSMSTRCT *adfsm)
/* scan string str to limit slen for hits;
return generated hit list */
{
CA_ADHIT *hlst;
char *sp;
int spos;
FS_RESELT *frp;
CA_ADHIT *hlp;
int mpos;

hlst = NULL;
sp = str;
spos = 1;
if (adfsm != NULL)
  fs_initrun(adfsm);
while ((*sp != '\0') && (spos <= slen))
  {
  if ((adfsm != NULL) && (spos > rpp->skipres))
    {
    if ((frp = fs_procchr(adfsm,*sp,tr_nares2int)) != NULL)
      {
      mpos = spos - ca_len4fs_reselt(frp);
      if (((hlp = ca_hit4pos5p(hlst,mpos)) != NULL) &&
           (ca_len4fs_reselt(frp) > ca_len4fs_reselt(hlp->resp)))
        hlp->resp = frp;    /* longest match */
      else
        (void) ca_insertnewhit(&hlst,mpos,frp,ca_inrankordr);
      }
    }
  sp++;
  spos++;
  }
return(hlst);
}

void ca_scanreads(CA_RUNPARS *rpp,
                  FS_FSMSTRCT *adfsm)
/* Generate line buffers, progressively get read1 (&read 2)
into buffers.  check for adaptors/length.  Write out 
compliant reads if single-end or both pass */
{
CA_FASTQLNTYPE fqlt;
int lcnt1;
int lcnt2;
CA_ADHIT *hlist1;
CA_ADHIT *hlist2;
int olen1;
int olen2;
int slen;
CA_ADHIT *p5mosthit;
int s1start;
int s2start;
int hcnt;
CA_ADAPTOR *ap;
FS_DATPRELT *pp;
CA_ADPTRELT *paep;
CA_ADHIT *ret_hlst;
char *rp;
FS_RESELT *reseltp;
FS_DATPRELT *datprp;
int scancnt;
int r1valid;
int r2valid;
/* int rcnt; */

/* generate buffers as required */
for (fqlt = CA_fqlt_hdrline; fqlt < CA_fqlt_unk; fqlt++)
  {
  rpp->r1bufs[fqlt] = (char *) getmemory(rpp->readflbuflen+1,"Read1buf");
  *(rpp->r1bufs[fqlt]) = '\0';
  if (rpp->pairedends)
    {
    rpp->r2bufs[fqlt] = (char *) getmemory(rpp->readflbuflen+1,"Read2buf");
    *(rpp->r2bufs[fqlt]) = '\0';
    }
  else
    rpp->r2bufs[fqlt] = NULL;
  }
lcnt1 = lcnt2 = 0;
/* rcnt = 0; */
while ((ca_readflread(rpp->read1fl,adfsm,rpp,&lcnt1,rpp->r1bufs) && (!rpp->pairedends)) ||
         (rpp->pairedends && ca_readflread(rpp->read2fl,adfsm,rpp,&lcnt2,rpp->r2bufs)))
  {
  switch (rpp->ostyle)
    {
    case CA_out_listhitreads:
      hlist1 = ca_scanstring(rpp,rpp->r1bufs[CA_fqlt_sqline],strlen(rpp->r1bufs[CA_fqlt_sqline]),
                               adfsm);
      if (rpp->pairedends)
        hlist2 = ca_scanstring(rpp,rpp->r2bufs[CA_fqlt_sqline],strlen(rpp->r2bufs[CA_fqlt_sqline]),
                                 adfsm);
      if (ca_ahitmeetscriteria(rpp->r1bufs[CA_fqlt_sqline],hlist1,rpp) &&
        (!rpp->pairedends || ca_ahitmeetscriteria(rpp->r2bufs[CA_fqlt_sqline],hlist2,rpp)))
        {
        ca_writelstreadout(rpp->dstfl1,rpp->r1bufs,hlist1,rpp);
        if (rpp->pairedends)
          ca_writelstreadout(rpp->dstfl2,rpp->r2bufs,hlist2,rpp);
        }
      break;
    case CA_out_listreads:
      hlist1 = ca_scanstring(rpp,rpp->r1bufs[CA_fqlt_sqline],strlen(rpp->r1bufs[CA_fqlt_sqline]),
                               adfsm);
      if (rpp->pairedends)
        hlist2 = ca_scanstring(rpp,rpp->r2bufs[CA_fqlt_sqline],strlen(rpp->r2bufs[CA_fqlt_sqline]),
                                 adfsm);
      ca_writelstreadout(rpp->dstfl1,rpp->r1bufs,hlist1,rpp);
      if (rpp->pairedends)
        ca_writelstreadout(rpp->dstfl2,rpp->r2bufs,hlist2,rpp);
      break;
    case CA_out_trimhits:
    case CA_out_nfill:
    default:
/* rcnt++;
fprintf(stdout,"Read# %d\n",rcnt); */
      rp = rpp->r1bufs[CA_fqlt_sqline];
      ret_hlst = hlist1 = NULL;
      scancnt = rpp->repeatscan;
      r1valid = 1;
      olen1 = strlen(rp);
      s1start = 0;
      while (r1valid && ((scancnt == -1) || (scancnt > 0)) &&
               ((ret_hlst = ca_scanstring(rpp,rp,strlen(rp),adfsm)) != NULL))
        {
        r1valid =  ca_readisvalid(rpp,rp,rpp->r1bufs[CA_fqlt_qualln],
                                   ret_hlst,&s1start,&olen1);
        *(rp + olen1) = '\0';
        ca_killhitlst(&hlist1);
        hlist1 = ret_hlst;
        scancnt--;
        }
      hlist2 = ret_hlst = NULL;
      if (r1valid && rpp->pairedends) 
        {
        rp = rpp->r2bufs[CA_fqlt_sqline];
        scancnt = rpp->repeatscan;
        r2valid = 1;
        olen2 = strlen(rp);
        s2start = 0;
        while (r2valid && ((scancnt == -1) || (scancnt > 0)) &&
                ((ret_hlst = ca_scanstring(rpp,rp,strlen(rp),adfsm)) != NULL))
          {
          r2valid = ca_readisvalid(rpp,rp,rpp->r2bufs[CA_fqlt_qualln],
                                    ret_hlst,&s1start,&olen2);
          *(rp + olen2) = '\0';
          ca_killhitlst(&hlist2);
          hlist2 = ret_hlst;
          scancnt--;
          }
        }
      if (r1valid &&
            (!rpp->pairedends || r2valid))
        {
        ca_writetrmreadout(rpp->dstfl1,rpp->r1bufs,rpp,s1start,olen1);
        if (rpp->pairedends)
          ca_writetrmreadout(rpp->dstfl2,rpp->r2bufs,rpp,s2start,olen2);
        }
      break;
    }    
  ca_killhitlst(&hlist1);
  if (rpp->pairedends)
    ca_killhitlst(&hlist2);
  }
/* dispose of buffers as required */
for (fqlt = CA_fqlt_hdrline; fqlt < CA_fqlt_unk; fqlt++)
  {
  memfree(rpp->r1bufs[fqlt]);
  if (rpp->pairedends)
    memfree(rpp->r2bufs[fqlt]);
  }
/* close dst files */
ca_closesrcdststrct(rpp->dstfl1);
if (rpp->pairedends)
  ca_closesrcdststrct(rpp->dstfl2);
}

void ca_lstrawadapt(CA_RUNPARS *rpp,
                    char *hdr)
/* list raw stuff preceded by hdr if
nonNULL */
{
CA_ADAPTOR *ap;

if (hdr != NULL)
  ca_puts(hdr,rpp->dstfl1);
ap = rpp->rawalist;
while (ap != NULL)
  {
#ifndef NO_ZLIB
  if (rpp->dstfl1->comprssd)
    gzprintf(rpp->dstfl1->gzfile,"%d:\"%s\":\"%s\":@%p\n",ap->asqno,ap->adaptstr,
             (ap->adnam!=NULL?(char *)ap->adnam:""),ap);
  else
#endif
    fprintf(rpp->dstfl1->txtfile,"%d:\"%s\":\"%s\":@%p\n",ap->asqno,ap->adaptstr,
             (ap->adnam!=NULL?(char *)ap->adnam:""),ap);
  ap = ap->nxtadaptor;
  }
}

void ca_lstadapts(CA_RUNPARS *rpp,
                  char *hdr)
/* as above for complete adaptor list */
{
CA_ADAPTOR *ap;
CA_ADPTRELT *aep;

if (hdr != NULL)
  ca_puts(hdr,rpp->dstfl1);
ap = rpp->adaptlst;
while (ap != NULL)
  {
#ifndef NO_ZLIB
  if (rpp->dstfl1->comprssd)
    gzprintf(rpp->dstfl1->gzfile,"%d:\"%s\":%p:",ap->asqno,ap->adaptstr,ap->adnam);
  else
#endif
    fprintf(rpp->dstfl1->txtfile,"%d:\"%s\":%p:",ap->asqno,ap->adaptstr,ap->adnam);
  aep = (CA_ADPTRELT *) ap->rlist;
  while (aep != NULL)
    {
#ifndef NO_ZLIB
    if (rpp->dstfl1->comprssd)
      gzprintf(rpp->dstfl1->gzfile,"%d%s",aep->aptr->asqno,(aep->nxtapelt!=NULL?"/":""));
    else
#endif
      fprintf(rpp->dstfl1->txtfile,"%d%s",aep->aptr->asqno,(aep->nxtapelt!=NULL?"/":""));
    aep = aep->nxtapelt;
    }
  ca_putc('\n',rpp->dstfl1);
  ap = ap->nxtadaptor;
  }
}

int main(int argc,
         char *argv[])
{
FILE *adsrcfl;
FS_FSMSTRCT *adfsm;
int adcnt;
int ap;
char op;
char *fendp;
int mismatches;
CA_ADAPTOR *adaptlst;
CA_RUNPARS runpars;
char *defadaptor;

bas_initmalfl("cleanadaptors_malloc");
if (argc <=1)
  {
  say_usage(stdout,argv[0]);
  exit(0);
  }
debuglevel = 1;
adsrcfl = NULL;
adfsm = NULL;
runpars.read1fl = ca_getsrcdststrct("R1_source");
runpars.read1fl->txtfile = stdin;
runpars.read2fl = ca_getsrcdststrct("R2_source");
runpars.dstfl1 = ca_getsrcdststrct("R1_dest");
runpars.dstfl1->txtfile = stdout;
runpars.dstfl2 = ca_getsrcdststrct("R2_dest");
runpars.readflbuflen = DEF_READBUFLEN;
runpars.minleadtrail = MIN_LEAD_TRAIL;
runpars.maxadlen = MAX_ADAPTOR_LEN;
runpars.margn3p = DEF_3PMARGIN;
runpars.mismatches = 0;
runpars.adaptlst = NULL;
runpars.trimcnt = 0;
runpars.ostyle = CA_out_listreads;
runpars.skipres = 0;
runpars.matcriterion = DEF_MATCHCRITERION;
runpars.mintrimlen = 1;
runpars.trimback = 0;
runpars.srcstyle = CA_src_fastq;
runpars.pairedends = 0;
runpars.phredqlimit = 0;
runpars.qualbase = (int) '!';
runpars.p5trunc = 0;
runpars.adflfmt = CA_adflfmt_list;
runpars.maxoutlen = 0;
runpars.rawalist = NULL;
runpars.compress = 0;
runpars.repeatscan = 2;
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 'L':     /* print/list fsm */
        runpars.ostyle = CA_out_fsmlist;
        break;
      case 'a':    /* adaptor max len */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.maxadlen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.maxadlen <= 0)
              err_msg_die("Invalid max adaptor length '%s'\n",argv[ap]);
          }
        break;
      case 'm':    /* min lead/trail */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.minleadtrail = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.minleadtrail <= 0)
              err_msg_die("Invalid min lead/trail length '%s'\n",argv[ap]);
          }
        break;
      case 'M':    /* 3p margin for trim */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.margn3p = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.margn3p <= 0)
              err_msg_die("Invalid 3' trim margin '%s'\n",argv[ap]);
          }
        break;
      case 'R':    /* read buff length */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.readflbuflen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.readflbuflen <= 0)
              err_msg_die("Invalid read file buffer length '%s'\n",argv[ap]);
          }
        break;
      case 'p':    /* match criterion percentage */
        if (++ap > argc)
          err_msg_die("-%c needs float\n",op);
        else
          {
          runpars.matcriterion = strtof(argv[ap],&fendp);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert float '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.matcriterion < 0.0)
              err_msg_die("Invalid match criterion value '%s'\n",argv[ap]);
          }
        break;
      case 's':    /* skip residues each read */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.skipres = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.skipres <= 0)
              err_msg_die("Invalid read skip length '%s'\n",argv[ap]);
          }
        break;
      case 'S':    /* disable single base mismatches */
        runpars.mismatches = 1;
        break;
      case 'i':    /* user name for adaptor file */
      case 'I':    /* fasta format */
        if (++ap > argc)
          err_msg_die("-%c needs adaptor seq file name\n",op);
        else
          if ((adsrcfl = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open adaptor seq file '%s'\n",argv[ap]);
          else
            if (op == 'i')
              runpars.adflfmt = CA_adflfmt_list;
            else
              runpars.adflfmt = CA_adflfmt_fasta;
        break;
      case 'H':    /* report hit matches only */
      case 'F':    /* Trim reads in file */
      case 'f':    /* Report all reads matches & other */
      case 'N':    /* N-fill lines */
        switch (op)
          {
          case 'H':
            runpars.ostyle = CA_out_listhitreads;
            break;
          case 'F':
            runpars.ostyle = CA_out_trimhits;
            break;
          case 'N':
            runpars.ostyle = CA_out_nfill;
            break;
          case 'f':
          default:
            runpars.ostyle = CA_out_listreads;
            break;
          }
        if (++ap > argc)
          err_msg_die("-%c needs read file name\n",op);
        else
          if (strcmp(argv[ap],"-") != 0)
            if (!ca_opensrcdststrct(runpars.read1fl,runpars.compress,argv[ap],"r"))
              err_msg_die("Can't open read file '%s'\n",argv[ap]);
        break;
      case 'x':    /* set minimum trimmed seq length */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.mintrimlen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 't':    /* trimback */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.trimback = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 'T':    /* 5' truncation of adaptor seqs */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.p5trunc = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 'q':    /* Qualitytrim */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.phredqlimit = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.phredqlimit < 0)
              err_msg_die("Invalid phred score '%s'\n",argv[ap]);
          }
        break;
      case 'Q':    /* Quality base */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.qualbase = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if ((runpars.qualbase < 0) || (runpars.qualbase > 127))
              err_msg_die("Invalid base for quality scores '%s' should be >0&<128\n",argv[ap]);
          }
        break;
      case 'G':   /* paired end read trim */
        if (++ap > argc)
          err_msg_die("-%c needs read 2 filename\n",op);
        else
          if (!ca_opensrcdststrct(runpars.read2fl,runpars.compress,argv[ap],"r"))
            err_msg_die("Can't open read 2 file '%s'\n",argv[ap]);
          else
            runpars.pairedends = 1;
        break;
      case 'o':   /* output file read 1 */
        if (++ap > argc)
          err_msg_die("-%c needs read 1 out name\n",op);
        else
          if (!ca_opensrcdststrct(runpars.dstfl1,runpars.compress,argv[ap],"w"))
            err_msg_die("Can't open read 1 out file '%s'\n",argv[ap]);
        break;
      case 'O':   /* output file read 2 */
        if (++ap > argc)
          err_msg_die("-%c needs read 2 out name\n",op);
        else
          if (!ca_opensrcdststrct(runpars.dstfl2,runpars.compress,argv[ap],"w"))
            err_msg_die("Can't open read 2 out file '%s'\n",argv[ap]);
        break;
      case 'l':    /* read max out len */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.maxoutlen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (runpars.maxoutlen <= 0)
              err_msg_die("Invalid max read output length '%s'\n",argv[ap]);
          }
        break;
      case 'd':   /* debug on */
        debuglevel = 1;
        break;
      case 'D':   /* debug^2 on */
        debuglevel = 2;
        break;
      case 'z':   /* compress on */
#ifdef NO_ZLIB
        err_msg_die("Compiled with NO_ZLIB set: compression/decompression unavailable\n");
#else
        runpars.compress = 1;
#endif
        break;
      case 'Z':   /* compress off */
        runpars.compress = 0;
        break;
      case 'r':   /* set repeat scan count */
        if (++ap > argc)
          err_msg_die("-%c needs integer\n",op);
        else
          {
          runpars.repeatscan = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
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
/* sanity checks: if 2nd source file, need 2nd dest file */
#ifdef NO_ZLIB
if ((runpars.pairedends && (runpars.dstfl2->txtfile == NULL)) ||
     (((runpars.dstfl2->txtfile != NULL)) && !runpars.pairedends))
#else
if ((runpars.pairedends && ((runpars.dstfl2->gzfile == NULL) && runpars.dstfl2->txtfile == NULL)) ||
     (((runpars.dstfl2->gzfile != NULL) || (runpars.dstfl2->txtfile != NULL)) && !runpars.pairedends))
#endif
  err_msg_die("-G <read2flin> & -O <read2flout> both needed\n");
/* if we are 5' truncating, then we must trim back by that much */
runpars.trimback += runpars.p5trunc;
if (adsrcfl != NULL)
  {
  adfsm = fs_initnewfsm(4,1,FS_inv_reset);
  adcnt = ca_readadaptorfile(adsrcfl,&runpars,adfsm);
  (void) fs_bldfsm(adfsm,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
  if (runpars.ostyle == CA_out_fsmlist)
    {
    ca_lstrawadapt(&runpars,"#Raw-SqNo.:Seq:Name:PtrAddress\n");
    ca_lstadapts(&runpars,"#Adaptors-ASqNo.:Seq:PtrAddress:RawSqNos_list\n");
/*    fputs("#fsm:data pairs\n",runpars.dstfl1);
    fs_lststrpr(runpars.dstfl1,adfsm->prlst,1,ca_prtlval,ca_nl); */
    ca_puts("#STT:\n",runpars.dstfl1);
    if (runpars.dstfl1->comprssd)
      exit(0);
    else
      {
      fs_lststttbl(runpars.dstfl1->txtfile,adfsm);
      if (runpars.read1fl->txtfile == stdin)  /* don't want to do anything else */
      exit(0);
      }
    }
  }
ca_scanreads(&runpars,adfsm);
exit(0);
}
