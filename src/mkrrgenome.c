/* mkrrgenome: to generate reduced representation genome files for all
or selected chromosomes */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/param.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "mrg_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"
#include "bin_cnts.h"
#include "fsm_ops.h"

/* Local defines */
/* version: first OK version: 1-Apr-2011 */
/* #define PROG_VERSION 1.0 */
/* generate CpG position list: Nov-2011 */
/* #define PROG_VERSION 1.1 */
/* Allow for non-human & non XY chromosomes: Feb-2012 */
/* #define PROG_VERSION 1.2 */
/* #define PROG_VERSION 1.21 */
/* cpg counts too high: Dec-2012 */
/* #define PROG_VERSION 1.22 */
/* allow chromosome info file -G (modified -H): Sep-2014 */
#define PROG_VERSION 1.23
/* additional site changes: Oct-2015 */

#define MRG_DEFMIN 40
#define MRG_DEFMAX 220

typedef enum MRG_outstyle
  {
  MRG_out_rrgenomes = 0,  /* usual mode: generate RR genome files */
  MRG_out_rrlist,         /* tabbed list of valid RR fragment positions */
  MRG_out_cpginrng,       /* tabbed list of CpGs in valid RR fragments */
  MRG_out_allcpg          /* tabbed list of all CpGs */
  }
MRG_OUTSTYLE;

typedef struct DM_chr_info     /* information about & bins for a chromosome */
  {        /* particularly the individual bin lst and the bin list elements */
/*   DM_CPG_BIN *binlst; */
/*   DM_BIN_LSTELT *beltlst; */
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

/* let's make a structure to hold the things we need */

typedef struct MRG_rundata
  {
  int maxchrno;            /* number of separate chromosomes */
/*  char **chromoseq;  */      /* to be allocated an array of seq pointers */
/*  int *chrlens;     */       /* list of their lengths */
  MRG_REGN_ELT **rrregions; /* region list list */
  MRG_OUTSTYLE ostyle;
  int havexy;              /* regard highest 2 chromosomes as X&Y, else numbered */
  FS_FSMSTRCT *cpgfsmp;    /* ptr to CpG searching fsm */
  int maxfrag;             /* size limits for fragments */
  int minfrag;
  char *chrinfoflnam;      /* name of chr info file */
  FILE *chrinfofl;         /* the chr info file handle */
  char *genomsqhdrstr;     /* hdr for numbered (+XY) chromosome files */
  char *destsqhdrstr;      /* hdr for output files */
  DM_CHR_INFO *chrinfo;    /* array of information on chromsomes */
  WRD_LUSTRCT *chridlu;    /* chromosome ID lookup */
  FILE *rstfile;        /* non-NULL, if we have opened a file of restriction sites */
  char *rstflnam;
  DM_SITES *sitelist;
  }
MRG_RUNDATA;

/* global debug & out style variables, for simplicity of access */
RBC_DBG debuglevel;

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
fprintf(fl,"%s (v%.2f): generate reduced representation genome files for MspI (or other) digests\n",
          pnam,PROG_VERSION);
fputs("Options:\n",fl);
fputs("     -g <genomehead> dir and file string to locate genomic seq files by adding n.fa\n",fl);
fputs(" or  -G <chr_info_file> file of chromosome IDs and filenames\n",fl);
fputs("     -k <MaxCNo> set maximum chromosome number (def = ChrY for Humans)\n",fl);
fputs("     -z Don't expect X/Y chromosomes, only numbered names (def=do)\n",fl);
fputs("     -c <C_no> restrict activity to chromosome C_no (1..22,X,Y), def = all\n",fl);
fprintf(fl,"     -M <j,k> scan for restricted rep fragment sizes between j & k residues. (def=%d..%d)\n",
          MRG_DEFMIN,MRG_DEFMAX);
fputs("     -m <j,k>   ditto, produce tab delimited list of positions to stdout\n",fl);
fputs("     -p <j,k> generate tabbed list of CpG positions for j-k fragments\n",fl);
fputs("     -P generate tabbed list of all CpG positions\n",fl);
fputs("     -H <desthead> dir & file string for output genome files (-M), completed with n.fa\n",fl);
fputs("     -n <restrictionfile> - use sites in <restrictionfile>, def=MspI\n",fl);
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

int bc_fragsizeok(int fragmin,
                  int fragmax,
                  int fragsize)
/* return if fragsize fragment is acceptable,
noting that zero limits will always return true. */
{
if ((fragmin != 0) && (fragmax != 0))
  return(rbc_intinrng(fragmin,fragsize,fragmax));
else
  return(1);
}

int rbc_getnscangenomesqs(MRG_RUNDATA *rdp,
                          char *hdrstr,
                          int uchrno)
/* use hdrstr to create a series of file names, one
for each chromosome, open each as a Fasta sequence
file and create a buffer for each.  Return the
number of chromosomes processed */
{
int chno;
int nblen;
SQFL_STRCT *chsqfl;
int rcnt;
char nsc;
int cpos;
FS_RESELT *frp;
FS_DATPRELT *dpep;
int prevmat;
MRG_REGN_ELT *lstend;
int thismat;
int cpgcnt;
MRG_REGN_ELT *cpglstp;
MRG_REGN_ELT *cpglstend;
MRG_REGN_ELT *cp;

rcnt = 0;
for (chno = 0; chno < rdp->maxchrno; chno++)
  if ((uchrno == 0) || (uchrno == chno+1))
    {
    if ((chsqfl = sqfl_opnsqstrct((rdp->chrinfo+chno)->chrflname,
                                   SFMT_fasta,"r")) != NULL)
      {
      (rdp->chrinfo+chno)->chrlen = readsrcsq(chsqfl,NULL);
      if (debuglevel > RBC_dbg_on)
        {
        fprintf(stdout,"%s: %d res,",(rdp->chrinfo+chno)->chrflname,(rdp->chrinfo+chno)->chrlen);
        fflush(stdout);
        }
      if (rdp->ostyle == MRG_out_rrgenomes)
        (rdp->chrinfo+chno)->chromseq =
           (char *) getmemory((rdp->chrinfo+chno)->chrlen+1,"Chr Buff");
      sqfl_rewind(chsqfl);
      cpos = 0;
      prevmat = 1;
      lstend = NULL;
      fs_initrun(rdp->cpgfsmp);
      cpgcnt = 0;
      (void) sqfl_skipsqflhdr(chsqfl);
      cpglstp = cpglstend = NULL;
      while ((nsc = sqfl_getnxtres(chsqfl)) != '\0')
        {
        if ((frp = fs_procchr(rdp->cpgfsmp,nsc,tr_nares2int)) != NULL)   /* have a match */
          {
          dpep = (FS_DATPRELT *) frp->action;
          switch (dpep->ldstr)
            {
            case 2:        /* CpG */
              cpgcnt++;
              switch (rdp->ostyle)
                {
                case MRG_out_allcpg:
                  fprintf(stdout,"%s\t%d\n",(rdp->chrinfo+chno)->chrid,cpos);
                  break;
                case MRG_out_cpginrng:
                  cpglstend = mrg_appndrgnelt(&cpglstend,cpos,cpos,1);
                  if (cpglstp == NULL)
                    cpglstp = cpglstend;
                  break;
                default:
                  break;
                }
              break;
            case 4:       /* CCGG */
              thismat = cpos - dpep->ldstr + ((int) dpep->stxt) + 1;
/*          if (debuglevel > RBC_dbg_none)
            fprintf(stdout,"%s@%d, Frag=%d\n",dpep->dstrng,thismat,thismat-prevmat+1); */
              if (bc_fragsizeok(rdp->minfrag,rdp->maxfrag,thismat - prevmat + 1))
                {
                switch (rdp->ostyle)
                  {
                  case MRG_out_cpginrng:
                    cp = cpglstp;
                    while (cp != NULL)
                      {
                      fprintf(stdout,"%s\t%d\n",(rdp->chrinfo+chno)->chrid,cp->rstart);
                      cp = cp->nxtregn;
                      }
                    break;
                  case MRG_out_rrgenomes:
                  case MRG_out_rrlist:
                    lstend = mrg_appndrgnelt(&lstend,prevmat,thismat,cpgcnt);
                    if (*(rdp->rrregions+chno) == NULL)
                      *(rdp->rrregions+chno) = lstend;
                    break;
                  default:
                    break;
                  }
                }
              mrg_clrallregnelts(&cpglstp);
              cpglstp = cpglstend = NULL;
              prevmat = thismat + 1;
              cpgcnt = 0;       /* is a cpg in CCGG */
              break;
            }
          }
        if (rdp->ostyle == MRG_out_rrgenomes)
          {
          *((rdp->chrinfo+chno)->chromseq + cpos) = nsc;
          *((rdp->chrinfo+chno)->chromseq + cpos + 1) = '\0';
          }
        cpos++;
        }
      if (bc_fragsizeok(rdp->minfrag,rdp->maxfrag,cpos - prevmat + 1))
        (void) mrg_appndrgnelt(&lstend,prevmat,cpos,cpgcnt);
      sqfl_clssqstrct(chsqfl);
      if ((debuglevel > RBC_dbg_on) && (rdp->ostyle == MRG_out_rrgenomes))
        {
        fprintf(stdout,"'%.10s...'\n",(rdp->chrinfo+chno)->chromseq);
        fflush(stdout);
        }
      rcnt++;
      mrg_clrallregnelts(&cpglstp);
      }
    else
      {
      err_msg("Can't open chromosome file %s\n",(rdp->chrinfo+chno)->chrflname);
      (rdp->chrinfo+chno)->chrlen = 0;
      }
    }
return(rcnt);
}

void mrg_putreducedrepseq(MRG_RUNDATA *rdp,
                          char *hdrstr,
                          char *seqbuf,
                          int chrlen,
                          MRG_REGN_ELT *rrfragp,
                          int chno)
/* use hdrstr to create a file name for thisp
chromosome, open as a Fasta output sequence
file and write regions in rrfragp list there to. */
{
char *sqfnam;
int nblen;
SQFL_STRCT *chsqfl;
int cpos;
MRG_REGN_ELT *rrfp;
int rc;
int lcnt;

sqfnam = (char *) getmemory((nblen = strlen(hdrstr) + 16),"Sq dest file name buf");
snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,(rdp->chrinfo+chno)->chrid);
if ((chsqfl = sqfl_opnsqstrct(sqfnam,SFMT_fasta,"w")) != NULL)
  {
  rrfp = rrfragp;
  fprintf(chsqfl->sfl,">%srr reduced repr %d..%d for MspI digest Chr%s %d/%dbp CpG: %d\n",
            (rdp->chrinfo+chno)->chrid,rdp->minfrag,rdp->maxfrag, 
            (rdp->chrinfo+chno)->chrid,mrg_sumregnelts(rrfragp),chrlen,
            mrg_sumcpgregnelts(rrfragp));
  rc = lcnt = 0;
  while (rrfp != NULL)
    {
    cpos = rrfp->rstart;
    while (cpos <= rrfp->rstop)
      {
      sqfl_putres(chsqfl,*(seqbuf + cpos - 1),rc,&lcnt);
      cpos++;
      }
    rrfp = rrfp->nxtregn;
    }
  sqfl_termsqfl(chsqfl,"",&lcnt);
  sqfl_clssqstrct(chsqfl);
  }
memfree(sqfnam);
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

int mrg_str2chrno(MRG_RUNDATA *rdp,
                  char *ustr)
/* see if ustr matches any thing on rdp-dependent
list.  Return that value, else 0 =(CHR_unk) */
{
return(any_str2chrno(rdp->maxchrno,rdp->havexy,ustr));
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
                        MRG_RUNDATA *rpp)
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

DM_SITES *dm_appndsite(DM_SITES **sitlst,
                       char *site,
                       int cutoffset)
/* dup site and append an element for it to the end of sitlst.
return address of new element in case useful */
{
DM_SITES *prev, *end_ptr;

if (sitlst != NULL)
  {
  prev = end_ptr = *sitlst;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtsite;
    }
  end_ptr = (DM_SITES *) getmemory(sizeof(DM_SITES),"Site_elt");
  end_ptr->nxtsite = NULL;
  end_ptr->sitestr = bas_strdup(site);
  end_ptr->offset = cutoffset;
  if (*sitlst == NULL)
    *sitlst = end_ptr;
  else
    prev->nxtsite = end_ptr;
  return(end_ptr);
  }
else
  return(NULL);
}

DM_SITES *dm_appndsitestr(DM_SITES **sitlst,
                          char *recogstr)
/* recogstr is a restriction recognition site
with (optional) cleavage position indicated with
'^'.  append the site/offset derived from this to
*sitlst */
{
char *scopy;
char *sp;
char *dp;
int offst;

sp = dp = scopy = bas_strdup(recogstr);
offst = 0;
while (*sp != '\0')
  {
  if (*sp == '^')
    offst = (int) (sp - scopy);
  else
    {
    *dp = *sp;
    dp++;
    }
  sp++;
  }
*dp = '\0';
return(dm_appndsite(sitlst,scopy,offst));
}

int dm_readrstsites(FILE *rfl,
                    DM_SITES **rstlst)
/* read successive lines from rfl, check
length > 0 and process as restriction sites.
return number of sites found */
{
char *lbuf;
char *lp;
int scnt;

lbuf = (char *) getmemory(FILENAME_MAX + 1,"site buff");
scnt = 0;
while ((lp = fgets(lbuf,FILENAME_MAX,rfl)) != NULL)
  if (strlen(lp) > 0)
    {
    while (*lp != '\0')
      {
      if (*lp == '\n')
        *lp = '\0';
      lp++;
      }
    (void) dm_appndsitestr(rstlst,lbuf);
    scnt++;
    }
memfree(lbuf);
return(scnt);
}

int main(int argc,
         char *argv[])
{
int ap;
char op;
SQFL_STRCT *srcsq;
int ecnts;
int modcnt;
int chrcnt;
int uchrno;
RBC_CHRNO chrp;
MRG_REGN_ELT *rlstp;
int maxchrno;
MRG_RUNDATA rdata;
char *fendp;
DM_SITES *sitp;

debuglevel = RBC_dbg_none;
modcnt = ecnts = 0;
uchrno = 0;
rdata.maxchrno = ChrY;
rdata.havexy = 1;
rdata.ostyle = MRG_out_rrgenomes;
rdata.cpgfsmp = NULL;
rdata.minfrag = MRG_DEFMIN;
rdata.maxfrag = MRG_DEFMAX;
rdata.genomsqhdrstr = NULL;
rdata.destsqhdrstr = NULL;
rdata.chrinfoflnam = NULL;
rdata.chrinfofl = NULL;
rdata.chridlu = NULL;
rdata.rstfile = NULL;
rdata.rstflnam = NULL;
rdata.sitelist = NULL;
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 'c':   /* a chromosome */
        if (++ap > argc)
          err_msg_die("-%c needs a chromosome identifier (1..20,X,Y)\n",op);
        else
          if ((uchrno = mrg_str2chrno(&rdata,argv[ap])) == Chr_unk)
            err_msg_die("Can't determine Chromosome '%s'\n",argv[ap]);
        break;
      case 'g':
        if (++ap > argc)
          err_msg_die("-%c needs genome files header string\n",op);
        else
          rdata.genomsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'G':   /* file name for chromosome spec file */
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          {
          rdata.chrinfoflnam = bas_strdup(argv[ap]);
          if ((rdata.chrinfofl = fopen(rdata.chrinfoflnam,"r")) == NULL)
            err_msg_die("Can't open chromosome information file '%s' for -%c\n",
                          rdata.chrinfoflnam,op);
          }
        break;
      case 'H':
        if (++ap > argc)
          err_msg_die("-%c needs header string\n",op);
        else
          rdata.destsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'M':
      case 'm':
      case 'p':
        switch (op)
          {
          case 'p':
            rdata.ostyle = MRG_out_cpginrng;
            break;
          case 'm':
            rdata.ostyle = MRG_out_rrlist;
            break;
          case 'M':
          default:
            rdata.ostyle = MRG_out_rrgenomes;
            break;
          }
        if (++ap < argc)
          if (*argv[ap] != '-')
            bc_commalst2ints(argv[ap],&rdata.minfrag,&rdata.maxfrag);
          else
            ap--;
        break;
      case 'P':   /* all CpG list */
        rdata.ostyle = MRG_out_allcpg;
        break;
      case 'n':   /* -n restriction site file */
        if (++ap > argc)
          err_msg_die("-%c needs nuclease site file name\n",op);
        else
          {
          rdata.rstflnam = bas_strdup(argv[ap]);
          if ((rdata.rstfile = fopen(rdata.rstflnam,"r")) == NULL)
            err_msg_die("Can't open nuclease site file '%s' for -%c\n",
                          rdata.rstflnam,op);
          }
        break;
      case 'd':   /* debug on */
        debuglevel = RBC_dbg_on;
        break;
      case 'D':   /* debug^2 on */
        debuglevel = RBC_dbg_serious;
        break;
      case 'z':   /* disable XY chromosome naming */
        rdata.havexy = 0;
        break;
      case 'k':   /* set max chromosome number */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          {
          rdata.maxchrno = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (rdata.maxchrno < 0)
              err_msg_die("Invalid chromosome number '%s'\n",argv[ap]);
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
/* set up fsm as required: */
switch (rdata.ostyle)
  {
  case MRG_out_cpginrng:
  case MRG_out_rrlist:
  case MRG_out_rrgenomes:
    rdata.cpgfsmp = fs_initnewfsm(4,1,FS_inv_ignor);
    if (rdata.rstfile != NULL)
      {
      (void) dm_readrstsites(rdata.rstfile,&rdata.sitelist);
      fclose(rdata.rstfile);
      }
    else
      dm_appndsitestr(&rdata.sitelist,"C^CGG");
    sitp = rdata.sitelist;
    while (sitp != NULL)
      {
      (void) fs_adddatprs(rdata.cpgfsmp,sitp->sitestr,(void *) sitp->offset);
      sitp = sitp->nxtsite;
      }
    fs_adddatprs(rdata.cpgfsmp,"CG",NULL);
    (void) fs_bldfsm(rdata.cpgfsmp,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
    break;
  case MRG_out_allcpg:
    rdata.cpgfsmp = fs_initnewfsm(4,1,FS_inv_ignor);
    fs_adddatprs(rdata.cpgfsmp,"CG","CG");
    (void) fs_bldfsm(rdata.cpgfsmp,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
    break;
  default:
    break;
  }    
if (rdata.chrinfofl != NULL)
  rdata.maxchrno = dm_scanchrinfofl(rdata.chrinfofl,NULL,NULL);
/* create chromosome information array */
rdata.chrinfo = (DM_CHR_INFO *) getmemory(sizeof(DM_CHR_INFO)*rdata.maxchrno,
                                         "Chr bin info");
/* chr ID lookup stuff */
rdata.chridlu = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"ChrIDlookup");
wlu_initlustrct(rdata.chridlu,WLU_CASEIND,-1);
rdata.rrregions = (MRG_REGN_ELT **) getmemory(rdata.maxchrno*sizeof(MRG_REGN_ELT *),
                     "Chr Region Elts");
/* init info array */
for (chrp = 0; chrp < rdata.maxchrno; chrp++)
  {
  (rdata.chrinfo+chrp)->chromseq = NULL;
  (rdata.chrinfo+chrp)->chrlen = 0;
  (rdata.chrinfo+chrp)->chrflname = NULL;
  (rdata.chrinfo+chrp)->chrid = NULL;
  *(rdata.rrregions+chrp) = NULL;
  }
/* establish file names, etc. either from params or the info file */
if ((rdata.chrinfofl == NULL) && (rdata.genomsqhdrstr != NULL))
  dm_bldchrfilenames(rdata.genomsqhdrstr,&rdata);
else
  {
  rewind(rdata.chrinfofl);
  (void) dm_scanchrinfofl(rdata.chrinfofl,rdata.chrinfo,rdata.chridlu);
  fclose(rdata.chrinfofl);
  }
if (rdata.maxchrno > 0)
  {
  chrcnt = rbc_getnscangenomesqs(&rdata,rdata.genomsqhdrstr,uchrno);
  if (debuglevel > RBC_dbg_none)
    fprintf(stdout,"%d chromosomes read\n",chrcnt);
  for (chrp = 0; chrp < rdata.maxchrno; chrp++)
    if ((uchrno == 0) || (uchrno == chrp+1))
      {
      rlstp = *(rdata.rrregions+chrp);
      switch (rdata.ostyle)
        {
        case MRG_out_rrlist:
          while (rlstp != NULL)
            {
            fprintf(stdout,"%s\t%d..%d (%d bp) CpG: %d\n",
                      (rdata.chrinfo+chrp)->chrid,
                      rlstp->rstart,
                      rlstp->rstop,(rlstp->rstop - rlstp->rstart + 1),rlstp->cpgcnt);
            rlstp = rlstp->nxtregn;
            }
          break;
        case MRG_out_rrgenomes:
          if (rdata.destsqhdrstr == NULL)
            (void) err_msg_die("Need -H <destinationheader>\n");
          else
            mrg_putreducedrepseq(&rdata,rdata.destsqhdrstr,(rdata.chrinfo+chrp)->chromseq,
                                   (rdata.chrinfo+chrp)->chrlen,rlstp,chrp);
          break;
        default:
          break;
        }
      }
  }
else
  {
  if (rdata.genomsqhdrstr == NULL)
    (void) err_msg("No genome file header string (-g option) used\n");
  if (rdata.cpgfsmp == NULL)
    (void) err_msg("No fragment size range specified (-M or -m options)\n");
  exit(1);
  }
exit(0);
}
