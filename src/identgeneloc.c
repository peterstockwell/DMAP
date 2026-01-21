/* identgeneloc.c: load up Genbank/EMBL/gff3/gtf feature table info.
  identifiy proximal genes for a table of:
ChrNo. start stop... */

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
/* #include <strings.h> */
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <sys/param.h>
/* #include <gdbm.h> */

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqtrans.h"
#include "cod_fns.h"
#include "sqmat_fns.h"
#include "dbpars.h"
/* #include "dblu_fns.h" */
#include "mrg_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"

/* local defines */

#define PROG_NAME "identgeneloc"
/* #define PROG_VERS 0.00 */
/* #define PROG_VERS 0.01 */
/* require associated CDS option (-P): 4-Apr-2012 */
/* #define PROG_VERS 0.02 */
/* allow for CpGI & TSS tags (seqmonk), note features starting in RR frags: May-2012 */
/* #define PROG_VERS 0.03 */
/* allow fragments internal to genes, classify wrt exon boundaries: Jun-2012 */
/* #define PROG_VERS 0.04 */
/* Frag-based CpGI & TSS distances: Nov-2012 */
/* #define PROG_VERS 0.05 */
/* relate frag/gene & CpGI positions: Mar-2013 */
/* #define PROG_VERS 0.06 */
/* correct CpGI classification error: Sep-2013 */
/* #define PROG_VERS 0.07 */
/* non directional CpG (-u) & -G chromosome info file: Jul/Aug-2014 */
/* #define PROG_VERS 0.08 */
/* count genes overlapping fragments/regions: Dec-2014 */
/* #define PROG_VERS 0.09 */
/* correct logic for -p & -s options and -c: Mar-2015 */
/* #define PROG_VERS 0.10 */
/* Modify GFF3 code to pick up Amel gene names: Jun-2015 */
/* #define PROG_VERS 0.11 */
/* accept variant GTF files: Jul-2015 */
/* #define PROG_VERS 0.12 */
/* updated dbpars functions for SeqMonk GRCh38-type ID lines: Oct-2016 */
/* #define PROG_VERS 0.13 */
/* -q option for CDS/protein_id name qualifiers for Genbank (esp) files: Apr-2017 */
/* #define PROG_VERS 0.14 */
/* more freedom for defining required gff/gtf attribute: Apr-2017 */
/* #define PROG_VERS 0.15 */
/* include 'gene' attribute values for GenBank, etc.: May-2017 */
/* #define PROG_VERS 0.16 */
/* sanity checks for ft file parameters: May-2019 */
/* #define PROG_VERS 0.17 */
/* optional output for null hitters: -W : Aug-2019 */
/* #define PROG_VERS 0.18 */
/* maintain consistent columns for cpgi NULLs: Aug-2019 */
/* #define PROG_VERS 0.19 */
/* check for enough tokens (file corruption): Oct-2019 */
/* #define PROG_VERS 0.20 */
/* fix .gtf logic: Jun-2020 */
/* #define PROG_VERS 0.21 */
/* increase GTF token limit: Jul-2020 */
/* #define PROG_VERS 0.22 */
/* allow multiple GTF attribute outputs: Aug-2020 */
/* #define PROG_VERS 0.23 */
/* allow multiple user feature types: Mar-2021 */
/* #define PROG_VERS 0.24 */
/* allow Genbank gene_ids to be found + variant Genbank LOCUS lines: Feb-2022 */
/* #define PROG_VERS 0.25 */
/* allow biotypes for GTF feature data: Jun-2023 */
/* #define PROG_VERS 0.26 */
/* append biotype to lines for SeqMonk and GTF annotations: Jul-2023 */
/* #define PROG_VERS 0.27 */
/* -o output file option: May-2024 */
/* #define PROG_VERS 0.28 */
/* correct header line field count anomaly for header lines: 14-Jan-2025 */
/* #define PROG_VERS 0.29 */
/* tidy previous modification: 15-Jan-2025 */
/* #define PROG_VERS 0.30 */
/* actually implement -b input buffer length option: 25-Sep-2025 */
/* #define PROG_VERS 0.31 */
/* allow extra length for input line buffer: 7-Nov-2025 */
/* #define PROG_VERS 0.32 */
/* correct gtf/gff3 srcbuflen issue: 6-Jan-2026 */
/* #define PROG_VERS 0.33 */
/* include dbpars updates for T2T gtf variations: 9-Jan-2026 */
#define PROG_VERS 0.34
/* further dbpars refinement for T2T gtf: 21-Jan-2026 */

#define DEF_SRCBUFLEN 2048

typedef enum IGL_runstyle
  {
  IGL_findproximal = 0,
  IGL_listfeatdata,
  IGL_listchrids
  }
IGL_RUNSTYLE;

typedef enum IGL_intrnl_strtgy   /* which strategy for internal fragments */
  {
  IGL_intrnl_none = 0,           /* don't permit fragments internal to genes */
  IGL_intrnl_allow,              /* do     ditto */
  IGL_intrnl_resolv              /* ditto, further resolve exon/intron locations */
  }
IGL_INTRNL_STRTGY;

typedef enum IGL_outmode    /* part of output style control */
  {
  IGL_omod_genebased = 0,
  IGL_omod_fragbased,
  IGL_omod_genecnt
  }
IGL_OUTMODE;

typedef enum IGL_CpG_dir   /* direction of CpG Is search */
  {
  IGL_CpG_none = 0,
  IGL_CpG_upstream,        /* upstream of gene */
  IGL_CpG_proximal         /* proximal to fragment */
  }
IGL_CpG_DIR;

typedef struct IGL_chr_info     /* information about & bins for a chromosome */
  {        /* particularly the individual bin lst and the bin list elements */
  char *chrflname;             /* file name for fa sequencs */
  char *chrid;                 /* id for this chr */
  }
IGL_CHR_INFO;

typedef struct IGL_runpars  /* runtime parameters/data */
  {
  char *ftflnam;            /* name of feature table file */
  char *ftprefix;           /* prefix for ft files */
  char *ftsuffix;           /* suffix for ft files */
  char *srcflnam;           /* name of source chr_pos_pos file */
  DBLU_DBFMT dfmt;          /* this db format */
  int srclno;               /* count source lines */
  int uchrno;               /* user-defined chromosome No, 0=>implies all */
  DBP_MULTI_ELEMNT *meplist;   /* list of chromosomal feature information */
  DB_FTYPELT *ftlist;       /* list of features wanted */
  DB_FTYPELT *tsslist;      /* tss feature special list */
  DB_FTYPELT *cpgilist;     /* CpGIs feature special list */
  WRD_LUSTRCT *ftkw;        /* useful keywords sfor db format */
  int maxchrno;             /* max different chr, including X,Y */
  int havexy;               /* note if X,Y or not */
  char *chrinfoflnam;       /* name of file of chromosomal information */
  FILE *chrinfofl;          /* file for it */
  WRD_LUSTRCT *chridlu;     /* lookup table for chromosome ids */
  IGL_CHR_INFO *chrinfo;    /* array of chromosomal information */
  IGL_RUNSTYLE rstyle;      /* control listing/vs proximal */
  int srcbuflen;            /* modifiable source buffer length */
  char *chrtokdelmtr;       /* string of token separators for chr pos file */
  char *outtokdelmtr;       /* what separates output elements */
  int srccollmt;            /* columns of chr/pos file to scan & list */
  int needcds;              /* 1 if we require genes to have a corresponding CDS */
  int classolap;            /* 1 to show feature relative to fragment */
  int seektss;              /* 1 if want first upstream tss (SeqMonk data) shown */
  IGL_CpG_DIR seekCpGis;    /* control CpG Island scan (SeqMonk data) */
  int distlimit;            /* limit how far out we look: 0=> nolimit */
  int displinfo;            /* 1 to list other gene info after name, def to 0 */
  int tssdistlmt;           /* limit how far from fragment we will look for TSS; 0=> no limit */
  IGL_INTRNL_STRTGY intrnlfrag; /* allow fragments internal to genes & classify wrt introns/exons */
  int geneexmrna;           /* use mRNA feature to define gene names/positions */
  DB_STR_ELT *wantedcodingtypes; /* list to restrict genes to /biotype=<type> entries, NULL to ignore */
  WRD_LUSTRCT *wantedtypewlu; /* from wantedcodingtypes list */
  int tss_cpgi_rng;         /* 1 if we are to have range positions for these entities */
  int relatefragcpgi;       /* 1 => show fragment relative to CpG Is */
  WRD_LUSTRCT *ftqul;       /* feature table qualifier look up */
  IGL_OUTMODE omode;        /* out style: normal/fragbased */
  int cds_prot_id;          /* 1 if gene identifier is from CDS/protein_id field */
  DB_STR_ELT *gtfgeneidattr;      /* if non-null, these are the wanted gtf/gff attributes */
  DBLU_FTKW ufeattype;      /* if we need a user-specified feature */
  DB_STR_ELT *ufeatlist;    /* list of values for ftlist */
  int printnullhits;        /* print lines for non-hitters (-W) */
  FILE *outfile;            /* user output file, defaults to stdout */
  }
IGL_RUNPARS;

typedef enum IGL_fraggen_loc   /* resolve where a region/fragment lies wrt gene info */
  {
  IGL_frgen_unknown = 0,
  IGL_frgen_upstream,
  IGL_frgen_on5p,
  IGL_frgen_inexon,
  IGL_frgen_inintron,
  IGL_frgen_onexonintron,
  IGL_frgen_onintronexon,
  IGL_frgen_includsexon,
  IGL_frgen_includsintron,
  IGL_frgen_on3p,
  IGL_frgen_genein
  }
IGL_FRAGGENE_LOC;

typedef enum IGL_cpgi_loc  /* relate position to CpGIs range */
  {
  IGL_cpgi_unknown = 0,    /* unknown or not on island,shore,shelf */
  IGL_cpgi_oncore,
  IGL_cpgi_coreshore,
  IGL_cpgi_onshore,
  IGL_cpgi_shoreshelf,
  IGL_cpgi_onshelf,
  IGL_cpgi_shelfolap
  }
IGL_CPGI_LOC;

/* typedef enum IGL_ */

char *oricpy;

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

void err_msg(char *fmt,
             ...)
/* write user error message but continue */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
fflush(stderr);
}

void igl_sayusage(FILE *fl,
                  char *pnam)
{
fprintf(fl,"%s v%.2f: proximal genes for chromosomal position table for a feature table\n",
          pnam,PROG_VERS);
fputs("Usage: -e|g|E|F|T|Q {-f <ftfile> | -p <ftprefix> [-s <ftsuffix>] | -G <infofile>} -r <chr_pos_file>\n",fl);
fputs("where: -e: EMBL feature format (New)\n",fl);
fputs("       -g: Genbank feature format\n",fl);
fputs("       -F: GFF3 feature format\n",fl);
fputs("       -T: GTF feature format\n",fl);
fputs("       -Q: SeqMonk feature format (for CpGIs & TSS)\n",fl);
fputs("       -q: take gene name from protein_id qualifier for CDS (esp. for Genbank, etc. def=gene)\n",fl);
fputs("       -f <ftfile> read feature table from <ftfile>\n",fl);
fputs("       -p <ftprefix> generate feature files from <ftprefix>+chrno+<ftsuffix>\n",fl);
fputs("       -s <ftsuffix> use <ftsuffix> for generated names (def=empty string)\n",fl); 
fputs("       -G <infofile> read chromosome ID & feature file info from <infofile>\n",fl);
fputs("       -r positions from chr_pos_file as Chr<ht>start<ht>stop...: file - => stdin\n",fl);
fputs("       -k <maxchrno>: allow <maxchrno> different chromosomes (including X,Y: def=human)\n",fl);
fputs("       -z: don't use X,Y chromosome names for two highest\n",fl);
fputs("       -c <userchrno>: confine activity to chromosome <userchrno>\n",fl);
fputs("       -l list feature data, don't process chromosome positions (def=do)\n",fl);
fputs("       -L list chromosome ids only (def=find proximal genes)\n",fl);
fputs("       -C <srccollmt> include <srccollmt> columns of chr/pos file (def=9,min=3)\n",fl);
fputs("       -P require corresponding Protein CDS for any gene (def=don't)\n",fl);
fputs("       -B <type> restrict to /biotype=\"<type>\" gene/mRNA features (SeqMonk features), multiples OK, (def=don't)\n",fl);
fputs("       -A <attribute> choose these attributes for gtf/gff3 files, multiples OK\n",fl);
fputs("       -a <FeatureType> specify the feature type: override defaults, multiples OK\n",fl);
fputs("       -m use mRNA to define gene/position (def=use 'gene')\n",fl);
fputs("       -I allow fragments internal to genes (def=don't)\n",fl);
fputs("       -i as for -I, but identify fragment relationships with exon boundaries\n",fl);
fputs("       -K show fragment vs feature info (def=don't)\n",fl);
fputs("       -t show closest upstream TSS (SeqMonk feature list) (def=don't)\n",fl);
fputs("       -U show closest upstream CpG Is (SeqMonk feature list) (def=don't)\n",fl);
fputs("       -u show closest CpG Is (SeqMonk feature list) in any direction (def=don't)\n",fl);
fputs("       -d <limit> don't look beyond <limit> for features (def=nolimit)\n",fl);
fputs("       -D <tsslimit> don't look for TSS beyond this distance(def=nolimit)\n",fl);
fputs("       -R show ranges for TSS & CpGi positions\n",fl);
fputs("       -X append Xref info after gene name (def=don't)\n",fl);
fputs("       -n show nearest TSS/CpGI to fragment rather than gene\n",fl);
fputs("       -N count qualifying genes in fragment range\n",fl);
fprintf(fl,"       -b <buflen> allow input lines up to <buflen> (def=%d)\n",DEF_SRCBUFLEN);
fputs("       -W print lines for non-hitters (def=don't)\n",fl);
fputs("       -o <outfile> - write output to <outfile> (def=stdout)\n",fl);
}

void tg_putgfffeatrecurs(FILE *fl,
                         DB_FEATSTRCT *fp,
                         int indent)
/* recursively write fp and parents to fl, increasing
indent each level */
{
int icnt;
DB_SEGELT *sp;
DB_TBLINFO *tbip;

if (fp != NULL)
  {
  icnt = indent;
  while (icnt-- > 0)
    fputc(' ',fl);
  fprintf(fl,"%d %s %s %s",fp->featno,
            db_ftkw2str(fp->featur),fp->idp,
            (fp->strctsens==DB_sens_norm?"+":"-"));
  sp = fp->fstseg;
  while (sp != NULL)
    {
    fprintf(fl," %d-%d",sp->sgstart,sp->sgstop);
    sp = sp->nextseg;
    }
  if (((tbip = db_tblent4udata(fp->infolist,FTQU_db_xref,NULL)) != NULL) &&
       (tbip->qval != NULL))
    fprintf(fl," DBXref=%s",tbip->qval);
  if (((tbip = db_tblent4udata(fp->infolist,FTQU_GO_info,NULL)) != NULL) &&
       (tbip->qval != NULL))
    fprintf(fl," GO=%s",tbip->qval);
  if (((tbip = db_tblent4udata(fp->infolist,FTQU_note,NULL)) != NULL) &&
       (tbip->qval != NULL))
    fprintf(fl," Note=%s",tbip->qval);
  if (((tbip = db_tblent4udata(fp->infolist,FTQU_gene,NULL)) != NULL) &&
       (tbip->qval != NULL))
    fprintf(fl," gene=%s",tbip->qval);
  if (((tbip = db_tblent4udata(fp->infolist,FTQU_biotype,NULL)) != NULL) &&
       (tbip->qval != NULL))
    fprintf(fl," biotype=%s",tbip->qval);
  fputc('\n',fl);
  tg_putgfffeatrecurs(fl,fp->parentfeat,indent+2);
  }
}

void igl_cleansrcline(char *ln)
/* scan ln for \n and replace with \0 */
{
char *lp;

lp = ln;
while ((*lp != '\0') && (*lp != '\n'))
  lp++;
if (*lp == '\n')
  *lp = '\0';
}

int igl_parseftsource(FILE *src,
                      IGL_RUNPARS *rpp)
/* depending on ft format, parse src and build
lists of features for chromosomes. return
number of features scanned OK */
{
DBP_MULTI_ELEMNT *curseq;
char *sline;
int ecnt;
DB_ENTSTRCT *estpt;
char *trnam;
int ccnt;

ecnt = 0;
switch (rpp->dfmt)
  {
  case DBFMT_gtf:
  case DBFMT_gff3:
    sline = (char *) getmemory(rpp->srcbuflen+8,"scanline");
    curseq = NULL;
    while (fgets(sline,rpp->srcbuflen+7,src) != NULL)
      {
      igl_cleansrcline(sline);
      switch (rpp->dfmt)
        {
        case DBFMT_gtf:
          ecnt += dbp_parse_gtf_ln(src,&curseq,sline,rpp->ftkw,rpp->ftlist,rpp->wantedtypewlu,&rpp->meplist,rpp->gtfgeneidattr);
          break;
        case DBFMT_gff3:
        default:
          ecnt += dbp_parse_gff3_ln(src,&curseq,sline,rpp->ftkw,rpp->ftlist,&rpp->meplist);
          break;
        }
      }
    memfree(sline);
    switch (rpp->dfmt)
      {
      case DBFMT_gtf:
        wlu_clrlustrct(rpp->chridlu);
        wlu_initlustrct(rpp->chridlu,WLU_CASEIND,-1);
        (void) dbp_bldgtfgenes(rpp->meplist);
        curseq = rpp->meplist;
        ccnt = 1;
        while (curseq != NULL)
          {
          wlu_addwrd(rpp->chridlu,curseq->me_id,ccnt,NULL);
          ccnt++;
          curseq = curseq->nxt_dme;
          }
        break;
      case DBFMT_gff3:
        curseq = rpp->meplist;
        while (curseq != NULL)
          {
          (void) dbp_findparentfeats(curseq->me_entstrptr,strcmp,0);
          curseq = curseq->nxt_dme;
          }
        break;
      default:
        break;
      }
    break;
  DBFMT_embl:
  DBFMT_genbank:
  DBFMT_emblold:
  DBFMT_sqmonk:
  default:
    while ((estpt = db_parsfl4feats(src,rpp->dfmt,rpp->ftlist)) != NULL)
      {
      if ((trnam = rindex(estpt->ename,':')) == NULL)
        trnam = estpt->ename;
      else
        trnam++;
      curseq = dbp_appnd_mult_elemnt(&rpp->meplist,trnam,estpt);
      ecnt += db_countfeats(estpt,FTKW_unknown);
      }
    break;
  }
return(ecnt);
}

int igl_parseftflnam(IGL_RUNPARS *rpp,
                     char *flfnam)
  /* using information in rpp, try to open a ft source
file flfnam and read contents.  Return the consequences */
{
FILE *src;
int ecnt;

src = NULL;
ecnt = 0;
if ((flfnam != NULL) && ((src = fopen(flfnam,"r")) != NULL))
  {
  ecnt = igl_parseftsource(src,rpp);
  fclose(src);
  }
else
  if (flfnam == NULL)
    err_msg_die("No ft source file name\n");
  else
    err_msg_die("Can't open ft source file '%s'\n",flfnam);
return(ecnt);
}

int igl_sensupstreamdist(DB_FEATSTRCT *fp,
                         DB_SENS sens,
                         int r5ppos,
                         int r3ppos)
/* return the distance from start of region r5ppos..r3ppos to
start of fp in sense sens */
{
if (sens == DB_sens_comp)
  return(db_5pextnt4feat(fp) - r3ppos + 1);
else
  return(r5ppos - db_3pextnt4feat(fp) + 1);
}

int igl_sens3pupstreamdist(DB_FEATSTRCT *fp,
                           DB_SENS sens,
                           int r5ppos,
                           int r3ppos)
/* return the distance from the 3' end of
region to start of fp, depending on sense
of fp.  Assume fp is non-NULL */
{
if (sens == DB_sens_comp)
  return(r5ppos - db_3pextnt4feat(fp) + 1);
else
  return(db_5pextnt4feat(fp) - r3ppos + 1);
}

int igl_sens5pupstreamdist(DB_FEATSTRCT *fp,
                           DB_SENS sens,
                           int r5ppos,
                           int r3ppos)
/* return the distance from the 5' end of
region to start of fp, depending on sense
of fp.  Assume fp is non-NULL */
{
if (sens == DB_sens_comp)
  return(db_3pextnt4feat(fp) - r3ppos + 1);
else
  return(r5ppos - db_5pextnt4feat(fp) + 1);
}

DB_FEATSTRCT *igl_locatnearestfeat(IGL_RUNPARS *rpp,
                                   DB_ENTSTRCT *entsp,
                                   DB_FTYPELT *wantftlist,
                                   int rstart,
                                   int rstop)
/* look thru features locating nearest feature of feattype
3' to rstart..rstop, in sense of feature.
(any if feattyp == FTKW_unknown) in entsp.  return 
pointer to feature structure, NULL if can't be done */
{
DB_FEATSTRCT *fp;
int dist;
DB_FEATSTRCT *proxfp;
int proxdist;
/* int dist5p; */
int dist3p;
DB_TBLINFO *ilstp;
DBP_RELPOSTYPE relpos;

fp = entsp->featlst;
if ((proxdist = entsp->sqlen) <= 0)
  proxdist = db_3pextnt4feat(entsp->lastfeat);
proxfp = NULL;
while (fp != NULL)
  {
  if ((wantftlist == NULL) || (db_ptr4felt(wantftlist,fp->featur) != NULL))
    {
    ilstp = db_tblent4udata(fp->infolist,FTQU_biotype,NULL);
    if ((rpp->wantedcodingtypes == NULL) ||
         ((ilstp != NULL) && (db_matstrelt(rpp->wantedcodingtypes,ilstp->qval,strcmp) != NULL)))
      {
      dist = dist3p = igl_sens3pupstreamdist(fp,fp->strctsens,rstart,rstop);
      if (rpp->intrnlfrag > IGL_intrnl_none)
        dist += db_rawfeatlength(fp);
      relpos = dbp_relatefeat2posns(fp,rstart,rstop);
/* if (relpos == DBP_relpos_3polap)
  fprintf(stderr,"%s: %d..%d: %s feat %d..%d:dist=%d:dist3p=%d\n",
            dbp_relpostype2str(DBP_relpos_3polap),
            rstart,rstop,(fp->strctsens==DB_sens_comp?"3'":"5'"),
            db_5pextnt4feat(fp),db_3pextnt4feat(fp),dist,dist3p); */
      if ((dist < proxdist) &&
             ((relpos <= DBP_relpos_5polap) || 
               ((relpos < DBP_relpos_3pto) && (rpp->intrnlfrag > IGL_intrnl_none))) &&
             ((!rpp->needcds) || (dbp_featincludesfeat(fp,fp,FTKW_CDS))) &&
          ((rpp->distlimit == 0) || (dist3p <= rpp->distlimit)))
        {
/* fprintf(stdout,"Chk: %d..%d dist=%d,%s for %d..%d\n",db_5pextnt4feat(fp),db_3pextnt4feat(fp),
          dist,(fp->strctsens==DB_sens_comp?"3'":"5'"),rstart,rstop); */
        proxdist = dist;
        proxfp = fp;
        }
      }
    }
  fp = fp-> nextst;
  }
return(proxfp);
}

int igl_olaplength(int r1start,
                   int r1stop,
                   int r2start,
                   int r2stop)
/* return the length overlapping between
r1start..r1stop & r2start..r2stop. 0 for
no overlap.  Assume that r1start<=r2stop
and r2start<=r2stop. */
{
if (rbc_intinrng(r1start,r2start,r1stop) &&
      rbc_intinrng(r1start,r2stop,r1stop))
  return(r2stop - r2start + 1);
else
  if (rbc_intinrng(r2start,r1start,r2stop) &&
        rbc_intinrng(r2start,r1stop,r2stop))
    return(r2stop - r2stop + 1);
  else
    if (rbc_intinrng(r1start,r2start,r1stop))
      return(r1stop - r2start + 1);
    else
      if (rbc_intinrng(r1start,r2stop,r1stop))
        return(r2stop - r1start + 1);
      else
        if (rbc_intinrng(r2start,r1start,r2stop))
          return(r2stop - r1start + 1);
        else
          if (rbc_intinrng(r2start,r1stop,r2stop))
            return(r1stop - r2start + 1);
          else
            return(0);
}

int igl_featcntinrng(IGL_RUNPARS *rpp,
                     DB_ENTSTRCT *entsp,
                     DB_FTYPELT *feattype,
                     int rstart,
                     int rstop)
/* look thru features to count the number of feattype
in range rstart..rstop (any if feattyp == FTKW_unknown) in entsp.
Overlapping features are counted if their majority overlaps
rstart..rstop. */
{
DB_FEATSTRCT *fp;
DB_TBLINFO *ilstp;
DBP_RELPOSTYPE relpos;
int ft5p;
int ft3p;
int featlen;
int gcnt;

fp = entsp->featlst;
gcnt = 0;
while (fp != NULL)
  {
  if ((feattype == NULL) || (db_ptr4felt(feattype,fp->featur) != NULL))
    {
    if ((rpp->wantedcodingtypes == NULL) ||
         (((ilstp = db_tblent4udata(fp->infolist,FTQU_biotype,NULL)) != NULL) &&
           (db_matstrelt(rpp->wantedcodingtypes,ilstp->qval,strcmp) != NULL)))
      {
      ft5p = db_5pextnt4feat(fp);
      ft3p = db_3pextnt4feat(fp);
      if (((featlen = db_rawfeatlength(fp)) > 0) &&
           (igl_olaplength(rstart,rstop,ft5p,ft3p)/featlen >= 0.5))
        gcnt++;
      }
    }
  fp = fp-> nextst;
  }
return(gcnt);
}
int igl_absdist2feat(DB_FEATSTRCT *fp,
                     int spos)
/* return the absolute distance of spos
to fp 5' absolute position */
{
return(abs(db_lmt4feat(fp,imin)-spos));
}

DB_FEATSTRCT *igl_locatproximalfeat(IGL_RUNPARS *rpp,
                                    DB_ENTSTRCT *entsp,
                                    DB_FTYPELT *wantftlist,
                                    int rstart)
/* look thru features locating nearest feature of feattype
3' to rstart..rstop, 
(any if feattyp == FTKW_unknown) in entsp.  return 
pointer to feature structure, NULL if can't be done */
{
DB_FEATSTRCT *fp;
int dist;
DB_FEATSTRCT *proxfp;
int proxdist;
DB_TBLINFO *ilstp;

fp = entsp->featlst;
if ((proxdist = entsp->sqlen) <= 0)
  proxdist = db_3pextnt4feat(entsp->lastfeat);
proxfp = NULL;
while (fp != NULL)
  {
  if ((wantftlist == NULL) || (db_ptr4felt(wantftlist,fp->featur) != NULL))
    {
    if ((rpp->wantedcodingtypes == NULL) ||
         (((ilstp = db_tblent4udata(fp->infolist,FTQU_biotype,NULL)) != NULL) &&
           (db_matstrelt(rpp->wantedcodingtypes,ilstp->qval,strcmp) != NULL)))
      {
      dist = igl_absdist2feat(fp,rstart);
      if (dist < proxdist)
        {
        proxdist = dist;
        proxfp = fp;
        }
      }
    }
  fp = fp-> nextst;
  }
return(proxfp);
}

void igl_headerappend(FILE *dst,
                      IGL_RUNPARS *rpp)
/* add to existing header as far as we can */
{
int fp;

switch (rpp->omode)
  {
  case IGL_omod_genecnt:
    fprintf(dst,"%sGene_count\n",rpp->outtokdelmtr);
    break;
  default:
    fprintf(dst,"%sEnddist",rpp->outtokdelmtr);
    if ((rpp->classolap) || (rpp->intrnlfrag > IGL_intrnl_none))
      fprintf(dst,"%sOverlap",rpp->outtokdelmtr);
    if (rpp->seektss)
      fprintf(dst,"%sTSS_%s%sTSS_dist",rpp->outtokdelmtr,(rpp->tss_cpgi_rng?"range":"pos"),rpp->outtokdelmtr);
    if (rpp->seekCpGis)
      fprintf(dst,"%sCpGIs_%s%sCpGIs_dist%sCpGI_relation",rpp->outtokdelmtr,(rpp->tss_cpgi_rng?"range":"pos"),
                rpp->outtokdelmtr,rpp->outtokdelmtr);
    fprintf(dst,"%sStrand%sGeneID%sGeneRange",rpp->outtokdelmtr,rpp->outtokdelmtr,
              rpp->outtokdelmtr);
    if (rpp->displinfo)
      fprintf(dst,"%sXref...",rpp->outtokdelmtr);
    if ((rpp->dfmt == DBFMT_gtf) || (rpp->dfmt == DBFMT_sqmonk))
      fprintf(dst,"%sBiotype",rpp->outtokdelmtr);
    fputc('\n',dst);
    break;
  }
}

DB_FEATSTRCT *igl_firstupstreamfeat(DB_FEATSTRCT *fp,
                                    DB_FEATSTRCT *flst,
                                    DBLU_FTKW kw,
                                    int (* usdistfun)(DB_FEATSTRCT *fp,
                                                      DB_SENS sens,
                                                      int r5ppos,
                                                      int r3ppos))
/* scan flst for the closest feature of type kw in sense
sens and return the lst pointer.  if sens is DB_sens_unk,
then return closest. */
{
int clsdist;
DB_FEATSTRCT *closfp;
DB_FEATSTRCT *flp;
int fdist;
int fp5p;
int fp3p;

clsdist = INT32_MAX;
flp = flst;
fp5p = db_5pextnt4feat(fp);
fp3p = db_3pextnt4feat(fp);
closfp = NULL;
while (flp != NULL)
  {
  if ((kw == FTKW_unknown) || (flp->featur == kw))
    {
    if (((fdist = (*usdistfun)(flp,fp->strctsens,fp5p,fp3p)) < clsdist) &&
         (fdist > 0))
      {
      clsdist = fdist;
      closfp = flp;
      }
    }
  flp = flp->nextst;
  }
return(closfp);
}

char *igl_exintloc2str(IGL_FRAGGENE_LOC eiloc)
  /* return a string appropriate to eiloc */
{
switch (eiloc)
  {
  case IGL_frgen_upstream:
    return("upstream");
    break;
  case IGL_frgen_on5p:
    return("5prime");
    break;
  case IGL_frgen_inexon:
    return("on_exon");
    break;
  case IGL_frgen_inintron:
    return("on_intron");
    break;
  case IGL_frgen_onexonintron:
    return("exon_intron_boundary");
    break;
  case IGL_frgen_onintronexon:
    return("intron_exon_boundary");
    break;
  case IGL_frgen_includsexon:
    return("exon_internal");
    break;
  case IGL_frgen_includsintron:
    return("intron_internal");
    break;
  case IGL_frgen_on3p:
    return("3prime");
    break;
  case IGL_frgen_genein:
    return("includes_gene");
    break;
  case IGL_frgen_unknown:
    return("Unknown");
    break;
  }
}

IGL_FRAGGENE_LOC igl_classfyexintloc(IGL_RUNPARS *rpp,
                                     DB_FEATSTRCT *fp,
                                     int rgnstrt,
                                     int rgnstop)
/* given rgnstrt & rgnstop define a sequence region wrt
gene boundaries.  If internal, then further wrt 
intron exon boundaries.  Assume fp is a gene that has
corresponding mRNA or CDS following it in order.  Locate
mRNA or CDS structure and scan the segment list thereof for
the section relating to rgnstrt..rgnstop.  Return a classification
of how rgnstrt..rgnstop relate.  This code has only been checked
for effectiveness with GenBank/EMBL/SeqMonk-based feature lists.
Attempt to make it work for gtf gene features with exon positions
in list */
{
DB_FEATSTRCT *mrnacdsp;
DB_SEGELT *rstrtsegp;
DB_SEGELT *rstopsegp;
int p5featextent;
int p3featextent;

if (rgnstop < rgnstrt)
  return(igl_classfyexintloc(rpp,fp,rgnstop,rgnstrt));
else
  switch (rpp->dfmt)
    {
    case DBFMT_gtf:
      p5featextent = db_5pextnt4feat(fp);
      p3featextent = db_3pextnt4feat(fp);
      if (db_rnginclds(p5featextent,rgnstrt,p3featextent) &&
            db_rnginclds(p5featextent,rgnstop,p3featextent))
        {
        rstrtsegp = db_featseg4pos(fp,rgnstrt);
        rstopsegp = db_featseg4pos(fp,rgnstop);
        if ((rstrtsegp != NULL) && (rstopsegp != NULL))
          if (rstrtsegp == rstopsegp)
            return(IGL_frgen_inexon);
          else
            return(IGL_frgen_includsintron);
        else
          if (rstrtsegp == NULL)
            if (rstopsegp != NULL)
              if (fp->strctsens == DB_sens_comp)
                return(IGL_frgen_onexonintron);
              else
                return(IGL_frgen_onintronexon);
            else /* need to see if we have a segment between */
              if (db_nxtfeatseg4pos(fp,rgnstrt) == db_nxtfeatseg4pos(fp,rgnstop))
                return(IGL_frgen_inintron);
              else
                return(IGL_frgen_includsexon);
          else  /* rstopsegp is NULL, but not rstrtsegp */
            if (fp->strctsens == DB_sens_comp)
              return(IGL_frgen_onintronexon);
            else
              return(IGL_frgen_onexonintron);
        }
      else
        {
        if (db_rnginclds(rgnstrt,p5featextent,rgnstop))
          return(IGL_frgen_on5p);
        if (db_rnginclds(rgnstrt,p3featextent,rgnstop))
          return(IGL_frgen_on3p);
        if (db_rnginclds(rgnstrt,p5featextent,rgnstop) &&
              db_rnginclds(rgnstrt,p3featextent,rgnstop))
          return(IGL_frgen_genein);
        if ((db_rnginclds(rgnstrt,rgnstop,p5featextent) && (fp->strctsens != DB_sens_comp)) ||
              (db_rnginclds(p3featextent,rgnstrt,rgnstop) && (fp->strctsens == DB_sens_comp)))
          return(IGL_frgen_upstream);
        return(IGL_frgen_unknown);
        }
      break;
    default:
      if (((mrnacdsp = dbp_featincludesfeat(fp,fp,FTKW_mRNA)) != NULL) ||
            ((mrnacdsp = dbp_featincludesfeat(fp,fp,FTKW_CDS)) != NULL))
        {
        p5featextent = db_5pextnt4feat(fp);
        p3featextent = db_3pextnt4feat(fp);
        if (db_rnginclds(p5featextent,rgnstrt,p3featextent) &&
              db_rnginclds(p5featextent,rgnstop,p3featextent))
          {
          rstrtsegp = db_featseg4pos(mrnacdsp,rgnstrt);
          rstopsegp = db_featseg4pos(mrnacdsp,rgnstop);
          if ((rstrtsegp != NULL) && (rstopsegp != NULL))
            if (rstrtsegp == rstopsegp)
              return(IGL_frgen_inexon);
            else
              return(IGL_frgen_includsintron);
          else
            if (rstrtsegp == NULL)
              if (rstopsegp != NULL)
                if (mrnacdsp->strctsens == DB_sens_comp)
                  return(IGL_frgen_onexonintron);
                else
                  return(IGL_frgen_onintronexon);
              else /* need to see if we have a segment between */
                if (db_nxtfeatseg4pos(mrnacdsp,rgnstrt) ==
                      db_nxtfeatseg4pos(mrnacdsp,rgnstop))
                  return(IGL_frgen_inintron);
                else
                  return(IGL_frgen_includsexon);
            else  /* rstopsegp is NULL, but not rstrtsegp */
              if (mrnacdsp->strctsens == DB_sens_comp)
                return(IGL_frgen_onintronexon);
              else
                return(IGL_frgen_onexonintron);
          }
        else
          {
          if (db_rnginclds(rgnstrt,p5featextent,rgnstop))
            return(IGL_frgen_on5p);
          if (db_rnginclds(rgnstrt,p3featextent,rgnstop))
            return(IGL_frgen_on3p);
          if (db_rnginclds(rgnstrt,p5featextent,rgnstop) &&
                db_rnginclds(rgnstrt,p3featextent,rgnstop))
            return(IGL_frgen_genein);
          if ((db_rnginclds(rgnstrt,rgnstop,p5featextent) && (fp->strctsens != DB_sens_comp)) ||
                (db_rnginclds(p3featextent,rgnstrt,rgnstop) && (fp->strctsens == DB_sens_comp)))
            return(IGL_frgen_upstream);
          return(IGL_frgen_unknown);
          }
        }
      break;
    }
return(IGL_frgen_unknown);
}

char *igl_cpgiloc2str(IGL_CPGI_LOC cpgiloc)
  /* return a string for a cpgi location */
{
switch (cpgiloc)
  {
  case IGL_cpgi_oncore:
    return("CpGI_core");
    break;
  case IGL_cpgi_coreshore:
    return("CpGI_coreshore");
    break;
  case IGL_cpgi_onshore:
    return("CpGI_shore");
    break;
  case IGL_cpgi_shoreshelf:
    return("CpGI_shoreshelf");
    break;
  case IGL_cpgi_onshelf:
    return("CpGI_shelf");
    break;
  case IGL_cpgi_shelfolap:
    return("CpGI_shelfoverlap");
    break;
  case IGL_cpgi_unknown:    /* unknown or not on island,shore,shelf */
  default:
    return("-");
    break;
  }
}

IGL_CPGI_LOC igl_classfragcpgi(int fstart,
                               int fstop,
                               int cpgistart,
                               int cpgistop)
/* compare the fragment fstart..fstop to CpGIs
cpgistart..cpgistop and return an appropriate
position */
{
int shelflow;
int shelfhigh;
int shorelow;
int shorehigh;

if (fstart > fstop)
  return(igl_classfragcpgi(fstop,fstart,cpgistart,cpgistop));
else
  if (cpgistart > cpgistop)
    return(igl_classfragcpgi(fstart,fstop,cpgistop,cpgistart));
  else
    {
    shorelow = imax(1,cpgistart-2000);
    shorehigh = cpgistop + 2000;
    shelflow = imax(1,shorelow-2000);
    shelfhigh = shorehigh + 2000;
    if (!db_rnginclds(shelflow,fstart,shelfhigh) && !db_rnginclds(shelflow,fstop,shelfhigh))
      return(IGL_cpgi_unknown);      /* bomb out early if not in range at all */
    else
      if (db_rnginclds(cpgistart,fstart,cpgistop) && db_rnginclds(cpgistart,fstop,cpgistop))
        return(IGL_cpgi_oncore);
      else
        if ((db_rnginclds(shorelow,fstart,cpgistart) && db_rnginclds(cpgistart,fstop,cpgistop)) ||
             (db_rnginclds(cpgistart,fstart,cpgistop) && db_rnginclds(cpgistop,fstop,shorehigh)))
          return(IGL_cpgi_coreshore);
        else
          if ((db_rnginclds(shorelow,fstart,cpgistart) && db_rnginclds(shorelow,fstop,cpgistart)) ||
               (db_rnginclds(cpgistop,fstart,shorehigh) && db_rnginclds(cpgistop,fstop,shorehigh)))
            return(IGL_cpgi_onshore);
          else
            if ((db_rnginclds(shelflow,fstart,shorelow) && db_rnginclds(shorelow,fstop,cpgistart)) ||
                  (db_rnginclds(cpgistop,fstart,shorehigh) && db_rnginclds(shorehigh,fstop,shelfhigh)))
              return(IGL_cpgi_shoreshelf);
            else
              if ((db_rnginclds(shelflow,fstart,shorelow) && db_rnginclds(shelflow,fstop,shorelow)) ||
                    (db_rnginclds(shorehigh,fstart,shelfhigh) && db_rnginclds(shorehigh,fstop,shelfhigh)))
                return(IGL_cpgi_onshelf);
              else
                if (db_rnginclds(shelflow,fstop,shelfhigh) || db_rnginclds(shelflow,fstart,shelfhigh))
                  return(IGL_cpgi_shelfolap);
    }
return(IGL_cpgi_unknown);    /* should have returned before here, but... */
}

void igl_procchrposfile(FILE *src,
                        FILE *dst,
                        IGL_RUNPARS *rpp)
/* read successivelines from src, expecting
Chr\tStart\tStop\t....\n.  If Chr matches
rpp-uchrno then try to find nearest gene
from feature table */
{
char *lbuf;
char *lcpy;
/* char *oricpy; */
char **lp;
char **tokns;
int tcnt;
int chrno;
DBP_MULTI_ELEMNT *elpt;
int istrt;
int istop;
DB_FEATSTRCT *proxftp;
DB_TBLINFO *infp;
int tpt;
DB_FEATSTRCT *upstrmftp;
int usd;
DBP_RELPOSTYPE rptype;
DB_FTYPELT *feattype;
DB_FEATSTRCT *proxcpgip;
char *gnamptr;
int infocnt;
DB_STR_ELT *gtptr;
DBLU_FTKW gtkw;

/* allow some extra length to encompass terminal '\n' */
lbuf = (char *) getmemory(rpp->srcbuflen+8,"srclinebuf");
tokns = (char **) getmemory(rpp->srccollmt*sizeof(char *),"tokenlist");
feattype = NULL;
switch (rpp->dfmt)
  {
  case DBFMT_gtf:
    if (rpp->dfmt != DBFMT_gtf)
      {
      if (rpp->ufeattype == FTKW_unknown)
        if (rpp->geneexmrna)
          (void) db_appnfelt(&feattype,FTKW_mRNA);
        else
          if (rpp->cds_prot_id)
            (void) db_appnfelt(&feattype,FTKW_CDS);
          else
            {
            (void) db_appnfelt(&feattype,FTKW_FT_gene);
	    (void) db_appnfelt(&feattype,FTKW_lncRNA);
	    }
      else
        (void) db_appnfelt(&feattype,rpp->ufeattype);
      }
    break;
  case DBFMT_genbank:
  case DBFMT_embl:
  case DBFMT_embl_old:
    if (rpp->ufeattype == FTKW_unknown)
      if (rpp->geneexmrna)
        (void) db_appnfelt(&feattype,FTKW_mRNA);
       else
        (void) db_appnfelt(&feattype,FTKW_CDS);  /* gene name seems to be put on CDS lines */
    else
      (void) db_appnfelt(&feattype,rpp->ufeattype);
    break;
  default:
    if (rpp->gtfgeneidattr == NULL)
      (void) db_appnfelt(&feattype,FTKW_FT_gene);
    else
      {
      gtptr = rpp->gtfgeneidattr;
      while (gtptr != NULL)
        {
        if ((gtkw = wlu_chkwrd(rpp->ftkw,gtptr->strval)) != FTKW_unknown)
          (void) db_appnfelt(&feattype,gtkw);
        gtptr = gtptr->nxtselt;
        }
      }
    break;
  }
while (fgets(lbuf,(rpp->srcbuflen+7),src) != NULL)
  {
  igl_cleansrcline(lbuf);
  rpp->srclno++;
/*   if (*lbuf != '#') */
  lcpy = oricpy = bas_strdup(lbuf);
  tcnt = 0;
  for (lp = tokns; (*lp = strsep(&lcpy,rpp->chrtokdelmtr)) != NULL;)
    {
    if (**lp != '\0')
      {
      tcnt++;
      if (++lp >= (tokns+rpp->srccollmt))
        break;
      }
    }
/* check for adequate no of tokens: avoid file corruption issues */
  if (tcnt < 3)
    err_msg_die("Position file corruption at line %d: inadequate tokens in line\n",rpp->srclno);
  if (*lbuf == '#')  /* comment or header line */
    {
    for (tpt = 0; tpt < rpp->srccollmt; tpt++)
      {
      if (tpt < tcnt)
        fputs(*(tokns+tpt),dst);
      else
        fputc('-',dst);
      if (tpt != (rpp->srccollmt-1))
        fputs(rpp->outtokdelmtr,dst);
      }
    if (rpp->srclno == 1)  /* is header line */
      igl_headerappend(dst,rpp);
    else
      fputc('\n',dst);
    }
  else
    {
    chrno = wlu_chkwrd(rpp->chridlu,*tokns);
    if ((rpp->uchrno <= 0) || ((rpp->uchrno-1) == chrno))
      {
      if ((elpt = dbp_melemnt4id(rpp->meplist,*tokns,strcmp)) != NULL)
        {
        istrt = (int) strtol(*(tokns+1),NULL,10);
        istop = (int) strtol(*(tokns+2),NULL,10);
        switch (rpp->omode)
          {
          case IGL_omod_fragbased:
            proxftp = proxcpgip = NULL;
            if (rpp->seektss)
              proxftp = igl_locatnearestfeat(rpp,elpt->me_entstrptr,rpp->tsslist,istrt,istop);
            if ((((rpp->seekCpGis == IGL_CpG_upstream)) &&
                 ((proxcpgip = igl_locatnearestfeat(rpp,elpt->me_entstrptr,rpp->cpgilist,
                                                      istrt,istop)) != NULL) ||
                 ((proxcpgip = igl_locatproximalfeat(rpp,elpt->me_entstrptr,rpp->cpgilist,
                                                      istrt)) != NULL)) &&
                 (proxftp != NULL))
              {
              fputs(*tokns,dst);
              tpt = 1;
              while (tpt < tcnt)
                {
                fputs(rpp->outtokdelmtr,dst);
                fputs(*(tokns+tpt),dst);
                tpt++;
                }
              while (tpt < rpp->srccollmt)
                {
                fputs(rpp->outtokdelmtr,dst);
                fputc('-',dst);
                tpt++;
                }
              if (proxftp != NULL)
                {
                fprintf(dst,"%s%d",rpp->outtokdelmtr,db_5pextnt4feat(proxftp));
                if (rpp->tss_cpgi_rng)
                  fprintf(dst,"-%d",db_3pextnt4feat(proxftp));
                fprintf(dst,"%s%d",rpp->outtokdelmtr,
                          igl_sens5pupstreamdist(proxftp,DB_sens_norm,istrt,istop));
                }
              else
                fprintf(dst,"%s0%s0",rpp->outtokdelmtr,rpp->outtokdelmtr);
              if (proxcpgip != NULL)
                {
                fprintf(dst,"%s%d",rpp->outtokdelmtr,db_5pextnt4feat(proxcpgip));
                if (rpp->tss_cpgi_rng)
                  fprintf(dst,"-%d",db_3pextnt4feat(proxcpgip));
                fprintf(dst,"%s%d",rpp->outtokdelmtr,
                          igl_sens5pupstreamdist(proxcpgip,DB_sens_norm,istrt,istop));
                }
              else
                fprintf(dst,"%s0%s0",rpp->outtokdelmtr,rpp->outtokdelmtr);
              fputc('\n',dst);
              }
            break;
          case IGL_omod_genecnt:
            fputs(*tokns,dst);
            tpt = 1;
            while (tpt < tcnt)
              {
              fputs(rpp->outtokdelmtr,dst);
              fputs(*(tokns+tpt),dst);
              tpt++;
              }
            while (tpt < rpp->srccollmt)
              {
              fputs(rpp->outtokdelmtr,dst);
              fputc('-',dst);
              tpt++;
              }
            fprintf(dst,"%s%d\n",rpp->outtokdelmtr,
                      igl_featcntinrng(rpp,elpt->me_entstrptr,feattype,istrt,istop));
            break;
          default:
          case IGL_omod_genebased:
            if ((proxftp =
                   igl_locatnearestfeat(rpp,elpt->me_entstrptr,feattype,istrt,istop)) != NULL)
              {
              fputs(*tokns,dst);
              tpt = 1;
              while (tpt < tcnt)
                {
                fputs(rpp->outtokdelmtr,dst);
                fputs(*(tokns+tpt),dst);
                tpt++;
                }
              while (tpt < rpp->srccollmt)
                {
                fputs(rpp->outtokdelmtr,dst);
                fputc('-',dst);
                tpt++;
                }
              fprintf(dst,"%s%d",rpp->outtokdelmtr,
                        igl_sens3pupstreamdist(proxftp,proxftp->strctsens,istrt,istop));
              if ((rpp->classolap) || (rpp->intrnlfrag > IGL_intrnl_none))
                {
                switch (rptype = dbp_relatefeat2posns(proxftp,istrt,istop))
                  {
                  case DBP_relpos_5pto:
                    if (rpp->classolap)
                      fprintf(dst,"%supstream",rpp->outtokdelmtr);
                    else
                      fprintf(dst,"%s-",rpp->outtokdelmtr);
                    break;
                  case DBP_relpos_included:
                    if (rpp->classolap)
                      fprintf(dst,"%s%s",rpp->outtokdelmtr,dbp_relpostype2str(rptype));
                    else
                      fprintf(dst,"%s-",rpp->outtokdelmtr);
                    break;
                  case DBP_relpos_5polap:
                  case DBP_relpos_internal:
                  case DBP_relpos_3polap:
                    if (rpp->intrnlfrag == IGL_intrnl_resolv)
                      fprintf(dst,"%s%s",rpp->outtokdelmtr,
                                ((rptype==DBP_relpos_internal)?
                                igl_exintloc2str(igl_classfyexintloc(rpp,proxftp,istrt,istop)):"-"));
                    else
                      fprintf(dst,"%s%s",rpp->outtokdelmtr,dbp_relpostype2str(rptype));
                    break;
                  case DBP_relpos_3pto:
                  case DBP_relpos_undef:
                  default:
                    fprintf(dst,"%s-",rpp->outtokdelmtr);
                    break;
                  }
                }
/*              if (rpp->intrnlfrag == IGL_intrnl_resolv)
                fprintf(dst,"%s%s",rpp->outtokdelmtr,
                          ((rptype==DBP_relpos_internal)?
                             igl_exintloc2str(igl_classfyexintloc(proxftp,istrt,istop)):"-")); */
              if (rpp->seektss)
                {
                upstrmftp = igl_firstupstreamfeat(proxftp,elpt->me_entstrptr->featlst,FTKW_TSS,
                                                    igl_sens3pupstreamdist);
                if ((upstrmftp != NULL) &&
                     (((usd = igl_sensupstreamdist(upstrmftp,proxftp->strctsens,istrt,istop))
                                <= rpp->tssdistlmt) || (rpp->tssdistlmt == 0)))
                  {
                  fprintf(dst,"%s%d",rpp->outtokdelmtr,db_5pextnt4feat(upstrmftp));
                  if (rpp->tss_cpgi_rng)
                    fprintf(dst,"-%d",db_3pextnt4feat(upstrmftp));
                  fprintf(dst,"%s%d",rpp->outtokdelmtr,usd);
                  }
                else
                  fprintf(dst,"%s0%s0",rpp->outtokdelmtr,rpp->outtokdelmtr);
                }
              if (((rpp->seekCpGis == IGL_CpG_upstream) &&
                   ((upstrmftp = igl_firstupstreamfeat(proxftp,elpt->me_entstrptr->featlst,FTKW_CpG_is,
                                                         igl_sens5pupstreamdist)) != NULL)) ||
                   ((rpp->seekCpGis == IGL_CpG_proximal) &&
                   ((upstrmftp = igl_locatproximalfeat(rpp,elpt->me_entstrptr,rpp->cpgilist,
                                                         istrt)) != NULL)))
                {
                fprintf(dst,"%s%d",rpp->outtokdelmtr,db_5pextnt4feat(upstrmftp));
                if (rpp->tss_cpgi_rng)
                  fprintf(dst,"-%d",db_3pextnt4feat(upstrmftp));
                fprintf(dst,"%s%d",rpp->outtokdelmtr,
                          igl_sens5pupstreamdist(upstrmftp,proxftp->strctsens,istrt,istop));
                fprintf(dst,"%s%s",rpp->outtokdelmtr,
                          igl_cpgiloc2str(igl_classfragcpgi(istrt,istop,db_5pextnt4feat(upstrmftp),
                                                              db_3pextnt4feat(upstrmftp))));
                }
              else
                {
                if (rpp->tss_cpgi_rng)
                  fprintf(dst,"%s-",rpp->outtokdelmtr);
                fprintf(dst,"%s-%s-",rpp->outtokdelmtr,rpp->outtokdelmtr);
                }
/*                fprintf(dst,"%s0%s0%s%s",rpp->outtokdelmtr,rpp->outtokdelmtr,rpp->outtokdelmtr,
                          igl_cpgiloc2str(IGL_cpgi_unknown)); */
              fprintf(dst,"%s%s'",rpp->outtokdelmtr,(proxftp->strctsens==DB_sens_comp?"3":"5"));
              if ((infp = db_tblent4udata(proxftp->infolist,FTQU_name,NULL)) != NULL)
                fprintf(dst,"%s%s%s%d-%d",rpp->outtokdelmtr,infp->qval,rpp->outtokdelmtr,
                          db_5pextnt4feat(proxftp),db_3pextnt4feat(proxftp));
              else
                {
                gnamptr = NULL;
                if (rpp->cds_prot_id)
                  infp = db_tblent4udata(proxftp->infolist,FTQU_protid,NULL);
                else
                  infp = db_tblent4udata(proxftp->infolist,FTQU_gene,NULL);
                if (infp != NULL)
                  gnamptr = infp->qval;
                else
                  gnamptr = proxftp->idp;
                fprintf(dst,"%s%s%s%d-%d",rpp->outtokdelmtr,(gnamptr!=NULL?gnamptr:"<No Name>"),
                            rpp->outtokdelmtr,db_5pextnt4feat(proxftp),db_3pextnt4feat(proxftp));
                }
              if (rpp->displinfo)
                {
                infocnt = 0;
                infp = db_tblent4udata(proxftp->infolist,FTQU_db_xref,NULL);
                if (infp != NULL)
                  fputs(rpp->outtokdelmtr,dst);
                while (infp != NULL)
                  {
                  if (infocnt > 0)
                    fputs("; ",dst);
                  fprintf(dst,"%s",infp->qval);
                  infocnt++;
                  infp = db_tblent4udata(infp->nxtielt,FTQU_db_xref,NULL);
                  }
                infp = db_tblent4udata(proxftp->infolist,FTQU_product,NULL);
                while (infp != NULL)
                  {
                  if (infocnt > 0)
                    fputs("; ",dst);
                  fprintf(dst,"%s",infp->qval);
                  infocnt++;
                  infp = db_tblent4udata(infp->nxtielt,FTQU_product,NULL);
                  }
                }
              if (((rpp->dfmt == DBFMT_gtf) || (rpp->dfmt == DBFMT_sqmonk)) &&
	        (infp = db_tblent4udata(proxftp->infolist,FTQU_biotype,NULL)))
/* want to put biotype output at end of line */
	       fprintf(dst,"%s%s",rpp->outtokdelmtr,infp->qval);
              fputc('\n',dst);
              }
            else
              if (rpp->printnullhits)
                {
                tpt = 0;
                while (tpt < tcnt)
                  {
                  fprintf(dst,"%s%s",(tpt>0?rpp->outtokdelmtr:""),*(tokns+tpt));
                  tpt++;
                  }
                while (tpt < rpp->srccollmt)
                  {
                  fprintf(dst,"%s-",rpp->outtokdelmtr);
                  tpt++;
                  }
                fputc('\n',dst);
                }
            break;
          }
        }
      }
    }
  memfree(oricpy);
  }
db_killfeltlst(&feattype);
memfree(lbuf);
}

int igl_scanchrinfofl(FILE *cifile,
                      IGL_CHR_INFO *cinfo,
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

int igl_maknamsnparsftfls(IGL_RUNPARS *rpp)
  /* either use given feature table file name
or generate using prefix/suffix and try to
read.  Die on failures with appropriate message.
Return total number of feat entries stored */
{
char *fnm;
int fcnt;
int chrno;

rpp->chridlu = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"ChrID lu struct");
wlu_initlustrct(rpp->chridlu,WLU_CASEIND,-1);
if (rpp->ftflnam != NULL)
  {
  for (chrno = 0; chrno <= rpp->maxchrno; chrno++)
    wlu_addwrd(rpp->chridlu,any_chrno2str(rpp->maxchrno,rpp->havexy,chrno,1),
                  (int) chrno,NULL);
  fcnt = igl_parseftflnam(rpp,rpp->ftflnam);
  }
else
  {
  if ((rpp->ftprefix == NULL) && (rpp->chrinfofl != NULL))
    { /* scan info file to find number chromosomes */
    rpp->maxchrno = igl_scanchrinfofl(rpp->chrinfofl,NULL,NULL);
    rewind(rpp->chrinfofl);
    }
  rpp->chrinfo = (IGL_CHR_INFO *) getmemory(sizeof(IGL_CHR_INFO)*rpp->maxchrno,
                                              "Chr bin info");
/* init these */
  for (chrno = 0; chrno < rpp->maxchrno; chrno++)
    {
    (rpp->chrinfo+chrno)->chrflname = NULL;
    (rpp->chrinfo+chrno)->chrid = NULL;
    }
  if (rpp->ftprefix != NULL)
    {
    for (chrno = 0; chrno < rpp->maxchrno; chrno++)
     if ((rpp->uchrno == 0) || ((rpp->uchrno-1) == chrno))
        if (asprintf(&fnm,"%s%s%s",rpp->ftprefix,
                       any_chrno2str(rpp->maxchrno,rpp->havexy,chrno+1,1),
                       rpp->ftsuffix) > 0)
          {
          (rpp->chrinfo+chrno)->chrflname = fnm;
          (rpp->chrinfo+chrno)->chrid = bas_strdup(any_chrno2str(rpp->maxchrno,
                                                                   rpp->havexy,(chrno+1),1));
          }
    }
  else
    {
    (void) igl_scanchrinfofl(rpp->chrinfofl,rpp->chrinfo,NULL);
    fclose(rpp->chrinfofl);
    }
  fcnt = 0;
  for (chrno = 0; chrno < rpp->maxchrno; chrno++)
    if ((rpp->chrinfo+chrno)->chrid != NULL)
      {
      wlu_addwrd(rpp->chridlu,(rpp->chrinfo+chrno)->chrid,chrno,NULL);
      fcnt += igl_parseftflnam(rpp,(rpp->chrinfo+chrno)->chrflname);
      }
  }
if (fcnt <= 0)
  {
  if ((rpp->ftflnam == NULL) && (rpp->ftprefix == NULL) &&
        (rpp->chrinfoflnam == NULL))
    err_msg_die("Error: no names given for -f, -p or -G\n");
  else
    err_msg_die("Error: no appropriate feature information recovered\n");
  }
return(fcnt);
}

int igl_chkcommadelstrs(DB_STR_ELT **slst,
                        char *astr)
/* astr may or may not be a comma-delimited
set of strings.  put what is there into *slst.
return number elements found if it is any use */
{
char *scopy;
char *scpyori;
int tcnt;
char *cp;
char **tokns;
char **lp;
int ttcnt;
int kcnt;

scopy = scpyori = bas_strdup(astr);
tcnt = 1;
cp = astr;
while (*cp != '\0')
  {
  if (*cp == ',')
    tcnt++;
  cp++;
  }
tokns = (char **) getmemory(tcnt*sizeof(char *),"tokenlist");
ttcnt = 0;
for (lp = tokns; (*lp = strsep(&scopy,",")) != NULL;)
  {
  if (**lp != '\0')
    {
    ttcnt++;
    if (++lp >= (tokns+tcnt))
      break;
    }
  }
kcnt = 0;
tcnt = 0;
while (kcnt < ttcnt)
  {
  (void) dbp_appstrelt(slst,*(tokns+kcnt));
  kcnt++;
  tcnt++;
  }
memfree(scpyori);
memfree(tokns);
return(tcnt);
}

int igl_chkvalidbiotypelist(WRD_LUSTRCT *validwrds,
                            DB_STR_ELT **slst,
                            char *errhdr,
                            WRD_LUSTRCT **wantlst)
/* Check *slst entries are in validwrds,
writing an error message if not (non-fatal).
Remove any invalid entries.
Return number of valid words and
generate *wantlst lookup list */
{
DB_STR_ELT *dseptr;
int scancnt;
DBP_BIOTYPE biotype;

dseptr = *slst;
scancnt = 0;
while (dseptr != NULL)
  {
  if ((biotype = wlu_chkwrd(validwrds,dseptr->strval)) == BIOTYPE_unknown)
    {
    err_msg("Unknown %s word '%s'\n",((errhdr==NULL)?"":errhdr),dseptr->strval);
    db_delstrelt(dseptr,slst,NULL);
    }
  else
    {
    wlu_addwrd(*wantlst,dseptr->strval,biotype,NULL);
    scancnt++;
    }
  dseptr = dseptr->nxtselt;
  }
return(scancnt);
}

int igl_chkbiotypes4fmt(DB_STR_ELT **slst,
                        DBLU_DBFMT dfmt,
                        WRD_LUSTRCT **wantlst)
/* works for seqmonk and gtf formats:
others default to seqmonk.
Ensure that astr entries are OK for
this format, return valid No. */
{
WRD_LUSTRCT *okwrds;
int ret;

okwrds = wlu_getlustrct(WLU_CASEIND,BIOTYPE_unknown);
switch (dfmt)
  {
  case DBFMT_gtf:
    wlu_addwrd(okwrds,"IG_C_gene",BIOTYPE_IG_C_gene,NULL);
    wlu_addwrd(okwrds,"IG_C_pseudogene",BIOTYPE_IG_C_pseudogene,NULL);
    wlu_addwrd(okwrds,"IG_D_gene",BIOTYPE_IG_D_gene,NULL);
    wlu_addwrd(okwrds,"IG_J_gene",BIOTYPE_IG_J_gene,NULL);
    wlu_addwrd(okwrds,"IG_J_pseudogene",BIOTYPE_IG_J_pseudogene,NULL);
    wlu_addwrd(okwrds,"IG_V_gene",BIOTYPE_IG_V_gene,NULL);
    wlu_addwrd(okwrds,"IG_V_pseudogene",BIOTYPE_IG_V_pseudogene,NULL);
    wlu_addwrd(okwrds,"IG_pseudogene",BIOTYPE_IG_pseudogene,NULL);
    wlu_addwrd(okwrds,"Mt_rRNA",BIOTYPE_Mt_rRNA,NULL);
    wlu_addwrd(okwrds,"Mt_tRNA",BIOTYPE_Mt_tRNA,NULL);
    wlu_addwrd(okwrds,"TEC",BIOTYPE_TEC,NULL);
    wlu_addwrd(okwrds,"TR_C_gene",BIOTYPE_TR_C_gene,NULL);
    wlu_addwrd(okwrds,"TR_D_gene",BIOTYPE_TR_D_gene,NULL);
    wlu_addwrd(okwrds,"TR_J_gene",BIOTYPE_TR_J_gene,NULL);
    wlu_addwrd(okwrds,"TR_J_pseudogene",BIOTYPE_TR_J_pseudogene,NULL);
    wlu_addwrd(okwrds,"TR_V_gene",BIOTYPE_TR_V_gene,NULL);
    wlu_addwrd(okwrds,"TR_V_pseudogene",BIOTYPE_TR_V_pseudogene,NULL);
    wlu_addwrd(okwrds,"antisense",BIOTYPE_antisense,NULL);
    wlu_addwrd(okwrds,"lincRNA",BIOTYPE_lincRNA,NULL);
    wlu_addwrd(okwrds,"lncRNA",BIOTYPE_lncRNA,NULL);
    wlu_addwrd(okwrds,"miRNA",BIOTYPE_miRNA,NULL);
    wlu_addwrd(okwrds,"misc_RNA",BIOTYPE_misc_RNA,NULL);
    wlu_addwrd(okwrds,"polymorphic_pseudogene",BIOTYPE_polymorphic_pseudogene,NULL);
    wlu_addwrd(okwrds,"processed_pseudogene",BIOTYPE_processed_pseudogene,NULL);
    wlu_addwrd(okwrds,"processed_transcript",BIOTYPE_processed_transcript,NULL);
    wlu_addwrd(okwrds,"protein_coding",BIOTYPE_protein_coding,NULL);
    wlu_addwrd(okwrds,"pseudogene",BIOTYPE_pseudogene,NULL);
    wlu_addwrd(okwrds,"rRNA",BIOTYPE_rRNA,NULL);
    wlu_addwrd(okwrds,"rRNA_pseudogene",BIOTYPE_rRNA_pseudogene,NULL);
    wlu_addwrd(okwrds,"scRNA",BIOTYPE_scRNA,NULL);
    wlu_addwrd(okwrds,"sense_intronic",BIOTYPE_sense_intronic,NULL);
    wlu_addwrd(okwrds,"sense_overlapping",BIOTYPE_sense_overlapping,NULL);
    wlu_addwrd(okwrds,"snRNA",BIOTYPE_snRNA,NULL);
    wlu_addwrd(okwrds,"snoRNA",BIOTYPE_snoRNA,NULL);
    wlu_addwrd(okwrds,"transcribed_processed_pseudogene",BIOTYPE_transcribed_processed_pseudogene,NULL);
    wlu_addwrd(okwrds,"transcribed_unitary_pseudogene",BIOTYPE_transcribed_unitary_pseudogene,NULL);
    wlu_addwrd(okwrds,"transcribed_unprocessed_pseudogene",BIOTYPE_transcribed_unprocessed_pseudogene,NULL);
    wlu_addwrd(okwrds,"translated_processed_pseudogene",BIOTYPE_translated_processed_pseudogene,NULL);
    wlu_addwrd(okwrds,"translated_unprocessed_pseudogene",BIOTYPE_translated_unprocessed_pseudogene,NULL);
    wlu_addwrd(okwrds,"unitary_pseudogene",BIOTYPE_unitary_pseudogene,NULL);
    wlu_addwrd(okwrds,"unprocessed_pseudogene",BIOTYPE_unprocessed_pseudogene,NULL);
    wlu_addwrd(okwrds,"vaultRNA",BIOTYPE_vaultRNA,NULL);
    *wantlst = wlu_getlustrct(WLU_CASEIND,BIOTYPE_unknown);
    ret = igl_chkvalidbiotypelist(okwrds,slst,"GTF gene_type",wantlst);
    break;
  case DBFMT_sqmonk:
  default:
    wlu_addwrd(okwrds,"Mt_tRNA_pseudogene",BIOTYPE_Mt_tRNA_pseudogene,NULL);
    wlu_addwrd(okwrds,"TEC",BIOTYPE_TEC,NULL);
    wlu_addwrd(okwrds,"ambiguous_orf",BIOTYPE_ambiguous_orf,NULL);
    wlu_addwrd(okwrds,"antisense",BIOTYPE_antisense,NULL);
    wlu_addwrd(okwrds,"lincRNA",BIOTYPE_lincRNA,NULL);
    wlu_addwrd(okwrds,"miRNA",BIOTYPE_miRNA,NULL);
    wlu_addwrd(okwrds,"miRNA_pseudogene",BIOTYPE_miRNA_pseudogene,NULL);
    wlu_addwrd(okwrds,"misc_RNA",BIOTYPE_misc_RNA,NULL);
    wlu_addwrd(okwrds,"non_coding",BIOTYPE_non_coding,NULL);
    wlu_addwrd(okwrds,"nonsense_mediated_decay",BIOTYPE_nonsense_mediated_decay,NULL);
    wlu_addwrd(okwrds,"polymorphic_pseudogene",BIOTYPE_polymorphic_pseudogene,NULL);
    wlu_addwrd(okwrds,"processed_transcript",BIOTYPE_processed_transcript,NULL);
    wlu_addwrd(okwrds,"protein_coding",BIOTYPE_protein_coding,NULL);
    wlu_addwrd(okwrds,"pseudogene",BIOTYPE_pseudogene,NULL);
    wlu_addwrd(okwrds,"rRNA",BIOTYPE_rRNA,NULL);
    wlu_addwrd(okwrds,"rRNA_pseudogene",BIOTYPE_rRNA_pseudogene,NULL);
    wlu_addwrd(okwrds,"retained_intron",BIOTYPE_retained_intron,NULL);
    wlu_addwrd(okwrds,"scRNA_pseudogene",BIOTYPE_scRNA_pseudogene,NULL);
    wlu_addwrd(okwrds,"sense_intronic",BIOTYPE_sense_intronic,NULL);
    wlu_addwrd(okwrds,"snRNA",BIOTYPE_snRNA,NULL);
    wlu_addwrd(okwrds,"snRNA_pseudogene",BIOTYPE_snRNA_pseudogene,NULL);
    wlu_addwrd(okwrds,"snoRNA",BIOTYPE_snoRNA,NULL);
    wlu_addwrd(okwrds,"snoRNA_pseudogene",BIOTYPE_snoRNA_pseudogene,NULL);
    wlu_addwrd(okwrds,"tRNA_pseudogene",BIOTYPE_tRNA_pseudogene,NULL);
    *wantlst = wlu_getlustrct(WLU_CASEIND,BIOTYPE_unknown);
    ret = igl_chkvalidbiotypelist(okwrds,slst,"SeqMonk biotype",wantlst);
    break;
  }    
wlu_clrlustrct(okwrds);
return(ret);
}

int main(int argc,
         char *argv[])
{
FILE *src;
DBP_MULTI_ELEMNT *curseq;
int modcnt;
int thiscnt;
DB_FEATSTRCT *fp;
int ap;
char uopt;
int lcnt;
DB_SEGELT *sp;
IGL_RUNPARS *rpp;
DB_FTYPELT *fep;
WRD_LUSTRCT *tmpkwlu;
DBLU_FTKW kwp;
DB_STR_ELT *sfp;
DBLU_FTKW feattype;

bas_initmalfl("igl_mallc.txt");
src = NULL;
rpp = (IGL_RUNPARS *) getmemory(sizeof(IGL_RUNPARS),"Runpar Struct");
rpp->srcflnam = rpp->ftflnam = rpp->ftprefix =rpp->ftsuffix = NULL;
rpp->srclno = 0;
rpp->dfmt = DBFMT_unknown;
rpp->uchrno = 0;
rpp->meplist = NULL;
rpp->ftlist = NULL;
rpp->ftkw = NULL;
rpp->havexy = 1;
rpp->maxchrno = ChrY;
rpp->rstyle = IGL_findproximal;
rpp->srcbuflen = DEF_SRCBUFLEN;
rpp->chrtokdelmtr = "\t ";
rpp->outtokdelmtr = "\t";
rpp->srccollmt = 9;
rpp->needcds = rpp->classolap = rpp->seektss =
  rpp->displinfo = rpp->distlimit = rpp->tssdistlmt = rpp->geneexmrna = 0;
rpp->seekCpGis = IGL_CpG_none;
rpp->wantedcodingtypes = NULL;
rpp->wantedtypewlu = NULL;
rpp->intrnlfrag = IGL_intrnl_none;
rpp->tss_cpgi_rng = 0;
rpp->omode = IGL_omod_genebased;
rpp->chridlu = NULL;
rpp->chrinfoflnam = NULL;
rpp->chrinfofl = NULL;
rpp->cds_prot_id = 0;
rpp->gtfgeneidattr = NULL;
rpp->ufeattype = FTKW_unknown;
rpp->ufeatlist = NULL;
rpp->tsslist = rpp->cpgilist = NULL;
rpp->printnullhits = 0;
rpp->outfile = stdout;
for (ap = 1; ap < argc; ap++)
  if (*argv[ap] == '-') /* option */
    switch (uopt = *(argv[ap]+1))
      {
      case 'e':
      case 'g':
      case 'E':
      case 'F':
      case 'T':
      case 'Q':
        rpp->dfmt = db_chr2dbfmt(uopt);
        rpp->ftkw = db_getkwstrct(rpp->dfmt);
        break;
      case 'f':       /* file name */
        if (++ap > argc)
          err_msg_die("%c needs feature file name\n",uopt);
        else
          {
          rpp->ftflnam = bas_strdup(argv[ap]);
          }
        break;
      case 'p':       /* file name prefix string */
        if (++ap > argc)
          err_msg_die("%c needs feature file prefix string\n",uopt);
        else
          rpp->ftprefix = bas_strdup(argv[ap]);
        break;
      case 's':       /* file name suffix string */
        if (++ap > argc)
          err_msg_die("%c needs feature file suffix string\n",uopt);
        else
          rpp->ftsuffix = bas_strdup(argv[ap]);
        break;
      case 'r':  /* region read/position info */
        if (++ap > argc)
          err_msg_die("%c needs chromosome region info file name\n",uopt);
        else
          if (strcmp(argv[ap],"-") == 0)      /* use stdin */
            {
            src = stdin;
            rpp->srcflnam = "stdin";
            }
          else
            if ((src = fopen(argv[ap],"r")) != NULL)
              rpp->srcflnam = bas_strdup(argv[ap]);
            else
              err_msg_die("Can't open CpG info file '%s'\n",argv[ap]);
        break;
      case 'k':     /* reset maxchrno */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",uopt);
        else
          rpp->maxchrno = (int) strtol(argv[ap],NULL,10);
        break;
      case 'G':   /* file name for chromosome spec file */
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",uopt);
        else
          {
          rpp->chrinfoflnam = bas_strdup(argv[ap]);
          if ((rpp->chrinfofl = fopen(rpp->chrinfoflnam,"r")) == NULL)
            err_msg_die("Can't open chromosome information file '%s' for -%c\n",
                          rpp->chrinfoflnam,uopt);
          }
        break;
      case 'z':       /* disable XY chromosome labelling */
        rpp->havexy = 0;
        break;
      case 'l':       /* list feature stuff (check) */
        rpp->rstyle = IGL_listfeatdata;
        break;
      case 'L':      /* list chromosome ids */
        rpp->rstyle = IGL_listchrids;
        break;
      case 'c':     /* user-defined chromosome */
        if (++ap > argc)
          err_msg_die("-%c needs a chromosome identifier (%s..%s,%s,%s)\n",uopt,
                        any_chrno2str(rpp->maxchrno,rpp->havexy,1,1),
                        any_chrno2str(rpp->maxchrno,rpp->havexy,rpp->maxchrno-2,1),
                        any_chrno2str(rpp->maxchrno,rpp->havexy,rpp->maxchrno-1,1),
                        any_chrno2str(rpp->maxchrno,rpp->havexy,rpp->maxchrno,1));
        else
          if ((rpp->uchrno = any_str2chrno(rpp->maxchrno,rpp->havexy,argv[ap])) == Chr_unk)
            err_msg_die("Can't determine Chromosome '%s'\n",argv[ap]);
        break;
      case 'C':     /* number columns of chr/pos file to scan/include */
        if (++ap > argc)
          err_msg_die("-%c needs column number (min 3)\n",uopt);
        else
          if ((rpp->srccollmt = (int) strtol(argv[ap],NULL,10)) < 3)
            err_msg_die("Can't work with less than 3 columns\n");
        break;
      case 'P':      /* require CDS for a gene to be considered */
        rpp->needcds = 1;
        break;
      case 'I':      /* allow internal fragments */
        rpp->intrnlfrag = IGL_intrnl_allow;
        break;
      case 'i':      /* internal fragments & relate to exon boundaries */
        rpp->intrnlfrag = IGL_intrnl_resolv;
        break;
      case 'K':      /* classify frag & feature positions */
        rpp->classolap = 1;
        break;
      case 't':      /* look for TSS */
        rpp->seektss = 1;
        (void) db_appnfelt(&rpp->tsslist,FTKW_TSS);
        break;
      case 'U':      /* look for CpG Is closest upstream */
        rpp->seekCpGis = IGL_CpG_upstream;
        (void) db_appnfelt(&rpp->tsslist,FTKW_TSS);
        break;
      case 'u':      /* CpG Is nearest to fragment */
        rpp->seekCpGis = IGL_CpG_proximal;
        (void) db_appnfelt(&rpp->tsslist,FTKW_CpG_is);
        break;
      case 'd':     /* limit distance */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",uopt);
        else
          rpp->distlimit = (int) strtol(argv[ap],NULL,10);
        break;
      case 'D':     /* limit TSS distance */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",uopt);
        else
          rpp->tssdistlmt = (int) strtol(argv[ap],NULL,10);
        break;
      case 'X':         /* include Xref & production info after name */
        rpp->displinfo = 1;
        break;
      case 'R':
        rpp->tss_cpgi_rng = 1;
        break;
      case 'm':         /* take gene name from mRNA */ 
        rpp->geneexmrna = 1;
        break;
      case 'B':         /* only use biotype:protein_coding features */
        if (++ap > argc)
          err_msg_die("-%c needs biotype string value\n",uopt);
        else
          (void) igl_chkcommadelstrs(&rpp->wantedcodingtypes,argv[ap]);
        break;
      case 'A':       /* gtf/gff2 wanted attributes */
        if (++ap > argc)
          err_msg_die("-%c needs attribute string value(s)\n",uopt);
        else
          (void) igl_chkcommadelstrs(&rpp->gtfgeneidattr,argv[ap]);
        break;
      case 'a':      /* user-specified feature type to override defaults */
        if (++ap > argc)
          err_msg_die("-%c needs feature type string value(s)\n",uopt);
        else
          (void) igl_chkcommadelstrs(&rpp->ufeatlist,argv[ap]);
        break;
      case 'n':       /* fragment based CpGI/TSS distances */
        rpp->omode = IGL_omod_fragbased;
        break;
      case 'N':      /* count genes/fragment */
        rpp->omode = IGL_omod_genecnt;
        break;
      case 'q':   /* gene id from CDS/protein_id field */
        rpp->cds_prot_id = 1;
        break;
      case 'W':
        rpp->printnullhits = 1;
        break;
      case 'o':         /* change output file */
        if (++ap >= argc)
          err_msg_die("-%c needs output file name\n",uopt);
        else
          {
          if ((rpp->outfile = fopen(argv[ap],"w")) == NULL)
            err_msg_die("Can't open output file '%s'\n",argv[ap]);
          }
        break;
      case 'b':      /* src file bufferlenGth */
        if (++ap > argc)
          err_msg_die("-%c needs buffer length value\n",uopt);
        else
          rpp->srcbuflen = (int) strtol(argv[ap],NULL,10);
        break;
      case 'h':
        igl_sayusage(rpp->outfile,argv[0]);
        exit(0);
        break;
      default:
        fprintf(stderr,"Unknown option: '-%c'\n",uopt);
        igl_sayusage(stderr,argv[0]);
        exit(1);
        break;
      }
/* sanity check for necessary parameters */
if ((rpp->ftflnam == NULL) && (rpp->ftprefix == NULL) && (rpp->ftsuffix == NULL)
      && (rpp->chrinfoflnam == NULL))
  err_msg_die("No file spec. for feature information, -f, -G or -p & -s values needed\n");
switch (rpp->dfmt)
  {
  case DBFMT_unknown:
    err_msg_die("No valid feature table style set (-e,-g,-F,-T or -Q)\n");
    break;
  case DBFMT_gtf:
    if (rpp->gtfgeneidattr == NULL)
      (void) igl_chkcommadelstrs(&rpp->gtfgeneidattr,"transcript_name,transcript_id,gene_name");
    if (rpp->ufeatlist == NULL)
      (void) db_scnfeatwrds("exon,CDS,start_codon,stop_codon,transcript,lncRNA",rpp->ftkw,&rpp->ftlist,0);
    else
      {
      sfp = rpp->ufeatlist;
      while (sfp != NULL)
        {
        if ((feattype = wlu_chkwrd(rpp->ftkw,sfp->strval)) == FTKW_unknown)
          err_msg_die("Error: %s not suitable for GTF feature type\n",sfp->strval);
        else
          db_appnfelt(&rpp->ftlist,feattype);
        sfp = sfp->nxtselt;
        }
      }
    if (rpp->wantedcodingtypes == NULL)  /* not specified, default to "protein_coding" for gtf */
      (void) igl_chkcommadelstrs(&rpp->wantedcodingtypes,"protein_coding");
    (void) igl_chkbiotypes4fmt(&rpp->wantedcodingtypes,rpp->dfmt,&rpp->wantedtypewlu);
    break;
  case DBFMT_embl:
  case DBFMT_genbank:
  case DBFMT_embl_old:
    (void) db_scnfeatwrds("gene,mRNA,CDS",rpp->ftkw,&rpp->ftlist,0);
    rpp->ftqul = db_getftqualstrct();
    fep = rpp->ftlist;
    while (fep != NULL)
      {
      (void) db_scnqulwrds("db_xref,name,gene",rpp->ftqul,&fep->wquals,0);
      if ((rpp->cds_prot_id) && (fep->fval == FTKW_CDS))
        (void) db_scnqulwrds("protein_id",rpp->ftqul,&fep->wquals,0);
      if ((rpp->displinfo) && (fep->fval == FTKW_CDS))
        (void) db_scnqulwrds("product",rpp->ftqul,&fep->wquals,0);
      fep = fep->nxtfelt;
      }
    break;
  case DBFMT_sqmonk:
    (void) db_scnfeatwrds("gene,mRNA,CDS,CpG,TSS",rpp->ftkw,&rpp->ftlist,0);
    rpp->ftqul = db_getftqualstrct();
    fep = rpp->ftlist;
    while (fep != NULL)
      {
      (void) db_scnqulwrds("db_xref,name,biotype",rpp->ftqul,&fep->wquals,0);
      fep = fep->nxtfelt;
      }
    (void) igl_chkbiotypes4fmt(&rpp->wantedcodingtypes,rpp->dfmt,&rpp->wantedtypewlu);
    break;
  case DBFMT_gff3:
    (void) db_scnfeatwrds("gene,mRNA,CDS",rpp->ftkw,&rpp->ftlist,0);
    rpp->ftqul = wlu_getlustrct(WLU_CASEIND,-1);
    wlu_addwrd(rpp->ftqul,"Dbxref",1,NULL);
    wlu_addwrd(rpp->ftqul,"ID",2,NULL);
    wlu_addwrd(rpp->ftqul,"Is_circular",3,NULL);
    wlu_addwrd(rpp->ftqul,"Name",4,NULL);
    wlu_addwrd(rpp->ftqul,"Note",5,NULL);
    wlu_addwrd(rpp->ftqul,"Parent",6,NULL);
    wlu_addwrd(rpp->ftqul,"Target",7,NULL);
    wlu_addwrd(rpp->ftqul,"anticodon",8,NULL);
    wlu_addwrd(rpp->ftqul,"chromosome",9,NULL);
    wlu_addwrd(rpp->ftqul,"codons",10,NULL);
    wlu_addwrd(rpp->ftqul,"description",11,NULL);
    wlu_addwrd(rpp->ftqul,"dev-stage",12,NULL);
    wlu_addwrd(rpp->ftqul,"exception",13,NULL);
    wlu_addwrd(rpp->ftqul,"exon_number",14,NULL);
    wlu_addwrd(rpp->ftqul,"gbkey",15,NULL);
    wlu_addwrd(rpp->ftqul,"gene",16,NULL);
    wlu_addwrd(rpp->ftqul,"gene_synonym",17,NULL);
    wlu_addwrd(rpp->ftqul,"genome",18,NULL);
    wlu_addwrd(rpp->ftqul,"linkage-group",19,NULL);
    wlu_addwrd(rpp->ftqul,"mol_type",20,NULL);
    wlu_addwrd(rpp->ftqul,"ncrna_class",21,NULL);
    wlu_addwrd(rpp->ftqul,"partial",22,NULL);
    wlu_addwrd(rpp->ftqul,"product",23,NULL);
    wlu_addwrd(rpp->ftqul,"protein_id",24,NULL);
    wlu_addwrd(rpp->ftqul,"pseudo",25,NULL);
    wlu_addwrd(rpp->ftqul,"transcript_id",26,NULL);
    fep = rpp->ftlist;
    while (fep != NULL)
      {
      (void) db_scnqulwrds("Dbxref,Name,gene",rpp->ftqul,&fep->wquals,0);
      fep = fep->nxtfelt;
      }
    break;
  default:
    (void) db_scnfeatwrds("gene,mRNA,CDS",rpp->ftkw,&rpp->ftlist,0);
    break;
  }
if (igl_maknamsnparsftfls(rpp) > 0)
  {
  switch (rpp->rstyle)
    {
    case IGL_listfeatdata:
      switch (rpp->dfmt)
        {
        case DBFMT_gtf:
          curseq = rpp->meplist;
          while (curseq != NULL)
            {
            fp = curseq->me_entstrptr->featlst;
            while (fp != NULL)
              {
              fputs(curseq->me_id,rpp->outfile); 
              tg_putgfffeatrecurs(rpp->outfile,fp,1);
              fp = fp->nextst;
              }
            curseq = curseq->nxt_dme;
            }      
          break;
        case DBFMT_gff3:
          curseq = rpp->meplist;
          while (curseq != NULL)
            {
            fp = curseq->me_entstrptr->featlst;
            while (fp != NULL)
              {
              if ((fp->featur == FTKW_CDS) || (fp->featur == FTKW_FT_gene))
                {
                fputs(curseq->me_id,rpp->outfile); 
                tg_putgfffeatrecurs(rpp->outfile,fp,1);
                }
              fp = fp->nextst;
              }
            curseq = curseq->nxt_dme;
            }      
          break;
        default:
          curseq = rpp->meplist;
          while (curseq != NULL)
            {
            fprintf(rpp->outfile,"%s: %d items from %d\n",curseq->me_id,curseq->me_entstrptr->nfeats,curseq->me_entstrptr->nseen);
            fp = curseq->me_entstrptr->featlst;
            while (fp != NULL)
              {
              fprintf(rpp->outfile,"  %d: %s\n",fp->savdno,db_ftkw2str(fp->featur));
              tg_putgfffeatrecurs(rpp->outfile,fp,1);
              fp = fp->nextst;
              }
            curseq = curseq->nxt_dme;
            }      
          break;
        }
      break;
    case IGL_listchrids:
      curseq = rpp->meplist;
      while (curseq != NULL)
        {
        fprintf(rpp->outfile,"%s\n",curseq->me_id);
        curseq = curseq->nxt_dme;
        }
      break;
    case IGL_findproximal:
      if (src == NULL)
        err_msg_die("proximal gene finding needs CpG position info\n");
     else
        igl_procchrposfile(src,rpp->outfile,rpp);
      break;
    }
  exit(0);
  }
else
  exit(1);
}
