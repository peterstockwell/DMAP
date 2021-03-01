/* sqfl_fns.h: datatype and routine definitions for c sequence file functions
 */

typedef enum SFMT_type        /* seq file format types */
  {
  SFMT_undefined = 0,
  SFMT_staden,
  SFMT_molgen,
  SFMT_nbrf,
  SFMT_fasta,
  SFMT_gcg,
  SFMT_raw
  }
SFMT_TYPE;

  /* max chars in header buffer */
#define HDRBUF_MAX 512

typedef enum sqfl_sqtopol  /* sequence topology */
  {
  SQTP_unknown = 0,
  SQTP_linear,
  SQTP_circular
  }
SQFL_SQTOPOL;

typedef enum sqfl_sqtos    /* type of sequence */
  {
  SQFL_tosunknown = 0,
  SQFL_dna,
  SQFL_rna,
  SQFL_peptide
  }
SQFL_SQTOS;

typedef enum SQ_restype    /* enumerated set of sequence residues */
  {
  RES_x = 0, /* unknown,unspecified */
  RES_a,     /* A or Ala */
  RES_c,     /* C or Cys */
  RES_g,     /* G or Gly */
  RES_t,     /* T,U or Thr */
  RES_b,     /* Asx */
  RES_d,     /* Asp */
  RES_e,     /* Glu */
  RES_f,     /* Phe */
  RES_h,     /* His */
  RES_i,     /* Ile */
  RES_k,     /* Lys */
  RES_l,     /* Leu */
  RES_m,     /* Met */
  RES_n,     /* Asn */
  RES_p,     /* Pro */
  RES_q,     /* Gln */
  RES_r,     /* Arg */
  RES_s,     /* Ser */
  RES_v,     /* Val */
  RES_w,     /* Trp */
  RES_y,     /* Tyr */
  RES_z,     /* Glx */
  RES_blk    /* null character */
  }
SQ_RESTYPE;

#define MAXSQCHRVAL 128    /* max byte value for any valid seq character */

typedef struct sqfl_strct  /* information relevant to sequence files */
  {
  SFMT_TYPE flfmt;
  FILE *sfl;               /* stream for sequence file */
  char *filnam;            /* sequence file name */
  char *seqnam;            /* sequence name */
  SQFL_SQTOPOL stopol;     
  SQFL_SQTOS sqtos;
  int slength;             /* length of sequence */
  int gcgchksum;           /* checksum for gcg output */
  int scnt;                /* count of written residues */
  char lstinchr;           /* cache the last char read */
  char *gnam;              /* name of auxfile for gcg temp storage */
  FILE *gfl;               /* aux file for gcg temporary storage */
  int lnlen;               /* length of output lines */
  char *annot;             /* annotation - use def if NULL */
  int abuflen;             /* length of allocated buffer */
  char *thisuse;           /* "r" or "w": fopen usage strings */
  int validlut[MAXSQCHRVAL]; /* lut for valid residues */
  }
SQFL_STRCT;

  /* chars in sequence name */
/* #define SQNAMEMAXLEN 12 */
#define SQNAMEMAXLEN 80

/* allow longer names to avoid downstream problems: Jan-2005 */

/* GCG uses this value to calculate checksum positions */

#define SQFL_GCGCHK 57

typedef enum sqfl_gcgwrds  /* useful header words for gcg files */
  {
  SQFL_gcgunknown = 0,
  SQFL_gcglength,
  SQFL_gcgtype,
  SQFL_gcgdots            /* .. data indicator */
  }
SQFL_GCGWRDS;

/* function headers */

SFMT_TYPE sqfl_chr2fmttp(char fchr);
  /* convert the letters s,q,f,n,u to corresponding file format type.
  currently case independent */

SFMT_TYPE sqfl_getdeffmt();
  /* use env variable SEQDEFFILEFMT to establish default file type */

char *sqfl_fmttp2strng(SFMT_TYPE sfmt);
  /* return the pointer to a string defining sfmt */

char sqfl_fmttp2chr(SFMT_TYPE sfmt);
  /* return a character sfmt */

char sqfl_restype2chr(SQ_RESTYPE rt);
  /* return the character corresponding to rt */

char *sqfl_restype2str(SQ_RESTYPE rt);
  /* return the address of a string corresponding to rt */

SQ_RESTYPE sqfl_chr2restype(char res);
  /* return an internal restype for this character */

SQ_RESTYPE sqfl_chr2narestype(char rs);
  /* return RES_x..RES_t only */

SQ_RESTYPE sqfl_chr2aarestype(char r);
  /* return the corresponding restype for r, assuming it is a valid AA code,
else return RES_x */

SQFL_SQTOS sqfl_chr2tos(SFMT_TYPE sfmt,
                        char tc);
/* return the tos for tc - case independent */

int sqfl_linelength(SFMT_TYPE sfmt);
  /* return normal linelength used by this format */

int sqfl_gcgchkinc(int sp,
                   char res);
/* return the incremental addition for the GCG check algorithm for res at
  position sp, (1..n) */

int sqfl_gcgsqchksum(char *sq,
                     int sqlen);
/* return gcg checksum for (1..sqlen) chars of sq.  Stop at '\0'
if encountered ealier */

SQFL_SQTOS sqfl_scantos(char *sqbuf,
                        int sqln);
/* scan sqbuf, looking for evidence of NA or peptide.
  Algorithm looks for the letters E,F,I,O,X,Z: none of which should occur in
  DNA sequences, even using IUB or LDNA redundant base codes.  A proportion
  1% of strangers is allowed */

SQFL_SQTOS sqfl_scanfltos(SQFL_STRCT *src);
  /* scan contents of previously opened src, looking for evidence of 
NA or peptide.  Rewind source file on completion.
  Algorithm looks for the letters E,F,I,O,X,Z: none of which should occur in
  DNA sequences, even using IUB or LDNA redundant base codes.  A proportion
  1% of strangers is allowed */

char sqfl_tos2chr(SFMT_TYPE sfmt,
                  SQFL_SQTOS st);
  /* return the character for sequence st - default to D */

char sqfl_topol2chr(SFMT_TYPE sfmt,
                    SQFL_SQTOPOL st);
/* return a character corresponding to st for sfmt */

void sqfl_setsqdetails(SQFL_STRCT *st,
                       char *sqbuf,
                       int sqlen);
/* set as many as possible of st details by examination of sqbuf */

void sqfl_headsfstr(SQFL_STRCT *sst,
                    char *sqname,
                    char *ann);
/* write sequence header info to previously opened sst->sfl, using string ann
 */

void sqfl_headsfstrct(SQFL_STRCT *sst,
                      char *sqname,
                      char *maker,
                      char *origin);
/* write sequence header info to previously opened sst->sfl */

void sqfl_headsfstrctann(SQFL_STRCT *sst,
                         char *sqname);
/* write sequence header info to previously opened sst->sfl. use annotation
  if set */

void sqfl_newlnstart(FILE *sfl,
                     SFMT_TYPE sfmt,
                     int sp);
/* if necessary, write this number as a line start */

void sqfl_putchr(SQFL_STRCT *sqs,
                 char chr,
                 int rp,
                 int *lcnt);
/* put chr out to file, ignoring validity. observe line count etc */

void sqfl_putres(SQFL_STRCT *sqs,
                 char res,
                 int rp,
                 int *lcnt);
/* put res out to file, keeping track of valid bases, line count etc */

void sqfl_termsqfl(SQFL_STRCT *sqs,
                   char *pname,
                   int *lcnt);
/* finish off sequence if necessary */

int sqfl_lookforchr(FILE *sfl,
                    char ec);
/* scan sfl, looking for ec.  stop and return 1.  Return 0 if EOF encountered
  before hand */

int sqfl_look4chrs(FILE *sfl,
                   char *ecs);
/* scan sfl, looking for anything in ecs.  stop and return that char.
  Return EOF if not encountered by EOF */

void sqfl_skipeol(SQFL_STRCT *sqs);
  /* try to skip to end of current line - check if just seen eol already */

void sqfl_forceeol(SQFL_STRCT *sqs);
  /* force skip to eol irrespective of what has been seen todate */

int sqfl_cpy2eol(SQFL_STRCT *sqs,
                 char *buf,
                 int bufmax,
                 int (*ufgetc)(FILE *fl));
/* copy to next line end to buf, up to bufmax, using ufgetc,
  return no of chars written */

int sqfl_skipsqflhdr(SQFL_STRCT *sqs);
  /* skip file header, getting to start of sequence.  Eventually this routine
should be able to extract data from the header.  Return 1 if nothing untoward
occured */

char *sqfl_valid4fmt(SFMT_TYPE fmt);
  /* return pointer to string of valid residues for fmt type of file.
  This is written as case-independent - uppercase only being shown */

char sqfl_normfillc(SFMT_TYPE sf);
  /* return the normal fill char used for this format */

char sqfl_getnxtres(SQFL_STRCT *sqs);
  /* return the next valid residue from sqs.  EOF or non-alphanumeric
   return NULL.  invalid chars are skipped */

int loadsrcrng(SQFL_STRCT *sqs,
               char *seqbuf,
               int fstrt,
               int fstop);
/* read successive chars from sfl, write valid bases in fstrt..fstop into
  seqbuf if it is not NULL.  fstrt&fstop are sequence positions (1..n)
  Don't put in terminal null char.
  return no of chars that would be inserted (whether stored or not) */

int loadsrcsq(SQFL_STRCT *sqs,
              char *seqbuf);
/* read successive chars from sfl, write valid bases into seqbuf if it is 
  not NULL.  Don't put in terminal null char.
  return highest position read (whether stored or not) */

SQFL_STRCT *sqfl_creatsqstrctann(FILE *ufl,
                                 char *fnam,
                                 SFMT_TYPE sfmt,
                                 char *myuse,
                                 char *annotat);
/* use an open sequence file fnam, returning pointer to a new structure
  if successful.  NULL if not.
  annotat is pointer to annotation string (most useful for output files) to
  be used if non-NULL.  myuse is as for fopen */

SQFL_STRCT *sqfl_opnsqstrctann(char *fnam,
                               SFMT_TYPE sfmt,
                               char *myuse,
                               char *annotat);
/* attempt to open a sequence file fnam, returning pointer to a new structure
  if successful.  NULL if not.  myuse is a normal parameter to fopen, though
  it is likely that "r" and "w" are the most useful of them,
  annotat is pointer to annotation string (most useful for output files) to
  be used if non-NULL */

SQFL_STRCT *sqfl_opnsqstrct(char *fnam,
                            SFMT_TYPE sfmt,
                            char *myuse);
/* attempt to open a sequence file fnam, returning pointer to a new structure
  if successful.  NULL if not.  myuse is a normal parameter to fopen, though
  it is likely that "r" and "w" are the most useful of them */

SQFL_STRCT *sqfl_opntoread(char *fnam,
                           SFMT_TYPE sfmt);
/* attempt to open a sequence file fnam for reading, returning pointer to a 
  new structure if successful, NULL if not.  Perform header skipping, so
  seq can be read directly */

void sqfl_rewind(SQFL_STRCT *sqs);

void sqfl_clssqstrct(SQFL_STRCT *sqs);
  /* close this sequence data structure, freeing the allocated memory */

void sqfl_clsstrctfil(SQFL_STRCT *sqs);
  /* close the file component of *sqs, retaining other information */

off_t sqfl_filelength(FILE *fl);
  /* use fstat to return the length of file fl in bytes */

int readsrcsq(SQFL_STRCT *sqs,
              char *seqbuf);
/* read successive chars from sqs, write valid bases into seqbuf if it is 
  not NULL.
  return highest position read (whether stored or not) */

char *sqfl_tos2strng(SQFL_SQTOS st);
  /* return pointer to character string defining st */

SFMT_TYPE sqfl_prmpt4sfmt(SFMT_TYPE def,
                          SFMT_TYPE rngex,
                          SFMT_TYPE rngto);
/* prompt for a sequence file format, offer default def, if defined */

char *sqfl_defextnsns(SFMT_TYPE sf);
  /* return a suggested extension for type sf */

void sqfl_sfmts2bf(char *ubuf,
                   int blim);
/* write sfmt details to ubuf */

int sqfl_wrtsq2fl(SQFL_STRCT *sfl,
                  char *ubuf,
                  int blen,
                  int *lcnt);
/* put blen chars of ubuf as sequence to sfl file.  Return res count */

int sqfl_wrtsqbuf2fl(SQFL_STRCT *sfl,
                     char *sqnam,
                     char *maker,
                     char *origin,
                     char *ubuf,
                     int blen);
/* write ubuf up to blen out to sfl, (previously opened).  Do headers & 
terminators.  Don't close sfl, leave to calling routine.  Return no res
written */

int sqfl_fgetc(FILE *fl);
  /* return character (as int) from fl.  ignore '\r' chars totally */
