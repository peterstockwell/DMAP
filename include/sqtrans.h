/* sqtrans.h: headers for sequence translation routines */

#define CODONLENGTH 3
#define MAXBASINT 4
#define MAXAAINT 21
#define SQT_TERMCHR '*'
#define MAX_CODON_NO 64

typedef enum sqt_gencod    /* different genetic codes - genbank numbering */
  {
  GC_undefined = 0,
  GC_universal,            /* standard */
  GC_vertmito,             /* vertebrate mitochondrial */
  GC_yeastmito,            /* yeast mitochondrial */
  GC_mycomito,      /* mold, protozoan, coelent, mycoplasma mitochondrial */
  GC_invertmito,           /* invertebrae mitochondrial */
  GC_ciliate,              /* Ciliate Dasycladacean Hexamita Nuclear */
  GC_echinmito = 9,        /* Echinoderm mitochondrial */
  GC_euplotid,             /* Euplotid nuclear */
  GC_bacterial,            /* Bacterial */
  GC_yeastalt,             /* Alternative Yeast Nuclear */
  GC_ascidmito,            /* Ascidian mitochondrial */
  GC_flatwmito,            /* Flatworm mitochrondiral */
  GC_blepharisma,          /* Blepharisma nuclear */
  GC_userdefined
  }
SQT_GENCOD;

typedef struct
  {
  char reschr;    /* character for translation */
  int trinit;     /* initiator */
  } TR_RESIDUE;

typedef union
  {
  TR_RESIDUE trres;
  int iresidue;
  float xresidue;
  } RESIDUEDATA;

typedef enum sqtrn_datype  /* type of data stored in matrix */
  {
  SQTDT_undefined = -1,
  SQTDT_char,
  SQTDT_int,
  SQTDT_float
  }
SQTRN_DATYPE;

typedef struct
  {
  SQTRN_DATYPE dtype;
  RESIDUEDATA resdata;
  } TRN_MATELT;

typedef TRN_MATELT TRANS_MATRX [MAXBASINT+1] [MAXBASINT+1] [MAXBASINT+1];

typedef int CODN_VECTR [CODONLENGTH];  /* array of integers for a codon */

typedef struct cod_vect_strct
  {
  CODN_VECTR codnvector;    /* vector for this codon */
  struct cod_vect_strct *nextcvec; /* next such struct */
  }
COD_VECT_STRCT;

/* typedef COD_VECT_STRCT *CVEC_PTRS [MAXAAINT + 1]; */

/* routine headers for sequence translation and genetic code 
     manipulation */

char int2nabas(int iv);
/* return base corresponding to iv: ?=0 A=1 C=2 G=3 T=4 */

SQT_GENCOD sqt_int2gbgencod(int iv);
  /* return the corresponding symbolic name for this integer */

char *sqt_gbgencd2strng(SQT_GENCOD gc);
  /* return a string pointer for this gc value */

void sqt_saygbgencodtbls(FILE *fl,
                         int from,
                         int to);
/* describe valid codes in range from..to at fl */

int sqt_nabas2int(char abase);
  /* return integer corresponding to abase: 0=? 1=A 2=C 3=G 4=T/U */

char sqt_int2nabas(int iv);
  /* return base corresponding to iv: ?=0 A=1 C=2 G=3 T=4 */

int remapbasptr(int bp);
  /* remap base pointers to return bases in order TCAG, as per conventional
  genetic code display */

int aares2int(char aares);
/* return an appropriate integer pointer for this residue */

char *aachr2str3(char aares);
/* return an appropriate 3 letter string for this residue */

int sqt_res2aaint(char ares);
  /* return an integer for this residue - same as aares2int, except that
     unknown = 0 and output is 0..MAXAAINT */

char int2aares(int aaint);
  /* return character for this integer pointer */

void fill_chr_trnarr(TRANS_MATRX trnsar,
                     char fillch);
/* fill trnsarr with fill character */

void fill_ini_trnarr(TRANS_MATRX trnsar,
                     int fillini);
/* fill trnsarr with fill character */

void fill_int_trnarr(TRANS_MATRX trnsar,
                     int fillval);
/* fill trnsarr with fill value */

void fill_flt_trnarr(TRANS_MATRX trnsar,
                     float fillval);
/* fill float trnsarr with fill value */


int set_chr_trnarr(TRANS_MATRX trnsar,
                   char *cdn,
                   char res);
/* put char into relevent codon positions data struct. Wild card *
  causes all nonzero values to be initialized */

int set_ini_trnarr(TRANS_MATRX trnsar,
                   char *cdn,
                   int rini);
/* put rini into relevent codon positions data struct. Wild card *
  causes all nonzero values to be initialized */

int set_int_trnarr(TRANS_MATRX trnsar,
                   char *cdn,
                   int ival);
/* put ival into relevent codon positions data struct. Wild card *
  causes all nonzero values to be initialized */

int set_flt_trnarr(TRANS_MATRX trnsar,
                   char *cdn,
                   float xval);
/* put xval into relevent codon positions data struct. Wild card *
  causes all nonzero values to be initialized */

TRN_MATELT *sqt_vect2eltaddr(TRANS_MATRX tm,
                             CODN_VECTR cv);
/* return the address of the array element corresponding to cv */

int sqt_cdn2codvect(char *cdn,
                    CODN_VECTR cv);
/* convert cdn to cv. return 1 if all bases OK, else 0 */

int sqt_codvect2cdn(CODN_VECTR cv,
                    char *cdn);
/* convert cv to corresponding codon in cdn (presumed to be big enough)
  return 1 if all of codon is valid, 0 otherwise */

int sqt_cdnno2cdn(int cno,
                  char *cdnbuf);
/* write to cdnbuf (assumed to be
CODONLENGTH+1 chars or longer) he codon corresponding to
enumerated codon cno (0..63). return 1 if OK, else 0 */

TRN_MATELT *sqt_cdnno2codvect(TRANS_MATRX tmat,
                              int cno,
                              char *ucdn);
/* if cno is valid, return the vector pointer
for it for translation table tmat.  NULL if
can't.  If ucdn is non-NULL, write the codon
to it (assumes enough room) */

char sqt_trnslate(TRANS_MATRX trnsar,
                  char *cdn);
/* translate cdn according to trnsar */

int sqt_get_int4tmatrix(TRANS_MATRX trnsar,
                        char *cdn);
/* return integer value for  cdn according to trnsar.
return 0 if not integer datum */

float sqt_get_flt4tmatrix(TRANS_MATRX trnsar,
                          char *cdn);
/* return the float value for translate cdn according to trnsar.
0.0 if not float datum */

void sqt_inc_cdncnt(TRANS_MATRX trnsar,
                    char *cdn);
/* increment iresidue for cdn */

void init_univ(TRANS_MATRX trnsar);
  /* fill trnsar with codes for universal genetic code */

void init_trnsar4gb(TRANS_MATRX trnsar,
                    SQT_GENCOD gencod);
/* initialize trnsar to any known genetic code - using Genbank trans-table
  values */

char init_trnsar4chr(TRANS_MATRX trnsar,
                     char gencod);
/* initialize trnsar to any known genetic code, according to following:
u = universal, m = mammalian mitochondrial, y = yeast mitochondrial,
a = aspergillus/Neurospora and Trypansoma mitochondrial,
d = drosophila mitochondrial

return character, if valid, else null. */

int uniquecdn(char *cdn);
  /* return a unique nonzero numeric score for this codon.
  only if all elements are A,C,G,T.  Else return 0 */

int validcdn(char *cdn);
  /* return true if this codon is valid.  only if all elements
   are A,C,G,T, or *. */

int sqt_getnxtcdn(char *src,
                  int slen,
                  int *sqpos,
                  char *cdn);
/* fill cdn from sqpos of src, return no of chars written.  
 Stop on null byte or when slen is reached.
 cdn should allow room for null byte.
 *sqpos is 0..n-1 */

int sqt_transbuffer(char *src,
                    int slen,
                    TRANS_MATRX gencod,
                    char *dst);
/* translate sequence in src to dst (which will be assumed to be long enough.
  stop on encountering null or slen has been translated */

void sqt_inivectarr(COD_VECT_STRCT **cvs);
  /* initialise all cvs aa vectors to null */

COD_VECT_STRCT *sqt_appnd_cvect(COD_VECT_STRCT **cvp,
                                CODN_VECTR cv);
/* append this vector to *cvp.  Return a pointer to the new component, NULL
 if problem */

COD_VECT_STRCT *sqt_appndcdn2vec(COD_VECT_STRCT **cvp,
                                 char *cdn);
/* append a vector for cdn to *cvp assuming it is a valid codon.
  Return a pointer to the new component, NULL if problem */

int sqt_bldvects4matrx(TRANS_MATRX tmtx,
                       COD_VECT_STRCT **cvecs);
/* build a set of vectors in cvecs for the translation matrix tmtx.
return the number of vectors added */

void sqt_clrvectarr(COD_VECT_STRCT **cvarr);
  /* clear all of vectors for all aas */

COD_VECT_STRCT *sqt_getvect4res(COD_VECT_STRCT **cvlist,
                                char aares,
                                COD_VECT_STRCT *cvp);
/* if cvp = NULL, then return the first vector in the linked list for residue
aares.  if cvp is non-null, return the next codon vector in the linked list */

COD_VECT_STRCT *sqt_ptr4cvec(COD_VECT_STRCT *clst,
                             CODN_VECTR cv);
/* return a pointer to the element of clst which matches cv.  NULL if none */

int sqt_codnbas4res(COD_VECT_STRCT **cvlist,
                    char aares,
                    int cpos,
                    char *cdn,
                    char *cstrng);
/* return a number of base choices for position cpos (1..3) of the
codons for aares.  Return bases in cstrng if non-NULL, 
assumed to be large enough for the purpose.
Error returns: -2 = cpos ivalid, -1 = invalid aares.
If the codon seen so far is present as cdn, then this information is used in
refining the output (Esp for Ser, etc in the standard genetic code */

COD_VECT_STRCT **sqt_getcvctarr(char *msg);
  /* allocate storage for a MAXAAINT-long array of pointers */

void sqt_sumintmatrx(TRANS_MATRX m1,
                     TRANS_MATRX m2,
                     TRANS_MATRX sum);
/* sum each element of m1 and m2 into sum */

int sqt_totvldicnts(TRANS_MATRX tmat);
  /* for all valid positions in tmat, sum the integer counts */

int sqt_gcodexfl(FILE *fl,
                 TRANS_MATRX tmat);
/* read CODONS/RES from fl and load into tmat, return no of valid codons read.
fields starting with ';' are regarded as comments */

void say_intcdn(FILE *strm,
                int p1,
                int p2,
                int p3);
/* tell strm the character representation of the codon (p1,p2,p3) */

void sqt_gcod2fl(FILE *fl,
                 char *desc,
                 TRANS_MATRX tmat,
                 int letter3);
/* write tmat contents to fl in approved form.  desc written as leading
comment if non-NULL, letter3 enables 3 letter AA ids */

COD_VECT_STRCT *sqt_maxintval4cvec(TRANS_MATRX tm,
                                   COD_VECT_STRCT *c1,
                                   COD_VECT_STRCT *c2);
/* compare integer value for tm values of c1,c2, returning
the larger.  if c1 or c2 is NULL, return the other. NULL
if cockup */

COD_VECT_STRCT *sqt_minintval4cvec(TRANS_MATRX tm,
                                   COD_VECT_STRCT *c1,
                                   COD_VECT_STRCT *c2);
/* compare integer value for tm values of c1,c2, returning
the larger.  if c1 or c2 is NULL, return the other. NULL
if cockup */

COD_VECT_STRCT *sqt_scnveclst(COD_VECT_STRCT *clst,
                              TRANS_MATRX tm,
                              COD_VECT_STRCT *(*cmpfun)(TRANS_MATRX tx,
                                                        COD_VECT_STRCT *c1,
                                                        COD_VECT_STRCT *c2));
/* Scan thru clst... using tm and cmpfun to decide on best
vector in list.  cmpfun will return which ever of c1,c2 is "better"
*/

int sqt_vectrsequal(COD_VECT_STRCT *c1,
                    COD_VECT_STRCT *c2);
/* return 1 if c1 & c2 are the same */

int sqt_tmatival4cvecp(TRANS_MATRX tm,
                       COD_VECT_STRCT *cvp,
                       SQTRN_DATYPE *dtp);
/* simply return the integer value for the
cvp postion of tm.  0 if not an integer value.
if dtp is non-NULL, then return the data type
to it */
