/* bam_fns.h: headers and data structures for a simple
callable set of routines to return BAM file material
in defined structures to the caller.  Intended use is
to permit reading of BAM file alignments within my
own SW without requiring bam->sam conversion of the
files first.

Peter Stockwell: Apr-2015 */
/* PAS: Correction for compression block alignment and read structures: Nov-2020 */

/* for a description of the BAM data file format see:

https://github.com/samtools/hts-specs

*/

#define BAM_FLAG_CMP 16

typedef struct bf_bgzf_hdr
  {
  uint8_t id1;
  uint8_t id2;
  uint8_t cm;
  uint8_t flg;
  uint32_t mtime;
  uint8_t xfl;
  uint8_t os;
  uint16_t xlen;
  }
BF_BGZF_HDR;

typedef struct bf_subfield
  {
  uint8_t si1;
  uint8_t si2;
  uint16_t slen;
  uint16_t bsize;
  }
BF_SUBFIELD;

typedef struct bf_datafield
  {
  int srcfdes;
  FILE *srcfl;        /* in case we need to close it this way */
  int blockno;
  z_stream *zstr;
  uint32_t crc32;
  uint32_t isize;
  char *cdata;
  int cdbuflen;
  int cdlen;
  int cdatp;
  char *udata;
  int udbuflen;
  int udatp;
  int blockremain;
  int isok;           /* for gdb watch checking */
  }
BF_DATAFIELD;

typedef struct bf_bamhdr_1
  {
  char magic[4];
  int32_t l_text;
  }
BF_BAMHDR_1;

typedef struct bf_binmqnl
  {
  uint8_t l_read_name: 8;
  uint8_t filler : 8;
  uint8_t MAPQ : 8;
  uint8_t bin : 8;
  }
BF_BINMQNL;

typedef struct bf_alignrec
  {
  int32_t block_size;
  int32_t refID;
  int32_t pos;
  uint8_t l_read_name: 8;
  uint8_t MAPQ : 8;
  uint16_t bin : 16;
/*  uint32_t bin_mq_nl; */
  uint16_t n_cigar_op : 16;
  uint16_t FLAG : 16;
/*   uint32_t flag_nc; */
  int32_t l_seq;
  int32_t next_refID;
  int32_t next_pos;
  int32_t tlen;
  }
BF_ALIGNREC;

typedef struct bf_cigar_op_type
  {
  int cigar_op : 4;
  int c_op_len : 28;
  }
BF_CIGAR_OP_TYPE;

typedef struct bf_aux_info       /* to contain data for auxilliary data */
  {
  char tag[2];
  char val_type;
  void *auxdat;
  int aival;
  float afval;
  char arrtype;  /* only needed for 'B' type arrays */
  int arrbtlen;  /* ditto: length of B arrays in bytes */
  struct bf_aux_info *nxtaux;
  struct bf_aux_info *prvaux;
  }       /* linked list for arbitrary expansion */
BF_AUX_INFO;

typedef struct bf_expndalignrec
  {
  BF_ALIGNREC *arec;
  char *read_name;
  BF_CIGAR_OP_TYPE *cigar;
  char *seq;
  char *qual;
  BF_AUX_INFO *auxinfo;
  }
BF_EXPNDALIGNREC;

typedef struct bf_refinfo
  {
  int l_name;
  char *refname;
  int l_ref;
  }
BF_REFINFO;

typedef struct bf_rundata /* complete kludge-up of everything */
  {
  BF_DATAFIELD *dfptr;
  int bhdrtxtlen;   /* length of next V */
  char *hdrtext;    /* malloc()ed text of header */
  int n_ref;        /* number ref sequences */
  BF_REFINFO *refinfo;  /* array of ref data */
  int bigendian;   /* 1 if we need to swap byte order */
  }
BF_RUNDATA;

/* function headers */

/* endian functions draw on samtools/bam definitions */

int bf_is_bigendian();
  /* return true if machine representation is BE */

uint16_t bf_endswap_uint16(uint16_t v);

uint32_t bf_endswap_uint32(uint32_t v);

void bf_err_msg_die(char *fmt,
                    ...);

void bf_err_msg(char *fmt,
                ...);
/* write user error message */

void bf_null_err_msg(char *fmt,
                     ...);
/* do nothing */

int bf_chk_bgzfhdr(BF_BGZF_HDR *hptr);
  /* check for valid  values in *hptr.

int bf_chk_subfield(BF_SUBFIELD *sfptr);
  /* check for expected fields in *sfptr
& if so return the buffer size, else 0 */

void bf_prt_bgzf_hdr(FILE *ofl,
                     BF_BGZF_HDR *hptr);

void bf_prt_subfield(FILE *ofl,
                     BF_SUBFIELD *sfp);


int bf_chkzlibstatus(char *op,
                     int z_ret,
                     void emsgfn(char *xfmt,
                                 ...));

int bf_chkbamfmt(char *ubuf);
/* check if ubuf contains the expected bam header info.
return the header text length, -1 if corrupt.  non-NULL
htextptr will be set to start of SAM header text */

char *bf_getrefinfo(char *buf,
                    int buflen,
                    int *refnamlen,
                    char **refnam,
                    int *refsqlen);
/* assume that buff points to the start
of a reference name block. checking that we
don't exceed the buffer size buflen, return
relevant information.  Return the final
position of the buffer after scanning the
data, NULL if can't do it */

void bf_fputanyc(char nc,
                 FILE *ofl);
/* try to meet C conventions for
printing any char */

void bf_strout(FILE *ofl,
               char *str,
               int slen);
/* write slen chars of str out to ofl,
mapping chars to octal, etc if necessary */

void bf_fputivals(uint32_t *intstream,
                  int n,
                  FILE *ofl);

int bf_getnewblock(BF_DATAFIELD *dfp);
  /* Assume that the file is ready to read the 
next bgzf block.  get and uncompress a new data block,
set parameters and return the uncompressed size */
  
void bf_initdatafield(BF_DATAFIELD *dfp,
                      FILE *srcfile);
/* make initial settings in dfp */

uint8_t bf_nxtuncompbyte(BF_DATAFIELD *dfp,
                         int *isok);
  /* return the next byte from the uncompressed
stream.  Redo block if necessary */

int bf_fill_structure(void *sptr,
                      size_t nbytes,
                      BF_DATAFIELD *dfp);
/* put nbytes of uncompressed data into sptr
from dfp.  Return number of bytes transferred */

int bf_fill_strct_blkrem(void *sptr,
                         size_t nbytes,
                         BF_DATAFIELD *dfp);
/* put nbytes of uncompressed data into sptr
from dfp while dfp->blockremain > 0.
Return number of bytes transferred */

char bf_int2sqchar(int ival);
  /* ival is a 4 bit compressed code
which should map to a residue.  Return
this */

char bf_int2cigarchar(int32_t ival);
  /* similar to above for CIGAR codes */

void bf_prtbfarec(BF_ALIGNREC *bfarecp,
                  FILE *ofl);
/* tell ofl about contents of bfarecp */

size_t bf_auxchr2byteno(char bc);
  /* number of bytes for char bc */

int bf_getbamauxint4hdr(BF_AUX_INFO *dhp,
                        BF_DATAFIELD *dfp,
                        int *isok);
/* depending on type of dhp->val_type, return
an integer value */

float bf_getbamauxflt(BF_DATAFIELD *dfp,
                      int *isok);
/* return a float: assumed that the
val_type char is f */

void bf_dispose_datafield(BF_DATAFIELD *dfp);
  /* lose all storage associated with dfp */

BF_DATAFIELD *bf_openinbamflmsg(char *bfname,
                                void err_fn(char *msg,
                                              ...));
  /* attempt to open bamfile bfname for input.
return the BF_DATAFIELD initialised if successful
otherwise NULL. err_fn used to communicate error
issues */

BF_DATAFIELD *bf_openinbamfile(char *bfname);
  /* attempt to open bamfile bfname for input.
return the BF_DATAFIELD initialised if successful
otherwise NULL */

void bf_initrundata(BF_RUNDATA *rdp);
  /* init *rdp to default settings */

void bf_disposerundata(BF_RUNDATA *rdp);
  /* do the obvious */

BF_RUNDATA *bf_opnchkbamflmsg(char *bfname,
                              void err_fn(char *msg,
                                             ...));
  /* attempt to open bfname as BAM file.  Read
header and validate it.  Return a run data
object with fields filled (including any header
text) if OK, else NULL */

BF_RUNDATA *bf_opencheckbamfile(char *bfname);
  /* attempt to open bfname as BAM file.  Read
header and validate it.  Return a run data
object with fields filled (including any header
text) if OK, else NULL */

BF_AUX_INFO *bf_appndauxinfo(BF_AUX_INFO **auxinfp,
                             char atag[2],
                             char valtype,
                             void *datp);
/* append an element atag & valtype to the end of auxinfp.
set datp as defined,
return address of new element in case useful */

int bf_cntauxinflst(BF_AUX_INFO *auxinflst);
  /* recursively count auxinflst */

void bf_delauxinfolst(BF_AUX_INFO **auxinfp);
/* iteratively delete elements from auxinfp */

void bf_disposeexpalignrec(BF_EXPNDALIGNREC *earp);
  /* recover storage of earp */

BF_EXPNDALIGNREC *bf_nxtalignrec(BF_RUNDATA *rdp,
                                 int skipaux,
                                 void err_fn(char *msg,
                                             ...));
/* read a new alignment record from rdp->dfptr, returning
the structure containing the data, use eff_fn to 
return any necessary error info, if skipaux just
lose any aux info */

void bf_expalignrec2sam(FILE *ofl,
                        BF_RUNDATA *rdp,
                        BF_EXPNDALIGNREC *earp);
/* write earp stuff out in sam type format */
