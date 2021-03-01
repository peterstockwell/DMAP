/* bam_fns.c: routines intended to provide a simple
callable set of routines to return BAM file material
in defined structures to the caller.  Intended use is
to permit reading of BAM file alignments within my
own SW without requiring bam->sam conversion of the
files first.

Peter Stockwell: Apr-2015 */
/* PAS: Correction for compression block alignment and read structures: Nov-2020 */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <fcntl.h>
#include <ctype.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <sys/param.h>
#include <zlib.h>
#include <sys/types.h>

#include "bas_fns.h"
#include "bam_fns.h"

/* endian functions draw on samtools/bam definitions */

int bf_is_bigendian()
  /* return true if machine representation is BE */
{
long one = 1;

return(!*((char *) &one));
}

uint16_t bf_endswap_uint16(uint16_t v)
{
return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}

uint32_t bf_endswap_uint32(uint32_t v)
{
v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}

void bf_err_msg_die(char *fmt,
                    ...)
/* write user error message then exit with error status */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
exit(1);
}

void bf_err_msg(char *fmt,
                ...)
/* write user error message */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
}

void bf_null_err_msg(char *fmt,
                     ...)
/* do nothing */
{
}

int bf_chk_bgzfhdr(BF_BGZF_HDR *hptr)
  /* check for valid  values in *hptr.
If OK, return the number of xlen bytes,
else 0 */
{
if ((hptr->id1 == 31) && (hptr->id2 == 139) &&
      (hptr->cm == 8) && (hptr->flg = 4))
  return(hptr->xlen);
else
  return(0);
}

int bf_chk_subfield(BF_SUBFIELD *sfptr)
  /* check for expected fields in *sfptr
& if so return the buffer size, else 0 */
{
if ((sfptr->si1 == 66) && (sfptr->si2 == 67) &&
     (sfptr->slen == 2))
  return(sfptr->bsize);
else
  return(0);
}

void bf_prt_bgzf_hdr(FILE *ofl,
                     BF_BGZF_HDR *hptr)
{
if (hptr != NULL)
  fprintf(ofl,"ID1: %hhu\tID2: %hhu\tCM: %hhu\tFLG: %hhu\tXFL: %hhu\tOS: %hhu\tXLEN: %hhu\n",
            hptr->id1,hptr->id2,hptr->cm,hptr->flg,hptr->xfl,hptr->os,hptr->xlen);
else
  fputs("<NULL>\n",ofl);
}

void bf_prt_subfield(FILE *ofl,
                     BF_SUBFIELD *sfp)
{
if (sfp != NULL)
  fprintf(ofl,"SI1: %hhu\t SI2: %hhu\tSLEN: %hhu\tBSIZE: %hu\n",
            sfp->si1,sfp->si2,sfp->slen,sfp->bsize);
else
  fputs("<NULL>\n",ofl);
}

void bf_chkbuffer(char **buf,
                  int reqlen,
                  int *blen)
/* check if *blen exceeds reqlen.  If not
free old buffer and malloc another. */
{
if (reqlen > *blen)
  {
  if (*buf != NULL)
    memfree(*buf);
  *buf = (char *) getmemzero(reqlen*sizeof(char),"buf");
  *blen = reqlen;
  }
}

int bf_chkzlibstatus(char *op,
                     int z_ret,
                     void emsgfn(char *xfmt,
                                 ...))
{
switch (z_ret)
  {
  case Z_MEM_ERROR:
    (* emsgfn)("zlib_%s: Insufficient memory\n",op);
    break;
  case Z_BUF_ERROR:
    (* emsgfn)("zlib_%s: Output buffer too small\n",op);
    break;
  case Z_DATA_ERROR:
    (* emsgfn)("zlib_%s: Corrupt/incomplete input data\n",op);
    break;
  case Z_STREAM_END:     /* OK, returned if stream end reached on success */
    z_ret = Z_OK;
    break;
  case Z_OK:
    break;
  default:
    (* emsgfn)("zlib_%s: Unknown err (%d)\n",op,z_ret);
    break;
   }
return(z_ret);
}

int bf_chkbamfmt(char *ubuf)
/* check if ubuf contains the expected bam header info.
return the header text length, -1 if corrupt.  non-NULL
htextptr will be set to start of SAM header text */
{
BF_BAMHDR_1 bf1;

if ((strncmp("BAM",ubuf,3) == 0) && (*(ubuf+3) == '\001'))
  {
  memcpy(&bf1,ubuf,sizeof(BF_BAMHDR_1));
  return((int) bf1.l_text);
  }
else
  return(-1);
}

char *bf_getrefinfo(char *buf,
                    int buflen,
                    int *refnamlen,
                    char **refnam,
                    int *refsqlen)
/* assume that buff points to the start
of a reference name block. checking that we
don't exceed the buffer size buflen, return
relevant information.  Return the final
position of the buffer after scanning the
data, NULL if can't do it */
{
int32_t rnlen;
int32_t rsqlen;
char *bp;

bp = buf;
if (buflen > sizeof(int32_t))
  {
  memcpy(&rnlen,bp,sizeof(int32_t));
  if (refnamlen != NULL)
    *refnamlen = (int) rnlen;
  bp += sizeof(int32_t);
  if ((bp + rnlen + sizeof(int32_t)) <= (buf + buflen))
    {
    if (refnam != NULL)
      *refnam = bp;
    bp += rnlen;
    memcpy(&rsqlen,bp,sizeof(int32_t));
    if (refsqlen != NULL)
      *refsqlen = (int) rsqlen;
    return(bp + sizeof(int32_t));
    }
  else
    return(NULL);
  }
else
  return(NULL);
}

void bf_fputanyc(char nc,
                 FILE *ofl)
/* try to meet C conventions for
printing any char */
{
if ((nc >= ' ') && (nc <= '~'))
  fputc(nc,ofl);
else
  switch (nc)
    {
    case '\t':
      fputs("\\t",ofl);
      break;
    case '\n':
      fputs("\\n",ofl);
      break;
    default:
      fprintf(ofl,"\\%.3o",nc);
      break;
    }
}

void bf_strout(FILE *ofl,
               char *str,
               int slen)
/* write slen chars of str out to ofl,
mapping chars to octal, etc if necessary */
{
char *sp;

sp = str;
while (slen-- > 0)
  {
  bf_fputanyc(*sp,ofl);
  sp++;
  }
}  

void bf_fputivals(uint32_t *intstream,
                  int n,
                  FILE *ofl)
{
int c;

c = 0;
while (c < n)
 {
 fprintf(ofl,"%s%d",(c==0?"":","),*(intstream+c));
 c++;
 }
}

int bf_getnewblock(BF_DATAFIELD *dfp)
  /* Assume that the file is ready to read the 
next bgzf block.  get and uncompress a new data block,
set parameters and return the uncompressed size */
{
BF_BGZF_HDR bgzfhdr;
BF_SUBFIELD subfld;
int rlen;
int sflen;
int buflen;
int r2len;

if (((rlen = read(dfp->srcfdes,&bgzfhdr,sizeof(BF_BGZF_HDR))) > 0) &&
       ((sflen = bf_chk_bgzfhdr(&bgzfhdr)) > 0))
  {
  if (((rlen = read(dfp->srcfdes,&subfld,(size_t) sflen)) > 0) &&
          ((buflen = bf_chk_subfield(&subfld)) > 0))
    {
    dfp->cdlen = buflen - bgzfhdr.xlen - 19;
    if (((rlen = read(dfp->srcfdes,dfp->cdata,dfp->cdlen)) > 0) &&
         ((r2len = read(dfp->srcfdes,&dfp->crc32,2*sizeof(uint32_t))) > 0))
      {
      dfp->zstr->next_in = dfp->cdata;
      dfp->zstr->avail_in = dfp->cdlen;
      dfp->zstr->next_out = dfp->udata;
      dfp->zstr->avail_out = dfp->udbuflen;
      if (dfp->isize > dfp->udbuflen)
        bf_err_msg_die("Buffer too small for decompress at block %d\n",dfp->blockno);
      else
        if ((bf_chkzlibstatus("init",inflateInit2(dfp->zstr,-15),bf_err_msg_die) == Z_OK) &&
               (bf_chkzlibstatus("inflate",inflate(dfp->zstr,Z_FINISH),bf_err_msg_die) == Z_OK) &&
               (bf_chkzlibstatus("end",inflateEnd(dfp->zstr),bf_err_msg_die) == Z_OK))
          {
          dfp->udatp = 0;
          dfp->blockno++;
/*           dfp->blockremain = 0; */  /* read data does not necessarily align with compression blocks */
          return(1);
          }
        else
          return(0);
      }
    else
      bf_err_msg_die("Unexpected EOF at block %d\n",dfp->blockno);
    }
  else
    {
    if (rlen <= 0)
      bf_err_msg_die("Unexpected EOF at block %d\n",dfp->blockno);
    else
      bf_err_msg_die("Invalid subfield header at block %d\n",dfp->blockno);
    }
  }
else
  if ((rlen > 0) && (sflen <= 0))
    bf_err_msg_die("Invalid BZGF header at block %d\n",dfp->blockno+1);
  else
    return(0);
}
  
void bf_initdatafield(BF_DATAFIELD *dfp,
                      FILE *srcfile)
/* make initial settings in dfp */
{
dfp->cdbuflen = dfp->udbuflen = 0x10000;
dfp->cdata = (char *) getmemzero(dfp->cdbuflen,"compr buf");
dfp->udata = (char *) getmemzero(dfp->udbuflen,"uncompr buf");
dfp->zstr = (z_stream *) getmemzero(sizeof(z_stream),"z_stream");
dfp->zstr->zalloc = NULL;
dfp->zstr->zfree = NULL;
dfp->udatp = dfp->udbuflen;  /* force new buffer read */
dfp->cdatp = dfp->cdbuflen;
dfp->blockno = 0;
dfp->srcfdes = fileno(srcfile);
dfp->srcfl = srcfile;
}

uint8_t bf_nxtuncompbyte(BF_DATAFIELD *dfp,
                         int *isok)
  /* return the next byte from the uncompressed
stream.  Redo block if necessary */
{
uint8_t nc;

if (dfp->udatp < dfp->isize)
  {
  dfp->isok = 1;
  dfp->blockremain--;
  return(*(dfp->udata + dfp->udatp++));
  }
else
  if (bf_getnewblock(dfp) > 0)
    return(bf_nxtuncompbyte(dfp,isok));
  else
    {
    dfp->isok = 0;
    return('\0');
    }
}

int bf_fill_structure(void *sptr,
                      size_t nbytes,
                      BF_DATAFIELD *dfp)
/* put nbytes of uncompressed data into sptr
from dfp.  Return number of bytes transferred */
{
int isok;
int bcnt;
uint8_t *dp;

bcnt = 0;
dp = (uint8_t *) sptr;
dfp->isok = 1;
while ((bcnt < nbytes) && dfp->isok)
  {
  *dp = bf_nxtuncompbyte(dfp,&dfp->isok);
  bcnt++;
  dp++;
  }
return(bcnt);
}

int bf_fill_strct_blkrem(void *sptr,
                         size_t nbytes,
                         BF_DATAFIELD *dfp)
/* put nbytes of uncompressed data into sptr
from dfp while dfp->blockremain > 0.
Return number of bytes transferred */
{
int isok;
int bcnt;
uint8_t *dp;

bcnt = 0;
dp = (uint8_t *) sptr;
dfp->isok = 1;
while ((bcnt < nbytes) && isok && (dfp->blockremain > 0))
  {
  *dp = bf_nxtuncompbyte(dfp,&dfp->isok);
  bcnt++;
  dp++;
  }
return(bcnt);
}

char bf_int2sqchar(int ival)
  /* ival is a 4 bit compressed code
which should map to a residue.  Return
this */
{
char bascodes[] = "=ACMGRSVTWYHKDBN";

if ((ival >= 0) && (ival <= 15))
  return(bascodes[ival]);
else
  return(bascodes[15]);
}

char bf_int2cigarchar(int32_t ival)
  /* similar to above for CIGAR codes */
{
char cigarcodes[] = "MIDNSHP=X";

if ((ival >=0) && (ival <= 8))
  return(cigarcodes[ival]);
else
  return('*');
}

void bf_prtbfarec(BF_ALIGNREC *bfarecp,
                  FILE *ofl)
/* tell ofl about contents of bfarecp */
{
fprintf(ofl,"block_size: %d refID: %d pos-1: %d bin: %o MAPQ: %d l_read_name: %d\n",
          bfarecp->block_size,bfarecp->refID,bfarecp->pos,bfarecp->bin,
          bfarecp->MAPQ,bfarecp->l_read_name);
fprintf(ofl," FLAG: %o n_cigar_op: %d l_seq: %d next_refID: %d next_pos: %d tlen: %d\n",
          bfarecp->FLAG,bfarecp->n_cigar_op,bfarecp->l_seq,bfarecp->next_refID,
          bfarecp->next_pos,bfarecp->tlen);
}

size_t bf_auxchr2byteno(char bc)
  /* number of bytes for char bc */
{
switch (bc)
  {
  case 'c':
  case 'C':
    return(sizeof(int8_t));
    break;
  case 's':
  case 'S':
    return(sizeof(int16_t));
    break;
  case 'i':
  case 'I':
    return(sizeof(int32_t));
    break;
  case 'f':
    return(sizeof(float));
    break;
  default:
    return(0);
    break;
  }
}

int bf_getbamauxint4hdr(BF_AUX_INFO *dhp,
                        BF_DATAFIELD *dfp,
                        int *isok)
/* depending on type of dhp->val_type, return
an integer value */
{
uint8_t ival8;
uint16_t ival16;
uint32_t ival32;

dfp->isok = 1;
switch (dhp->val_type)
  {
  case 'c':
  case 'C':
    if (bf_fill_strct_blkrem(&ival8,bf_auxchr2byteno(dhp->val_type),dfp) ==
          bf_auxchr2byteno(dhp->val_type))
      return((int)ival8);
    break;
  case 's':
  case 'S':
    if (bf_fill_strct_blkrem(&ival16,bf_auxchr2byteno(dhp->val_type),dfp) ==
          bf_auxchr2byteno(dhp->val_type))
      return((int)ival16);
    break;
  case 'i':
  case 'I':
    if (bf_fill_strct_blkrem(&ival32,bf_auxchr2byteno(dhp->val_type),dfp) ==
          bf_auxchr2byteno(dhp->val_type))
      return((int)ival32);
    break;
  default:
    break;
  }
dfp->isok = 0;
return(0);
}

float bf_getbamauxflt(BF_DATAFIELD *dfp,
                      int *isok)
/* return a float: assumed that the
val_type char is f */
{
float fval;

dfp->isok = 1;
if (bf_fill_strct_blkrem(&fval,bf_auxchr2byteno('f'),dfp) == bf_auxchr2byteno('f'))
  return(fval);
else
  {
  dfp->isok = 0;
  return(0);
  }
}

BF_DATAFIELD *bf_openinbamflmsg(char *bfname,
                                void err_fn(char *msg,
                                              ...))
  /* attempt to open bamfile bfname for input.
return the BF_DATAFIELD initialised if successful
otherwise NULL. err_fn used to communicate error
issues */
{
BF_DATAFIELD *bfp;
FILE *srcfl;

if ((srcfl = fopen(bfname,"r")) != NULL)
  {
  bfp = (BF_DATAFIELD *) getmemzero(sizeof(BF_DATAFIELD),"bamdatafield");
  bf_initdatafield(bfp,srcfl);
  return(bfp);
  }
else
  {
  (* err_fn)("Can't open BAM input file '%s'\n",bfname);
  return(NULL);
  }
}

BF_DATAFIELD *bf_openinbamfile(char *bfname)
  /* attempt to open bamfile bfname for input.
return the BF_DATAFIELD initialised if successful
otherwise NULL */
{
return(bf_openinbamflmsg(bfname,bf_null_err_msg));
}

void bf_dispose_datafield(BF_DATAFIELD *dfp)
  /* lose all storage associated with dfp */
{
if (dfp != NULL)
  {
  memfree(dfp->cdata);
  memfree(dfp->udata);
  memfree(dfp->zstr);
/*   close(dfp->srcfdes); */
  fclose(dfp->srcfl);
  dfp->srcfl = NULL;
  memfree(dfp);
  }
}

void bf_initrundata(BF_RUNDATA *rdp)
  /* init *rdp to default settings */
{
rdp->dfptr = NULL;
rdp->bhdrtxtlen = 0;
rdp->hdrtext = NULL;
rdp->n_ref = 0;
rdp->refinfo = NULL;
rdp->bigendian = bf_is_bigendian();
}

void bf_disposerundata(BF_RUNDATA *rdp)
  /* do the obvious */
{
int i;

if (rdp != NULL)
  {
  bf_dispose_datafield(rdp->dfptr);
  if (rdp->hdrtext != NULL)
    memfree(rdp->hdrtext);
  for (i = 0; i < rdp->n_ref; i++)
    if ((rdp->refinfo + i)->refname != NULL)
      memfree((rdp->refinfo+i)->refname);
  if (rdp->refinfo != NULL)
    memfree(rdp->refinfo);
  memfree(rdp);
  }
}

BF_RUNDATA *bf_opnchkbamflmsg(char *bfname,
                              void err_fn(char *msg,
                                            ...))
  /* attempt to open bfname as BAM file.  Read
header and validate it.  Return a run data
object with fields filled (including any header
text) if OK, else NULL */
{
BF_DATAFIELD *dfp;
int bhlen;
BF_RUNDATA *rdp;
BF_BAMHDR_1 bhdr1;
int i;
int fillen;

if ((dfp = bf_openinbamflmsg(bfname,err_fn)) != NULL)
  {
  if (((fillen = bf_fill_structure(&bhdr1,sizeof(BF_BAMHDR_1),dfp)) == sizeof(BF_BAMHDR_1))
        && ((bhlen = bf_chkbamfmt((char *)&bhdr1)) >= 0))
    {
    rdp = (BF_RUNDATA *)getmemzero(sizeof(BF_RUNDATA),"rundata");
    bf_initrundata(rdp);
    rdp->dfptr = dfp;
    if (rdp->bigendian)
      rdp->bhdrtxtlen = bf_endswap_uint32(bhlen);
    else
      rdp->bhdrtxtlen = bhlen;
    rdp->hdrtext = (char *)getmemzero(rdp->bhdrtxtlen,"headertext");
    if ((bf_fill_structure(rdp->hdrtext,rdp->bhdrtxtlen,rdp->dfptr) == rdp->bhdrtxtlen)
           && (bf_fill_structure(&rdp->n_ref,sizeof(int32_t),rdp->dfptr) == sizeof(int32_t)))
      {
      if (rdp->bigendian)
        rdp->n_ref = bf_endswap_uint32(rdp->n_ref);
      rdp->refinfo = (BF_REFINFO *) getmemzero(sizeof(BF_REFINFO)*rdp->n_ref,"refseqdata");
      bzero(rdp->refinfo,sizeof(BF_REFINFO)*rdp->n_ref);
      i = 0;
      while (i < rdp->n_ref)
        if (bf_fill_structure(&(rdp->refinfo + i)->l_name,sizeof(int32_t),rdp->dfptr) !=
              sizeof(int32_t))
          {
          bf_disposerundata(rdp);
          (* err_fn)("Unexpected EOF at ref seq info for #%d in '%s'\n",i+1,bfname);
          return(NULL);
          }
        else
          {
          if (rdp->bigendian)
            (rdp->refinfo+i)->l_name = bf_endswap_uint32((rdp->refinfo+i)->l_name);
          (rdp->refinfo + i)->refname = (char *)getmemzero((rdp->refinfo+i)->l_name,"refname");
          if (bf_fill_structure((rdp->refinfo+i)->refname,(rdp->refinfo+i)->l_name,rdp->dfptr)
                != (rdp->refinfo+i)->l_name)
            {
            bf_disposerundata(rdp);
            (* err_fn)("Unexpected EOF at ref seq info for #%d in '%s'\n",i+1,bfname);
            return(NULL);
            }
          else
            {
            if (bf_fill_structure(&(rdp->refinfo+i)->l_ref,sizeof(int32_t),rdp->dfptr)
                 != sizeof(int32_t))
              {
              bf_disposerundata(rdp);
              (* err_fn)("Unexpected EOF at ref seq info for #%d in '%s'\n",i+1,bfname);
              return(NULL);
              }
            else
              {
              if (rdp->bigendian)
                (rdp->refinfo+i)->l_ref = bf_endswap_uint32((rdp->refinfo+i)->l_ref);
              i++;
              }
            }
          }
      return(rdp);
      }
    else
      {
      bf_disposerundata(rdp);
      (* err_fn)("Unexpected EOF in '%s'\n",bfname);
      return(NULL);
      }
    }
  else
    {
    bf_dispose_datafield(dfp);
    (* err_fn)("Invalid BAM header in '%s'\n",bfname);
    return(NULL);
    }
  }
else
  {
  (* err_fn)("Can't open BAM input file '%s'\n",bfname);
  return(NULL);
  }
}

BF_RUNDATA *bf_opencheckbamfile(char *bfname)
  /* attempt to open bfname as BAM file.  Read
header and validate it.  Return a run data
object with fields filled (including any header
text) if OK, else NULL */
{
return(bf_opnchkbamflmsg(bfname,bf_null_err_msg));
}

BF_AUX_INFO *bf_appndauxinfo(BF_AUX_INFO **auxinfp,
                             char atag[2],
                             char valtype,
                             void *datp)
/* append an element atag & valtype to the end of auxinfp.
set datp as defined,
return address of new element in case useful */
{
BF_AUX_INFO *prev, *end_ptr;

if (auxinfp != NULL)
  {
  prev = end_ptr = *auxinfp;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtaux;
    }
  end_ptr = (BF_AUX_INFO *) getmemzero(sizeof(BF_AUX_INFO),"auxinfo elt");
  end_ptr->nxtaux = NULL;
  memcpy(&end_ptr->tag[0],&atag[0],2);
  end_ptr->val_type = valtype;
  end_ptr->auxdat = datp;
  end_ptr->arrtype = '\0';    /* set separately */
  end_ptr->arrbtlen = 0;        /* ditto */
  if (*auxinfp == NULL)
    {    
    *auxinfp = end_ptr;
    end_ptr->prvaux = NULL;
    }
  else
    {
    prev->nxtaux = end_ptr;
    end_ptr->prvaux = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

int bf_cntauxinflst(BF_AUX_INFO *auxinflst)
  /* recursively count auxinflst */
{
if (auxinflst == NULL)
  return(0);
else
  return(bf_cntauxinflst(auxinflst->nxtaux) + 1);
}

void bf_delauxinfolst(BF_AUX_INFO **auxinfp)
/* iteratively delete elements from auxinfp */
{
BF_AUX_INFO *fp;

while (*auxinfp != NULL)
  {
  fp = *auxinfp;
  *auxinfp = fp->nxtaux;
  if (fp->auxdat != NULL)
    memfree(fp->auxdat);
  memfree(fp);
  }
}

void bf_disposeexpalignrec(BF_EXPNDALIGNREC *earp)
  /* recover storage of earp */
{
if (earp != NULL)
  {
  if (earp->read_name != NULL)
    memfree(earp->read_name);
  if (earp->cigar != NULL)
    memfree(earp->cigar);
  if (earp->seq != NULL)
    memfree(earp->seq);
  if (earp->qual != NULL)
    memfree(earp->qual);
  if (earp->arec != NULL)
    memfree(earp->arec);
  bf_delauxinfolst(&earp->auxinfo);
  memfree(earp);
  }
}

BF_EXPNDALIGNREC *bf_nxtalignrec(BF_RUNDATA *rdp,
                                 int skipaux,
                                 void err_fn(char *msg,
                                             ...))
/* read a new alignment record from dfp, returning
the structure containing the data, use eff_fn to 
return any necessary error info, if skipaux just
lose any aux info */
{
BF_ALIGNREC *arp;
int bfareclen;
BF_EXPNDALIGNREC *earp;
char *rp;
int i;
int isok;
uint8_t nc;
BF_AUX_INFO ainfo;
int auxarrlen;
char auxarrtype;
char auxstr[0x10000];
BF_AUX_INFO *newap;
int abytelen;

arp = (BF_ALIGNREC *) getmemzero(sizeof(BF_ALIGNREC),"alignrec");
bfareclen = sizeof(BF_ALIGNREC);
if (bf_fill_structure(arp,bfareclen,rdp->dfptr) == bfareclen)
  {
  if (rdp->bigendian)
    arp->block_size = bf_endswap_uint32(arp->block_size);
  rdp->dfptr->blockremain = arp->block_size - bfareclen + sizeof(int32_t);
  earp = (BF_EXPNDALIGNREC *) getmemzero(sizeof(BF_EXPNDALIGNREC),"expalignrec");
  earp->read_name = earp->seq = earp->qual = NULL;
  earp->arec = arp;
  earp->read_name = (char *) getmemzero((size_t) arp->l_read_name,"readname");
  if (rdp->bigendian)
    {
    arp->n_cigar_op = bf_endswap_uint16(arp->n_cigar_op);
    arp->refID = bf_endswap_uint32(arp->refID);
    arp->l_seq = bf_endswap_uint32(arp->l_seq);
    arp->pos = bf_endswap_uint32(arp->pos);
    arp->next_refID = bf_endswap_uint32(arp->next_refID);
    arp->next_pos = bf_endswap_uint32(arp->next_pos);
    arp->tlen = bf_endswap_uint32(arp->tlen);
    }
  earp->cigar = (BF_CIGAR_OP_TYPE *) getmemzero(sizeof(BF_CIGAR_OP_TYPE)*arp->n_cigar_op,"cigar");
  earp->seq = (char *) getmemzero((size_t)arp->l_seq,"seq");
  earp->qual = (char *) getmemzero((size_t)arp->l_seq,"qual");
  earp->auxinfo = NULL;
  i = arp->l_read_name;
  rdp->dfptr->isok = 1;
  rp = earp->read_name;
  while ((i > 0) && rdp->dfptr->isok)
    {
    *rp = bf_nxtuncompbyte(rdp->dfptr,&rdp->dfptr->isok);
    rp++;
    i--;
    }
  i = 0;
  while ((i < arp->n_cigar_op) && rdp->dfptr->isok)
    {
    rdp->dfptr->isok = bf_fill_structure((earp->cigar+i),sizeof(BF_CIGAR_OP_TYPE),rdp->dfptr)
              == sizeof(BF_CIGAR_OP_TYPE);
    i++;
    }
  i = (arp->l_seq + 1)/2;
  rp = earp->seq;
  while ((i > 0) && rdp->dfptr->isok)
    {
    nc = bf_nxtuncompbyte(rdp->dfptr,&rdp->dfptr->isok);
    *rp = bf_int2sqchar(nc >> 4);
    rp++;
    *rp = bf_int2sqchar(nc & 0x0f);
    rp++;
    i--;
    }
  i = arp->l_seq;
  rp = earp->qual;
  while ((i > 0) && rdp->dfptr->isok)
    {
    *rp = bf_nxtuncompbyte(rdp->dfptr,&rdp->dfptr->isok) + 33;
    rp++;
    i--;
    }
  if (skipaux)
    {
    while ((rdp->dfptr->blockremain > 0) && rdp->dfptr->isok)
      nc = bf_nxtuncompbyte(rdp->dfptr,&rdp->dfptr->isok);
    }
  else
    {
    memset(auxstr,'\0',0x10000);
    while ((rdp->dfptr->blockremain > 0) && rdp->dfptr->isok)
      {
      rdp->dfptr->isok = bf_fill_strct_blkrem(&ainfo,3,rdp->dfptr) == 3;
      if (rdp->dfptr->isok)
        switch (ainfo.val_type)
          {
          case 'c':
          case 'C':
          case 's':
          case 'S':
          case 'i':
          case 'I':
            newap = bf_appndauxinfo(&earp->auxinfo,&ainfo.tag[0],ainfo.val_type,NULL);
            newap->aival = bf_getbamauxint4hdr(&ainfo,rdp->dfptr,&rdp->dfptr->isok);
/*              fprintf(stdout,"%.2s:%c:%d\n",aux_dathdr.tag,aux_dathdr.val_type,auxival); */
            break;
          case 'f':
            newap = bf_appndauxinfo(&earp->auxinfo,&ainfo.tag[0],ainfo.val_type,NULL);
            newap->afval = bf_getbamauxflt(rdp->dfptr,&rdp->dfptr->isok);
/*              fprintf(stdout,"%.2s:%c:%f\n",aux_dathdr.tag,aux_dathdr.val_type,auxfval); */
            break;
          case 'Z':
            rp = &auxstr[0];
            do
              {
              *rp = nc = bf_nxtuncompbyte(rdp->dfptr,&rdp->dfptr->isok);
              rp++;
              }
            while ((nc != '\0') && rdp->dfptr->isok);
            newap = bf_appndauxinfo(&earp->auxinfo,&ainfo.tag[0],ainfo.val_type,
                                      bas_strdup(&auxstr[0]));
/*            fprintf(stdout,"%.2s:%c:",aux_dathdr.tag,aux_dathdr.val_type); */
            break;
          case 'B':
            auxarrtype = bf_nxtuncompbyte(rdp->dfptr,&rdp->dfptr->isok);
            if ((rdp->dfptr->blockremain > 0) && rdp->dfptr->isok &&
                  (bf_fill_strct_blkrem(&auxarrlen,sizeof(int32_t),rdp->dfptr) == sizeof(int32_t)))
              {
              rdp->dfptr->isok = 1;
              if (rdp->bigendian)
                auxarrlen = bf_endswap_uint32(auxarrlen);
              abytelen = auxarrlen * bf_auxchr2byteno(auxarrtype);
              if (bf_fill_structure(&auxstr[0],abytelen,rdp->dfptr) != abytelen)
                (* err_fn)("Unexpected EOF for aux data '%.2s:%c' at block %d\n",
                             &ainfo.tag[0],ainfo.val_type,rdp->dfptr->blockno);
              else
                {
                newap = bf_appndauxinfo(&earp->auxinfo,&ainfo.tag[0],ainfo.val_type,NULL);
                newap->arrtype = auxarrtype;
                newap->arrbtlen = abytelen;
                newap->auxdat = (char *) getmemzero(abytelen,"B auxdata");
                memcpy(newap->auxdat,&auxstr[0],abytelen);
                }
/*             fprintf(stdout,"%.2s:%c:%dvalues@%dbytes\n",aux_dathdr.tag,aux_dathdr.val_type,
                      auxarrlen,bf_auxchr2byteno(auxarrtype)); */
              }
            break;
          default:
            while ((rdp->dfptr->blockremain > 0) && rdp->dfptr->isok)
              nc = bf_nxtuncompbyte(rdp->dfptr,&rdp->dfptr->isok);
            break;
          }
      }
    }
  if (rdp->dfptr->isok)
    return(earp);
  else
    {
    (* err_fn)("Unexpected EOF at bgzf block %d\n",rdp->dfptr->blockno);
    bf_disposeexpalignrec(earp);
    return(NULL);
    }
  }
else
  return(NULL);
}

void bf_expalignrec2sam(FILE *ofl,
                        BF_RUNDATA *rdp,
                        BF_EXPNDALIGNREC *earp)
/* write earp stuff out in sam type format */
{
int i;
char *sp;
char mq;
BF_AUX_INFO *aip;

if (earp != NULL)
  {
  fprintf(ofl,"%s\t%d\t%s\t%d\t%hhu\t",
            earp->read_name,earp->arec->FLAG,
            (earp->arec->refID>=0?(rdp->refinfo+earp->arec->refID)->refname:"*"),
            earp->arec->pos+1,earp->arec->MAPQ);
  if (earp->arec->n_cigar_op <= 0)
    fputc('*',ofl);
  else
    {
    i = 0;
    while (i < earp->arec->n_cigar_op)
      {
      fprintf(ofl,"%d%c",(earp->cigar+i)->c_op_len,
                bf_int2cigarchar((int) (earp->cigar+i)->cigar_op));
      i++;
      }
    }
  fprintf(ofl,"\t%s\t%d\t%d\t",
            (earp->arec->next_refID==-1?"*":(rdp->refinfo+earp->arec->next_refID)->refname),
            earp->arec->next_pos+1,earp->arec->tlen);
  i = earp->arec->l_seq;
  sp = earp->seq;
  while (i > 0)
    {
    fputc(*sp,ofl);
    sp++;
    i--;
    }
  fputc('\t',ofl);
  i = earp->arec->l_seq;
  sp = earp->qual;
  while (i > 0)
    {
    fputc(*sp,ofl);
    sp++;
    i--;
    }
  aip = earp->auxinfo;
  while (aip != NULL)
    {
    fputc('\t',ofl);
    switch (aip->val_type)
      {
      case 'c':
      case 'C':
      case 's':
      case 'S':
      case 'i':
      case 'I':
        fprintf(stdout,"%.2s:%c:%d",aip->tag,aip->val_type,aip->aival);
        break;
      case 'f':
        fprintf(stdout,"%.2s:%c:%f",aip->tag,aip->val_type,aip->afval);
        break;
      case 'Z':
        fprintf(stdout,"%.2s:%c:%s",aip->tag,aip->val_type,(char *)aip->auxdat);
        break;
      case 'B':  /* not writing this type out at present */
        fprintf(stdout,"%.2s:%c:%dvalues@%dbytes\n",aip->tag,aip->val_type,
                      aip->arrbtlen,bf_auxchr2byteno(aip->arrtype));
        break;
      default:
        break;
      }
    aip = aip->nxtaux;
    }
  fputc('\n',ofl);
  }
else
  fputs("<NULL>\n",ofl);
}
