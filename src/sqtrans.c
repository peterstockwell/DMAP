/* sqtrans.c : routines for sequence translation and genetic code 
     manipulation */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "bas_fns.h"
#include "sqtrans.h"

#ifdef MALLOCDBG
int sqt_chkvectarr(COD_VECT_STRCT **cvarr)
  /* traverse all vectors for surveilance purposes.
  return number of valid nodes - die otherwise */
{
int ap;
COD_VECT_STRCT *cvp;
int ncnt;

ncnt = 0;
for (ap = 0; ap <= MAXAAINT; ap++)
  {
  cvp = *(cvarr + ap);
  while (cvp != NULL)
    {
    if (scan4free(cvp) != NULL)
      ncnt++;
    else
      {
      fprintf(stderr,"\nCvec pointer error: %lx\n",(long) cvp);
      exit(1);
      }
    cvp = cvp->nextcvec;
    }
  }
return(ncnt);
}
#endif

char int2nabas(int iv)
/* return base corresponding to iv: ?=0 A=1 C=2 G=3 T=4 */
{
switch (iv)
  {
  case 1:
    return('A');
    break;
  case 2:
    return('C');
    break;
  case 3:
    return('G');
    break;
  case 4:
    return('T');
    break;
  case 0:
  default:
    return('?');
    break;
  }
}

SQT_GENCOD sqt_int2gbgencod(int iv)
  /* return the corresponding symbolic name for this integer */
{
switch (iv)
  {
  case 1:
    return(GC_universal);
    break;
  case 2:
    return(GC_vertmito);
    break;
  case 3:
    return(GC_yeastmito);
    break;
  case 4:
    return(GC_mycomito);
    break;
  case 5:
    return(GC_invertmito);
    break;
  case 6:
    return(GC_ciliate);
    break;
  case 7:
  case 8:
    return(GC_undefined);
    break;
  case 9:
    return(GC_echinmito);
    break;
  case 10:
    return(GC_euplotid);
    break;
  case 11:
    return(GC_bacterial);
    break;
  case 12:
    return(GC_yeastalt);
    break;
  case 13:
    return(GC_ascidmito);
    break;
  case 14:
    return(GC_flatwmito);
    break;
  case 15:
    return(GC_blepharisma);
    break;
  case 16:
    return(GC_userdefined);
    break;
  default:
    return(GC_undefined);
    break;
  }
}

char *sqt_gbgencd2strng(SQT_GENCOD gc)
  /* return a string pointer for this gc value */
{
switch (gc)
  {
  case GC_universal:
    return("Standard/Universal");
    break;
  case GC_vertmito:
    return("Vertebrate mitochondrial");
    break;
  case GC_yeastmito:
    return("Yeast mitochondrial");
    break;
  case GC_mycomito:
    return("Mold/protozoan/coelent/mycoplasma mitochondrial");
    break;
  case GC_invertmito:
    return("Invertebrate mitochondrial");
    break;
  case GC_ciliate:
    return("Ciliate/Dasycladacean Hexamita nuclear");
    break;
  case GC_echinmito:
    return("Echinoderm mitochondrial");
    break;
  case GC_euplotid:
    return("Euplotid nuclear");
    break;
  case GC_bacterial:
    return("Prokaryotic");
    break;
  case GC_yeastalt:
    return("Alternative Yeast nuclear");
    break;
  case GC_ascidmito:
    return("Ascidian mitochondrial");
    break;
  case GC_flatwmito:
    return("Flatworm mitochrondrial");
    break;
  case GC_blepharisma:
    return("Blepharisma nuclear");
    break;
  case GC_userdefined:
    return("User defined");
    break;
  case GC_undefined:
  default:
    return("Undefined");
    break;
  }
}

void sqt_saygbgencodtbls(FILE *fl,
                         int from,
                         int to)
/* describe valid codes in range from..to at fl */
{
int gp;
SQT_GENCOD gc;
char *strng;

for (gp = from; gp <= to; gp++)
  if ((gc = sqt_int2gbgencod(gp)) != GC_undefined)
    if ((strng = sqt_gbgencd2strng(gc)) != NULL)
      fprintf(fl,"%2d: %s code\n",gp,strng);
    else
      fprintf(fl,"%2d: Invalid code",gp);
}    

int sqt_nabas2int(char abase)
/* return integer corresponding to abase: 0=? 1=A 2=C 3=G 4=T/U */
{
switch (tolower(abase))
  {
  case 'a':
    return(1);
    break;
  case 'c':
    return(2);
    break;
  case 'g':
    return(3);
    break;
  case 't':
  case 'u':
    return(4);
    break;
  default:
    return(0);
    break;
  }
}

char sqt_int2nabas(int iv)
/* return base corresponding to iv: ?=0 A=1 C=2 G=3 T=4 */
{
switch (iv)
  {
  case 1:
    return('A');
    break;
  case 2:
    return('C');
    break;
  case 3:
    return('G');
    break;
  case 4:
    return('T');
    break;
  case 0:
  default:
    return('?');
    break;
  }
}

int remapbasptr(int bp)
/* remap base pointers to return bases in order TCAG, as per conventional
  genetic code display */
{
switch (bp)
  {
  case 1:
    return(4);
    break;
  case 2:
    return(2);
    break;
  case 3:
    return(1);
    break;
  case 4:
    return(3);
    break;
  case 0:
  default:
    return(0);
    break;
  }
}

int aares2int(char aares)
/* return an appropriate integer pointer for this residue */
{
switch (tolower(aares))
  {
  case 'a':
    return(1);
    break;
  case 'r':
    return(2);
    break;
  case 'n':
    return(3);
    break;
  case 'd':
    return(4);
    break;
  case 'c':
    return(5);
    break;
  case 'q':
    return(6);
    break;
  case 'e':
    return(7);
    break;
  case 'g':
    return(8);
    break;
  case 'h':
    return(9);
    break;
  case 'i':
    return(10);
    break;
  case 'l':
    return(11);
    break;
  case 'k':
    return(12);
    break;
  case 'm':
    return(13);
    break;
  case 'f':
    return(14);
    break;
  case 'p':
    return(15);
    break;
  case 's':
    return(16);
    break;
  case 't':
    return(17);
    break;
  case 'w':
    return(18);
    break;
  case 'y':
    return(19);
    break;
  case 'v':
    return(20);
    break;
  case SQT_TERMCHR:
    return(21);
    break;
  default:
    return(0);
    break;
  }
}

char *aachr2str3(char aares)
/* return an appropriate 3 letter string for this residue */
{
switch (tolower(aares))
  {
  case 'a':
    return("Ala");
    break;
  case 'b':
    return("Asx");
    break;
  case 'r':
    return("Arg");
    break;
  case 'n':
    return("Asn");
    break;
  case 'd':
    return("Asp");
    break;
  case 'c':
    return("Cys");
    break;
  case 'q':
    return("Gln");
    break;
  case 'e':
    return("Glu");
    break;
  case 'g':
    return("Gly");
    break;
  case 'h':
    return("His");
    break;
  case 'i':
    return("Ile");
    break;
  case 'l':
    return("Leu");
    break;
  case 'k':
    return("Lys");
    break;
  case 'm':
    return("Met");
    break;
  case 'f':
    return("Phe");
    break;
  case 'p':
    return("Pro");
    break;
  case 's':
    return("Ser");
    break;
  case 't':
    return("Thr");
    break;
  case 'w':
    return("Trp");
    break;
  case 'y':
    return("Tyr");
    break;
  case 'v':
    return("Val");
    break;
  case SQT_TERMCHR:
    return("Ter");
    break;
  case 'z':
    return("Glx");
    break;
  default:
    return("???");
    break;
  }
}

int sqt_res2aaint(char ares)
  /* return an integer for this residue - same as aares2int, except that
     unknown = 0 and output is 0..MAXAAINT */
{
int rv;

if (((rv = aares2int(ares)) <= 0) || (rv >= MAXAAINT))
  rv = 0;
return(rv);
}

char int2aares(int aaint)
/* return character for this integer pointer */
{
switch (aaint)
  {
  case 1:
    return('A');
    break;
  case 2:
    return('R');
    break;
  case 3:
    return('N');
    break;
  case 4:
    return('D');
    break;
  case 5:
    return('C');
    break;
  case 6:
    return('Q');
    break;
  case 7:
    return('E');
    break;
  case 8:
    return('G');
    break;
  case 9:
    return('H');
    break;
  case 10:
    return('I');
    break;
  case 11:
    return('L');
    break;
  case 12:
    return('K');
    break;
  case 13:
    return('M');
    break;
  case 14:
    return('F');
    break;
  case 15:
    return('P');
    break;
  case 16:
    return('S');
    break;
  case 17:
    return('T');
    break;
  case 18:
    return('W');
    break;
  case 19:
    return('Y');
    break;
  case 20:
    return('V');
    break;
  case 21:
    return(SQT_TERMCHR);
    break;
  default:
    return('?');
    break;
  }
}

int set_chr_datum(TRN_MATELT *tdatum,
                  RESIDUEDATA *newdat)
/* fill this tdatum element with char data type */
{
tdatum->dtype = SQTDT_char;
tdatum->resdata.trres.reschr = newdat->trres.reschr;
return(1);
}

int set_ini_datum(TRN_MATELT *tdatum,
                  RESIDUEDATA *newdat)
/* fill this tdatum element with initiation status type */
{
tdatum->dtype = SQTDT_char;
tdatum->resdata.trres.trinit = newdat->trres.trinit;
return(1);
}

int set_int_datum(TRN_MATELT *tdatum,
                  RESIDUEDATA *newdat)
/* fill this tdatum element with integer data type */
{
tdatum->dtype = SQTDT_int;
tdatum->resdata.iresidue = newdat->iresidue;
return(1);
}

void inc_int_datum(TRN_MATELT *tdatum,
                   RESIDUEDATA *dummy)
/* increment this tdatum element */
{
tdatum->resdata.iresidue += 1;
}

int set_flt_datum(TRN_MATELT *tdatum,
                  RESIDUEDATA *newdat)
/* fill this tdatum element with float data type */
{
tdatum->dtype = SQTDT_float;
tdatum->resdata.xresidue = newdat->xresidue;
return(1);
}

TRN_MATELT *sqt_vect2eltaddr(TRANS_MATRX tm,
                             CODN_VECTR cv)
/* return the address of the array element corresponding to cv */
{
return(&tm[cv[0]][cv[1]][cv[2]]);
}

int sqt_cdn2codvect(char *cdn,
                    CODN_VECTR cv)
/* convert cdn to cv. return 1 if all bases OK, else 0 */
{
int cp;
int cc;

cc = 0;
for (cp = 0; cp < CODONLENGTH; cp++)
  if ((cv[cp] = sqt_nabas2int(*(cdn+cp))) > 0)
    cc++;
return(cc == CODONLENGTH);
}

int sqt_cdnno2cdn(int cno,
                  char *cdnbuf)
/* write to cdnbuf (assumed to be
CODONLENGTH+1 chars or longer) he codon corresponding to
enumerated codon cno (0..63). return 1 if OK, else 0 */
{
int cp;
int acp;

*(cdnbuf + CODONLENGTH) = '\0';
if ((cno >= 0) && (cno < MAX_CODON_NO))
  {
  cp = CODONLENGTH - 1;
  while (cp >= 0)
    {
    switch (cp)
      {
      case 2:
        acp = 2;
        break;
      case 1:
        acp = 0;
        break;
      case 0:
        acp = 1;
        break;
      default:
        acp = cp;
        break;
      }
    *(cdnbuf + acp) = int2nabas(remapbasptr((cno % MAXBASINT) + 1));
    cno = (int) cno / MAXBASINT;
    cp--;
    }
  return(1);
  }
else
  return(0);
}

TRN_MATELT *sqt_cdnno2codvect(TRANS_MATRX tmat,
                              int cno,
                              char *ucdn)
/* if cno is valid, return the vector pointer
for it for translation table tmat.  NULL if
can't.  If ucdn is non-NULL, write the codon
to it (assumes enough room) */
{
int ccp;
char cdn[CODONLENGTH + 1];
CODN_VECTR cv;

if (sqt_cdnno2cdn(cno,&cdn[0]))
  {
  for (ccp = 0; ccp < CODONLENGTH; ccp++)
    cv[ccp] = sqt_nabas2int(cdn[ccp]);
  if (ucdn != NULL)
    memcpy(ucdn,&cdn[0],(size_t) (CODONLENGTH + 1));
  return(sqt_vect2eltaddr(tmat,cv));
  }
else
  return(NULL);
}

int sqt_codvect2cdn(CODN_VECTR cv,
                    char *cdn)
/* convert cv to corresponding codon in cdn (presumed to be big enough)
  return 1 if all of codon is valid, 0 otherwise */
{
int cp;
char *sp;
char nb;
int cc;

cc = 0;
sp = cdn;
for (cp = 0; cp < CODONLENGTH; cp++)
  if ((nb = sqt_int2nabas(cv[cp])) != '?')
    {
    bas_appchr(cdn,&sp,nb,CODONLENGTH);
    cc++;
    }
bas_appchr(cdn,&sp,'\0',CODONLENGTH);
return(cc == CODONLENGTH);
}

void fill_chr_trnarr(TRANS_MATRX trnsar,
                     char fillch)
/* fill trnsarr with fill character */
{
int p1, p2, p3;
RESIDUEDATA rdat;

rdat.trres.reschr = fillch;
for (p1 = 0; p1 <= MAXBASINT; p1++)
  for (p2 = 0; p2 <= MAXBASINT; p2++)
    for (p3 = 0; p3 <= MAXBASINT; p3++)
      (void) set_chr_datum(&trnsar[p1][p2][p3],&rdat);
}

void fill_ini_trnarr(TRANS_MATRX trnsar,
                     int fillini)
/* fill trnsarr with fill character */
{
int p1, p2, p3;
RESIDUEDATA rdat;

rdat.trres.trinit = fillini;
for (p1 = 0; p1 <= MAXBASINT; p1++)
  for (p2 = 0; p2 <= MAXBASINT; p2++)
    for (p3 = 0; p3 <= MAXBASINT; p3++)
      (void) set_ini_datum(&trnsar[p1][p2][p3],&rdat);
}

void fill_int_trnarr(TRANS_MATRX trnsar,
                     int fillval)
/* fill trnsarr with fill value */
{
int p1, p2, p3;
RESIDUEDATA rdat;

rdat.iresidue = fillval;
for (p1 = 0; p1 <= MAXBASINT; p1++)
  for (p2 = 0; p2 <= MAXBASINT; p2++)
    for (p3 = 0; p3 <= MAXBASINT; p3++)
      (void) set_int_datum(&trnsar[p1][p2][p3],&rdat);
}

void fill_flt_trnarr(TRANS_MATRX trnsar,
                     float fillval)
/* fill float trnsarr with fill value */
{
int p1, p2, p3;
RESIDUEDATA rdat;

rdat.xresidue = fillval;
for (p1 = 0; p1 <= MAXBASINT; p1++)
  for (p2 = 0; p2 <= MAXBASINT; p2++)
    for (p3 = 0; p3 <= MAXBASINT; p3++)
      (void) set_flt_datum(&trnsar[p1][p2][p3],&rdat);
}

int set_cdnpos3(TRANS_MATRX trnsar,
                int b1,
                int b2,
                char *cdn,
                RESIDUEDATA *dat_ptr,
                int (*fill_proc)(TRN_MATELT *td,
                                 RESIDUEDATA *rd))
/* consider 3rd codon position, set rest of codon for this, b1 and b2.  
  Wild card * causes all nonzero values to be initialized */
{
int b3pt;
int fcnt;

fcnt = 0;
if (*cdn == '*')
  for (b3pt = 1; b3pt <= 4; b3pt++)
    fcnt += (*fill_proc)(&trnsar[b1][b2][b3pt],dat_ptr);
else
  fcnt += (*fill_proc)(&trnsar[b1][b2][sqt_nabas2int(*cdn)],dat_ptr);
return(fcnt);
}

int set_cdnpos2(TRANS_MATRX trnsar,
                int b1,
                char *cdn,
                RESIDUEDATA *dat_ptr,
                int (*fill_proc)(TRN_MATELT *td,
                                 RESIDUEDATA *rd))
/* consider second codon position, set rest of codon for this and b1.  
  Wild card * causes all nonzero values to be initialized */
{
int b2pt;
char bas;
int fcnt;

fcnt = 0;
if ((bas = *cdn++) == '*')
  for (b2pt = 1; b2pt <= 4; b2pt++)
    fcnt += set_cdnpos3(trnsar,b1,b2pt,cdn,dat_ptr,fill_proc);
else
  fcnt += set_cdnpos3(trnsar,b1,sqt_nabas2int(bas),cdn,dat_ptr,fill_proc);
return(fcnt);
}

int set_trnarr(TRANS_MATRX trnsar,
               char *cdn,
               RESIDUEDATA *dat_ptr,
               int (*fill_proc)(TRN_MATELT *td,
                                RESIDUEDATA *rd))
/* consider first codon position, set rest of codon for this.  Wild card *
  causes all nonzero values to be initialized */
{
int bpt;
char bas;
int fcnt;

fcnt = 0;
if ((bas = *cdn++) == '*')
  for (bpt = 1; bpt <= 4; bpt++)
    fcnt += set_cdnpos2(trnsar,bpt,cdn,dat_ptr,fill_proc);
else
  fcnt += set_cdnpos2(trnsar,sqt_nabas2int(bas),cdn,dat_ptr,fill_proc);
return(fcnt);
}

int set_chr_trnarr(TRANS_MATRX trnsar,
                   char *cdn,
                   char res)
/* put char into relevent codon positions data struct. Wild card *
  causes all nonzero values to be initialized */
{
RESIDUEDATA rdat;

rdat.trres.reschr = res;
return(set_trnarr(trnsar,cdn,&rdat,set_chr_datum));
}

int set_ini_trnarr(TRANS_MATRX trnsar,
                   char *cdn,
                   int rini)
/* put rini into relevent codon positions data struct. Wild card *
  causes all nonzero values to be initialized */
{
RESIDUEDATA rdat;

rdat.trres.trinit = rini;
return(set_trnarr(trnsar,cdn,&rdat,set_ini_datum));
}

int set_int_trnarr(TRANS_MATRX trnsar,
                   char *cdn,
                   int ival)
/* put ival into relevent codon positions data struct. Wild card *
  causes all nonzero values to be initialized */
{
RESIDUEDATA rdat;

rdat.iresidue = ival;
return(set_trnarr(trnsar,cdn,&rdat,set_int_datum));
}

int set_flt_trnarr(TRANS_MATRX trnsar,
                   char *cdn,
                   float xval)
/* put xval into relevent codon positions data struct. Wild card *
  causes all nonzero values to be initialized */
{
RESIDUEDATA rdat;

rdat.xresidue = xval;
return(set_trnarr(trnsar,cdn,&rdat,set_flt_datum));
}

char sqt_trnslate(TRANS_MATRX trnsar,
                  char *cdn)
/* translate cdn according to trnsar */
{
int p1,p2,p3;

p1 = sqt_nabas2int(*cdn++);
p2 = sqt_nabas2int(*cdn++);
p3 = sqt_nabas2int(*cdn);
return(trnsar[p1][p2][p3].resdata.trres.reschr);
}

int sqt_get_int4tmatrix(TRANS_MATRX trnsar,
                        char *cdn)
/* return integer value for  cdn according to trnsar.
return 0 if not integer datum */
{
int p1,p2,p3;

p1 = sqt_nabas2int(*cdn++);
p2 = sqt_nabas2int(*cdn++);
p3 = sqt_nabas2int(*cdn);
if (trnsar[p1][p2][p3].dtype == SQTDT_int)
  return(trnsar[p1][p2][p3].resdata.iresidue);
else
  return(0);
}

float sqt_get_flt4tmatrix(TRANS_MATRX trnsar,
                          char *cdn)
/* return the float value for translate cdn according to trnsar.
0.0 if not float datum */
{
int p1,p2,p3;

p1 = sqt_nabas2int(*cdn++);
p2 = sqt_nabas2int(*cdn++);
p3 = sqt_nabas2int(*cdn);
if (trnsar[p1][p2][p3].dtype == SQTDT_float)
  return(trnsar[p1][p2][p3].resdata.xresidue);
else
  return(0.0);
}

void sqt_inc_cdncnt(TRANS_MATRX trnsar,
                    char *cdn)
/* increment iresidue for cdn */
{
int p1,p2,p3;

p1 = sqt_nabas2int(*cdn++);
p2 = sqt_nabas2int(*cdn++);
p3 = sqt_nabas2int(*cdn);
trnsar[p1][p2][p3].resdata.iresidue += 1;
}

void init_univ(TRANS_MATRX trnsar)
  /* fill trnsar with codes for universal genetic code */
{
fill_chr_trnarr(trnsar,'?');
fill_ini_trnarr(trnsar,0);
set_chr_trnarr(trnsar,"TTT",'F');
set_chr_trnarr(trnsar,"TTC",'F');
set_chr_trnarr(trnsar,"TTA",'L');
set_chr_trnarr(trnsar,"TTG",'L');
set_chr_trnarr(trnsar,"TC*",'S');
set_chr_trnarr(trnsar,"TAT",'Y');
set_chr_trnarr(trnsar,"TAC",'Y');
set_chr_trnarr(trnsar,"TAA",SQT_TERMCHR);
set_chr_trnarr(trnsar,"TAG",SQT_TERMCHR);
set_chr_trnarr(trnsar,"TGT",'C');
set_chr_trnarr(trnsar,"TGC",'C');
set_chr_trnarr(trnsar,"TGA",SQT_TERMCHR);
set_chr_trnarr(trnsar,"TGG",'W');
set_chr_trnarr(trnsar,"CT*",'L');
set_chr_trnarr(trnsar,"CC*",'P');
set_chr_trnarr(trnsar,"CAT",'H');
set_chr_trnarr(trnsar,"CAC",'H');
set_chr_trnarr(trnsar,"CAA",'Q');
set_chr_trnarr(trnsar,"CAG",'Q');
set_chr_trnarr(trnsar,"CG*",'R');
set_chr_trnarr(trnsar,"ATT",'I');
set_chr_trnarr(trnsar,"ATC",'I');
set_chr_trnarr(trnsar,"ATA",'I');
set_chr_trnarr(trnsar,"ATG",'M');
set_chr_trnarr(trnsar,"AC*",'T');
set_chr_trnarr(trnsar,"AAT",'N');
set_chr_trnarr(trnsar,"AAC",'N');
set_chr_trnarr(trnsar,"AAA",'K');
set_chr_trnarr(trnsar,"AAG",'K');
set_chr_trnarr(trnsar,"AGT",'S');
set_chr_trnarr(trnsar,"AGC",'S');
set_chr_trnarr(trnsar,"AGA",'R');
set_chr_trnarr(trnsar,"AGG",'R');
set_chr_trnarr(trnsar,"GT*",'V');
set_chr_trnarr(trnsar,"GC*",'A');
set_chr_trnarr(trnsar,"GAT",'D');
set_chr_trnarr(trnsar,"GAC",'D');
set_chr_trnarr(trnsar,"GAA",'E');
set_chr_trnarr(trnsar,"GAG",'E');
set_chr_trnarr(trnsar,"GG*",'G');
set_ini_trnarr(trnsar,"ATG",1);
}

void init_trnsar4gb(TRANS_MATRX trnsar,
                    SQT_GENCOD gencod)
/* initialize trnsar to any known genetic code - using Genbank trans-table
  values */
{
init_univ(trnsar);
switch (gencod)
  {
  case GC_vertmito:             /* vertebrate mitochondrial */
    set_chr_trnarr(trnsar,"TGA",'W');
    set_chr_trnarr(trnsar,"ATA",'M');
    set_chr_trnarr(trnsar,"AGA",SQT_TERMCHR);
    set_chr_trnarr(trnsar,"AGG",SQT_TERMCHR);
    set_ini_trnarr(trnsar,"AT*",1);
    break;
  case GC_yeastmito:            /* yeast mitochondrial */
    set_chr_trnarr(trnsar,"TGA",'W');
    set_chr_trnarr(trnsar,"CT*",'T');
    set_chr_trnarr(trnsar,"ATA",'M');
    set_chr_trnarr(trnsar,"CGA",'-');   /* absent */
    set_ini_trnarr(trnsar,"ATG",1);
    set_ini_trnarr(trnsar,"GTG",1);
    break;
  case GC_mycomito:   /* mold, protozoan, coelent, mycoplasma mitochondrial */
    set_chr_trnarr(trnsar,"TGA",'W');
    set_ini_trnarr(trnsar,"ATG",1);
    set_ini_trnarr(trnsar,"GTG",1);
    break;
  case GC_invertmito:           /* invertebrae mitochondrial */
    set_chr_trnarr(trnsar,"TGA",'W');
    set_chr_trnarr(trnsar,"ATA",'M');
    set_chr_trnarr(trnsar,"AGA",'S');
    set_chr_trnarr(trnsar,"AGG",'S');
    set_ini_trnarr(trnsar,"ATG",1);
    set_ini_trnarr(trnsar,"GTG",1);
    set_ini_trnarr(trnsar,"ATT",1);
    break;
  case GC_ciliate:         /* Ciliate Dasycladacean Hexamita Nuclear */
    set_chr_trnarr(trnsar,"TAA",'Q');
    set_chr_trnarr(trnsar,"TAG",'Q');
    set_chr_trnarr(trnsar,"ATG",1);
    break;
  case GC_echinmito:       /* Echinoderm mitochondrial */
    set_chr_trnarr(trnsar,"AAA",'N');
    set_chr_trnarr(trnsar,"AGA",'S');
    set_chr_trnarr(trnsar,"AGG",'S');
    set_chr_trnarr(trnsar,"TGA",'W');
    set_ini_trnarr(trnsar,"ATG",1);
    break;
  case GC_euplotid:        /* Euplotid nuclear */
    set_chr_trnarr(trnsar,"TGA",'C');
    set_ini_trnarr(trnsar,"ATG",1);
    break;
  case GC_bacterial:       /* Bacterial */
    set_ini_trnarr(trnsar,"ATG",1);
    set_ini_trnarr(trnsar,"GTG",1);
    break;
  case GC_yeastalt:        /* Alternative Yeast Nuclear */
    set_chr_trnarr(trnsar,"CTG",'S');
    set_ini_trnarr(trnsar,"ATG",1);
    break;
  case GC_ascidmito:       /* Ascidian mitochondrial */
    set_chr_trnarr(trnsar,"AGA",'G');
    set_chr_trnarr(trnsar,"AGG",'G');
    set_chr_trnarr(trnsar,"ATA",'M');
    set_chr_trnarr(trnsar,"TGA",'W');
    set_ini_trnarr(trnsar,"ATG",1);
    break;
  case GC_flatwmito:       /* Flatworm mitochrondiral */
    set_chr_trnarr(trnsar,"AAA",'N');
    set_chr_trnarr(trnsar,"AGA",'S');
    set_chr_trnarr(trnsar,"AGG",'S');
    set_chr_trnarr(trnsar,"TAA",'Y');
    set_chr_trnarr(trnsar,"TGA",'W');
    set_ini_trnarr(trnsar,"ATG",1);
    break;
  case GC_blepharisma:          /* Blepharisma nuclear */
    set_chr_trnarr(trnsar,"TAG",'Q');
    set_ini_trnarr(trnsar,"ATG",1);
    break;
  case GC_undefined:
  case GC_universal:            /* standard */
  default:
    break;
  }
}

char init_trnsar4chr(TRANS_MATRX trnsar,
                     char gencod)
/* initialize trnsar to any known genetic code, according to following:
u = universal, m = mammalian mitochondrial, y = yeast mitochondrial,
a = aspergillus/Neurospora and Trypansoma mitochondrial,
d = drosophila mitochondrial

return character, if valid, else null. */
{
init_univ(trnsar);
switch (tolower(gencod))
  {
  case 'u':
    return(gencod);
    break;
  case 'm':
    init_trnsar4gb(trnsar,GC_vertmito);
/*    set_chr_trnarr(trnsar,"TGA",'W');
    set_chr_trnarr(trnsar,"ATA",'M');
    set_chr_trnarr(trnsar,"AGG",SQT_TERMCHR);
    set_chr_trnarr(trnsar,"AGG",SQT_TERMCHR); */
    return(gencod);
    break;
  case 'y':
    init_trnsar4gb(trnsar,GC_yeastmito);
/*    set_chr_trnarr(trnsar,"TGA",'W');
    set_chr_trnarr(trnsar,"CTT",'T');
    set_chr_trnarr(trnsar,"CTC",'T');
    set_chr_trnarr(trnsar,"CTA",'T');
    set_chr_trnarr(trnsar,"CTG",'T');
    set_chr_trnarr(trnsar,"ATA",'M'); */
    return(gencod);
    break;
  case 'a':
    set_chr_trnarr(trnsar,"TGA",'W');
    return(gencod);
    break;
  case 'd':
    set_chr_trnarr(trnsar,"TGA",'W');
    set_chr_trnarr(trnsar,"ATA",'M');
    set_chr_trnarr(trnsar,"AGA",'S');
    set_chr_trnarr(trnsar,"AGG",'S');
    return(gencod);
    break;
  default:
    return('\0');
  }
}

int uniquecdn(char *cdn)
  /* return a unique nonzero numeric score for this codon.
  only if all elements are A,C,G,T.  Else return 0 */
{
int ccnt;
int scr;
int bs;

scr = 0;
ccnt = CODONLENGTH;
while (--ccnt)
  if ((bs = sqt_nabas2int(*cdn++)))
    scr = scr * 5 + bs;
  else
    return(0);
return (scr);
}

int validcdn(char *cdn)
  /* return true if this codon is valid.  only if all elements
   are A,C,G,T, or *. */
{
int ccnt;
char bas;

ccnt = CODONLENGTH;
while (--ccnt)
  {
  bas = *cdn++;
  if (! sqt_nabas2int(bas))
    if (bas != '*')
      return(0);
  }
return(1);
}

int sqt_getnxtcdn(char *src,
                  int slen,
                  int *sqpos,
                  char *cdn)
/* fill cdn from sqpos of src, return no of chars written.  
 Stop on null byte or when slen is reached.
 cdn should allow room for null byte.
 *sqpos is 0..n-1 */
{
int cp;

cp = 0;
while ((cp < CODONLENGTH) && (*sqpos < slen) && (*(src+*sqpos) != '\0'))
  {
  *(cdn + cp++) = *(src + *sqpos);
  (*sqpos)++;
  }
*(cdn + cp) = '\0';
return(cp);
}

int sqt_transbuffer(char *src,
                    int slen,
                    TRANS_MATRX gencod,
                    char *dst)
/* translate sequence in src to dst (which will be assumed to be long enough.
  stop on encountering null or slen has been translated */
{
char cdn[CODONLENGTH+1];
int sp;
char *dp;

sp = 0;
dp = dst;
while (sqt_getnxtcdn(src,slen,&sp,&cdn[0]))
  *dp++ = sqt_trnslate(gencod,&cdn[0]);
*dp = '\0';
return((int) (dp - dst));
}

void sqt_inivectarr(COD_VECT_STRCT **cvs)
  /* initialise all cvs aa vectors to null */
{
int vp;

for (vp = 0; vp <= MAXAAINT; vp++)
  *(cvs+vp) = NULL;
}

COD_VECT_STRCT *sqt_appnd_cvect(COD_VECT_STRCT **cvp,
                                CODN_VECTR cv)
/* append this vector to *cvp.  Return a pointer to the new component, NULL
 if problem */
{
COD_VECT_STRCT *prev, *end_ptr;

if (cvp != NULL)
  {
  prev = end_ptr = *cvp;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nextcvec;
    }
  if ((end_ptr = (COD_VECT_STRCT *) getmemory(sizeof(COD_VECT_STRCT),
                                      "Codon vector element")) != NULL)
    {
    end_ptr->nextcvec = NULL;
    bcopy(cv,end_ptr->codnvector,(size_t) sizeof(CODN_VECTR));
    if (*cvp == NULL)
      *cvp = end_ptr;
    else
      prev->nextcvec = end_ptr;
    }
  return(end_ptr);
  }
return(NULL);
}

COD_VECT_STRCT *sqt_appndcdn2vec(COD_VECT_STRCT **cvp,
                                 char *cdn)
/* append a vector for cdn to *cvp assuming it is a valid codon.
  Return a pointer to the new component, NULL if problem */
{
CODN_VECTR cv;

if (sqt_cdn2codvect(cdn,cv))
  return(sqt_appnd_cvect(cvp,cv));
else
  return(NULL);
}

int sqt_bldvects4matrx(TRANS_MATRX tmtx,
                       COD_VECT_STRCT **cvecs)
/* build a set of vectors in cvecs for the translation matrix tmtx.
return the number of vectors added */
{
int vp;
int p1,p2,p3;
CODN_VECTR cv;
int cnt;

cnt = 0;
sqt_inivectarr(cvecs);
for (p1 = 1; p1 <= MAXBASINT; p1++)
  {
  cv[0] = p1;
  for (p2 = 1; p2 <= MAXBASINT; p2++)
    {
    cv[1] = p2;
    for (p3 = 1; p3 <= MAXBASINT; p3++)
      {
      cv[2] = p3;
      vp = aares2int(tmtx[p1][p2][p3].resdata.trres.reschr);
      if (sqt_appnd_cvect(cvecs+vp,cv))
        cnt++;
      }
    }
  }
return(cnt);
}

void sqt_behedvect(COD_VECT_STRCT **cvp)
  /* eat Ist elt of *cvp */
{
COD_VECT_STRCT *lp;

if ((lp = *cvp) != NULL)
  {
  *cvp = lp->nextcvec;
  memfree(lp);
  }
}

void sqt_clrvect(COD_VECT_STRCT **cvp)
  /* clear and free entire linked list *cvp */
{
while (*cvp != NULL)
  sqt_behedvect(cvp);
}

void sqt_clrvectarr(COD_VECT_STRCT **cvarr)
  /* clear all of vectors for all aas */
{
int ap;

for (ap = 0; ap <= MAXAAINT; ap++)
  sqt_clrvect(cvarr+ap);
}

COD_VECT_STRCT *sqt_getvect4res(COD_VECT_STRCT **cvlist,
                                char aares,
                                COD_VECT_STRCT *cvp)
/* if cvp = NULL, then return the first vector in the linked list for residue
aares.  if cvp is non-null, return the next codon vector in the linked list */
{
int rp;

if (cvp == NULL)
  {
  rp = aares2int(aares);
  return(*(cvlist+rp));
  }
else
  return(cvp->nextcvec);
}

COD_VECT_STRCT *sqt_ptr4cvec(COD_VECT_STRCT *clst,
                             CODN_VECTR cv)
/* return a pointer to the element of clst which matches cv.  NULL if none */
{
int cp;
COD_VECT_STRCT *lp;

lp = clst;
while (lp != NULL)
  {
  cp = 0;
  while (cp < CODONLENGTH)
    if (cv[cp] == lp->codnvector[cp])
      {
      if (++cp >= CODONLENGTH) /* have a match, return lp */
        return(lp);
      }
    else       /* mismatch */
      {
      cp = CODONLENGTH;
      lp = lp->nextcvec;
      }
  }
return(NULL);
}

int sqt_tmatival4cvecp(TRANS_MATRX tm,
                       COD_VECT_STRCT *cvp,
                       SQTRN_DATYPE *dtp)
/* simply return the integer value for the
cvp postion of tm.  0 if not an integer value.
if dtp is non-NULL, then return the data type
to it */
{
TRN_MATELT *mep;

if (cvp != NULL)
  {
  mep = sqt_vect2eltaddr(tm,cvp->codnvector);
  if (mep->dtype == SQTDT_int)
    {
    if (dtp != NULL)
      *dtp = SQTDT_int;
    return(mep->resdata.iresidue);
    }
  }
if (dtp != NULL)
  *dtp = SQTDT_undefined;
return(0);
}

COD_VECT_STRCT *sqt_maxintval4cvec(TRANS_MATRX tm,
                                   COD_VECT_STRCT *c1,
                                   COD_VECT_STRCT *c2)
/* compare integer value for tm values of c1,c2, returning
the larger.  if c1 or c2 is NULL, return the other. NULL
if cockup */
{
SQTRN_DATYPE dt1;
SQTRN_DATYPE dt2;

if (c1 == NULL)
  return(c2);
else
  if (c2 == NULL)
    return(c1);
  else
    if ((sqt_tmatival4cvecp(tm,c1,&dt1) >= sqt_tmatival4cvecp(tm,c2,&dt2)) &&
          (dt1 == dt2) && (dt1 == SQTDT_int))
      return(c1);
    else
      if (dt2 == SQTDT_int)
        return(c2);
      else
        return(NULL);
}

COD_VECT_STRCT *sqt_minintval4cvec(TRANS_MATRX tm,
                                   COD_VECT_STRCT *c1,
                                   COD_VECT_STRCT *c2)
/* compare integer value for tm values of c1,c2, returning
the larger.  if c1 or c2 is NULL, return the other. NULL
if cockup */
{
SQTRN_DATYPE dt1;
SQTRN_DATYPE dt2;

if (c1 == NULL)
  return(c2);
else
  if (c2 == NULL)
    return(c1);
  else
    if ((sqt_tmatival4cvecp(tm,c1,&dt1) <= sqt_tmatival4cvecp(tm,c2,&dt2)) &&
          (dt1 == dt2) && (dt1 == SQTDT_int))
      return(c1);
    else
      if (dt2 == SQTDT_int)
        return(c2);
      else
        return(NULL);
}

int sqt_vectrsequal(COD_VECT_STRCT *c1,
                    COD_VECT_STRCT *c2)
/* return 1 if c1 & c2 are the same */
{
int cp;

cp = 0;
while (cp < CODONLENGTH)
  if (c1->codnvector[cp] != c2->codnvector[cp])
    return(0);
return(1);
}

COD_VECT_STRCT *sqt_scnveclst(COD_VECT_STRCT *clst,
                              TRANS_MATRX tm,
                              COD_VECT_STRCT *(*cmpfun)(TRANS_MATRX tx,
                                                        COD_VECT_STRCT *c1,
                                                        COD_VECT_STRCT *c2))
/* Scan thru clst... using tm and cmpfun to decide on best
vector in list.  cmpfun will return which ever of c1,c2 is "better"
*/
{
COD_VECT_STRCT *clp;
COD_VECT_STRCT *best;

clp = clst;
best = clp;
while (clp != NULL)
  {
  best = (*cmpfun)(tm,best,clp->nextcvec);
  clp = clp->nextcvec;
  }
return(best);
}

int sqt_cdnmatb4pos(int cp,
                    char *cdn,
                    COD_VECT_STRCT *cvp)
/* cp is codon position (1..3), cdn is string containing codon so far
 (or NULL) and cvp is pointer to codon structure.
  return true if: cp <= 1 || cdn == NULL
                    || chars of cdn prior to cp match those of *cvp */
{
if (cdn == NULL) 
  return(1);
else
  {
  while (--cp >= 0)
    if (toupper(cdn[cp]) != toupper(sqt_int2nabas(cvp->codnvector[cp])))
      return(0);
  return(1);
  }
}

int sqt_codnbas4res(COD_VECT_STRCT **cvlist,
                    char aares,
                    int cpos,
                    char *cdn,
                    char *cstrng)
/* return a number of base choices for position cpos (1..3) of the
codons for aares.  Return bases in cstrng if non-NULL, 
assumed to be large enough for the purpose.
Error returns: -2 = cpos ivalid, -1 = invalid aares.
If the codon seen so far is present as cdn, then this information is used in
refining the output (Esp for Ser, etc in the standard genetic code */
{
char *sp;
COD_VECT_STRCT *cvp;
char nr;
int bcnt;
char lclstr[MAXBASINT+1];

bcnt = 0;
sp = &lclstr[0];
*sp = '\0';
cpos--;
if ((cpos < 0) || (cpos >= CODONLENGTH))
  return(-2);
else
  if ((cvp = sqt_getvect4res(cvlist,aares,NULL)) == NULL)
    return(-1);    /* couldn't find any codons */
  else
    {
    while (cvp != NULL)
      {
      if (sqt_cdnmatb4pos(cpos,cdn,cvp))
        {
        nr = sqt_int2nabas(cvp->codnvector[cpos]);
        if (index(lclstr,nr) == NULL)
          {
          bcnt++;
          *sp++ = nr;
          *sp = '\0';
          }
        }
      cvp = sqt_getvect4res(cvlist,aares,cvp);
      }
    if (cstrng != NULL)
      (void) strcpy(cstrng,lclstr);
    return(bcnt);
    }
}

COD_VECT_STRCT **sqt_getcvctarr(char *msg)
  /* allocate storage for a MAXAAINT-long array of pointers */
{
return((COD_VECT_STRCT**) getmemory((sizeof(COD_VECT_STRCT*)*(MAXAAINT+1)),
                                      msg));
}

void sqt_sumintmatrx(TRANS_MATRX m1,
                     TRANS_MATRX m2,
                     TRANS_MATRX sum)
/* sum each element of m1 and m2 into sum */
{
TRN_MATELT *p1;
TRN_MATELT *p2;
TRN_MATELT *sp;
CODN_VECTR cv;

for (cv[0] = 0; cv[0] <= MAXBASINT; cv[0]++)
  for (cv[1] = 0; cv[1] <= MAXBASINT; cv[1]++)
    for (cv[2] = 0; cv[2] <= MAXBASINT; cv[2]++)
      {
      p1 = sqt_vect2eltaddr(m1,cv);
      p2 = sqt_vect2eltaddr(m2,cv);
      sp = sqt_vect2eltaddr(sum,cv);
      sp->resdata.iresidue = p1->resdata.iresidue + p2->resdata.iresidue;
      }
}

int sqt_totvldicnts(TRANS_MATRX tmat)
  /* for all valid positions in tmat, sum the integer counts */
{
int tsum;
TRN_MATELT *p1;
CODN_VECTR cv;

tsum = 0;
for (cv[0] = 0; cv[0] <= MAXBASINT; cv[0]++)
  for (cv[1] = 0; cv[1] <= MAXBASINT; cv[1]++)
    for (cv[2] = 0; cv[2] <= MAXBASINT; cv[2]++)
      {
      p1 = sqt_vect2eltaddr(tmat,cv);
      tsum += p1->resdata.iresidue;
      }
return(tsum);
}

int sqt_gcodexfl(FILE *fl,
                 TRANS_MATRX tmat)
/* read CODONS/RES from fl and load into tmat, return no of valid codons read.
fields starting with ';' are regarded as comments */
{
char cdn[33];
char ares[9];
int ccnt;

ccnt = 0;
fill_chr_trnarr(tmat,'?');
fill_ini_trnarr(tmat,0);
while (bas_fgetatm(fl,&cdn[0],32," \n\t") > 0)
  if (cdn[0] == ';')       /* comment line, jump to end of it */ 
    bas_skipeol(fl,NULL);
  else
    {
    if (bas_fgetatm(fl,&ares[0],8," \n\t") > 0)
      if (ares[0] == '+')      /* must note this as an init */
        {
        (void) set_ini_trnarr(tmat,&cdn[0],1);
        ccnt += set_chr_trnarr(tmat,&cdn[0],ares[1]);
        }
      else
        ccnt += set_chr_trnarr(tmat,&cdn[0],ares[0]);
    }
return(ccnt);
}

void say_intcdn(FILE *strm,
                int p1,
                int p2,
                int p3)
/* tell strm the character representation of the codon (p1,p2,p3) */
{
fputc(int2nabas(p1),strm);
fputc(int2nabas(p2),strm);
fputc(int2nabas(p3),strm);
}

void sqt_gcod2fl(FILE *fl,
                 char *desc,
                 TRANS_MATRX tmat,
                 int letter3)
/* write tmat contents to fl in approved form.  desc written as leading
comment if non-NULL, letter3 enables 3 letter AA ids */
{
int p1,p2,p3;
CODN_VECTR cv;
TRN_MATELT *dp;

if (desc != NULL)
  fprintf(fl,"; %s\n\n",desc);
for (p1 = 1; p1 <= MAXBASINT; p1++)
  {
  cv[0] = remapbasptr(p1);
  for (p3 = 1; p3 <= MAXBASINT; p3++)
    {
    cv[2] = remapbasptr(p3);
    for (p2 = 1; p2 <= MAXBASINT; p2++)
      {
      if (p2 > 1)
        fputs("  ",fl);
      cv[1] = remapbasptr(p2);
      say_intcdn(fl,cv[0],cv[1],cv[2]);
      dp = sqt_vect2eltaddr(tmat,cv);
      switch (dp->dtype)
        {
        case SQTDT_char:
          if (dp->resdata.trres.trinit)
            fputs(" +",fl);
          else
            fputs("  ",fl);
          if (letter3)
            fputs(aachr2str3(dp->resdata.trres.reschr),fl);
          else
            fputc(dp->resdata.trres.reschr,fl);
          break;
        case SQTDT_int:
        case SQTDT_float:
        case SQTDT_undefined:
        default:
          break;
        }
      }
    fputc('\n',fl);
    }
  fputc('\n',fl);
  }
}
