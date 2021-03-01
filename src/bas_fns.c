/* bas_fns.c: basic library functions written in c */

#include <stdlib.h>
#include <time.h>
#include <locale.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "bas_fns.h"

#ifdef MALLOCDBG

FILE *malfl = NULL;

MALCHK *mfirst;
MALCHK *mlast;
int mlccnt;

MALCHK *scan4free(void *pt)
  /* return ptr to structure containing address pt, NULL if not found */
{
MALCHK *lp;

lp = mfirst;
while (lp != NULL)
  if (lp->ptr == pt)
    return(lp);
  else
    lp = lp->nextmelt;
return(NULL);
}

void bas_initmalfl(char *nm)
/* initialise mf to write to nm for tracing purposes */
{
if ((malfl = fopen(nm,"w")) == NULL)
  {
  fprintf(stderr,"Can't open malloc stream file %s\n",nm);
  exit(1);
  }
mfirst = mlast = NULL;
mlccnt = 0;
}

void mchk_addpt(void *pt)
  /* add pt to list of stored items */
{
MALCHK *newp;

if ((newp = (MALCHK *) malloc(sizeof(MALCHK))) == NULL)
  {
  fprintf(stderr,"check element malloc failed\n");
  if (malfl)
    fflush(malfl);
  exit(1);
  }
if (mlast != NULL)
  {
  mlast->nextmelt = newp;
  newp->prvmelt = mlast;
  }
else
  {
  mfirst = newp;
  newp->prvmelt = NULL;
  }
newp->ptr = pt;
newp->nextmelt = NULL;
mlccnt++;
mlast = newp;
}

#else

void bas_initmalfl(char *nm)
  /* dummy under these conditions */
{
}

#endif

char *getmemory(int nbytes,
                char *msg)
/* malloc nbytes, and return address: die with msg if problem */
{
char *pt;

if ((pt = malloc(nbytes)) == NULL)
  {
  fprintf(stderr,"Malloc failed for %d: %s\n",nbytes,msg);
  exit(1);
  }
#ifdef MALLOCDBG
if (malfl)
  {
  fprintf(malfl,"%d @ %lx: %s\n",nbytes,(long) pt,msg);
  fflush(malfl);
  }
mchk_addpt(pt);
#endif
return(pt);
}

char *getmemzero(int nbytes,
                 char *msg)
/* malloc requested memory & zero it */
{
char *mp;

mp = getmemory(nbytes,msg);
bzero(mp,(size_t) nbytes);
return(mp);
}


char *bas_strdup(char *str)
  /* local strdup, intended to be trackable if needed */
{
#ifdef MALLOCDBG
int slen;
char *nmem;

if (str == NULL)
  return(NULL);
else
  {
  slen = strlen(str) + 1;
  nmem = getmemory(slen,"bas_strdup");
  bcopy(str,nmem,(size_t) slen);
  return(nmem);
  }
#else
if (str == NULL)
  return(NULL);
else
  return(strdup(str));
#endif
}

void memfree(void *pt)
  /* free malloced pt, tell malfl if appropriate */
{
#ifdef MALLOCDBG
MALCHK *pts;
MALCHK *adjp;

if ((pts = scan4free(pt)) == NULL)
  {
  fprintf(stderr,"free nonalloced memory %lx\n",(long) pt);
  fflush(stderr);
  if (malfl)
    {
    fprintf(malfl,"Free nonalloced memory %lx\n",(long) pt);
    fflush(malfl);
    }
/*  free((void *) -1); */  /* try to engineer a traceable dump */
  exit(1);
  }
else
  {
  if ((adjp = pts->nextmelt) == NULL) /* last element */
    {
    mlast = pts->prvmelt;
    if (mlast != NULL)
      mlast->nextmelt = NULL;
    }
  else
    adjp->prvmelt = pts->prvmelt;
  if ((adjp = pts->prvmelt) == NULL) /* first element */
    {
    mfirst = pts->nextmelt;
    if (mfirst != NULL)
      mfirst->prvmelt = NULL;
    }
  else
    adjp->nextmelt = pts->nextmelt;
  free(pts);
  mlccnt--;
  if (malfl)
    {
    fprintf(malfl,"Free %lx\n",(long) pt);
    fflush(malfl);
    }
  }

#endif

free(pt);
}

void bas_tellunfreed()
  /* report presently unfreed memory - write to malfl if appropriate */
{
#ifdef MALLOCDBG
MALCHK *mp;
FILE *fl;

mp = mfirst;
if (malfl != NULL)
  fl = malfl;
else
  fl = stderr;
while (mp != NULL)
  {
  fprintf(fl,"Unfreed: %lx\n",(long) mp->ptr);
  mp = mp->nextmelt;
  }
#else
/* do nothing */
#endif
}

void bas_rptchrout(FILE *fl,
                   int cc,
                   char c)
/* put cc repetitions of c to fl */
{
while (cc-- > 0)
  fputc(c,fl);
}

char *bas_int2strmnth(int imnth,
                      int shrtfm)
/* return a pointer to an array of month names, 3 letter form for shrtfm */
{
if (shrtfm)
  switch (imnth)
    {
    case 1:
      return("Jan");
      break;
    case 2:
      return("Feb");
      break;
    case 3:
      return("Mar");
      break;
    case 4:
      return("Apr");
      break;
    case 5:
      return("May");
      break;
    case 6:
      return("Jun");
      break;
    case 7:
      return("Jul");
      break;
    case 8:
      return("Aug");
      break;
    case 9:
      return("Sep");
      break;
    case 10:
      return("Oct");
      break;
    case 11:
      return("Nov");
      break;
    case 12:
      return("Dec");
      break;
    default:
      return("?Mnth?");
      break;
    }
else
  switch (imnth)
    {
    case 1:
      return("January");
      break;
    case 2:
      return("February");
      break;
    case 3:
      return("March");
      break;
    case 4:
      return("April");
      break;
    case 5:
      return("May");
      break;
    case 6:
      return("June");
      break;
    case 7:
      return("July");
      break;
    case 8:
      return("August");
      break;
    case 9:
      return("Sepember");
      break;
    case 10:
      return("October");
      break;
    case 11:
      return("November");
      break;
    case 12:
      return("December");
      break;
    default:
      return("?Month?");
      break;
    }
}

void bas_scatprintf(char *buf,
                    int blen,
                    char *fmt,
                    ...)
/* varags method of appending (concatenating more material to
buf */
{
va_list args;
char *fstrng;

va_start(args,fmt);
if (vasprintf(&fstrng,fmt,args) > 0)
  {
  (void) strncat(buf,fstrng,(size_t) (blen - strlen(buf)));
  free(fstrng);
  }
va_end(args);
}

void bas_appdate2str(char *dbuf,
                     char **bp,
                     int blen)
/* append system date & time to dbuf if room at *bp, update *bp
   - code nicked from man strftime pages */
{
char dtstr[80];
time_t dtbin;
struct tm *dtstrct;

(void) setlocale(LC_ALL, "");
if (time(&dtbin) == (time_t) -1)
  bas_appstr(dbuf,bp,"time call failed",blen);
else
  {
  dtstrct = localtime(&dtbin);
  if (strftime(dtstr,80, "%a %e-%b-%Y %T",dtstrct) == (size_t) 0)
    bas_appstr(dbuf,bp,"time call failed",blen);
  else
    bas_appstr(dbuf,bp,&dtstr[0],blen);
  }
} 

void sqfl_date2file(FILE *sfl)
  /* write system date & time to sfl */
{
char dtstr[80];
char *bp;

bp = &dtstr[0],
bas_appdate2str(&dtstr[0],&bp,79);
fprintf(sfl,"%s",&dtstr[0]);
} 

int imin(int i,
         int j)
/* return the minimum of i,j */
{
if (i <= j)
  return(i);
else
  return(j);
}

int imax(int i,
         int j)
/* return the larger of i,j */
{
return(i>j ? i:j);
}

void bas_skipeol(FILE *fl,
                 int *lno)
/* skip past next \n char or \r, increment lno if non-null */
{
int nc;

while (((nc = fgetc(fl)) != EOF) && ((char) nc != '\n') &&
          ((char) nc != '\r'));
  if (lno != NULL)
    (*lno)++;
}

char *bas_fgets(char *buf,
                int blen,
                FILE *fl,
                int *lcnt)
/* perform fgets of fl, remove \n from string */
{
char *lp;
char *slim;
char *rslt;

rslt = fgets(buf,(blen+1),fl);
if (rslt != NULL)
  {
  lp = buf;
  slim = buf + blen;
  while ((lp <= slim) && (*lp != '\0'))
    {
    if (*lp == '\n')
      {
      (*lcnt)++;
      *lp = '\0';
      lp = slim;
      }
    lp++;
    }
  if (*buf == '\0') /* have read an empty line in effect, try again */
    return(bas_fgets(buf,blen,fl,lcnt));
  else
    return(rslt);
  }
else
  return(NULL);
}

int bas_chr2ubuf(char *lbuf,
                 char **bp,
                 char nc,
                 int blen)
/* append nc to lbuf at *bp.  return 1 if room */
{
if (*bp - lbuf <= blen)
  {
  **bp = nc;
  (*bp)++;
  if (*bp - lbuf <= blen)  /* try to ensure there is a terminal null chr */
    **bp = '\0';
  return(1);
  }
else
  return(0);
}

void bas_appchr(char *lbuf,
                char **bp,
                char nc,
                int blen)
/* append nc to lbuf at *bp, if *bp does not exceed blen chars in total */
{
(void) bas_chr2ubuf(lbuf,bp,nc,blen);
}

void bas_catchr(char *lbuf,
                char nc,
                int blen)
/* append nc to lbuf end of current string, check won't exceed blen chars in 
total */
{
char *bp;

bp = lbuf;
while (*bp != '\0')
  bp++;
bas_appchr(lbuf,&bp,nc,blen);
}

int bas_str2ubuf(char *lbuf,
                 char **bp,
                 char *str,
                 int blen)
/* append str to lbuf at *bp, if room.  return no of chars inserted */
{
char *sp;
int ccnt;

ccnt = 0;
sp = str;
while ((sp != NULL) && (*sp != '\0'))
  if (bas_chr2ubuf(lbuf,bp,*sp++,blen))
    ccnt++;
return(ccnt);
}

void bas_appstr(char *lbuf,
                char **bp,
                char *str,
                int blen)
/* append str to lbuf at *bp, if *bp does not exceed blen chars in total */
{
(void) bas_str2ubuf(lbuf,bp,str,blen);
}

int bas_int2ubuf(char *lbuf,
                 char **bp,
                 int ival,
                 char *fmt,
                 int blen)
/* append ival to lbuf at *bp, return 1 if there was sufficient room */
{
char ibuf[33];
int slen;

sprintf(&ibuf[0],fmt,ival);
slen = strlen(&ibuf[0]);
if ((*bp - lbuf) <= (blen - slen))
  return(bas_str2ubuf(lbuf,bp,&ibuf[0],blen) > 0);
else
  return(0);
}

int bas_fgetatm_ufn(FILE *fl,
                    char *lbuf,
                    int blen,
                    char *brkchrs,
                    int (*ufgetc)(FILE *f))
/* read characters from fl using ufgetc, skipping any in brkchrs.
  write to lbuf, if not overlength, stopping at EOF or next 
occurrence of brkchrs.  Return number of chrs read. 
NULL terminate string if underlength */
{
int nc;
char *bp;
int rv;

bp = lbuf;
while (((nc = (*ufgetc)(fl)) != EOF) && (index(brkchrs,(char)nc) != NULL));
if (nc == EOF)
  return(-1);
else
  {
  bas_appchr(lbuf,&bp,(char) nc,blen);
  while (((nc = (*ufgetc)(fl)) != EOF) && (index(brkchrs,(char) nc) == NULL))
    bas_appchr(lbuf,&bp,(char) nc,blen);
  rv = (int) (bp - lbuf);
  bas_appchr(lbuf,&bp,'\0',blen);
  return(rv);
  }
}

int bas_fgetatm(FILE *fl,
                char *lbuf,
                int blen,
                char *brkchrs)
/* read characters from fl, skipping any in brkchrs.  write to lbuf, if
not overlength, stopping at EOF or next occurrence of brkchrs.  Return
number of chrs read. NULL terminate string if underlength */
{
return(bas_fgetatm_ufn(fl,lbuf,blen,brkchrs,fgetc));
}

int bas_fgetlin(FILE *fl,
                char *lbuf,
                int blen,
                int *lcnt)
/* read characters from fl to next eoln, \r or EOF.  Return
number of chrs read, even if 0. NULL terminate string if underlength
return -1 for EOF */
{
int nc;
char *bp;
int rv;

bp = lbuf;
while (((nc = fgetc(fl)) != EOF) && ((char) nc != '\n') && ((char) nc != '\r'))
  bas_appchr(lbuf,&bp,(char) nc,blen);
rv = (int)(bp - lbuf);
bas_appchr(lbuf,&bp,'\0',blen);
(*lcnt)++;
if (nc == EOF)
  return(-1);
else
  return(rv);
}

int bas_ptrgetatm(char *sbuf,
                  char **sp,
                  char *lbuf,
                  int blen,
                  char *brkchrs)
/* scan characters from sbuf, starting at *bp, skipping any in brkchrs.  
 update *bp to final position.  write to lbuf, if
not overlength, stopping at '\0' or next occurrence of brkchrs.  Return
number of chrs written. NULL terminate string if underlength */
{
char nc;
char *bp;
int rv;

bp = lbuf;
while (((nc = *(*sp)++) != '\0') && (index(brkchrs,nc) != NULL));
if (nc == '\0')
  return(0);
else
  {
  bas_appchr(lbuf,&bp,nc,blen);
  while (((nc = *(*sp)++) != '\0') && (index(brkchrs,nc) == NULL))
    bas_appchr(lbuf,&bp,nc,blen);
  rv = (int) (bp - lbuf);
  bas_appchr(lbuf,&bp,'\0',blen);
  return(rv);
  }
}

int bas_sgetatm(char *sbuf,
                char *lbuf,
                int blen,
                char *brkchrs)
/* scan characters from sbuf, skipping any in brkchrs.  write to lbuf, if
not overlength, stopping at '\0' or next occurrence of brkchrs.  Return
number of chrs written. NULL terminate string if underlength */
{
char *sp;

sp = sbuf;
return(bas_ptrgetatm(sbuf,&sp,lbuf,blen,brkchrs));
}

int bas_digitsin(int nmb)
  /* return the number of decimal digits to represent nmb */
{
int cnt;

if (nmb == 0)
  return(1);
else
  if (nmb < 0)
    return(bas_digitsin(-nmb) + 1);
  else
    {
    cnt = 0;
    while (nmb > 0)
      {
      nmb = (int) nmb/10;
      cnt++;
      }
    return(cnt);
    }
}

char *bas_strcpyskip(char *lin,
                     int smax,
                     char *dst,
                     int dmax,
                     char *skipset)
/* copy lin up to smax chars to dst, up to dmax chars.  skip leading chars
in skipset, stop on first skipset char.  Return pointer to remainder of string
 - NULL if nothing was found */
{
char *sp;
char *dp;
int ocnt;

sp = lin;
while (index(skipset,*sp) != NULL)
  {
  sp++;
  smax--;
  }
dp = dst;
ocnt = 0;
while ((smax-- > 0) && (*sp != '\0') && (dmax-- > 0) && 
        (index(skipset,*sp) == NULL))
  {
  if (dp != NULL)
    *dp++ = *sp++;
  else
    sp++;
  ocnt++;
  }
if ((dp != NULL) && (dmax > 0))
  *dp = '\0';
if (ocnt > 0)
  return(sp);
else
  return(NULL);
}

char *bas_modstrng(char *str,
                   int (*chrmod)(int nc))
/* replace each character of str by chrmod()ed equivalent.  return str */
{
char *sp;

sp = str;
while (*sp != '\0')
  {
  *sp = (char)(*chrmod)((int) *sp);
  sp++;
  }
return(str);
}

int bas_null4nl(int nc)
  /* return nc, unless it is a nl, in which case return '\0' */
{
if (nc == '\n')
  return('\0');
else
  return(nc);
}

char *bas_sharchrs(char *s1,
                   char *s2)
/* return a pointer to any part of s1 that appears in s2, NULL if none */
{
char *p1;

p1 = s1;
while (*p1 != '\0')
  if (index(s2,*p1) != NULL)   /* have a match, exit here */
    return(p1);
  else
    p1++;
return(NULL);                  /* fell off end, no match, return NULL */
}

int bas_cntfputs(char *str,
                 FILE *fl)
/* return the length of this string, and put to fl */
{
fputs(str,fl);
return(strlen(str));
}

int bas_coutputchr(FILE *fl,
                   char c)
/* write c to fl, in "approved form. return no of chars to do this */
{
switch (c)
  {
  case '\0':
    return(bas_cntfputs("\\0",fl));
    break;
  case '\n':
    return(bas_cntfputs("\\n",fl));
    break;
  case '\r':
    return(bas_cntfputs("\\r",fl));
    break;
  case '\t':
    return(bas_cntfputs("\\t",fl));
    break;
  default:
    if ((c >= ' ') && (c <= '~'))
      {
      fputc(c,fl);
      return(1);
      }
    else
      {
      fprintf(fl,"\\%03o",(int) c);
      return(4);
      }
    break;
  }
}

int bas_putochr(FILE *fl,
                char c)
/* write c to fl, in octal form. return no of chars to do this */
{
fprintf(fl,"\\%03o",(int) 0xff&c);
return(4);
}

int bas_coutputstr(FILE *fl,
                   char *str,
                   char dlmt,
                   int ccnt,
                   int (*coutfun)(FILE *xfl,
                                  char c))
/* put out ccnt characters from str, irrespective of their value.  Return
actual no of char positions consumed */
{
char *sp;
int cc;

cc = 0;
sp = str;
if (dlmt != '\0')
  {
  fputc(dlmt,fl);
  cc++;
  }
while (ccnt-- > 0)
  cc += (*coutfun)(fl,*sp++);
if (dlmt != '\0')
  {
  fputc(dlmt,fl);
  cc++;
  }
return(cc);
}

int bas_strncmp(char *s1,
                char *s2,
                int n)
/* compare s1 & s2 byte by byte up to n bytes or till either differ.
  return 0 if same, -1 if s1 precedes s2 in collating order, +1 if succeeds.
Differs from library strncmp() in that this will continue beyond null chars */
{
char *p1;
char *p2;

p1 = s1;
p2 = s2;
while (n-- > 0)
  {
  if (*p1 != *p2)
    if (*p1 > *p2)
      return(1);
    else
      return(-1);
  else
    {
    p1++;
    p2++;
    }
  }
return(0);
}

int bas_cmatcasdep(char c1,
                   char c2)
/* return 1 if c1 == c2, case dependent */
{
return(c1 == c2);
}

int bas_cmatnocas(char c1,
                  char c2)
/* return 1 if c1 == c2, case independent */
{
return(toupper(c1) == toupper(c2));
}

int bas_wldstrcmp(char *strng,
                  char *patn,
                  int (*chrcmpfn)(char c1,
                                  char c2))
/* simple minded string matching function which will match chars using 
(*chrcmpfn)() (returns 1 for match, 0 otherwise), and '*'s.
  return 1 for match 0 otherwise */
{
char *pp;
char *sp;
char *wp;

pp = patn;
sp = strng;
if (*pp == '\0')       /* at end of pattern - success if at end of string */
  return(*sp == '\0');
else
  if (*pp != '*')
    if ((*chrcmpfn)(*pp,*sp))
      return(bas_wldstrcmp(++sp,++pp,chrcmpfn));
    else
      return(0);
  else   /* on '*' in pattern */
    {
    while (*pp == '*')
      wp = pp++;    /* skip ahead of '*' */
    if (*pp == '\0')     /* we've won */
      return(1);
    else
      {
      while ((*sp != '\0') && (!(*chrcmpfn)(*sp,*pp)))
        sp++;
      if (*sp == '\0')
        return(0);         /* next c doesn't match */
      else
        return(bas_wldstrcmp(++sp,++pp,chrcmpfn) ||
                 bas_wldstrcmp(sp,wp,chrcmpfn));
      }
    }
}

int bas_getustr(FILE *fl,
                char *ubuf,
                int ublen)
/* get input from fl - read chars, don't overflow ubuf.
return length of string returned, -1 for EOF */
{
char *bp;
int nc;

bp = ubuf;
while (((nc = fgetc(fl)) != EOF) && ((char) nc != '\n'))
  bas_appchr(ubuf,&bp,(char) nc,ublen);
if (nc == EOF)
  return(-1);
else
  return((int) (bp - ubuf));
}

int bas_ugetstr(char *prmpt,
                char *ubuf,
                int ublen)
/* prompt for user input with prmpt - read chars, don't overflow ubuf.
return length of string returned */
{
char *bp;
int nc;

if (prmpt != NULL)
  fputs(prmpt,stdout);
fflush(stdout);
bp = ubuf;
while (((nc = fgetc(stdin)) != EOF) && ((char) nc != '\n'))
  bas_appchr(ubuf,&bp,(char) nc,ublen);
return((int)(bp - ubuf));
}

void bas_pause()
  /* request user <CR> */
{
int nc;

fprintf(stdout,"  Press <Return> to continue:");
fflush(stdout);
while (((nc = fgetc(stdin)) != EOF) && ((char) nc != '\n'));
}

int bas_getuint(char *prmpt,
                int min,
                int max,
                int def,
                int *uval)
/* prompt stdout/stdin for an integer value generating a prompt from prmpt
range and default values.  return returned value in uval and 1 if decoded OK,
else return 0 */
{
char *pbuf;
char *ubuf;
int iv;
char *ep;
int rv;

pbuf = (char *) getmemory((strlen(prmpt)+40),"prompt buffer");
ubuf = (char *) getmemory(129,"user buff");
sprintf(pbuf,prmpt,min,max,def);
if (bas_ugetstr(pbuf,ubuf,128) == 0)
  iv = def;
else
  iv = (int) strtol(ubuf,&ep,10);
if ((iv >= min) && (iv <= max))
  {
  *uval = iv;
  rv = 1;
  }
else
  rv = 0;
memfree(ubuf);
memfree(pbuf);
return(rv);
}

int bas_uconfirm(char *prmpt,
                 char def)
/* prompt for a yes/no response */
{
char *ubuf;
char *pbuf;
int rv;

pbuf = (char *) getmemory((strlen(prmpt)+20),"prompt buff");
ubuf = (char *) getmemory(129,"user buf");
sprintf(pbuf,"%s (Y/N) [%c] > ",prmpt,def);
do
  if (bas_ugetstr(pbuf,ubuf,128) == 0)
    *ubuf = def;
while ((toupper(*ubuf) != 'Y') && (toupper(*ubuf) != 'N'));
rv = toupper(*ubuf) == 'Y';
memfree(pbuf);
memfree(ubuf);
return(rv);
}
