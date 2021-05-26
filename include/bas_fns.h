/* bas_fns.h: headers for basic library functions in c */

#ifdef MALLOCDBG

typedef struct Malchk     /* for checking mallocs */
  {
  void *ptr;
  struct Malchk *nextmelt;
  struct Malchk *prvmelt;
  }
MALCHK;

MALCHK *scan4free(void *pt);
  /* return ptr to structure containing address pt, NULL if not found */

#endif

void bas_initmalfl(char *nm);
/* initialise mf to write to nm for tracing purposes */

char *getmemory(int nbytes,
                char *msg);
/* malloc nbytes, and return address: die with msg if problem */

char *getmemzero(int nbytes,
                 char *msg);
/* malloc requested memory & zero it */

char *bas_strdup(char *str);
  /* local strdup, intended to be trackable if needed */

void memfree(void *pt);
  /* free malloced pt, tell malfl if appropriate */

void bas_tellunfreed();
  /* report presently unfreed memory - write to malfl if appropriate */

void bas_rptchrout(FILE *fl,
                   int cc,
                   char c);
/* put cc repetitions of c to fl */

char *bas_int2strmnth(int imnth,
                      int shrtfm);
/* return a pointer to an array of month names, 3 letter form for shrtfm */

void bas_scatprintf(char *buf,
                    int blen,
                    char *fmt,
                    ...);
/* varags method of appending (concatenating more material to
buf */

void bas_appdate2str(char *dbuf,
                     char **bp,
                     int blen);
/* append system date & time to dbuf if room at *bp, update *bp
   - code nicked from man strftime pages */

void sqfl_date2file(FILE *sfl);
  /* write system date & time to sfl */

int imin(int i,
         int j);
/* return the minimum of i,j */

int imax(int i,
         int j);
/* return the larger of i,j */

void bas_skipeol(FILE *fl,
                 int *lno);
/* skip past next \n char, increment lno if non-null */

char *bas_fgets(char *buf,
                int blen,
                FILE *fl,
                int *lcnt);
/* perform fgets of fl, remove \n from string */

int bas_chr2ubuf(char *lbuf,
                 char **bp,
                 char nc,
                 int blen);
/* append nc to lbuf at *bp.  return 1 if room */

void bas_appchr(char *lbuf,
                char **bp,
                char nc,
                int blen);
/* append nc to lbuf at *bp, if *bp does not exceed blen chars in total */

void bas_catchr(char *lbuf,
                char nc,
                int blen);
/* append nc to lbuf end of current string, check won't exceed blen chars in 
total */

int bas_str2ubuf(char *lbuf,
                 char **bp,
                 char *str,
                 int blen);
/* append str to lbuf at *bp, if room.  return no of chars inserted */

void bas_appstr(char *lbuf,
                char **bp,
                char *str,
                int blen);
/* append str to lbuf at *bp, if *bp does not exceed blen chars in total */

int bas_int2ubuf(char *lbuf,
                char **bp,
                int ival,
                char *fmt,
                int blen);
/* append ival to lbuf at *bp, return 1 if there was sufficient room */

int bas_fgetatm_ufn(FILE *fl,
                    char *lbuf,
                    int blen,
                    char *brkchrs,
                    int (*ufgetc)(FILE *f));
/* read characters from fl using ufgetc, skipping any in brkchrs.
  write to lbuf, if not overlength, stopping at EOF or next 
occurrence of brkchrs.  Return number of chrs read. 
NULL terminate string if underlength */

int bas_fgetatm(FILE *fl,
                char *lbuf,
                int blen,
                char *brkchrs);
/* read characters from fl, skipping any in brkchrs.  write to lbuf, if
not overlength, stopping at EOF or next occurrence of brkchrs.  Return
number of chrs read. NULL terminate string if underlength */

int bas_fgetlin(FILE *fl,
                char *lbuf,
                int blen,
                int *lcnt);
/* read characters from fl to next eoln or EOF.  Return
number of chrs read, even if 0. NULL terminate string if underlength
return -1 for EOF */

int bas_ptrgetatm(char *sbuf,
                  char **sp,
                  char *lbuf,
                  int blen,
                  char *brkchrs);
/* scan characters from sbuf, starting at *bp, skipping any in brkchrs.  
 update *bp to final position.  write to lbuf, if
not overlength, stopping at '\0' or next occurrence of brkchrs.  Return
number of chrs written. NULL terminate string if underlength */

int bas_sgetatm(char *sbuf,
                char *lbuf,
                int blen,
                char *brkchrs);
/* scan characters from sbuf, skipping any in brkchrs.  write to lbuf, if
not overlength, stopping at '\0' or next occurrence of brkchrs.  Return
number of chrs written. NULL terminate string if underlength */

int bas_digitsin(int nmb);
  /* return the number of decimal digits to represent nmb */

char *bas_strcpyskip(char *lin,
                     int smax,
                     char *dst,
                     int dmax,
                     char *skipset);
/* copy lin up to smax chars to dst, up to dmax chars.  skip leading chars
in skipset, stop on first skipset char.  Return pointer to remainder of string
 - NULL if nothing was found */

char *bas_modstrng(char *str,
                   int (*chrmod)(int nc));
/* replace each character of str by chrmod()ed equivalent.  return str */

int bas_null4nl(int nc);
  /* return nc, unless it is a nl, in which case return '\0' */

char *bas_sharchrs(char *s1,
                   char *s2);
/* return a pointer to any part of s1 that appears in s2, NULL if none */

int bas_cntfputs(char *str,
                 FILE *fl);
/* return the length of this string, and put to fl */

int bas_coutputchr(FILE *fl,
                   char c);
/* write c to fl, in "approved form. return no of chars to do this */

int bas_putochr(FILE *fl,
                char c);
/* write c to fl, in octal form. return no of chars to do this */

int bas_coutputstr(FILE *fl,
                   char *str,
                   char dlmt,
                   int ccnt,
                   int (*coutfun)(FILE *xfl,
                                  char c));
/* put out ccnt characters from str, irrespective of their value.  Return
actual no of char positions consumed */

int bas_strncmp(char *s1,
                char *s2,
                int n);
/* compare s1 & s2 byte by byte up to n bytes or till either differ.
  return 0 if same, -1 if s1 precedes s2 in collating order, +1 if succeeds.
Differs from library strncmp() in that this will continue beyond null chars */

int bas_cmatcasdep(char c1,
                   char c2);
/* return 1 if c1 == c2, case dependent */

int bas_cmatnocas(char c1,
                  char c2);
/* return 1 if c1 == c2, case independent */

int bas_wldstrcmp(char *strng,
                  char *patn,
                  int (*chrcmpfn)(char c1,
                                  char c2));
/* simple minded string matching function which will match chars using 
(*chrcmpfn)() (returns 1 for match, 0 otherwise), and '*'s.
  return 1 for match 0 otherwise */

int bas_getustr(FILE *fl,
                char *ubuf,
                int ublen);
/* get user input - read chars, don't overflow ubuf.
return length of string returned, -1 for EOF */

int bas_ugetstr(char *prmpt,
                char *ubuf,
                int ublen);
/* prompt for user input with prmpt - read chars, don't overflow ubuf.
return length of string returned */

void bas_pause();
  /* request user <CR> */

int bas_getuint(char *prmpt,
                int min,
                int max,
                int def,
                int *uval);
/* prompt stdout/stdin for an integer value generating a prompt from prmpt
range and default values.  return returned value in uval and 1 if decoded OK,
else return 0 */

int bas_uconfirm(char *prmpt,
                 char def);
/* prompt for a yes/no response */
