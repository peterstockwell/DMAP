/* sqmat_fns.c: sequence matching functions */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "bas_fns.h"
#include "sqmat_fns.h"

char case_mirror(char nc,
                 char casepat)
/* return nc in same case as casepat */
{
if (isupper((int) casepat))
  return((char) toupper((int) nc));
else
  if (islower((int) casepat))
    return((char) tolower((int) nc));
  else
    return(nc);
}

char ssd_bascmplmnt(char bs,
                    BAS_MATMODE mmod)
/* return the complementary base to bs, depending on mmod */
{
char cb;

switch (mmod)
  {
  case BAS_iub:
    if ((cb = ssd_bascmplmnt(bs,BAS_exact)) != '-')
      return(cb);
    else
      switch (tolower(bs))
        {
        case 'r':
          return(case_mirror('Y',bs));
          break;
        case 'y':
          return(case_mirror('r',bs));
          break;
        case 'm':
          return(case_mirror('k',bs));
          break;
        case 'k':
          return(case_mirror('m',bs));
          break;
        case 'n':
        case 's':
        case 'w':
        case '-':
          return(bs);
          break;
        case 'h':
          return(case_mirror('d',bs));
        case 'b':
          return(case_mirror('v',bs));
        case 'v':
          return(case_mirror('b',bs));
        case 'd':
          return(case_mirror('h',bs));
        default:
          return('-');
        }
    break;
  case BAS_ldna:
    if ((cb = ssd_bascmplmnt(bs,BAS_exact)) != '-')
      return(cb);
    else
      switch (tolower(bs))
        {
        case 'r':
          return(case_mirror('Y',bs));
          break;
        case 'y':
          return(case_mirror('r',bs));
          break;
        case 'h':
          return(case_mirror('q',bs));
          break;
        case 'q':
          return(case_mirror('h',bs));
          break;
        case 'p':
        case 's':
        case '-':
          return(bs);
          break;
        case 'k':
          return(case_mirror('l',bs));
          break;
        case 'l':
          return(case_mirror('k',bs));
          break;
        case 'm':
          return(case_mirror('j',bs));
          break;
        case 'j':
          return(case_mirror('m',bs));
          break;
        default:
          return('-');
          break;
        }
    break;
  case BAS_exact:
  default:
    switch (tolower(bs))
      {
      case 'a':
        return(case_mirror('T',bs));
        break;
      case 'c':
        return(case_mirror('G',bs));
        break;
      case 'g':
        return(case_mirror('C',bs));
        break;
      case 't':
      case 'u':
        return(case_mirror('A',bs));
        break;
      default:
        return('-');
        break;
      }
    break;
  }
}

char same_residue(char bas,
                  BAS_MATMODE dummy)
/* return bas */
{
return(bas);
}

void reverse_seq(char *sq,
                 int sqlen,
                 char (* cmpfun)(char xres,
                                 BAS_MATMODE xmod),
                 BAS_MATMODE mmod)
/* reverse sq character order, apply
cmpfun to each transfer */
{
char *hipt;
char *lopt;
char schr;

lopt = sq;
hipt = sq + sqlen - 1;
while (hipt >= lopt)
  {
  schr = *hipt;
  *hipt-- = (* cmpfun)(*lopt,mmod);
  *lopt++ = (* cmpfun)(schr,mmod);
  }
}

void complmnt_seq(char *sq,
                  int sqlen,
                  BAS_MATMODE mmod)
/* reverse&complement the contents of sq, depending on mmod */
{
reverse_seq(sq,sqlen,ssd_bascmplmnt,mmod);
}

int ssd_bas2bit(char bas)
/* return a bit map representing bas with A=bit 1, C=2, G=3 & T=4,
 otherwise 0 */
{
switch ((char) tolower((int) bas))
  {
  case 'a':
    return(1);
    break;
  case 'c':
    return(2);
    break;
  case 'g':
    return(4);
    break;
  case 't':
  case 'u':
    return(8);
    break;
  default:
    return(0);
    break;
  }
}

int ssd_bits4basechr(char bas,
                     BAS_MATMODE mmod)
/* return a bit map representing bas with A=bit 1, C=2, G=3 & T=4 and
  bits set for redundant specs, as for mmod */
{
int bret;

if ((bret = ssd_bas2bit(bas)))  /* if A,C,G,T/U return value right away */
  return(bret);
else
  switch (mmod)               /* not A,C,B,T/U - combine bits */
    {
    case BAS_exact:
      return(0);
      break;
    case BAS_iub:
      switch (toupper(bas))
        {
        case 'R':
          return(ssd_bas2bit('a') | ssd_bas2bit('g'));
          break;
        case 'Y':
          return(ssd_bas2bit('c') | ssd_bas2bit('t'));
          break;
        case 'M':
          return(ssd_bas2bit('a') | ssd_bas2bit('c'));
          break;
        case 'K':
          return(ssd_bas2bit('g') | ssd_bas2bit('t'));
          break;
        case 'S':
          return(ssd_bas2bit('c') | ssd_bas2bit('g'));
          break;
        case 'W':
          return(ssd_bas2bit('a') | ssd_bas2bit('t'));
          break;
        case 'H':
          return(ssd_bas2bit('a') | ssd_bas2bit('c') | ssd_bas2bit('t'));
          break;
        case 'B':
          return(ssd_bas2bit('c') | ssd_bas2bit('g') | ssd_bas2bit('t'));
          break;
        case 'V':
          return(ssd_bas2bit('a') | ssd_bas2bit('c') | ssd_bas2bit('g'));
          break;
        case 'D':
          return(ssd_bas2bit('a') | ssd_bas2bit('g') | ssd_bas2bit('t'));
          break;
        case 'N':
        case '-':
          return(ssd_bas2bit('a') | ssd_bas2bit('g') | ssd_bas2bit('t') |
                   ssd_bas2bit('c'));
          break;
        default:
          return(0);
          break;
        }
      break;
    case BAS_ldna:
      switch (toupper(bas))
        {
        case 'R':
          return(ssd_bas2bit('a') | ssd_bas2bit('g'));
          break;
        case 'Y':
          return(ssd_bas2bit('c') | ssd_bas2bit('t'));
          break;
        case 'S':
          return(ssd_bas2bit('a') | ssd_bas2bit('c'));
          break;
        case 'H':
          return(ssd_bas2bit('g') | ssd_bas2bit('t'));
          break;
        case 'Q':
          return(ssd_bas2bit('c') | ssd_bas2bit('g'));
          break;
        case 'P':
          return(ssd_bas2bit('a') | ssd_bas2bit('t'));
          break;
        case 'L':
          return(ssd_bas2bit('a') | ssd_bas2bit('c') | ssd_bas2bit('t'));
          break;
        case 'J':
          return(ssd_bas2bit('c') | ssd_bas2bit('g') | ssd_bas2bit('t'));
          break;
        case 'M':
          return(ssd_bas2bit('a') | ssd_bas2bit('c') | ssd_bas2bit('g'));
          break;
        case 'K':
          return(ssd_bas2bit('a') | ssd_bas2bit('g') | ssd_bas2bit('t'));
          break;
        case 'N':
        case '-':
          return(ssd_bas2bit('a') | ssd_bas2bit('g') | ssd_bas2bit('t') |
                   ssd_bas2bit('c'));
          break;
        default:
          return(0);
          break;
        }
      break;
    default:
      return(0);
      break;
    }
}

int ssd_basmatch(char b1,
                 char b2,
                 BAS_MATMODE mmod)
/* return 1 if b1 matches b2, given match strategy mmod. case-independent */
{
char lb1;
char lb2;

if ((lb1 = (char) tolower((int) b1)) == (lb2 = (char) tolower((int) b2)))
  return(1);
else
  switch (mmod)
    {
    case BAS_iub:
    case BAS_ldna:
      return((ssd_bits4basechr(lb1,mmod) & ssd_bits4basechr(lb2,mmod)) != 0);
      break;
    case BAS_exact:
    default:
      return(0);
      break;
    }
}

char *ssd_nxtbasmatch(char *seq,
                      char bas,
                      BAS_MATMODE mmod)
/* return a pointer to the next occurrence of bas in seq.  NULL if none */
{
char *sp;

if ((sp = seq))
  {
  while (*sp != '\0')
    if (ssd_basmatch(*sp,bas,mmod))
      return(sp);
    else
      sp++;
  return(NULL);
  }
else
  return(NULL);
}

char *ssd_prvbasmatch(char *seq,
                      char *spos,
                      char bas,
                      BAS_MATMODE mmod)
/* return a pointer to the previous occurrence of bas in seq, starting 
from spos.  NULL if none */
{
char *sp;

if ((sp = spos))
  {
  while (sp >= seq)
    if (ssd_basmatch(*sp,bas,mmod))
      return(sp);
    else
      sp--;
  return(NULL);
  }
else
  return(NULL);
}

char *ssd_nxtstrmatch(char *seq,
                      char *pat,
                      BAS_MATMODE mmod)
/* return the position of the next match of pat in seq. NULL if not found */
{
char *sp;
char *pp;
char *s;
char *sstr;

pp = pat;
sp = seq;
while ((sstr = s = ssd_nxtbasmatch(sp,*pp++,mmod)))
  {
  s++;
  while ((*pp != '\0') && (*s != '\0') && (ssd_basmatch(*s,*pp,mmod)))
    {
    pp++;
    s++;
    }
  if (*pp == '\0') /* found a match, return sp */
    return(sstr);
  else
    if (*s == '\0') /* ran out of seq, return NULL */
      return(NULL);
    else            /* base mismatch, reset and continue */
      {
      sp = sstr + 1;
      pp = pat;
      }
  }
return(NULL);  /* no match if got here, return NULL */
}

char *ssd_prvstrmatch(char *seq,
                      int sqlen,
                      int startpos,
                      char *pat,
                      BAS_MATMODE mmod)
/* return the position of the next match of pat in seq. NULL if not found */
{
char *sp;
char *pp;
char *s;
char *sstr;

pp = pat;
sp = seq + startpos - 1;
while ((sstr = s = ssd_prvbasmatch(seq,sp,*pp++,mmod)))
  {
  s++;
  while ((*pp != '\0') && (*s != '\0') && (ssd_basmatch(*s,*pp,mmod)))
    {
    pp++;
    s++;
    }
  if (*pp == '\0') /* found a match, return sp */
    return(sstr);
  else
    if (*s == '\0') /* ran out of seq, reset and continue, if start in seq */
      if (s >= seq)
        {
        sp = sstr - 1;
        pp = pat;
        }
      else
        return(NULL);
    else            /* base mismatch, reset and continue */
      {
      sp = sstr - 1;
      if (sp < seq)
        return(NULL);
      pp = pat;
      }
  }
return(NULL);  /* no match if got here, return NULL */
}

int ssd_nxtmatchpos(char *seq,
                    int sqlen,
                    int curpos,
                    char *pat,
                    BAS_MATMODE mmod,
                    int srchfwd)
/* return the sequence position of next match after curpos of pat in seq.
  0 if not found.  If string starts with '#' then try to read rest of
string as an integer value and return that */
{
char *nxt;
int npos;

if (curpos <= sqlen)
  {
  if (*pat == '#')   /* numerical argument, decode it and return if valid */
    {
    npos = atoi(pat+1);
    if ((npos <= sqlen) && (npos > 0))
      return(npos);
    else
      return(0);
    }
  else
    {
    if (srchfwd)
      nxt = ssd_nxtstrmatch((seq+curpos-1),pat,mmod);
    else
      nxt = ssd_prvstrmatch(seq,sqlen,curpos,pat,mmod);
    if (nxt)
      if (((npos = (int) (nxt - seq + 1)) <= sqlen) && (npos > 0))
        return(npos);
      else
        return(0);      /* spurious match outside seq */
    else
      return(0);
    }
  }
else
  return(0);
}

int ssd_bitcnt4bas(char ares,
                   BAS_MATMODE mmod)
/* count up the bits set for ares (i.e. no bases it represents) */
{
int bv;
int ccnt;
int msk;
int bc;

bv = ssd_bits4basechr(ares,mmod);
ccnt = 0;
msk = 1;
bc = 0;
while (bc <= 3)
  {
  if (bv & msk)
    ccnt++;
  msk = msk << 1;
  bc++;
  }
return(ccnt);
}
