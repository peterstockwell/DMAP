/* wlu_fns.c: fast word lookup functions for various purposes */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "bas_fns.h"
#include "sqfl_fns.h"
#include "wlu_fns.h"

int wlu_init2offst(int cased,
                   char initc)
/* return a numerical offset for initc, depending on
 cased:  0 => case depend, 1 case independ, -1 wild card */
{
switch (cased)
  {
  case WLU_CASEIND:
  case WLU_CASEWILD:
    if (islower(initc))         /* lowercase, upcase it */
      return(wlu_init2offst(cased,(char) toupper(initc)));
    else
       if (isupper(initc))
         return((int) initc - 'A' + 1);
       else   /* other */
         return(0);
    break;
  case WLU_CASEDEP:       /* 0,A..Z,a..z counted */
  default:
    if (isalpha(initc))
      if (isupper(initc))
        return((int) initc - 'A' + 1);
      else
        return(wlu_init2offst(cased,(char) toupper(initc)) +
                 wlu_init2offst(cased,'Z'));
    else
      return(0);
    break;
  }
}

int wlu_no_inits(int cased)
  /* return the number of initials for
  cased: 0 => case depend, 1 case independ, -1 wild card */
{
return(wlu_init2offst(cased,'z')+1);
}

void wlu_initlustrctptr(WRD_LUSTRCT *wlus,
                        int cased,
                        void *failval)
/*  init wlus values prior to loading words
  cased: 0 => case depend, 1 case independ, -1 wild card.
This is for initing pointer based lists */
{
int ninits;
int wp;
WRD_LU_REC **lp;
WRD_LU_REC **fp;

wlus->casedep = cased;
wlus->failret = failval;
ninits = wlu_no_inits(cased);
fp = wlus->firstlet = (WRD_LU_REC **)
                           getmemory(((ninits+1)*sizeof(WRD_LU_REC *)),
                                       "Word first letter look up table");
lp = wlus->lastlet = (WRD_LU_REC **)
                           getmemory(((ninits+1)*sizeof(WRD_LU_REC *)),
                                       "Word last letter look up table");
for (wp = 0; wp < ninits; wp++)
  {
  *(fp++) = NULL;
  *(lp++) = NULL;
  }
}

void wlu_initlustrct(WRD_LUSTRCT *wlus,
                     int cased,
                     int failval)
/*  init wlus values prior to loading words
  cased: 0 => case depend, 1 case independ, -1 wild card */
{
int ninits;
int wp;
WRD_LU_REC **lp;
WRD_LU_REC **fp;

wlus->casedep = cased;
wlus->failret = NULL;
wlus->failret = (void *) ((long) failval);
ninits = wlu_no_inits(cased);
fp = wlus->firstlet = (WRD_LU_REC **)
                           getmemory(((ninits+1)*sizeof(WRD_LU_REC *)),
                                       "Word first letter look up table");
lp = wlus->lastlet = (WRD_LU_REC **)
                           getmemory(((ninits+1)*sizeof(WRD_LU_REC *)),
                                       "Word last letter look up table");
for (wp = 0; wp < ninits; wp++)
  {
  *(fp++) = NULL;
  *(lp++) = NULL;
  }
}

WRD_LUSTRCT *wlu_getlustrct(int cased,
                            int failval)
/* allocate a wlu strcture and initialise it
with case & failval */
{
WRD_LUSTRCT *newwlu;

newwlu = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"wordlookup");
wlu_initlustrct(newwlu,cased,failval);
return(newwlu);
}

void wlu_killwelt(WRD_LU_REC *rp)
  /* kill this element, lose storage associated with it.  This does not
change forward links - that must be handled by calling routine */
{
memfree(rp->wrd);
if (rp->helpline != NULL)
  memfree(rp->helpline);
memfree(rp);
}

void wlu_clrlustrct(WRD_LUSTRCT *wlu)
  /* relinquish malloced memory for wlu lists */
{
int lp;
int lmax;
WRD_LU_REC *rp;
WRD_LU_REC *np;

if (wlu != NULL)
  {                   /* free all allocated word records */
  lmax = wlu_no_inits(wlu->casedep);
  for (lp = 0; lp < lmax; lp++)
    {
    rp = *(wlu->firstlet+lp);
    while (rp != NULL)
      {
      np = rp->nextiwrd;
      wlu_killwelt(rp);
      rp = np;
      }
    }
  memfree(wlu->firstlet);       /* free arrays of starts */
  wlu->firstlet = NULL;
  memfree(wlu->lastlet);        /* free arrays of stops */
  wlu->lastlet = NULL;
  }
}

void wlu_clralllustrct(WRD_LUSTRCT *wlu)
  /* relinquish malloced memory for wlu lists and wlu itself */
{
if (wlu != NULL)
  {                   /* free all allocated word records */
  wlu_clrlustrct(wlu);
  memfree(wlu);                 /* free the thing itself - late change!! */
  }
}

void wlu_addwrdptr(WRD_LUSTRCT *wlus,
                   char *newwrd,
                   void *newval,
                   char *hlpline)
/* add newwrd to lookup structure.  No check is performed for multiple
  insertions */
{
int offst;
WRD_LU_REC *newrec;

offst = wlu_init2offst(wlus->casedep,*newwrd);
newrec = (WRD_LU_REC *) getmemory(sizeof(WRD_LU_REC),
                                    "New word lookup element");
newrec->wrd = bas_strdup(newwrd);
newrec->retval = newval;
newrec->nextiwrd = NULL;
if (hlpline != NULL)
  newrec->helpline = bas_strdup(hlpline);
else
  newrec->helpline = NULL;
if ((*(wlus->lastlet+offst)) == NULL) /* haven't seen this one before */
  (*(wlus->firstlet+offst)) = (*(wlus->lastlet+offst)) = newrec;
else
  {
  (*(wlus->lastlet+offst))->nextiwrd = newrec;
  *(wlus->lastlet+offst) = newrec;
  }
}

void wlu_addwrd(WRD_LUSTRCT *wlus,
                char *newwrd,
                int newval,
                char *hlpline)
/* add newwrd to lookup structure.  No check is performed for multiple
  insertions */
{
wlu_addwrdptr(wlus,newwrd,(void *) ((long) newval),hlpline);
}

void wlu_addintgrs(WRD_LUSTRCT *wlus,
                   int newint,
                   int newval,
                   char *hlpline)
/* add newint as a string to lookup structure, returning newval.
  No check is performed for multiple insertions */
{
char nbuf[33];

sprintf(&nbuf[0],"%d",newint);
wlu_addwrd(wlus,&nbuf[0],newval,hlpline);
}

void wlu_addintgr(WRD_LUSTRCT *wlus,
                  int newval,
                  char *hlpline)
/* add newval as a string to lookup structure, returning newval.
  No check is performed for multiple insertions */
{
wlu_addintgrs(wlus,newval,newval,hlpline);
}

int wlu_cntwldmats(WRD_LU_REC *wrpt,
                   char *uwrd,
                   int wlen)
/* return count of words which match this one upto wlen chars in linked
  list wrpt */
{
int wc;

wc = 0;
while (wrpt != NULL)
  {
  if (strncasecmp(uwrd,wrpt->wrd,(size_t) wlen) == 0)
    wc++;
  wrpt = wrpt->nextiwrd;
  }
return(wc);
}

WRD_LU_REC *wlu_lulst4wrd(WRD_LUSTRCT *wlus,
                          char *uwrd)
/* return the start of the list containing this word */
{
return(*(wlus->firstlet+wlu_init2offst(wlus->casedep,*uwrd)));
}

WRD_LU_REC *wlu_lurec4wrd(WRD_LUSTRCT *wlus,
                          WRD_LU_REC *wlst,
                          char *uwrd)
/* scan wlst for uwrd, return the first match found, NULL if not */
{
WRD_LU_REC *wp;
int ulen;

wp = wlst;
if (wlus->casedep == WLU_CASEWILD)
  ulen = strlen(uwrd);
while (wp != NULL)
  {
  switch (wlus->casedep)
    {
    case WLU_CASEIND:
      if (strcasecmp(uwrd,wp->wrd) == 0)    /* found it */
        return(wp);
      break;
    case WLU_CASEWILD:         /* wild - check if unique */
      if (strncasecmp(uwrd,wp->wrd,(size_t) ulen) == 0)
        if (wlu_cntwldmats(wp,uwrd,ulen) > 1)
          return(NULL);
        else
          return(wp);
      break;        
    case WLU_CASEDEP:
    default:
      if (strcmp(uwrd,wp->wrd) == 0)      /* found it */
        return(wp);
      break;
    }
  wp = wp->nextiwrd;
  }
return(NULL);
}

WRD_LU_REC *wlu_wrd2lurec(WRD_LUSTRCT *wlus,
                          char *uwrd)
/* return the pointer to the lu_rec for uwrd, NULL if non-existent */
{
return(wlu_lurec4wrd(wlus,wlu_lulst4wrd(wlus,uwrd),uwrd));
}

void *wlu_chkwrdptr(WRD_LUSTRCT *wlus,
                    char *uwrd)
/* check for uwrd in wlus word look up structure.  if not found return NULL
value else return the set parameter.  word comparisons are based on case
  dependency setting */
{
WRD_LU_REC *wrpt;

if ((wrpt = wlu_wrd2lurec(wlus,uwrd)) != NULL)
  return(wrpt->retval);
else
  return(wlus->failret);  /* not found */
}

int wlu_chkwrd(WRD_LUSTRCT *wlus,
               char *uwrd)
/* check for uwrd in wlus word look up structure.  if not found return failret
value else return the set parameter.  word comparisons are based on case
  dependency setting */
{
WRD_LU_REC *wrpt;

if ((wrpt = wlu_wrd2lurec(wlus,uwrd)) != NULL)
  return((int)((long) wrpt->retval));
else
  return((int)((long) wlus->failret));  /* not found */
}

int wlu_chkint(WRD_LUSTRCT *wlus,
               int uvl)
/* check for uvl in wlus word look up structure. */
{
char nbuf[33];

sprintf(&nbuf[0],"%d",uvl);
return(wlu_chkwrd(wlus,&nbuf[0]));
}

void wlu_maktoklu(WLU_CHRLUTBL utbl,
                  char *luchrs)
  /* allocate and create a character look up table (7 bit ascii only),
 with chars in luchrs */
{
int cp;

for (cp = 0; cp <= MAXASCIIVAL; cp++)
  utbl[cp] = (index(luchrs,(char) cp) != NULL);
}

void wlu_makcmptoklu(WLU_CHRLUTBL utbl,
                     char *luchrs)
/* create a character look up table (7 bit ascii only),
for chars NOT in luchrs */
{
int cp;

wlu_maktoklu(utbl,luchrs);
for (cp = 0; cp <= MAXASCIIVAL; cp++)
  utbl[cp] = !utbl[cp];
}

int wlu_gettokensep(FILE *fl,
                    char *ubuf,
                    int blen,
                    WLU_CHRLUTBL stbl,
                    char *schr)
/* pull next token from fl, write into ubuf, upto blen.  sbtl is a 
table of ascii values of token separators. if *schr non-NULL, then set it
to the separating char */
{
int nc;
char *tp;
int bl;

tp = ubuf;
while ((nc = fgetc(fl)) != EOF)
  if (stbl[nc])     /* have a token separator */
    {
    if ((bl = (int) (tp - ubuf)) > 0)    /* does separate something */
      {
      if (bl < blen)
        *tp = '\0';         /* room for null, put it in */
      if (schr != NULL)
        *schr = (char) nc;
      return(bl);
      }
    }
  else                      /* wasn't separator, try to fit it in */
    if ((bl = (int)(tp - ubuf)) >= blen)  /* full, return anyway */
      {
      if (schr != NULL)
        *schr = '\0';
      return(bl);
      }
    else
      *tp++ = (char) nc;
if ((bl = (int)(tp - ubuf)) < blen)
  *tp = '\0';
if (schr != NULL)
  *schr = '\0';
return(bl);
}

int wlu_gettoken(FILE *fl,
                 char *ubuf,
                 int blen,
                 WLU_CHRLUTBL stbl)
/* pull next token from fl, write into ubuf, upto blen.  sbtl is a 
table of ascii values of token separators */
{
return(wlu_gettokensep(fl,ubuf,blen,stbl,NULL));
}

int wlu_sgettokensep(char *src,
                     int *sp,
                     char *ubuf,
                     int blen,
                     WLU_CHRLUTBL stbl,
                     char *schr)
/* pull next token from src at *sp, write into ubuf, upto blen.  *sp is
incremented.  sbtl is a 
table of ascii values of token separators. if *schr non-NULL, then set it
to the separating char */
{
char nc;
char *tp;
int bl;

tp = ubuf;
while ((nc = *(src+*sp)) != '\0')
  {
  (*sp)++;
  if (stbl[(int) nc])     /* have a token separator */
    {
    if ((bl = (int)(tp - ubuf)) > 0)    /* does separate something */
      {
      if (bl < blen)
        *tp = '\0';         /* room for null, put it in */
      if (schr != NULL)
        *schr = nc;
      return(bl);
      }
    }
  else                      /* wasn't separator, try to fit it in */
    if ((bl = (int)(tp - ubuf)) >= blen)  /* full, return anyway */
      {
      if (schr != NULL)
        *schr = '\0';
      return(bl);
      }
    else
      *tp++ = nc;
  }
if ((bl = (int)(tp - ubuf)) < blen)
  *tp = '\0';
if (schr != NULL)
  *schr = '\0';
return(bl);
}

int wlu_sgettoken(char *src,
                  int *sp,
                  char *ubuf,
                  int blen,
                  WLU_CHRLUTBL stbl)
/* pull next token from src at *sp (incremented), write into ubuf, 
  upto blen.  sbtl is a 
  table of ascii values of token separators */
{
return(wlu_sgettokensep(src,sp,ubuf,blen,stbl,NULL));
}

int wlu_newstrng(WRD_LUSTRCT *wlus,
                 char *swrd,
                 char *newlin)
/* insert newlin in place of existing (if any) text line for swrd.  Return
1 if OK, 0 if swrd is not in wlus */
{
WRD_LU_REC *wrpt;

if ((wrpt = wlu_wrd2lurec(wlus,swrd)) == NULL)
  return(0);
else
  {
  if (wrpt->helpline != NULL)
    memfree(wrpt->helpline);
  wrpt->helpline = bas_strdup(newlin);
  return(1);
  }
}

int wlu_newstr4int(WRD_LUSTRCT *wlus,
                   int sint,
                   char *newlin)
/* insert newlin in place of existing (if any) text line for sint.  Return
1 if OK, 0 if swrd is not in wlus */
{
char is[33];

sprintf(&is[0],"%d",sint);
return(wlu_newstrng(wlus,&is[0],newlin));
}

int wlu_delwrd(WRD_LUSTRCT *wlus,
               char *swrd)
/* delete the first entry related to swrd from *wlus, return 1 if found &
deleted */
{
WRD_LU_REC *wrpt;
WRD_LU_REC *prvp;
WRD_LU_REC *lpt;

if ((wrpt = wlu_wrd2lurec(wlus,swrd)) == NULL)
  return(0);
else
  {
  lpt = wlu_lulst4wrd(wlus,swrd);
  if (wrpt == lpt)              /* is first element */
    {
    if ((*(wlus->firstlet+wlu_init2offst(wlus->casedep,*swrd)) =
                                              wrpt->nextiwrd) == NULL)
      *(wlus->lastlet+wlu_init2offst(wlus->casedep,*swrd)) = NULL;
    wlu_killwelt(wrpt);
    return(1);
    }
  else
    {
    prvp = lpt;
    while (prvp != NULL)
      if (prvp->nextiwrd == wrpt)         /* have prv */
        {
        if ((prvp->nextiwrd = wrpt->nextiwrd) == NULL)  /* is last elt */
          *(wlus->lastlet+wlu_init2offst(wlus->casedep,*swrd)) = prvp;
        wlu_killwelt(wrpt);
        return(1);
        }
      else
        prvp = prvp->nextiwrd;
    return(0);            /* got to end of list without finding previous!! */
    }  
  }
}

int wlu_delint(WRD_LUSTRCT *wlus,
               int di)
/* delete the first entry related to int di from *wlus, return 1 if found &
deleted */
{
char nbuf[33];

sprintf(&nbuf[0],"%d",di);
return(wlu_delwrd(wlus,&nbuf[0]));
}

WRD_LUSTRCT *wlu_getmnthwrds()
  /* return a wlu table that will convert 3 letter month names into
corresponding integers in case-independent manner */
{
WRD_LUSTRCT *ws;
int mpt;

ws = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"Month word lu struct");
wlu_initlustrct(ws,WLU_CASEIND,0);
for (mpt = 1; mpt <= 12; mpt++)
  wlu_addwrd(ws,bas_int2strmnth(mpt,1),1,NULL);
return(ws);
}

void wlu_dsplymenu(FILE *fl,
                   WRD_LUSTRCT *menu)
/* display the existing contents of *menu */
{
int ap;  /* alpha pointer */
WRD_LU_REC *wlpt;
int ninits;

ap = 0;
ninits = wlu_no_inits(menu->casedep);
while (ap < ninits)
  {
  wlpt = *(menu->firstlet+ap);
  while (wlpt != NULL)
    {
    if (wlpt->helpline != NULL)
      fprintf(fl,"    %s - %s\n",wlpt->wrd,wlpt->helpline);
    wlpt = wlpt->nextiwrd;
    }
  ap++;
  }
}

char *wlu_intval2str(WRD_LUSTRCT *ws,
                     int qval)
/* scan through ws, returning pointer to
the first corresponding string found,
NULL if none */
{
int ap;  /* alpha pointer */
WRD_LU_REC *wlpt;
int ninits;

ap = 0;
ninits = wlu_no_inits(ws->casedep);
while (ap < ninits)
  {
  wlpt = *(ws->firstlet+ap);
  while (wlpt != NULL)
    {
    if (wlpt->retval == (void *) ((long) qval))
      return(wlpt->wrd);
    else
      wlpt = wlpt->nextiwrd;
    }
  ap++;
  }
/* fell off end, didn't find value */
return(NULL);
}

void wlu_dbgwlustrct(FILE *fl,
                     WRD_LUSTRCT *menu,
                     char *msg)
/* tell fl about contents of *menu for diagnostic purposes */
{
int ap;  /* alpha pointer */
WRD_LU_REC *wlpt;
int ninits;

ap = 0;
fprintf(fl,"%s:@%lx:",(msg==NULL?"":msg),(long) menu);
fflush(fl);
ninits = wlu_no_inits(menu->casedep);
fprintf(fl,"casedep=%d,fail=%lx,initials=%d\n",menu->casedep,
            (long) menu->failret,ninits);
fflush(fl);
while (ap < ninits)
  {
  wlpt = *(menu->firstlet+ap);
  fprintf(fl,"ap=%d:@%lx:",ap,(long) wlpt);
  fflush(fl);
  while (wlpt != NULL)
    {
    fprintf(fl,"@%lx:",(long) wlpt);
    fflush(fl);
    fprintf(fl,"%lx\"%s\"",(long) wlpt->wrd,(wlpt->wrd==NULL?"<NULL>":wlpt->wrd));
    fflush(fl);
    if (wlpt->helpline != NULL)
      fputs("...",fl);
    if ((wlpt = wlpt->nextiwrd) != NULL)
      fputc(',',fl);
    fflush(fl);
    }
  ap++;
  fputc('\n',fl);
  fflush(fl);
  }
}

int wlu_maxwrdlen(WRD_LUSTRCT *wlus)
  /* return the length of the longest string in wlus */
{
int mxlen;
int ap;
int ninits;
WRD_LU_REC *wp;
int slen;

mxlen = 0;
ap = 0;
ninits = wlu_no_inits(wlus->casedep);
while (ap < ninits)
  {
  wp = *(wlus->firstlet+ap++);
  while (wp != NULL)
    {
    if ((slen = strlen(wp->wrd)) > mxlen)
      mxlen = slen;
    wp = wp->nextiwrd;
    }
  }
return(mxlen);
}

int wlu_initwrdscan(WRD_LUSTRCT *wls,
                    int *ap,
                    WRD_LU_REC **rp)
/* set initial values of ap & *rp, return true if it worked */
{
if (wls != NULL)
  {
  *ap = 0;
  while (*ap < wlu_no_inits(wls->casedep))
    if ((*rp = *(wls->firstlet + *ap)) != NULL)
      return(1);
    else
      (*ap)++;  /* try next init */
/* fell right off end, */
  return(0);
  }
else
  return(0);
}

WRD_LU_REC *wlu_nxtwrd(WRD_LUSTRCT *wls,
                       int *ap,
                       WRD_LU_REC *rprv)
/* given a pointer rprv, return the next element, if necessary, *ap is
incremented and the element taken from there.  NULL for no more */
{
if (rprv != NULL)
  rprv = rprv->nextiwrd;
if (rprv != NULL)
  return(rprv);
else
  {
  while (*ap < wlu_no_inits(wls->casedep))
    {
    (*ap)++;
    if ((rprv = *(wls->firstlet + *ap)) != NULL)
      return(rprv);
    }
  return(NULL);
  }
}
