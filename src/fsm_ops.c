/* fsm_ops.c: common routines for general FSM operations */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "fsm_ops.h"

void fs_initfsmstrct(FS_FSMSTRCT *fs,
                     int nstates,
                     int dsubs,
                     FS_ONINVLD inv)
/* assign initial (NULL largely) values to fs components.  dsubs determines if
fsm will completely handle substrings - The method of Smith (1988) CABIOS, 4,
459 does not correctly handle substrings.  It is made optional here in order
to allow some uses to work properly.  dsubs=1 - build fsm for substrings,
dsubs=0 - build fsm as in Smith, 1988 */
{
fs->instrngs = fs->cis = NULL;
fs->firststt = fs->laststt = fs->currstt = NULL;
fs->nstats = nstates;
fs->null_reslt = NULL;
fs->cisbuf = NULL;
fs->strmax = 0;
fs->sttcnt = 0;
fs->prlst = fs->lstpr = NULL;
fs->invact = inv;
fs->dosubs = dsubs;
}

FS_FSMSTRCT *fs_initnewfsm(int nstates,
                           int dsubs,
                           FS_ONINVLD inv)
/* malloc() storage for a new fsm & assign initial (NULL largely) values to
fs components, return address. dsubs is described in fs_initfsmstrct() */
{
FS_FSMSTRCT *fss;

fss = (FS_FSMSTRCT *) getmemory(sizeof(FS_FSMSTRCT),"FSM structure");
fs_initfsmstrct(fss,nstates,dsubs,inv);
return(fss);
}

FS_RESELT *fs_appndreselt(FS_RESELT **rlst,
                          void *rsp)
/* append a new value to *rlst, return address of new element */
{
FS_RESELT *prvp, *endp;

if (rlst != NULL)
  {
  prvp = endp = *rlst;
  while (endp != NULL)
    {
    prvp = endp;
    endp = endp->nxtreselt;
    }
  endp = (FS_RESELT *) getmemory(sizeof(FS_RESELT),"result element");
  endp->nxtreselt = NULL;
  endp->locatd = 0;
  endp->action = rsp;
  endp->nxtstate = NULL;
  if (*rlst == NULL)
    {
    *rlst = endp;
    endp->prvreselt = NULL;
    }
  else
    {
    prvp->nxtreselt = endp;
    endp->prvreselt = prvp;
    }
  return(endp);
  }
else
  return(NULL);
}

void fs_clrreslst(FS_RESELT **rlst)
  /* Chain through *rlst, freeing all malloc()ed memory */
{
FS_RESELT *nxt;

while (*rlst != NULL)
  {
  nxt = (*rlst)->nxtreselt;
  memfree(*rlst);
  *rlst = nxt;
  }
}

FS_STTELT *fs_appndsttelt(FS_STTELT **lstrt,
                          int nstats)
/* append a new element to *lstrt, and making nstats result slots for it.
  return address of new element */
{
FS_STTELT *prvp,*endp;

if (lstrt != NULL)
  {
  prvp = endp = *lstrt;
  while (endp != NULL)
    {
    prvp = endp;
    endp = endp->nxtsttelt;
    }
  endp = (FS_STTELT *) getmemory(sizeof(FS_STTELT),"STT element");
  endp->nxtsttelt = NULL;
  endp->s_cis = NULL;
#ifdef U_TRACK
  endp->axcnt = 0;
#endif
  endp->reslts = (FS_RESELT *) getmemory(nstats*sizeof(FS_RESELT),
                                           "result array for SST");
  if (*lstrt == NULL)
    *lstrt = endp;
  else
    prvp->nxtsttelt = endp;
  return(endp);
  }
else
  return(NULL);
}

void fs_initreselt(FS_RESELT *rep)
  /* put NULL values into rep */
{
rep->action = NULL;
rep->nxtstate = NULL;
rep->locatd = 0;
rep->nxtreselt = rep->prvreselt = NULL;
}

FS_STTELT *fs_newsttelt(FS_FSMSTRCT *fsp,
                        char *ccis)
  /* get a new stt element and clear it.  Return the address of the new
element.  Append it to the SST list */
{
FS_STTELT *newelt;
int sp;

newelt = fsp->laststt = fs_appndsttelt(&fsp->firststt,fsp->nstats);
newelt->eltno = fsp->sttcnt++;
newelt->achk = 0;
/* fprintf(stdout,"stt:%d=%s\n",newelt->eltno,(ccis?ccis:"-")); */
if (ccis != NULL)
  {
  newelt->s_cis = bas_strdup(ccis);
  wlu_addwrdptr(fsp->cis,ccis,(void *) newelt,NULL);  /* add to cis lu list */
  }
else
  newelt->s_cis = NULL;
for (sp = 0; sp < fsp->nstats; sp++)
  fs_initreselt(newelt->reslts + sp);
return(newelt);
}                        

int fs_sttonlst(FS_STTELT *plst,
                FS_STTELT *pp)
/* return 1 if pp is on list plst */
{
FS_STTELT *pt;

if (pp != NULL)
  {
  pt = plst;
  while (pt != NULL)
    if (pt == pp)
      return(1);
    else
      pt = pt->nxtsttelt;
  }
/* fell off end, return 0 */
return(0);
}

void fs_initfsmbld(FS_FSMSTRCT *fsm,
                   int casedep)
/* init data structures prior to building fsm */
{
fsm->cis = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"cis word lu struct");
wlu_initlustrct(fsm->cis,casedep,0);
fsm->null_reslt = NULL;
fsm->currstt = NULL;
fsm->strmax = (wlu_maxwrdlen(fsm->instrngs));
fsm->cisbuf = getmemory((fsm->strmax+2),"new cis buf");
}

WRD_LU_REC *fs_shed2lurec(FS_FSMSTRCT *fss,
                          char *cstrng,
                          size_t lcstr)
/* return a pointer to the lu_rec for which cstrng matches the head.
NULL if non-existent */
{
WRD_LU_REC *wrpt;

wrpt = wlu_lulst4wrd(fss->instrngs,cstrng);
while (wrpt != NULL)
  switch (fss->instrngs->casedep)
    {
    case WLU_CASEIND:
    case WLU_CASEWILD:
      if (strncasecmp(cstrng,wrpt->wrd,lcstr) == 0)    /* found it */
        return(wrpt->retval);
      else
        wrpt = wrpt->nextiwrd;
      break;
    case WLU_CASEDEP:
    default:
      if (strncmp(cstrng,wrpt->wrd,lcstr) == 0)      /* found it */
        return(wrpt->retval);
      else
        wrpt = wrpt->nextiwrd;
      break;
    }
return(NULL);  /* not found */
}

FS_DATPRELT *fs_appnddatpr(FS_DATPRELT **dlst,
                           char *sstr,
                           void *sstxt,
                           int mset)
/* append a new value to *dlst, return address of new element.
sstr values are strdup()ed locally, mset notes if the sstxt quantities
were malloc()ed and can thus be free()ed */
{
FS_DATPRELT *prvp, *endp;

if (dlst != NULL)
  {
  prvp = endp = *dlst;
  while (endp != NULL)
    {
    prvp = endp;
    endp = endp->nxtpelt;
    }
  endp = (FS_DATPRELT *) getmemory(sizeof(FS_DATPRELT),"data pair element");
  endp->nxtpelt = NULL;
  endp->dstrng = bas_strdup(sstr);
  endp->ldstr = strlen(sstr);
  endp->stxt = sstxt;
  endp->plocat = 0;
  endp->mallced = mset;
  if (*dlst == NULL)
    {
    *dlst = endp;
    endp->prvpelt = NULL;
    }
  else
    {
    prvp->nxtpelt = endp;
    endp->prvpelt = prvp;
    }
  return(endp);
  }
else
  return(NULL);
}

void fs_clrprlst(FS_DATPRELT **rlst)
  /* Chain through *rlst, freeing all malloc()ed memory */
{
FS_DATPRELT *nxt;

while (*rlst != NULL)
  {
  nxt = (*rlst)->nxtpelt;
  if ((*rlst)->mallced && ((*rlst)->stxt != NULL))
    memfree((*rlst)->stxt);
  if ((*rlst)->dstrng != NULL)
    memfree((*rlst)->dstrng);
  memfree(*rlst);
  *rlst = nxt;
  }
}

int fs_cntdatprs(FS_DATPRELT *pp)
  /* return the number of datalist items in *pp */
{
if (pp == NULL)
  return(0);
else
  return(fs_cntdatprs(pp->nxtpelt)+1);
}

int fs_inprcnt(FS_DATPRELT *plst,
               FS_DATPRELT *pp)
/* count pair elements in plst till pp, 0 for none */
{
FS_DATPRELT *pt;
int cnt;

if (pp == NULL)
  return(0);
else
  {
  cnt = 1;
  pt = plst;
  while (pt != NULL)
    if (pt == pp)
      return(cnt);
    else
      {
      cnt++;
      pt = pt->nxtpelt;
      }
/* fell off end, return 0 */
  return(0);
  }
}

void fs_adddatprs(FS_FSMSTRCT *fss,
                  char *sstrng,
                  void *lval)
/* put a data pair (sstrng associated with lval) into fss. */
{
fss->lstpr = fs_appnddatpr(&fss->prlst,sstrng,lval,0);
}

void fs_addstrprs(FS_FSMSTRCT *fss,
                  char *sstrng,
                  char *stxt)
/* put a data pair (sstrng associated with stxt) into fss, stxt is
strdup()ed locally. */
{
fss->lstpr = fs_appnddatpr(&fss->prlst,sstrng,bas_strdup(stxt),1);
}

int fs_saycchr(FILE *fl,
               char c)
/* write character c to fl in "c" style, return no of chars actually used */
{
if ((c < ' ') || (c > '~'))
  {
  fprintf(fl,"\\%03o",(int) c);
  return(4);
  }
else
  {
  fputc(c,fl);
  return(1);
  }
}

int fs_lstcstr(FILE *fl,
               char *str)
/* display characters of str to fl, converting non-printable chars to standard
C representation, return no chars actually written */
{
char *sp;
int ccnt;

ccnt = 0;
if (str != NULL)
  {
  sp = str;
  while (*sp != '\0')
    ccnt += fs_saycchr(fl,*sp++);
  }
return(ccnt);
}

void fs_lstdatpr(FILE *fl,
                 FS_DATPRELT *pp,
                 int cnt)
/* recursively list *pp */
{
if (pp != NULL)
  {
  fprintf(fl,"%d: \"",cnt);
  (void) fs_lstcstr(fl,pp->dstrng);
  fprintf(fl,"\" - %lx\n",(long) pp->stxt);
  fs_lstdatpr(fl,pp->nxtpelt,++cnt);
  }
}

void fs_showcstr(FILE *fl,
                 void *p)
/* let p be pointer to a char string, print it c-style */
{
if (p)
  {
  fputc('"',fl);
  (void) fs_lstcstr(fl,p);
  fputc('"',fl);
  }
else
  fputs("<NULL>",fl);
}

void fs_showpval(FILE *fl,
                 void *p)
/* let p be an arbitrary pointer, print it */
{
fprintf(fl,"%lx",(long) p);
}

void fs_lststrpr(FILE *fl,
                 FS_DATPRELT *pp,
                 int cnt,
                 void (*slstfn)(FILE *bfl,
                                void *p),
                 void (*newlnfn)(FILE *afl))
/* recursively list *pp */
{
if (pp != NULL)
  {
  fprintf(fl,"%d:\"",cnt);
  (void) fs_lstcstr(fl,pp->dstrng);
  fputs("\"-",fl);
  (* slstfn)(fl,(char *) pp->stxt);
  (*newlnfn)(fl);
  fs_lststrpr(fl,pp->nxtpelt,++cnt,slstfn,newlnfn);
  }
}

void fs_itlststrpr(FILE *fl,
                   FS_DATPRELT *plst,
                   int cstrt,
                   void (*slstfn)(FILE *bfl,
                                  void *p),
                   void (*newlnfn)(FILE *afl))
/* interatively list *pp */
{
FS_DATPRELT *pp;
int cnt;

pp = plst;
cnt = cstrt;
while (pp != NULL)
  {
  fprintf(fl,"%d:\"",cnt);
  (void) fs_lstcstr(fl,pp->dstrng);
  fputs("\"-",fl);
  (* slstfn)(fl,(char *) pp->stxt);
  (*newlnfn)(fl);
  pp = pp->nxtpelt;
  cnt++;
  }
}

int fs_dprs2wlu(FS_FSMSTRCT *fss,
                int casedep)
  /* malloc()s storage for instrng wlu table, inits and loads data pairs
into it.  Return count of included values */
{
FS_DATPRELT *pp;
int cnt;

fss->instrngs = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),
                                            "in word lu struct");
wlu_initlustrct(fss->instrngs,casedep,0);
pp = fss->prlst;
cnt = 0;
while (pp != NULL)
  {
  if (pp->mallced)
    wlu_addwrdptr(fss->instrngs,pp->dstrng,(void *) pp,pp->stxt);
  else
    wlu_addwrdptr(fss->instrngs,pp->dstrng,(void *) pp,NULL);
  pp = pp->nxtpelt;
  cnt++;
  }
return(cnt);
}

void fs_lstprwlureclst(FILE *fl,
                       WRD_LU_REC *wp)
/* recursively list elements of wp */
{
if (wp != NULL)
  {
  fprintf(fl,"\"%s\" - \"%s\" %lx\n",wp->wrd,wp->helpline,(long) wp->retval);
  fs_lstprwlureclst(fl,wp->nextiwrd);
  }
}

void fs_lstwluprs(FILE *fl,
                  WRD_LUSTRCT *wlus)
/* list the contents of wlus at fl */
{
int ap;
int ninits;

ap = 0;
ninits = wlu_no_inits(wlus->casedep);
while (ap < ninits)
  fs_lstprwlureclst(fl,*(wlus->firstlet+ap++));
}

FS_STTELT *fs_nxtincmplt(FS_FSMSTRCT *fss,
                         int *ii)
/* return the address of the next incomplete element in stt.
ii, if non-NULL will be set to the alphabetic value for the first non-completed position of this stt element.  NULL if none */
{
FS_STTELT *sp;
int ap;

if (fss->currstt == NULL)
  fss->currstt = fss->firststt;
sp = fss->currstt;
while (sp != NULL)
  {
  ap = 0;
  while (ap < fss->nstats)
    if ((sp->reslts+ap)->nxtstate == NULL)
      {
      if (ii != NULL)
        *ii = ap;
      fss->currstt = sp;       /* note this state, to reduce scan time */
      return(sp);
      }
    else
      ap++;
  sp = sp->nxtsttelt;
  }
return(fss->currstt = NULL);
}

size_t fs_cislen(FS_STTELT *stt)
  /* return the length of the current cis */
{
if ((stt == NULL) || (stt->s_cis == NULL))
  return(0);
else
  return(strlen(stt->s_cis));
}

char *fs_bldstat2str(FS_BLDSTATUS bs)
  /* a string for bs */
{
switch (bs)
  {
  case FS_bldnoint:
    return("FS_bldnoint");
    break;
  case FS_bldmatcis:
    return("FS_bldmatcis");
    break;
  case FS_bldmatinpt:
    return("FS_bldmatinpt");
    break;
  case FS_bldintrst:
    return("FS_bldintrst");
    break;
  default:
    return("FS_bld???");
    break;
  }
}

int fs_chkinpstr(FS_FSMSTRCT *fss,
                 FS_STTELT *sttp,
                 int ap,
                 char *acis)
/* check if acis matches any input strings & if so, note them as outputs.
Effectively implements Rule 3, return count of matching strings */
{
FS_RESELT *pt;
WRD_LU_REC *wrp;
FS_DATPRELT *dp;
int mcnt;

mcnt = 0;
wrp = wlu_wrd2lurec(fss->instrngs,acis);
while (wrp != NULL)
  { /* it did, note this as a result */
  mcnt++;
  dp = (FS_DATPRELT *) wrp->retval;
  pt = sttp->reslts + ap;
  (void) fs_appndreselt(&pt,(void *) dp);
  dp->plocat = 1;
  wrp = wlu_lurec4wrd(fss->instrngs,wrp->nextiwrd,acis);
  }
return(mcnt);
}

int fs_chksubstrs(FS_FSMSTRCT *fss,
                  FS_STTELT *sttp,
                  int ap,
                  char *acis)
/* check if beheaded substrings of acis match any output states & if so,
note them as outputs. return count of any matches */
{
char *sp;
int mcnt;

mcnt = 0;
if (*acis != '\0')
  {
  sp = acis + 1;
  while (*sp != '\0')
    {
    mcnt += fs_chkinpstr(fss,sttp,ap,sp);
    sp++;
    }
  }
return(mcnt);
}

void fs_trynewcis(FS_FSMSTRCT *fss,
                  FS_STTELT *sttp,
                  int ap,
                  char *nwcis,
                  int clen,
                  int (*chkinpstrfn)(FS_FSMSTRCT *xfss,
                                     FS_STTELT *xsttp,
                                     int xap,
                                     char *xacis),
                  WRD_LU_REC *(*shed2lurecfn)(FS_FSMSTRCT *zfss,
                                              char *zcstrng,
                                              size_t zlcstr))
/* perform Rules 3,4,5 on acis, return information on the status */
{
FS_STTELT *cp;
char *acis;            /* "start" of nwcis */
int chkintcis;

acis = nwcis;
chkintcis = 1;
(sttp->reslts+ap)->nxtstate = fss->firststt;
while (*acis != '\0')   /* while we still have some new cis to look for */
  {  /* rule 3: look for a match in input strings */
  (void) (*chkinpstrfn)(fss,sttp,ap,acis);
/* Rule 4: Does this match an existing cis?? */
  if (chkintcis)
    if ((cp = (FS_STTELT *) wlu_chkwrdptr(fss->cis,acis)) != NULL)
       /* yep, note the newstate as the next state for this result */
      {
      (sttp->reslts+ap)->nxtstate = cp;
      chkintcis = 0;
      }
/* Rule 5: is this of interest?? */
    else
      {
      if ((*shed2lurecfn)(fss,acis,clen) != NULL)
        {  /* is interesting - create new state for it */
        fss->laststt = fs_newsttelt(fss,acis);      /* Rule 5 */
        (sttp->reslts+ap)->nxtstate = fss->laststt;
        chkintcis = 0;
        }
      }
  acis++;
  clen--;
  }
}

int fs_trycisnosub(FS_FSMSTRCT *fss,
                   FS_STTELT *sttp,
                   int ap,
                   char *nwcis,
                   int clen,
                   int (*chkinpstrfn)(FS_FSMSTRCT *xfss,
                                      FS_STTELT *xsttp,
                                      int xap,
                                      char *xacis),
                   WRD_LU_REC *(*shed2lurecfn)(FS_FSMSTRCT *zfss,
                                               char *zcstrng,
                                               size_t zlcstr))
/* perform Rules 3,4,5 on acis, return information on the status, method of
Smith (1988), Cabios, 4, 459 */
{
FS_DATPRELT *dp;
FS_STTELT *cp;
char *acis;            /* "start" of nwcis */
int inlutmat;
int addcnt;
int cismat;

acis = nwcis;
addcnt = 0;
(sttp->reslts+ap)->nxtstate = fss->firststt;
/* is this cis an output? if so note it - Rule 3 */
inlutmat = (*chkinpstrfn)(fss,sttp,ap,acis) > 0;
while (clen > 0)
  {    /* is this cis already present as a cis? - Rule 4 */
  if ((cismat = ((cp = (FS_STTELT *) wlu_chkwrdptr(fss->cis,acis)) != NULL)))
    { /* yes, note new state as next state for this result */
    (sttp->reslts+ap)->nxtstate = cp;
    return(addcnt);        /* and stop looking */
    }
      /* does acis match start of any site - Rule 5 */
  if ((dp = (FS_DATPRELT *) (*shed2lurecfn)(fss,acis,clen))
         != NULL)
    {     /* is interesting - create new state for it, if not located */
    if (dp->plocat == 0)
      {
      fss->laststt = fs_newsttelt(fss,acis);      /* Rule 5 */
      (sttp->reslts+ap)->nxtstate = fss->laststt;
      return(++addcnt);
      }
    if (inlutmat || cismat)   /* bale out if matched input or cis */
      return(addcnt);
    }
  acis++;
  clen--;
  }  
return(addcnt);
}

int fs_bldstt4prs(FS_FSMSTRCT *fss,
                  char (* int2chrfn)(int ci),
                  int (*chkinpstrfn)(FS_FSMSTRCT *xfss,
                                     FS_STTELT *xsttp,
                                     int xap,
                                     char *xacis),
                  WRD_LU_REC *(*shed2lurecfn)(FS_FSMSTRCT *zfss,
                                              char *zcstrng,
                                              size_t zlcstr))
/* go thru logic for building stt for loaded data, return no of states
required.  (*int2chrfn)() is used to convert an alpha pointer to a char */
{
int ap;
FS_STTELT *nxtinc;      /* next incomplete element */
size_t clen;
char *sp;

fss->laststt = fs_newsttelt(fss,NULL);      /* Rule 1 */
while ((nxtinc = fs_nxtincmplt(fss,&ap)) != NULL)
  {
  clen = fs_cislen(nxtinc);  /* make the new cis: Rule 2 */
  if (clen > fss->strmax)
    {
    fprintf(stderr,"Cis length (%d) exceeds buffer\n",(int) clen);
    exit(1);
    }
  *fss->cisbuf = '\0';
  if (nxtinc->s_cis != NULL)
    strncat(fss->cisbuf,nxtinc->s_cis,clen);
  sp = fss->cisbuf + clen;
  *sp++ = (*int2chrfn)(ap);
  *sp = '\0';
               /* is this interesting?? */
  if (fss->dosubs)
    fs_trynewcis(fss,nxtinc,ap,fss->cisbuf,(int)(clen+1),chkinpstrfn,
                   shed2lurecfn);
  else
    (void) fs_trycisnosub(fss,nxtinc,ap,fss->cisbuf,(int)(clen+1),chkinpstrfn,
                            shed2lurecfn);
  }
return(fss->sttcnt);
}

void fs_purgecis(FS_FSMSTRCT *fss)
  /* remove cis storage from fss */
{
FS_STTELT *sp;

sp = fss->firststt;
while (sp != NULL)
  {
  if (sp->s_cis != NULL)
    {
    memfree(sp->s_cis);
    sp->s_cis = NULL;
    }
  sp = sp->nxtsttelt;
  }
}

void fs_purgelustrct(FS_FSMSTRCT *fss)
  /* remove lookup structures and cis buffer from fss */
{
wlu_clrlustrct(fss->cis);
fss->cis = NULL;
wlu_clrlustrct(fss->instrngs);
fss->instrngs = NULL;
if (fss->cisbuf != NULL)   /* make sure OK for possible repeat purges */
  memfree(fss->cisbuf);
fss->cisbuf = NULL;
}

void fs_purgefsm(FS_FSMSTRCT *fss)
  /* remove unnecessary storage from fss */
{
fs_purgecis(fss);
fs_purgelustrct(fss);
}

void fs_clrfsm(FS_FSMSTRCT *fss)
  /* clear storage from fss, leave basic fss structure, but denude it of 
any malloc()ed stuff */
{
FS_STTELT *sp;
FS_STTELT *nxt;
FS_RESELT *rp;
int ap;

fs_purgefsm(fss);
sp = fss->firststt;
while (sp != NULL)      /* for each stt table entry */
  {   /*  interate thru each reslt, clearing any results */
  for (ap = 0; ap < fss->nstats; ap++)
    {
    rp = (sp->reslts+ap)->nxtreselt;
    fs_clrreslst(&rp);
    }
  memfree(sp->reslts);    /* de_malloc basic reslts array */
  nxt = sp->nxtsttelt;   /* check nxt element */
  memfree(sp);           /* kill this element */
  sp = nxt;
  }
fss->firststt = fss->laststt = NULL;
fs_clrprlst(&fss->prlst);
fss->lstpr = NULL;
}

void fs_killfsm(FS_FSMSTRCT **fss)
  /* completely free() *fss storage, base fss structure included */
{
if (*fss != NULL)
  {
  fs_clrfsm(*fss);
  memfree(*fss);
  *fss = NULL;
  }
}

int fs_chkfsm(FS_FSMSTRCT *fss,
              int erprt,
              char (*int2chr)(int ic))
/* check fss to ensure:
  1: all data pairs are located
  2: all states are valid & exist
  3: all stt entries will be accessed, somehow (except state 0)
  4: anything else I think of
return 1 if OK.  erprt controls error reporting: 0 means say nothing, nonzero
mention each error */
{
FS_DATPRELT *dp;
int icnt;
int pcnt;
int stat;
FS_STTELT *sp;
int ap;
FS_STTELT *np;

/* Test 1: all input states are located */
stat = 1;
pcnt = icnt = 0;
dp = fss->prlst;
if (fss->dosubs)
  {
  while (dp != NULL)
    {
    pcnt++;
    if (!dp->plocat)
      {
      icnt++;
      if (erprt)
        {
        fprintf(stderr,"string %d \"",fs_inprcnt(fss->prlst,dp));
        (void) fs_lstcstr(stderr,dp->dstrng);
        fputs("\" not located\r\n",stderr);
        }
      }
    dp = dp->nxtpelt;
    }
  if ((pcnt > 0) && (icnt > 0))
    stat = 0;
  }
/* tests 2 & 3: all stt pointers are valid and all stt states are accessible */
sp = fss->firststt;
while (sp != NULL)
  {
  for (ap = 0; ap < fss->nstats; ap++)
    {
    np = (sp->reslts+ap)->nxtstate;
    if (!fs_sttonlst(fss->firststt,np))
      {
      stat = 0;
      if (erprt)
        fprintf(stderr,"STT state %d has invalid nextstate (%lx) for '%c'\n",
                  sp->eltno,(long) (sp->reslts+ap)->nxtstate,(*int2chr)(ap));
      }
    else   /* test 3 - mark next state as accessible */
      np->achk = 1;
    }
  sp = sp->nxtsttelt;
  }
sp = fss->firststt;
while (sp != NULL)
  {
  if ((!sp->achk) && (sp->eltno > 0))   /* state 0 may not be addressed */
    {
    stat = 0;
    if (erprt)
      fprintf(stderr,"STT state %d inaccessible\n",sp->eltno);
    }    
  sp = sp->nxtsttelt;
  }
return(stat);
}

int fs_bldfsm(FS_FSMSTRCT *fss,
              int casedep,
              int errpt,
              int purge,
              char (*int2chr)(int ic),
              int (*chkinpstrfn)(FS_FSMSTRCT *xfss,
                                 FS_STTELT *xsttp,
                                 int xap,
                                 char *xacis),
              WRD_LU_REC *(*shed2lurecfn)(FS_FSMSTRCT *zfss,
                                          char *zcstrng,
                                          size_t zlcstr))
/* combine a series of fsm building steps, return no of states if checks OK,
else -1.  errpt controls error reporting: 1 = mention each error, 0 = don't,
  -1 = check, but ignore result  */
{
int nstate;

if (fs_dprs2wlu(fss,casedep) > 0)
  {
  fs_initfsmbld(fss,casedep);
  if ((nstate = fs_bldstt4prs(fss,int2chr,chkinpstrfn,shed2lurecfn)) > 0)
    {
    if (purge)
      fs_purgefsm(fss);
    if ((fs_chkfsm(fss,errpt,int2chr)) || (errpt == -1))
      return(nstate);
    }
  }
return(-1);
}

void fs_initrun(FS_FSMSTRCT *fss)
/* initialise the fsm for running */
{
fss->currstt = fss->firststt;
}

FS_RESELT *fs_procchr(FS_FSMSTRCT *fss,
                      char nc,
                      int (*chr2int)(char c))
/* feed nc into fsm and return any result list.  if nc is invalid, behavior
is determined by fss->invact. The returned value is the heading to a linked 
list of results */
{
int ap;
FS_STTELT *sp;
FS_RESELT *rtp;

if (((ap = (*chr2int)(nc)) < 0) || (ap >= fss->nstats))
  {
  switch (fss->invact)
    {
    case FS_inv_ressay:
      fputs("Invalid char '",stderr);
      (void) fs_saycchr(stderr,nc);
      fputs("' - FSM is reset\n",stderr);
    case FS_inv_reset:
      fs_initrun(fss);
      break;
    case FS_inv_ignsay:
      fputs("Invalid char '",stderr);
      (void) fs_saycchr(stderr,nc);
      fputs("' - FSM not reset\n",stderr);
    case FS_inv_ignor:
    default:
      break;
    }
  return(NULL);
  }
else
  {
  rtp = (fss->currstt->reslts+ap)->nxtreselt;
  sp = (fss->currstt->reslts+ap)->nxtstate;
#ifdef U_TRACK
  sp->axcnt++;
#endif
  fss->currstt = sp;
  return(rtp);
  }
}

int fs_sttsmatch(FS_FSMSTRCT *ffst,
                 FS_STTELT *stt1,
                 FS_STTELT *stt2)
/* return 1 if stt1 and stt2 match: same nextstate for all positions and
NULL action */
{
int ep;
FS_RESELT *rs1;
FS_RESELT *rs2;

if ((stt1 == NULL) || (stt2 == NULL))  /* one/both doesn't exist, ret 0 */
  return(0);
else
  {
  for (ep = 0; ep < ffst->nstats; ep++)
    {
    rs1 = stt1->reslts + ep;
    rs2 = stt2->reslts + ep;
    if (rs1->nxtstate != rs2->nxtstate)
      return(0);
    else
      if (!((rs1->nxtreselt == NULL) && (rs2->nxtreselt == NULL)))
        return(0);
    }
/* fell off end: match - return 1 */
  return(1);
  }
}

int fs_modnxtsttelts(FS_FSMSTRCT *ffst,
                     FS_STTELT *oldstt,
                     FS_STTELT *newstt)
/* traverse all of stt in ffst, replacing every reference to oldstt with
one to newstt.  Return the number so modified */
{
FS_STTELT *stp;
int ep;
int mdcnt;
FS_RESELT *rp;

stp = ffst->firststt;
mdcnt = 0;
while (stp != NULL)
  {
  for (ep = 0; ep < ffst->nstats; ep++)
    {
    rp = stp->reslts + ep;
    if (rp->nxtstate == oldstt)
      {
      rp->nxtstate = newstt;
      mdcnt++;
      }
    }
  stp = stp->nxtsttelt;
  }
return(mdcnt);
}

int fs_rmdupstts(FS_FSMSTRCT *ffst)
  /* scan stts, looking for duplicates which have no actions: combine into
first occurrence, and delete the duplicate.  Return number of dups so
removed */
{
FS_STTELT *sp1;
FS_STTELT *sp2;
int kcnt;
FS_STTELT *prvs;
int ep;
FS_RESELT *rp;

kcnt = 0;
sp1 = ffst->firststt;
while (sp1 != NULL)
  {
  prvs = sp1;
  sp2 = sp1->nxtsttelt;
  while (sp2 != NULL)
    {
    if (fs_sttsmatch(ffst,sp1,sp2))
      {
      kcnt++;
      if (fs_modnxtsttelts(ffst,sp2,sp1) > 0)
        {
        for (ep = 0; ep < ffst->nstats; ep++)
          {
          rp = (sp2->reslts+ep)->nxtreselt;
          fs_clrreslst(&rp);
          }
        memfree(sp2->reslts);
        prvs->nxtsttelt = sp2->nxtsttelt;
        memfree(sp2);
        }
      sp2 = prvs->nxtsttelt;
      }
    prvs = sp2;
    if (sp2 != NULL)
      sp2 = sp2->nxtsttelt;
    }
  sp1 = sp1->nxtsttelt;
  }
return(kcnt++);
}

void fs_lststtelt(FILE *fl,
                  FS_FSMSTRCT *fss,
                  FS_STTELT *stp)
/* tell fl about stp */
{
int ap;
FS_RESELT *nxt;
FS_STTELT *snx;
FS_RESELT *dp;

fprintf(fl,"%d:%s",stp->eltno,(stp->s_cis?stp->s_cis:"-"));
for (ap = 0; ap < fss->nstats; ap++)
  {
  nxt = stp->reslts+ap;
  if (nxt == NULL)
    fputs(" ?/?",fl);
  else
    if ((snx = (FS_STTELT *) nxt->nxtstate) == NULL)
      fprintf(fl," Null/*");
    else
      {
      fprintf(fl," %d",snx->eltno);
      if (nxt->nxtreselt == NULL)
        fputs("/-",fl);
      else
        {
        dp = (FS_RESELT *) nxt->nxtreselt;
        while (dp != NULL)
          {
          fprintf(fl,"/%d",fs_inprcnt(fss->prlst,dp->action));
          dp = dp->nxtreselt;
          }
        }
      }
  }
fputc('\n',fl);
fflush(fl);
}

void fs_lststttbl(FILE *fl,
                  FS_FSMSTRCT *fss)
/* list info in fss */
{
FS_STTELT *stp;

stp = fss->firststt;
while (stp != NULL)
  {
  fs_lststtelt(fl,fss,stp);
  stp = stp->nxtsttelt;
  }
}

void fs_tellusage(FILE *fl,
                  FS_FSMSTRCT *fss)
/* tell fl about the usage of each state of fss */
{
#ifdef U_TRACK

FS_STTELT *sp;
int cttl;
int scnt;
int ocnt;

ocnt = scnt = cttl = 0;
sp = fss->firststt;
while (sp != NULL)
  {
  fprintf(fl,"%d\t%d\t%s\n",sp->eltno,sp->axcnt,(sp->s_cis?sp->s_cis:""));
  scnt++;
  cttl += sp->axcnt;
  if (sp->axcnt > 0)
    ocnt++;
  sp = sp->nxtsttelt;
  }
fprintf(fl,"\n%d States (%d used) with total %d counts\n",scnt,ocnt,cttl);

#else

fputs("STT usage is not enabled: compile with -DU_TRACK\n",fl);

#endif
}
