/* wlu_fns.h: headers for fast word lookup routines (for command processing )
 */

#define WLU_CASEDEP 0
#define WLU_CASEIND 1
#define WLU_CASEWILD -1

typedef struct Wrd_lu_rec
  {    /* basic structure for looking up words quickly */
  char *wrd;         /* this word */
  void *retval;      /* value to be returned */
  struct Wrd_lu_rec *nextiwrd;  /* next word same initial */
  char *helpline;    /* line of help info: can be NULL */
  }
WRD_LU_REC;

typedef struct Wrd_lustruct /* info for building word lookup */
  {
  int casedep;          /* 0 => case depend, 1 case independ, -1 wild card */
  void *failret;         /* value returned on failure */
  WRD_LU_REC **firstlet; /* array of first occurence of init letters */
  WRD_LU_REC **lastlet;  /* array of last occurrence of first letters */
  }
WRD_LUSTRCT;

#define MAXASCIIVAL 127  /* max 7 bit ascii value */

typedef int WLU_CHRLUTBL[MAXASCIIVAL+1]; 

/* function headers */

int wlu_no_inits(int cased);
  /* return the number of initials for
  cased: 0 => case depend, 1 case independ, -1 wild card */

void wlu_initlustrctptr(WRD_LUSTRCT *wlus,
                        int cased,
                        void *failval);
/*  init wlus values prior to loading words
  cased: 0 => case depend, 1 case independ, -1 wild card.
This is for initing pointer based lists */

void wlu_initlustrct(WRD_LUSTRCT *wlus,
                     int cased,
                     int failval);
/*  init wlus values prior to loading words
  cased: 0 => case depend, 1 case independ, -1 wild card */

WRD_LUSTRCT *wlu_getlustrct(int cased,
                            int failval);
/* allocate a wlu strcture and initialise it
with case & failval */

int wlu_init2offst(int cased,
                   char initc);
/* return a numerical offset for initc, depending on
 cased:  0 => case depend, 1 case independ, -1 wild card */

void wlu_addwrdptr(WRD_LUSTRCT *wlus,
                   char *newwrd,
                   void *newval,
                   char *hlpline);
/* add newwrd to lookup structure.  No check is performed for multiple
  insertions */

void wlu_addwrd(WRD_LUSTRCT *wlus,
                char *newwrd,
                int newval,
                char *hlpline);
/* add newwrd to lookup structure.  No check is performed for multiple
  insertions */

void wlu_addintgrs(WRD_LUSTRCT *wlus,
                   int newint,
                   int newval,
                   char *hlpline);
/* add newint as a string to lookup structure, returning newval.
  No check is performed for multiple insertions */

void wlu_addintgr(WRD_LUSTRCT *wlus,
                  int newval,
                  char *hlpline);
/* add newval as a string to lookup structure, returning newval.
  No check is performed for multiple insertions */

int wlu_delwrd(WRD_LUSTRCT *wlus,
               char *swrd);
/* delete the first entry related to swrd from *wlus, return 1 if found &
deleted */

int wlu_delint(WRD_LUSTRCT *wlus,
               int di);
/* delete the first entry related to int di from *wlus, return 1 if found &
deleted */

int wlu_cntwldmats(WRD_LU_REC *wrpt,
                   char *uwrd,
                   int wlen);
/* return count of words which match this one upto wlen chars in linked
  list wrpt */

WRD_LU_REC *wlu_lulst4wrd(WRD_LUSTRCT *wlus,
                          char *uwrd);
/* return the start of the list containing this word */

WRD_LU_REC *wlu_lurec4wrd(WRD_LUSTRCT *wlus,
                          WRD_LU_REC *wlst,
                          char *uwrd);
/* scan wlst for uwrd, return the first match found, NULL if not */

WRD_LU_REC *wlu_wrd2lurec(WRD_LUSTRCT *wlus,
                          char *uwrd);
/* return the pointer to the lu_rec for uwrd, NULL if non-existent */

void *wlu_chkwrdptr(WRD_LUSTRCT *wlus,
                    char *uwrd);
/* check for uwrd in wlus word look up structure.  if not found return NULL
value else return the set parameter.  word comparisons are based on case
  dependency setting */

int wlu_chkwrd(WRD_LUSTRCT *wlus,
               char *uwrd);
/* check for uwrd in wlus word look up structure.  if not found return failret
value else return the set parameter.  word comparisons are based on case
  dependency setting */

int wlu_chkint(WRD_LUSTRCT *wlus,
               int uvl);
/* check for uvl in wlus word look up structure. */

void wlu_clrlustrct(WRD_LUSTRCT *wlu);
  /* relinquish malloced memory for wlu lists */

void wlu_clralllustrct(WRD_LUSTRCT *wlu);
  /* relinquish malloced memory for wlu lists and wlu itself */

void wlu_maktoklu(WLU_CHRLUTBL utbl,
                  char *luchrs);
  /* allocate and create a character look up table (7 bit ascii only),
 with chars in luchrs */

void wlu_makcmptoklu(WLU_CHRLUTBL utbl,
                     char *luchrs);
/* create a character look up table (7 bit ascii only),
for chars NOT in luchrs */

int wlu_gettokensep(FILE *fl,
                    char *ubuf,
                    int blen,
                    WLU_CHRLUTBL stbl,
                    char *schr);
/* pull next token from fl, write into ubuf, upto blen.  sbtl is a 
table of ascii values of token separators. if *schr non-NULL, then set it
to the separating char */

int wlu_gettoken(FILE *fl,
                 char *ubuf,
                 int blen,
                 WLU_CHRLUTBL stbl);
/* pull next token from fl, write into ubuf, upto blen.  sbtl is a 
table of ascii values of token separators */

int wlu_sgettokensep(char *src,
                     int *sp,
                     char *ubuf,
                     int blen,
                     WLU_CHRLUTBL stbl,
                     char *schr);
/* pull next token from src at *sp, write into ubuf, upto blen.  *sp is
incremented.  sbtl is a 
table of ascii values of token separators. if *schr non-NULL, then set it
to the separating char */

int wlu_sgettoken(char *src,
                  int *sp,
                  char *ubuf,
                  int blen,
                  WLU_CHRLUTBL stbl);
/* pull next token from src at *sp (incremented), write into ubuf, 
  upto blen.  sbtl is a 
  table of ascii values of token separators */

int wlu_newstrng(WRD_LUSTRCT *wlus,
                 char *swrd,
                 char *newlin);
/* insert newlin in place of existing (if any) text line for swrd.  Return
1 if OK, 0 if swrd is not in wlus */

int wlu_newstr4int(WRD_LUSTRCT *wlus,
                   int sint,
                   char *newlin);
/* insert newlin in place of existing (if any) text line for sint.  Return
1 if OK, 0 if swrd is not in wlus */

WRD_LUSTRCT *wlu_getmnthwrds();
  /* return a wlu table that will convert 3 letter month names into
corresponding integers in case-independent manner */

void wlu_dsplymenu(FILE *fl,
                   WRD_LUSTRCT *menu);
/* display the existing contents of *menu */

char *wlu_intval2str(WRD_LUSTRCT *ws,
                     int qval);
/* scan through ws, returning pointer to
the first corresponding string found,
NULL if none */

void wlu_dbgwlustrct(FILE *fl,
                     WRD_LUSTRCT *menu,
                     char *msg);
/* tell fl about contents of *menu for diagnostic purposes */

int wlu_maxwrdlen(WRD_LUSTRCT *wlus);
  /* return the length of the longest string in wlus */

int wlu_initwrdscan(WRD_LUSTRCT *wls,
                    int *ap,
                    WRD_LU_REC **rp);
/* set initial values of ap & *rp, return true if it worked */

WRD_LU_REC *wlu_nxtwrd(WRD_LUSTRCT *wls,
                       int *ap,
                       WRD_LU_REC *rprv);
/* given a pointer rprv, return the next element, if necessary, *ap is
incremented and the element taken from there.  NULL for no more */
