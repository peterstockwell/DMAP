/* fsm_ops.h: type definitions for general FSM operations */

typedef struct FS_reselt  /* a result element - part of l-list of actions */
  {
  void *action;           /* what to return when this result is attained */
  struct FS_sttelt *nxtstate;  /* next stt state - actually *FS_STTELT */
  int locatd;             /* true when this result is found */
  struct FS_reselt *nxtreselt;  /* forward link */
  struct FS_reselt *prvreselt;  /* back link (needed?) */
  }
FS_RESELT;

typedef struct FS_sttelt  /* element in stt table */
  {
  int eltno;                 /* number of this element */
  FS_RESELT *reslts;      /* pointer to malloced array of result elements */
  char *s_cis;                 /* Current Interesting Sequence */
  int achk;                     /* for access checking */
  struct FS_sttelt *nxtsttelt;  /* forward list of these for housekeeping */
#ifdef U_TRACK
  int axcnt;                    /* count the frequency of access */
#endif
  }
FS_STTELT;

typedef struct FS_datprelt       /* for linked list of data elements */
  {
  char *dstrng;                  /* data string */
  int ldstr;                     /* length of dstrng */
  void *stxt;                    /* accompanying text/etc. */
  int mallced;                   /* set if stxt can be free()ed */
  int plocat;                    /* 0 if not located */
  struct FS_datprelt *nxtpelt;   /* forward link */ 
  struct FS_datprelt *prvpelt;   /* back link */
  }
FS_DATPRELT;

typedef enum FS_oninvld   /* options for invalid input characters */
  {
  FS_inv_reset = 0,       /* reset fsm to init state */
  FS_inv_ressay,          /* ditto, report error to stderr */
  FS_inv_ignor,           /* ignore */
  FS_inv_ignsay           /* ignore, but mention to stderr */
  }
FS_ONINVLD;

typedef struct FS_fsmstrct       /* fundamental info */
  {
  WRD_LUSTRCT *cis;
  WRD_LUSTRCT *instrngs;
  FS_STTELT *firststt;             /* state transition table */
  FS_STTELT *laststt;              /* max current state */
  FS_STTELT *currstt;              /* set to current state */
  int nstats;                      /* No states for input */
  int dosubs;                      /* 1=complete substrings, 0=don't */
  void *null_reslt;                /* what to return for noscore */
  char *cisbuf;                    /* buffer malloc()ed for cis */
  int strmax;                      /* longest string in input set */
  int sttcnt;                      /* current no of stt elements */
  FS_DATPRELT *prlst;              /* start of paired element list */
  FS_DATPRELT *lstpr;              /* last pair element */
  FS_ONINVLD invact;               /* define behavior for invalid in chars */
  int dbgbld;                      /* trace behaviour, defaults to 0=don't */
  }
FS_FSMSTRCT;

typedef enum FS_bldstatus          /* values when fsm is under construction */
  {
  FS_bldnoint = 0,      /* no interest */
  FS_bldmatcis,         /* matches existing stt cis */
  FS_bldmatinpt,        /* matches an input string */
  FS_bldintrst          /* matches head of input string */
  }
FS_BLDSTATUS;

/* function headers */

WRD_LU_REC *fs_shed2lurec(FS_FSMSTRCT *fss,
                          char *cstrng,
                          size_t lcstr);
/* return a pointer to the lu_rec for which cstrng matches the head.
NULL if non-existent */

void fs_initfsmstrct(FS_FSMSTRCT *fs,
                     int nstates,
                     int dsubs,
                     FS_ONINVLD inv);
/* assign initial (NULL largely) values to fs components.  dsubs determines if
fsm will completely handle substrings - The method of Smith (1988) CABIOS, 4,
459 does not correctly handle substrings.  It is made optional here in order
to allow some uses to work properly.  dsubs=1 - build fsm for substrings,
dsubs=0 - build fsm as in Smith, 1988 */

FS_FSMSTRCT *fs_initnewfsm(int nstates,
                           int dsubs,
                           FS_ONINVLD inv);
/* malloc() storage for a new fsm & assign initial (NULL largely) values to
fs components, return address. dsubs is described in fs_initfsmstrct() */

FS_RESELT *fs_appndreselt(FS_RESELT **rlst,
                          void *rsp);
/* append a new value to *rlst, return address of new element */

void fs_clrreslst(FS_RESELT **rlst);
  /* Chain through *rlst, freeing all malloc()ed memory */

FS_STTELT *fs_appndsttelt(FS_STTELT **lstrt,
                          int nstats);
/* append a new element to *lstrt, and making nstats result slots for it.
  return address of new element */

void fs_initreselt(FS_RESELT *rep);
  /* put NULL values into rep */

void fs_initfsmbld(FS_FSMSTRCT *fsm,
                   int casedep);
/* init data structures prior to building fsm */

void fs_clrprlst(FS_DATPRELT **rlst);
  /* Chain through *rlst, freeing all malloc()ed memory */

int fs_cntdatprs(FS_DATPRELT *pp);
  /* return the number of datalist items in *pp */

void fs_adddatprs(FS_FSMSTRCT *fss,
                  char *sstrng,
                  void *lval);
/* put a data pair (sstrng associated with lval) into fss. */

void fs_addstrprs(FS_FSMSTRCT *fss,
                  char *sstrng,
                  char *stxt);
/* put a data pair (sstrng associated with stxt) into fss, stxt is
strdup()ed locally. */

int fs_lstcstr(FILE *fl,
               char *str);
/* display characters of str to fl, converting non-printable chars to standard
C representation, return no chars actually written */

void fs_showpval(FILE *fl,
                 void *p);
/* let p be an arbitrary pointer, print it */

int fs_inprcnt(FS_DATPRELT *plst,
               FS_DATPRELT *pp);
/* count pair elements in plst till pp, 0 for none */

void fs_lststrpr(FILE *fl,
                 FS_DATPRELT *pp,
                 int cnt,
                 void (*slstfn)(FILE *bfl,
                                void *p),
                 void (*newlnfn)(FILE *afl));
/* recursively list *pp */

void fs_itlststrpr(FILE *fl,
                   FS_DATPRELT *plst,
                   int cstrt,
                   void (*slstfn)(FILE *bfl,
                                  void *p),
                   void (*newlnfn)(FILE *afl));
/* interatively list *pp */

int fs_dprs2wlu(FS_FSMSTRCT *fss,
                int casedep);
  /* malloc()s storage for instrng wlu table, inits and loads data pairs
into it.  Return count of included values */

void fs_lstwluprs(FILE *fl,
                  WRD_LUSTRCT *wlus);
/* list the contents of wlus at fl */

int fs_chkinpstr(FS_FSMSTRCT *fss,
                 FS_STTELT *sttp,
                 int ap,
                 char *acis);
/* check if acis matches any input strings & if so, note them as outputs.
Effectively implements Rule 3, return count of matching strings */

WRD_LU_REC *shed2lureclin(FS_FSMSTRCT *fss,
                          char *cstrng,
                          size_t lcstr);
/* return a pointer to the lu_rec for which cstrng matches the head.
NULL if non-existent */

int fs_bldstt4prs(FS_FSMSTRCT *fss,
                  char (* int2chrfn)(int ci),
                  int (*chkinpstrfn)(FS_FSMSTRCT *xfss,
                                     FS_STTELT *xsttp,
                                     int xap,
                                     char *xacis),
                  WRD_LU_REC *(*shed2lurecfn)(FS_FSMSTRCT *zfss,
                                              char *zcstrng,
                                              size_t zlcstr));
/* go thru logic for building stt for loaded data, return no of states
required.  (*int2chrfn)() is used to convert an alpha pointer to a char */

void fs_purgecis(FS_FSMSTRCT *fss);
  /* remove cis storage from fss */

void fs_purgelustrct(FS_FSMSTRCT *fss);
  /* remove lookup structures and cis buffer from fss */

void fs_purgefsm(FS_FSMSTRCT *fss);
  /* remove unnecessary storage from fss */

void fs_clrfsm(FS_FSMSTRCT *fss);
  /* clear storage from fss, leave basic fss structure, but denude it of 
any malloc()ed stuff */

void fs_killfsm(FS_FSMSTRCT **fss);
  /* completely free() *fss storage, base fss structure included */

int fs_chkfsm(FS_FSMSTRCT *fss,
              int erprt,
              char (*int2chr)(int ic));
/* check fss to ensure:
  1: all states are valid & exist
  2: all data pairs are located
  3: anything else I think of
return 1 if OK.  erprt controls error reporting: 0 means say nothing, 1
mention each error */

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
                                          size_t zlcstr));
/* combine a series of fsm building steps, return no of states if checks OK,
else -1 */

void fs_initrun(FS_FSMSTRCT *fss);
/* initialise the fsm for running */

FS_RESELT *fs_procchr(FS_FSMSTRCT *fss,
                      char nc,
                      int (*chr2int)(char c));
/* feed nc into fsm and return any result list.  if nc is invalid, behavior
is determined by fss->invact. The returned value is the heading to a linked 
list of results */

int fs_rmdupstts(FS_FSMSTRCT *ffst);
  /* scan stts, looking for duplicates which have no actions: combine into
first occurrence, and delete the duplicate.  Return number of dups so
removed */

void fs_lststtelt(FILE *fl,
                  FS_FSMSTRCT *fss,
                  FS_STTELT *stp);
/* tell fl about stp */

void fs_lststttbl(FILE *fl,
                  FS_FSMSTRCT *fss);
/* list info in fss */

void fs_tellusage(FILE *fl,
                  FS_FSMSTRCT *fss);
/* tell fl about the usage of each state of fss */
