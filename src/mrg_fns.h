typedef struct MRG_regn_elt   /* for linked list of regions */
  {
  int rstart;                /* start this region */
  int rstop;                 /* end this region (inclusive) */
  int cpgcnt;                  /* count of CpGs this frag */
  struct MRG_regn_elt *nxtregn; /* forward link */
  struct MRG_regn_elt *prvregn; /* back link */
  }
MRG_REGN_ELT;

MRG_REGN_ELT *mrg_appndrgnelt(MRG_REGN_ELT **lstrt,
                             int rgstart,
                             int rgstop,
                             int cpgs);
/* create and append a new element to *lstrt,
 Return address of new element */

void mrg_delregnelt(MRG_REGN_ELT *ep,
                    MRG_REGN_ELT **lstrt);
/* delete ep from list *lstrt */

void mrg_clrallregnelts(MRG_REGN_ELT **lstrt);
  /* iteratively delete all of lstrt */

int mrg_cntregnelts(MRG_REGN_ELT *clst);
  /* recursively count list elements */

int mrg_sumregnelts(MRG_REGN_ELT *rlst);
  /* recursively sum lengths of elements */

int mrg_sumcpgregnelts(MRG_REGN_ELT *rlst);
  /* recursively ... */

MRG_REGN_ELT *mrg_regnelt4pos(MRG_REGN_ELT *rlst,
			      int pos);
/* return the region element corresponding to pos, if any */
