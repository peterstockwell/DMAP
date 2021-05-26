/* sqmat_fns.h: headers for sequence matching functions */

typedef enum Bas_MatMode    /* match strategies for bases */
  {
  BAS_exact = 0,
  BAS_iub,
  BAS_ldna
  }
BAS_MATMODE;

/* routine headers */

char case_mirror(char nc,
                 char casepat);
/* return nc in same case as casepat */

char ssd_bascmplmnt(char bs,
                    BAS_MATMODE mmod);
/* return the complementary base to bs, depending on mmod */

char same_residue(char bas,
                  BAS_MATMODE dummy);
/* return bas */

void reverse_seq(char *sq,
                 int sqlen,
                 char (* cmpfun)(char xres,
                                 BAS_MATMODE xmod),
                 BAS_MATMODE mmod);
/* reverse sq character order, apply
cmpfun to each transfer */

void complmnt_seq(char *sq,
                  int sqlen,
                  BAS_MATMODE mmod);
/* reverse&complement the contents of sq, depending on mmod */

int ssd_bas2bit(char bas);
/* return a bit map representing bas with A=bit 1, C=2, G=3 & T=4,
 otherwise 0 */

int ssd_bits4basechr(char bas,
                     BAS_MATMODE mmod);
/* return a bit map representing bas with A=bit 1, C=2, G=3 & T=4 and
  bits set for redundant specs, as for mmod */

int ssd_basmatch(char b1,
                 char b2,
                 BAS_MATMODE mmod);

char *ssd_nxtbasmatch(char *seq,
                      char bas,
                      BAS_MATMODE mmod);
/* return a pointer to the next occurrence of bas in seq.  NULL if none */

char *ssd_prvbasmatch(char *seq,
                      char *spos,
                      char bas,
                      BAS_MATMODE mmod);
/* return a pointer to the previous occurrence of bas in seq, starting 
from spos.  NULL if none */

char *ssd_nxtstrmatch(char *seq,
                      char *pat,
                      BAS_MATMODE mmod);
/* return the position of the next match of pat in seq. NULL if not found */

char *ssd_prvstrmatch(char *seq,
                      int sqlen,
                      int startpos,
                      char *pat,
                      BAS_MATMODE mmod);
/* return the position of the next match of pat in seq. NULL if not found */

int ssd_nxtmatchpos(char *seq,
                    int sqlen,
                    int curpos,
                    char *pat,
                    BAS_MATMODE mmod,
                    int srchfwd);
/* return the sequence position of next match after curpos of pat in seq.
  0 if not found */

int ssd_bitcnt4bas(char ares,
                   BAS_MATMODE mmod);
/* count up the bits set for ares (i.e. no bases it represents) */
