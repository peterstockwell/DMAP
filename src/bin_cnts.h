/* bin_cnts.h: defs for count-binning stuff */

#define BC_DEF_BINLEN 1000

typedef struct BC_chr_bin       /* bin for chromosomal position counts */
  {
  int spos;                      /* start position */
  int binlen;                    /* length this bin */
  int cpgcnt;                    /* count CpG's this bin */
  int bcntfwd;                   /* fwd (meth) counts this bin */
  int bcntrev;                   /* rev (non-meth) counts this bin */
  struct BC_chr_bin *nxtbin;    /* forward/backward links */
  struct BC_chr_bin *prvbin;
  }
BC_CHR_BIN;

typedef enum BC_outmode         /* output mode */
  {
  BC_out_none = 0,           /* really to let debugging info be shown by self */
  BC_out_list,               /* list of counts+/- */
  BC_out_listnz,             /* list nonzero bins */
  BC_out_dmthrstrep,         /* output for differentially methylated regions */
  BC_out_rstrrepbins,        /* scan for CCGG and generate bins based on positions & size */
  BC_out_rstrreplst,         /* ditto, but list bins to stdout */
  BC_out_rstrrepmiss,        /* show reads which miss reduced rep bins */
  BC_out_rrmisstots,         /* just show totals for reads which miss RR bins */
  BC_out_sqmonkdat           /* generate a SeqMonk dat file for RR genome */
  }
BC_OUTMODE;

typedef struct BC_regn_elt   /* for linked list of regions */
  {
  int rstart;                /* start this region */
  int rstop;                 /* end this region (inclusive) */
  struct BC_regn_elt *nxtregn; /* forward link */
  struct BC_regn_elt *prvregn; /* back link */
  }
BC_REGN_ELT;

