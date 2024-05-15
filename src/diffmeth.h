/* diffmeth.h: defs for diffmeth count-binning stuff */

#define PROG_VERSION 1.83
/* -o option for user output file: May-2024 */
/* #define PROG_VERSION 1.82 */
/* count CpG details for WGBS: Nov-2022 */
/* #define PROG_VERSION 1.81 */
/* correct WGBS count check: Nov-2022 */
/* #define PROG_VERSION 1.80 */
/* correct final CpG for fixed windows: Jun-2022 */
/* #define PROG_VERSION 1.79 */
/* improve meth proportion output: Feb-2022 */
/* #define PROG_VERSION 1.78 */
/* tidy some details for DMAP2: Aug-2021 */
/* #define PROG_VERSION 1.77 */
/* user-defined fragments: Jan/Feb-2018 */
/* #define PROG_VERSION 1.76 */
/* allow sam/bam FLAG values other than 0 & 16; Nov-2017 */
/* #define PROG_VERSION 1.75 */
/* manage bam files with '*' genome hits: Oct-2017 */
/* #define PROG_VERSION 1.74 */
/* use fstat in sqfl_filelength() to get seq buff length: Aug-2017 */
/* #define PROG_VERSION 1.73 */
/* correct anova R & S methylation figures with valid check: Aug-2017 */
/*  #define PROG_VERSION 1.72 */
/* minor changes for NO_ZLIB: Apr-2017 */
/* #define PROG_VERSION 1.71 */
/* multi treatment group ANOVA: Nov-2016 */
/* #define PROG_VERSION 1.70 */
/* allow for non-CpG methylation: Nov-2016 */
/* #define PROG_VERSION 1.63 *
/* bam_fns bzero() z_stream field; Dec-2015 */
/* #define PROG_VERSION 1.62 */
/* -n other restriction sites: Oct-2015 */
/* #define PROG_VERSION 1.61 */
/* close chromosome source files: May-2015 */
/* #define PROG_VERSION 1.60 */
/* bam file input, change in -z option: Apr-2015 */
/* #define PROG_VERSION 1.55 */
/* tidy valid output checking for -L CpG listing: Mar-2015 */
/* #define PROG_VERSION 1.54 */
/* check valid chromosomal positions for SAM input: Dec-2014 */
/* #define PROG_VERSION 1.53 */
/* correct fixed bin width issue for end of chromosomes: Dec-2014 */
/* #define PROG_VERSION 1.52 */
/* check that at least -g or -G used: Nov-2014 */
/* #define PROG_VERSION 1.51 */
/* include fold difference column for pairwise -R & -S samples: Sep-2014 */
/* #define PROG_VERSION 1.50 */
/* -G file to specify chromosomes&files: Jul-2014 */
/* #define PROG_VERSION 1.49 */
/* correct fragment/sample checking for stddev output: Apr-2014 */
/* #define PROG_VERSION 1.48 */
/* reinstall std dev (-D) output: Mar-2014 */
/* #define PROG_VERSION 1.47 */
/* add meth proportions to -A anova output: Feb-2014 */
/* #define PROG_VERSION 1.46 */
/* other group option for Fisher's exact: Feb-2014 */
/* #define PROG_VERSION 1.45 */
/* add -A Anova group methylation info value: Jan-2014 */
/* #define PROG_VERSION 1.44 */
/* tidy up some of F-ratio logic: make consistent with double vs float: Jan-2014 */
/* #define PROG_VERSION 1.43 */
/* options to compare control/treatment, lose old -S-s options: Nov-2013 */
/* #define PROG_VERSION 1.42 */
/* correct -f logic (removed by mistake): Oct-2013 */
/* #define PROG_VERSION 1.41 */
/* allow cpg count listing on one sample: Jul-2013 */
/* #define PROG_VERSION 1.40 */
/* correct Fisher's Exact criterion recursive test error: June 2013 */
/* #define PROG_VERSION 1.39 */
/* -N CpG correction for non-adjacent fragments: Mar-2013 */
/* #define PROG_VERSION 1.38 */
/* -N CpG count increment for chi, etc.: Feb-2013 */
/* #define PROG_VERSION 1.37 */
/* correct dm_cntsarevalid threshold comparison: Feb-2013 */
/* #define PROG_VERSION 1.36 */
/* correct error with binp in dm_notemethstatinbin: Feb-2013 */
/* #define PROG_VERSION 1.35 */
/* correct per CpG count 3'SAM problem: Jan-2013 */
/* #define PROG_VERSION 1.34 */
/* allow modified coverge criteria, density/2 leading CpGs: Nov-2012 */
/* #define PROG_VERSION 1.33 */
/* Allow fixed-width bins (Bock, et al, 2012 Molecular Cell, 47): Oct-2012 */
/* #define PROG_VERSION 1.32 */
/* manage 3' reads, leading CpG position: Oct-2012 */
/* #define PROG_VERSION 1.31 */
/* read .sam files for input: Sep-2012 */
/* #define PROG_VERSION 1.30 */
/* join adjacent RRBS fragments, implement bin listing: Sep-2012 */
/* #define PROG_VERSION 1.23 */
/* implement hypomethylation criterion: Aug-2012 */
/* #define PROG_VERSION 1.22 */
/* correct FE pairwise test failing to write CpG counts: Jul-2012 */
/* #define PROG_VERSION 1.21 */
/* -U option & per-CpG counts: May-2012 */
/* #define PROG_VERSION 1.20 */
/* allow individual CpG mapping: Apr-2012 */
/*#define PROG_VERSION 1.10 */
/* various incremental changes to output, +&-cnts/CpG: Mar-2012 */
/* #define PROG_VERSION 1.00 */
/* faster lookup, modified meth proportions: Dec-2011 */
/* first hack, ex bin_cnts: Nov-2011 */
/* #define PROG_VERSION 0.0 */

#define DM_LIST_FREQ_DEF 1000
#define DEF_SAMBUFLEN 150
#define SAM_CMT_CHR '@'

typedef struct dm_cnts4cpg      /* to contain counts for an individual CpG */
  {
  int cpgpos;
  int metcnt;
  int unmetcnt;
  }
DM_CNTS4CPG;

typedef struct dm_meth_cnts     /* linked list element for meth/unmeth counts */
  {
  int methcnt;
  int unmethcnt;
  int cntused;                    /* for listing relevant counts */
  int smplgroup;                  /* integer to specify which group (control/treat) 1-based */
  DM_CNTS4CPG *cpgcntarr;         /* malloced() array of cpg count elements, NULL if unwanted */
  struct dm_meth_cnts *nxtelt;
  struct dm_meth_cnts *prvelt;
  }
DM_METH_CNTS;

typedef struct dm_cpgpos_elt     /* for a linked list of CpG positions in bin */
  {
  int cpgpos;
  struct dm_cpgpos_elt *nxtcpgelt;
  struct dm_cpgpos_elt *prvcpgelt;
  }
DM_CPGPOS_ELT;

typedef struct dm_cpg_bin       /* bin for chromosomal position counts */
  {
  int spos;                      /* start position */
  int binlen;                    /* length this bin */
  int cpgcnt;                    /* count CpG's this bin */
  int maxccnt;                   /* cpgcnt if treated as 1, else 2xcpgcnt */
  DM_METH_CNTS *cntlst;          /* linked list of count elements */
  DM_CPGPOS_ELT *cpgposlist;    /* linked list of cpg positions */
  struct dm_cpg_bin *nxtbin;    /* forward/backward links */
  struct dm_cpg_bin *prvbin;
  }
DM_CPG_BIN;

typedef struct DM_bin_lstelt    /* to implement a way of looking up bins more rapidly */
  {
  int binsstrt;                 /* first position represented by this bin series */
  int binsstop;                 /* last position  ditto */
  DM_CPG_BIN *binptr;           /* first bin in for above start */
  struct DM_bin_lstelt *nxtbinelt;
  struct DM_bin_lstelt *prvbinelt;
  }
DM_BIN_LSTELT;

typedef struct DM_chr_info     /* information about & bins for a chromosome */
  {        /* particularly the individual bin lst and the bin list elements */
  DM_CPG_BIN *binlst;
  DM_BIN_LSTELT *beltlst;
  char *chromseq;              /* for if we need sam file input */
  int chrlen;
  char *chrflname;             /* file name for fa sequencs */
  char *chrid;                 /* id for this chr */
  }
DM_CHR_INFO;

typedef enum dm_outmode         /* output mode */
  {
  dm_out_none = 0,           /* really to let debugging info be shown by self */
  dm_out_list,               /* list of counts+/- */
  dm_out_listnz,             /* list nonzero bins */
  dm_out_pairwise,           /* pairwise output for differentially methylated regions */
  dm_out_chisq,              /* chisquare counts for as many individuals as possible */
  dm_out_fisherforce,        /* fisher's exact obligatory */
  dm_out_chiforce,           /* chisq obligatory */
  dm_out_bestpairchi,        /* try to run with best of FE or chi, giving min pairwise PR or chiPR */
  dm_out_allpairchi,         /* list all prs (for paired FE) */
  dm_out_anova,              /* do analysis of variance */
  dm_out_anova_gtr,          /* Anova: show greater meth group */
  dm_out_anova_more,         /* ANOVA, even more detail */
  dm_out_sdev,               /* std dev of proportions */
  dm_out_cpg_detail,         /* list counts on perCpG basis */
  dm_out_cpg_nzdetl          /* list nonZero counts on perCpG basis */
  }
DM_OUTMODE;

typedef enum DM_cntchk_type    /* type of check for bin count validity */
  {
  dm_chk_none = 0,
  dm_chk_cntcrit,              /* -t */
  dm_chk_cpghits,              /* -T */
  dm_chk_maxcpghits,           /* -u */
  dm_chk_minmeth,              /* -M */
  dm_chk_cpgcnt,               /* -C */
  dm_chk_indcnt,               /* -I */
  dm_chk_mincpgcrit,           /* -F */
  dm_chk_folddiff              /* -f */
  }
DM_CNTCHK_TYPE;

typedef enum DM_src_mode
  {
  DM_src_cpglist = 0,
  DM_src_samfile,
  DM_src_bamfile,
  DM_src_ctlist
  }
DM_SRC_MODE;

typedef struct dm_fnam_elt   /* for linked list of file names */
  {
  char *flnam;
  DM_SRC_MODE srcmode;
  int sgroup;                /* number of sample group (ctrl/treatment) */
  char *groupid;             /* cmdline identifier this sample */
  struct dm_fnam_elt *nxtfelt; /* forward link */
  }
DM_FNAM_ELT;

#define ANOVA_GROUPS 2

typedef struct dm_anova_dat  /* data storage for anova - to avoid remallocing */
  {
  int maxgroup;              /* how many groups we have; 1 based */
  double *ssx;               /* array of ss */
  double *sx;                /* array of sums */
  int *grpcnts;              /* array of counts in each group */
  double grandtotal;         /* total of all observations */
  }
DM_ANOVA_DAT;

typedef struct DM_sites
  {
  char *sitestr;
  int offset;
  struct DM_sites *nxtsite;
  struct DM_sites *prvsite;
  }
DM_SITES;

typedef enum DM_bin_style
  {
  DM_bin_restenz = 0,
  DM_bin_fixed,
  DM_bin_userlist
  }
DM_BIN_STYLE;

typedef struct DM_runpars
  {
  DM_OUTMODE omode;
  int maxchrno;       /* changeable limit for various species */
  int havexy;         /* for species which do/don't have XY */
  char *uchrno;       /* user-specified chromosome: NULL=>all */
  float folddiff;     /* fold difference for Fisher's exact */
  int cntthreshold;   /* must exceed this number of counts */
  float minhitspercpg; /* must exceed this number hits/CpG, 0.0 to ignore */
  float maxhitspercpg; /* must not exceed this No hits/CpG, 0.0 to ignore */
  int displycnts;     /* control listing of counts for each individual */
  int binclusterfreq; /* frequency with which bins are clustered */
  char *listdelmtr;   /* what delimits output fields */
  char *chrprefix;    /* write this prior to chromosome ID */
  int cntcpgsas1;     /* whether we want complementary CpGs to map to 5' C or not */
  float prthreshold;  /* pr threshold for comparisons: 1.0 to disable */
  int minsamples;     /* require a minimum of minsamples for a fragment */
  int mincpgs;        /* minimum CpGs for a fragment to be valid */
  float minmethprop;  /* minimum methylation for valid frag: restrict hypometh: 0.0 to ignore */
  int joinfrags;      /* 1 to join adjacent frags */
  int needseq;        /* 1 if we need to store the sequence (sam file input) */
  int sambuflen;      /* bufferlength for reads in sam files */
  DM_CPG_BIN *currbin;  /* current bin, to reduce lookups */
  int cachechrno;     /* as above */
  int samfst3p2prv;   /* set => count first 3' sam CpG to preceding frag */
  int samsrccnt;      /* count of sam source files */
  int binwidth;       /* how to construct bins 0=>MspI, else (Bock, et al, 2012 Molecular Cell, 47) */
  DM_BIN_STYLE binstyle;
  FILE *binfile;      /* for user file of bins (chr start end) */
  int binchkmask;     /* bit mask for bin validity checks to perform */
  int cpgdetails;     /* set if we need to note individual cpg info */
  int cpgcntscrit;    /* no of CpGs for count criteria: 0 for all */
  double (*binprfun)(struct DM_runpars *rpxx,
                     DM_CPG_BIN *binxx,
                     int *dfxx);  /* function to return probability & df */
  DM_CHR_INFO *chrbininfo;
  DM_ANOVA_DAT *anova_dat;
  int maxgrp;        /* actual number of input groups: NOTE 0-based */
/*  int srclno; */
  char *genomsqhdrstr; /* header string for genomic sequences automatic name generation */
  char *chrinfoflnam;  /* name of a file containing chromosome information -G option */
  FILE *chrinfofl;     /* file of chromosome information -G option */
  WRD_LUSTRCT *chridlu; /* for chromosome ID lookup */
  DM_SRC_MODE srcmode;  /* mode for position or ct list files */
  int compressed;       /* 1 if -Z, else 0 (-z) */
#ifndef NO_ZLIB
  BF_RUNDATA *bamrdp;   /* pointer to bam fns rundata */
#endif
  FILE *rstfile;        /* non-NULL, if we have opened a file of restriction sites */
  char *rstflnam;
  DM_SITES *sitelist;
  int allcs;            /* 1 if we want all Cs ie non-CpG & CpG methylation */
  char **groupidlist;   /* to contain a list of group ids */
  FILE *outfile;        /* user output file, defaults to stdout */
  }
DM_RUNPARS;

typedef enum SAM_bsmk_fld
  {
  SAM_fld_header = 0,
  SAM_fld_sense,
  SAM_fld_chromosome,
  SAM_fld_positn,
  SAM_fld_flag,
  SAM_fld_cigar,
  SAM_fld_rnext,
  SAM_fld_pnext,
  SAM_fld_tlen,
  SAM_fld_seq,
  SAM_fld_qual,
  SAM_fld_nmtag,
  SAM_fld_xxtag,
  SAM_fld_methcall,
  SAM_fld_xrtag,
  SAM_fld_xgtag
  }
SAM_BSMK_FLD;

typedef enum DM_read_sens
  {
  DM_rdsens_5p = 0,
  DM_rdsens_3p
  }
DM_READ_SENS;

#ifndef BAM_FLAG_CMP
#define BAM_FLAG_CMP 16
#endif
/* Should be defined in bam_fns.h, where it should have been
in the first place, have it here in case */
