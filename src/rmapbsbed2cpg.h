/* headers for rmapbsbed2cpg: generate list of CpG positions from rmapbs BED file */

typedef enum RBC_flowcellversn
  {
  RBC_fcv_2 = 0,     /* GAII */
  RBC_fcv_3_3,       /* HiSeq V3 chemistry */
  RBC_fcv_3_2        /* HiSeq V2 chemistry */
  }
RBC_FLOWCELLVERSN;

typedef enum RBC_readflform
  {
  RBC_readflfm_fasta = 0, /* header is >s<laneno>_<tileno>_<xpixel>_<ypixel> */
  RBC_readflfm_fastq      /* header is @...:<laneno>:<tileno>:<xpixel>:<ypixe> where ... may contain ':' */
  }
RBC_READFLFORM;

#define MAXTILENO 120

typedef enum RBC_Chrno
  {
  Chr_unk = 0,
  Chr1,
  Chr2,
  Chr3,
  Chr4,
  Chr5,
  Chr6,
  Chr7,
  Chr8,
  Chr9,
  Chr10,
  Chr11,
  Chr12,
  Chr13,
  Chr14,
  Chr15,
  Chr16,
  Chr17,
  Chr18,
  Chr19,
  Chr20,
  Chr21,
  Chr22,
  ChrX,
  ChrY
  }
RBC_CHRNO;

typedef struct RBC_read_elt   /* element containing read info */
  {
  int xpix;                   /* x,y pixel positions on tile */
  int ypix;
  int chrno;                  /* chromosome No. */
  int pstart;                 /* start/stop positions from .BED */
  int pstop;
  int readfwd;                /* 1=> forward read, 0=> reverse */
  char *readseq;              /* this read seq */
  struct RBC_read_elt *nxtrdelt;
  struct RBC_read_elt *prvrdelt;
  }
RBC_READ_ELT;

typedef enum RBC_dbg     /* verbosity level */
  {
  RBC_dbg_none = 0,
  RBC_dbg_on,
  RBC_dbg_serious
  }
RBC_DBG;

typedef enum RBC_out_style
  {
  RBC_out_cpg = 0
  }
RBC_OUT_STYLE;

typedef enum RBC_src_style
  {
  RBC_src_rmapbs = 0,
  RBC_src_bsmap102,     /* v1.02 or earlier? bsmap */
  RBC_src_bsmap12       /* v1.2 bsmap */
  }
RBC_SRC_STYLE;

typedef enum RBC_rrread_form  /* to enumerate rr read maps */
  {
  RBC_rrread_nomap = 0,
  RBC_rrread_5p,             /* maps at 5' end of read */
  RBC_rrread_3p,             /* ditto   3' */
  RBC_rrread_internal,       /* maps in middle of rrbs fragment */
  RBC_rrread_onjoin,         /* maps over rrbs fragment join */
  RBC_rrread_is5pto,         /* read is 5' to rrbs fragment */
  RBC_rrread_is3pto          /* read is 3' to rrbs fragment */
  }
RBC_RRREAD_FORM;
