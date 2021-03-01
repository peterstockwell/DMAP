/* dbpars.h: header datatypes and declarations for C db parsing routines */

#define DBLUVERSION 5.36
/* add enamelist element to DB_FEATSTRCT type to handle multiple values: Aug-2020 */
/* #define DBLUVERSIONNO 5.35 */
/* increase GTFMAXTOKCNT for gencode gtf to 256: Jul-2020 */
/* #define DBLUVERSIONNO 5.34 */
/* attribute value table for gtf/gff2 values: Apr-2017 */
/* #define DBLUVERSIONNO 5.33 */
/* allow GRCh38-type SeqMonk ID lines: Oct-2016 */
/* #define DBLUVERSIONNO 5.32 */
/* allow SeqMonk features, TSS & CpG: May-2012 */
/* #define DBLUVERSIONNO 5.31 */
/* attampt to parse gft files: Feb-2012 */
/* #define DBLUVERSIONNO 5.3 */
/* attempt to parse gff3 files: Jan-2012 */
/* #define DBLUVERSIONNO 5.2 */
/* to accommodate changes in EMBL ID lines: Sep-2011 */
/* #define DBLUVERSIONNO 5.1 */
/* add lclass parameter to db_parsecore to retrieve specified line text: Dec-2008*/
/* #define DBLUVERSIONNO 5.0 */
/* revise structures for feature/qualifier specification */
/* #define DBLUVERSIONNO 4.04 */
/* correct segment matching in db_featsplicsame(): Jul-2005 */
/* #define DBLUVERSIONNO 4.03 */
/* longer ids for refseqs, revise db_scanidlin(): Jun-2003 */
/* #define DBLUVERSIONNO 4.02 */
/* revised genomic locus line */
/* #define DBLUVERSIONNO 4.01 */
/* try to manage mixed entry joins: Nov-2002 */
/* #define DBLUVERSIONNO 4.00 */
/* increased field for file number */
/* #define DBLUVERSIONNO 3.02 */
/* check for too many source files */
/* #define DBLUVERSIONNO 3.01 */

#define DBLU_MAXIDLEN 15
#define DBLU_MBLHDR 2
#define DBLU_GENBHDR 10
#define DBLU_MAXHDR 10
#define DBLU_MBLFLDSTRT 6
#define DBLU_GBFLDSTRT 13
#define DBLU_MAXSRCLEN 256
#define DBLU_FKWSTART 6
#define DBLU_FLCSTART 22
#define DBLU_GNBIDFLDSTRT 13
#define DBLU_MBLIDFLDSTRT 6
#define DBLU_GNBLENFLDSTRT 23
#define DBLU_GNBNATPFLDSTRT 37
#define DBLU_GNBDATFLDSTRT 63
#define DBLU_GNBTOPFLDSTRT 43
#define DBLU_GNBSTRNDFLDSTRT 34
#define DB_FEATINSET 6
#define DB_FEATLOCATN 22
#define DB_GBSQINSET 11
#define DB_MBLSQINSET 6
#define DB_MAXKEYLN 32
#define DB_SWSFKWCOL 6
#define DB_SWFTPOS1COL 15
#define DB_SWFTPOS2COL 22
#define DB_SWFTQUALCOL 35
#define DB_FNODIGITS 2
#define DBLU_NGNBNATPFLDSTRT 48
#define DBLU_NGNBDATFLDSTRT 69
#define DBLU_NGNBTOPFLDSTRT 56
#define DBLU_NGNBSTRNDFLDSTRT 45
#define DBLU_NGNBLENFLDSTRT 30
#define GFF3TOKENCNT 9
#define GTFTOKENCNT 10
#define GTFMAXTOKCNT 256     /* arbitrary - to manage lots of tabbed attribute values */
#define GTFMAXQUAL 32
#define MAX_SWFTTOKENS 64   /* arbitrary, really */

typedef enum DBLU_dbfmt
  {
  DBFMT_unknown = 0,
  DBFMT_genbank,
  DBFMT_embl,
  DBFMT_embl_old,
  DBFMT_swiss,
  DBFMT_gff3,
  DBFMT_gtf,
  DBFMT_sqmonk
  }
DBLU_DBFMT;

typedef enum DBLU_lineclass
  {
  DBLC_unk = 0,
  DBLC_ident,
  DBLC_access,
  DBLC_naident,
  DBLC_datestamp,
  DBLC_description,
  DBLC_keyword,
  DBLC_source,
  DBLC_taxonomy,
  DBLC_reference,
  DBLC_author,
  DBLC_title,
  DBLC_journal,
  DBLC_remark,
  DBLC_feature,
  DBLC_basecount,
  DBLC_sequence,
  DBLC_ignore,
  DBLC_seqversn,
  DBLC_refnumbr,
  DBLC_refpositn,
  DBLC_xref,
  DBLC_dbxref,
  DBLC_terminator
  }
DBLU_LINECLASS;

typedef enum DBLU_ftkw
  {
  FTKW_unknown = 0,
  FTKW_ft_error,
  FTKW_allele,
  FTKW_attenuator,
  FTKW_CAAT_signal,
  FTKW_CDS,
  FTKW_cellular,
  FTKW_conflict,
  FTKW_C_region,
  FTKW_D_loop,
  FTKW_D_region,
  FTKW_enhancer,
  FTKW_exon,
  FTKW_GC_signal,
  FTKW_FT_gene,
  FTKW_iDNA,
  FTKW_insertion_seq,
  FTKW_intron,
  FTKW_J_region,
  FTKW_LTR,
  FTKW_mat_peptide,
  FTKW_misc_binding,
  FTKW_misc_difference,
  FTKW_misc_feature,
  FTKW_misc_recomb,
  FTKW_misc_RNA,
  FTKW_misc_signal,
  FTKW_misc_structure,
  FTKW_modified_base,
  FTKW_mRNA,
  FTKW_mutation,
  FTKW_N_region,
  FTKW_old_sequence,
  FTKW_polyA_signal,
  FTKW_polyA_site,
  FTKW_precursor_RNA,
  FTKW_prim_transcript,
  FTKW_primer,
  FTKW_primer_bind,
  FTKW_promoter,
  FTKW_protein_bind,
  FTKW_provirus,
  FTKW_RBS,
  FTKW_rep_origin,
  FTKW_repeat_region,
  FTKW_repeat_unit,
  FTKW_rRNA,
  FTKW_satellite,
  FTKW_scRNA,
  FTKW_sig_peptide,
  FTKW_snRNA,
  FTKW_stem_loop,
  FTKW_S_region,
  FTKW_STS,
  FTKW_TATA_signal,
  FTKW_terminator_sq,
  FTKW_transit_peptide,
  FTKW_transposon,
  FTKW_tRNA,
  FTKW_unsure,
  FTKW_variation,
  FTKW_virion,
  FTKW_V_region,
  FTKW_hyphen,
  FTKW_m10_signal,
  FTKW_m35_signal,
  FTKW_clip3p,
  FTKW_UTR3p,
  FTKW_clip5p,
  FTKW_UTR5p,
  FTKW_startcodon,
  FTKW_CNS,           /* conserved non-coding sequence */
  FTKW_intergenic,
  FTKW_source_ftkw,
  FTKW_CpG_is,
  FTKW_TSS,
  FTKW_SWP_act_site,
  FTKW_SWP_binding,
  FTKW_SWP_carbohyd,
  FTKW_SWP_ca_bind,
  FTKW_SWP_chain,
  FTKW_SWP_conflict,
  FTKW_SWP_disulfid,
  FTKW_SWP_dna_bind,
  FTKW_SWP_domain,
  FTKW_SWP_helix,
  FTKW_SWP_init_met,
  FTKW_SWP_lipid,
  FTKW_SWP_metal,
  FTKW_SWP_mod_res,
  FTKW_SWP_mutagen,
  FTKW_SWP_non_cons,
  FTKW_SWP_non_ter,
  FTKW_SWP_np_bind,
  FTKW_SWP_peptide,
  FTKW_SWP_propep,
  FTKW_SWP_repeat,
  FTKW_SWP_signal,
  FTKW_SWP_similar,
  FTKW_SWP_site,
  FTKW_SWP_strand,
  FTKW_SWP_thioeth,
  FTKW_SWP_thiolest,
  FTKW_SWP_transit,
  FTKW_SWP_transmem,
  FTKW_SWP_turn,
  FTKW_SWP_unsure,
  FTKW_SWP_variant,
  FTKW_SWP_varsplic,
  FTKW_SWP_zn_fing,
  FTKW_continuation
  }
DBLU_FTKW;

typedef enum DBLU_floctype
  {
  FLOC_unknown = 0,
  FLOC_complement,        /* "complement" */
  FLOC_join,              /* "join" */
  FLOC_order,             /* "order" */
  FLOC_group,             /* "group" */
  FLOC_one_of,            /* "one_of" */
  FLOC_replace,           /* "replace */
  FLOC_lpar,              /* '(' */
  FLOC_rpar,              /* ')' */
  FLOC_drange,            /* .. defined range */
  FLOC_irange,            /* in range: n.mm dot */
  FLOC_btween,            /* n^mm */
  FLOC_lthan,             /* < */
  FLOC_gthan,             /* > */
  FLOC_comma,             /* , */
  FLOC_nvalue,            /* a numeric value */
  FLOC_xentry,            /* id: cross reference */
  FLOC_ateol              /* run out of stuff */
  }
DBLU_FLOCTYPE;

typedef enum DB_ctype     /* classify characters */
  {
  DB_unkchr = 0,
  DB_alpha,               /* alphabetic */
  DB_numeric,             /* a decimal digit */
  DB_sep_scan,            /* a separator which can be scanned in wlu_chkwrd */
  DB_sep_immed,           /* a separator for immediate ident */
  DB_eos                  /* end of string */
  }
DB_CTYPE;

typedef enum DB_toktype   /* classify tokens depending on appearance */
  {
  DB_toknum = 0,          /* all numerical digits */
  DB_tokalph,             /* all alphabetic */
  DB_tokmixd,             /* mixed */
  DB_tok_mt,              /* 0 length */
  DB_tokunk               /* ??? */
  }
DB_TOKTYPE;

typedef enum DB_dlmttype  /* type of token delimiter */
  {
  DB_dlmt_lpar = 0,       /* '(' */
  DB_dlmt_rpar,           /* ')' */
  DB_dlmt_comma,          /* ',' */
  DB_dlmt_2dots,          /* '..' */
  DB_dlmt_colon,          /* ':' */
  DB_dlmt_not,            /* Not a recognised delimiter */
  DB_dlmt_eol             /* EOL */
  }
DB_DLMTTYPE;

typedef enum DB_locqual   /* note if gt or lt applies to a location */
  {
  DB_rq_none = 0,
  DB_rq_lthan,
  DB_rq_gthan,
  DB_rq_mt                /* null string */
  }
DB_LOCQUAL;

typedef enum DBLU_ftqual
  {
  FTQU_unknown = 0,
  FTQU_anticodon,
  FTQU_bound_moiety,
  FTQU_citation,
  FTQU_codon,
  FTQU_codon_start,
  FTQU_cons_splice,
  FTQU_db_xref,
  FTQU_direction,
  FTQU_EC_number,
  FTQU_evidence,
  FTQU_frequency,
  FTQU_function_qu,
  FTQU_gene,
  FTQU_label_qu,
  FTQU_mod_base,
  FTQU_note,
  FTQU_number,
  FTQU_organism,
  FTQU_partial,
  FTQU_phenotype,
  FTQU_product,
  FTQU_protid,
  FTQU_pseudo,
  FTQU_rpt_family,
  FTQU_rpt_type,
  FTQU_rpt_unit,
  FTQU_standard_name,
  FTQU_transl_except,
  FTQU_transl_table,
  FTQU_type_qu,
  FTQU_usedin,
  FTQU_locus_tag,
  FTQU_GO_info,
  FTQU_name,
  FTQU_biotype,
  FTQU_phosphothre,
  FTQU_phosphoser
  }
DBLU_FTQUAL;

typedef enum DB_rngqual
  {
  RQ_bothdefind = 0,
  RQ_p5undefind,
  RQ_p3undefind,
  RQ_p5p3undefind,
  RQ_p5patchd,
  RQ_unknown
  }
DB_RNGQUAL;

typedef enum DB_relpostype     /* relative position of items */
  {
  DB_rp_same = 0,              /* same position */
  DB_rp_3p,                    /* 3' WRT x */
  DB_rp_5p,                    /* 5' WRT x */
  DB_rp_not                    /* not related */
  }
DB_RELPOSTYPE;

typedef enum DB_lbltype
  {
  DBLBL_unk = 0,
  DBLBL_boffst,
  DBLBL_lblrec,
  DBLBL_notrec
  }
DB_LBLTYPE;

typedef enum DB_sptstat
  {
  DBPT_undef = 0,
  DBPT_5ptocds,
  DBPT_incds,
  DBPT_3ptocds
  }
DB_SPTSTAT;

typedef enum DB_natype
  {
  DBNA_unk = 0,
  DBNA_dna,
  DBNA_rna,
  DBNA_trna,
  DBNA_rrna,
  DBNA_mrna,
  DBNA_urna,
  DBNA_prot
  }
DB_NATYPE;

typedef enum DB_nastrand
  {
  DB_unkstrnd = 0,
  DB_snglstrnd,
  DB_dblstrnd,
  DB_mxdstrnd
  }
DB_NASTRAND;

typedef enum DB_oddwrds
  {
  DB_unknown = 0,
  DB_std,
  DB_unrev,
  DB_prelim,
  DB_unannot,
  DB_backbone,
  DB_circlr,
  DB_linear
  }
DB_ODDWRDS;

typedef enum DB_parstypes
  {
  PARS_numeric = 0,
  PARS_alpha,
  PARS_noise,
  PARS_separator,
  PARS_rangemod,
  PARS_ateoln
  }
DB_PARSTYPE;

typedef enum DB_sens      /* sense of a coding region */
  {
  DB_sens_unk = 0,
  DB_sens_norm,
  DB_sens_comp,
  DB_sens_mixd
  }
DB_SENS;

typedef struct db_str_elt   /* for linked list of strings */
  {
  char *strval;
  struct db_str_elt *nxtselt;  /* forward link */
  struct db_str_elt *prvselt;  /* reverse link */
  }
DB_STR_ELT;

typedef struct DB_segelt      /* segment of a coding seq */
  {
  int sgstart;
  int sgstop;
  int sgstart_2;              /* alternate start */
  DB_SENS sgsens;             /* sense of this segment */
  int lthan;                  /* 1 if <sgstart */
  int gthan;                  /* 1 if >sgstop */
  char *segid;                /* null if not mixed */
  struct DB_segelt *nextseg;  /* forward link */
  struct DB_segelt *prevseg;  /* back link */
  }
DB_SEGELT;

typedef struct DB_tblinfo    /* actually becomes a linked list */
  {
  DBLU_FTQUAL fqual;         /* this qualifier */
  char *qval;                /* its text */
  struct DB_tblinfo *nxtielt; /* forward link */
  struct DB_tblinfo *prvielt; /* back link */
  }
DB_TBLINFO;

typedef struct DB_lcinfo       /* like tblinfo, but for line class data */
  {
  DBLU_LINECLASS lclss;        /* class of this line */
  char *cinfo;                 /* strdup()ed info on that line */
  struct DB_lcinfo *nxtlcinf;  /* forward pointer */
  struct DB_lcinfo *prvlcinf;  /* back pointer */
  }
DB_LCINFO;

typedef enum DB_olaptyp      /* classify overlaps */
  {
  OLAP_none = 0,             /* no overlap of interest */
  OLAP_incl,                 /* one object is included in another */
  OLAP_5pabut,               /* 5' ends abutt, 3' included */
  OLAP_3pabut,               /* 3' ends abutt, 5' included */
  OLAP_abuts                 /* 5' & 3' ends abut */
  }
DB_OLAPTYP;

typedef enum DB_gbgene         /* "gene" definition varies */
  {
  DB_gene_unknown = 0,         /* undefined */
  DB_gene_genbank,             /* trad "genbank": /gene="xyz" qualifier */
  DB_gene_genomic,             /* genomic-style gene ft entries */
  DB_gene_refseq               /* refseq-style */
  }
DB_GBGENE;

typedef struct DB_featstrct  /* information about feature - incl. splicing */
  {
  DB_SEGELT *fstseg;          /* first segment */
  DB_SEGELT *lstseg;          /* last */
  char *idp;                  /* ptr to ent ename */
  DB_STR_ELT *enamelist;      /* list of ename values if multiples required */
  DB_STR_ELT *enamlstend;     /* end of that list */
  DBLU_FTKW featur;           /* what this feature is */
  int featno;                 /* number of this feature wrt whole entry */
  int savdno;                 /* number of feature wrt to "ofinterest" feats */
  DB_SENS strctsens;          /* sens of this whole feature */
  int cod_start;              /* codon start */
  SQT_GENCOD gencod;          /* genetic code table */
  int pseudogn;               /* 1 if pseudo gene qualifier found */
  int partcds;                /* 1 if qualified as partial CDS */
  DB_TBLINFO *infolist;       /* list of ft qualifiers */
  DB_TBLINFO *infoend;        /* end of list of ft qualifiers */
/* splicing-related details */
  struct DB_featstrct *mRNA;  /* points to any mRNA if defined, else NULL */
/* More stuff esp. for gff3 handling */
  char *parentname;          /* name from "Parent=";  Prob. malloc()ed, need free()ing if not NULL */
  struct DB_featstrct *parentfeat;  /* pointer to it */
/* list linkage info */
  struct DB_featstrct *nextst; /* forward link */
  struct DB_featstrct *prevst; /* back link */
  }
DB_FEATSTRCT;

typedef struct DB_accelt       /* accession no list elt */
  {
  char *accno;                 /* acc no itself */
  struct DB_accelt *nxtacc;    /* forward link */
  }
DB_ACCELT;

typedef struct DB_wantqu       /* structure for linked list of qual wants */
  {
  DBLU_FTQUAL fqu;             /* desired feat qualifier */
  struct DB_wantqu *nxtwntqu;  /* forward link */
  struct DB_wantqu *prvwntqu;  /* back link */
  }
DB_WANTQU;

typedef struct DB_ftypelt     /* for a linked list of feature types */
  {
  DBLU_FTKW fval;             /* feature keyword */
  DB_WANTQU *wquals;          /* qualifiers wanted for this feature */
  struct DB_ftypelt *nxtfelt; /* forward link */
  struct DB_ftypelt *prvfelt; /* back link */
  }
DB_FTYPELT;

typedef enum DB_spliclmt       /* control splicing behaviour */
  {
  DBS_allres = 0,              /* return any possible residue */
  DBS_segonly,                 /* return base if in actual or mrna segment */
  DBS_featlimit,               /* limit to actual feature */
  DBS_mrnalimit                /* limit to a linked mRNA, if exists */
  }
DB_SPLICLMT;

typedef struct DB_entstrct     /* info for complete entry */
  {
  char *ename;                 /* entry name */
  unsigned long eoffst;        /* offset of this entry start */
  int nfeats;                  /* no of features stored */
  int nseen;                   /* no of features seen */
  DB_ACCELT *acclst;           /* list of accession numbers */
  int sqlen;                   /* length of this sequence */
  char *seq;                   /* sequence data malloc()ed as needed */
  char *sfpt;                  /* current length of sequence (fill pointer) */
  DB_FEATSTRCT *featlst;       /* list of features */
  DB_FEATSTRCT *lastfeat;      /* end of feature list */
  DB_NATYPE natype;            /* type of this sequence */
  DB_NASTRAND strnd;           /* strandedness of entry */
  char *descr;                 /* description line */
  char *date;                  /* date of entry */
  char *nid;                   /* nid if present */
  int circ;                    /* 1 if circular sequence */
  DB_FTYPELT *intrst;          /* list of feature types of interest */
/*  DB_FTYPELT *fquals;       */   /* list of feature qualifiers to save */
  DBLU_DBFMT efmt;             /* format of this entry */
/* stuff relating to source and scanning buffers */
  FILE *sfl;                   /* the source file handle */
  char *slbuf;                 /* source line buf */
  int slcnt;                   /* source line count */
  int slblen;                  /* length of slbuf */
  char *cp;                    /* current position in line buf */
  DBLU_LINECLASS curltyp;      /* type of current line */
  DBLU_LINECLASS prvltyp;      /* type of previous identifiable line */
  int onqual;                  /* true once quals are found - stop odd nos */
/* for parsing and splicing operations */
  DB_SEGELT *curseg;           /* current segment */
  DB_SENS cursens;             /* current sens */
/* for splicing operations */
  DB_FEATSTRCT *splcft;        /* struct for splice operations */
  int sqpt;                    /* present position */
  int p5limt;                  /* set to 5' limit of actual feature */
  int p3limt;                  /* set to 3' limit of actual feature */
/* for translating/checking */
  SQT_GENCOD curgcod;          /* current genetic code */
  TRANS_MATRX *trnsmatrx;      /* current gcode matrix */
  COD_VECT_STRCT **cvecs;      /* allocated set of vectors for translation */
/* parsing stuff for convenience */
  WRD_LUSTRCT *lclss;          /* raw line classes */
  WRD_LUSTRCT *nawlu;          /* identify NA type */
  WRD_LUSTRCT *oddwlu;         /* identify odd words */
  WRD_LUSTRCT *strndwlu;       /* identify na strandedness */
  WRD_LUSTRCT *fkwlu;          /* identify feature keywords */
  WRD_LUSTRCT *flocwlu;        /* identify location keywords */
  WRD_LUSTRCT *fqulwlu;        /* identify feat qualifiers */
/* housekeeping */
  int freefeatids;             /* 1 if these are malloc()ed and need freeing */
/* lineclass information */
  DB_LCINFO *lcinfo;           /* linked list of line classes & assoc values */
  int seqversno;               /* integer, hopefully */
  }
DB_ENTSTRCT;

typedef struct DB_splicinfo    /* necessary information for a splice position */
  {
  DB_ENTSTRCT *cacheentp;      /* may as well cache pointer to original entry */
  DB_FEATSTRCT *cachefp;       /* current feature pointer */
  DB_SEGELT *cachesegp;        /* current segment */
  int cachesqpt;               /* current seq pointer */
  int lim5p;                   /* 5' limit of feature */
  int lim3p;                   /* 3' limit */
  DB_SENS cachesens;           /* current sense setting */
  }
DB_SPLICINFO;

typedef enum DBT_dbmtype  /* file/id/acc */
  {
  DBT_db_file = 0,
  DBT_db_id,
  DBT_db_acc,
  DBT_db_versn
  }
DBT_DBMTYPE;

typedef struct DBP_multi_elemnt  /* data for a multi entry scan */
  {
  char *me_id;   /* id for this element */
  DB_ENTSTRCT *me_entstrptr;   /* ptr to entry data structure */
  WRD_LUSTRCT *genids;               /* quick look up for gene names */
  struct DBP_multi_elemnt *nxt_dme;  /* fwd ptr for list */
  struct DBP_multi_elemnt *prv_dme;  /* back ptr for list */
  }
DBP_MULTI_ELEMNT;

typedef enum DBP_relpostype  /* ways that entities can be positionally related */
  {
  DBP_relpos_undef  = 0,
  DBP_relpos_5pto,
  DBP_relpos_5polap,
  DBP_relpos_internal,
  DBP_relpos_included,
  DBP_relpos_3polap,
  DBP_relpos_3pto
  }
DBP_RELPOSTYPE;

/* function headers */

DBLU_DBFMT db_chr2dbfmt(char c);
  /* return the database format corresponding to this character */

char db_dbfmt2chr(DBLU_DBFMT fmt);
  /* return a letter for database format fmt */

int db_logerr(char *msg,
              int lno,
              char *lin);
  /* at present just write to stderr, return 0 to allow logical test */

WRD_LUSTRCT *db_getclasstrct(DBLU_DBFMT dfmt);
  /* depending on dfmt, create and fill the datastructure for lineclass 
parsing, returning the address of the structure created.  Invalid dfmt
will return NULL.  malloc failures will die */
  
char *db_ftkw2str(DBLU_FTKW kw);
  /* return a string corresponding to the kw given */

char *db_ftkw2text(DBLU_FTKW kw);
  /* return a long string corresponding to the kw given */

WRD_LUSTRCT *db_getkwstrct(DBLU_DBFMT dfmt);
  /* create and fill datastructure for feature table key word parsing,
depending on dfmt.  return pointer to that structure, NULL for invalid
dfmt  */

WRD_LUSTRCT *db_getftlocstrct(void);
  /* create and fill the datastructure for location word 
parsing, returning the address of the structure created. */

WRD_LUSTRCT *db_getftqualstrct(void);
  /* create and fill the datastructure for location qualifier word 
parsing, returning the address of the structure created. */

WRD_LUSTRCT *db_getnatypestrct(void);
  /* create and fill the datastructure for NA_TYPE words, return
 the address of the structure created. */

WRD_LUSTRCT *db_getoddwrdsdstrct(void);
  /* create and fill the datastructure for odd words, return
 the address of the structure created. */

WRD_LUSTRCT *db_getstrndstrct(void);
  /* create and fill the datastructure for strandedness words, return
 the address of the structure created. */

DBLU_LINECLASS db_classfylin(char *lin,
                             DBLU_DBFMT dfmt,
                             WRD_LUSTRCT *ws);
/* return the class of this line */

DBLU_LINECLASS db_classfyentln(DB_ENTSTRCT *es);
  /* return the line classification of es line */

char *db_linclass2str(DBLU_LINECLASS lc);
  /* return a string for the linetype - (don't take account of database
format at this point */

char *db_oddwrd2str(DB_ODDWRDS ow);
  /* return a string for the odd word ow */

char *db_natyp2str(DB_NATYPE nat);
  /* return a string for the type of this na */

char *db_strnd2str(DB_NASTRAND nas);
  /* return a string for nas */

void db_initfeatstrct(DB_FEATSTRCT *fs,
                      char *id);
  /* set initial values for this feature structure */

DB_STR_ELT *dbp_appstrelt(DB_STR_ELT **spt,
                          char *str);
/* malloc and append a string element to list *spt.
return the address of the new element.
The new value is strduped, to avoid unexpected
changes/behavior */

void db_delstrelt(DB_STR_ELT *sp,
                  DB_STR_ELT **lstrt,
                  DB_STR_ELT **lend);
/* delete sp from list *lstrt..*lend */

void db_killstrlst(DB_STR_ELT **lstrt,
                   DB_STR_ELT **lend);
/* iteratively remove the header from *lstrt */

DB_STR_ELT *db_matstrelt(DB_STR_ELT *strlst,
                         char *qstr,
                         int (* scmpfun)(const char *x1,
                                         const char *x2));
/* return an element in strlst matching qstr, using
scmpfun for comparison.  return NULL if none */

int db_cntstrlst(DB_STR_ELT *strlst);
  /* recursively count strlst elements */

int db_cntstrlst(DB_STR_ELT *strlst);
  /* recursively count strlst elements */

DB_SEGELT *db_appnseg(DB_SEGELT **spt,
                      int strt,
                      int stop,
                      DB_SENS sens,
                      char *xidnm);
/* malloc and append a new segment to list *spt.  return the address of the new
element */

void db_delsegelt(DB_SEGELT *sp,
                  DB_SEGELT **lstrt,
                  DB_SEGELT **lend);
/* delete sp from list *lstrt..*lend */

void db_killsegs4featstrct(DB_FEATSTRCT *fs);
  /* iteratively delete the segment list for fs - memory is freed */

DB_FEATSTRCT *db_appnfeat(DB_FEATSTRCT **spt,
                          DBLU_FTKW feat,
                          char *id);
/* malloc and append a new featstruct to list *spt.  return the address of 
the new element */

DB_FEATSTRCT *dbp_prependfeatstrct(DB_FEATSTRCT **flststrtp,
                                   DB_FEATSTRCT *spt,
                                   DBLU_FTKW feat,
                                   char *id);
/* malloc and prepend a new featstruct before element *spt.  return the address of 
the new element */

void db_killsegs4featstrct(DB_FEATSTRCT *fs);
  /* iteratively delete the segment list for fs - memory is freed */

void db_delfeatelt(DB_FEATSTRCT *sp,
                   DB_FEATSTRCT **lstrt,
                   DB_FEATSTRCT **lend);
/* delete sp from list *lstrt..*lend */

void db_killallfeats(DB_FEATSTRCT **ffst,
                     DB_FEATSTRCT **flst);
  /* iteratively delete all of *flst */

DB_TBLINFO *db_appnielt(DB_TBLINFO **spt,
                        DBLU_FTQUAL fql,
                        char *info);
/* malloc and append a information element to list *spt.
  return the address of the new element */

void db_delielt(DB_TBLINFO *sp,
                DB_TBLINFO **lstrt,
                DB_TBLINFO **lend);
/* delete sp from list *lstrt..*lend */

void db_killtblinfolst(DB_TBLINFO **lstrt,
                       DB_TBLINFO **lend);
/* iteratively remove the header from *lstrt */

DB_TBLINFO *db_tblent4udata(DB_TBLINFO *tlst,
                            DBLU_FTQUAL uqul,
                            char *ustr);
/* scan tlst for the first element which matches uqul and (if non-NULL) 
contains ustr.  Return a pointer to that element, NULL for no match */

DB_LCINFO *db_appnlcelt(DB_LCINFO **spt,
                        DBLU_LINECLASS lcl,
                        char *info);
/* malloc and append a information element to list *spt.
  return the address of the new element */

void db_dellcelt(DB_LCINFO *sp,
                 DB_LCINFO **lstrt);
/* delete sp from list *lstrt */

void db_killlclsslst(DB_LCINFO **lstrt);
/* iteratively remove the header from *lstrt */

DB_LCINFO *db_lcent4udata(DB_LCINFO *tlst,
                          DBLU_LINECLASS uqul,
                          char *ustr);
/* scan tlst for the first element which matches uqul and (if non-NULL) 
contains ustr.  Return a pointer to that element, NULL for no match */

void db_endfeatsplic(DB_ENTSTRCT *es);
  /* revert es splice components to an un-initiallised state */

DB_ACCELT *db_fndaccno(DB_ACCELT *alst,
                       char *acc);
/* return a pointer to any element of alst which matches acc.  Don't do
any case conversion */

void db_initentstrct(DB_ENTSTRCT *es,
                     DBLU_DBFMT fmt,
                     WRD_LUSTRCT *ftkwlu);
/* initialise contents of an entry structure */

char *db_strdupskip(char *lin,
                    int smax,
                    char *skipset);
/* produce a strdup-ed copy of lin up to smax chars, stopping at first non 
skipset char after object.  If leading skipset chars, then skip them also.
  return new copy of string */

void db_mltdup_quotd(DB_ENTSTRCT *es,
                     DB_FEATSTRCT *fp,
                     DBLU_FTQUAL qul);
/* look on current line for a quoted string, append each such line to
fp info list. */

void db_unpickidline(char *slin,
                     DB_ENTSTRCT *es,
                     DBLU_DBFMT fmt);
/* examine slin as a locus/ident line and extract useful info into es */

void db_killacclst(DB_ACCELT **alst);
  /* iteratively delete the accno list alst - memory is freed */

void db_dispos_ent(DB_ENTSTRCT **es);
  /* lose all allocated storage associated with *es and set to NULL */

void db_unpickaccno(char *lin,
                    DB_ACCELT **alst);
/* scan lin for a list of accession nos, append strdup-ed copy to *alst */

void db_tellacclst(FILE *fl,
                   DB_ACCELT *alp,
                   char *sep);
/* iteratively write alp names to fl */

DB_FTYPELT *db_appnfelt(DB_FTYPELT **spt,
                        DBLU_FTKW eval);
/* malloc and append a new value to list *spt.  return the address of the new
element */

void db_addcntfelts(DB_FTYPELT **clst,
                    int lcls,
                    int cnt);
/* append lcls records cntX onto *clst */

int db_delfelt(DB_FTYPELT **fstrt,
               DB_FTYPELT *fp);
/* fp is an element of *fstrt, delete it, correct pointers */

void db_delfstfelt(DB_FTYPELT **lstrt);
/* delete *lstrt */

void db_killfeltlst(DB_FTYPELT **alst);
  /* iteratively delete the feature list alst - memory is freed */

DB_FTYPELT *db_ptr4felt(DB_FTYPELT *flst,
                        DBLU_FTKW eval);
/* return the (first) address of any element of flst which matches eval */

int db_delmatfelt(DB_FTYPELT **fstrt,
                  DBLU_FTKW eval);
/* delete all elements of *fstrt which match eval. return no deleted */

DB_WANTQU *db_appnquelt(DB_WANTQU **spt,
                        DBLU_FTQUAL qval);
/* malloc and append a new value to list *spt.  return the address of the new
element */

DB_WANTQU *db_appnquelt2felt(DB_FTYPELT *fep,
                             DBLU_FTQUAL qval);
/* malloc and append a new qualifier value for feature element fep.
return the address of the new element */

DB_WANTQU *db_appnquelt4feat(DB_FTYPELT *flst,
                             DBLU_FTKW fkw,
                             DBLU_FTQUAL qval);
/* scan flst for a value matching fkw, malloc and append a new qualifier value
for it. return the address of the new element */

int db_appnquelt4all(DB_FTYPELT *flst,
                     DBLU_FTQUAL qval);
/* malloc and append a new qualifier value for all flst elts.
return number of elements entered */

DB_WANTQU *db_ptr4qualval(DB_WANTQU *qlst,
                          DBLU_FTQUAL qv);
/* if an element matching qv exists on qlst, return its address */

void db_delwantquelt(DB_WANTQU *sp,
                     DB_WANTQU **lstrt);
/* delete sp from list *lstrt */

void db_beheadquallst(DB_WANTQU **qlst);
  /* remove the element *qlst */

void db_killwntqulst(DB_WANTQU **lstrt);
/* iteratively remove the header from *lstrt */

DB_WANTQU *db_wntquptr4featnqul(DB_FTYPELT *flst,
                                DBLU_FTKW ftkw,
                                DBLU_FTQUAL fq);
/* return a pointer to the element matched by ftkw & fq, NULL otherwise */

char *db_fgets(char *buf,
               int blen,
               FILE *fl,
               int *lcnt);
/* perform fgets of fl, remove \n from string */

DB_CTYPE db_classchr(char c);
  /* return a character classification of c */

DBLU_FTKW db_clssfeatln(DB_ENTSTRCT *es);
  /* return (if any) the feature type of es */

int db_rnginclds(int p1,
                 int p2,
                 int p3);
/* return 1 if p1..p3 includes p2 (p1 <= p2 <= p3) or (p3 <= p2 <= p1) */

DBP_RELPOSTYPE dbp_relatepos2feat(DB_FEATSTRCT *fp,
                                  int qpos);
/* indicate the relative position of qpos wrt feat fp */

DBP_RELPOSTYPE dbp_relatefeat2posns(DB_FEATSTRCT *fp,
                                    int p5p,
                                    int p3p);
/* p5p and p3p are 5' & 3' positions in absolute sense.
return the relative position of that feature wrt
p5p..p3p */

char *dbp_relpostype2str(DBP_RELPOSTYPE rpt);
  /* text for rpt */

int db_savdno4feat(DB_FEATSTRCT *fp);
  /* return the saved number for this feature - should correspond to old-style
feature numbering for fish_term */

int db_featno4feat(DB_FEATSTRCT *fp);
  /* return the feature number for this feature - relative to whole feature
table */

DB_FEATSTRCT *db_nthfeat(DB_ENTSTRCT *es,
                         int nth,
                         int (*numbfn)(DB_FEATSTRCT *fs));
/* return ptr to the nth feature.  (*numbfn)() relates which type of numbering
to use.  NULL if not found */

DB_SEGELT *db_5pseg4feat(DB_FEATSTRCT *fp);
  /* return the 5' segment for this feature, NULL if non-existent */

DB_SEGELT *db_3pseg4feat(DB_FEATSTRCT *fp);
  /* return the 3' segment for this feature, NULL if non-existent */

int db_5ppos4seg(DB_SEGELT *sp);
  /* return the 5' position for this sequence segment */

int db_3ppos4seg(DB_SEGELT *sp);
  /* return the 3' position for this sequence segment */

int db_5ppos4feat(DB_FEATSTRCT *fp);
  /* return the 5' position for this feature */

int db_3ppos4feat(DB_FEATSTRCT *fp);
  /* return the 3' position for this feature */

int db_5pextnt4feat(DB_FEATSTRCT *fp);
  /* return the 5' extent of feature fp (i.e. its minimum position) */

int db_3pextnt4feat(DB_FEATSTRCT *fp);
  /* return the 3' extent of feature fp (i.e. its maximum position) */

void db_featextremes(DB_ENTSTRCT *es,
                     DB_FEATSTRCT *fp,
                     int *e5p,
                     int *e3p);
/* return to e5p & e3p (if non-null) the extreme 
extent of fp, noting that mixed splice features
misfire the usual checks */

DB_SEGELT *db_nxtseg4feat(DB_SEGELT *cseg,
                          DB_FEATSTRCT *fp);
/* return the next segment in order for this feature */

DB_SEGELT *db_prvseg4feat(DB_SEGELT *cseg,
                          DB_FEATSTRCT *fp);
/* return the previous segment in order for this feature */

DB_RELPOSTYPE db_segorder(DB_SEGELT *s1,
                          DB_SEGELT *s2,
                          DB_FEATSTRCT *fp);
/* return relationship of s1 & s2.  NOTE that if s1 is 3' to s2, then
DB_rp_3p is returned */

int db_is5pto(int p1,
              int p2,
              DB_SENS sens);
/* return 1 if p1 is 5' to p2, depending on sense.
  NOTE: returns 0 for p1=p2 */

int db_segis5pto(DB_SEGELT *sp1,
                 DB_SEGELT *sp2);
/* return 1 if sp1 is 5' wrt sp2 */

void db_saysegposns(FILE *fl,
                    DB_SEGELT *sp);
/* tell fl about all of this segment */

int db_featlength(DB_FEATSTRCT *fs,
                  int *segcnt);
/* scan this feature for length, if segcnt is non-NULL, return the count */

int db_rawfeatlength(DB_FEATSTRCT *fp);
  /* return the raw length of fp, irrespective of orientation.  ignores
introns */

int db_segcnt4feat(DB_FEATSTRCT *fs);
  /* return the number of segments in fs */

DB_OLAPTYP db_chkfeatolap(DB_FEATSTRCT *f1,
                          DB_FEATSTRCT *f2);
/* return an assessment of the feature overlap for f2 included in f1 */

DB_SEGELT *db_seglstelt4pos(DB_SEGELT *seglst,
                            int pos);
/* return the seg element pointer which contains pos, NULL
if none */

DB_SEGELT *db_featseg4pos(DB_FEATSTRCT *fs,
                          int pos);
/* return the feature segment of fs which contains pos.  NULL if none */

DB_SEGELT *db_nxtfeatseg4pos(DB_FEATSTRCT *fs,
                             int pos);
/* if pos is on a segment, return that, else traverse
the segment list, returning the next segment after
pos in the sense direction.  NULL if isn't in fs range */

DB_FEATSTRCT *db_nxtfeat4pos(DB_FEATSTRCT *flst,
                             DBLU_FTKW ftype,
                             int gpos);
/* return pointer to the next feature including present
that embraces gpos.  NULL for none.  Match on ftype unless
it is FTKW_unknown. */

DB_FEATSTRCT *db_feat4pos(DB_ENTSTRCT *ep,
                          DBLU_FTKW ftype,
                          int gpos);
/* return the first feature of *ep that includes
position gpos.  Match on ftype unless
it is FTKW_unknown. */

DB_FEATSTRCT *db_nxtfeat4typenid(DB_FEATSTRCT *ftlst,
                                 DBLU_FTKW fttype,
                                 char *queryid,
                                 int (*scmpfn)(const char *x1,
                                               const char *x2));
/* look thru ftlst for next feature with an id
matching queryid as indicated by (scmpfn*)() and of 
type fttype (if this is not unknown) */

DB_FEATSTRCT *db_nxtfeat4id(DB_FEATSTRCT *ftlst,
                            char *queryid,
                            int (*scmpfn)(const char *x1,
                                          const char *x2));
/* look thru ftlst for next feature with an id
matching querid as indicated by (*scmpfn)() */

int db_segsoverlap(DB_SEGELT *seglst);
  /* check if any segments in seglst overlap */

DB_FEATSTRCT *db_fndcdsandmrna4pos(DB_FEATSTRCT *fstrt,
                                   int gpos,
                                   DB_SPTSTAT *postype);
/* from feature list fstrt try to find a CDS with or without mrna which includes
gpos.  Return the cds feature pointer for this.  if
postype is non-NULL, return a position status value to it */

DB_FEATSTRCT *dbp_featincludesfeat(DB_FEATSTRCT *fp,
                                   DB_FEATSTRCT *flst,
                                   DBLU_FTKW ftype);
/* return the structure address of the first feature in flst which
is of type ftype (if not FTKW_unknown) and is included in feature
fp */

char *db_sptstat2str(DB_SPTSTAT sstat);
  /* return a relevant string for sstat */

int db_segcnt4featseg(DB_FEATSTRCT *fp,
                      DB_SEGELT *sep);
/* count through fp segments, returning ordinal
position of sep in elements.  0 if not present */

int db_segincldseg(DB_SEGELT *sout,
                   DB_SEGELT *sin);
/* return true if sin is included in sout */

int db_restsplicsame(DB_FEATSTRCT *fout,
                     DB_FEATSTRCT *fin,
                     DB_SEGELT *oseg,
                     DB_SEGELT *iseg,
                     int scnt);
/* recursively check if oseg.. is a valid identical splice for iseg */

int db_featsplicsame(DB_FEATSTRCT *fout,
                     DB_FEATSTRCT *fin);
/* check that fin and fout use the same splicing pattern over the range of
fin */

int db_intrncnt4seg(DB_FEATSTRCT *fs,
                    DB_SEGELT *sp,
                    int to3p);
/* return count of introns from sp to 3' (to3p) or 5' (!to3p) direction */

int db_intrncnt4pos(DB_FEATSTRCT *fs,
                    int posn,
                    int to3p);
/* return count of the introns from posn to end of fs in 3' (to3p) or 
5' (!to3p) direction */

int db_intrnsinrng(DB_FEATSTRCT *fs,
                   int pstrt,
                   int pstop);
/* return count of the fs introns from pstrt to pstop */

DB_FEATSTRCT *db_gene4mrna(DB_ENTSTRCT *es,
                           DB_FEATSTRCT *mfs);
/* look through es feat list and return the first gene feature that 
encompasses mfs (must have same sense and start&end match) */

DB_FEATSTRCT *db_mrna4feat(DB_ENTSTRCT *es,
                           DB_FEATSTRCT *fs);
/* look through es feat list and return the first mRNA that encompasses fs
 (must have same sense and the splicing pattern must match) */

int db_scnmrnas4feats(DB_ENTSTRCT *es,
                      DBLU_FTKW feat);
/* scan through feature list, for any of type feat, link up any mRNA's which
correspond to them and return count of No so matched */

int db_countfeats(DB_ENTSTRCT *es,
                  DBLU_FTKW feat);
/* scan through feature list, for any of type feat - return count.
if feat == FTKW_unknown, return count of all features */

char *db_olaptyp2str(DB_OLAPTYP ol);
  /* return a string for ol */

void db_initftsplic4pos(DB_ENTSTRCT *es,
                        DB_FEATSTRCT *fs,
                        int spos);
/* try to initiallise seq splicing for fs at position spos.  Set relevant
es parameters, even for invalid positions */

void db_chktransmatrx(DB_ENTSTRCT *es,
                      DB_FEATSTRCT *fs);
/* check (and if necessary init) that a translation
matrix exists */

void db_init5psplic4feat(DB_ENTSTRCT *es,
                         DB_FEATSTRCT *fs);
/* initialise extracting the spliced seq data for fs alone */

void db_init3psplic4feat(DB_ENTSTRCT *es,
                         DB_FEATSTRCT *fs);
/* initialise splicing for 3' end of fs */

void db_nxtesptr(DB_ENTSTRCT *es);
  /* increment the splice pointer in es to the next 5'->3' position */

void db_prvesptr(DB_ENTSTRCT *es);
  /* increment the splice pointer in es to the next 3'->5' position */

char db_getesresnoinc(DB_ENTSTRCT *es,
                      DB_SPLICLMT segsonly,
                      char (* bascmpfun)(char yb,
                                           BAS_MATMODE xmmod));
/* guts of the process to return the appropriate residue, irrespective any any complex links,
but don't increment/decrement seq pointer.
Depending on segsonly, return
null if there is not a valid segment or if pointer is out of feature limit.
(* bascmpfun)() is used to create base complement */

char db_getesrescore(DB_ENTSTRCT *es,
                     DB_SPLICLMT segsonly,
                     void (* esincdecfun)(DB_ENTSTRCT *xs),
                     char (* bascmpfun)(char yb,
                                          BAS_MATMODE xmmod));
/* guts of the process to return the appropriate residue, irrespective any any complex links,
  use esincdecfun() to adjust the seq pointer. depending on segsonly, return
  null if there is not a valid segment or if pointer is out of feature limit.
  (* bascmpfun)() is used to create base complement */

char db_getesres(DB_ENTSTRCT *es,
                 DB_SPLICLMT segsonly,
                 void (* esincdecfun)(DB_ENTSTRCT *xs));
/* return the appropriate residue, irrespective any any complex links,
  use esincdecfun() to adjust the seq pointer. if segsonly, then return
  null if there is not a valid segment */

char db_getnxtesres(DB_ENTSTRCT *es,
                    DB_SPLICLMT segsonly,
                    char (* bascmpfun)(char yb,
                                         BAS_MATMODE xmmod));
/* return the appropriate residue, irrespective any any complex links,
  use esincdecfun() to adjust the seq pointer. if segsonly, then return
  null if there is not a valid segment. use (*bascmpfun)() for reverse
  strand residues */

char db_getprvesres(DB_ENTSTRCT *es,
                     DB_SPLICLMT segsonly,
                     char (* bascmpfun)(char yb,
                                          BAS_MATMODE xmmod));
/* return the previous residue, irrespective any any complex links,
  use esincdecfun() to adjust the seq pointer. if segsonly, then return
  null if there is not a valid segment. use (*bascmpfun)() for reverse
  strand residues */

int db_getescodn(DB_ENTSTRCT *es,
                 DB_SPLICLMT splimt,
                 char *codn,
                 void (* esincdecfun)(DB_ENTSTRCT *xs));
/* use db_getesres() to fill a codon (codn assumed to be long enough) -
return the number of valid bases inserted */

char db_nxtsplcres(DB_ENTSTRCT *es,
                   DB_SPLICLMT segsonly);
/* return the next base in sequence for presently-init-ed feature.
  if segsonly, then return '\0' if the position is not on segments */

char db_prvsplcres(DB_ENTSTRCT *es,
                   DB_SPLICLMT segsonly);
/* return the previous base in sequence for presently-init-ed feature.
  if segsonly, then return '\0' if the position is not on segments */

int db_segnoinfeat(DB_FEATSTRCT *fp,
                   DB_SEGELT *sep);
/* return the segment number from start of feature of sep.
0 if not observed in fp.  */

void db_cntshift(DB_ENTSTRCT *es,
                 int nshft,
                 char (* getbasfn)(DB_ENTSTRCT *xs,
                                   DB_SPLICLMT sgsonly));
/* move nshft residues in direction dicated by (* getbasfn)() - presently by 
iteratively calling it */

int db_cntlmtshift(DB_ENTSTRCT *es,
                   int nshft,
                   DB_SPLICLMT slmt,
                   char (* getbasfn)(DB_ENTSTRCT *xs,
                                       DB_SPLICLMT sgsonly));
/* move nshft residues in direction dicated by (* getbasfn)() - presently by 
iteratively calling it, Uses slmt to determine splice behaviour.  return
1 if shift was completed OK. */

void db_cntnolmtshift(DB_ENTSTRCT *es,
                      int nshft,
                      DB_SPLICLMT slmt,
                      char (* getbasfn)(DB_ENTSTRCT *xs,
                                          DB_SPLICLMT sgsonly));
/* move nshft residues in direction dicated by (* getbasfn)() - presently by 
iteratively calling it. Does not limit pointer adjustment at seq limits.
Uses slmt to determine splice behaviour. */

char db_prvsplcresnolmt(DB_ENTSTRCT *es,
                        DB_SPLICLMT segsonly);
/* return the previous base in sequence for presently-init-ed feature.
  if segsonly, then return '\0' if the position is not on segments.
The seq pointer will be decremented whether we are on sequence or
not */

void db_fwdshift(DB_ENTSTRCT *es,
                 int nshft);
/* iteratively call db_nxtsplcres to move nshft positions forward (5'->3') */

void db_bakshift(DB_ENTSTRCT *es,
                 int nshft);
/* iteratively call db_prvsplcres to move nshft positions back (3'->5') */


int db_splcft2cdn(DB_ENTSTRCT *es,
                  DB_SPLICLMT segsonly,
                  char *cdn);
/* try to fill cdn with 3 successive calls to db_nxtsplcres().  Return
the number of non-null characters returned */

int db_rvsplcft2cdn(DB_ENTSTRCT *es,
                    DB_SPLICLMT segsonly,
                    char *cdn);
/* try to fill cdn with 3 successive calls to db_prvsplcres().  Return
the number of non-null characters returned */

void db_marksplicposn(DB_ENTSTRCT *esp,
                      DB_SPLICINFO *spinfo);
/* cache various esp values into spinfo to store
the present splice information for later restoration */

int db_recovrsplicposn(DB_ENTSTRCT *esp,
                       DB_SPLICINFO *spinfo);
/* restore stored splice parameters from *spinfo to
esp.  Include sanity check that esp matches stored
entry pointer */

int db_chkft4term(DB_ENTSTRCT *es,
                  DB_FEATSTRCT *fs,
                  int assrtnxt);
/* establish if the final codon of fs is a terminator: if not then on basis of
assrtnxt check the next codon and if it is, then correct fs value */

int db_ft2imatrx(DB_ENTSTRCT *es,
                 DB_FEATSTRCT *fs,
                 DB_SPLICLMT slmt,
                 TRANS_MATRX cnts);
/* count all codons for a feature into cnts,  return the number of valid
 codons recovered */

int db_ft2imatrxrev(DB_ENTSTRCT *es,
                    DB_FEATSTRCT *fs,
                    DB_SPLICLMT slmt,
                    TRANS_MATRX cnts);
/* count all codons for a feature into cnts,  return the number of valid
 codons recovered */

int db_rescnt4imatrx(DB_ENTSTRCT *es,
                     char pres,
                     TRANS_MATRX cnts);
/* return the total number of pres codes for cnts, using the 
pre-built codon vector table of es */

int db_termcnt4imatrx(DB_ENTSTRCT *es,
                      TRANS_MATRX cnts);
/* return the total number of terminator codes for cnts, using the 
pre-built codon vector table of es */

int db_totimtrx(DB_ENTSTRCT *es,
                TRANS_MATRX ccnts);
/* total the count of valid residues in ccnts, using genetic code defined
for es */

void prt_mxd_tdata(FILE *strm,
                   TRANS_MATRX trnsar,
                   TRANS_MATRX datarr);
/* tell strm about trnsar contents in structured table for valid codons */

char *db_rngqual2str(DB_RNGQUAL rq);
  /* return a string for rq */

DB_RNGQUAL db_feat2rngqual(DB_FEATSTRCT *fp);
  /* look at range qualifiers of fp segment ends and return the classification
for them */

int db_lmt4feat(DB_FEATSTRCT *fp,
                int (*maxminfn)(int x,
                                int y));
/* return the extent of fp as established by maxminfn(). */

int db_directeddst(DB_FEATSTRCT *f1,
                   DB_FEATSTRCT *f2,
                   int p3dst);
/* determine f1..f2 distance, depending on p3dst.  Return -ve for non-sensible
  configurations */

DB_FEATSTRCT *db_nxtfeat(DB_ENTSTRCT *es,
                         DB_FEATSTRCT *fp,
                         DB_FTYPELT *flst,
                         int p3dst);
/* return the next 3' or 5' feature to fp.  return NULL if no such feature 
 exists.  Scan entire featlst, so we don't make any assumptions about feature
 ordering. */

int db_dst2nxtfeat(DB_ENTSTRCT *es,
                   DB_FEATSTRCT *fp,
                   DB_FTYPELT *flst,
                   int p3dst);
/* return the distance to the next 3' or 5' feature to fp.  return remaining 
 free sequence if no such feature  exists.  Scan entire featlst, so we 
 don't make any assumptions about feature ordering. */

int db_dist4pos2segend(DB_SEGELT *sep,
                       int pos,
                       int p3dist);
/* assuming pos is in range of sep, return
either the distance to 5' or 3' extremity
depending on p3dist.  return -1 if
pos is not in seg */

int db_splicdist4feat2pos(DB_FEATSTRCT *fp,
                          int spos,
                          int p3dist);
/* return spliced distance from 5' or 3' end of 
feat fp to spos.  -1 if spos not on segment */

void db_sayafeatcore(FILE *fl,
                     DB_ENTSTRCT *es,
                     DB_FEATSTRCT *fs,
                     int showseq);
/* tell fl about fs contents. showseq controls
whether any sequence is to be shown */

void db_sayafeat(FILE *fl,
                 DB_ENTSTRCT *es,
                 DB_FEATSTRCT *fs);
/* tell fl about fs contents */

DB_FTYPELT *db_str2lst(char *ustr);
  /* scan ustr for a list of seq no integers - return pointer to a linked list
of values specified therein.
  Input can be a string of integers, separated by non-digit chars.
  n..m or n-m is interpreted as all ints from n to m inclusive */

DBLU_FTKW db_scnnafeat(DB_ENTSTRCT *es);
/* attempt to identify this line as a NA feature descriptor.  if wanted, append
it to es->featlst.  In any event, return the feature identified */

DBLU_FTKW db_scnswfeat(DB_ENTSTRCT *es);
/* attempt to identify this line as a SwissProt feature descriptor.
if wanted, append it to es->featlst.  In any event, return the feature 
identified */

int db_crossplccnt(DB_FEATSTRCT *fs);
  /* return the number of cross entry splices for this feature */

int db_scnfeatwrds(char *intwrds,
                   WRD_LUSTRCT *kwlu,
                   DB_FTYPELT **flst,
                   int wldok);
/* scan intwrds for space or comma separated set of words, try to parse by
kwlu and create linked list of entries.  Allow and respond to wildcard strings
if wldok */

int db_scnqulwrds(char *intwrds,
                  WRD_LUSTRCT *kwlu,
                  DB_WANTQU **flst,
                  int wldok);
/* scan intwrds for space or comma separated set of words, try to parse by
kwlu and create linked list of entries.  Allow and respond to wildcard strings
if wldok */

DB_ENTSTRCT *db_parsemain(FILE *df,
                          WRD_LUSTRCT *lclass,
                          DBLU_DBFMT fmt,
                          DB_FTYPELT *ufwrds,
                          int ndxrun,
                          DB_LCINFO *linflst);
/* now the real guts of the parsing process.  This version
accepts a linked list of lineclass data which, if non-NULL, will
be used to capture information from relevant lines */
/* scan df from current position to next sq terminator.  wlus is used for
key words.  Return ptr to entry structure, NULL if no ident found.
scan string intrst and note what is to be kept from feature table.
ndxrun flags if it is only an indexing run or whether it is necessary to
create storage for and read the sequence data  */

DB_ENTSTRCT *db_parsecore(FILE *df,
                          WRD_LUSTRCT *lclass,
                          DBLU_DBFMT fmt,
                          DB_FTYPELT *ufwrds,
                          int ndxrun);
/* scan df from current position to next sq terminator.  wlus is used for
key words.  Return ptr to entry structure, NULL if no ident found.
scan string intrst and note what is to be kept from feature table.
ndxrun flags if it is only an indexing run or whether it is necessary to
create storage for and read the sequence data */

DB_ENTSTRCT *db_parseflent(FILE *df,
                           WRD_LUSTRCT *lclass,
                           DBLU_DBFMT fmt,
                           char *ofintrst,
                           char *fqwants,
                           int ndxrun);
/* scan df from current position to next sq terminator.  wlus is used for
key words.  Return ptr to entry structure, NULL if no ident found.
scan string ofintrst and note what is to be kept from feature table.
ndxrun flags if it is only an indexing run or whether it is necessary to
create storage for and read the sequence
data
*/

DB_ENTSTRCT *db_parsfl4feats(FILE *df,
                             DBLU_DBFMT fmt,
                             DB_FTYPELT *ftwrds);
/* scan df from current position to next sq terminator.  Note features
in ftwrds and qualifiers from fqwrds.
Return ptr to entry structure, NULL if no ident found.*/

DBLU_LINECLASS db_netlclass(DBLU_LINECLASS curlt,
                            DBLU_LINECLASS prvlt,
                            DBLU_DBFMT fmt);
/* depending on fmt, return the current database file category for these
settings */

int db_lstparse(FILE *sf,
                WRD_LUSTRCT *lclass,
                DBLU_DBFMT fmt,
                DB_FTYPELT *unlmted,
                DB_FTYPELT **lmted,
                FILE *df);
/* scan sf from current position to next sq terminator.  wlus is used for
key words.  list lines to df according to contents of unlmted & lmtd lists.
return no of lines printed */

int db_readsqlen(DB_ENTSTRCT *es);
  /* return the actual read sequence length */

char *db_strdupn(char *str,
                 size_t sln);
/* malloc() sln+1 bytes, copy sln chars of str, null complete & return */

void dbt_putidstr2fl(FILE *ofl,
                     int fwidth,
                     int *cc,
                     char *idst,
                     int idln);
/* attempt to put idstr contents to ofl. Keep track of filewidth and cc */

void dbt_putstr2fl(FILE *ofl,
                   int fwidth,
                   int *cc,
                   char *idst);
/* attempt to put idstr contents to ofl. Keep track of filewidth and cc */

DBLU_DBFMT db_chkflfmt(FILE *fl);
  /* check for the database type of the file open on fl by checking each line
for Genbank or Swiss/EMBL id.  Determine latter by checking for NAtype on
id line.  Return DBFMT_unknown for indeterminate file. Rewind file. */

int db_offst2str(unsigned long offst,
                 char *str,
                 int sln);
/* attempt to encode offst into str in sln chars.  return 1 if succeeds */

unsigned long db_str2offst(char *str,
                           int sln);
/* scan str for offset value */

char *db_dbfmt2str(DBLU_DBFMT fmt);
  /* return ptr to string for fmt */

SQFL_SQTOS db_natype2sqtos(DB_NATYPE nat);
  /* remapping function: return a sqfl_fns.h sequence type respresentation
  of dbpars.h nat */ 

char *db_sens2str(DB_SENS sns);
  /* return a string value for sns */

char *db_ftqual2str(DBLU_FTQUAL fqul);
  /* return a string value for fqul */

DBP_MULTI_ELEMNT *dbp_appnd_mult_elemnt(DBP_MULTI_ELEMNT **melst,
                                        char *meid,
                                        DB_ENTSTRCT *entptr);
/* append a new element to *melst, return address of new element if useful.
meid is strduped, entptr is likely to be NULL at this stage */

DBP_MULTI_ELEMNT *dbp_melemnt4id(DBP_MULTI_ELEMNT *melst,
                                 char *queryid,
                                 int (*cmpfn)(const char *x1,
                                              const char *x2));
/* use (*cmpfn)() to locate the first element matching
queryid, return its address, else NULL */

void dbp_delmultelt(DBP_MULTI_ELEMNT *dep,
                    DBP_MULTI_ELEMNT **lstrt,
                    DBP_MULTI_ELEMNT **lend);
/* dep is an element of *lstrt: delete it from the list 
*lstrt..*lend,  lend can be NULL */

void dbp_killmultieltlist(DBP_MULTI_ELEMNT **lstrt,
                          DBP_MULTI_ELEMNT **lend);
/* iteratively remove all of *lstrt..*lend.  lend
can be NULL */

char *dbp_cpystrmodurlesc(char *srcstr,
                          char *dststr);
/* scan srcstr, looking for %xx where (xx is hexadecimal)
and generate a copy of the string in dststr with
escaped URL chars substituted.  return address of
dststr for convenience.  Invalid xx are simply 
copied verbatim. */

char *dbp_cpynreplurlescchars(char *ustr);
  /* malloc() ustr and return a version with
any URL-escaped chars replaced.  Strings
created by this function must be free()ed */

int dbp_fndnappgff3info(char *toks[],
                        int tmax,
                        char *leadstr,
                        DBLU_FTQUAL ftq,
                        DB_FEATSTRCT *ftp);
/* scan thru toks looking for a string headed with leadstr.
take the remainder, replacing any URL-escaped chars and append
with qualifier-type ftq to infolist of ftp.  return 1 if
succeeded */

int dbp_parse_gff3_ln(FILE *srcfl, 
                      DBP_MULTI_ELEMNT **curseq,
                      char *ln,
                      WRD_LUSTRCT *ftkw,
                      DB_FTYPELT *ftwrds,
                      DBP_MULTI_ELEMNT **elmntlst);
/* ln contains text: check for strlen() > 0;  if so, 
check leading '#' and return OK but do nothing if found.
otherwise break into ht separated tokens and check
if we've seen this.
Assume line is '\0' terminated.
if so, insert relevant information to that elemnt,
else append a new element and init it.
return 1 if something was parsed */

int dbp_parse_gtf_ln(FILE *srcfl, 
                     DBP_MULTI_ELEMNT **curseq,
                     char *ln,
                     WRD_LUSTRCT *ftkw,
                     DB_FTYPELT *ftwrds,
                     DBP_MULTI_ELEMNT **elmntlst,
                     DB_STR_ELT *attr);
/* ln contains text: check for strlen() > 0;  if so, 
check leading '#' and return OK but do nothing if found.
otherwise break into ht separated tokens and check
if we've seen this.
Assume line is '\0' terminated.
if so, insert relevant information to that elemnt,
else append a new element and init it.
if attr is non-NULL, then use it to select the
attribute values for the last field.
return 1 if something was parsed */

int dbp_findparentfeats(DB_ENTSTRCT *esp,
                        int (*strcmpfn)(const char *s1,
                                        const char *s2),
                        int lnkgene2chr);
/* using (*strcmp)(), scan all features of esp, looking for
any which can be linked as parents.  Note that the
chromosome line (FTKW_source_ftkw) is a parent of gene features
de facto: make this assertion if lnkgene2chr.
When a parentfeat is assigned, the parentname can
be free()ed and NULLed to avoid later complications.
Return number of assignments made.
Also set any feature mRNA links.
This process must be done separately,  by scanning for parent
 names since it doesn't seem
that gff3 files have any enforced requirement for parent
features to precede their children. */

int dbp_gtfgenes4entstrct(DB_ENTSTRCT *esptr,
                          WRD_LUSTRCT *idlut);
  /* scan ent structure esptr, looking
for exon groups which can be regarded as an mRNA.
reset the feature type for these and create a
corresponding gene entry.  Create other relevant
links for related CDS feature.
return number of mRNAs converted */

int dbp_bldgtfgenes(DBP_MULTI_ELEMNT *melstp);
  /* progressively scan melstp ent structures, looking
for exon groups which can be regarded as an mRNA.
reset the feature type for these and create a
corresponding gene entry.  Create other relevant
links for related CDS feature.
return number of mRNAs converted */

/* EOF */
