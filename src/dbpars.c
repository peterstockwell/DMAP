/* dbpars.c: seq database parsing functions in c */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
/* #include <sys/mode.h> */
#include <sys/types.h>
#include <sys/stat.h>
/* #include <gdbm.h> */

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqtrans.h"
#include "sqmat_fns.h"
#include "dbpars.h"
/* #include "dblu_fns.h" */

/* global variables: where reference to
a defined location is appropriate */

/* Dequote string used in deauotbiotypestr() */

char btypecopy[MAX_BIOTYPE_LENGTH+1];

int db_datoffst(DBLU_DBFMT fmt)
  /* return the start of normal data */
{
switch (fmt)
  {
  case DBFMT_embl:
  case DBFMT_swiss:
  case DBFMT_embl_old:
  case DBFMT_sqmonk:
    return(DBLU_MBLFLDSTRT - 1);
    break;
  case DBFMT_genbank:
  case DBFMT_unknown:
  default:
    return(DBLU_GBFLDSTRT - 1);
    break;
  }
}

DBLU_DBFMT db_chr2dbfmt(char c)
  /* return the database format corresponding to this character */
{
switch (toupper(c))
  {
  case 'E':
    return(DBFMT_embl);
    break;
  case 'G':
    return(DBFMT_genbank);
    break;
  case 'S':
    return(DBFMT_swiss);
    break;
  case 'O':
    return(DBFMT_embl_old);
    break;
  case 'F':
    return(DBFMT_gff3);
    break;
  case 'T':
    return(DBFMT_gtf);
    break;
  case 'Q':
    return(DBFMT_sqmonk);
    break;
  default:
    return(DBFMT_unknown);
    break;
  }
}

char db_dbfmt2chr(DBLU_DBFMT fmt)
  /* return a letter for database format fmt */
{
switch (fmt)
  {
  case DBFMT_embl:
    return('e');
    break;
  case DBFMT_genbank:
    return('g');
    break;
  case DBFMT_swiss:
    return('s');
    break;
  case DBFMT_embl_old:
    return('E');
    break;
  case DBFMT_gff3:
    return('F');
    break;
  case DBFMT_gtf:
    return('T');
    break;
  case DBFMT_sqmonk:
    return('Q');
    break;
  case DBFMT_unknown:
  default:
    return('u');
    break;
  }
}

int db_logerr(char *msg,
              int lno,
              char *lin)
  /* at present just write to stderr, return 0 to allow logical test */
{
fprintf(stderr,"Error: %s\n",msg);
if (lin != NULL)
  fprintf(stderr,"Line %d: %s\n",lno,lin);
return(0);
}

WRD_LUSTRCT *db_getclasstrct(DBLU_DBFMT dfmt)
  /* depending on dfmt, create and fill the datastructure for lineclass 
parsing, returning the address of the structure created.  Invalid dfmt
will return NULL.  malloc failures will die */
{
WRD_LUSTRCT *ws;

switch (dfmt) 
  {
  case DBFMT_swiss:
  case DBFMT_embl:
  case DBFMT_embl_old:
  case DBFMT_sqmonk:
    ws = wlu_getlustrct(1,DBLC_unk);
    wlu_addwrd(ws,"ID",DBLC_ident,NULL);
    wlu_addwrd(ws,"AC",DBLC_access,NULL);
    wlu_addwrd(ws,"NI",DBLC_naident,NULL);
    wlu_addwrd(ws,"DT",DBLC_datestamp,NULL);
    wlu_addwrd(ws,"DE",DBLC_description,NULL);
    wlu_addwrd(ws,"KW",DBLC_keyword,NULL);
    wlu_addwrd(ws,"OS",DBLC_source,NULL);
    wlu_addwrd(ws,"OC",DBLC_taxonomy,NULL);
    wlu_addwrd(ws,"RA",DBLC_author,NULL);
    wlu_addwrd(ws,"FT",DBLC_feature,NULL);
    wlu_addwrd(ws,"SQ",DBLC_sequence,NULL);
    wlu_addwrd(ws,"RT",DBLC_title,NULL);
    wlu_addwrd(ws,"RL",DBLC_journal,NULL);
    wlu_addwrd(ws,"RM",DBLC_remark,NULL);
    wlu_addwrd(ws,"XX",DBLC_ignore,NULL);
    wlu_addwrd(ws,"CC",DBLC_ignore,NULL);
    wlu_addwrd(ws,"FH",DBLC_ignore,NULL);
    wlu_addwrd(ws,"SV",DBLC_seqversn,NULL);
    wlu_addwrd(ws,"RN",DBLC_refnumbr,NULL);
    wlu_addwrd(ws,"RP",DBLC_refpositn,NULL);
    wlu_addwrd(ws,"RX",DBLC_xref,NULL);
    wlu_addwrd(ws,"DR",DBLC_dbxref,NULL);
    wlu_addwrd(ws,"//",DBLC_terminator,NULL);
    return(ws);
    break;
  case DBFMT_genbank:
    ws = wlu_getlustrct(1,DBLC_unk);
    wlu_addwrd(ws,"LOCUS",DBLC_ident,NULL);
    wlu_addwrd(ws,"ACCESSION",DBLC_access,NULL);
    wlu_addwrd(ws,"NID",DBLC_naident,NULL);
    wlu_addwrd(ws,"DATE",DBLC_datestamp,NULL);
    wlu_addwrd(ws,"DEFINITION",DBLC_description,NULL);
    wlu_addwrd(ws,"KEYWORDS",DBLC_keyword,NULL);
    wlu_addwrd(ws,"SOURCE",DBLC_source,NULL);
    wlu_addwrd(ws,"  ORGANISM",DBLC_taxonomy,NULL);
    wlu_addwrd(ws,"Reference",DBLC_reference,NULL);
    wlu_addwrd(ws,"  AUTHORS",DBLC_author,NULL);
    wlu_addwrd(ws,"FEATURES",DBLC_feature,NULL);
    wlu_addwrd(ws,"BASE COUNT",DBLC_basecount,NULL);
    wlu_addwrd(ws,"ORIGIN",DBLC_sequence,NULL);
    wlu_addwrd(ws,"  TITLE",DBLC_title,NULL);
    wlu_addwrd(ws,"  JOURNAL",DBLC_journal,NULL);
    wlu_addwrd(ws,"  REMARK",DBLC_remark,NULL);
    wlu_addwrd(ws,"VERSION",DBLC_seqversn,NULL);
    wlu_addwrd(ws,"//",DBLC_terminator,NULL);
    return(ws);
    break;
  case DBFMT_unknown:
  case DBFMT_gff3:
  case DBFMT_gtf:
  default:
    return(NULL);
    break;
  }
}
  
char *db_ftkw2str(DBLU_FTKW kw)
  /* return a string corresponding to the kw given */
{
switch (kw)
  {
  case FTKW_allele:
    return("allele");
    break;
  case FTKW_attenuator:
    return("attenuator");
    break;
  case FTKW_C_region:
    return("C_region");
    break;
  case FTKW_CAAT_signal:
    return("CAAT_signal");
    break;
  case FTKW_CDS:
    return("CDS");
    break;
  case FTKW_cellular:
    return("cellular");
    break;
  case FTKW_conflict:
    return("conflict");
    break;
  case FTKW_D_loop:
    return("D-loop");
    break;
  case FTKW_D_region:
    return("D_region");
    break;
  case FTKW_enhancer:
    return("enhancer");
    break;
  case FTKW_exon:
    return("exon");
    break;
  case FTKW_FT_gene:
    return("gene");
    break;
  case FTKW_GC_signal:
    return("GC_signal");
    break;
  case FTKW_iDNA:
    return("iDNA");
    break;
  case FTKW_insertion_seq:
    return("insertion_seq");
    break;
  case FTKW_intron:
    return("intron");
    break;
  case FTKW_J_region:
    return("J_region");
    break;
  case FTKW_LTR:
    return("LTR");
    break;
  case FTKW_lncRNA:
    return("lncRNA");
    break;
  case FTKW_mat_peptide:
    return("mat_peptide");
    break;
  case FTKW_misc_binding:
    return("misc_binding");
    break;
  case FTKW_misc_difference:
    return("misc_difference");
    break;
  case FTKW_misc_feature:
    return("misc_feature");
    break;
  case FTKW_misc_recomb:
    return("misc_recomb");
    break;
  case FTKW_misc_RNA:
    return("misc_RNA");
    break;
  case FTKW_misc_signal:
    return("misc_signal");
    break;
  case FTKW_misc_structure:
    return("misc_structure");
    break;
  case FTKW_modified_base:
    return("modified_base");
    break;
  case FTKW_mRNA:
    return("mRNA");
    break;
  case FTKW_mutation:
    return("mutation");
    break;
  case FTKW_N_region:
    return("N_region");
    break;
  case FTKW_old_sequence:
    return("old_sequence");
    break;
  case FTKW_polyA_signal:
    return("polyA_signal");
    break;
  case FTKW_polyA_site:
    return("polyA_site");
    break;
  case FTKW_precursor_RNA:
    return("precursor_RNA");
    break;
  case FTKW_prim_transcript:
    return("prim_transcript");
    break;
  case FTKW_primer:
    return("primer");
    break;
  case FTKW_primer_bind:
    return("primer_bind");
    break;
  case FTKW_promoter:
    return("promoter");
    break;
  case FTKW_protein_bind:
    return("protein_bind");
    break;
  case FTKW_provirus:
    return("provirus");
    break;
  case FTKW_RBS:
    return("RBS");
    break;
  case FTKW_rep_origin:
    return("rep_origin");
    break;
  case FTKW_repeat_region:
    return("repeat_region");
    break;
  case FTKW_repeat_unit:
    return("repeat_unit");
    break;
  case FTKW_rRNA:
    return("rRNA");
    break;
  case FTKW_S_region:
    return("S_region");
    break;
  case FTKW_satellite:
    return("satellite");
    break;
  case FTKW_scRNA:
    return("scRNA");
    break;
  case FTKW_sig_peptide:
    return("sig_peptide");
    break;
  case FTKW_snRNA:
    return("snRNA");
    break;
  case FTKW_stem_loop:
    return("stem_loop");
    break;
  case FTKW_STS:
    return("STS");
    break;
  case FTKW_TATA_signal:
    return("TATA_signal");
    break;
  case FTKW_terminator_sq:
    return("terminator");
    break;
  case FTKW_transit_peptide:
    return("transit_peptide");
    break;
  case FTKW_transposon:
    return("transposon");
    break;
  case FTKW_tRNA:
    return("tRNA");
    break;
  case FTKW_unsure:
    return("unsure");
    break;
  case FTKW_V_region:
    return("V_region");
    break;
  case FTKW_variation:
    return("variation");
    break;
  case FTKW_virion:
    return("virion");
    break;
  case FTKW_hyphen:
    return("-");
    break;
  case FTKW_m10_signal:
    return("-10_signal");
    break;
  case FTKW_m35_signal:
    return("-35_signal");
    break;
  case FTKW_clip3p:
    return("3'clip");
    break;
  case FTKW_UTR3p:
    return("3'UTR");
    break;
  case FTKW_clip5p:
    return("5'clip");
    break;
  case FTKW_UTR5p:
    return("5'UTR");
    break;
  case FTKW_source_ftkw:
    return("source");
    break;
  case FTKW_startcodon:
    return("start_codon");
    break;
  case FTKW_unknown:
    return("unknown_ftkw");
    break;
  case FTKW_ft_error:
    return("FT_Error");
    break;
  case FTKW_SWP_act_site:
    return("ACT_SITE");
    break;
  case FTKW_SWP_binding:
    return("BINDING");
    break;
  case FTKW_SWP_carbohyd:
    return("CARBOHYD");
    break;
  case FTKW_SWP_ca_bind:
    return("CA_BIND");
    break;
  case FTKW_SWP_chain:
    return("CHAIN");
    break;
  case FTKW_SWP_conflict:
    return("CONFLICT");
    break;
  case FTKW_SWP_disulfid:
    return("DISULFID");
    break;
  case FTKW_SWP_dna_bind:
    return("DNA_BIND");
    break;
  case FTKW_SWP_domain:
    return("DOMAIN");
    break;
  case FTKW_SWP_helix:
    return("HELIX");
    break;
  case FTKW_SWP_init_met:
    return("INIT_MET");
    break;
  case FTKW_SWP_lipid:
    return("LIPID");
    break;
  case FTKW_SWP_metal:
    return("METAL");
    break;
  case FTKW_SWP_mod_res:
    return("MOD_RES");
    break;
  case FTKW_SWP_mutagen:
    return("MUTAGEN");
    break;
  case FTKW_SWP_non_cons:
    return("NON_CONS");
    break;
  case FTKW_SWP_non_ter:
    return("NON_TER");
    break;
  case FTKW_SWP_np_bind:
    return("NP_BIND");
    break;
  case FTKW_SWP_peptide:
    return("PEPTIDE");
    break;
  case FTKW_SWP_propep:
    return("PROPEP");
    break;
  case FTKW_SWP_repeat:
    return("REPEAT");
    break;
  case FTKW_SWP_signal:
    return("SIGNAL");
    break;
  case FTKW_SWP_similar:
    return("SIMILAR");
    break;
  case FTKW_SWP_site:
    return("SITE");
    break;
  case FTKW_SWP_strand:
    return("STRAND");
    break;
  case FTKW_SWP_thioeth:
    return("THIOETH");
    break;
  case FTKW_SWP_thiolest:
    return("THIOLEST");
    break;
  case FTKW_SWP_transit:
    return("TRANSIT");
    break;
  case FTKW_SWP_transmem:
    return("TRANSMEM");
    break;
  case FTKW_SWP_turn:
    return("TURN");
    break;
  case FTKW_SWP_unsure:
    return("UNSURE");
    break;
  case FTKW_SWP_variant:
    return("VARIANT");
    break;
  case FTKW_SWP_varsplic:
    return("VARSPLIC");
    break;
  case FTKW_SWP_zn_fing:
    return("ZN_FING");
    break;
  default:
    return("????");
    break;
  }
}

char *db_ftkw2text(DBLU_FTKW kw)
  /* return a long string corresponding to the kw given */
{
switch (kw)
  {
  case FTKW_allele:
    return("Related strain contains alternative gene form");
    break;
  case FTKW_attenuator:
    return("Sequence related to transcription termination");
    break;
  case FTKW_C_region:
    return("Span of the C immunological feature");
    break;
  case FTKW_CAAT_signal:
    return("`CAAT box' in eukaryotic promoters");
    break;
  case FTKW_CDS:
    return("Sequence coding for amino acids in protein");
    break;
  case FTKW_cellular:
    return("Region of cellular DNA");
    break;
  case FTKW_conflict:
    return("Independent determinations differ");
    break;
  case FTKW_D_loop:
    return("Displacement loop");
    break;
  case FTKW_D_region:
    return("Span of the D immunological feature");
    break;
  case FTKW_enhancer:
    return("Cis-acting enhancer of promoter function");
    break;
  case FTKW_exon:
    return("Region that codes for part of spliced mRNA");
    break;
  case FTKW_FT_gene:
    return("gene");
    break;
  case FTKW_GC_signal:
    return("`GC box' in eukaryotic promoters");
    break;
  case FTKW_iDNA:
    return("Intervening DNA eliminated by recombination");
    break;
  case FTKW_insertion_seq:
    return("Insertion sequence (IS), a small transposon");
    break;
  case FTKW_intron:
    return("Transcribed region excised by mRNA splicing");
    break;
  case FTKW_J_region:
    return("Span of the J immunological feature");
    break;
  case FTKW_LTR:
    return("Long Terminal Repeat");
    break;
  case FTKW_lncRNA:
    return("long noncoding RNA");
    break;
  case FTKW_mat_peptide:
    return("Mature peptide coding region (does not include stop codon)");
    break;
  case FTKW_misc_binding:
    return("Miscellaneous binding site");
    break;
  case FTKW_misc_difference:
    return("Miscellaneous difference feature");
    break;
  case FTKW_misc_feature:
    return("Region of biological significance not otherwise described");
    break;
  case FTKW_misc_recomb:
    return("Miscellaneous recombination feature");
    break;
  case FTKW_misc_RNA:
    return("Miscellaneous transcript feature not defined by other RNA keys");
    break;
  case FTKW_misc_signal:
    return("Miscellaneous signal");
    break;
  case FTKW_misc_structure:
    return("Miscellaneous DNA or RNA structure");
    break;
  case FTKW_modified_base:
    return("The indicated base is a modified nucleotide");
    break;
  case FTKW_mRNA:
    return("Messenger RNA");
    break;
  case FTKW_mutation:
    return("A mutation alters the sequence here");
    break;
  case FTKW_N_region:
    return("Span of the N immunological feature");
    break;
  case FTKW_old_sequence:
    return("Presented sequence revises a previous version");
    break;
  case FTKW_polyA_signal:
    return("Signal for cleavage & polyadenylation");
    break;
  case FTKW_polyA_site:
    return("Site at which polyadenine is added to mRNA");
    break;
  case FTKW_precursor_RNA:
    return("Any RNA species that is not yet the mature RNA product");
    break;
  case FTKW_prim_transcript:
    return("Primary (unprocessed) transcript");
    break;
  case FTKW_primer:
    return("Primer binding region used with PCR");
    break;
  case FTKW_primer_bind:
    return("Non-covalent primer binding site");
    break;
  case FTKW_promoter:
    return("A region involved in transcription initiation");
    break;
  case FTKW_protein_bind:
    return("Non-covalent protein binding site on DNA or RNA");
    break;
  case FTKW_provirus:
    return("Proviral sequence");
    break;
  case FTKW_RBS:
    return("Ribosome binding site");
    break;
  case FTKW_rep_origin:
    return("Replication origin for duplex DNA");
    break;
  case FTKW_repeat_region:
    return("Sequence containing repeated subsequences");
    break;
  case FTKW_repeat_unit:
    return("One repeated unit of a repeat_region");
    break;
  case FTKW_rRNA:
    return("Ribosomal RNA");
    break;
  case FTKW_S_region:
    return("Span of the S immunological feature");
    break;
  case FTKW_satellite:
    return("Satellite repeated sequence");
    break;
  case FTKW_scRNA:
    return("Small cytoplasmic RNA");
    break;
  case FTKW_sig_peptide:
    return("Signal peptide coding region");
    break;
  case FTKW_snRNA:
    return("Small nuclear RNA");
    break;
  case FTKW_stem_loop:
    return("Hair-pin loop structure in DNA or RNA");
    break;
  case FTKW_STS:
    return("Sequence Tagged Site");
    break;
  case FTKW_TATA_signal:
    return("`TATA box' in eukaryotic promoters");
    break;
  case FTKW_terminator_sq:
    return("Sequence causing transcription termination");
    break;
  case FTKW_transit_peptide:
    return("Transit peptide coding region");
    break;
  case FTKW_transposon:
    return("Transposable element (TN)");
    break;
  case FTKW_tRNA:
    return("Transfer RNA");
    break;
  case FTKW_unsure:
    return("Authors are unsure about the sequence in this region");
    break;
  case FTKW_V_region:
    return("Span of the V immunological feature");
    break;
  case FTKW_variation:
    return("A related population contains stable mutation");
    break;
  case FTKW_virion:
    return("Virion (encapsidated) viral sequence");
    break;
  case FTKW_hyphen:
    return("hyphen");
    break;
  case FTKW_m10_signal:
    return("`Pribnow box' in prokaryotic promoters");
    break;
  case FTKW_m35_signal:
    return("`-35 box' in prokaryotic promoters");
    break;
  case FTKW_clip3p:
    return("3'-most region of a precursor transcript removed in processing");
    break;
  case FTKW_UTR3p:
    return("3' untranslated region (trailer)");
    break;
  case FTKW_clip5p:
    return("5'-most region of a precursor transcript removed in processing");
    break;
  case FTKW_UTR5p:
    return("5' untranslated region (leader)");
    break;
  case FTKW_source_ftkw:
    return("sequence source");
    break;
  case FTKW_unknown:
    return("Unknown Feature Keyword");
    break;
  case FTKW_ft_error:
    return("FT_Error");
    break;
  case FTKW_SWP_act_site:
    return("ACT_SITE");
    break;
  case FTKW_SWP_binding:
    return("BINDING");
    break;
  case FTKW_SWP_carbohyd:
    return("CARBOHYD");
    break;
  case FTKW_SWP_ca_bind:
    return("CA_BIND");
    break;
  case FTKW_SWP_chain:
    return("CHAIN");
    break;
  case FTKW_SWP_conflict:
    return("CONFLICT");
    break;
  case FTKW_SWP_disulfid:
    return("DISULFID");
    break;
  case FTKW_SWP_dna_bind:
    return("DNA_BIND");
    break;
  case FTKW_SWP_domain:
    return("DOMAIN");
    break;
  case FTKW_SWP_helix:
    return("HELIX");
    break;
  case FTKW_SWP_init_met:
    return("INIT_MET");
    break;
  case FTKW_SWP_lipid:
    return("LIPID");
    break;
  case FTKW_SWP_metal:
    return("METAL");
    break;
  case FTKW_SWP_mod_res:
    return("MOD_RES");
    break;
  case FTKW_SWP_mutagen:
    return("MUTAGEN");
    break;
  case FTKW_SWP_non_cons:
    return("NON_CONS");
    break;
  case FTKW_SWP_non_ter:
    return("NON_TER");
    break;
  case FTKW_SWP_np_bind:
    return("NP_BIND");
    break;
  case FTKW_SWP_peptide:
    return("PEPTIDE");
    break;
  case FTKW_SWP_propep:
    return("PROPEP");
    break;
  case FTKW_SWP_repeat:
    return("REPEAT");
    break;
  case FTKW_SWP_signal:
    return("SIGNAL");
    break;
  case FTKW_SWP_similar:
    return("SIMILAR");
    break;
  case FTKW_SWP_site:
    return("SITE");
    break;
  case FTKW_SWP_strand:
    return("STRAND");
    break;
  case FTKW_SWP_thioeth:
    return("THIOETH");
    break;
  case FTKW_SWP_thiolest:
    return("THIOLEST");
    break;
  case FTKW_SWP_transit:
    return("TRANSIT");
    break;
  case FTKW_SWP_transmem:
    return("TRANSMEM");
    break;
  case FTKW_SWP_turn:
    return("TURN");
    break;
  case FTKW_SWP_unsure:
    return("UNSURE");
    break;
  case FTKW_SWP_variant:
    return("VARIANT");
    break;
  case FTKW_SWP_varsplic:
    return("VARSPLIC");
    break;
  case FTKW_SWP_zn_fing:
    return("ZN_FING");
    break;
  default:
    return("????");
    break;
  }
}

char *db_ftqu2str(DBLU_FTQUAL fqual)
  /* return a string for fqual */
{
switch (fqual)
  {
  case FTQU_anticodon:
    return("anticodon");
    break;
  case FTQU_bound_moiety:
    return("bound_moiety");
    break;
  case FTQU_citation:
    return("citation");
    break;
  case FTQU_codon:
    return("codon");
    break;
  case FTQU_codon_start:
    return("codon_start");
    break;
  case FTQU_cons_splice:
    return("cons_splice");
    break;
  case FTQU_db_xref:
    return("db_xref");
    break;
  case FTQU_direction:
    return("direction");
    break;
  case FTQU_EC_number:
    return("EC_number");
    break;
  case FTQU_evidence:
    return("evidence");
    break;
  case FTQU_frequency:
    return("frequency");
    break;
  case FTQU_function_qu:
    return("function_qu");
    break;
  case FTQU_gene:
    return("gene");
    break;
  case FTQU_label_qu:
    return("label_qu");
    break;
  case FTQU_mod_base:
    return("mod_base");
    break;
  case FTQU_note:
    return("note");
    break;
  case FTQU_number:
    return("number");
    break;
  case FTQU_organism:
    return("organism");
    break;
  case FTQU_partial:
    return("partial");
    break;
  case FTQU_phenotype:
    return("phenotype");
    break;
  case FTQU_product:
    return("product");
    break;
  case FTQU_protid:
    return("protid");
    break;
  case FTQU_pseudo:
    return("pseudo");
    break;
  case FTQU_rpt_family:
    return("rpt_family");
    break;
  case FTQU_rpt_type:
    return("rpt_type");
    break;
  case FTQU_rpt_unit:
    return("rpt_unit");
    break;
  case FTQU_standard_name:
    return("standard_name");
    break;
  case FTQU_transl_except:
    return("transl_except");
    break;
  case FTQU_transl_table:
    return("transl_table");
    break;
  case FTQU_type_qu:
    return("type_qu");
    break;
  case FTQU_usedin:
    return("usedin");
    break;
  case FTQU_locus_tag:
    return("locus_tag");
    break;
  case FTQU_GO_info:
    return("GO_info");
    break;
  case FTQU_name:
    return("name");
    break;
  case FTQU_biotype:
    return("biotype");
    break;
  case FTQU_phosphothre:
    return("phosphothre");
    break;
  case FTQU_phosphoser:
    return("phsophoser");
    break;
  FTQU_unknown:
  default:
    return("Unknown");
    break;
  }
}

WRD_LUSTRCT *db_getkwstrct(DBLU_DBFMT dfmt)
  /* create and fill datastructure for feature table key word parsing:
    return("phosphoser

  case FTQU_unknown:
  default:
    return("Unknown");
    break;
  }
}

WRD_LUSTRCT *db_getkwstrct(DBLU_DBFMT dfmt)
  /* create and fill datastructure for feature table key word parsing");
break;
depending on dfmt.  return pointer to that structure, NULL for invalid
dfmt  */
{
WRD_LUSTRCT *ws;
DBLU_FTKW wp;

switch (dfmt)
  {
  case DBFMT_genbank:
  case DBFMT_embl:
  case DBFMT_embl_old:
    ws = wlu_getlustrct(WLU_CASEIND,FTKW_unknown);
    for (wp = FTKW_allele; wp <= FTKW_source_ftkw; wp++)
      wlu_addwrd(ws,db_ftkw2str(wp),(int) wp,NULL);
    return(ws);
    break;
  case DBFMT_swiss:
    ws = wlu_getlustrct(WLU_CASEIND,FTKW_unknown);
    for (wp = FTKW_SWP_act_site; wp <= FTKW_SWP_zn_fing; wp++)
      wlu_addwrd(ws,db_ftkw2str(wp),wp,NULL);
    return(ws);
    break;
  case DBFMT_sqmonk:
    ws = db_getkwstrct(DBFMT_embl);
    wlu_addwrd(ws,"CpG",(int) FTKW_CpG_is,NULL);
    wlu_addwrd(ws,"TSS",(int) FTKW_TSS,NULL);
    return(ws);
    break;
  case DBFMT_gff3:
    ws = wlu_getlustrct(WLU_CASEIND,FTKW_unknown);
    wlu_addwrd(ws,"chromosome",(int) FTKW_source_ftkw,NULL);
    wlu_addwrd(ws,"interior_coding_exon",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"satellite_DNA",(int) FTKW_satellite,NULL);
    wlu_addwrd(ws,"transposable_element",(int) FTKW_transposon,NULL);
    wlu_addwrd(ws,"polypeptide",(int) FTKW_mat_peptide,NULL);
    wlu_addwrd(ws,"sequence_variant_obs",(int) FTKW_variation,NULL);
    wlu_addwrd(ws,"sequence_feature",(int) FTKW_misc_feature,NULL);
    wlu_addwrd(ws,"primer",(int) FTKW_primer,NULL);
    wlu_addwrd(ws,"proviral_region",(int) FTKW_provirus,NULL);
    wlu_addwrd(ws,"protein_coding_primary_transcript",(int) FTKW_FT_gene,NULL);
    wlu_addwrd(ws,"ribosome_entry_site",(int) FTKW_RBS,NULL);
    wlu_addwrd(ws,"attenuator",(int) FTKW_attenuator,NULL);
    wlu_addwrd(ws,"terminator",(int) FTKW_terminator_sq,NULL);
    wlu_addwrd(ws,"exon",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"methylated_C",(int) FTKW_modified_base,NULL);
    wlu_addwrd(ws,"methylated_A",(int) FTKW_modified_base,NULL);
    wlu_addwrd(ws,"enhancer",(int) FTKW_enhancer,NULL);
    wlu_addwrd(ws,"promoter",(int) FTKW_promoter,NULL);
    wlu_addwrd(ws,"primary_transcript",(int) FTKW_prim_transcript,NULL);
    wlu_addwrd(ws,"intron",(int) FTKW_intron,NULL);
    wlu_addwrd(ws,"coding_exon",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"five_prime_coding_exon_coding_region",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"three_prime_coding_exon_coding_region",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"noncoding_exon",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"five_prime_coding_exon",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"five_prime_UTR",(int) FTKW_UTR5p,NULL);
    wlu_addwrd(ws,"three_prime_UTR",(int) FTKW_UTR3p,NULL);
    wlu_addwrd(ws,"rRNA_primary_transcript",(int) FTKW_precursor_RNA,NULL);
    wlu_addwrd(ws,"mature_transcript",(int) FTKW_mRNA,NULL);
    wlu_addwrd(ws,"mRNA",(int) FTKW_mRNA,NULL);
    wlu_addwrd(ws,"ORF",(int) FTKW_CDS,NULL);
    wlu_addwrd(ws,"rRNA",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"tRNA",(int) FTKW_tRNA,NULL);
    wlu_addwrd(ws,"snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"snoRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"microsatellite",(int) FTKW_satellite,NULL);
    wlu_addwrd(ws,"inverted_repeat",(int) FTKW_stem_loop,NULL);
    wlu_addwrd(ws,"modified_base",(int) FTKW_modified_base,NULL);
    wlu_addwrd(ws,"methylated_base_feature",(int) FTKW_modified_base,NULL);
    wlu_addwrd(ws,"direct_repeat",(int) FTKW_repeat_region,NULL);
    wlu_addwrd(ws,"CDS",(int) FTKW_CDS,NULL);
    wlu_addwrd(ws,"rRNA_large_subunit_primary_transcript",(int) FTKW_prim_transcript,NULL);
    wlu_addwrd(ws,"conserved_region",(int) FTKW_C_region,NULL);
    wlu_addwrd(ws,"STS",(int) FTKW_STS,NULL);
    wlu_addwrd(ws,"origin_of_replication",(int) FTKW_rep_origin,NULL);
    wlu_addwrd(ws,"insertion_site",(int) FTKW_insertion_seq,NULL);
    wlu_addwrd(ws,"enzymatic_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"ribozyme",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"rRNA_5_8S",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"hammerhead_ribozyme",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"RNase_MRP_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"RNase_P_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"telomerase_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"U1_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U2_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U4_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U4atac_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U5_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U6_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U6atac_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U11_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U12_snRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"U14_snoRNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"vault_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"Y_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"rRNA_18S",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"protein_binding_site",(int) FTKW_protein_bind,NULL);
    wlu_addwrd(ws,"sequence_difference",(int) FTKW_conflict,NULL);
    wlu_addwrd(ws,"signal_peptide",(int) FTKW_sig_peptide,NULL);
    wlu_addwrd(ws,"mature_protein_region",(int) FTKW_mat_peptide,NULL);
    wlu_addwrd(ws,"rasiRNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"nc_primary_transcript",(int) FTKW_prim_transcript,NULL);
    wlu_addwrd(ws,"three_prime_coding_exon_noncoding_region",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"five_prime_coding_exon_noncoding_region",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"polyA_signal_sequence",(int) FTKW_polyA_signal,NULL);
    wlu_addwrd(ws,"polyA_site",(int) FTKW_polyA_site,NULL);
    wlu_addwrd(ws,"group_I_intron",(int) FTKW_intron,NULL);
    wlu_addwrd(ws,"autocatalytically_spliced_intron",(int) FTKW_intron,NULL);
    wlu_addwrd(ws,"SRP_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"C_D_box_snoRNA",(int) FTKW_snRNA,NULL);
    wlu_addwrd(ws,"guide_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"group_II_intron",(int) FTKW_intron,NULL);
    wlu_addwrd(ws,"minisatellite",(int) FTKW_satellite,NULL);
    wlu_addwrd(ws,"chromosomal_structural_element",(int) FTKW_misc_structure,NULL);
    wlu_addwrd(ws,"antisense_RNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"antisense_primary_transcript",(int) FTKW_prim_transcript,NULL);
    wlu_addwrd(ws,"siRNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"stRNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"small_subunit_rRNA",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"large_subunit_rRNA",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"rRNA_5S",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"rRNA_28S",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"ncRNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"repeat_region",(int) FTKW_repeat_region,NULL);
    wlu_addwrd(ws,"dispersed_repeat",(int) FTKW_repeat_unit,NULL);
    wlu_addwrd(ws,"spliceosomal_intron",(int) FTKW_intron,NULL);
    wlu_addwrd(ws,"insertion",(int) FTKW_insertion_seq,NULL);
    wlu_addwrd(ws,"transcript",(int) FTKW_prim_transcript,NULL);
    wlu_addwrd(ws,"SNP",(int) FTKW_variation,NULL);
    wlu_addwrd(ws,"possible_base_call_error",(int) FTKW_unsure,NULL);
    wlu_addwrd(ws,"possible_assembly_error",(int) FTKW_unsure,NULL);
    wlu_addwrd(ws,"gene",(int) FTKW_FT_gene,NULL);
    wlu_addwrd(ws,"tandem_repeat",(int) FTKW_repeat_region,NULL);
    wlu_addwrd(ws,"nucleotide_motif",(int) FTKW_misc_structure,NULL);
    wlu_addwrd(ws,"RNA_motif",(int) FTKW_misc_feature,NULL);
    wlu_addwrd(ws,"transit_peptide",(int) FTKW_transit_peptide,NULL);
    wlu_addwrd(ws,"pseudogenic_rRNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"pseudogenic_tRNA",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"primary_transcript_region",(int) FTKW_prim_transcript,NULL);
    wlu_addwrd(ws,"mRNA_region",(int) FTKW_misc_RNA,NULL);
    wlu_addwrd(ws,"polypeptide_region",(int) FTKW_mat_peptide,NULL);
    wlu_addwrd(ws,"CDS_region",(int) FTKW_CDS,NULL);
    wlu_addwrd(ws,"exon_region",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"rRNA_16S",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"rRNA_23S",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"rRNA_25S",(int) FTKW_rRNA,NULL);
    wlu_addwrd(ws,"copy_number_variation",(int) FTKW_variation,NULL);
    wlu_addwrd(ws,"sequence_alteration",(int) FTKW_variation,NULL);
    wlu_addwrd(ws,"immature_peptide_region",(int) FTKW_FT_gene,NULL);
    wlu_addwrd(ws,"noncoding_region_of_exon",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"coding_region_of_exon",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"replicon",(int) FTKW_rep_origin,NULL);
    wlu_addwrd(ws,"experimental_feature",(int) FTKW_misc_feature,NULL);
    wlu_addwrd(ws,"biological_region",(int) FTKW_misc_feature,NULL);
    wlu_addwrd(ws,"SNV",(int) FTKW_variation,NULL);
    wlu_addwrd(ws,"peptide_localization_signal",(int) FTKW_sig_peptide,NULL);
    wlu_addwrd(ws,"nucleotide_to_protein_binding_site",(int) FTKW_protein_bind,NULL);
    wlu_addwrd(ws,"sequence_motif",(int) FTKW_misc_feature,NULL);
    wlu_addwrd(ws,"epigenetically_modified_region",(int) FTKW_modified_base,NULL);
    wlu_addwrd(ws,"substitution",(int) FTKW_variation,NULL);
    wlu_addwrd(ws,"complex_substitution",(int) FTKW_conflict,NULL);
    wlu_addwrd(ws,"point_mutation",(int) FTKW_mutation,NULL);
    wlu_addwrd(ws,"edited_to",(int) FTKW_conflict,NULL);
    wlu_addwrd(ws,"exemplar_of",(int) FTKW_misc_feature,NULL);
    return(ws);
    break;
  case DBFMT_gtf:
    ws = wlu_getlustrct(WLU_CASEIND,FTKW_unknown);
    wlu_addwrd(ws,"CDS",(int) FTKW_CDS,NULL);
    wlu_addwrd(ws,"start_codon",(int) FTKW_startcodon,NULL);
    wlu_addwrd(ws,"stop_codon",(int) FTKW_terminator_sq,NULL);
    wlu_addwrd(ws,"exon",(int) FTKW_exon,NULL);
    wlu_addwrd(ws,"intron_CNS",(int) FTKW_intron,NULL);
    wlu_addwrd(ws,"intron",(int) FTKW_intron,NULL);
    wlu_addwrd(ws,"5UTR",(int) FTKW_UTR5p,NULL);
    wlu_addwrd(ws,"3UTR",(int) FTKW_UTR3p,NULL);
    wlu_addwrd(ws,"CNS",(int) FTKW_CNS,NULL);
    wlu_addwrd(ws,"inter",(int) FTKW_intergenic,NULL);
    wlu_addwrd(ws,"transcript",(int) FTKW_prim_transcript,NULL);
    wlu_addwrd(ws,"gene",(int) FTKW_FT_gene,NULL);
    wlu_addwrd(ws,"lncRNA",(int) FTKW_lncRNA,NULL);
    return(ws);
    break;
  case DBFMT_unknown:
  default:
    return(NULL);
    break;
  }
}

WRD_LUSTRCT *db_getftlocstrct(void)
  /* create and fill the datastructure for location word 
parsing, returning the address of the structure created. */
{
WRD_LUSTRCT *ws;

ws = wlu_getlustrct(1,FLOC_unknown);
wlu_addwrd(ws,"complement",FLOC_complement,NULL);
wlu_addwrd(ws,"join",FLOC_join,NULL);
wlu_addwrd(ws,"order",FLOC_order,NULL);
wlu_addwrd(ws,"group",FLOC_group,NULL);
wlu_addwrd(ws,"one_of",FLOC_one_of,NULL);
wlu_addwrd(ws,"replace",FLOC_replace,NULL);
wlu_addwrd(ws,"..",FLOC_drange,NULL);
wlu_addwrd(ws,".",FLOC_irange,NULL);
wlu_addwrd(ws,"^",FLOC_btween,NULL);
return(ws);
}

WRD_LUSTRCT *db_getftqualstrct(void)
  /* create and fill the datastructure for location qualifier word 
parsing, returning the address of the structure created. */
{
WRD_LUSTRCT *ws;

ws = wlu_getlustrct(1,FTQU_unknown);
wlu_addwrd(ws,"anticodon",FTQU_anticodon,NULL);
wlu_addwrd(ws,"bound_moiety",FTQU_bound_moiety,NULL);
wlu_addwrd(ws,"citation",FTQU_citation,NULL);
wlu_addwrd(ws,"codon",FTQU_codon,NULL);
wlu_addwrd(ws,"codon_start",FTQU_codon_start,NULL);
wlu_addwrd(ws,"cons_splice",FTQU_cons_splice,NULL);
wlu_addwrd(ws,"direction",FTQU_direction,NULL);
wlu_addwrd(ws,"EC_number",FTQU_EC_number,NULL);
wlu_addwrd(ws,"evidence",FTQU_evidence,NULL);
wlu_addwrd(ws,"frequency",FTQU_frequency,NULL);
wlu_addwrd(ws,"function",FTQU_function_qu,NULL);
wlu_addwrd(ws,"gene",FTQU_gene,NULL);
wlu_addwrd(ws,"label",FTQU_label_qu,NULL);
wlu_addwrd(ws,"mod_base",FTQU_mod_base,NULL);
wlu_addwrd(ws,"note",FTQU_note,NULL);
wlu_addwrd(ws,"number",FTQU_number,NULL);
wlu_addwrd(ws,"organism",FTQU_organism,NULL);
wlu_addwrd(ws,"partial",FTQU_partial,NULL);
wlu_addwrd(ws,"phenotype",FTQU_phenotype,NULL);
wlu_addwrd(ws,"product",FTQU_product,NULL);
wlu_addwrd(ws,"protein_id",FTQU_protid,NULL);
wlu_addwrd(ws,"pseudo",FTQU_pseudo,NULL);
wlu_addwrd(ws,"rpt_family",FTQU_rpt_family,NULL);
wlu_addwrd(ws,"rpt_type",FTQU_rpt_type,NULL);
wlu_addwrd(ws,"rpt_unit",FTQU_rpt_unit,NULL);
wlu_addwrd(ws,"standard_name",FTQU_standard_name,NULL);
wlu_addwrd(ws,"transl_except",FTQU_transl_except,NULL);
wlu_addwrd(ws,"type",FTQU_type_qu,NULL);
wlu_addwrd(ws,"usedin",FTQU_usedin,NULL);
wlu_addwrd(ws,"transl_table",FTQU_transl_table,NULL);
wlu_addwrd(ws,"db_xref",FTQU_db_xref,NULL);
wlu_addwrd(ws,"locus_tag",FTQU_locus_tag,NULL);
wlu_addwrd(ws,"name",FTQU_name,NULL);
wlu_addwrd(ws,"biotype",FTQU_biotype,NULL);
wlu_addwrd(ws,"Phosphothreonine",FTQU_phosphothre,NULL);
wlu_addwrd(ws,"Phosphoserine",FTQU_phosphoser,NULL);
return(ws);
}

WRD_LUSTRCT *db_getnatypestrct(void)
  /* create and fill the datastructure for NA_TYPE words, return
 the address of the structure created. */
{
WRD_LUSTRCT *ws;

ws = wlu_getlustrct(1,DBNA_unk);
wlu_addwrd(ws,"DNA",DBNA_dna,NULL);
wlu_addwrd(ws,"RNA",DBNA_rna,NULL);
wlu_addwrd(ws,"tRNA",DBNA_trna,NULL);
wlu_addwrd(ws,"rRNA",DBNA_rrna,NULL);
wlu_addwrd(ws,"mRNA",DBNA_mrna,NULL);
wlu_addwrd(ws,"uRNA",DBNA_urna,NULL);
wlu_addwrd(ws,"PRT",DBNA_prot,NULL);
return(ws);
}

WRD_LUSTRCT *db_getoddwrdsdstrct(void)
  /* create and fill the datastructure for odd words, return
 the address of the structure created. */
{
WRD_LUSTRCT *ws;

ws = wlu_getlustrct(1,DB_unknown);
wlu_addwrd(ws,"standard",DB_std,NULL);
wlu_addwrd(ws,"unreviewed",DB_unrev,NULL);
wlu_addwrd(ws,"preliminary",DB_prelim,NULL);
wlu_addwrd(ws,"unannotated",DB_unannot,NULL);
wlu_addwrd(ws,"backbone",DB_backbone,NULL);
wlu_addwrd(ws,"circular",DB_circlr,NULL);
wlu_addwrd(ws,"linear",DB_linear,NULL);
return(ws);
}

WRD_LUSTRCT *db_getstrndstrct(void)
  /* create and fill the datastructure for strandedness words, return
 the address of the structure created. */
{
WRD_LUSTRCT *ws;

ws = wlu_getlustrct(1,DB_unkstrnd);
wlu_addwrd(ws,"ss-",DB_snglstrnd,NULL);
wlu_addwrd(ws,"ds-",DB_dblstrnd,NULL);
wlu_addwrd(ws,"ms-",DB_mxdstrnd,NULL);
return(ws);
}

DBLU_LINECLASS db_classfylin(char *lin,
                             DBLU_DBFMT dfmt,
                             WRD_LUSTRCT *ws)
/* return the class of this line */
{
char tok[DBLU_GENBHDR+1];
int cp;
int lim;
char *dp;
char nc;

if (lin == NULL)
  return(DBLC_terminator);
else
  {
  switch (dfmt)
    {
    case DBFMT_genbank:
      lim = DBLU_GENBHDR;
      break;
    case DBFMT_embl:
    case DBFMT_embl_old:
    case DBFMT_swiss:
    case DBFMT_sqmonk:
      lim = DBLU_MBLHDR;
      break;
    case DBFMT_unknown:
    default:
      return(DBLC_terminator);
      break;
    }
  cp = 0;
  dp = &tok[0];
  while (cp < lim)
    {
    if (((nc = *(lin+cp)) == '\n') || (nc == '\0'))
      cp = lim;
    else
      bas_appchr(&tok[0],&dp,nc,DBLU_GENBHDR+1);
    cp++;
    }
  while ((dp >= &tok[0]) && ((*dp == ' ') || (*dp == '\0')))
    *dp-- = '\0';
  return(wlu_chkwrd(ws,&tok[0]));
  }
}

DBLU_LINECLASS db_classfyentln(DB_ENTSTRCT *es)
  /* return the line classification of es line */
{
return(db_classfylin(es->slbuf,es->efmt,es->lclss));
}

DBLU_LINECLASS db_actlintyp(DB_ENTSTRCT *es,
                            DBLU_LINECLASS fndlt)
/* return the actual line type, depending on format and previous linetype */
{
switch (es->efmt)
  {
  case DBFMT_embl:
  case DBFMT_embl_old:
  case DBFMT_swiss:
  case DBFMT_sqmonk:
    return(fndlt);
    break;
  case DBFMT_genbank:
  default:
    if (fndlt == DBLC_unk)
      return(es->prvltyp);
    else
      return(fndlt);
    break;
  }
}

char *db_linclass2str(DBLU_LINECLASS lc)
  /* return a string for the linetype - (don't take account of database
format at this point */
{
switch (lc)
  {
  case DBLC_ident:
    return("ident");
    break;
  case DBLC_access:
    return("access");
    break;
  case DBLC_naident:
    return("naident");
    break;
  case DBLC_datestamp:
    return("datestamp");
    break;
  case DBLC_description:
    return("description");
    break;
  case DBLC_keyword:
    return("keyword");
    break;
  case DBLC_source:
    return("source");
    break;
  case DBLC_taxonomy:
    return("taxonomy");
    break;
  case DBLC_reference:
    return("reference");
    break;
  case DBLC_author:
    return("author");
    break;
  case DBLC_title:
    return("title");
    break;
  case DBLC_journal:
    return("journal");
    break;
  case DBLC_remark:
    return("remark");
    break;
  case DBLC_feature:
    return("feature");
    break;
  case DBLC_basecount:
    return("basecount");
    break;
  case DBLC_sequence:
    return("sequence");
    break;
  case DBLC_terminator:
    return("terminator");
    break;
  case DBLC_unk:
  default:
    return("unknown");
    break;
  }
}

char *db_oddwrd2str(DB_ODDWRDS ow)
  /* return a string for the odd word ow */
{
switch (ow)
  {
  case DB_std:
    return("standard");
    break;
  case DB_unrev:
    return("unreviewed");
    break;
  case DB_prelim:
    return("preliminary");
    break;
  case DB_unannot:
    return("unannotated");
    break;
  case DB_backbone:
    return("backbone");
    break;
  case DB_circlr:
    return("circular");
    break;
  case DB_linear:
    return("linear");
    break;
  default:
    return("");
    break;
  }
}  

char *db_natyp2str(DB_NATYPE nat)
  /* return a string for the type of this na */
{
switch (nat)
  {
  case DBNA_dna:
    return("DNA");
    break;
  case DBNA_rna:
    return("RNA");
    break;
  case DBNA_trna:
    return("tRNA");
    break;
  case DBNA_rrna:
    return("rRNA");
    break;
  case DBNA_mrna:
    return("mRNA");
    break;
  case DBNA_urna:
    return("uRNA");
    break;
  case DBNA_prot:
    return("Protein");
    break;
  case DBNA_unk:
    return("");
    break;
  default:
    return("?NA?");
    break;
  }
}

char *db_strnd2str(DB_NASTRAND nas)
  /* return a string for nas */
{
switch (nas)
  {
  case DB_snglstrnd:
    return("ss-");
    break;
  case DB_dblstrnd:
    return("ds-");
    break;
  case DB_mxdstrnd:
    return("mixedstrnd-");
    break;
  case DB_unkstrnd:
    return("");
    break;
  default:
    return("?strd-");
    break;
  }
}

DB_STR_ELT *dbp_appstrelt(DB_STR_ELT **spt,
                          char *str)
/* malloc and append a string element to list *spt.
return the address of the new element.
The new value is strduped, to avoid unexpected
changes/behavior */
{
DB_STR_ELT *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nxtselt;
    }
  endp = (DB_STR_ELT *) getmemory(sizeof(DB_STR_ELT),"str elt");
  endp->nxtselt = NULL;
  endp->strval = bas_strdup(str);;
  if (*spt == NULL)
    {
    *spt = endp;
    endp->prvselt = NULL;
    }
  else
    {
    prev->nxtselt = endp;
    endp->prvselt = prev;
    }
  }
return(endp);
}

void db_delstrelt(DB_STR_ELT *sp,
                  DB_STR_ELT **lstrt,
                  DB_STR_ELT **lend)
/* delete sp from list *lstrt..*lend */
{
DB_STR_ELT *pt;

if (sp != NULL)
  {
  if ((pt = sp->prvselt) == NULL)
    *lstrt = sp->nxtselt;
  else
    pt->nxtselt = sp->prvselt;
  if ((pt = sp->nxtselt) != NULL)
    pt->prvselt = sp->prvselt;
  else
    if (lend != NULL)
      *lend = sp->prvselt;
  memfree(sp->strval);
  memfree(sp);
  }
}

void db_killstrlst(DB_STR_ELT **lstrt,
                   DB_STR_ELT **lend)
/* iteratively remove the header from *lstrt */
{
while (*lstrt != NULL)
  db_delstrelt(*lstrt,lstrt,lend);
}

void db_initfeatstrct(DB_FEATSTRCT *fs,
                      char *id)
/* set initial values for this feature structure */
{
fs->fstseg = fs->lstseg = NULL;
fs->seg_append_ok = 0;
fs->featur = FTKW_unknown;
fs->savdno = fs->featno = 0;
fs->strctsens = DB_sens_unk;
fs->cod_start = 0;
fs->gencod = GC_undefined;
fs->pseudogn = 0;
fs->partcds = 0;
fs->infolist = fs->infoend = NULL;
fs->mRNA = NULL;
fs->nextst = fs->prevst = 0;
fs->idp = id;
fs->enamelist = fs->enamlstend = NULL;
fs->enamlstend =  dbp_appstrelt(&fs->enamelist,id);
fs->parentname = NULL;
fs->parentfeat = NULL;
}

DB_SEGELT *db_appnseg(DB_SEGELT **spt,
                      int strt,
                      int stop,
                      DB_SENS sens,
                      char *xidnm)
/* malloc and append a new segment to list *spt.  return the address of the new
element */
{
DB_SEGELT *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nextseg;
    }
  endp = (DB_SEGELT *) getmemory(sizeof(DB_SEGELT),"Segment element");
  endp->nextseg = NULL;
  endp->segid = NULL;
  endp->sgstart = strt;
  endp->sgstop = stop;
  endp->lthan = 0;
  endp->gthan = 0;
  endp->sgsens = sens;
  endp->segid = xidnm;  /* assumed to have been malloc()ed already */
  if (*spt == NULL)
    {
    *spt = endp;
    endp->prevseg = NULL;
    }
  else
    {
    prev->nextseg = endp;
    endp->prevseg = prev;
    }
  }
return(endp);
}

DB_SEGELT *db_prepndseg(DB_SEGELT **spt,
                        int strt,
                        int stop,
                        DB_SENS sens,
                        char *xidnm)
/* malloc and prepend a new segment to list *spt.  return the address of the new
element */
{
DB_SEGELT *prev, *endp;

endp = (DB_SEGELT *) getmemory(sizeof(DB_SEGELT),"Segment element");
endp->nextseg = NULL;
endp->segid = NULL;
endp->sgstart = strt;
endp->sgstop = stop;
endp->lthan = 0;
endp->gthan = 0;
endp->sgsens = sens;
endp->segid = xidnm;  /* assumed to have been malloc()ed already */
if (*spt == NULL)
  endp->nextseg = endp->prevseg = NULL;
else
  {
  prev = *spt;
  prev->prevseg = endp;
  endp->nextseg = prev;
  endp->prevseg = NULL;
  }
*spt = endp;
return(endp);
}

void db_delsegelt(DB_SEGELT *sp,
                  DB_SEGELT **lstrt,
                  DB_SEGELT **lend)
/* delete sp from list *lstrt..*lend */
{
DB_SEGELT *pt;

if (sp != NULL)
  {
  if ((pt = sp->prevseg) == NULL)
    *lstrt = sp->nextseg;
  else
    pt->nextseg = sp->prevseg;
  if ((pt = sp->nextseg) != NULL)
    pt->prevseg = sp->prevseg;
  else
    *lend = sp->prevseg;
  if (sp->segid != NULL)
    memfree(sp->segid);
  memfree(sp);
  }
}

DB_STR_ELT *db_matstrelt(DB_STR_ELT *strlst,
                         char *qstr,
                         int (* scmpfun)(const char *x1,
                                         const char *x2))
/* return an element in strlst matching qstr, using
scmpfun for comparison.  return NULL if none */
{
DB_STR_ELT *sp;

if (qstr != NULL)
  {
  sp = strlst;
  while (sp != NULL)
    if ((*scmpfun)(sp->strval,qstr) == 0)
      return(sp);
    else
      sp = sp->nxtselt;
  return(NULL);
  }
else
  return(NULL);
}  

int db_cntstrlst(DB_STR_ELT *strlst)
  /* recursively count strlst elements */
{
if (strlst == NULL)
  return(0);
else
  return(db_cntstrlst(strlst->nxtselt) + 1);
}

DB_TBLINFO *db_appnielt(DB_TBLINFO **spt,
                        DBLU_FTQUAL fql,
                        char *info)
/* malloc and append a information element to list *spt.
  return the address of the new element. info may have
  needed to be malloced.*/
{
DB_TBLINFO *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
/* printf("App: '%s' = %s\n",db_ftqu2str(fql),info); */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nxtielt;
    }
  endp = (DB_TBLINFO *) getmemory(sizeof(DB_TBLINFO),"Qual info element");
  endp->nxtielt = NULL;
  endp->qval = info;
  endp->fqual = fql;
  if (*spt == NULL)
    {
    *spt = endp;
    endp->prvielt = NULL;
    }
  else
    {
    prev->nxtielt = endp;
    endp->prvielt = prev;
    }
  }
return(endp);
}

void db_delielt(DB_TBLINFO *sp,
                DB_TBLINFO **lstrt,
                DB_TBLINFO **lend)
/* delete sp from list *lstrt..*lend */
{
DB_TBLINFO *pt;

if (sp != NULL)
  {
  if ((pt = sp->prvielt) == NULL)
    *lstrt = sp->nxtielt;
  else
    pt->nxtielt = sp->prvielt;
  if ((pt = sp->nxtielt) != NULL)
    pt->prvielt = sp->prvielt;
  else
    *lend = sp->prvielt;
  memfree(sp->qval);
  memfree(sp);
  }
}

void db_killtblinfolst(DB_TBLINFO **lstrt,
                       DB_TBLINFO **lend)
/* iteratively remove the header from *lstrt */
{
while (*lstrt != NULL)
  db_delielt(*lstrt,lstrt,lend);
}

DB_TBLINFO *db_tblent4udata(DB_TBLINFO *tlst,
                            DBLU_FTQUAL uqul,
                            char *ustr)
/* scan tlst for the first element which matches uqul and (if non-NULL) 
contains ustr.  Return a pointer to that element, NULL for no match */
{
DB_TBLINFO *tp;

tp = tlst;
while (tp != NULL)
  {
  if (tp->fqual == uqul)
    if (ustr == NULL)
      return(tp);
    else
      if (strstr(tp->qval,ustr) != NULL)  /* found a string match */
        return(tp);
  tp = tp->nxtielt;
  }
return(NULL);
}

DB_LCINFO *db_appnlcelt(DB_LCINFO **spt,
                        DBLU_LINECLASS lcl,
                        char *info)
/* malloc and append a information element to list *spt.
  return the address of the new element */
{
DB_LCINFO *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nxtlcinf;
    }
  endp = (DB_LCINFO *) getmemory(sizeof(DB_LCINFO),"Qual info element");
  endp->nxtlcinf = NULL;
  endp->cinfo = info;
  endp->lclss = lcl;
  if (*spt == NULL)
    {
    *spt = endp;
    endp->prvlcinf = NULL;
    }
  else
    {
    prev->nxtlcinf = endp;
    endp->prvlcinf = prev;
    }
  }
return(endp);
}

void db_dellcelt(DB_LCINFO *sp,
                 DB_LCINFO **lstrt)
/* delete sp from list *lstrt */
{
DB_LCINFO *pt;

if (sp != NULL)
  {
  if ((pt = sp->prvlcinf) == NULL)
    *lstrt = sp->nxtlcinf;
  else
    pt->nxtlcinf = sp->prvlcinf;
  if ((pt = sp->nxtlcinf) != NULL)
    pt->prvlcinf = sp->prvlcinf;
  if (sp->cinfo != NULL)
    memfree(sp->cinfo);
  memfree(sp);
  }
}

void db_killlclsslst(DB_LCINFO **lstrt)
/* iteratively remove the header from *lstrt */
{
while (*lstrt != NULL)
  db_dellcelt(*lstrt,lstrt);
}

DB_LCINFO *db_lcent4udata(DB_LCINFO *tlst,
                          DBLU_LINECLASS uqul,
                          char *ustr)
/* scan tlst for the first element which matches uqul and (if non-NULL) 
contains ustr.  Return a pointer to that element, NULL for no match */
{
DB_LCINFO *tp;

tp = tlst;
while (tp != NULL)
  {
  if (tp->lclss == uqul)
    if (ustr == NULL)
      return(tp);
    else
      if (strstr(tp->cinfo,ustr) != NULL)  /* found a string match */
        return(tp);
  tp = tp->nxtlcinf;
  }
return(NULL);
}

DB_FEATSTRCT *db_appnfeat(DB_FEATSTRCT **spt,
                          DBLU_FTKW feat,
                          char *id)
/* malloc and append a new featstruct to list *spt.  return the address of 
the new element */
{
DB_FEATSTRCT *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nextst;
    }
  endp = (DB_FEATSTRCT *) getmemory(sizeof(DB_FEATSTRCT),"Feature element");
  db_initfeatstrct(endp,id);
  endp->featur = feat;
  endp->nextst = NULL;
  if (*spt == NULL)
    {
    *spt = endp;
    endp->prevst = NULL;
    }
  else
    {
    prev->nextst = endp;
    endp->prevst = prev;
    }
  }
return(endp);
}

DB_FEATSTRCT *dbp_prependfeatstrct(DB_FEATSTRCT **flststrtp,
                                   DB_FEATSTRCT *spt,
                                   DBLU_FTKW feat,
                                   char *id)
/* malloc and prepend a new featstruct before element *spt.  return the address of 
the new element */
{
DB_FEATSTRCT *prefp;
DB_FEATSTRCT *newelt;

if (spt == NULL)
  return(NULL);
else
  {
  newelt = NULL;
  newelt = db_appnfeat(&newelt,feat,id);
  if (spt == NULL)       /* haven't got anything in list yet, doubtful, though */
    spt = newelt;
  else
    {
    newelt->nextst = spt;
    prefp = spt->prevst;
    spt->prevst = newelt;
    if (prefp == NULL) /* was at head of list */
      *flststrtp = newelt;
    else
      prefp->nextst = newelt;
    newelt->prevst = prefp;
    }
  return(newelt);
  }
}

void db_killsegs4featstrct(DB_FEATSTRCT *fs)
  /* iteratively delete the segment list for fs - memory is freed */
{
while (fs->fstseg != NULL)
  db_delsegelt(fs->fstseg,&fs->fstseg,&fs->lstseg);
}

void db_freeifok(char *pt)
  /* free memory pt if non-NULL */
{
if (pt != NULL)
  memfree(pt);
}

DB_ACCELT *db_appnacc(DB_ACCELT **spt,
                      char *acnam)
/* malloc and append a new segment to list *spt.  return the address of the new
element */
{
DB_ACCELT *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nxtacc;
    }
  endp = (DB_ACCELT *) getmemory(sizeof(DB_ACCELT),"Segment element");
  endp->nxtacc = NULL;
  endp->accno = acnam;
  if (*spt == NULL)
    *spt = endp;
  else
    prev->nxtacc = endp;
  }
return(endp);
}

void db_delfstaccelt(DB_ACCELT **lstrt)
/* delete *lstrt */
{
DB_ACCELT *pt;

if (*lstrt != NULL)
  {
  pt = *lstrt;
  *lstrt = pt->nxtacc;
  db_freeifok(pt->accno);
  memfree(pt);
  }
}

void db_killacclst(DB_ACCELT **alst)
  /* iteratively delete the accno list alst - memory is freed */
{
while (*alst != NULL)
  db_delfstaccelt(alst);
}

DB_ACCELT *db_fndaccno(DB_ACCELT *alst,
                       char *acc)
/* return a pointer to any element of alst which matches acc.  Don't do
any case conversion */
{
DB_ACCELT *lp;

lp = alst;
while (lp != NULL)
  if (strcmp(acc,lp->accno) == 0)
    return(lp);
  else
    lp = lp->nxtacc;
return(NULL);
}

DB_WANTQU *db_appnquelt(DB_WANTQU **spt,
                        DBLU_FTQUAL qval)
/* malloc and append a new value to list *spt.  return the address of the new
element */
{
DB_WANTQU *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nxtwntqu;
    }
  endp = (DB_WANTQU *) getmemory(sizeof(DB_WANTQU),"want element");
  endp->nxtwntqu = NULL;
  endp->fqu = qval;
  if (*spt == NULL)
    {
    *spt = endp;
    endp->prvwntqu = NULL;
    }
  else
    {
    prev->nxtwntqu = endp;
    endp->prvwntqu = prev;
    }
  }
return(endp);
}

DB_WANTQU *db_appnquelt2felt(DB_FTYPELT *fep,
                             DBLU_FTQUAL qval)
/* malloc and append a new qualifier value for feature element fep.
return the address of the new element */
{
if (fep != NULL)
  return(db_appnquelt(&fep->wquals,qval));
else
  return(NULL);
}

DB_WANTQU *db_appnquelt4feat(DB_FTYPELT *flst,
                             DBLU_FTKW fkw,
                             DBLU_FTQUAL qval)
/* scan flst for a value matching fkw, malloc and append a new qualifier value
for it. return the address of the new element */
{
DB_FTYPELT *fep;

if ((fep = db_ptr4felt(flst,fkw)) != NULL)
  return(db_appnquelt2felt(fep,qval));
else
  return(NULL);
}

int db_appnquelt4all(DB_FTYPELT *flst,
                     DBLU_FTQUAL qval)
/* malloc and append a new qualifier value for all flst elts.
return number of elements entered */
{
DB_FTYPELT *fep;
int cnt;

cnt = 0;
fep = flst;
while (fep != NULL)
  {
  if (db_appnquelt2felt(fep,qval) != NULL)
    cnt++;
  fep = fep->nxtfelt;
  }
return(cnt);
}

DB_WANTQU *db_ptr4qualval(DB_WANTQU *qlst,
                          DBLU_FTQUAL qv)
/* if an element matching qv exists on qlst, return its address */
{
DB_WANTQU *qp;

qp = qlst;
while (qp != NULL)
  if (qp->fqu == qv)
    return(qp);
  else
    qp = qp->nxtwntqu;
return(NULL);
}

void db_delwantquelt(DB_WANTQU *sp,
                     DB_WANTQU **lstrt)
/* delete sp from list *lstrt */
{
DB_WANTQU *pt;

if (sp != NULL)
  {
  if ((pt = sp->prvwntqu) == NULL)
    *lstrt = sp->nxtwntqu;
  else
    pt->nxtwntqu = sp->prvwntqu;
  if ((pt = sp->nxtwntqu) != NULL)
    pt->prvwntqu = sp->prvwntqu;
  memfree(sp);
  }
}

void db_beheadquallst(DB_WANTQU **qlst)
  /* remove the element *qlst */
{
db_delwantquelt(*qlst,qlst);
}

void db_killwntqulst(DB_WANTQU **lstrt)
/* iteratively remove the header from *lstrt */
{
while (*lstrt != NULL)
  db_beheadquallst(lstrt);
}

DB_FTYPELT *db_appnfelt(DB_FTYPELT **spt,
                        DBLU_FTKW eval)
/* malloc and append a new value to list *spt.  return the address of the new
element */
{
DB_FTYPELT *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nxtfelt;
    }
  endp = (DB_FTYPELT *) getmemory(sizeof(DB_FTYPELT),"Feature element");
  endp->nxtfelt = NULL;
  endp->fval = eval;
  endp->wquals = NULL;
  if (*spt == NULL)
    {
    *spt = endp;
    endp->prvfelt = NULL;
    }
  else
    {
    prev->nxtfelt = endp;
    endp->prvfelt = prev;
    }
  }
return(endp);
}

void db_addcntfelts(DB_FTYPELT **clst,
                    int lcls,
                    int cnt)
/* append lcls records cntX onto *clst */
{
while (cnt-- > 0)
  (void) db_appnfelt(clst,lcls);
}

int db_delfelt(DB_FTYPELT **fstrt,
               DB_FTYPELT *fp)
/* fp is an element of *fstrt, delete it, correct pointers */
{
DB_FTYPELT *prv;
DB_FTYPELT *nxt;

if ((*fstrt!= NULL) && (fp != NULL))
  {
  prv = *fstrt;
  while ((prv != NULL) && (prv->nxtfelt != fp))
    prv = prv->nxtfelt;
  if (prv != NULL)
    {
    prv->nxtfelt = fp->nxtfelt;
    prv->prvfelt = *fstrt;
    }
  if (*fstrt == fp)
    *fstrt = fp->nxtfelt;
  db_killwntqulst(&fp->wquals);
  memfree(fp);
  return(1);
  }      
else
  return(0);
}

void db_killfeltlst(DB_FTYPELT **alst)
  /* iteratively delete the feature list alst - memory is freed */
{
while (*alst != NULL)
  db_delfelt(alst,*alst);
}

DB_FTYPELT *db_ptr4felt(DB_FTYPELT *flst,
                        DBLU_FTKW eval)
/* return the (first) address of any element of flst which matches eval */
{
DB_FTYPELT *fp;

fp = flst;
while (fp != NULL)
  if (fp->fval == eval)
    return(fp);
  else
    fp = fp->nxtfelt;
return(NULL);
}

int db_delmatfelt(DB_FTYPELT **fstrt,
                  DBLU_FTKW eval)
/* delete all elements of *fstrt which match eval. return no deleted */
{
DB_FTYPELT *fp;
int dcnt;

dcnt = 0;
while ((fp = db_ptr4felt(*fstrt,eval)) != NULL)
  dcnt += db_delfelt(fstrt,fp);
return(dcnt);
}

void db_endfeatsplic(DB_ENTSTRCT *es)
  /* revert es splice components to an un-initiallised state */
{
es->splcft = NULL;
es->curseg = NULL;
es->sqpt = 0;
es->cursens = DB_sens_norm;
}

void db_initentstrct(DB_ENTSTRCT *es,
                     DBLU_DBFMT fmt,
                     WRD_LUSTRCT *ftkwlu)
/* initialise contents of an entry structure */
{
es->ename = NULL;
es->nseen = es->nfeats = 0;
es->acclst = NULL;
es->sqlen = 0;
es->seq = NULL;
es->featlst = es->lastfeat = NULL;
es->natype = DBNA_unk;
es->descr = es->date = es->nid = NULL;
es->circ = 0;
es->nawlu = db_getnatypestrct();
es->oddwlu = db_getoddwrdsdstrct();
es->strndwlu = db_getstrndstrct();
es->flocwlu = db_getftlocstrct();
es->fqulwlu = db_getftqualstrct();
es->fkwlu = ftkwlu;
es->intrst = NULL;
es->lclss = NULL;
/* es->fquals = NULL; */
es->efmt = fmt;
es->sfl = NULL;
es->cp = es->slbuf = NULL;
es->slblen = 0;
es->slcnt = 0;
es->curltyp = es->prvltyp = DBLC_unk;
es->onqual = 0;
es->curgcod = GC_undefined;
es->trnsmatrx = NULL;
es->cvecs = NULL;
db_endfeatsplic(es);
es->p5limt = es->p3limt = 0;
es->lcinfo = NULL;
es->seqversno = 0;
es->freefeatids = 0;
}

char *db_strdupskip(char *lin,
                    int smax,
                    char *skipset)
/* produce a strdup-ed copy of lin up to smax chars, stopping at first non 
skipset char after object.  If leading skipset chars, then skip them also.
  return new copy of string */
{
char cpy[DBLU_MAXSRCLEN+1];

(void) bas_strcpyskip(lin,smax,&cpy[0],DBLU_MAXSRCLEN,skipset);
return(bas_strdup(&cpy[0]));
}

DB_NATYPE db_scanidlin(char *slin,
                       DBLU_DBFMT fmt,
                       WRD_LUSTRCT *nawlu,
                       WRD_LUSTRCT *owlu,
                       WRD_LUSTRCT *swlu,
                       char **enm,
                       int *sqln,
                       char **dat,
                       int *crc,
                       DB_NASTRAND *nas,
                       int *sqvers)
/* do the raw work of pulling slin apart */
{
char token[DBLU_MAXSRCLEN+1];
char *lp;
DB_NATYPE nat;
char *nm;
SQFL_SQTOPOL stopol;
char *lincpy;
char *lincpori;
char **clp;
int tcnt;
char *tokns[DB_TOKEN_MAX];
int svers;
int topotok;
char *chrname;
char *cnameori;
char *chrtokns[DB_TOKEN_MAX];
int ctcnt;
char *pcolon;
char *strtok_cntxt;
char *newgbtok;
int i;

switch (fmt)
  {
  case DBFMT_swiss:  /* not properly implemented, yet */
  case DBFMT_embl_old:
  case DBFMT_embl:
    lincpy = lincpori = bas_strdup(slin);
    tcnt = 0;
    for (clp = tokns; (*clp = strsep(&lincpy," ;\t")) != NULL;)
      {
      if (**clp != '\0')
        {
        tcnt++;
        if (++clp >= &tokns[DB_TOKEN_MAX])
          break;
        }
      }
    if (tcnt >= 2)
      {
      nm = bas_strdup(tokns[1]);
      if (strcmp(tokns[2],"SV") == 0)
        {
        svers = (int) strtol(tokns[3],NULL,10);
        topotok = 4;
        }
      else
        {
        svers = 0;
        topotok = 2;
        }
      if (sqvers != NULL)
        *sqvers = svers;
      if ((strcasecmp(tokns[topotok],"circular") == 0) &&
            (fmt != DBFMT_swiss))
        stopol = SQTP_circular;
      else
        stopol = SQTP_linear;
      nat = (DB_NATYPE) wlu_chkwrd(nawlu,tokns[topotok+1]);
      if (sqln != NULL)
        *sqln = (int) strtol(tokns[tcnt-2],NULL,10);
      if (crc != NULL)
        *crc = 0;
      if (nas != NULL)
        *nas = DB_unkstrnd;
      memfree(lincpori);
      }
    else
      nm = bas_strdup("");
    break;
  case DBFMT_sqmonk:
    lincpy = lincpori = bas_strdup(slin);
    tcnt = 0;
    for (clp = tokns; (*clp = strsep(&lincpy," ;\t")) != NULL;)
      {
      if (**clp != '\0')
        {
        tcnt++;
        if (++clp >= &tokns[DB_TOKEN_MAX])
          break;
        }
      }
    switch (tcnt)   /* distinguish between GRCh37 & GRCh38 style */
      {
      case 7:       /* GRCh37-style */
        if ((pcolon = rindex(tokns[1],':')) == NULL)
          nm = bas_strdup(tokns[1]);
        else
          nm = bas_strdup(pcolon+1);
        stopol = SQTP_linear;
        svers = 0;
        nat = (DB_NATYPE) wlu_chkwrd(nawlu,tokns[4]);
        break;
      case 10:
      default:
        chrname = cnameori = bas_strdup(tokns[1]);
        ctcnt = 0;
        for ( clp = chrtokns; (*clp = strsep(&chrname,":")) != NULL;)
          {
          ctcnt++;
          if (++clp >= &chrtokns[DB_TOKEN_MAX])
            break;
          }
        nm = bas_strdup(chrtokns[2]);
        if (strcmp(tokns[2],"SV") == 0)
          svers = (int) strtol(tokns[3],NULL,10);
        if (strcmp(tokns[4],"circular") == 0)
          stopol = SQTP_circular;
        else
          stopol = SQTP_linear;
        nat = (DB_NATYPE) wlu_chkwrd(nawlu,tokns[6]);
        memfree(cnameori);
        break;
      }      
    if (sqvers != NULL)
      *sqvers = svers;
    if (sqln != NULL)
      *sqln = (int) strtol(tokns[tcnt-2],NULL,10);
    if (crc != NULL)
      *crc = 0;
    if (nas != NULL)
      *nas = DB_unkstrnd;
    memfree(lincpori);
    break;
  case DBFMT_genbank:
  case DBFMT_unknown:
  default:
/* pull line apart by tokens */
    lincpy = lincpori = bas_strdup(slin);
    tcnt = 0;
    for (newgbtok = strtok_r(lincpy," \t",&strtok_cntxt);
         ((newgbtok != NULL) && (tcnt < DB_TOKEN_MAX));
         newgbtok = strtok_r(NULL," \t",&strtok_cntxt))
    tokns[tcnt++] = bas_strdup(newgbtok);
    nm = bas_strdup(tokns[1]);
    if ((nat = (DB_NATYPE) wlu_chkwrd(nawlu,tokns[4])) > DBNA_unk)
      {
      if (dat != NULL)
        *dat = bas_strdup(tokns[tcnt-1]);
      if (crc != NULL)
        *crc = wlu_chkwrd(owlu,tokns[5]) == (int) DB_circlr;
      if (nas != NULL)
        *nas = wlu_chkwrd(swlu,tokns[3]);
      if (sqln != NULL)
        *sqln = atoi(tokns[2]);
      }	      
/* clean up malloc()ed memory */
    for (i = 0; i < tcnt; i++)
      memfree(tokns[i]);
    memfree(lincpori);
    break;
  }
if (enm != NULL)
  *enm = nm;
return(nat);
}

void db_unpickidline(char *slin,
                     DB_ENTSTRCT *es,
                     DBLU_DBFMT fmt)
/* examine slin as a locus/ident line and extract useful info into es */
{
es->natype = db_scanidlin(slin,fmt,es->nawlu,es->oddwlu,es->strndwlu,
                            &es->ename,&es->sqlen,&es->date,
                            &es->circ,&es->strnd,&es->seqversno);
}

char *db_fgets(char *buf,
               int blen,
               FILE *fl,
               int *lcnt)
/* perform fgets of fl, remove \n from string */
{
char *lp;
char *slim;
char *rslt;

rslt = fgets(buf,(blen+1),fl);
if (rslt != NULL)
  {
  lp = buf;
  slim = buf + blen;
  while ((lp <= slim) && (*lp != '\0'))
    {
    if (*lp == '\n')
      {
      (*lcnt)++;
      *lp = '\0';
      lp = slim;
      }
    lp++;
    }
  if (*buf == '\0') /* have read an empty line in effect, try again */
    return(db_fgets(buf,blen,fl,lcnt));
  else
    return(rslt);
  }
else
  return(NULL);
}

char *db_getesln(DB_ENTSTRCT *es)
  /* check if a new line is needed, if so get it */
{
if (es->cp == NULL)
  es->cp = db_fgets(es->slbuf,es->slblen,es->sfl,&es->slcnt);
else
  es->cp = es->slbuf;
return(es->cp);
}

void db_unpickaccno(char *lin,
                    DB_ACCELT **alst)
/* scan lin for a list of accession nos, append strdup-ed copy to *alst */
{
char buf[DBLU_MAXSRCLEN+1];
char *pt;
int len;

pt = lin;
while ((pt != NULL) && (*pt != '\0') && ((len = strlen(pt)) > 0))
  if ((pt = bas_strcpyskip(pt,len,&buf[0],DBLU_MAXSRCLEN,", \t;")) != NULL)
    if (db_fndaccno(*alst,&buf[0]) == NULL)  /* avoid duplicates */
      (void) db_appnacc(alst,bas_strdup(&buf[0]));
}

void db_tellacclst(FILE *fl,
                   DB_ACCELT *alp,
                   char *sep)
/* iteratively write alp names to fl */
{
DB_ACCELT *lp;

lp = alp;
while (lp != NULL)
  {
  fputs(lp->accno,fl);
  if ((lp = lp->nxtacc) != NULL)
    fputs(sep,fl);
  }
}

int db_scnfeatwrds(char *intwrds,
                   WRD_LUSTRCT *kwlu,
                   DB_FTYPELT **flst,
                   int wldok)
/* scan intwrds for space or comma separated set of words, try to parse by
kwlu and create linked list of entries.  Allow and respond to wildcard strings
if wldok */
{
char buf[DBLU_MAXSRCLEN+1];
char *lp;
int len;
DBLU_FTKW feat;
int fcnt;
int ap;
WRD_LU_REC *rep;

fcnt = 0;
if ((lp = intwrds) != NULL)
  while ((*lp != '\0') && ((len = strlen(lp)) > 0))
    {
    lp = bas_strcpyskip(lp,len,&buf[0],DBLU_MAXSRCLEN,", \t");
    if ((feat = wlu_chkwrd(kwlu,&buf[0])) != FTKW_unknown)
      {
      (void) db_appnfelt(flst,(int) feat);
      fcnt++;
      }
    else       /* does it contain a * and is wldok?? */
      if ((wldok || (kwlu->casedep == WLU_CASEWILD)) &&
            (index(&buf[0],'*') != NULL) && wlu_initwrdscan(kwlu,&ap,&rep))
        {
        do
          if (((kwlu->casedep == WLU_CASEIND) &&
                   bas_wldstrcmp(rep->wrd,&buf[0],bas_cmatnocas)) ||
                 ((kwlu->casedep == WLU_CASEDEP) &&
                   (bas_wldstrcmp(rep->wrd,&buf[0],bas_cmatcasdep))))
            {
            (void) db_appnfelt(flst,(DBLU_FTKW) rep->retval);
            fcnt++;
            }
        while ((rep = wlu_nxtwrd(kwlu,&ap,rep)) != NULL);
        }
    }  
return(fcnt);
}

int db_scnqulwrds(char *intwrds,
                  WRD_LUSTRCT *kwlu,
                  DB_WANTQU **flst,
                  int wldok)
/* scan intwrds for space or comma separated set of words, try to parse by
kwlu and create linked list of entries.  Allow and respond to wildcard strings
if wldok */
{
char buf[DBLU_MAXSRCLEN+1];
char *lp;
int len;
DBLU_FTQUAL feat;
int fcnt;
int ap;
WRD_LU_REC *rep;

fcnt = 0;
if ((lp = intwrds) != NULL)
  while ((*lp != '\0') && ((len = strlen(lp)) > 0))
    {
    lp = bas_strcpyskip(lp,len,&buf[0],DBLU_MAXSRCLEN,", \t");
    if ((feat = wlu_chkwrd(kwlu,&buf[0])) != FTQU_unknown)
      {
      (void) db_appnquelt(flst,feat);
      fcnt++;
      }
    else       /* does it contain a * and is wldok?? */
      if ((wldok || (kwlu->casedep == WLU_CASEWILD)) &&
            (index(&buf[0],'*') != NULL) && wlu_initwrdscan(kwlu,&ap,&rep))
        {
        do
          if (((kwlu->casedep == WLU_CASEIND) &&
                   bas_wldstrcmp(rep->wrd,&buf[0],bas_cmatnocas)) ||
                 ((kwlu->casedep == WLU_CASEDEP) &&
                   (bas_wldstrcmp(rep->wrd,&buf[0],bas_cmatcasdep))))
            {
            (void) db_appnquelt(flst,(DBLU_FTQUAL) rep->retval);
            fcnt++;
            }
        while ((rep = wlu_nxtwrd(kwlu,&ap,rep)) != NULL);
        }
    }  
return(fcnt);
}

int db_classfixdfld(char *buf,
                    int lbuf,
                    char *spos,
                    WRD_LUSTRCT *fclass,
                    char *skipst)
/* scan the fixed field, starting at *spos return any type found there */
{
char tbuf[DBLU_MAXSRCLEN+1];
char *sp;
char *dp;

sp = spos;
dp = &tbuf[0];
while (((sp - buf) < lbuf) && (index(skipst,*sp) == NULL) && (*sp != '\0'))
  *dp++ = *sp++;
*dp = '\0';
return(wlu_chkwrd(fclass,&tbuf[0]));
}

DBLU_FTQUAL db_isqual(char *buf,
                      int blen,
                      char *lpos,
                      WRD_LUSTRCT *qclass)
/* see if this is a qualifier line (contains / at lpos), if so return
  the type */
{
if (((lpos - buf) < blen) && (*lpos++ == '/'))
  return(db_classfixdfld(buf,blen,lpos,qclass," =\t"));
else
  return(FTQU_unknown);
}

char *db_ftqual2str(DBLU_FTQUAL fqul)
  /* return a string value for fqul */
{
switch (fqul)
  {
  case FTQU_anticodon:
    return("anticodon");
    break;
  case FTQU_bound_moiety:
    return("bound_moiety");
    break;
  case FTQU_citation:
    return("citation");
    break;
  case FTQU_codon:
    return("codon");
    break;
  case FTQU_codon_start:
    return("codon_start");
    break;
  case FTQU_cons_splice:
    return("cons_splice");
    break;
  case FTQU_direction:
    return("direction");
    break;
  case FTQU_EC_number:
    return("EC_number");
    break;
  case FTQU_evidence:
    return("evidence");
    break;
  case FTQU_frequency:
    return("frequency");
    break;
  case FTQU_function_qu:
    return("function");
    break;
  case FTQU_gene:
    return("gene");
    break;
  case FTQU_label_qu:
    return("label");
    break;
  case FTQU_mod_base:
    return("mod_base");
    break;
  case FTQU_note:
    return("note");
    break;
  case FTQU_number:
    return("number");
    break;
  case FTQU_organism:
    return("organism");
    break;
  case FTQU_partial:
    return("partial");
    break;
  case FTQU_phenotype:
    return("phenotype");
    break;
  case FTQU_product:
    return("product");
    break;
  case FTQU_protid:
    return("protein_id");
    break;
  case FTQU_pseudo:
    return("pseudo");
    break;
  case FTQU_rpt_family:
    return("rpt_family");
    break;
  case FTQU_rpt_type:
    return("rpt_type");
    break;
  case FTQU_rpt_unit:
    return("rpt_unit");
    break;
  case FTQU_standard_name:
    return("standard_name");
    break;
  case FTQU_transl_except:
    return("transl_except");
    break;
  case FTQU_type_qu:
    return("type");
    break;
  case FTQU_usedin:
    return("usedin");
    break;
  case FTQU_transl_table:
    return("transl_table");
    break;
  case FTQU_db_xref:
    return("db_xref");
    break;
  case FTQU_biotype:
    return("biotype");
    break;
  case FTQU_phosphothre:
    return("Phosphothreonine");
    break;
  case FTQU_phosphoser:
    return("Phosphoserine");
    break;
  case FTQU_unknown:
  default:
    return("??FT_qual??");
    break;
  }
}

DB_CTYPE db_classchr(char c)
  /* return a character classification of c */
{
switch (c)
  {
  case '\0':
    return(DB_eos);
    break;
  case '.':
  case '^':
    return(DB_sep_scan);
    break;
  case '(':
  case ')':
  case '<':
  case '>':
  case ',':
  case ':':
    return(DB_sep_immed);
    break;
  default:
    if (isalpha(c))
      return(DB_alpha);
   else
      if (isdigit(c))
        return(DB_numeric);
      else
        return(DB_unkchr);
  }
}

char *db_floctp2str(DBLU_FLOCTYPE flt)
  /* return a string for flt */
{
switch (flt)
  {
  case FLOC_complement:
    return("complement");
    break;
  case FLOC_join:
    return("join");
    break;
  case FLOC_order:
    return("order");
    break;
  case FLOC_group:
    return("group");
    break;
  case FLOC_one_of:
    return("one_of");
    break;
  case FLOC_replace:
    return("replace");
    break;
  case FLOC_lpar:
    return("'('-Left par");
    break;
  case FLOC_rpar:
    return("')'-Right par");
    break;
  case FLOC_drange:
    return("..-defined range");
    break;
  case FLOC_irange:
    return(".-inrange");
    break;
  case FLOC_btween:
    return("^-between");
    break;
  case FLOC_lthan:
    return("less_than");
    break;
  case FLOC_gthan:
    return("gtr_than");
    break;
  case FLOC_comma:
    return("','-comma");
    break;
  case FLOC_nvalue:
    return("numeric");
    break;
  case FLOC_xentry:
    return("IDcrossref");
    break;
  case FLOC_ateol:
    return("EOL");
    break;
  case FLOC_unknown:
  default:
    return("Unknown");
    break;
  }
}

DB_TOKTYPE db_tokclass(char *tbuf)
  /* establish what tbuf is */
{
char *p;
int acnt;
int dcnt;
int ocnt;

p = tbuf;
acnt = dcnt = ocnt = 0;
while (*p != '\0')
  {
  if (isalpha(*p))
    acnt++;
  else
    if (index("0123456789<>",(int) *p) != NULL)  /* allow <> with positions */
      dcnt++;
    else
      ocnt++;
  p++;
  }
if (ocnt > 0)
  return(DB_tokmixd);
else
  if ((dcnt > 0) && (acnt == 0))     /* number */
    return(DB_toknum);
  else
    if ((dcnt == 0) && (acnt > 0))   /* pure alpha */
      return(DB_tokalph);
    else
      if ((dcnt == 0) && (acnt == 0) && (ocnt == 0))
        return(DB_tok_mt);
      else
        return(DB_tokmixd);
}

DB_DLMTTYPE db_typedlmter(char *buf)
  /* return the type of delimiter currently at start of buf */
{
if (buf == NULL)
  return(DB_dlmt_not);
else
  switch (*buf)
    {
    case '\0':
    case '\n':
    case '\r':         /* shouldn't see this, but... */
      return(DB_dlmt_eol);
      break;
    case '(':
      return(DB_dlmt_lpar);
      break;
    case ')':
      return(DB_dlmt_rpar);
      break;
    case ',':
      return(DB_dlmt_comma);
      break;
    case ':':
      return(DB_dlmt_colon);
      break;
    default:              /* not recognised so far is this a ".."?? */
      if (strncmp("..",buf,2) == 0)
        return(DB_dlmt_2dots);
      else
        return(DB_dlmt_not);
      break;
    }
}

DB_DLMTTYPE db_gettok(char *dbuf,
                      char **dp,
                      char *tbuf,
                      int tlen)
/* fish a token out of dbuf into tbuf, setting *dp past end of token &
 delimiter.  Return type of delimiter */
{
DB_DLMTTYPE dt;
char *tp;

tp = tbuf;
while ((dt = db_typedlmter(*dp)) != DB_dlmt_eol)
  switch (dt)
    {
    case DB_dlmt_lpar:
    case DB_dlmt_rpar:
    case DB_dlmt_comma:
    case DB_dlmt_colon:
      (void) bas_chr2ubuf(tbuf,&tp,'\0',tlen);
      (*dp)++;     /* skip past delimiter */
      return(dt);
      break;
    case DB_dlmt_2dots:
      (void) bas_chr2ubuf(tbuf,&tp,'\0',tlen);
      (*dp) += 2;     /* skip past delimiter */
      return(dt);
      break;
    case DB_dlmt_eol:    /* EOL */
      return(DB_dlmt_eol);
      break;
    case DB_dlmt_not:   /* Not a recognised delimiter, append to tbuf */
    default:
      (void) bas_chr2ubuf(tbuf,&tp,**dp,tlen);
      (*dp)++;
      break;
    }
/* fell off end: return eol */
return(DB_dlmt_eol);
}

DBLU_FTKW db_clssfeatln(DB_ENTSTRCT *es)
  /* return (if any) the feature type of es */
{
return(db_classfixdfld(es->slbuf,es->slblen,(es->slbuf+DB_FEATINSET-1),
         es->fkwlu," \t"));
}

void db_ftok2str(FILE *ofl,
                 DBLU_FLOCTYPE floc,
                 char *tbuf)
/* for checking purposes, list token and floc value for it */
{
#ifdef DB_DEBUG

if (floc != FLOC_nvalue)
  fprintf(ofl,"%s-%s\n",tbuf,db_floctp2str(floc));
else
  fprintf(ofl,"%s-%s=%d\n",tbuf,db_floctp2str(floc),atoi(tbuf));

#endif
}

char *db_sens2str(DB_SENS sns)
  /* return a string value for sns */
{
switch (sns)
  {
  case DB_sens_norm:
    return("5'");
    break;
  case DB_sens_comp:
    return("3'");
    break;
  case DB_sens_mixd:
    return("3&5'");
    break;
  case DB_sens_unk:
  default:
    return("?Sens?");
    break;
  }
}

void db_sayqulval(FILE *fl,
                  DB_TBLINFO *inf,
                  char *spacer)
/* tell fl about fqul and inf, if inf is non-NULL */
{
if (inf != NULL)
  fprintf(fl,"%s: %s%s",db_ftqual2str(inf->fqual),inf->qval,
            (spacer==NULL?"":spacer));
}

int db_savdno4feat(DB_FEATSTRCT *fp)
  /* return the saved number for this feature - should correspond to old-style
feature numbering for fish_term */
{
if (fp != NULL)
  return(fp->savdno);
else
  return(0);
}

int db_featno4feat(DB_FEATSTRCT *fp)
  /* return the feature number for this feature - relative to whole feature
table */
{
if (fp != NULL)
  return(fp->featno);
else
  return(0);
}

DB_FEATSTRCT *db_nthfeat(DB_ENTSTRCT *es,
                         int nth,
                         int (*numbfn)(DB_FEATSTRCT *fs))
/* return ptr to the nth feature.  (*numbfn)() relates which type of numbering
to use.  NULL if not found */
{
DB_FEATSTRCT *fp;

if (es == NULL)
  return(NULL);
else
  {
  fp = es->featlst;
  while (fp != NULL)
    if ((*numbfn)(fp) == nth)
      return(fp);
    else
      fp = fp->nextst;
  return(NULL);
  }
}

DB_SEGELT *db_5pseg4feat(DB_FEATSTRCT *fp)
  /* return the 5' segment for this feature, NULL if non-existent */
{
if (fp == NULL)
  return(NULL);
else
  if (fp->strctsens == DB_sens_comp)
    return(fp->lstseg);
  else
    return(fp->fstseg);
}

DB_SEGELT *db_3pseg4feat(DB_FEATSTRCT *fp)
  /* return the 3' segment for this feature, NULL if non-existent */
{
if (fp == NULL)
  return(NULL);
else
  if (fp->strctsens == DB_sens_comp)
    return(fp->fstseg);
  else
    return(fp->lstseg);
}

int db_5ppos4seg(DB_SEGELT *sp)
  /* return the 5' position for this sequence segment */
{
if (sp == NULL)
  return(0);
else
  if (sp->sgsens == DB_sens_comp)
    return(sp->sgstop);
  else
    return(sp->sgstart);
}

int db_3ppos4seg(DB_SEGELT *sp)
  /* return the 3' position for this sequence segment */
{
if (sp == NULL)
  return(0);
else
  if (sp->sgsens == DB_sens_comp)
    return(sp->sgstart);
  else
    return(sp->sgstop);
}

int db_5ppos4feat(DB_FEATSTRCT *fp)
  /* return the 5' position for this feature */
{
return(db_5ppos4seg(db_5pseg4feat(fp)));
}

int db_3ppos4feat(DB_FEATSTRCT *fp)
  /* return the 3' position for this feature */
{
return(db_3ppos4seg(db_3pseg4feat(fp)));
}

int db_5pextnt4feat(DB_FEATSTRCT *fp)
  /* return the 5' extent of feature fp (i.e. its minimum position) */
{
return(imin(db_5ppos4feat(fp),db_3ppos4feat(fp)));
}

int db_3pextnt4feat(DB_FEATSTRCT *fp)
  /* return the 3' extent of feature fp (i.e. its maximum position) */
{
return(imax(db_5ppos4feat(fp),db_3ppos4feat(fp)));
}

void db_featextremes(DB_ENTSTRCT *es,
                     DB_FEATSTRCT *fp,
                     int *e5p,
                     int *e3p)
/* return to e5p & e3p (if non-null) the extreme 
extent of fp, noting that mixed splice features
misfire the usual checks */
{
DB_SEGELT *segp;
int p5min;
int p3max;

p3max = 0;
p5min = es->sqlen;
segp = fp->fstseg;
while (segp != NULL)
  {
  p3max = imax(p3max,imax(segp->sgstart,segp->sgstop));
  p5min = imin(p5min,imin(segp->sgstart,segp->sgstop));
  segp = segp->nextseg;
  }
if (e5p != NULL)
  *e5p = p5min;
if (e3p != NULL)
  *e3p = p3max;
}

DB_SEGELT *db_nxtseg4feat(DB_SEGELT *cseg,
                          DB_FEATSTRCT *fp)
/* return the next segment in order for this feature */
{
if ((fp == NULL) || (cseg == NULL))
  return(NULL);
else
  if (fp->strctsens == DB_sens_comp)
    return(cseg->prevseg);
  else
    return(cseg->nextseg);
}

DB_SEGELT *db_prvseg4feat(DB_SEGELT *cseg,
                          DB_FEATSTRCT *fp)
/* return the previous segment in order for this feature */
{
if ((fp == NULL) || (cseg == NULL))
  return(NULL);
else
  if (fp->strctsens == DB_sens_comp)
    return(cseg->nextseg);
  else
    return(cseg->prevseg);
}

DB_RELPOSTYPE db_segorder(DB_SEGELT *s1,
                          DB_SEGELT *s2,
                          DB_FEATSTRCT *fp)
/* return relationship of s1 & s2.  NOTE that if s1 is 3' to s2, then
DB_rp_3p is returned */
{
DB_SEGELT *sp;

if (s1 == s2)
  return(DB_rp_same);
else
  {
  sp = s1;
  while ((sp = db_nxtseg4feat(sp,fp)) != NULL)
    if (sp == s2)         /* yep, s2 is 5p wrt s1, so... */
      return(DB_rp_3p);
  sp = s1;                /* look in opposite direction */
  while ((sp = db_prvseg4feat(sp,fp)) != NULL)
    if (sp == s2)         /* yep, s2 is 3' wrt s1, so... */
      return(DB_rp_5p);
  return(DB_rp_not);      /* weren't on same feature */
  }
}

int db_is5pto(int p1,
              int p2,
              DB_SENS sens)
/* return 1 if p1 is 5' to p2, depending on sense.
  NOTE: returns 0 for p1=p2 */
{
if (sens == DB_sens_comp)
  return(p1 > p2);
else
  return(p1 < p2);
}

int db_segis5pto(DB_SEGELT *sp1,
                 DB_SEGELT *sp2)
/* return 1 if sp1 is 5' wrt sp2 */
{
return(db_is5pto(db_5ppos4seg(sp1),db_5ppos4seg(sp2),sp1->sgsens));
}

DBP_RELPOSTYPE dbp_relatepos2feat(DB_FEATSTRCT *fp,
                                  int qpos)
/* indicate the relative position of qpos wrt feat fp */
{
int ft5p;

ft5p = db_5ppos4feat(fp);
if (db_rnginclds(ft5p,qpos,db_3ppos4feat(fp)))
  return(DBP_relpos_internal);
else
  if (db_is5pto(qpos,ft5p,fp->strctsens))
    return(DBP_relpos_5pto);
  else
    return(DBP_relpos_3pto);
}

DBP_RELPOSTYPE dbp_relatefeat2posns(DB_FEATSTRCT *fp,
                                    int p5p,
                                    int p3p)
/* p5p and p3p are 5' & 3' positions in absolute sense.
return the relative position of that feature wrt
p5p..p3p */
{
DBP_RELPOSTYPE p5relp;
DBP_RELPOSTYPE p3relp;

p5relp = dbp_relatepos2feat(fp,p5p);
p3relp = dbp_relatepos2feat(fp,p3p);
switch (p5relp)
  {
  case DBP_relpos_internal:
    switch (p3relp)
      {
      case DBP_relpos_internal:
        return(DBP_relpos_internal);
        break;
      case DBP_relpos_3pto:
        return(DBP_relpos_3polap);
      case DBP_relpos_5pto:
        if (fp->strctsens == DB_sens_comp)
          return(DBP_relpos_5polap);
        else
          return(DBP_relpos_undef);
        break;
      default:
        return(DBP_relpos_undef);
        break;
      }
    break;
  case DBP_relpos_5pto:
    switch (p3relp)
      {
      case DBP_relpos_5pto:
        return(DBP_relpos_5pto);
        break;
      case DBP_relpos_internal:
        return(DBP_relpos_5polap);
        break;
      case DBP_relpos_3pto:
        return(DBP_relpos_included);
        break;
      default:
        return(DBP_relpos_undef);
        break;
      }
    break;
  case DBP_relpos_3pto:
    switch (p3relp)
      {
      case DBP_relpos_3pto:
        return(DBP_relpos_3pto);
        break;
      case DBP_relpos_5pto:
        return(DBP_relpos_included);
        break;
      case DBP_relpos_internal:
        return(DBP_relpos_3polap);
        break;
      default:
        return(DBP_relpos_undef);
        break;
      }
    break;
  default:
    return(DBP_relpos_undef);
    break;
  }
}

char *dbp_relpostype2str(DBP_RELPOSTYPE rpt)
  /* text for rpt */
{
switch (rpt)
  {
  case DBP_relpos_5pto:
    return("5'_to");
    break;
  case DBP_relpos_5polap:
    return("5'_overlap");
    break;
  case DBP_relpos_internal:
    return("internal");
    break;
  case DBP_relpos_included:
    return("included");
    break;
  case DBP_relpos_3polap:
    return("3'_overlap");
    break;
  case DBP_relpos_3pto:
    return("3'_to");
    break;
  case DBP_relpos_undef:
  default:
    return("Undefined");
    break;
  }
}

int db_crossplccnt(DB_FEATSTRCT *fs)
  /* return the number of cross entry splices for this feature */
{
DB_SEGELT *sp;
int xcnt;

xcnt = 0;
sp = fs->fstseg;
while (sp != NULL)
  {
  if (sp->segid != NULL)
    xcnt++;
  sp = sp->nextseg;
  }
return(xcnt);
}

void db_saysegposns(FILE *fl,
                    DB_SEGELT *sp)
/* tell fl about all of this segment */
{
int p5;
int p3;

p5 = db_5ppos4seg(sp);
p3 = db_3ppos4seg(sp);
if (sp->sgsens == DB_sens_comp)
  {
  if (sp->gthan)
    fputc('>',fl);
  if (sp->segid != NULL)
    fprintf(fl,"%s:",sp->segid);
  if (p5 == p3)
    fprintf(fl,"%d",p5);
  else
    if (sp->lthan)
      fprintf(fl,"%d..<%d",p5,p3);
    else
      fprintf(fl,"%d..%d",p5,p3);
  }
else
  {
  if (sp->lthan)
    fputc('<',fl);
  if (sp->segid != NULL)
    fprintf(fl,"%s:",sp->segid);
  if (p5 == p3)
    fprintf(fl,"%d",p5);
  else
    if (sp->gthan)
      fprintf(fl,"%d..>%d",p5,p3);
    else
      fprintf(fl,"%d..%d",p5,p3);
  }  
}

int db_featlength(DB_FEATSTRCT *fp,
                  int *segcnt)
/* scan this feature for length, if segcnt is non-NULL, return the count */
{
int scnt;
int ltot;
DB_SEGELT *sp;

if (fp == NULL)
  return(0);
else
  {
  scnt = ltot = 0;
  sp = fp->fstseg;
  while (sp != NULL)
    {
    scnt++;
    ltot += sp->sgstop - sp->sgstart + 1;
    sp = sp->nextseg;
    }
  if (segcnt != NULL)
    *segcnt = scnt;
  return(ltot);
  }
}

int db_rawfeatlength(DB_FEATSTRCT *fp)
  /* return the raw length of fp, irrespective of orientation.  ignores
introns */
{
if (fp != NULL)
  if (fp->strctsens == DB_sens_comp)
    return(db_5ppos4feat(fp) - db_3ppos4feat(fp) + 1);
  else
    return(db_3ppos4feat(fp) - db_5ppos4feat(fp) + 1);
else
  return(0);
}

int db_rnginclds(int p1,
                 int p2,
                 int p3)
/* return 1 if p1..p3 includes p2 (p1 <= p2 <= p3) or (p3 <= p2 <= p1) */
{
if (p1 > p3)
  return(db_rnginclds(p3,p2,p1));
else
  return((p1 <= p2) && (p2 <= p3));
}

int db_segcnt4feat(DB_FEATSTRCT *fs)
  /* return the number of segments in fs */
{
int scnt;

if (fs == NULL)
  return(0);
else
  {
  (void) db_featlength(fs,&scnt);
  return(scnt);
  }
}

DB_SEGELT *db_seglstelt4pos(DB_SEGELT *seglst,
                            int pos)
/* return the seg element pointer which contains pos, NULL
if none */
{
DB_SEGELT *sp;

sp = seglst;
while (sp != NULL)
  if (db_rnginclds(sp->sgstart,pos,sp->sgstop))
    return(sp);
  else
    sp = sp->nextseg;
return(NULL);
}

DB_SEGELT *db_featseg4pos(DB_FEATSTRCT *fs,
                          int pos)
/* return the feature segment of fs which contains pos.  NULL if none */
{
DB_SEGELT *sp;

if (fs == NULL)
  return(NULL);
else
  return(db_seglstelt4pos(fs->fstseg,pos));
}

DB_SEGELT *db_nxtfeatseg4pos(DB_FEATSTRCT *fs,
                             int pos)
/* if pos is on a segment, return that, else traverse
the segment list, returning the next segment after
pos in the sense direction.  NULL if isn't in fs range */
{
DB_SEGELT *segp;

if (!db_rnginclds(db_5pextnt4feat(fs),pos,db_3pextnt4feat(fs)))
  return(NULL);
else
  if ((segp = db_featseg4pos(fs,pos)) != NULL)
    return(segp);
  else
    if (fs->strctsens == DB_sens_comp)
      {
      segp = fs->lstseg;
      while (segp != NULL)
        if (pos > segp->sgstop)
          return(segp);
        else
          segp = segp->prevseg;
      return(NULL);
      }
    else
      {
      segp = fs->fstseg;
      while (segp != NULL)
        if (pos < segp->sgstart)
          return(segp);
        else
        segp = segp->nextseg;
      return(NULL);
      }
}

DB_FEATSTRCT *db_nxtfeat4pos(DB_FEATSTRCT *flst,
                             DBLU_FTKW ftype,
                             int gpos)
/* return pointer to the next feature including present
that embraces gpos.  NULL for none.  Match on ftype unless
it is FTKW_unknown. */
{
DB_FEATSTRCT *fp;

fp = flst;
while (fp != NULL)
  {
  if ((fp->featur == ftype) || (ftype == FTKW_unknown))
    if (db_rnginclds(db_5ppos4feat(fp),gpos,db_3ppos4feat(fp)))
      return(fp);
  fp = fp->nextst;
  }
return(NULL);
}

DB_FEATSTRCT *db_feat4pos(DB_ENTSTRCT *ep,
                          DBLU_FTKW ftype,
                          int gpos)
/* return the first feature of *ep that includes
position gpos.  Match on ftype unless
it is FTKW_unknown. */
{
return(db_nxtfeat4pos(ep->featlst,ftype,gpos));
}

DB_FEATSTRCT *db_nxtfeat4typenid(DB_FEATSTRCT *ftlst,
                                 DBLU_FTKW fttype,
                                 char *queryid,
                                 int (*scmpfn)(const char *x1,
                                               const char *x2))
/* look thru ftlst for next feature with an id
matching queryid as indicated by (scmpfn*)() and of 
type fttype (if this is not unknown) */
{
DB_FEATSTRCT *fp;

fp = ftlst;
while (fp != NULL)
  if (((fttype == fp->featur) || (fttype == FTKW_unknown))
        && ((*scmpfn)(queryid,fp->idp) == 0) ||
        (db_matstrelt(fp->enamelist,queryid,scmpfn) != NULL))
    return(fp);
  else
    fp = fp->nextst;
return(NULL);
}

DB_FEATSTRCT *db_nxtfeat4id(DB_FEATSTRCT *ftlst,
                            char *queryid,
                            int (*scmpfn)(const char *x1,
                                          const char *x2))
/* look thru ftlst for next feature with an id
matching querid as indicated by (scmpfn*)() */
{
return(db_nxtfeat4typenid(ftlst,FTKW_unknown,queryid,scmpfn));
}

DB_FEATSTRCT *db_fndcdsandmrna4pos(DB_FEATSTRCT *fstrt,
                                   int gpos,
                                   DB_SPTSTAT *postype)
/* from feature list fstrt try to find a CDS with or without mrna which includes
gpos.  Return the cds feature pointer for this.  if
postype is non-NULL, return a position status value to it */
{
DB_FEATSTRCT *fp;
DB_FEATSTRCT *mrnap;

/* hierarchically, look for CDS first */
fp = fstrt;
while (fp != NULL)
  {
  if (fp->featur == FTKW_CDS)
    if (db_rnginclds(db_5ppos4feat(fp),gpos,db_3ppos4feat(fp)))
      {
      if (postype != NULL)
        *postype = DBPT_incds;
      return(fp);
      }
    else
      if ((fp->mRNA != NULL) && (db_rnginclds(db_5ppos4feat(fp->mRNA),gpos,db_3ppos4feat(fp->mRNA))))
        {
        if (postype != NULL)
          {
          if (db_rnginclds(db_5ppos4feat(fp->mRNA),gpos,db_5ppos4feat(fp)))
            *postype = DBPT_5ptocds;
          else
            if (db_rnginclds(db_3ppos4feat(fp),gpos,db_3ppos4feat(fp->mRNA)))
              *postype = DBPT_3ptocds;
            else
              *postype = DBPT_undef;
          }
        return(fp);
        }
  fp = fp->nextst;
  }
/* fell off end */
return(NULL);
}

DB_FEATSTRCT *dbp_featincludesfeat(DB_FEATSTRCT *fp,
                                   DB_FEATSTRCT *flst,
                                   DBLU_FTKW ftype)
/* return the structure address of the first feature in flst which
is of type ftype (if not FTKW_unknown) and is included in feature
fp */
{
DB_FEATSTRCT *fpt;
int p5ext;
int p3ext;

p5ext = db_5pextnt4feat(fp);
p3ext = db_3pextnt4feat(fp);
fpt = flst;
while (fpt != NULL)
  if (((fp != fpt) && ((ftype == FTKW_unknown) || (fpt->featur == ftype))) &&
         (db_rnginclds(p5ext,db_5pextnt4feat(fpt),p3ext) &&
         db_rnginclds(p5ext,db_3pextnt4feat(fpt),p3ext)))
    return(fpt);
  else
    fpt = fpt->nextst;
return(NULL);
}

char *db_sptstat2str(DB_SPTSTAT sstat)
  /* return a relevant string for sstat */
{
switch (sstat)
  {
  case DBPT_incds:
    return("CDS");
    break;
  case DBPT_5ptocds:
    return("5'UTR");
    break;
  case DBPT_3ptocds:
    return("3'UTR");
    break;
  default:
    return("Unknown");
    break;
  }
}

int db_segsoverlap(DB_SEGELT *seglst)
  /* check if any segments in seglst overlap */
{
if (seglst == NULL)
  return(0);
else
  return((db_seglstelt4pos(seglst->nextseg,seglst->sgstart) != NULL) ||
            (db_seglstelt4pos(seglst->nextseg,seglst->sgstop) != NULL) ||
            db_segsoverlap(seglst->nextseg));
}

int db_segcnt4featseg(DB_FEATSTRCT *fp,
                      DB_SEGELT *sep)
/* count through fp segments, returning ordinal
position of sep in elements.  0 if not present */
{
int scnt;
DB_SEGELT *sp;

scnt = 1;
if (fp->strctsens == DB_sens_comp)
  {
  sp = fp->lstseg;
  while (sp != NULL)
    if (sp == sep)
      return(scnt);
    else
      {
      scnt++;
      sp = sp->prevseg;
      }
/* fell off end: not found */
  return(0);
  }
else
  {
  sp = fp->fstseg;
  while (sp != NULL)
    if (sp == sep)
      return(scnt);
    else
      {
      scnt++;
      sp = sp->nextseg;
      }
/* fell off end: not found */
  return(0);
  }
}

DB_OLAPTYP db_chkfeatolap(DB_FEATSTRCT *f1,
                          DB_FEATSTRCT *f2)
/* return an assessment of the feature overlap for f2 included in f1 */
{
DB_SEGELT *p5p;
DB_SEGELT *p3p;
int in5p;
int in3p;
int out5p;
int out3p;

if ((f1 == NULL) || (f2 == NULL))
  return(OLAP_none);
else
  {
  p5p = db_featseg4pos(f1,(in5p = db_5ppos4feat(f2)));
  p3p = db_featseg4pos(f1,(in3p = db_3ppos4feat(f2)));
  if ((p5p == NULL) || (p3p == NULL))    /* didn't overlap */
    return(OLAP_none);
  else
    {
    out5p = db_5ppos4feat(f1);
    out3p = db_3ppos4feat(f1);
    if (db_rnginclds(out5p,in5p,out3p) && db_rnginclds(out5p,in3p,out3p))
      if ((out5p == in5p) && (out3p == in3p))
        return(OLAP_abuts);
      else
        if (out5p == in5p)
          return(OLAP_5pabut);
        else
          if (out3p == in3p)
            return(OLAP_3pabut);
          else
            return(OLAP_incl);
    else
      return(OLAP_none);
    }
  }
}

int db_segissame(DB_SEGELT *s1,
                 DB_SEGELT *s2)
/* return 1 if s1 & s2 are identical in sense and position */
{
return((s1 != NULL) && (s2 != NULL) && (s1->sgsens == s2->sgsens) &&
         (s1->sgstart == s2->sgstart) && (s1->sgstop == s2->sgstop));
}

int db_segincldseg(DB_SEGELT *sout,
                   DB_SEGELT *sin)
/* return true if sin is included in sout */
{
int olm5;
int olm3;

olm5 = db_5ppos4seg(sout);
olm3 = db_3ppos4seg(sout);
return(db_rnginclds(olm5,db_5ppos4seg(sin),olm3) &&
          db_rnginclds(olm5,db_3ppos4seg(sin),olm3));
}

int db_restsplicsame(DB_FEATSTRCT *fout,
                     DB_FEATSTRCT *fin,
                     DB_SEGELT *oseg,
                     DB_SEGELT *iseg,
                     int scnt)
/* recursively check if oseg.. is a valid identical splice for iseg */
{
if (db_segissame(oseg,iseg))  /* yep these match: look for rest */
  return(db_restsplicsame(fout,fin,db_nxtseg4feat(oseg,fout),
                            db_nxtseg4feat(iseg,fin),++scnt));
else
  if (iseg == NULL)   /* run off end */
    return(1);
  else
    if ((scnt == 1) && db_segincldseg(oseg,iseg)) /* 5' end */
      return(db_restsplicsame(fout,fin,db_nxtseg4feat(oseg,fout),
                                db_nxtseg4feat(iseg,fin),++scnt));
    else
      if (db_nxtseg4feat(iseg,fin) == NULL)  /* at 3' end */
        return((oseg != NULL) && (db_5ppos4seg(iseg) == db_5ppos4seg(oseg)));
      else
        return(0);
}  

int db_featsplicsame(DB_FEATSTRCT *fout,
                     DB_FEATSTRCT *fin)
/* check that fin and fout use the same splicing pattern over the range of
fin */
{
return(db_restsplicsame(fout,fin,db_featseg4pos(fout,db_5ppos4feat(fin)),
                          db_featseg4pos(fin,db_5ppos4feat(fin)),1));
}

int db_intrncnt4seg(DB_FEATSTRCT *fs,
                    DB_SEGELT *sp,
                    int to3p)
/* return count of introns from sp to 3' (to3p) or 5' (!to3p) direction */
{
int icnt;

icnt = 0;
while (sp != NULL)
  {
  if (to3p)
    sp = db_nxtseg4feat(sp,fs);
  else
    sp = db_prvseg4feat(sp,fs);
  if (sp != NULL)
    icnt++;
  }
return(icnt);
}

int db_intrncnt4pos(DB_FEATSTRCT *fs,
                    int posn,
                    int to3p)
/* return count of the introns from posn to end of fs in 3' (to3p) or 
5' (!to3p) direction */
{
return(db_intrncnt4seg(fs,db_featseg4pos(fs,posn),to3p));
}

int db_intrnsinrng(DB_FEATSTRCT *fs,
                  int pstrt,
                  int pstop)
/* return count of the fs introns from pstrt to pstop */
{
int icnt;
DB_SEGELT *spstrt;
int to5p;
DB_SEGELT *spstop;

if (fs != NULL)
  {
  to5p = db_is5pto(pstop,pstrt,fs->strctsens);
  icnt = 0;
  spstop = db_featseg4pos(fs,pstop);
  spstrt = db_featseg4pos(fs,pstrt);
  while ((spstrt != NULL) && (spstop != NULL))
    if (spstrt == spstop)
      return(icnt);
   else
      {
      icnt++;
      if (to5p)
        spstrt = db_prvseg4feat(spstrt,fs);
      else
        spstrt = db_nxtseg4feat(spstrt,fs);
      }
  return(0);
  }
else
  return(0);
}

DB_FEATSTRCT *db_gene4mrna(DB_ENTSTRCT *es,
                           DB_FEATSTRCT *mfs)
/* look through es feat list and return the first gene feature that 
encompasses mfs (must have same sense and start&end match) */
{
DB_FEATSTRCT *fp;

if (mfs != NULL)
  {
  fp = es->featlst;
  while (fp != NULL)
    if ((fp->featur == FTKW_FT_gene) && (db_chkfeatolap(fp,mfs) != OLAP_none)
         && (db_5ppos4feat(fp) == db_5ppos4feat(mfs))
         && (db_3ppos4feat(fp) == db_3ppos4feat(mfs)))
      return(fp);
    else
      fp = fp->nextst;
  }
return(NULL);
}

DB_FEATSTRCT *db_mrna4feat(DB_ENTSTRCT *es,
                           DB_FEATSTRCT *fs)
/* look through es feat list and return the first mRNA that encompasses fs
 (must have same sense and the splicing pattern must match) */
{
DB_FEATSTRCT *fp;

fp = es->featlst;
while (fp != NULL)
  if ((fp->featur == FTKW_mRNA) && (db_chkfeatolap(fp,fs) != OLAP_none) &&
        db_featsplicsame(fp,fs))
    return(fp);
  else
    fp = fp->nextst;
return(NULL);
}

int db_scnmrnas4feats(DB_ENTSTRCT *es,
                      DBLU_FTKW feat)
/* scan through feature list, for any of type feat, link up any mRNA's which
correspond to them and return count of No so matched */
{
int mcnt;
DB_FEATSTRCT *fs;

mcnt = 0;
fs = es->featlst;
while (fs != NULL)
  {
  if (fs->featur == feat)
    if ((fs->mRNA = db_mrna4feat(es,fs)) != NULL)
      mcnt++;
  fs = fs->nextst;
  }
return(mcnt);
}

int db_countfeats(DB_ENTSTRCT *es,
                  DBLU_FTKW feat)
/* scan through feature list, for any of type feat - return count.
if feat == FTKW_unknown, return count of all features */
{
int mcnt;
DB_FEATSTRCT *fs;

mcnt = 0;
fs = es->featlst;
while (fs != NULL)
  {
  if ((feat == FTKW_unknown) || (fs->featur == feat))
    mcnt++;
  fs = fs->nextst;
  }
return(mcnt);
}

void db_chktransmatrx(DB_ENTSTRCT *es,
                      DB_FEATSTRCT *fs)
/* check (and if necessary init) that a translation
matrix exists */
{
if (es->trnsmatrx == NULL)   /* must allocate a translation matrix */
  {
  es->trnsmatrx = (TRANS_MATRX *) getmemory(sizeof(TRANS_MATRX),
                                              "Genetic code table");
  if ((fs != NULL) &&(fs->gencod != GC_undefined))
    es->curgcod = fs->gencod;
  if (es->curgcod == GC_undefined)
    es->curgcod = GC_universal;
  init_trnsar4gb(*(es->trnsmatrx),es->curgcod);
  }
if ((fs != NULL) && (fs->gencod != GC_undefined) &&
       (es->curgcod != fs->gencod))
  init_trnsar4gb(*(es->trnsmatrx),(es->curgcod = fs->gencod));
if (es->cvecs == NULL)
  {
  es->cvecs = sqt_getcvctarr("Entry Struct codon vector table");
  (void) sqt_bldvects4matrx(*(es->trnsmatrx),es->cvecs);
  }
}

void db_initftsplic4pos(DB_ENTSTRCT *es,
                        DB_FEATSTRCT *fs,
                        int spos)
/* try to initiallise seq splicing for fs at position spos.  Set relevant
es parameters, even for invalid positions */
{
SQT_GENCOD gcod;

if ((fs != NULL) && (fs->mRNA != NULL))
  {
  db_initftsplic4pos(es,fs->mRNA,spos);
  es->p5limt = db_5ppos4feat(fs);
  es->p3limt = db_3ppos4feat(fs);
  }
else
  {
  es->splcft = fs;
  es->sqpt = spos;
  if (fs == NULL)
    {
    es->curseg = NULL;
    es->cursens = DB_sens_norm;
    es->p5limt = spos;
    es->p3limt = es->sqlen;
    }
  else
    {
    es->curseg = db_featseg4pos(fs,spos);
    es->cursens = es->curseg->sgsens; /* fs->strctsens; */
    es->p5limt = db_5ppos4feat(fs);
    es->p3limt = db_3ppos4feat(fs);
    }
  db_chktransmatrx(es,fs);
  }
}

void db_init5psplic4feat(DB_ENTSTRCT *es,
                         DB_FEATSTRCT *fs)
/* initialise extracting the spliced
 seq data for fs alone */
{
db_initftsplic4pos(es,fs,db_5ppos4feat(fs));
}

void db_init3psplic4feat(DB_ENTSTRCT *es,
                         DB_FEATSTRCT *fs)
/* initialise splicing for 3' end of fs */
{
db_initftsplic4pos(es,fs,db_3ppos4feat(fs));
}

int db_nxtsqptr(int sp,
                DB_SENS sns)
/* return incremented/decremented sp as required */
{
if (sns == DB_sens_comp)
  return(sp - 1);
else
  return(sp + 1);
}

void db_nxtesptr(DB_ENTSTRCT *es)
  /* increment the splice pointer in es to the next 5'->3' position */
{
if (es->splcft != NULL)
  switch (es->splcft->strctsens)
    {
    case DB_sens_comp:
      es->sqpt = db_nxtsqptr(es->sqpt,DB_sens_comp);
      break;
    case DB_sens_mixd:
      if (es->curseg != NULL)       /* then can respond intelligently */
        es->sqpt = db_nxtsqptr(es->sqpt,es->curseg->sgsens);
      else
        es->sqpt = db_nxtsqptr(es->sqpt,DB_sens_norm);
      break;
    case DB_sens_norm:
    default:
      es->sqpt = db_nxtsqptr(es->sqpt,DB_sens_norm);
      break;
    }
else
  es->sqpt = db_nxtsqptr(es->sqpt,DB_sens_norm);
}

int db_prvsqptr(int sp,
                DB_SENS sns)
/* return incremented/decremented sp as required */
{
if (sns == DB_sens_comp)
  return(sp + 1);
else
  return(sp - 1);
}

void db_prvesptr(DB_ENTSTRCT *es)
  /* increment the splice pointer in es to the next 3'->5' position */
{
if (es->splcft != NULL)
  switch (es->splcft->strctsens)
    {
    case DB_sens_comp:
      es->sqpt = db_prvsqptr(es->sqpt,DB_sens_comp);
      break;
    case DB_sens_mixd:
      if (es->curseg != NULL)       /* then can respond intelligently */
        es->sqpt = db_prvsqptr(es->sqpt,es->curseg->sgsens);
      else
        es->sqpt = db_prvsqptr(es->sqpt,DB_sens_norm);
      break;
    case DB_sens_norm:
    default:
      es->sqpt = db_prvsqptr(es->sqpt,DB_sens_norm);
      break;
    }
}

char db_getesresnoinc(DB_ENTSTRCT *es,
                      DB_SPLICLMT segsonly,
                      char (* bascmpfun)(char yb,
                                           BAS_MATMODE xmmod))
/* guts of the process to return the appropriate residue, irrespective any any complex links,
but don't increment/decrement seq pointer.
Depending on segsonly, return
null if there is not a valid segment or if pointer is out of feature limit.
(* bascmpfun)() is used to create base complement */
{
char nb;
DB_SEGELT *segp;

if (db_rnginclds(1,es->sqpt,es->sqlen))
  switch (segsonly)
    {
    case DBS_segonly:
      if (es->curseg == NULL)
        nb = '\0';
      else
        if (db_rnginclds(es->curseg->sgstart,es->sqpt,es->curseg->sgstop))
          {
          nb = *(es->seq + es->sqpt - 1);
          if (es->cursens == DB_sens_comp)
            nb = (* bascmpfun)(nb,BAS_exact);
          }
        else
          if ((es->curseg = db_nxtseg4feat(es->curseg,es->splcft)) != NULL)
            {
            es->sqpt = db_5ppos4seg(es->curseg);
            es->cursens = es->curseg->sgsens;
            return(db_getesresnoinc(es,segsonly,bascmpfun));
            }
          else
            return('\0');
      break;
    case DBS_featlimit:
      if (db_featseg4pos(es->splcft,es->sqpt) != NULL)
/*      if (db_rnginclds(es->p5limt,es->sqpt,es->p3limt)) */ /* don't work with circular seqs */
        {
        nb = *(es->seq + es->sqpt - 1);
        if (es->cursens == DB_sens_comp)
          nb = (* bascmpfun)(nb,BAS_exact);
        }
      else
        nb = '\0';
      break;
    case DBS_mrnalimit:
      if ((es->splcft != NULL) &&
           (db_featseg4pos(es->splcft,es->sqpt) != NULL))
/*     db_rnginclds(db_5ppos4feat(es->splcft),es->sqpt,db_3ppos4feat(es->splcft))) */
/* don't work with circular seqs */
        {
        nb = *(es->seq + es->sqpt - 1);
        if (es->cursens == DB_sens_comp)
          nb = (* bascmpfun)(nb,BAS_exact);
        }
      else
        nb = '\0';
      break;
    case DBS_allres:
    default:
      nb = *(es->seq + es->sqpt - 1);
      if (es->cursens == DB_sens_comp)
        nb = (* bascmpfun)(nb,BAS_exact);
      break;
    }     
else   /* check in case we have a circular seq */
  if ((segp = db_nxtseg4feat(es->curseg,es->splcft)) != NULL)
    {
    es->curseg = segp;
    es->cursens = es->curseg->sgsens;
    es->sqpt = db_5ppos4seg(es->curseg);
    return(db_getesresnoinc(es,segsonly,bascmpfun));
    }
  else
    nb = '\0';
return(nb);
}

char db_getesrescore(DB_ENTSTRCT *es,
                     DB_SPLICLMT segsonly,
                     void (* esincdecfun)(DB_ENTSTRCT *xs),
                     char (* bascmpfun)(char yb,
                                          BAS_MATMODE xmmod))
/* guts of the process to return the appropriate residue, irrespective any any complex links,
  use esincdecfun() to adjust the seq pointer. depending on segsonly, return
  null if there is not a valid segment or if pointer is out of feature limit.
  (* bascmpfun)() is used to create base complement */
{
char nb;

if ((nb = db_getesresnoinc(es,segsonly,bascmpfun)) != '\0')
  (* esincdecfun)(es);
return(nb);
}

char db_getesres(DB_ENTSTRCT *es,
                 DB_SPLICLMT segsonly,
                 void (* esincdecfun)(DB_ENTSTRCT *xs))
/* return the appropriate residue, irrespective any any complex links,
  use esincdecfun() to adjust the seq pointer. depending on segsonly, return
  null if there is not a valid segment or if pointer is out of feature limit */
{
return(db_getesrescore(es,segsonly,esincdecfun,ssd_bascmplmnt));
}

char db_getnxtesres(DB_ENTSTRCT *es,
                    DB_SPLICLMT segsonly,
                    char (* bascmpfun)(char yb,
                                         BAS_MATMODE xmmod))
/* return the appropriate residue, irrespective any any complex links,
  use esincdecfun() to adjust the seq pointer. if segsonly, then return
  null if there is not a valid segment. use (*bascmpfun)() for reverse
  strand residues */
{
return(db_getesrescore(es,segsonly,db_nxtesptr,bascmpfun));
}

char db_getprvesres(DB_ENTSTRCT *es,
                     DB_SPLICLMT segsonly,
                     char (* bascmpfun)(char yb,
                                          BAS_MATMODE xmmod))
/* return the previous residue, irrespective any any complex links,
  use esincdecfun() to adjust the seq pointer. if segsonly, then return
  null if there is not a valid segment. use (*bascmpfun)() for reverse
  strand residues */
{
return(db_getesrescore(es,segsonly,db_prvesptr,bascmpfun));
}

int db_getescodn(DB_ENTSTRCT *es,
                 DB_SPLICLMT splimt,
                 char *codn,
                 void (* esincdecfun)(DB_ENTSTRCT *xs))
/* use db_getesres() to fill a codon (codn assumed to be long enough) -
return the number of valid bases inserted */
{
int vcnt;
int cc;
char nr;
char *cp;

cc = vcnt = 0;
cp = codn;
while (cc++ < CODONLENGTH)
  {
  if ((nr = db_getesres(es,splimt,esincdecfun)) != '\0')
    vcnt++;
  bas_appchr(codn,&cp,nr,(CODONLENGTH+1));
  }
return(vcnt);
}

char db_nxtsplcres(DB_ENTSTRCT *es,
                   DB_SPLICLMT segsonly)
/* return the next base in sequence for presently-init-ed feature.
  if segsonly, then return '\0' if the position is not on segments */
{
char nb;
DB_SEGELT *segp;

if (es->splcft == NULL) /* not init-ed just return '\0' & do nothing else */
  return('\0');
else                    /* is position in feat range?? */
/*  if (!db_rnginclds(db_5ppos4feat(es->splcft),es->sqpt,
                      db_3ppos4feat(es->splcft))) */ /* don't work with circular seqs */
  {
  if (es->curseg != NULL)            /* have a defined segment, is it valid? */
    {
    if (db_rnginclds(db_5ppos4seg(es->curseg),es->sqpt,
                       db_3ppos4seg(es->curseg)))
    /* yep, in range, return this residue */
      return(db_getesres(es,segsonly,db_nxtesptr));
    else  /* not in range, try for another segment */
      {
      if ((segp = db_featseg4pos(es->splcft,es->sqpt)) != NULL)
        {
        es->curseg = segp;
        return(db_getesres(es,segsonly,db_nxtesptr));
        }
      else    /* not on a segment, is there another segment? */
        if ((segp = db_nxtseg4feat(es->curseg,es->splcft)) != NULL)
          {  /* yes, make it current and use it */
          es->curseg = segp;
          return(db_getesres(es,segsonly,db_nxtesptr));
          }
        else  /* no, set curseg NULL and return a base if OK */
          {
          es->curseg = NULL;
          return(db_getesres(es,segsonly,db_nxtesptr));
          }
      }
    }
  else  /* no current segment, try to find one */
    {
    es->curseg = db_featseg4pos(es->splcft,es->sqpt);
    if ((nb = db_getesres(es,segsonly,db_nxtesptr)) == '\0')
       db_nxtesptr(es);
    return(nb);
    }
  }
}

char db_prvsplcres(DB_ENTSTRCT *es,
                   DB_SPLICLMT segsonly)
/* return the previous base in sequence for presently-init-ed feature.
  if segsonly, then return '\0' if the position is not on segments */
{
char nb;
int p5extnt;
int p3extnt;

if (es->splcft == NULL) /* not init-ed just return '\0' & do nothing else */
  return('\0');
else                    /* is position in feat range?? */
  {
  db_featextremes(es,es->splcft,&p5extnt,&p3extnt);
  if (!db_rnginclds(p5extnt,es->sqpt,p3extnt))
/* no, return raw base, inc/dec ptr */
    {
    es->curseg = NULL;
    return(db_getesres(es,segsonly,db_prvesptr));
    }
  else               /* yes, is seg defined?? */
    if (es->curseg == NULL)
/* no, see if it should be */
      {
      if ((es->curseg = db_featseg4pos(es->splcft,es->sqpt)) != NULL)
        {
        es->cursens = es->curseg->sgsens;        /* make current sens that of segment */
        es->sqpt = db_3ppos4seg(es->curseg);     /* init pointer to 3' end */
        return(db_getesres(es,segsonly,db_prvesptr));
        }
      else
        return('\0');
      }
    else   /* seg defined, check we are in seg range  */
      if (!db_rnginclds(db_5ppos4seg(es->curseg),es->sqpt,
                          db_3ppos4seg(es->curseg)))
        {      /* out of this seg, get next seg, & init ptr */
        if ((es->curseg = db_prvseg4feat(es->curseg,es->splcft)) != NULL)
          es->sqpt = db_3ppos4seg(es->curseg);
        return(db_prvsplcres(es,segsonly));
        }
      else     /* in range, return this residue */
        return(db_getesres(es,segsonly,db_prvesptr));
  }
}

char db_prvsplcresnolmt(DB_ENTSTRCT *es,
                        DB_SPLICLMT segsonly)
/* return the previous base in sequence for presently-init-ed feature.
  if segsonly, then return '\0' if the position is not on segments.
The seq pointer will be decremented whether we are on sequence or
not */
{
char nb;

if (es->splcft == NULL) /* not init-ed just return '\0' & do nothing else */
  return('\0');
else                    /* is position in feat range?? */
  {
  if ((nb = db_prvsplcres(es,segsonly)) == '\0')
    es->sqpt = db_prvsqptr(es->sqpt,es->splcft->strctsens);
  return(nb);
  }
}

int db_segnoinfeat(DB_FEATSTRCT *fp,
                   DB_SEGELT *sep)
/* return the segment number from start of feature of sep.
0 if not observed in fp.  */
{
DB_SEGELT *sp;
int scnt;

sp = db_featseg4pos(fp,db_5ppos4feat(fp));
scnt = 1;
while (sp != NULL)
  if (sp == sep)
    return(scnt);
  else
    {
    sp = db_nxtseg4feat(sp,fp);
    scnt++;
    }
/* fell off end, didn'tt find sep */
return(0);
}

void db_cntshift(DB_ENTSTRCT *es,
                 int nshft,
                 char (* getbasfn)(DB_ENTSTRCT *xs,
                                   DB_SPLICLMT sgsonly))
/* move nshft residues in direction dicated by (* getbasfn)() - presently by 
iteratively calling it */
{
while (nshft-- > 0)
  (void) (*getbasfn)(es,DBS_segonly);
}

int db_cntlmtshift(DB_ENTSTRCT *es,
                   int nshft,
                   DB_SPLICLMT slmt,
                   char (* getbasfn)(DB_ENTSTRCT *xs,
                                       DB_SPLICLMT sgsonly))
/* move nshft residues in direction dicated by (* getbasfn)() - presently by 
iteratively calling it, Uses slmt to determine splice behaviour.  return
1 if shift was completed OK. */
{
int gotc;

gotc = 1;
while (gotc && (nshft-- > 0))
   gotc = (*getbasfn)(es,slmt) != '\0';
return(gotc);;
}

void db_cntnolmtshift(DB_ENTSTRCT *es,
                      int nshft,
                      DB_SPLICLMT slmt,
                      char (* getbasfn)(DB_ENTSTRCT *xs,
                                          DB_SPLICLMT sgsonly))
/* move nshft residues in direction dicated by (* getbasfn)() - presently by 
iteratively calling it. Does not limit pointer adjustment at seq limits.
Uses slmt to determine splice behaviour. */
{
while (nshft-- > 0)
   (void) (*getbasfn)(es,slmt);
}

void db_fwdshift(DB_ENTSTRCT *es,
                 int nshft)
/* iteratively call db_nxtsplcres to move nshft positions forward (5'->3') */
{
db_cntshift(es,nshft,db_nxtsplcres);
}

void db_bakshift(DB_ENTSTRCT *es,
                 int nshft)
/* iteratively call db_prvsplcres to move nshft positions back (3'->5') */
{
db_cntshift(es,nshft,db_prvsplcres);
}

int db_splcft2cdn(DB_ENTSTRCT *es,
                  DB_SPLICLMT segsonly,
                  char *cdn)
/* try to fill cdn with 3 successive calls to db_nxtsplcres().  Return
the number of non-null characters returned */
{
int vcnt;
char *cp;
int cc;
char nr;

cc = vcnt = 0;
cp = cdn;
while (cc++ < CODONLENGTH)
  {
  if ((nr = db_nxtsplcres(es,segsonly)) != '\0')
    vcnt++;
  bas_appchr(cdn,&cp,nr,CODONLENGTH+1);
  }
return(vcnt);
}

int db_rvsplcft2cdn(DB_ENTSTRCT *es,
                    DB_SPLICLMT segsonly,
                    char *cdn)
/* try to fill cdn with 3 successive calls to db_prvsplcres().  Return
the number of non-null characters returned */
{
int vcnt;
char *cp;
int cc;
char nr;

cc = vcnt = 0;
cp = cdn + CODONLENGTH - 1;
*(cdn + CODONLENGTH) = '\0';
while (cc++ < CODONLENGTH)
  {
  if ((nr = db_prvsplcres(es,segsonly)) != '\0')
    {
    vcnt++;
    *cp-- = nr;
    }
  else
    *cp-- = '?';
  }
return(vcnt);
}

void db_marksplicposn(DB_ENTSTRCT *esp,
                      DB_SPLICINFO *spinfo)
/* cache various esp values into spinfo to store
the present splice information for later restoration */
{
if ((esp != NULL) && (spinfo != NULL))
  {
  spinfo->cacheentp = esp;
  spinfo->cachefp = esp->splcft;
  spinfo->cachesegp = esp->curseg;
  spinfo->cachesqpt = esp->sqpt;
  spinfo->lim5p = esp->p5limt;
  spinfo->lim3p = esp->p3limt;
/*   spinfo->cachesens = esp->cursens; */
  }
}

int db_recovrsplicposn(DB_ENTSTRCT *esp,
                       DB_SPLICINFO *spinfo)
/* restore stored splice parameters from *spinfo to
esp.  Include sanity check that esp matches stored
entry pointer */
{
if ((esp != NULL) && (spinfo != NULL) && (esp == spinfo->cacheentp))
  {
  esp->splcft = spinfo->cachefp;
  esp->curseg = spinfo->cachesegp;
  esp->sqpt = spinfo->cachesqpt;
  esp->p5limt = spinfo->lim5p;
  esp->p3limt = spinfo->lim3p;
/*  esp->cursens = spinfo->cachesens; */
  return(1);
  }
else
  return(0);
}

int db_chkft4term(DB_ENTSTRCT *es,
                  DB_FEATSTRCT *fs,
                  int assrtnxt)
/* establish if the final codon of fs is a terminator: if not then on basis of
assrtnxt check the next codon and if it is, then correct fs value */
{
char codn[CODONLENGTH+1];
DB_SEGELT *finlseg;

db_init3psplic4feat(es,fs);
db_bakshift(es,(CODONLENGTH-1));
if (db_splcft2cdn(es,DBS_featlimit,&codn[0]) >= CODONLENGTH)
  if (sqt_trnslate(*(es->trnsmatrx),&codn[0]) == SQT_TERMCHR)    /* yep */
    {
    db_endfeatsplic(es);
    return(1);
    }
if ((assrtnxt) && (db_splcft2cdn(es,0,&codn[0]) >= CODONLENGTH) &&
      (sqt_trnslate(*(es->trnsmatrx),&codn[0]) == SQT_TERMCHR))   /* next is */
  {
  finlseg = db_3pseg4feat(fs);
  finlseg->sgstop = es->sqpt - 1;
  return(1);
  }
else
  return(0);
}

int db_ft2imatrx(DB_ENTSTRCT *es,
                 DB_FEATSTRCT *fs,
                 DB_SPLICLMT slmt,
                 TRANS_MATRX cnts)
/* count all codons for a feature into cnts,  return the number of valid
 codons recovered */
{
char cdn[CODONLENGTH+1];
int ccnt;
int cln;

ccnt = 0;
db_init5psplic4feat(es,fs);
fill_int_trnarr(cnts,0);
while ((cln = db_splcft2cdn(es,slmt,&cdn[0])) > 0)
  if (cln >= CODONLENGTH)
    {
    ccnt++;
    sqt_inc_cdncnt(cnts,&cdn[0]);
    }
return(ccnt);
}

int db_ft2imatrxrev(DB_ENTSTRCT *es,
                    DB_FEATSTRCT *fs,
                    DB_SPLICLMT slmt,
                    TRANS_MATRX cnts)
/* count all codons for a feature into cnts,  return the number of valid
 codons recovered */
{
char cdn[CODONLENGTH+1];
int ccnt;
int cln;

ccnt = 0;
db_init3psplic4feat(es,fs);
fill_int_trnarr(cnts,0);
while ((cln = db_rvsplcft2cdn(es,slmt,&cdn[0])) > 0)
  if (cln >= CODONLENGTH)
    {
    ccnt++;
    sqt_inc_cdncnt(cnts,&cdn[0]);
    }
return(ccnt);
}

int db_rescnt4imatrx(DB_ENTSTRCT *es,
                     char pres,
                     TRANS_MATRX cnts)
/* return the total number of pres codes for cnts, using the 
pre-built codon vector table of es */
{
int tcnt;
COD_VECT_STRCT *cvp;

tcnt = 0;
if ((cvp = sqt_getvect4res(es->cvecs,pres,NULL)) != NULL)
  while (cvp != NULL)
    {
    tcnt += cnts[cvp->codnvector[0]]
                  [cvp->codnvector[1]]
                  [cvp->codnvector[2]].resdata.iresidue;
    cvp = sqt_getvect4res(es->cvecs,pres,cvp);
    }
return(tcnt);
}

int db_termcnt4imatrx(DB_ENTSTRCT *es,
                      TRANS_MATRX cnts)
/* return the total number of terminator codes for cnts, using the 
pre-built codon vector table of es */
{
return(db_rescnt4imatrx(es,SQT_TERMCHR,cnts));
}

int db_tot4ress4imatrx(DB_ENTSTRCT *es,
                       char *rstrng,
                       TRANS_MATRX cnts)
/* total the counts for each valid residue in rstrng */
{
char *rp;
int tot;

tot = 0;
rp = rstrng;
while (*rp != '\0')
  {
  if (aares2int(*rp) > 0)
    tot += db_rescnt4imatrx(es,*rp,cnts);
  rp++;
  }
return(tot);
}

int db_totimtrx(DB_ENTSTRCT *es,
                TRANS_MATRX ccnts)
/* total the count of valid residues in ccnts, using genetic code defined
for es */
{
return(db_tot4ress4imatrx(es,"acdefghiklmnpqrstvwy*",ccnts));
}

char *db_olaptyp2str(DB_OLAPTYP ol)
  /* return a string for ol */
{
switch (ol)
  {
  case OLAP_incl:
    return("included");
    break;
  case OLAP_5pabut:
    return("5'abutted");
    break;
  case OLAP_3pabut:
    return("3'abutted");
    break;
  case OLAP_abuts:
    return("abutted");
    break;
  case OLAP_none:
  default:
    return("not-included");
    break;
  }
}

void db_allfeat2fl(FILE *fl,
                   DB_ENTSTRCT *es,
                   DB_FEATSTRCT *fs)
/* write all residues of fs out to fl. */
{
int nc;
int rcnt;

db_init5psplic4feat(es,fs);
rcnt = 0;
while ((nc = db_nxtsplcres(es,1)) != '\0')
  {
  fputc(nc,fl);
  rcnt++;
  }
fprintf(fl,"/%s#%d-%dres\n",db_ftkw2str(fs->featur),fs->featno,rcnt);
db_endfeatsplic(es);
}

void prt_datum(FILE *strm,
               TRN_MATELT *dat_ptr)
/* tell strm the value of *dat_ptr */
{
switch (dat_ptr->dtype)
  {
  case SQTDT_char:
    fputc(dat_ptr->resdata.trres.reschr,strm);
    break;
  case SQTDT_int:
    fprintf(strm,"%3d",dat_ptr->resdata.iresidue);
    break;
  case SQTDT_float:
    fprintf(strm,"%5.3f",dat_ptr->resdata.xresidue);
    break;
  case SQTDT_undefined:
  default:
    fprintf(stderr,"Invalid data type code: %d\n",dat_ptr->dtype);
    exit(1);
    break;
  }
}

void prt_mxd_tdata(FILE *strm,
                   TRANS_MATRX trnsar,
                   TRANS_MATRX datarr)
/* tell strm about trnsar contents in structured table for valid codons */
{
int p1,p2,p3;
int rp1,rp2,rp3;

for (p1 = 1; p1 <= MAXBASINT; p1++)
  for (p3 = 1; p3 <= MAXBASINT; p3++)
    {
    for (p2 = 1; p2 <= MAXBASINT; p2++)
      {
      if (p2 > 1)
        fprintf(strm,"   ");
      say_intcdn(strm,(rp1 = remapbasptr(p1)),(rp2 = remapbasptr(p2)),
                  (rp3 = remapbasptr(p3)));
      fprintf(strm," ");
      prt_datum(strm,&trnsar[rp1][rp2][rp3]);
      fprintf(strm," ");
      prt_datum(strm,&datarr[rp1][rp2][rp3]);
      }
    fputc('\n',strm);
    }
}

char *db_rngqual2str(DB_RNGQUAL rq)
  /* return a string for rq */
{
switch (rq)
  {
  case RQ_p5undefind:
    return("5'undefined");
    break;
  case RQ_p3undefind:
    return("3'undefined");
    break;
  case RQ_p5p3undefind:
    return("5'&3'undefined");
    break;
  case RQ_p5patchd:
    return("5'patched");
    break;
  case RQ_bothdefind:
    return("");          /* don't need to say anything */
    break;
  case RQ_unknown:
  default:
    return("Range Undefined");
    break;
  }
}

DB_RNGQUAL db_feat2rngqual(DB_FEATSTRCT *fp)
  /* look at range qualifiers of fp segment ends and return the classification
for them */
{
DB_SEGELT *p3sp;
DB_SEGELT *p5sp;

if (((p3sp = db_3pseg4feat(fp)) == NULL) ||
      ((p5sp = db_5pseg4feat(fp)) == NULL))
  return(RQ_unknown);
else
  if (((p3sp->gthan) && (p5sp->lthan)) || ((p3sp->lthan) && (p5sp->gthan)))
    return(RQ_p5p3undefind);
  else
    if (fp->strctsens == DB_sens_comp)
      if ((!p3sp->lthan) && (!p5sp->gthan))
        return(RQ_bothdefind);
      else
        if (p3sp->lthan)
          return(RQ_p3undefind);
        else
          return(RQ_p5undefind);
    else
      if ((!p3sp->gthan) && (!p5sp->lthan))
        return(RQ_bothdefind);
       else
         if (p5sp->lthan)
           return(RQ_p5undefind);
         else
           return(RQ_p3undefind);
}

int db_lmt4feat(DB_FEATSTRCT *fp,
                int (*maxminfn)(int x,
                                int y))
/* return the extent of fp as established by maxminfn(). */
{
return((*maxminfn)(db_5ppos4feat(fp),db_3ppos4feat(fp)));
}

int db_directeddst(DB_FEATSTRCT *f1,
                   DB_FEATSTRCT *f2,
                   int p3dst)
/* determine f1..f2 distance, depending on p3dst.  Return -ve for non-sensible
  configurations */
{
if (p3dst)
  switch (f1->strctsens)
    {
    case DB_sens_comp:
      return(db_3ppos4feat(f1) - db_lmt4feat(f2,imax));
      break;
    case DB_sens_norm:
    default:
      return(db_lmt4feat(f2,imin) - db_3ppos4feat(f1));
      break;
    }
else
  switch (f1->strctsens)
    {
    case DB_sens_comp:
      return(db_lmt4feat(f2,imin) - db_5ppos4feat(f1));
      break;
    case DB_sens_norm:
    default:
      return(db_5ppos4feat(f1) - db_lmt4feat(f2,imax));
      break;
    }
}

DB_FEATSTRCT *db_nxtfeat(DB_ENTSTRCT *es,
                         DB_FEATSTRCT *fp,
                         DB_FTYPELT *flst,
                         int p3dst)
/* return the next 3' or 5' feature to fp.  return NULL if no such feature 
 exists.  Scan entire featlst, so we don't make any assumptions about feature
 ordering. */
{
DB_FEATSTRCT *ff;
DB_FEATSTRCT *nxt;
int mindst;
int curdst;

nxt = NULL;
mindst = es->sqlen;
ff = es->featlst;
while (ff != NULL)
  {
  if (ff != fp)
    if (db_ptr4felt(flst,(int) ff->featur) != NULL)  /* is wanted sort */
      if (((curdst = db_directeddst(fp,ff,p3dst)) < mindst) && (curdst >= 0))
        {
        mindst = curdst;
        nxt = ff;
        }
  ff = ff->nextst;
  }
return(nxt);
}

int db_dst2nxtfeat(DB_ENTSTRCT *es,
                   DB_FEATSTRCT *fp,
                   DB_FTYPELT *flst,
                   int p3dst)
/* return the distance to the next 3' or 5' feature to fp.  return remaining 
 free sequence if no such feature  exists.  Scan entire featlst, so we 
 don't make any assumptions about feature ordering. */
{
DB_FEATSTRCT *nxt;

if ((nxt = db_nxtfeat(es,fp,flst,p3dst)) != NULL)
  return(db_directeddst(fp,nxt,p3dst));
else
  if (((p3dst) && (fp->strctsens == DB_sens_comp)) ||
       ((!p3dst) && (fp->strctsens != DB_sens_comp)))
    return(db_lmt4feat(fp,imin));
  else
    return(es->sqlen - db_lmt4feat(fp,imax));
}

int db_dist4pos2segend(DB_SEGELT *sep,
                       int pos,
                       int p3dist)
/* assuming pos is in range of sep, return
either the distance to 5' or 3' extremity
depending on p3dist inclusively.  return -1 if
pos is not in seg */
{
int p3end;
int p5end;

p3end = db_3ppos4seg(sep);
p5end = db_5ppos4seg(sep);
if (db_rnginclds(p5end,pos,p3end))
  if (sep->sgsens == DB_sens_comp)
    if (p3dist)
      return(pos - p3end + 1);
    else
      return(p5end - pos + 1);
  else
     if (p3dist)
       return(p3end - pos + 1);
     else
       return(pos - p5end + 1);
else
  return(-1);
}

int db_splicdist4feat2pos(DB_FEATSTRCT *fp,
                          int spos,
                          int p3dist)
/* return spliced distance from 5' or 3' end of 
feat fp to spos.  -1 if spos not on segment */
{
DB_SEGELT *poseg;
DB_SEGELT *sp;
int dist;

if ((poseg = db_featseg4pos(fp,spos)) == NULL)
  return(-1);
else
  {
  dist = 0;
  if (p3dist)
    sp = db_featseg4pos(fp,db_3ppos4feat(fp));
  else
    sp = db_featseg4pos(fp,db_5ppos4feat(fp));
  while (sp != NULL)
    {
    if (sp == poseg)
      return(dist + db_dist4pos2segend(sp,spos,p3dist));
    else
      dist += sp->sgstop - sp->sgstart + 1;
    if (p3dist)
      sp = db_prvseg4feat(sp,fp);
    else
      sp = db_nxtseg4feat(sp,fp);
    }
/* fell off end, return -1 */
  return(-1);
  }
}

void db_sayafeatcore(FILE *fl,
                     DB_ENTSTRCT *es,
                     DB_FEATSTRCT *fs,
                     int showseq)
/* tell fl about fs contents. showseq controls
whether any sequence is to be shown */
{
DB_SEGELT *sp;
int scnt;
DB_TBLINFO *ip;
int ltot;
int rcnt;
char rs;

if (fs == NULL)
  fprintf(fl,"NULL feat\n");
else
  {
  fprintf(fl,"%s#%d: ",db_ftkw2str(fs->featur),fs->featno);
  if (db_crossplccnt(fs) > 0)
    fputs("\"Cross Entry Splice\"",fl);
  else
    {
    if (showseq)
      {
      db_init5psplic4feat(es,fs);
      fputc('"',fl);
      for (rcnt = 0; rcnt < 30; rcnt++)
        if ((rs = db_nxtsplcres(es,1)) == '\0')
          fputc('?',fl);
        else
          fputc(rs,fl);
      fputc('"',fl);
      db_endfeatsplic(es);
      }
    }
  sp = db_5pseg4feat(fs);
  while (sp != NULL)
    {
    db_saysegposns(fl,sp);
    if ((sp = db_nxtseg4feat(sp,fs)) != NULL)
      fputc(',',fl);
    }
  ltot = db_featlength(fs,&scnt);
  fprintf(fl," %s; %d%s %d seg%s; ",
            (es->efmt==DBFMT_swiss?"":db_sens2str(fs->strctsens)),
            ltot,(es->efmt==DBFMT_swiss?"res":"bp"),scnt,(scnt!=1?"s":""));
  ip = fs->infolist;
  while (ip != NULL)
    {
    db_sayqulval(fl,ip,";");
    ip = ip->nxtielt;
    }
  if (fs->mRNA != NULL)
    {
    fprintf(fl,"-%s by:\n  ",db_olaptyp2str(db_chkfeatolap(fs->mRNA,fs)));
    db_sayafeatcore(fl,es,fs->mRNA,showseq);
    }
  else
    fputc('\n',fl);
  }
}  

void db_sayafeat(FILE *fl,
                 DB_ENTSTRCT *es,
                 DB_FEATSTRCT *fs)
/* tell fl about fs contents. */
{
db_sayafeatcore(fl,es,fs,1);
}

DB_SENS db_chkfeatsens(DB_FEATSTRCT *fs)
  /* scan *fs for the features, return the value */
{
int scnt[DB_sens_mixd+1];
DB_SEGELT *ep;
DB_SENS sp;
int big;
DB_SENS sbig;

if (fs == NULL)
  return(DB_sens_unk);
else
  {
  for (sp = DB_sens_unk; sp <= DB_sens_mixd; sp++)
    scnt[sp] = 0;
  ep = fs->fstseg;
  while (ep != NULL)
    {
    scnt[ep->sgsens]++;
    ep = ep->nextseg;
    }
  if ((scnt[DB_sens_norm] > 0) && (scnt[DB_sens_comp] > 0))
    return(DB_sens_mixd);
  else
    {
    big = 0;
    sbig = DB_sens_unk;
    for (sp = DB_sens_unk; sp <= DB_sens_mixd; sp++)
      if (scnt[sp] > big)
        {
        big = scnt[sp];
        sbig = sp;
        }
     return(sbig);
     }
  }
}

char *db_strafter(char *slin,
                  char dchr)
/* check for dchr in slin, return the next character, NULL if dchr is not
found */
{
char *sp;

if ((sp = index(slin,dchr)) != NULL)
  return(sp+1);
else
  return(NULL);
}

char *db_strdup_quotd(DB_ENTSTRCT *es)
  /* look on current line for a quoted string, produce a strdup-ed copy of it
     up to DBLU_MAXSRCLEN chars.  Return NULL if no '"' occurred on 
     current line  */
{
char *np;
char cbuf[DBLU_MAXSRCLEN+1];
char *dp;

if ((np = db_strafter(es->cp,'"')) == NULL)
  return(NULL);
else
  {
  dp = &cbuf[0];
  while (*np != '"')
    {
    if (*np == '\0')        /* end of line */
      {
      es->cp = NULL;
      if (db_getesln(es) == NULL)    /* at end of seq, return what we have */
        {
        bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
        return(bas_strdup(&cbuf[0]));
        }
      else
        if ((db_actlintyp(es,db_classfyentln(es)) != DBLC_feature) ||
              (db_clssfeatln(es) != FTKW_unknown))
          {
          bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
          es->cp = es->slbuf;
          return(bas_strdup(&cbuf[0]));
          }
        else
          {
          np = es->cp = es->slbuf + DB_FEATLOCATN - 1;
          if (db_isqual(es->slbuf,es->slblen,es->cp,es->fqulwlu) != 
                FTQU_unknown)   /* new qualifier, so stop */
            {
            bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
            es->cp = es->slbuf;
            return(bas_strdup(&cbuf[0]));
            }
          }
      }
    else
      bas_appchr(&cbuf[0],&dp,*np++,DBLU_MAXSRCLEN);
    }
  bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
  es->cp = NULL;
  return(bas_strdup(&cbuf[0]));
  }
}

void db_mltdup_quotd(DB_ENTSTRCT *es,
                     DB_FEATSTRCT *fp,
                     DBLU_FTQUAL qul)
/* look on current line for a quoted string, append each such line to
fp info list. */
{
char *np;
char cbuf[DBLU_MAXSRCLEN+1];
char *dp;
int scn;

if ((np = db_strafter(es->cp,'"')) != NULL)
  {
  dp = &cbuf[0];
  scn = 1;
  while ((scn) &&(*np != '"'))
    {
    if (*np == '\0')        /* end of line */
      {
      bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
      fp->infoend = db_appnielt(&fp->infolist,qul,bas_strdup(&cbuf[0]));
      dp = &cbuf[0];
      es->cp = NULL;
      if (db_getesln(es) == NULL)
        {
        scn = 0;
        if (dp > &cbuf[0])     /* at end of seq, store what we have */
          {
          bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
          fp->infoend = db_appnielt(&fp->infolist,qul,bas_strdup(&cbuf[0]));
          }
        }
      else
        if ((db_actlintyp(es,db_classfyentln(es)) != DBLC_feature) ||
              (db_clssfeatln(es) != FTKW_unknown))
          {
          scn = 0;
          if (dp > &cbuf[0])
            {
            bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
            fp->infoend = db_appnielt(&fp->infolist,qul,bas_strdup(&cbuf[0]));
            }
          es->cp = es->slbuf;
          }
        else                /* continue */
          {
          np = es->cp = es->slbuf + DB_FEATLOCATN - 1;
          if (db_isqual(es->slbuf,es->slblen,es->cp,es->fqulwlu) != 
                FTQU_unknown)   /* new qualifier, so stop */
            {
            scn = 0;
            if (dp > &cbuf[0])
              {
              bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
              fp->infoend = db_appnielt(&fp->infolist,qul,
                                          bas_strdup(&cbuf[0]));
              }
            es->cp = es->slbuf;
            }
          }
      }
    else
      bas_appchr(&cbuf[0],&dp,*np++,DBLU_MAXSRCLEN);
    }
  if (dp > &cbuf[0])
    {
    bas_appchr(&cbuf[0],&dp,'\0',DBLU_MAXSRCLEN);
    fp->infoend = db_appnielt(&fp->infolist,qul,bas_strdup(&cbuf[0]));
    }    
  es->cp = NULL;
  scn = 0;
  }
}

DB_WANTQU *db_wntquptr4featnqul(DB_FTYPELT *flst,
                                DBLU_FTKW ftkw,
                                DBLU_FTQUAL fq)
/* return a pointer to the element matched by ftkw & fq, NULL otherwise */
{
DB_FTYPELT *fep;

if ((fep = db_ptr4felt(flst,ftkw)) != NULL)
  return(db_ptr4qualval(fep->wquals,fq));
else
  return(NULL);
}

void db_chkfqual(DB_ENTSTRCT *es,
                 DB_FEATSTRCT *fp,
                 DBLU_FTQUAL qul)
/* check this qul value out - save data or otherwise as required */
{
char *np;
int tvl;

#ifdef DB_DEBUG

fprintf(stdout,"qual: %s: %s\n",db_ftqual2str(qul),es->cp);

#endif
es->onqual = 1;
switch (qul)
  {
  case FTQU_codon_start:
    if ((np = db_strafter(es->cp,'=')) != NULL)
      fp->cod_start = atoi(np);
    es->cp = NULL;
    break;
  case FTQU_pseudo:
    fp->pseudogn = 1;
    es->cp = NULL;
    break;
  case FTQU_partial:
    fp->partcds = 1;
    es->cp = NULL;
    break;
  case FTQU_transl_table:
    if ((np = db_strafter(es->cp,'=')) != NULL)
      {
      tvl = atoi(np);
      fp->gencod = sqt_int2gbgencod(tvl);
      }
    es->cp = NULL;
    break;
  case FTQU_note:
  case FTQU_product:
    if (db_wntquptr4featnqul(es->intrst,fp->featur,qul) != NULL)
      db_mltdup_quotd(es,fp,qul);
    break;
  case FTQU_db_xref:
  case FTQU_EC_number:
  case FTQU_gene:
  case FTQU_function_qu:
  case FTQU_anticodon:
  case FTQU_bound_moiety:
  case FTQU_citation:
  case FTQU_codon:
  case FTQU_cons_splice:
  case FTQU_direction:
  case FTQU_evidence:
  case FTQU_frequency:
  case FTQU_label_qu:
  case FTQU_mod_base:
  case FTQU_number:
  case FTQU_organism:
  case FTQU_phenotype:
  case FTQU_protid:
  case FTQU_rpt_family:
  case FTQU_rpt_type:
  case FTQU_rpt_unit:
  case FTQU_standard_name:
  case FTQU_transl_except:
  case FTQU_type_qu:
  case FTQU_usedin:
  case FTQU_unknown:
  case FTQU_name:
  default:
    if (db_wntquptr4featnqul(es->intrst,fp->featur,qul) != NULL)
      fp->infoend = db_appnielt(&fp->infolist,qul,db_strdup_quotd(es));
    break;
  }
}

DB_FTYPELT *db_str2lst(char *ustr)
  /* scan ustr for a list of seq no integers - return pointer to a linked list
of values specified therein.
  Input can be a string of integers, separated by non-digit chars.
  n..m or n-m is interpreted as all ints from n to m inclusive */
{
DB_FTYPELT *ilst;
char *sp;
int newvl;
DB_FTYPELT *nxt;
char *tp;
int dcnt;
int prval;
int rpt;

ilst = NULL;
sp = ustr;
dcnt = 0;
while (*sp != '\0')
  {
  newvl = (int) strtol(sp,&tp,10);
  if (dcnt >= 2)       /* range of values, add each in turn */
    {
    if (newvl > (prval+1))
      for (rpt = (prval+1); rpt <= newvl; rpt++)
         nxt = db_appnfelt(&ilst,rpt);
    else
      for (rpt = newvl; rpt <= (prval-1); rpt++)
         nxt = db_appnfelt(&ilst,rpt);
    dcnt = 0;
    }
  else  /* single value */
    nxt = db_appnfelt(&ilst,newvl);
  prval = newvl;
  sp = tp;
  while ((*sp != '\0') && (index("0123456789",*sp) == NULL))
    {
    if (*sp == '.')
      dcnt++;
    else
      if (*sp == '-')
        dcnt = 2;
      else
        dcnt = 0;
    sp++;
    }
  }
return(ilst);
}  

int db_str2locval(char *str,
                  DB_LOCQUAL *lq)
/* read an integer value from str.  If starts or ends with < or >, set *lq
  accordingly */
{
if (str == NULL)
  {
  if (lq != NULL)
    *lq = DB_rq_mt;
  return(0);
  }
else
  switch (*str)
    {
    case '>':
      if (lq != NULL)
        *lq = DB_rq_gthan;
      return(db_str2locval((str+1),NULL));
      break;
    case '<':
      if (lq != NULL)
        *lq = DB_rq_lthan;
      return(db_str2locval((str+1),NULL));
      break;
    default:
      if (lq != NULL)
        *lq = DB_rq_none;
      return(atoi(str));
      break;
    }
}

void db_chklocnval(DB_ENTSTRCT *es,
                   char *tkn,
                   DB_FEATSTRCT *flst,
                   char *xidnam)
/* tkn should contain a position specifier: decode it and set appropriately
in es fields */
{
int vl;
DB_LOCQUAL lq;

switch (db_tokclass(tkn))
  {
  case DB_toknum:  /* location, unpick */
    if (!es->onqual)         /* ignore once qualifiers start */
      {
      db_ftok2str(stdout,FLOC_nvalue,tkn);
      vl = db_str2locval(tkn,&lq);
      if (es->curseg == NULL)      /* start of new seg */
        {
        es->curseg = flst->lstseg =
          db_appnseg(&flst->fstseg,vl,vl,es->cursens,xidnam);
        es->curseg->lthan = lq == DB_rq_lthan;
        }
      else
        {
        es->curseg->sgstop = vl;
        es->curseg->gthan = lq == DB_rq_gthan;
        es->curseg = NULL;       /* finished this value */
        }
      }
    *tkn = '\0';         /* force that this is empty */
    break;
  case DB_tokalph:
  case DB_tokmixd:
  case DB_tokunk:
  default:        /* other, ignore */
    break;
  }
}

DBLU_FTKW db_scnnafeat(DB_ENTSTRCT *es)
/* attempt to identify this line as a NA feature descriptor.  if wanted, append
it to es->featlst.  In any event, return the feature identified */
{
DB_FTYPELT *lte;
DBLU_FTKW feat;
DB_FEATSTRCT *flst;
char tbuf[DBLU_MAXSRCLEN+1];
DBLU_FLOCTYPE floc;
int plevl;       /* parenthesis level */
DBLU_FTQUAL qul; /* qualifier type of this line */
int scn;
int vl;
char *xidnam;    /* if find an id xref, store the name here */
DB_DLMTTYPE toktp;
DB_TOKTYPE tokcls;

if ((feat = db_clssfeatln(es)) != FTKW_unknown)
  es->nseen++;
lte = db_ptr4felt(es->intrst,feat);
if (lte != NULL)    /* did want this one */
  {
  flst = es->lastfeat = db_appnfeat(&es->featlst,feat,es->ename);
  flst->featno = es->nseen;
  flst->savdno = ++es->nfeats;
  es->cp = es->slbuf + DB_FEATLOCATN - 1;
  plevl = 0;
  es->curseg = NULL;
  es->cursens = DB_sens_norm;
  es->onqual = 0;
  scn = 1;
  xidnam = NULL;
  tbuf[0] = '\0';
  while (scn)
    {
    if (!es->onqual)
      {
      toktp = db_gettok(es->slbuf,&es->cp,&tbuf[0],DBLU_MAXSRCLEN);
      if (strlen(&tbuf[0]) > 0)  /* did get a token: check it */
        switch (tokcls = db_tokclass(&tbuf[0]))
          {
          case DB_toknum:
            floc = FLOC_nvalue;
            break;
          case DB_tokalph:       /* join, complement, etc. */
            if ((floc = wlu_chkwrd(es->flocwlu,&tbuf[0])) == FLOC_complement)
              {
              db_ftok2str(stdout,floc,&tbuf[0]);
              es->curseg = NULL;       /* finished this value */
              switch (es->cursens)
                {
                case DB_sens_comp:
                  flst->strctsens = es->cursens = DB_sens_norm;
                  break;
                case DB_sens_norm:
                case DB_sens_mixd:
                case DB_sens_unk:
                default:
                  flst->strctsens = es->cursens = DB_sens_comp;
                  break;
                }
              }
            break;
          case DB_tok_mt:
            break;
          case DB_tokmixd:
          case DB_tokunk:
          default:
            db_ftok2str(stdout,FLOC_unknown,&tbuf[0]);
            break;
          }
      }
    else
      {
      tokcls = DB_tok_mt;
      toktp = DB_dlmt_eol;
      }
    switch (toktp)
      {
      case DB_dlmt_eol:             /* EOL */
        if (tokcls == DB_toknum) /* have a no, finish processing it */
          db_chklocnval(es,&tbuf[0],flst,xidnam);
        es->cp = NULL;       /* force new line */
        if (plevl <= 0)      /* should have finished this location scan */
          flst->strctsens = db_chkfeatsens(flst);
        if (db_getesln(es) == NULL)
          scn = 0;
        else
          if ((db_actlintyp(es,db_classfyentln(es)) != DBLC_feature) ||
                (db_clssfeatln(es) != FTKW_unknown))
            {   /* on a new feature or linetype:
                 allow to be processed elsewhere */
            es->cp = es->slbuf;  /* can process all this line elsewhere */
            scn = 0;
            }
          else      /* not a new feat - check this line... */
            {
            es->cp = es->slbuf + DB_FEATLOCATN - 1;
            if ((qul = db_isqual(es->slbuf,es->slblen,es->cp,es->fqulwlu))
                   != FTQU_unknown)
              db_chkfqual(es,flst,qul);
            }
        break;
      case DB_dlmt_lpar:       /* '(' */
        es->curseg = NULL;       /* finished this value */
        db_ftok2str(stdout,floc,&tbuf[0]);
        plevl++;
        break;
      case DB_dlmt_rpar:           /* ')' */
        db_ftok2str(stdout,floc,&tbuf[0]);
        db_chklocnval(es,&tbuf[0],flst,xidnam);
        es->curseg = NULL;       /* finished this value */
        xidnam = NULL;
        plevl--;
        break;
      case DB_dlmt_comma:          /* ',' */
        db_chklocnval(es,&tbuf[0],flst,xidnam);
        es->curseg = NULL;       /* finished this value */
        xidnam = NULL;
        break;
      case DB_dlmt_colon:          /* ':' */
        xidnam = bas_strdup(&tbuf[0]); /* malloc() mem for this id here */
        break;
      case DB_dlmt_2dots:          /* '..' */
        db_ftok2str(stdout,floc,&tbuf[0]);
        db_chklocnval(es,&tbuf[0],flst,xidnam);
        break;
      case DB_dlmt_not:            /* Not a recognised delimiter */
      default:
        db_chklocnval(es,&tbuf[0],flst,xidnam);
        break;
      }
    }
  }
else
  flst = NULL;
if (flst == NULL)
  es->cp = NULL;
return(feat);
}

DBLU_FTKW db_classnchkftline(DB_ENTSTRCT *es,
                             char *line,
                             int *rstart,
                             int *rstop,
                             DBLU_FTQUAL *qul,
                             char **cmt)
/* using parameters embedded in es, examine line
as a feature table line (Uniprot fmt), returning
rstart & rstop for the positions, and checking that
this feature is of interest;
qualifier values.  The entire comment/qualifier
field is returned.  If the line appears to be a
continuation, then the value FTKW_continuation
is returned. */
{
char qualstr[DB_SWFTPOS1COL - DB_SWSFKWCOL + 2];
int qslmt;
DBLU_FTKW feat;
DBLU_FTQUAL qual;
int slen;
DB_FTYPELT *elt;

qslmt = DB_SWFTPOS1COL - DB_SWSFKWCOL + 2;
bzero(&qualstr[0],qslmt);
if (strncmp("         ",line+DB_FEATINSET-1,9) == 0) /* continuation?? */
  {
  if (cmt != NULL)
    *cmt = line + DB_SWFTQUALCOL - 1;
  return(FTKW_continuation);
  }
else
  {
  slen = (int) strlen(line);
  feat = db_classfixdfld(line,slen,line+DB_FEATINSET-1,es->fkwlu," \t");
  if ((elt = db_ptr4felt(es->intrst,(int) feat)) != NULL) /* want this sort of feature */
    {
    qual = (DBLU_FTQUAL) db_classfixdfld(line,slen,line+DB_SWFTQUALCOL-1,es->fqulwlu," \t.;:");
    if ((elt->wquals == NULL) || (qual != FTQU_unknown))  /* want this one */
      {
      if (rstart != NULL)
        *rstart = atoi(line+DB_SWFTPOS1COL-1);
      if (rstop != NULL)
        *rstop = atoi(line+DB_SWFTPOS2COL-1);
      if (qul != NULL)
        *qul = qual;
      if (cmt != NULL)
        *cmt = line + DB_SWFTQUALCOL - 1;
      return(feat);
      }
    }
  }
return(FTKW_unknown);    
}
                         

DBLU_FTKW db_scnswfeat(DB_ENTSTRCT *es)
/* attempt to identify this line as a SwissProt feature descriptor.
if wanted, append it to es->featlst.  In any event, return the feature 
identified */
{
DBLU_FTKW feat;
DB_FEATSTRCT *flst;
int scn;
int vl;
int fndlt;
int frompos;
int topos;
DBLU_FTQUAL qual;
char *cmt;
DBLU_FTQUAL newqual;
DBLU_FTKW newfeat;

scn = 1;
flst = NULL;
feat = FTKW_unknown;
while (scn)
  {
  if ((newfeat = db_classnchkftline(es,es->slbuf,&frompos,&topos,&qual,&cmt))
         == FTKW_unknown)
    {
    scn = 0;
    es->cp = NULL;
    }
  else
    {
    if (newfeat == FTKW_continuation)
      {
      if (flst != NULL)  /* appnd cont info */
        {
        flst->infoend = db_appnielt(&flst->infolist,qual,bas_strdup(cmt));
        es->cp = NULL;
        }
      else   /* cont line, but not one we want */
        {
        scn = 0;
        es->cp = NULL;
        }
      }
    else /* new wanted feat */
      {
      es->nseen++;
      flst = es->lastfeat = db_appnfeat(&es->featlst,feat,es->ename);
      flst->featno = es->nseen;
      flst->savdno = es->nfeats;
      flst->strctsens = DB_sens_norm;
      flst->lstseg = db_appnseg(&flst->fstseg,frompos,topos,DB_sens_norm,NULL);
      flst->infoend = db_appnielt(&flst->infolist,qual,bas_strdup(cmt));
      feat = newfeat;
      es->cp = NULL;
      }
    if (db_getesln(es) == NULL)
      scn = 0;
    }
  }
return(feat);
}

void db_sqlin2sq(DB_ENTSTRCT *es)
  /* read contents of current line as a sequence line, append to sq buffer */
{
char nc;

switch (es->efmt)
  {
  case DBFMT_embl:
  case DBFMT_swiss:
  case DBFMT_embl_old:
  case DBFMT_sqmonk:
    es->cp = es->slbuf + DB_MBLSQINSET - 1;
    break;
  case DBFMT_genbank:
  default:
    es->cp = es->slbuf + DB_GBSQINSET - 1;
    break;
  }
while ((nc = *(es->cp)) != '\0')
  {
  if (index(" 0123456789",nc) == NULL)  /* reject blanks,digits */
    bas_appchr(es->seq,&es->sfpt,nc,es->sqlen);
  es->cp++;
  }
es->cp = NULL;
}

void db_dispos_feat(DB_FEATSTRCT **fp)
  /* lose storage associated with *fp */
{
if (*fp != NULL)
  {
  db_killsegs4featstrct(*fp);
  db_killtblinfolst(&(*fp)->infolist,&(*fp)->infoend);
  db_killstrlst(&(*fp)->enamelist,&(*fp)->enamlstend);
  if ((*fp)->parentname != NULL)
    memfree((*fp)->parentname);
  }
}

void db_delfeatelt(DB_FEATSTRCT *sp,
                   DB_FEATSTRCT **lstrt,
                   DB_FEATSTRCT **lend)
/* delete sp from list *lstrt..*lend */
{
DB_FEATSTRCT *pt;

if (sp != NULL)
  {
  if ((pt = sp->prevst) == NULL)
    *lstrt = sp->nextst;
  else
    pt->nextst = sp->prevst;
  if ((pt = sp->nextst) != NULL)
    pt->prevst = sp->prevst;
  else
    *lend = sp->prevst;
  db_dispos_feat(&sp);
  memfree(sp);
  }
}

void db_killfeatst4entst(DB_ENTSTRCT *fs)
  /* iteratively delete the segment list for fs - memory is freed */
{
while (fs->featlst != NULL)
  db_delfeatelt(fs->featlst,&fs->featlst,&fs->lastfeat);
}

void db_dispos_ent(DB_ENTSTRCT **es)
  /* lose all allocated storage associated with *es and set to NULL */
{
DB_FEATSTRCT *fp;

if (*es != NULL)
  {
  db_freeifok((*es)->ename);
  db_freeifok((*es)->seq);
  db_freeifok((*es)->date);
  db_freeifok((*es)->nid);
  db_freeifok((*es)->descr);
  wlu_clralllustrct((*es)->nawlu);
  wlu_clralllustrct((*es)->oddwlu);
  wlu_clralllustrct((*es)->strndwlu);
  wlu_clralllustrct((*es)->fkwlu);
  wlu_clralllustrct((*es)->flocwlu);
  wlu_clralllustrct((*es)->fqulwlu);
  db_killacclst(&(*es)->acclst);
/*  db_killfeltlst(&(*es)->intrst); */ /* don't - this is now imported */
  if ((*es)->freefeatids)
    {
    fp = (*es)->featlst;
    while (fp != NULL)
      {
      if (fp->idp != NULL)
        memfree(fp->idp);
      fp = fp->nextst;
      }
    }
  db_killfeatst4entst(*es);
  db_freeifok((*es)->slbuf);
  db_freeifok((char *)(*es)->trnsmatrx);
  db_killlclsslst(&(*es)->lcinfo);
  if ((*es)->cvecs != NULL)
    {
    sqt_clrvectarr((*es)->cvecs);
    db_freeifok((char *)(*es)->cvecs);
    }
  memfree(*es);
  *es = NULL;
  }
}

int db_readsqlen(DB_ENTSTRCT *es)
  /* return the actual read sequence length */
{
return(es->sfpt - es->seq);
}

size_t db_strlen(char *cp)
  /* return strlen for cp, 0 for NULL cp */
{
return((cp == NULL)?(size_t) 0:strlen(cp));
}

char *db_linfo_offset(DB_ENTSTRCT *esp,
                      DBLU_DBFMT fmt)
/* return the offset of the actual information in the new
line buffer, depending on fmt */
{
return(esp->slbuf + db_datoffst(fmt));
}

DB_ENTSTRCT *db_parsemain(FILE *df,
                          WRD_LUSTRCT *lclass,
                          DBLU_DBFMT fmt,
                          DB_FTYPELT *ufwrds,
                          int ndxrun,
                          DB_LCINFO *linflst)
/* now the real guts of the parsing process.  This version
accepts a linked list of lineclass data which, if non-NULL, will
be used to capture information from relevant lines */
/* scan df from current position to next sq terminator.  wlus is used for
key words.  Return ptr to entry structure, NULL if no ident found.
scan string intrst and note what is to be kept from feature table.
ndxrun flags if it is only an indexing run or whether it is necessary to
create storage for and read the sequence data  */
{
int scanon;
DB_FEATSTRCT *newft;      /* new feature being handled */
DB_ENTSTRCT *newent;      /* new entry record */
WRD_LUSTRCT *ftkw;        /* feature table keywords */
long int coffs;           /* offset of start of this file */
int needoffs;             /* flag to note that offset is needed */
char *bp;                 /* buffer pointer */
int gotid;                /* set if id seen */
DBLU_FTKW feat;           /* "gene" feature if required */

needoffs = ndxrun;
ftkw = db_getkwstrct(fmt);
newent = (DB_ENTSTRCT *) getmemory(sizeof(DB_ENTSTRCT),"Entry structure");
db_initentstrct(newent,fmt,ftkw);
newent->intrst = ufwrds;
newent->lclss = lclass;
newent->sfl = df;
newent->slbuf = (char *) getmemory((DBLU_MAXSRCLEN+1),"Source buffer");
newent->slblen = DBLU_MAXSRCLEN;
newent->cp = NULL;
newent->slcnt = 1;
scanon = 1;
gotid = 0;
bp = NULL;
while (scanon)
  {
  if (needoffs)
    coffs = ftell(df);
  if (db_getesln(newent) == NULL)      /* EOF */
    {    
    scanon = 0;
    if (newent->prvltyp != DBLC_sequence)
      {
      if (newent->prvltyp != DBLC_unk)
        db_logerr("Unexpected EOF",newent->slcnt,NULL);
      db_dispos_ent(&newent);
      return(NULL);
      }
    }
  else
    {
/* fprintf(stdout,"%d\t%s\n",newent->slcnt,newent->slbuf); */
    switch (newent->curltyp = db_classfyentln(newent))
      {
      case DBLC_ident:
        db_unpickidline(newent->slbuf,newent,fmt);
        if (newent->sqlen <= 0)
          scanon = db_logerr("Invalid length in ID line",newent->slcnt,
                                newent->slbuf);
        else
          if (!ndxrun)
            newent->sfpt = newent->seq =
              (char *) getmemory((newent->sqlen+1),"Entry Seq buffer");
        newent->eoffst = (unsigned long) coffs;
        needoffs = 0;            /* don't keep calling ftell() */
        newent->prvltyp = newent->curltyp;
        newent->cp = NULL;
        gotid = 1;
        break;
      case DBLC_access:
        newent->cp = db_linfo_offset(newent,newent->efmt);
        db_unpickaccno(newent->cp,&newent->acclst);
        newent->prvltyp = newent->curltyp;
        newent->cp = NULL;
        break;
      case DBLC_naident:
      case DBLC_keyword:
      case DBLC_source:
      case DBLC_taxonomy:
      case DBLC_reference:
      case DBLC_author:
      case DBLC_title:
      case DBLC_journal:
      case DBLC_remark:
      case DBLC_basecount:
      case DBLC_ignore:
      case DBLC_seqversn:
      case DBLC_refnumbr:
      case DBLC_refpositn:
      case DBLC_xref:
      case DBLC_dbxref:
        if (db_lcent4udata(linflst,newent->curltyp,NULL) != NULL)
          db_appnlcelt(&newent->lcinfo,newent->curltyp,bas_strdup(db_linfo_offset(newent,fmt)));
        newent->cp = NULL;
        newent->prvltyp = newent->curltyp;
        break;
      case DBLC_description:
        if (newent->descr == NULL)
          bp = db_linfo_offset(newent,fmt);
        newent->descr = db_strdupskip(bp,db_strlen(bp),"");
        newent->prvltyp = newent->curltyp;
        newent->cp = NULL;
        break;
      case DBLC_datestamp:
        if (((fmt == DBFMT_embl) || (fmt == DBFMT_swiss)
                || (fmt == DBFMT_embl_old) || (fmt == DBFMT_sqmonk)) &&
             (newent->date == NULL))  /* only take first date entry */
          newent->date = db_strdupskip(db_linfo_offset(newent,fmt),
                                         db_strlen(bp)," \t");
        newent->prvltyp = newent->curltyp;
        newent->cp = NULL;
        break;
      case DBLC_sequence:
        newent->cp = NULL;
        newent->prvltyp = newent->curltyp;
        break;
      case DBLC_feature:
        switch (fmt)
          {
          case DBFMT_embl:
          case DBFMT_embl_old:
          case DBFMT_sqmonk:
            if (ndxrun)
              newent->cp = NULL;
            else
              (void) db_scnnafeat(newent);
            break;
          case DBFMT_swiss:
            if (ndxrun)
              newent->cp = NULL;
            else
              (void) db_scnswfeat(newent);
            break;
          case DBFMT_genbank: /* genbank doesn't have info on "feature" line */
          default:
            newent->cp = NULL;
          }
        newent->prvltyp = newent->curltyp;
        break;
      case DBLC_terminator:
        newent->cp = NULL;
        scanon = 0;
        break;
      case DBLC_unk:
      default:
        switch (newent->efmt)
          {
          case DBFMT_embl:
          case DBFMT_embl_old:
          case DBFMT_swiss:
          case DBFMT_sqmonk:
            switch (newent->prvltyp)
              {
              case DBLC_sequence:
                if (!ndxrun)
                  db_sqlin2sq(newent);
                else
                  newent->cp = NULL;
                break;
              default:                 /* ignore */
                newent->cp = NULL;
                break;
              } 
            break;
          case DBFMT_genbank:
          default:
            switch (newent->prvltyp)
              {
              case DBLC_access:
                newent->cp = db_linfo_offset(newent,newent->efmt);
                db_unpickaccno(newent->cp,&newent->acclst);
                newent->cp = NULL;
                break;
              case DBLC_ident:
              case DBLC_naident:
              case DBLC_datestamp:
              case DBLC_description:
              case DBLC_keyword:
              case DBLC_source:
              case DBLC_taxonomy:
              case DBLC_reference:
              case DBLC_author:
              case DBLC_title:
              case DBLC_journal:
              case DBLC_remark:
              case DBLC_basecount:
                newent->cp = NULL;
                break;
              case DBLC_sequence:
                if (!ndxrun)
                  db_sqlin2sq(newent);
                else
                  newent->cp = NULL;
                break;
              case DBLC_feature:
                if (ndxrun)
                  newent->cp = NULL;
                else
                  (void) db_scnnafeat(newent);
                break;
              case DBLC_terminator:
                newent->cp = NULL;
                break;
              case DBLC_unk:
              default:
                newent->cp = NULL;
                break;
              }
            break;
          }
        break;
      }
    }
  }
if (!gotid)
  db_dispos_ent(&newent);
return(newent);
}

DB_ENTSTRCT *db_parsecore(FILE *df,
                          WRD_LUSTRCT *lclass,
                          DBLU_DBFMT fmt,
                          DB_FTYPELT *ufwrds,
                          int ndxrun)
/* scan df from current position to next sq terminator.  wlus is used for
key words.  Return ptr to entry structure, NULL if no ident found.
scan string intrst and note what is to be kept from feature table.
ndxrun flags if it is only an indexing run or whether it is necessary to
create storage for and read the sequence data  */
{
return(db_parsemain(df,lclass,fmt,ufwrds,ndxrun,NULL));
}

DB_ENTSTRCT *db_parseflent(FILE *df,
                           WRD_LUSTRCT *lclass,
                           DBLU_DBFMT fmt,
                           char *ofintrst,
                           char *fqwants,
                           int ndxrun)
/* scan df from current position to next sq terminator.  wlus is used for
key words.  Return ptr to entry structure, NULL if no ident found.
scan string ofintrst and note what is to be kept from feature table.
ndxrun flags if it is only an indexing run or whether it is necessary to
create storage for and read the sequence data */
{
DB_ENTSTRCT *newent;      /* new entry record */
DB_FTYPELT *ufwrds;       /* list of user feature words */
WRD_LUSTRCT *ftkw;        /* feature table keywords */
WRD_LUSTRCT *fqkw;        /* feat qual keywords */
DB_FTYPELT *fep;          /* wanted feat element pointer */

ftkw = db_getkwstrct(fmt);
ufwrds = NULL;
(void) db_scnfeatwrds(ofintrst,ftkw,&ufwrds,0);
fep = ufwrds;
fqkw = db_getftqualstrct();
while (fep != NULL)  /* put all fqwants for each ufwrd */
  {
  (void) db_scnqulwrds(fqwants,fqkw,&fep->wquals,0);
  fep = fep->nxtfelt;
  }
newent = db_parsecore(df,lclass,fmt,ufwrds,ndxrun);
wlu_clralllustrct(fqkw);
wlu_clralllustrct(ftkw);
db_killfeltlst(&ufwrds);
return(newent);
}

DB_ENTSTRCT *db_parsfl4feats(FILE *df,
                             DBLU_DBFMT fmt,
                             DB_FTYPELT *ftwrds)
/* scan df from current position to next sq terminator.  Note features
in ftwrds and qualifiers from fqwrds.
Return ptr to entry structure, NULL if no ident found.*/
{
DB_ENTSTRCT *newent;      /* new entry record */
WRD_LUSTRCT *lclass;      /* line class tables */

lclass = db_getclasstrct(fmt);
newent = db_parsecore(df,lclass,fmt,ftwrds,0);
wlu_clralllustrct(lclass);
return(newent);
}

DBLU_LINECLASS db_netlclass(DBLU_LINECLASS curlt,
                            DBLU_LINECLASS prvlt,
                            DBLU_DBFMT fmt)
/* depending on fmt, return the current database file category for these
settings */
{
switch (fmt)
  {
  case DBFMT_embl:
  case DBFMT_embl_old:
  case DBFMT_swiss:
  case DBFMT_sqmonk:
    return(curlt);
    break;
  case DBFMT_genbank:
  default:
    if (curlt != DBLC_unk)
      return(curlt);
    else
      if (prvlt == DBLC_sequence)
        return(DBLC_unk);
      else
        return(prvlt);
    break;
  }
}

int db_lstparse(FILE *sf,
                WRD_LUSTRCT *lclass,
                DBLU_DBFMT fmt,
                DB_FTYPELT *unlmted,
                DB_FTYPELT **lmted,
                FILE *df)
/* scan sf from current position to next sq terminator.  wlus is used for
key words.  list lines to df according to contents of unlmted & lmtd lists.
return no of lines printed */
{
int scanon;
DBLU_LINECLASS prvclss;
DBLU_LINECLASS curclss;
char *slbuf;
char *cp;
int slcnt;
DB_FTYPELT *flt;
int pcnt;

slbuf = (char *) getmemory((DBLU_MAXSRCLEN+1),"Source buffer");
cp = NULL;
slcnt = 1;
scanon = 1;
pcnt = 0;
while (scanon)
  {
  if (cp == NULL)
    cp = db_fgets(slbuf,DBLU_MAXSRCLEN,sf,&slcnt);
  if (cp == NULL)      /* EOF */
    {    
    scanon = 0;
    if (prvclss != DBLC_sequence)
      db_logerr("Unexpected EOF",slcnt,NULL);
    }
  else
    {
    curclss = db_classfylin(slbuf,fmt,lclass);
    if (db_ptr4felt(unlmted,(int) db_netlclass(curclss,prvclss,fmt)) != NULL)
      {
      fputs(slbuf,df);
      fputc('\n',df);
      pcnt++;
      }
    else
      if ((flt = db_ptr4felt(*lmted,(int) db_netlclass(curclss,prvclss,fmt)))
             != NULL)  /* limited list, delete entry & print */
        {
        fputs(slbuf,df);
        pcnt++;
        (void) db_delfelt(lmted,flt);
        if (db_ptr4felt(*lmted,(int) db_netlclass(curclss,prvclss,fmt))
              == NULL)
          fputs("...",df);
        fputc('\n',df);
        }
    switch (curclss)
      {
      case DBLC_terminator:
        cp = NULL;
        scanon = 0;
        break;
      case DBLC_ident:
      case DBLC_access:
      case DBLC_naident:
      case DBLC_description:
      case DBLC_keyword:
      case DBLC_source:
      case DBLC_taxonomy:
      case DBLC_reference:
      case DBLC_author:
      case DBLC_title:
      case DBLC_journal:
      case DBLC_remark:
      case DBLC_basecount:
      case DBLC_ignore:
      case DBLC_seqversn:
      case DBLC_refnumbr:
      case DBLC_refpositn:
      case DBLC_xref:
      case DBLC_dbxref:
      case DBLC_datestamp:
      case DBLC_sequence:
      case DBLC_feature:
      case DBLC_unk:
      default:
        prvclss = db_netlclass(curclss,prvclss,fmt);
        cp = NULL;
        break;
      }
    }
  }
return(pcnt);
}

char *db_strdupn(char *str,
                 size_t sln)
/* malloc() sln+1 bytes, copy sln chars of str, null complete & return */
{
char *nbuf;

nbuf = (char *) getmemory((sln+1),"strdup buffer");
bcopy(str,nbuf,sln);
*(nbuf+sln) = '\0';
return(nbuf);
}
                      
void dbt_putidstr2fl(FILE *ofl,
                     int fwidth,
                     int *cc,
                     char *idst,
                     int idln)
/* attempt to put idstr contents to ofl. Keep track of filewidth and cc */
{
if ((fwidth > 0) && ((idln + *cc) >= fwidth))
  {
  fputc('\n',ofl);
  *cc = 0;
  }
if (*cc > 0)
  *cc += bas_coutputchr(ofl,' ');
*cc += bas_coutputstr(ofl,idst,'\0',idln,bas_coutputchr);
}

void dbt_putstr2fl(FILE *ofl,
                   int fwidth,
                   int *cc,
                   char *idst)
/* attempt to put idstr contents to ofl. Keep track of filewidth and cc */
{
dbt_putidstr2fl(ofl,fwidth,cc,idst,strlen(idst));
}

DBLU_DBFMT db_chkflfmt(FILE *fl)
  /* check for the database type of the file open on fl by checking each line
for Genbank or Swiss/EMBL id.  Determine latter by checking for NAtype on
id line.  Return DBFMT_unknown for indeterminate file. Rewind file. */
{
WRD_LUSTRCT *gb_lclss;
WRD_LUSTRCT *mbl_lclss;
DBLU_DBFMT ffmt;
char flin[DBLU_MAXSRCLEN+1];
int lcnt;
WRD_LUSTRCT *nwlu;

ffmt = DBFMT_unknown;
lcnt = 0;
gb_lclss = db_getclasstrct(DBFMT_genbank);
mbl_lclss = db_getclasstrct(DBFMT_embl);
nwlu = db_getnatypestrct();
while ((ffmt == DBFMT_unknown) && 
         (db_fgets(&flin[0],DBLU_MAXSRCLEN,fl,&lcnt) != NULL))
  if (db_classfylin(&flin[0],DBFMT_genbank,gb_lclss) == DBLC_ident)
      ffmt = DBFMT_genbank;
  else
    if (db_classfylin(&flin[0],DBFMT_embl,mbl_lclss) == DBLC_ident)
      if (db_scanidlin(&flin[0],DBFMT_embl,nwlu,NULL,NULL,NULL,NULL,NULL,NULL,
                         NULL,NULL) == DBNA_prot)
        ffmt = DBFMT_swiss;
      else
        ffmt = DBFMT_embl;
wlu_clralllustrct(gb_lclss);
wlu_clralllustrct(mbl_lclss);
wlu_clralllustrct(nwlu);
rewind(fl);
return(ffmt);
}

char *db_dbfmt2str(DBLU_DBFMT fmt)
  /* return ptr to string for fmt */
{
switch (fmt)
  {
  case DBFMT_embl:
    return("EMBL");
    break;
  case DBFMT_embl_old:
    return("EMBLold");
    break;
  case DBFMT_swiss:
    return("SwissProt");
    break;
  case DBFMT_genbank:
    return("Genbank");
    break;
  case DBFMT_sqmonk:
    return("SeqMonk_features");
    break;
  case DBFMT_gff3:
    return("GFF3");
    break;
  case DBFMT_gtf:
    return("GTF");
    break;
  case DBFMT_unknown:
  default:
    return("Unknown format");
    break;
  }
}

SQFL_SQTOS db_natype2sqtos(DB_NATYPE nat)
  /* remapping function: return a sqfl_fns.h sequence type respresentation
  of dbpars.h nat */ 
{
switch (nat)
  {
  case DBNA_prot:
    return(SQFL_peptide);
    break;
  case DBNA_dna:
    return(SQFL_dna);
    break;
  case DBNA_rna:
  case DBNA_trna:
  case DBNA_rrna:
  case DBNA_mrna:
  case DBNA_urna:
    return(SQFL_rna);
    break;
  case DBNA_unk:
  default:
    return(SQFL_tosunknown);
    break;
  }
}

DBP_MULTI_ELEMNT *dbp_appnd_mult_elemnt(DBP_MULTI_ELEMNT **melst,
                                        char *meid,
                                        DB_ENTSTRCT *entptr)
/* append a new element to *melst, return address of new element if useful.
meid is strduped, entptr is likely to be NULL at this stage */
{
DBP_MULTI_ELEMNT *prev, *endp;

if (melst != NULL)
  {                     /* chain to end of list */
  prev = endp = *melst;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nxt_dme;
    }
  endp = (DBP_MULTI_ELEMNT *) getmemory(sizeof(DBP_MULTI_ELEMNT),"multiscan element");
  endp->nxt_dme = NULL;
  endp->me_id = bas_strdup(meid);
  endp->me_entstrptr = entptr;
  endp->flistchunkarray = NULL;
  endp->flistchunkcount = 0;
  if (*melst == NULL)
    {
    *melst = endp;
    endp->prv_dme = NULL;
    }
  else
    {
    prev->nxt_dme = endp;
    endp->prv_dme = prev;
    }
  }
return(endp);
}

DBP_MULTI_ELEMNT *dbp_melemnt4id(DBP_MULTI_ELEMNT *melst,
                                 char *queryid,
                                 int (*cmpfn)(const char *x1,
                                              const char *x2))
/* use cmpfn() to locate the first element matching
queryid, return its address, else NULL */
{
DBP_MULTI_ELEMNT *lp;

lp = melst;
while (lp != NULL)
  if ((*cmpfn)(queryid,lp->me_id) == 0)
    return(lp);
  else
    lp = lp->nxt_dme;
return(NULL);
}

void dbp_delmultelt(DBP_MULTI_ELEMNT *dep,
                    DBP_MULTI_ELEMNT **lstrt,
                    DBP_MULTI_ELEMNT **lend)
/* dep is an element of *lstrt: delete it from the list 
*lstrt..*lend.  lend can be NULL.  */
{
DBP_MULTI_ELEMNT *pt;

if (dep != NULL)
  {
  if ((pt = dep->prv_dme) == NULL)
    *lstrt = dep->nxt_dme;
  else
    pt->nxt_dme = dep->prv_dme;
  if ((pt = dep->nxt_dme) != NULL)
    pt->prv_dme = dep->prv_dme;
  else
    if (lend != NULL)
      *lend = dep->prv_dme;
  db_dispos_ent(&dep->me_entstrptr);
  memfree(dep->me_id);   /* was strdup()-ed */
  memfree(dep);
  }
}

void dbp_killmultieltlist(DBP_MULTI_ELEMNT **lstrt,
                          DBP_MULTI_ELEMNT **lend)
/* iteratively remove all of *lstrt..*lend.  lend
can be NULL */
{
while (*lstrt != NULL)
  dbp_delmultelt(*lstrt,lstrt,lend);
}

char *dbp_strngafterlead(char *sstrng,
                         char *leadstr,
                         int (*strncmpfn)(const char *s1,
                                          const char *s2,
                                          size_t l1))
/* use (*strcmpfn)() to locate leadstr at the head of
sstrng.  Return the address of the following portion of sstrng,
NULL if leadstr not found */
{
size_t lslen;

if (sstrng != NULL)
  {
  lslen = strlen(leadstr);
  if ((lslen > 0) && ((*strncmpfn)(leadstr,sstrng,lslen) == 0))
    return(sstrng+lslen);
  else
    return(NULL);
  }
else
  return(NULL);
}

char *dbp_cpystrmodurlesc(char *srcstr,
                          char *dststr)
/* scan srcstr, looking for %xx where (xx is hexadecimal)
and generate a copy of the string in dststr with
escaped URL chars substituted.  return address of
dststr for convenience.  Invalid xx are simply 
copied verbatim. */
{
char *sp;
char *dp;
char hcod[3];
char nc;

sp = srcstr;
dp = dststr;
while (*sp != '\0')
  {
  if (*sp == '%')
    {
    strncpy(&hcod[0],sp+1,2);
    hcod[2] = '\0';
    nc = (char) strtol(&hcod[0],NULL,16);
    *dp = nc;
    sp += 3;
    }
  else
    {
    *dp = *sp;
    sp++;
    }
  dp++;
  }
*dp = '\0';
return(dststr);
}

char *dbp_cpynreplurlescchars(char *ustr)
  /* malloc() ustr and return a version with
any URL-escaped chars replaced.  Strings
created by this function must be free()ed */
{
char *ucpy;

ucpy = (char *) getmemory(strlen(ustr)+1,"urlreplstring");
return(dbp_cpystrmodurlesc(ustr,ucpy));
}

int dbp_fndnappgff3info(char *toks[],
                        int tmax,
                        char *leadstr,
                        DBLU_FTQUAL ftq,
                        DB_FEATSTRCT *ftp)
/* scan thru toks looking for a string headed with leadstr.
take the remainder, replacing any URL-escaped chars and append
with qualifier-type ftq to infolist of ftp.  return 1 if
succeeded */
{
char *ustr;
int tpt;

tpt = 0;
while (tpt < tmax)
  if ((ustr = dbp_strngafterlead(toks[tpt],leadstr,strncmp)) != NULL)
    {
    ftp->infoend = db_appnielt(&ftp->infolist,ftq,
                                 dbp_cpynreplurlescchars(ustr));
    return(1);
    }
  else
   tpt++;
return(0);
}

int dbp_parse_gff3_ln(FILE *srcfl, 
                      DBP_MULTI_ELEMNT **curseq,
                      char *ln,
                      WRD_LUSTRCT *ftkw,
                      DB_FTYPELT *ftwrds,
                      DBP_MULTI_ELEMNT **elmntlst)
/* ln contains text: check for strlen() > 0;  if so, 
check leading '#' and return OK but do nothing if found.
otherwise break into ht separated tokens and check
if we've seen this.
Assume line is '\0' terminated.
if so, insert relevant information to that elemnt,
else append a new element and init it.
return 1 if something was parsed */
{
char **lp;
char *tokens[GFF3TOKENCNT];
char *strcpy;
char *oricpy;
int tcnt;
DBP_MULTI_ELEMNT *eltp;
DBLU_FTKW thisfeat;
char *thisid;
char *t9cpy;
char *t9cpori;
char *t9toks[GFF3TOKENCNT];
int t9cnt;
int tpt;
DB_FEATSTRCT *ftptr;
int fstrt;
int fstop;
int appcnt;
int slen;
char *gennotecpy;

if (strlen(ln) > 0)
  {
  switch (*ln)
    {
    case '#':       /* comment, ignore */
      break;
    case '>':       /* sequence header - find which id */
      *curseq = dbp_melemnt4id(*elmntlst,ln+1,strcmp);
      if ((*curseq != NULL) &&
            ((ftptr = db_feat4pos((*curseq)->me_entstrptr,FTKW_source_ftkw,1)) != NULL))
        {
        (*curseq)->me_entstrptr->sqlen = db_3pextnt4feat(ftptr);
        (*curseq)->me_entstrptr->sfpt = (*curseq)->me_entstrptr->seq
           = (char *) getmemory((*curseq)->me_entstrptr->sqlen+1,"Seq buff");
        *(*curseq)->me_entstrptr->seq = '\0';
        }
      else
        *curseq = NULL;
      break;
    default:        /* everything else */
      if (*curseq != NULL)
        {
        slen = (*curseq)->me_entstrptr->sfpt - (*curseq)->me_entstrptr->seq;
        appcnt = imin(strlen(ln),(*curseq)->me_entstrptr->sqlen-slen);
        (void) strncat((*curseq)->me_entstrptr->sfpt,ln,appcnt);
        (*curseq)->me_entstrptr->sfpt += appcnt;
        }
      else
{  /* work on copy, since contents get altered */
        {
        strcpy = oricpy = bas_strdup(ln);
        tcnt = 0;
        for (lp = tokens; (*lp = strsep(&strcpy,"\t")) != NULL;)
          {
          tcnt++;
          if (**lp != '\0')
            if (++lp >= &tokens[GFF3TOKENCNT])
              break;
          }
        if ((eltp = dbp_melemnt4id(*elmntlst,tokens[0],strcmp)) == NULL)
          {
          eltp = dbp_appnd_mult_elemnt(elmntlst,tokens[0],NULL);
          eltp->genids = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"genids lu");
          wlu_initlustrctptr(eltp->genids,WLU_CASEIND,NULL);
          eltp->me_entstrptr = (DB_ENTSTRCT *) getmemory(sizeof(DB_ENTSTRCT),"Entry structure");
          db_initentstrct(eltp->me_entstrptr,DBFMT_gff3,ftkw);
          eltp->me_entstrptr->freefeatids = 1;
          }
        thisfeat = wlu_chkwrd(eltp->me_entstrptr->fkwlu,tokens[2]);
    /* do we want this sort of feature?? */
        eltp->me_entstrptr->nseen++;
        if (db_ptr4felt(ftwrds,thisfeat) != NULL)
          {
    /* have we already seen this feature?? */
          t9cpy = t9cpori = bas_strdup(tokens[8]);
          t9cnt = 0;
          for (lp = t9toks;(*lp = strsep(&t9cpy,";")) != NULL;)
            {
            t9cnt++;
            if (**lp != '\0')
              if (++lp >= &t9toks[GFF3TOKENCNT])
                break;
            }
          tpt = 0;
          thisid = NULL;
          while (tpt < t9cnt)
            if ((thisid = dbp_strngafterlead(t9toks[tpt],"Name=",strncmp)) != NULL)
              tpt = t9cnt;
            else
              tpt++;
          if (thisid != NULL)
            {
            if ((ftptr = db_nxtfeat4id(eltp->me_entstrptr->featlst,thisid,strcmp)) == NULL)
              {
              ftptr = eltp->me_entstrptr->lastfeat
                = db_appnfeat(&eltp->me_entstrptr->featlst,thisfeat,bas_strdup(thisid));
/* look for any "Parent=" token */
              ftptr->featno = eltp->me_entstrptr->nseen;
              ftptr->savdno = ++eltp->me_entstrptr->nfeats;
              tpt = 0;
              thisid = NULL;
              while (tpt < t9cnt)
                if ((thisid = dbp_strngafterlead(t9toks[tpt],"Parent=",strncmp)) != NULL)
                  {
                  ftptr->parentname = bas_strdup(thisid);
                  tpt = t9cnt;
                  }
                else
                  tpt++;
              if (strcmp("+",tokens[6]) == 0)
                ftptr->strctsens = DB_sens_norm;
              else
                ftptr->strctsens = DB_sens_comp;
              }
    /* need to append details for this ... */
            if (((fstrt = (int) strtol(tokens[3],NULL,10)) > 0) &&
                 ((fstop = (int) strtol(tokens[4],NULL,10)) > 0))  /* decoded OK */
              {
              eltp->me_entstrptr->curseg = ftptr->lstseg =
                db_appnseg(&ftptr->fstseg,fstrt,fstop,ftptr->strctsens,NULL);
              }
            if (thisfeat == FTKW_FT_gene)
              {
              (void) dbp_fndnappgff3info(&t9toks[0],t9cnt,"Note=",FTQU_note,ftptr);
              (void) dbp_fndnappgff3info(&t9toks[0],t9cnt,"dbxref=",FTQU_db_xref,ftptr);
              (void) dbp_fndnappgff3info(&t9toks[0],t9cnt,"Ontology_term=",FTQU_GO_info,ftptr);
              (void) dbp_fndnappgff3info(&t9toks[0],t9cnt,"gene=",FTQU_gene,ftptr);
              }
            }
          memfree(t9cpori);
          }
        memfree(oricpy);
        }
      return(1);
      }
    }
  return(0);
  }
else
  return(0);

}

int dbp_findparentfeats(DB_ENTSTRCT *esp,
                        int (*strcmpfn)(const char *s1,
                                        const char *s2),
                        int lnkgene2chr)
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
{
DB_FEATSTRCT *ftpt;
DB_FEATSTRCT *fp2;
DB_FEATSTRCT *srcftp;
int acnt;

acnt = 0;
if (esp != NULL)
  {
  srcftp = db_feat4pos(esp,FTKW_source_ftkw,1);
  ftpt = esp->featlst;
  while (ftpt != NULL)
    {
    if (lnkgene2chr && (ftpt->featur == FTKW_FT_gene)) /* check if current feature has source parent */
      {
      ftpt->parentfeat = srcftp;
      if (ftpt->parentname != NULL)
        memfree(ftpt->parentname);
      ftpt->parentname = NULL;
      acnt++;
      }
    else
      {
      if (ftpt->parentname != NULL)
        {
        fp2 = esp->featlst;       /* must traverse entire list */
        while (fp2 != NULL)
          if ((fp2 != ftpt) && (fp2->idp != NULL))
     /* don't waste time comparing with self */
            if (((*strcmpfn)(ftpt->parentname,fp2->idp) == 0) ||
                 (db_matstrelt(fp2->enamelist,ftpt->parentname,strcmpfn) != NULL))
              {
              ftpt->parentfeat = fp2;
              memfree(ftpt->parentname);
              ftpt->parentname = NULL;
              acnt++;
              if ((ftpt->featur == FTKW_CDS) && (fp2->featur == FTKW_mRNA))
                ftpt->mRNA = fp2;
              fp2 = NULL;     /* should only have one parent, so stop now */
              }
            else
              fp2 = fp2->nextst;
          else
            fp2 = fp2->nextst;
        }
      }    
    ftpt = ftpt->nextst;
    }
  }
return(acnt);
}

char *dbp_fndndupquotdgtfinfo(char *toks[],
                              int tmax,
                              char *leadstr)
/* scan thru toks looking for a string headed with leadstr, 
allowing for a possible leading blank.  strdup() the quoted
string following a match and return the address of that.
Else return NULL.  The output from this function must
be free()ed after use, since it is malloc()ed */
{
char *ustr;
int tpt;
char *up;
char *retstr;

tpt = 0;
while (tpt < tmax)
  if ((toks[tpt] != NULL) && (((ustr = dbp_strngafterlead(toks[tpt],leadstr,strncmp)) != NULL) ||
      ((*toks[tpt] == ' ') && ((ustr = dbp_strngafterlead(toks[tpt]+1,leadstr,strncmp)) != NULL))))
    {
/* ustr should now have trailing "String", so scan till after first " and work from there */
    up = ustr;
    while ((*up != '\0') && (*up != '"'))
      up++;
    if (*up == '"')
      {
      retstr = bas_strdup(up+1);
      up = retstr;
      while (index("\"",*up) == NULL)
        up++;
      *up = '\0';
      return(retstr);
      }
    }
  else
    tpt++;
return(NULL);
}

char *dequotebiotypestr(char *str)
  /* return ptr to a string copy with
double quotes removed. */
{
int clen;

clen = strlen(str) - 2;
(void) strncpy(&btypecopy[0],(str+1),clen);
btypecopy[clen] = '\0';
return(&btypecopy[0]);
}

int db_gtf_attr_match(char *tokens[],
                      int ntokns,
                      WRD_LUSTRCT *entftkwds,
                      DB_FTYPELT *ftwrds,
                      char *attr2fnd)
/* scan gtf line attributes from start,
looking for attr2fnd.  tokens ends up with
attribute types in index 8,10,12..., with attribute
values in 9,11,13.... Return true if its value
is feat2fnd.  tokens contain '"', so need to
remove them, probably copy.
Return token count of matching string value, 0 for
not found */
{
int tpt;

tpt = 8;
while (tpt < ntokns)
  {
  if (strcmp(tokens[tpt],attr2fnd) == 0)
    {
/* attribute values contain '"', so must lose them */
    if (wlu_chkwrd(entftkwds,dequotebiotypestr(tokens[tpt+1])) != FTKW_unknown)
      return(tpt+1);
    }
  tpt++;
  }
return(0);
}

int dbp_parse_gtf_ln(FILE *srcfl, 
                     DBP_MULTI_ELEMNT **curseq,
                     char *ln,
                     WRD_LUSTRCT *ftkw,
                     DB_FTYPELT *ftwrds,
                     WRD_LUSTRCT *wantedgtypes,
                     DBP_MULTI_ELEMNT **elmntlst,
                     DB_STR_ELT *attr)
/* ln contains text: check for strlen() > 0;  if so, 
check leading '#' and return OK but do nothing if found.
otherwise break into ht separated tokens and check
if we've seen this.
Assume line is '\0' terminated.
if so, insert relevant information to that elemnt,
else append a new element and init it.
if attr is non-NULL, then use it to select the
attribute values for the last field.
return 1 if something was parsed.
There is significant variation in how
the attribute field is used.  This version
allows for '\t', ' ', '=' & ';' to separate
the attribute pairs.  Attributes can then
be examined by skipping through list in pairs.
If wantedgtypes is non-NULL, then check that the
gene_type attribute is valid. */
{
char **lp;
char *tokens[GTFMAXTOKCNT];
char *strcpy;
char *oricpy;
int tcnt;
DBP_MULTI_ELEMNT *eltp;
DBLU_FTKW thisfeat;
char *thisid;
int tpt;
DB_FEATSTRCT *ftptr;
int fstrt;
int fstop;
int appcnt;
int slen;
char *gennotecpy;
DB_FEATSTRCT *firstidfeatp;
DB_STR_ELT *strp;
int p5val;
DBLU_FTQUAL ftqual;
int gtf_attr_mat_no;

/* check line for leading '#' */
if (*ln != '#')
  {
/* work on copy, since contents get altered */
  strcpy = oricpy = bas_strdup(ln);
  tcnt = 0;
  for (lp = tokens; (*lp = strsep(&strcpy," \t=;")) != NULL;)
    if (**lp != '\0')
      {
      tcnt++;
      if (++lp >= &tokens[GTFMAXTOKCNT])
        break;
      }
  if ((eltp = dbp_melemnt4id(*elmntlst,tokens[0],strcmp)) == NULL)
    {
    eltp = dbp_appnd_mult_elemnt(elmntlst,tokens[0],NULL);
    eltp->genids = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"genids lu");
    wlu_initlustrctptr(eltp->genids,WLU_CASEIND,NULL);
    eltp->me_entstrptr = (DB_ENTSTRCT *) getmemory(sizeof(DB_ENTSTRCT),"Entry structure");
    db_initentstrct(eltp->me_entstrptr,DBFMT_gtf,ftkw);
    eltp->me_entstrptr->freefeatids = 1;
    eltp->me_entstrptr->ename = bas_strdup(tokens[0]);
    }
  thisfeat = wlu_chkwrd(eltp->me_entstrptr->fkwlu,tokens[2]);
  /* do we want this sort of feature?? */
  eltp->me_entstrptr->nseen++;
  if (((gtf_attr_mat_no = db_gtf_attr_match(tokens,tcnt,wantedgtypes,ftwrds,"gene_type")) > 0) ||
       ((gtf_attr_mat_no = db_gtf_attr_match(tokens,tcnt,wantedgtypes,ftwrds,"gene_biotype")) > 0))  /* T2T gtf uses this form of attribute */
    ftqual = FTQU_biotype;
  else
    ftqual = FTQU_unknown;
  if ((db_ptr4felt(ftwrds,thisfeat) != NULL) && ((wantedgtypes == NULL) || (ftqual != FTQU_unknown)))
    {
/*  have we already seen this feature?? */
    tpt = 8;
    thisid = NULL;
    while ((tpt < tcnt) && (tpt < GTFMAXTOKCNT - 2))
      if (db_matstrelt(attr,tokens[tpt],strcmp) != NULL)
        {
        thisid = tokens[tpt+1];
        tpt = tcnt;
        }
      else
        tpt++;
    if (thisid != NULL)
      {
      if (((firstidfeatp = (DB_FEATSTRCT *) wlu_chkwrdptr(eltp->genids,thisid)) == NULL) ||
           ((ftptr = db_nxtfeat4typenid(firstidfeatp,thisfeat,
                                          thisid,strcmp)) == NULL))
  /* haven't seen this id and feat type before */
        {
        ftptr = eltp->me_entstrptr->lastfeat
          = db_appnfeat(&eltp->me_entstrptr->featlst,thisfeat,thisid);
        ftptr->featno = eltp->me_entstrptr->nseen;
        ftptr->savdno = ++eltp->me_entstrptr->nfeats;
        switch (*tokens[6])
          {
          case '-':
            ftptr->strctsens = DB_sens_comp;
            break;
          case '+':
          default:
            ftptr->strctsens = DB_sens_norm;
          break;
          }
        ftptr->lstseg = db_appnseg(&ftptr->fstseg,(int) strtol(tokens[3],NULL,10),
                                     (int) strtol(tokens[4],NULL,10),
                                     ftptr->strctsens,NULL);
/* fputs("appending: ",stdout);
db_sayafeatcore(stdout,eltp->me_entstrptr,ftptr,0); */
        if (firstidfeatp == NULL)  /* haven't even seen the id */
          wlu_addwrdptr(eltp->genids,thisid,ftptr,NULL);
        ftptr->featno = eltp->me_entstrptr->nseen;
        ftptr->savdno = ++eltp->me_entstrptr->nfeats;
        if (strcmp("+",tokens[6]) == 0)
          ftptr->strctsens = DB_sens_norm;
        else
          ftptr->strctsens = DB_sens_comp;
        ftptr->idp = bas_strdup(thisid);
        ftptr->infoend = db_appnielt(&ftptr->infolist,FTQU_biotype,bas_strdup(dequotebiotypestr(tokens[gtf_attr_mat_no])));
        }
      else /* does this new feature qualify an existing one? (e.g. exons in gene) */
        {
        if ((firstidfeatp != NULL) && (firstidfeatp->featur != FTKW_exon) &&
             (thisfeat == FTKW_exon))
          {
          if (!firstidfeatp->seg_append_ok)
            {
            db_killsegs4featstrct(firstidfeatp);
            firstidfeatp->seg_append_ok = 1;
            }
          p5val = (int) strtol(tokens[3],NULL,10);
          if ((firstidfeatp->fstseg != NULL) &&
                 (p5val < db_5pextnt4feat(firstidfeatp)))
            (void) db_prepndseg(&firstidfeatp->fstseg,
                                  p5val,(int) strtol(tokens[4],NULL,10),
                                  firstidfeatp->strctsens,NULL);
          else
            firstidfeatp->lstseg = db_appnseg(&firstidfeatp->fstseg,
                                                p5val,(int) strtol(tokens[4],NULL,10),
                                                firstidfeatp->strctsens,NULL);
          }
        }
      }
    }
  memfree(oricpy);
  return(1);
  }
else
  return(0);
}

int dbp_gtfgenes4entstrct(DB_ENTSTRCT *esptr,
                          WRD_LUSTRCT *idlut)
  /* scan ent structure esptr, looking
for exon groups which can be regarded as an mRNA.
reset the feature type for these and create a
corresponding gene entry.  Create other relevant
links for related CDS feature.
return number of mRNAs converted */
{
DB_FEATSTRCT *fp;
int rcnt;
DB_FEATSTRCT *wkingfp;
DB_FEATSTRCT *idfeatp;
WRD_LU_REC *idwrdrec;

rcnt = 0;
if (esptr != NULL)
  {
  fp = esptr->featlst;
  while (fp != NULL)
    {
    if (fp->featur == FTKW_exon)
      {
      fp->featur = FTKW_mRNA;
      if (((idwrdrec = wlu_wrd2lurec(idlut,fp->idp)) != NULL) &&
           ((wkingfp = db_nxtfeat4typenid((DB_FEATSTRCT *) idwrdrec->retval,
                                            FTKW_CDS,fp->idp,strcmp)) != NULL))
        {
        wkingfp->mRNA = fp;
        wkingfp->parentfeat = fp;
        wkingfp->parentname = bas_strdup(fp->idp);
        idwrdrec->retval = wkingfp;
        }
/* now create a gene entry */
      wkingfp = dbp_prependfeatstrct(&esptr->featlst,fp,FTKW_FT_gene,bas_strdup(fp->idp));
/*      wkingfp = esptr->lastfeat = db_appnfeat(&esptr->featlst,FTKW_FT_gene,bas_strdup(fp->idp)); */
      fp->parentfeat = wkingfp;
      fp->parentname = bas_strdup(fp->idp);
      wkingfp->strctsens = fp->strctsens;
      if (fp->strctsens != DB_sens_norm)
        wkingfp->lstseg = db_appnseg(&wkingfp->fstseg,db_3ppos4feat(fp),
                                       db_5ppos4feat(fp),fp->strctsens,NULL);
      else
        wkingfp->lstseg = db_appnseg(&wkingfp->fstseg,db_5ppos4feat(fp),
                                       db_3ppos4feat(fp),fp->strctsens,NULL);
      rcnt++;
      }
    fp = fp->nextst;
    }
  }
return(rcnt);
}

int dbp_bldgtfgenes(DBP_MULTI_ELEMNT *melstp)
  /* progressively scan melstp ent structures, looking
for exon groups which can be regarded as an mRNA.
reset the feature type for these and create a
corresponding gene entry.  Create other relevant
links for related CDS feature.
return number of mRNAs converted */
{
DBP_MULTI_ELEMNT *mep;
int rcnt;

mep = melstp;
rcnt = 0;
while (mep != NULL)
  {
  rcnt += dbp_gtfgenes4entstrct(mep->me_entstrptr,mep->genids);
  mep = mep->nxt_dme;
  }
return(rcnt);
}

