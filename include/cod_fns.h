/* cod_fns.h: headers for codon-related functions in c */

#define CF_NC_RARE 1.0e-6
#define CF_NON_RESULT -999.99
#define CF_NC_CUTOFF 100
#define CF_NC_MAX 61.0
#define CF_RSCUATOM 24
#define CF_MAXSRCLEN 80

typedef enum CF_rscuptype     /* parsing RSCU files */
  {
  CF_rscu_res = 0,
  CF_rscu_cdn,
  CF_rscu_real
  }
CF_RSCUPTYPE;

/* routine headers */

int cf_codno4cvptr(COD_VECT_STRCT *cvp);
  /* iterate thru cvp list counting no of elements */

int cf_codno4res(COD_VECT_STRCT **cvlist,
                 char aares);
/* return the number of codons for this residue - by iterating thru linked
list */

int cf_totimtrx(TRANS_MATRX ccnts);
  /* return the sum of valid ccnts */

int cf_cnts4veclst(TRANS_MATRX cnts,
                   COD_VECT_STRCT *cvec);
/* return the total cnts entries for this cvec list */

int cf_totcnts4res(TRANS_MATRX cnts,
                   COD_VECT_STRCT **cvs,
                   char ares);
/* sum up the total counts in cnts for residue ares */

char *cf_gc3_ref();

double cf_cnts2gc3(TRANS_MATRX cnts,
                   COD_VECT_STRCT **cvec);
/* return the value of gc3 statistic on cnts.  gc3 = proportion of codons
with c/g in 3rd position, discounting singleton codons and terms.
Ikemura (1985), Mol Biol Evol, 2, 13-34. */

double cf_hzygot(COD_VECT_STRCT *cvlst,
                 TRANS_MATRX cnts);
/* return homozygosity value for cvlst.  Return CF_NON_RESULT if "rare" aa
according to Wright's criterion */

char *cf_nc_ref();

double cf_cnts2nc(TRANS_MATRX ccnts,
                  COD_VECT_STRCT **cvecs);
/* perform the calculation for "effective number" of codons used in a gene;
  for ccounts.  Wright (1990), Gene, 87, 23-29 */

int cf_chkrscutbl(FILE *lfl,
                  TRANS_MATRX rscu,
                  COD_VECT_STRCT **cvecs,
                  float prec);
/* check that values in rscu table conform to requirement that they should
all add to n for each aa where n is the no of codons.  Allow prec*n deviation.
log errors if lfl is non-NULL */

int cf_rscuexfile(FILE *fl,
                  TRANS_MATRX tm);
/* attempt to read a rscu matrix from fl - assume "fish_term" std format:
     ================================================
     F TTT 1.143 S TCT 1.467 Y TAT 0.500 C TGT 2.000
     F TTC 0.857 S TCC 2.133 Y TAC 1.500 C TGC 0.000
     L TTA 0.067 S TCA 0.000 * TAA 0.000 * TGA 0.000
     L TTG 0.400 S TCG 0.133 * TAG 0.000 W TGG 1.000
     ================================================
     L CTT 0.000 P CCT 2.483 H CAT 0.872 R CGT 0.000
...

& drop lines starting with '#'

Return 1 for successful read */

void cf_rscu2reladpt(TRANS_MATRX rscu,
                     COD_VECT_STRCT **cvecs,
                     TRANS_MATRX rad);
/* relative adaptiveness is ration of given rscu to the max for that residue.
*/

char *cf_cai_ref();

double cf_cnts2cai(TRANS_MATRX cnts,
                   COD_VECT_STRCT **cvec,
                   TRANS_MATRX rscu);
/* use the rscu table and cvec to calculate the Codon Adaptive Index
for cnts.  Sharp & Li, (1987), NAR, 15, 1281). */

void cf_veccnts2rscutbl(COD_VECT_STRCT **cvs,
                        TRANS_MATRX cnts,
                        TRANS_MATRX rscu);
/* use codon vector set cvs to calculate rscu table for cnts */

void cf_cnts2rscutbl(TRANS_MATRX tmat,
                     TRANS_MATRX cnts,
                     TRANS_MATRX rscu);
/* make a codon vector set from tmat, then use this to calculate RSCU table */
