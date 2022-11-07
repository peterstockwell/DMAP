/* diffmeth: generate bins based on MspI (or other) fragments or fixed width, looking for
CpGs.  Scan different "Cno position +/-" lines for various
individuals/treatments and compare positions & run statistics.
Or read SAM/BAM files from bismark or other  mapping runs */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <sys/param.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "cmaths.h"
#include "mrg_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"
#ifndef NO_ZLIB

#include <zlib.h>
#include <sys/types.h>
#include "bam_fns.h"

#endif
#include "diffmeth.h"
#include "fsm_ops.h"

/* global debug & out style variables, for simplicity of access */
RBC_DBG debuglevel;
RBC_SRC_STYLE srcstyle;
int glblreadlen;          /* global value for readlength */
char istring[20];

int err_msg(char *fmt,
            ...)
/* write user error message.  Return 0 for err return status */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
return(0);
}

void err_msg_die(char *fmt,
                 ...)
/* write user error message then exit with error status */
{
va_list args;

va_start(args,fmt);
(void) vfprintf(stderr,fmt,args);
va_end(args);
exit(1);
}

void say_usage(FILE *fl,
               char *pnam)
{
fprintf(fl,"%s v%.2f: make MspI (or other) fragment or fixed bins, note meth/unmeth positions for individuals\n",
          pnam,PROG_VERSION);
fputs("Options:\n",fl);
fputs("     -r <posfile> read <posfile> as set of chr posit strand/meth (multiples allowed)\n",fl);
#ifdef NO_ZLIB
fputs("     -R <samfile> read info from .sam file (multiples allowed)\n",fl);
#else
fputs("     -R <samorbamfile> read info from .sam/.bam file (multiples allowed)\n",fl);
#endif
fputs("     -s <posfile> as for -r, second group data\n",fl);
#ifdef NO_ZLIB
fputs("     -S <samfile> as for -R, second group data\n",fl);
#else
fputs("     -S <samorbamfile> as for -R, second group data\n",fl);
fputs("     -z/-Z switch between sam (-Z) & bam (-z) for -R/-S input - positional (def=-Z)\n",fl);
#endif
fputs("     -g <genomehead> dir and file string to locate genomic seq files by adding n.fa\n",fl);
fputs(" or  -G <chr_info_file> file of chromosome IDs and filenames\n",fl);
fputs(
"     -p <j,k> pairwise scan meth variance of RR frags j..k (Li, et al. (2010) PLOSBiology,11,e1000533)\n"
,fl);
fputs(
"     -P <j,k> pairwise scan force Fisher's exact frags j..k (Li, et al. (2010) PLOSBiology,11,e1000533)\n"
,fl);
fputs("     -x <j,k> do Chi Square if possible, else pairwise FE\n",fl);
fputs("     -X <j,k> as -x but force Chi Square for valid count fragments\n",fl);
fputs("     -q <j,k> as -x but choose lower Pr of Fisher's Exact, or Chi if valid\n",fl);
fputs("     -Q <j,k> as -q but show all paired FE Prs, or Chi if valid\n",fl);
fputs("     -a/-A/-B <j,k> ANOVA on %meth for groups; -A->show >meth group & sample counts; -B->more detail\n",
  fl);
fputs("     -e/-E <j,k> list each CpG for bins; -E=> only nonzero bins\n",fl);
fputs("     -D <j,k> derive mean & std deviation for valid bins\n",fl);
fputs("     -l/-L <j,k> list bins (-L=>only nonzero bins)\n",fl);
fputs("     -m treat complementary Cs of CpGs as separate (def=don't)\n",fl);
fputs("     -f <fold> require fold methylation difference (def=don't)\n",fl);
fputs("     -t <threshold> ignore CpGs with fewer than threshold counts (def=1)\n",fl);
fputs("     -T <hitthreshold> ignore bins with less than hitthreshold/CpG counts (def=0.0)\n",fl);
fputs("     -F <cntscrit> No. CpGs that must meet -t count criterion (def=all)\n",fl);
fputs("     -u <maxhitspercpg> ignore bins with more counts/CpG (def=0.0=ignore)\n",fl);
fputs("     -M <minmethprop> ignore bins with methylation below this (def=0.0=ignore)\n",fl);
fputs("     -U <prthreshold> ignore comparisons with Pr > prthreshold (def=don't=>-U 1.0)\n",fl);
fputs("     -k <maxchrno> allow upto maxchrno different chromosomes (def = ChrY, 24)\n",fl);
fputs("     -Y don't expect XY chromosomes, just numbers (def=do)\n",fl);
fprintf(fl,"     -b <sambuflen> max read length for .sam files (def=%d)\n",DEF_SAMBUFLEN);
fputs("     -d display counts for each individual (def=don't)\n",fl);
fputs("     -c <n> restrict effort to Chromosome <n> (def = all chromosomes)\n",fl);
fputs("     -C <n> restrict bins to those with <n> or more CpGs\n",fl);
fputs("     -I <n> minimum of <n> individuals for a fragment (def=2)\n",fl);
fputs("     -j join adjacent RRBS fragments (def=don't)\n",fl);
fprintf(fl,"     -K <clustersize> cluster bins to clustersize groups for quick lookup (def=%d)\n",
          DM_LIST_FREQ_DEF);
fputs("     -H <Chdr> use Chdr as prefix for each chromosome (def=none)\n",fl);
#ifdef NO_ZLIB
fputs("     -N for 3' SAM RRBS reads map leading CpG to prev fragment (def=don't)\n",fl);
#else
fputs("     -N for 3' SAM/BAM RRBS reads map leading CpG to prev fragment (def=don't)\n",fl);
#endif
fputs("     -W <binwidth> make fixed width bins (def=rest. enz by size)\n",fl);
fputs("     -y <binfilename> read bin info from binfilename (chr start stop) (def=rest. enz)\n",fl);
fputs("     -n <restrictionfile> - use sites in <restrictionfile>, def=MspI\n",fl);
/* fputs("      -B read source in 'chr pos strand C T' form\n",fl); */
fputs("     -J include non-CpG methylation\n",fl);
}

DM_FNAM_ELT *dm_appndfnam(DM_FNAM_ELT **fnmlst,
                          char *fnam,
                          DM_SRC_MODE srcmod,
                          int smplgroup,
                          char *grpid)
/* dup fnam and append an element for it to the end of fnmlst.
return address of new element in case useful */
{
DM_FNAM_ELT *prev, *end_ptr;

if (fnmlst != NULL)
  {
  prev = end_ptr = *fnmlst;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtfelt;
    }
  end_ptr = (DM_FNAM_ELT *) getmemory(sizeof(DM_FNAM_ELT),"fnam elt");
  end_ptr->nxtfelt = NULL;
  end_ptr->flnam = bas_strdup(fnam);
  end_ptr->srcmode = srcmod;
  end_ptr->sgroup = smplgroup;
  end_ptr->groupid = grpid;
  if (*fnmlst == NULL)
    *fnmlst = end_ptr;
  else
    prev->nxtfelt = end_ptr;
  return(end_ptr);
  }
else
  return(NULL);
}

int dm_cntfnamlst(DM_FNAM_ELT *fnamlst)
  /* recursively count fnamlst */
{
if (fnamlst == NULL)
  return(0);
else
  return(dm_cntfnamlst(fnamlst->nxtfelt) + 1);
}

void dm_delfnamlst(DM_FNAM_ELT **fnmlst)
/* iteratively elements from fnamlst */
{
DM_FNAM_ELT *fp;

while (*fnmlst != NULL)
  {
  fp = *fnmlst;
  *fnmlst = fp->nxtfelt;
  if (fp->flnam != NULL)
    memfree(fp->flnam);
  memfree(fp);
  }
}

DM_METH_CNTS *dm_appndcntelt(DM_METH_CNTS **clststrt,
                             int pcnt,
                             int mcnt,
                             int sgroup)
/* append a new meth/unmeth cout element to clststrt.
return address of the new element if useful */
{
DM_METH_CNTS *prev, *end_ptr;

if (clststrt != NULL)
  {
  prev = end_ptr = *clststrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtelt;
    }
  end_ptr = (DM_METH_CNTS *) getmemory(sizeof(DM_METH_CNTS),"cnt elt");
  end_ptr->nxtelt = NULL;
  end_ptr->methcnt = pcnt;
  end_ptr->unmethcnt = mcnt;
  end_ptr->smplgroup = sgroup;
  end_ptr->cntused = 0;
  if (*clststrt == NULL)
    {
    *clststrt = end_ptr;
    end_ptr->prvelt = NULL;
    }
  else
    {
    prev->nxtelt = end_ptr;
    end_ptr->prvelt = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

int dm_cntmethcntelt(DM_METH_CNTS *clst)
  /* recursively count clst */
{
if (clst == NULL)
  return(0);
else
  return(dm_cntmethcntelt(clst->nxtelt)+1);
}

int dm_cntnzmethcntelt(DM_METH_CNTS *clst)
  /* recursively count non-zero elements in clst */
{
if (clst == NULL)
  return(0);
else
  if ((clst->methcnt > 0) || (clst->unmethcnt > 0))
    return(dm_cntnzmethcntelt(clst->nxtelt)+1);
  else
    return(dm_cntnzmethcntelt(clst->nxtelt));
}

void dm_delmethcntelt(DM_METH_CNTS *ep,
                      DM_METH_CNTS **lstrt)
/* delete ep from list *lstrt */
{
DM_METH_CNTS *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvelt) == NULL)
    *lstrt = ep->nxtelt;
  else
    pt->nxtelt = ep->nxtelt;
  if ((pt = ep->nxtelt) != NULL)
    pt->prvelt = ep->prvelt;
  memfree(ep);
  }
}

void dm_clrallmethcntelts(DM_METH_CNTS **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  dm_delmethcntelt(*lstrt,lstrt);
}

DM_CPGPOS_ELT *dm_appndcpgpos(DM_CPGPOS_ELT **cpgposlst,
                              int cgpos,
                              int do3pc)
/* append a new cpgpos element to end of cpgposlst.
return address of new element in case useful.  if do3pc
then append a count on higher, being the C position in
the complementary strand.  Reject 0 or -ve cpgpos entries
since these are a flag to ignore the values. */
{
DM_CPGPOS_ELT *prev, *end_ptr;

if ((cpgposlst != NULL) && (cgpos > 0))
  {
  prev = end_ptr = *cpgposlst;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtcpgelt;
    }
  end_ptr = (DM_CPGPOS_ELT *) getmemory(sizeof(DM_CPGPOS_ELT),"cpgpos elt");
  end_ptr->nxtcpgelt = NULL;
  end_ptr->cpgpos = cgpos;
  if (*cpgposlst == NULL)
    {
    *cpgposlst = end_ptr;
    end_ptr->prvcpgelt = NULL;
    }
  else
    {
    prev->nxtcpgelt = end_ptr;
    end_ptr->prvcpgelt = prev;
    }
  if (do3pc)
    return(dm_appndcpgpos(cpgposlst,cgpos+1,0));
  else
    return(end_ptr);
  }
else
  return(NULL);
}

void dm_delcpgposelt(DM_CPGPOS_ELT *ep,
                     DM_CPGPOS_ELT **lstrt)
/* delete ep from list *lstrt */
{
DM_CPGPOS_ELT *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvcpgelt) == NULL)
    *lstrt = ep->nxtcpgelt;
  else
    pt->nxtcpgelt = ep->nxtcpgelt;
  if ((pt = ep->nxtcpgelt) != NULL)
    pt->prvcpgelt = ep->prvcpgelt;
  memfree(ep);
  }
}

void dm_clrallcpgposelts(DM_CPGPOS_ELT **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  dm_delcpgposelt(*lstrt,lstrt);
}

DM_CPG_BIN *dm_appndbin(DM_CPG_BIN **lstrt,
                        int bstart,
                        int len,
                        int cpgs,
                        DM_RUNPARS *rpp)
/* create and append a new element to *lstrt,
init count to 0.
 Return address of new element */
{
DM_CPG_BIN *prev, *end_ptr;

if (lstrt != NULL)
  {
/*  if (debuglevel > RBC_dbg_none)
    fprintf(stdout,"Bin %d..%d, len=%d\n",bstart,bstart+len-1,len); */
  prev = end_ptr = *lstrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtbin;
    }
  end_ptr = (DM_CPG_BIN *) getmemory(sizeof(DM_CPG_BIN),"bin elt");
  end_ptr->nxtbin = NULL;
  end_ptr->spos = bstart;
  end_ptr->binlen = len;
  end_ptr->cpgcnt = cpgs;
  if (rpp->cntcpgsas1)
    end_ptr->maxccnt = cpgs;
  else
    end_ptr->maxccnt = 2*cpgs;
  end_ptr->cntlst = NULL;
  end_ptr->cpgposlist = NULL;
  if (*lstrt == NULL)
    {
    *lstrt = end_ptr;
    end_ptr->prvbin = NULL;
    }
  else
    {
    prev->nxtbin = end_ptr;
    end_ptr->prvbin = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void rbc_delcntbin(DM_CPG_BIN *ep,
                   DM_CPG_BIN **lstrt)
/* delete ep from list *lstrt */
{
DM_CPG_BIN *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvbin) == NULL)
    *lstrt = ep->nxtbin;
  else
    pt->nxtbin = ep->nxtbin;
  if ((pt = ep->nxtbin) != NULL)
    pt->prvbin = ep->prvbin;
  dm_clrallmethcntelts(&ep->cntlst);
  dm_clrallcpgposelts(&ep->cpgposlist);
  memfree(ep);
  }
}

void rbc_clrallcntbins(DM_CPG_BIN **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  rbc_delcntbin(*lstrt,lstrt);
}

int bc_cntcntbins(DM_CPG_BIN *clst)
  /* recursively read list elements */
{
if (clst == NULL)
  return(0);
else
  return(bc_cntcntbins(clst->nxtbin) + 1);
}

DM_CPG_BIN *cpg_lastcntbin(DM_CPG_BIN *clst)
  /* iterate thru clst, returning
last element, NULL if none */
{
DM_CPG_BIN *ep;

if ((ep = clst) == NULL)
  return(NULL);
else
  {
  while (ep->nxtbin != NULL)
    ep = ep->nxtbin;
  return(ep);
  }
}

int bc_adjacentbins(DM_CPG_BIN *b1,
                    DM_CPG_BIN *b2)
/* return true if b1 & b2 are adjacent with
no intervening gap */
{
if ((b1 == NULL) || (b2 == NULL))
  return(0);
else
  if (b1->spos > b2->spos)
    return(bc_adjacentbins(b2,b1));
  else
    return((b1->spos + b1->binlen) == b2->spos);
}

int bc_sumbinlens(DM_CPG_BIN *blst)
  /* recursively sum blst */
{
if (blst == NULL)
  return(0);
else
  return(bc_sumbinlens(blst->nxtbin) + blst->binlen);
}

int int_in_rng(int b1,
               int x,
               int b2)
/* return 1 if b1 <= x <= b2 or
b2 <= x <= b1 */
{
if (b1 <= b2)
  return((b1 <= x) && (x <= b2));
else
  return(int_in_rng(b2,x,b1));
}

DM_BIN_LSTELT *dm_appndbinlstelt(DM_BIN_LSTELT **lstrt,
                                 int bstart,
                                 DM_CPG_BIN *binp)
/* create and append a new element to *lstrt,
 Return address of new element */
{
DM_BIN_LSTELT *prev, *end_ptr;

if (lstrt != NULL)
  {
  prev = end_ptr = *lstrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtbinelt;
    }
  end_ptr = (DM_BIN_LSTELT *) getmemory(sizeof(DM_BIN_LSTELT),"bin lstelt");
  end_ptr->nxtbinelt = NULL;
  end_ptr->binsstrt = bstart;
  end_ptr->binsstop = 0;
  end_ptr->binptr = binp;
  if (*lstrt == NULL)
    {
    *lstrt = end_ptr;
    end_ptr->prvbinelt = NULL;
    }
  else
    {
    prev->nxtbinelt = end_ptr;
    end_ptr->prvbinelt = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

int dm_scanbinlsts2eltlsts(DM_BIN_LSTELT **beltlst,
                           DM_CPG_BIN *blst,
                           int mkfreq)
/* scan blst, making a new *beltlst element for every
mkfreq-th entry.  Predicated on the basis that blst
is in ascending order. return the number of beltlst
elements created */
{
DM_BIN_LSTELT *newblelt;
DM_CPG_BIN *bp;
int ecnt;
int bcnt;
DM_BIN_LSTELT *prvblelt;
DM_CPG_BIN *prvbin;

ecnt = bcnt = 0;
if (mkfreq > 1)
  {
  bp = blst;
  prvblelt = *beltlst;
  prvbin = NULL;
  while (bp != NULL)
    {
    if (bcnt <= 0)
      {
      newblelt = dm_appndbinlstelt(beltlst,bp->spos,bp);
      if (prvblelt == NULL)
        *beltlst = newblelt;
      else
        prvblelt->binsstop = bp->spos + bp->binlen - 1;
      prvblelt = newblelt;
      bcnt = mkfreq;
      ecnt++;
      }
    bcnt--;
    prvbin = bp;
    bp = bp->nxtbin;
    }
  if ((prvbin != NULL) && (prvblelt != NULL))
    prvblelt->binsstop = prvbin->spos + prvbin->binlen - 1;
  }
return(ecnt);
}

DM_CPG_BIN *bc_cntbin4pos(DM_CPG_BIN *blst,
                          int posn)
/* return a pointer to the bin that contains posn.
NULL if none */
{
DM_CPG_BIN *bp;

bp = blst;
while (bp != NULL)
  if (int_in_rng(bp->spos,posn,(bp->spos + bp->binlen -1)))
    return(bp);
  else
    if (bp->spos > posn)  /* die now: we have passed the thing */
      bp = NULL;
    else
      bp = bp->nxtbin;
/* fell off end, return NULL */
return(NULL);
}

DM_CPG_BIN *dm_binforpositn(DM_CHR_INFO *cbinfop,
                            int posn)
/* use the information in *cbinfop to
return the count bin for posn.  NULL if not */
{
DM_BIN_LSTELT *lep;

if ((lep = cbinfop->beltlst) != NULL)
  {
  while (lep != NULL)
    if (int_in_rng(lep->binsstrt,posn,lep->binsstop))
      return(bc_cntbin4pos(lep->binptr,posn));
    else
      if (lep->binsstrt > posn)
        lep = NULL;
      else
        lep = lep->nxtbinelt;
  return(NULL);
  }
else
  return(bc_cntbin4pos(cbinfop->binlst,posn));
}

int dm_cpgposmatch(DM_RUNPARS *rpp,
                   DM_CNTS4CPG *cntp,
                   int cpos)
/* according to rpp->cntcpgsas1, return 1 if
cpos corresponds to *cntp */
{
if (cntp->cpgpos == cpos)
  return(1);
else
  return((rpp->cntcpgsas1) && (cntp->cpgpos + 1 == cpos));
}

int dm_totbincnts(DM_METH_CNTS *cp,
                  int allcnts)
  /* iteratively sum meth and unmeth counts for cp... */
{
int total;
DM_METH_CNTS *mcp;

total = 0;
mcp = cp;
while (mcp != NULL)
  {
  if (allcnts || mcp->cntused)
    total += mcp->methcnt + mcp->unmethcnt;
  mcp = mcp->nxtelt;
  }
return(total);
}


void dm_startoputline(FILE *ofl,
                      DM_RUNPARS *rpp,
                      DM_CPG_BIN *binp,
                      int chrno)
/* write the first section of an output line to ofl */
{
int totcnt;
float cntspercpg;

totcnt = dm_totbincnts(binp->cntlst,rpp->cpgdetails||(rpp->omode==dm_out_fisherforce));
if (binp->cpgcnt > 0)
  cntspercpg = (float) totcnt/binp->cpgcnt;
else
  cntspercpg = -1.0;
switch (rpp->omode)
  {
  case dm_out_anova:
  case dm_out_anova_gtr:
  case dm_out_anova_more:
    fprintf(ofl,"%s%s%s%d%s%d%s%d%s%d",rpp->chrprefix,
              wlu_intval2str(rpp->chridlu,chrno),rpp->listdelmtr,
              binp->spos,rpp->listdelmtr,
              (binp->spos+binp->binlen-1),rpp->listdelmtr,binp->binlen,
              rpp->listdelmtr,binp->cpgcnt);
    break;
  default:
    fprintf(ofl,"%s%s%s%d%s%d%s%d%s%d%s%d%s%.2f",rpp->chrprefix,
              wlu_intval2str(rpp->chridlu,chrno),rpp->listdelmtr,
              binp->spos,rpp->listdelmtr,
              (binp->spos+binp->binlen-1),rpp->listdelmtr,binp->binlen,
              rpp->listdelmtr,binp->cpgcnt,rpp->listdelmtr,totcnt,
              rpp->listdelmtr,cntspercpg);
    break;
  }
}

int dm_notemethstatinbin(DM_RUNPARS *rpp,
                         DM_CPG_BIN *binp,
                         int chrno,
                         int sqpos,
                         int methstat)
/* note methstat in binp.  return 1 if ok matched to a
real position & details needed */
{
int cpepfnd;
int cpep;

rpp->cachechrno = chrno;
if (methstat)
  binp->cntlst->methcnt++;
else
  binp->cntlst->unmethcnt++;
if (rpp->cpgdetails)
  {
  cpep = 1;
  cpepfnd = 0;
  while (!cpepfnd && (cpep <= binp->maxccnt))
    if (dm_cpgposmatch(rpp,(binp->cntlst->cpgcntarr+cpep),sqpos))
      {
      if (methstat)
        (binp->cntlst->cpgcntarr+cpep)->metcnt++;
      else
        (binp->cntlst->cpgcntarr+cpep)->unmetcnt++;
      cpepfnd = 1;
      }
    else
      cpep++;
  if (cpepfnd)
    return(1);
  else
    {
    if (methstat)
      binp->cntlst->cpgcntarr->metcnt++;
    else
      binp->cntlst->cpgcntarr->unmetcnt++;
/* fprintf(stdout,"Stray at Chr %s %d\n",any_chrno2str(rpp->maxchrno,rpp->havexy,chrno+1,1),sqpos); */
    return(0);
    }
  }
else
  return(1);
}

int dm_notemethstat4chrnpos(DM_RUNPARS *rpp,
                            int chrno,
                            int sqpos,
                            int methstat,
                            int fst3pcpg)
/* the guts of noting methylation status.  Note that
chrno is expected from 0 */
{
DM_CPG_BIN *prevbin;

if (rpp->cachechrno != chrno)
  {
  rpp->currbin = NULL;
  rpp->cachechrno = -1;
  }
if ((rpp->currbin != NULL) &&
     !int_in_rng(rpp->currbin->spos,sqpos,
                   (rpp->currbin->spos + rpp->currbin->binlen -1)))
 {
 rpp->currbin = NULL;
 rpp->cachechrno = -1;
 }
if (rpp->currbin == NULL)
  rpp->currbin = dm_binforpositn((rpp->chrbininfo+chrno),sqpos);
if (rpp->currbin != NULL)
  {
/* dm_startoputline(stdout,rpp,rpp->currbin,chrno); */
  if (fst3pcpg && rpp->samfst3p2prv &&
        (sqpos == rpp->currbin->spos)) /* should accumulate in prev if any */
    {
    if ((prevbin = dm_binforpositn((rpp->chrbininfo+chrno),sqpos-1)) != NULL)
      return(dm_notemethstatinbin(rpp,prevbin,chrno,sqpos,methstat));
    else
      return(0);
    }
  else
    return(dm_notemethstatinbin(rpp,rpp->currbin,chrno,sqpos,methstat));
  }
else
  if (fst3pcpg && rpp->samfst3p2prv)
    {
    if ((prevbin = dm_binforpositn((rpp->chrbininfo+chrno),sqpos-1)) != NULL)
      return(dm_notemethstatinbin(rpp,prevbin,chrno,sqpos,methstat));
    else
      return(0);
    }
  else    
    return(0);
}

int dm_readcpgfl(FILE *sfl,
                 DM_RUNPARS *rpp)
/* use fscanf to read successive lines from sfl.
Expected to contain <Chr>\t<Pos>\t<+/-> */
{
char chrstr[5];
int sqpos;
char sensestr[5];
int matcnt;
int chrno;  /* now 0 based: -1 is invalid */
int scnt;

matcnt = 0;
rpp->cachechrno = -1;
rpp->currbin = NULL;
while ((scnt = fscanf(sfl,"%s %d %s",&chrstr[0],&sqpos,&sensestr[0])) != EOF)
  if ((scnt == 3) &&
        ((chrno = wlu_chkwrd(rpp->chridlu,&chrstr[0])) > -1) &&
        (chrno <= rpp->maxchrno) && ((rpp->uchrno == NULL) || (chrno == wlu_chkwrd(rpp->chridlu,rpp->uchrno))))
      {
      matcnt += dm_notemethstat4chrnpos(rpp,chrno,sqpos,sensestr[0]=='+',0);
      }
return(matcnt);
}

int dm_readctlstfl(FILE *sfl,
                   DM_RUNPARS *rpp)
/* use fscanf to read successive lines from sfl.
Expected to contain

<Chr>\t<Pos>\t<+/->\tCG\t<Ccnt>\t<Tcnt>

disregard lines that don't have CG in 4th field */
{
char chrstr[5];
int sqpos;
char sensestr[5];
int matcnt;
int chrno;  /* now 0 based: -1 is invalid */
int scnt;
int ccnt;
int tcnt;
char basefield[5];

matcnt = 0;
rpp->cachechrno = -1;
rpp->currbin = NULL;
while ((scnt = fscanf(sfl,"%s %d %s %s %d %d",&chrstr[0],&sqpos,&sensestr[0],&basefield[0],
                        &ccnt,&tcnt)) != EOF)
  if ((scnt == 6) && (strncmp(&basefield[0],"CG",2) == 0) &&
        ((chrno = wlu_chkwrd(rpp->chridlu,&chrstr[0])) > -1) &&
        (chrno <= rpp->maxchrno) && ((rpp->uchrno == NULL) || (chrno == wlu_chkwrd(rpp->chridlu,rpp->uchrno))))
      {
      if (sensestr[0] == '-')
        sqpos--;
      while (ccnt > 0)
        {
        matcnt += dm_notemethstat4chrnpos(rpp,chrno,sqpos,1,0);
        ccnt--;
        }
      while (tcnt > 0)
        {
        matcnt += dm_notemethstat4chrnpos(rpp,chrno,sqpos,0,0);
        tcnt--;
        }
      }
return(matcnt);
}

int dm_prntcpginfo(DM_RUNPARS *rpp,
                   int chrno,
                   int sqpos,
                   int methstat,
                   int fst3pcpg)
/* simply write info to stdout in conventional form */
{
fprintf(stdout,"%s%s%d%s%s",wlu_intval2str(rpp->chridlu,chrno),
                  rpp->listdelmtr,sqpos,rpp->listdelmtr,(methstat?"+":"-"));
fputc('\n',stdout);
return(1);
}

int dm_cmpsamread(DM_RUNPARS *rpp,
                  char *sqread,
                  int readlen,
                  int rdpos,
                  int chrno,
                  DM_READ_SENS rdsens,
                  int (*proccpgpos)(DM_RUNPARS *xrpp,
                                    int xchrno,
                                    int xsqpos,
                                    int xmetstat,
                                    int fst3pcpg))
/* compare this sequence read with rdsens, return
cpg count.  chrno now 0 based */
{
char *cspt;
char *rdpt;
int cpgcnt;
int rdlen;
int cpgpos;
int climit;

cpgcnt = 0;
if (((rpp->uchrno == NULL) || (chrno == wlu_chkwrd(rpp->chridlu,rpp->uchrno))) && (chrno >= 0))
  {  /* reject MT */
  rdlen = readlen;
  if (rdsens == DM_rdsens_5p)
    {
    cpgpos = rdpos;
    rdpt = sqread;
    cspt = ((rpp->chrbininfo+chrno)->chromseq+rdpos-1);
    climit = (rpp->chrbininfo+chrno)->chrlen;
    while ((rdlen > 0) && rbc_intinrng(1,cpgpos,climit))
      {
      if ((toupper(*cspt) == 'C') && (rpp->allcs || (toupper(*(cspt+1)) =='G')))
        switch (toupper(*rdpt))
          {
          case 'T':      /* unmeth */
            if((* proccpgpos)(rpp,chrno,cpgpos,0,0))
              cpgcnt++;
            break;
          case 'C':      /* meth */
            if ((* proccpgpos)(rpp,chrno,cpgpos,1,0))
              cpgcnt++;
            break;
          default:
            break;
          }
      cpgpos++;
      rdpt++;
      cspt++;
      rdlen--;
      }
    }
  else
    {
    rdpt = sqread + rdlen - 1;
    if (rpp->allcs)
      {
      cpgpos = rdpos + rdlen - 1;
      cspt = ((rpp->chrbininfo+chrno)->chromseq+cpgpos-1);
      climit = (rpp->chrbininfo+chrno)->chrlen;
      while ((rdpt >= sqread) && rbc_intinrng(1,cpgpos,climit))
        {  /* working in base complement space this time */
        if (toupper(*cspt) == 'C')
          {
          switch (toupper(*rdpt))
            {
            case 'A':   /* unmeth */
              if ((* proccpgpos)(rpp,chrno,cpgpos,0,cpgcnt==0))
                cpgcnt++;
              break;
            case 'G':
              if ((* proccpgpos)(rpp,chrno,cpgpos,1,cpgcnt==0))
                cpgcnt++;
              break;
            default:
              break;
            }
          }
        cpgpos--;
        rdpt--;
        cspt--;
        }
      }
    else
      {
      cpgpos = rdpos + rdlen - 2;
      cspt = ((rpp->chrbininfo+chrno)->chromseq+cpgpos);
      climit = (rpp->chrbininfo+chrno)->chrlen;
      while ((rdpt >= sqread) && rbc_intinrng(1,cpgpos,climit))
        {  /* working in base complement space this time */
        if ((toupper(*cspt) == 'G') && (toupper(*(cspt-1)) =='C'))
          {
          switch (toupper(*rdpt))
            {
            case 'A':   /* unmeth */
              if ((* proccpgpos)(rpp,chrno,cpgpos,0,cpgcnt==0))
                cpgcnt++;
              break;
            case 'G':
              if ((* proccpgpos)(rpp,chrno,cpgpos,1,cpgcnt==0))
                cpgcnt++;
              break;
            default:
              break;
            }
          }
        cpgpos--;
        rdpt--;
        cspt--;
        }
      }
    }
  }
return(cpgcnt);
}

int dm_readsamfl(FILE *sfl,
                 DM_RUNPARS *rpp)
/* sfl is a bismark .sam file: use fscanf to read 
successive lines from it. */
{
char *readbuf;
SAM_BSMK_FLD atmno;
char nc;
char *bp;
DM_READ_SENS rdsens;
int sqpos;
int chrno;
int matcnt;
int lno;

matcnt = 0;
readbuf = (char *) getmemory(rpp->sambuflen+1,"Sam Read buf");
atmno = SAM_fld_header;
bp = readbuf;
/* rpp->srclno = 0; */
lno = 0;
while ((nc = fgetc(sfl)) != EOF)
  switch (nc)
    {
    case '\n':    /* end of line, just reset things */
      lno++;
/*      rpp->srclno++; */
      matcnt += dm_cmpsamread(rpp,readbuf,strlen(readbuf),sqpos,chrno,rdsens,
                                dm_notemethstat4chrnpos);
      atmno = SAM_fld_header;
      bp = readbuf;
      break;
    case '\t':    /* end of atom: decide what we should do next */
      switch (atmno)
        {
        case SAM_fld_sense:
          if ((int)strtol(readbuf,NULL,10) & BAM_FLAG_CMP)
            rdsens = DM_rdsens_3p;
          else
            rdsens = DM_rdsens_5p;
/*          switch ((int)strtol(readbuf,NULL,10))
            {
            case 16:
              rdsens = DM_rdsens_3p;
              break;
            case 0:
            default:
              rdsens = DM_rdsens_5p;
              break;
            } */
          bp = readbuf;
          atmno++;
          break;
        case SAM_fld_chromosome:
          chrno = wlu_chkwrd(rpp->chridlu,readbuf);
/*          chrno = any_str2chrno(rpp->maxchrno,rpp->havexy,readbuf); */
          bp = readbuf;
          atmno++;
          break;
        case SAM_fld_positn:
          sqpos = (int) strtol(readbuf,NULL,10);
          bp = readbuf;
          atmno++;
          break;
        case SAM_fld_header:
          bp = readbuf;
          atmno++;
          break;
        case SAM_fld_seq:     /* got what we want for now, process & bail out this line */
          lno++;
/*          rpp->srclno++; */
          matcnt += dm_cmpsamread(rpp,readbuf,strlen(readbuf),sqpos,chrno,rdsens,
                                    dm_notemethstat4chrnpos);
          bas_skipeol(sfl,NULL);
          atmno = SAM_fld_header;
          bp = readbuf;
          break;
        default:
          bp = readbuf;
          atmno++;
          break;
        }
      break;
    case '@':
      if ((atmno == SAM_fld_header) && (bp == readbuf)) /* still on start of line, bail out */
        {
        lno++;
/*        rpp->srclno++; */
        bas_skipeol(sfl,NULL);
        bp = readbuf;
        break;
        }
    default:
      bas_appchr(readbuf,&bp,nc,rpp->sambuflen);
      break;
    }
memfree(readbuf);
return(matcnt);
}

#ifndef NO_ZLIB

int dm_readbamfl(DM_RUNPARS *rpp)
/* read from previously opened .bam file at rpp->bamrpd.
read successive successive records, return number of matches
*/
{
DM_READ_SENS rdsens;
int chrno;
int matcnt;
BF_EXPNDALIGNREC *earp;
int reccnt;

reccnt = matcnt = 0;
while ((earp = bf_nxtalignrec(rpp->bamrdp,1,err_msg_die)) != NULL)
  {
  reccnt++;
  if (earp->arec->refID >= 0)  /* was this a valid match (not *) */
    {
    chrno = wlu_chkwrd(rpp->chridlu,(rpp->bamrdp->refinfo+earp->arec->refID)->refname);
    if (earp->arec->FLAG & BAM_FLAG_CMP)
      rdsens = DM_rdsens_3p;
    else
      rdsens = DM_rdsens_5p;
/*   switch (earp->arec->FLAG)
      {
      case 16:
        rdsens = DM_rdsens_3p;
        break;
      case 0:
      default:
        rdsens = DM_rdsens_5p;
        break;
      } */
    matcnt += dm_cmpsamread(rpp,earp->seq,earp->arec->l_seq,(earp->arec->pos+1),
                              chrno,rdsens,dm_notemethstat4chrnpos);
    }
  bf_disposeexpalignrec(earp);
  }
return(matcnt);
}

#endif

void dm_resetcntlststohead(DM_CHR_INFO *bnlsts,
                           DM_RUNPARS *rpp)
/* for all used chromosomes, traverse the bin lists, 
reseting the cntlst pointer to the head of the list. */
{
RBC_CHRNO chrno;
DM_CPG_BIN *bp;
DM_METH_CNTS *clp;

for (chrno = 0; chrno < rpp->maxchrno; chrno++)
  if ((rpp->uchrno == NULL) || (wlu_chkwrd(rpp->chridlu,rpp->uchrno) == chrno))
    {
    bp = (bnlsts+chrno)->binlst;
    while (bp != NULL)
      {
      clp = bp->cntlst;
      while ((clp != NULL) && (clp->prvelt != NULL))
        clp = clp->prvelt;
      bp->cntlst = clp;
      bp = bp->nxtbin;
      }
    }
}

int dm_readnxtsrcfl(DM_CHR_INFO *bnlsts,
                    FILE *sfl,
                    DM_RUNPARS *rpp,
                    DM_SRC_MODE srcmod,
                    int sgroup)
/* prepare to read a new source file.  Go through
entire bin set and append new element, making that
the one for counting into.  The forward & reverse
lists heads can be reset later.
If rpp->uchrno is nonNULL, then only do that chromosome */
{
RBC_CHRNO chrno;
DM_CPG_BIN *binptr;
DM_CPGPOS_ELT *cpeptr;
int cp;

for (chrno = 0; chrno < rpp->maxchrno; chrno++)
  if ((rpp->uchrno == NULL) || (wlu_chkwrd(rpp->chridlu,rpp->uchrno) == chrno))
    {
    binptr = (bnlsts+chrno)->binlst;
    while (binptr != NULL)
      {
      binptr->cntlst = dm_appndcntelt(&binptr->cntlst,0,0,sgroup);
      if (rpp->cpgdetails)
        {
/* create an appropriately sized array to store CpG positions & counts.
Create a 0 position for non-matching positions */
        binptr->cntlst->cpgcntarr = (DM_CNTS4CPG *)
          getmemory((binptr->maxccnt+1)*sizeof(DM_CNTS4CPG),"bincpgcntarray");
        cpeptr = binptr->cpgposlist;
        cp = 1;
        binptr->cntlst->cpgcntarr->cpgpos = binptr->cntlst->cpgcntarr->metcnt =
          binptr->cntlst->cpgcntarr->unmetcnt = 0;
        while ((cpeptr != NULL) && (cp <= binptr->maxccnt))
          {
          (binptr->cntlst->cpgcntarr+cp)->cpgpos = cpeptr->cpgpos;
          (binptr->cntlst->cpgcntarr+cp)->metcnt =
            (binptr->cntlst->cpgcntarr+cp)->unmetcnt = 0;
          cp++;
          cpeptr = cpeptr->nxtcpgelt;
          }
        }
      binptr = binptr->nxtbin;
      }
    }
/* rpp->srclno = 0; */
switch (srcmod)
  {
  case DM_src_cpglist:
    return(dm_readcpgfl(sfl,rpp));
    break;
  case DM_src_ctlist:
    return(dm_readctlstfl(sfl,rpp));
    break;
#ifndef NO_ZLIB
  case DM_src_bamfile:
    return(dm_readbamfl(rpp));
    break;
#endif
  case DM_src_samfile:
  default:
    return(dm_readsamfl(sfl,rpp));
    break;
  }
}

char tr_int2nares(int iv)
  /* return a nucleic acid residue for iv */
{
switch (iv)
  {
  case 0:
    return('a');
    break;
  case 1:
    return('c');
    break;
  case 2:
    return('g');
    break;
  case 3:
    return('t');
    break;
  default:
    return('?');
    break;
  }
}

int tr_nares2int(char res)
  /* return an int value 0..3 for res, -1 for unknown */
{
switch (toupper(res))
  {
  case 'A':
    return(0);
    break;
  case 'C':
    return(1);
    break;
  case 'G':
    return(2);
    break;
  case 'T':
  case 'U':
    return(3);
    break;
  default:
    return(-1);
    break;
  }
}

int dm_cntcpgposlst(DM_CPGPOS_ELT *cpglst)
{
if (cpglst == NULL)
  return(0);
else
  return(1+dm_cntcpgposlst(cpglst->nxtcpgelt));
}

void dm_printcpglst(FILE *ofl,
                    DM_CPGPOS_ELT *cpglst)
{
DM_CPGPOS_ELT *cp;

fprintf(ofl,"%d: ",dm_cntcpgposlst(cpglst));
cp = cpglst;
while (cp != NULL)
  {
  fprintf(ofl,"%s%d",(cp==cpglst?"":","),cp->cpgpos);
  cp = cp->nxtcpgelt;
  }
fputc('\n',ofl);
}

DM_CPG_BIN *dm_completefrag(DM_RUNPARS *rpp,
                            DM_CPG_BIN **binlst,
                            RBC_CHRNO chrno,
                            int binstrt,
                            int binlength,
                            int cpgs,
                            DM_CPGPOS_ELT **cpgposlst)
/* append this fragment to binlst, return the new
element address chrno is 0..chrmax-1 - not checked
here */
{
DM_CPG_BIN *newelt;

switch (rpp->binstyle)
  {
  case DM_bin_userlist:
    if ((newelt = dm_binforpositn((rpp->chrbininfo+chrno),binstrt)) != NULL)
      {
      newelt->cpgcnt = cpgs;
      if (rpp->cntcpgsas1)
        newelt->maxccnt = cpgs;
      else
        newelt->maxccnt = 2*cpgs;
      }
    if (rpp->cpgdetails)
      {
      if (rpp->samfst3p2prv)
        (void) dm_appndcpgpos(cpgposlst,(binstrt+binlength),!rpp->cntcpgsas1);
      if (newelt != NULL)
        newelt->cpgposlist = *cpgposlst;
      *cpgposlst = NULL;
      }
    break;
  case DM_bin_restenz:
  case DM_bin_fixed:
  default:
    newelt = dm_appndbin(binlst,binstrt,binlength,cpgs,rpp);
    if (rpp->cpgdetails)
      {
      if (rpp->samfst3p2prv)
        (void) dm_appndcpgpos(cpgposlst,(binstrt+binlength),!rpp->cntcpgsas1);
      if (cpgposlst != NULL)
        {
        newelt->cpgposlist = *cpgposlst;
        *cpgposlst = NULL;
	}
      }
    if ((rpp->chrbininfo+chrno)->binlst == NULL)
      (rpp->chrbininfo+chrno)->binlst = newelt;
    break;
  }
return(newelt);
}

int dm_redrepgenomsqs(DM_RUNPARS *rpp,
                      FS_FSMSTRCT *cpgfsmp)
/* for each chromosome, open each as a Fasta sequence
file.  Scan with cpgfsmp (CCGG) and create bins, giving
a count of CpG found with CG.  Return the
number of chromosomes processed. */
{
int chno;
SQFL_STRCT *chsqfl;
int rcnt;
char nxtres;
int cpos;
FS_RESELT *frp;
FS_DATPRELT *dpep;
int thismat;
int prvmat;
int thislen;
DM_CPG_BIN *binsend;
int cpgs;
DM_CPGPOS_ELT *cpgposlst;
int prvcpgpos;
char *sqp;
off_t filelen;
DM_CPG_BIN *fxdbinp;
DM_CPG_BIN *fxdbinlst;

rcnt = 0;
for (chno = 0; chno < rpp->maxchrno; chno++)
  if ((rpp->uchrno == NULL) || (wlu_chkwrd(rpp->chridlu,rpp->uchrno) == chno))
    {
    if ((chsqfl = sqfl_opnsqstrct((rpp->chrbininfo+chno)->chrflname,SFMT_fasta,"r")) != NULL)
      {
      filelen = sqfl_filelength(chsqfl->sfl);
      if (rpp->needseq)
        sqp = (rpp->chrbininfo+chno)->chromseq =
          (char *) getmemory((size_t) filelen+1,"Chromosomebuff");
      sqfl_rewind(chsqfl);
      (void) sqfl_skipsqflhdr(chsqfl);
      fs_initrun(cpgfsmp);
      cpos = 0;
      prvmat = 1;
      binsend = NULL;
      if (rpp->binstyle == DM_bin_restenz)
        cpgs = -1;     /* avoid distortion since fsm finds CG before CCGG */
      else
        cpgs = 0;
      cpgposlst = NULL;
      prvcpgpos = 0;
/* make bin set for fixed windows, if this is bin style */
      fxdbinlst = NULL;
      if ((rpp->binstyle == DM_bin_fixed) && (rpp->binwidth > 0))
        {
	cpos = 1;
	while (cpos <= filelen)
	  {
	  binsend = dm_completefrag(rpp,&binsend,chno,cpos,rpp->binwidth,0,NULL);
	  if (fxdbinlst == NULL)
	    fxdbinlst = binsend;
	  cpos += rpp->binwidth;
	  }
        thislen = 0;
	}
      cpos = prvcpgpos = 0;
      fxdbinp = NULL;
      while ((nxtres = sqfl_getnxtres(chsqfl)) != '\0')
        {
        cpos++;
        thislen++;
        (rpp->chrbininfo+chno)->chrlen = cpos;
        if (rpp->needseq)
          {
          *sqp = nxtres;
          sqp++;
          *sqp = '\0';
          }
        if ((frp = fs_procchr(cpgfsmp,nxtres,tr_nares2int)) != NULL)
          {
          dpep = (FS_DATPRELT *) frp->action;
          if (debuglevel > RBC_dbg_none)
             fprintf(stdout,"Match: %s @ %d (prev=%d)\n",dpep->dstrng,cpos,prvcpgpos);
          switch (dpep->ldstr)
            {
            case 2: /* CpG */
              if (rpp->cpgdetails)
                {
                if (rpp->binstyle == DM_bin_fixed)
		  {  /* need to get pointer to correct bin */
                  if ((fxdbinp == NULL) || !int_in_rng(fxdbinp->spos,prvcpgpos,(fxdbinp->spos+fxdbinp->binlen-1)))
                    {
                    if (fxdbinp == NULL)
		      fxdbinp = bc_cntbin4pos(fxdbinlst,prvcpgpos);
		    else
		      fxdbinp = bc_cntbin4pos(fxdbinp,prvcpgpos);
		    }
                  if (fxdbinp != NULL)
                    fxdbinp->maxccnt++;
		  (void) dm_appndcpgpos(&fxdbinp->cpgposlist,prvcpgpos,!rpp->cntcpgsas1);
                  }
		else
                  (void) dm_appndcpgpos(&cpgposlst,prvcpgpos,!rpp->cntcpgsas1);
                }
              prvcpgpos = cpos-1;
              cpgs++;
              break;
            case 1: /* any C */
              if (rpp->cpgdetails)
                {
                if (rpp->binstyle == DM_bin_fixed)
		  {  /* need to get pointer to correct bin */
                  if ((fxdbinp == NULL) || !int_in_rng(fxdbinp->spos,prvcpgpos,(fxdbinp->spos+fxdbinp->binlen-1)))
                    {
                    if (fxdbinp == NULL)
		      fxdbinp = bc_cntbin4pos(fxdbinlst,prvcpgpos);
		    else
		      fxdbinp = bc_cntbin4pos(fxdbinp,prvcpgpos);
		    }
                  if (fxdbinp != NULL)
                    fxdbinp->maxccnt++;
		  (void) dm_appndcpgpos(&fxdbinp->cpgposlist,prvcpgpos,!rpp->cntcpgsas1);
                  }
		else
                  (void) dm_appndcpgpos(&cpgposlst,prvcpgpos,!rpp->cntcpgsas1);
                }
              prvcpgpos = cpos;
              cpgs++;
              break;
            default: /* some other site */
              thismat = cpos-dpep->ldstr + ((int) dpep->stxt) + 1;
              if (rpp->binwidth == 0)
                thislen = thismat - prvmat;
              if (rpp->binstyle != DM_bin_fixed)
                {
                binsend = dm_completefrag(rpp,&binsend,chno,prvmat,thislen,cpgs,&cpgposlst);
                if (rpp->cpgdetails)
                  prvcpgpos = thismat;
                if (rpp->cpgdetails && rpp->samfst3p2prv)
                  cpgs = 1;
                else
                  cpgs = 0;
		}
              prvmat = thismat;
              break;
            }
          }
        }
      if (rpp->binstyle != DM_bin_fixed)
        binsend = dm_completefrag(rpp,&binsend,chno,prvmat,thislen,cpgs,&cpgposlst);
      else
        {
	fxdbinp = (rpp->chrbininfo+chno)->binlst;
	while (fxdbinp != NULL)
	  {
	  fxdbinp->cpgcnt = fxdbinp->maxccnt;
	  fxdbinp = fxdbinp->nxtbin;
	  }
	}
      sqfl_clssqstrct(chsqfl);
      rcnt++;
      }
    else
      err_msg("Can't open chromosome file %s\n",(rpp->chrbininfo+chno)->chrflname);
    }
return(rcnt);
}

double bin_coeff(int n,
                 int k)
/* return as double the binomial
coefficient of n,k.  Multiplicative
formula */
{
double prodct;
int i;

prodct = 1.0;
for (i = 1; i <= k; i++)
  prodct *= (double) (n - k + i)/((double) i);
return(prodct);
}

double fishers_exact_t1(int a,
                        int b,
                        int c,
                        int d)
/* return a value for Fisher's exact
statistic - one tailed */
{
double rval;
int tot;

tot = a + b + c + d;
rval = (bin_coeff(a+b,a) / bin_coeff(tot,a+c)) * bin_coeff(c+d,c);
return(rval);
}

void bc_commalst2ints(char *str,
                      int *v1,
                      int *v2)
/* expect str to contain two base10 integers separated by a comma.
attempt to read these and return as v1 & v2.  zero return on
failures */
{
char *ep;
char *fp;

*v1 = *v2 = 0;
*v1 = (int) strtol(str,&ep,10);
if (ep != str)
  *v2 = strtol(ep+1,&fp,10);
}

int bc_fragsizeok(int fragmin,
                  int fragmax,
                  int fragsize)
/* return if fragsize fragment is acceptable,
noting that zero limits will always return true. */
{
if ((fragmin != 0) && (fragmax != 0))
  return(int_in_rng(fragmin,fragsize,fragmax));
else
  return(1);
}

int bc_cpglmttest(int cpgmin,
                  int cpgcnt)
/* return 1 if cpgcnt passes cpgmin, noting
that 0 cpgmin will always return true */
{
if (cpgmin <= 0)
  return(1);
else
  return(cpgcnt >= cpgmin);
}

int bc_binundercpglmt(DM_CPG_BIN *bp,
                      int cpgmin)
/* apply minimum CpG count test to bin, noting
that cgpmin <= 0 will always fail */
{
if (bp != NULL)
  return(!bc_cpglmttest(cpgmin,bp->cpgcnt));
else
  return(1);
}

int bc_prunebins4len(DM_CPG_BIN **binlsts,
                     int minfrag,
                     int maxfrag)
/* scan all or selected chromosome bin lists, removing all
bins which don't conform to minfrag..maxfrag (inclusive).
Return count of deleted bins */
{
DM_CPG_BIN *bp;
int delcnt;
DM_CPG_BIN *nxt;

delcnt = 0;
bp = *binlsts;
while (bp != NULL)
  {
  nxt = bp->nxtbin;
  if (!bc_fragsizeok(minfrag,maxfrag,bp->binlen))
    {
/* fprintf(stdout,"Del: %d..%d (%d), CpGs=%d\n",bp->spos,bp->spos + bp->binlen -1,bp->binlen,bp->cpgcnt); */
    rbc_delcntbin(bp,binlsts);
    delcnt++;
    }
  bp = nxt;
  }
return(delcnt);
}

int bc_prunebinslenok(DM_CPG_BIN **binlsts,
                      int minfrag,
                      int maxfrag)
/* scan all or selected chromosome bin lists, removing all
bins which do conform to minfrag..maxfrag (inclusive).
Return count of deleted bins */
{
DM_CPG_BIN *bp;
int delcnt;
DM_CPG_BIN *nxt;

delcnt = 0;
bp = *binlsts;
while (bp != NULL)
  {
  nxt = bp->nxtbin;
  if (bc_fragsizeok(minfrag,maxfrag,bp->binlen))
    {
/*      if (debuglevel > RBC_dbg_none)
        fprintf(stdout,"Del: %d..%d (%d), CpGs=%d\n",bp->spos,bp->spos + bp->binlen -1,bp->binlen,bp->cpgcnt); */
    rbc_delcntbin(bp,binlsts);
    delcnt++;
    }
  bp = nxt;
  }
return(delcnt);
}

int bc_prunebins4cpgmin(DM_CPG_BIN **binlsts,
                        int mincpg)
/* scan all or selected chromosome bin lists, removing all
bins which don't meet mincpg criterion.
Return count of deleted bins */
{
DM_CPG_BIN *bp;
int delcnt;
DM_CPG_BIN *nxt;

delcnt = 0;
bp = *binlsts;
while (bp != NULL)
  {
  nxt = bp->nxtbin;
  if (bc_binundercpglmt(bp,mincpg))
    {
/*      if (debuglevel > RBC_dbg_none)
        fprintf(stdout,"Del: %d..%d (%d), CpGs=%d\n",bp->spos,bp->spos + bp->binlen -1,bp->binlen,bp->cpgcnt); */
    rbc_delcntbin(bp,binlsts);
    delcnt++;
    }
  bp = nxt;
  }
return(delcnt);
}

int bc_losefailedbins(DM_RUNPARS *rpp,
                      DM_CHR_INFO *binlsts,
                      RBC_CHRNO maxchr,
                      int minfrag,
                      int maxfrag,
                      int cpgmin)
/* scan all or selected chromosome bin lists, removing all
bins which do conform to minfrag..maxfrag (inclusive) and
fail any cpgmin criterion.
Return count of deleted bins */
{
RBC_CHRNO cp;
int retv;

retv = 0;
for (cp = Chr1; cp <= maxchr; cp++)
  {
  if (rpp->binwidth == 0)
    retv += bc_prunebins4len(&(binlsts+cp-1)->binlst,minfrag,maxfrag);
  retv += bc_prunebins4cpgmin(&(binlsts+cp-1)->binlst,cpgmin);
  }
return(retv);
}

int bc_hypometh(DM_METH_CNTS *cp)
  /* return 1 if this is not hypomethylated:
The requirement is stated by Li, et al. (2010) PLOSBiology,11,e1000533,
but is not actually defined by them.  Let's try some simple scheme
for now */
{
return(cp->methcnt <= 0);
}

float dm_methpropnoftot(DM_METH_CNTS *bp)
  /* work out methylation proportion of
total counts this count element.  -1.0 if
total is zero */
{
int ctot;

if ((bp != 0) && ((ctot = bp->methcnt + bp->unmethcnt) > 0))
  return((float)bp->methcnt/ctot);
else
  return(-1.0);
}

float dm_methfraction(DM_METH_CNTS *bp)
  /* methylation fraction - return -1.0 if
can't compute */
{
if ((bp != NULL) && (bp->unmethcnt > 0))
  return((float)bp->methcnt/bp->unmethcnt);
else
  return(-1.0);
}

int dm_foldmethdiff(DM_METH_CNTS *b1p,
                    DM_METH_CNTS *b2p,
                    float fold,
                    float (*methpropfn)(DM_METH_CNTS *x1p))
/* test if the methylation difference between
b1p & b2p is >= fold-fold.  (*methpropfn)() is used
to return the methylation proportion, -1.0 indicating that
the quantity couldn't be computed */
{
float b1prop;
float b2prop;
float flddiff;

if (fold == 0.0)
  return(1);
else
  if (((b1prop = (*methpropfn)(b1p)) >= 0.0) &&
      ((b2prop = (*methpropfn)(b2p)) > 0.0))
    {
    flddiff = b1prop/b2prop;
    return((flddiff >= fold) || (flddiff <= 1.0/fold));
    }
  else
    return(0);
}

int dm_chkmincntcriteria(DM_RUNPARS *rpars,
                         DM_METH_CNTS *c1p)
/* see if c1p counts meet any criteria
set in rpars */
{
if (c1p == NULL)
  return(0);
else
  return((rpars->cntthreshold == 0) ||
          ((c1p->methcnt+c1p->unmethcnt) >= rpars->cntthreshold));
}

int dm_pairwisecnt(int i)
  /* return the number of pairwise comparisons for i items */
{
return((int) i*(i-1)/2);
}

double dm_ascendsort(double fvals[],
                     int itemcnt)
/* do a simple bubble sort of fvals, leaving in ascending order.
return first value (least) */
{
int p1;
int p2;
int scnt;
double tmp;

scnt = itemcnt;
while (scnt > 0)
  {
  p1 = 0;
  p2 = 1;
  scnt = 0;
  while (p2 < itemcnt)
    {
    if (fvals[p2] < fvals[p1])
      {
      tmp = fvals[p1];
      fvals[p1] = fvals[p2];
      fvals[p2] = tmp;
      scnt++;
      }
    p1++;
    p2++;
    }
  }
return(fvals[0]);
}

int dm_summethnunmeth(DM_METH_CNTS *clstp,
                      int *methtot,
                      int *unmethtot)
/* traverse clstp, summing meth & unmeth counts.
return number of elements summed */
{
DM_METH_CNTS *cp;
int ecnt;

ecnt = *methtot = *unmethtot = 0;
cp = clstp;
while (cp != NULL)
  {
  *methtot += cp->methcnt;
  *unmethtot += cp->unmethcnt;
  ecnt++;
  cp = cp->nxtelt;
  }
return(ecnt);
}

int dm_sumusedmethnunmeth(DM_METH_CNTS *clstp,
                          int *methtot,
                          int *unmethtot)
/* traverse clstp, summing meth & unmeth counts
for used bins.
return number of elements summed */
{
DM_METH_CNTS *cp;
int ecnt;

ecnt = *methtot = *unmethtot = 0;
cp = clstp;
while (cp != NULL)
  {
  if (cp->cntused)
    {
    *methtot += cp->methcnt;
    *unmethtot += cp->unmethcnt;
    ecnt++;
    }
  cp = cp->nxtelt;
  }
return(ecnt);
}

float dm_minlstcontingency(DM_METH_CNTS *clstp)
  /* traverse clstp, returning minimum contingency
expected count */
{
int methtot;
int unmethtot;
float cmin;
DM_METH_CNTS *cp;
float methprop;
float unmethprop;
float exp;
float thistot;

(void) dm_summethnunmeth(clstp,&methtot,&unmethtot);
cmin = FLT_MAX;
methprop = (float) methtot/(methtot+unmethtot);
unmethprop = 1.0 - methprop;
cp = clstp;
while (cp != NULL)
  {
  thistot = (float) (cp->methcnt + cp->unmethcnt);
  if ((exp = (methprop * thistot)) < cmin)
    cmin = exp;
  if ((exp = (unmethprop * thistot)) < cmin)
    cmin = exp;
  cp = cp->nxtelt;
  }
return(cmin);
}

int dm_min2x2contingency(DM_METH_CNTS *c1p,
                         DM_METH_CNTS *c2p)
  /* return minimum contingency
expected count for c1p & c2p */
{
int methtot;
int unmethtot;
float cmin;
float methprop;
float unmethprop;
float exp;
float thistot;

methtot = c1p->methcnt + c2p->methcnt;
unmethtot = c1p->unmethcnt + c2p->unmethcnt;
cmin = FLT_MAX;
methprop = (float) methtot/(methtot+unmethtot);
unmethprop = 1.0 - methprop;
thistot = (float) (c1p->methcnt + c2p->unmethcnt);
if ((exp = (methprop * thistot)) < cmin)
  cmin = exp;
if ((exp = (unmethprop * thistot)) < cmin)
  cmin = exp;
thistot = (float) (c2p->methcnt + c2p->unmethcnt);
if ((exp = (methprop * thistot)) < cmin)
  cmin = exp;
if ((exp = (unmethprop * thistot)) < cmin)
  cmin = exp;
return(cmin);
}

int dm_chkhypomethcriteria(DM_RUNPARS *rpp,
                           DM_METH_CNTS *clstp)
/* */
{
int totc;

totc = clstp->methcnt + clstp->unmethcnt;
return((rpp->minmethprop <= 0.0) || ((totc > 0) && ((float)clstp->methcnt/totc >= rpp->minmethprop)));
}

int dm_chkcnts4critrn(DM_RUNPARS *rpp,
                      DM_CNTCHK_TYPE chk,
                      DM_CPG_BIN *binp,
                      DM_METH_CNTS *cntp)
/* apply indicated test to cntp
returning  0 if failed */
{
float hitspercpg;
int indcnt;
int passcpgcrit;
int cpgno;
DM_CNTS4CPG *cpgcp;
DM_METH_CNTS *cp;
float minmethprop;
float maxmethprop;
int totcnt;
float mprop;

hitspercpg = 0.0;
switch (chk)
  {
  case dm_chk_cntcrit:
    if (!dm_chkmincntcriteria(rpp,cntp))
      return(0);
    break;
  case dm_chk_cpghits:
  case dm_chk_maxcpghits:
    if (binp->cpgcnt <= 0)
      return(0);
    else
      {
      if (hitspercpg <= 0.0)
        hitspercpg = (float) (cntp->methcnt + cntp->unmethcnt)/binp->cpgcnt;
      if ((chk == dm_chk_cpghits) && (rpp->minhitspercpg > 0.0) &&
           (hitspercpg < rpp->minhitspercpg))
        return(0);
      if ((chk == dm_chk_maxcpghits) && (rpp->maxhitspercpg > 0.0) &&
            (hitspercpg > rpp->maxhitspercpg))
        return(0);
      }
    break;
  case dm_chk_minmeth:
    if (!dm_chkhypomethcriteria(rpp,cntp))
      return(0);
    break;
  case dm_chk_cpgcnt:
    if ((rpp->mincpgs > 0) && (binp->cpgcnt <= rpp->mincpgs))
      return(0);
    break;
  case dm_chk_indcnt:
/*     if (rpp->binprfun != NULL)
      (void) (*rpp->binprfun)(rpp,binp,&indcnt);
    else */
    indcnt = dm_cntnzmethcntelt(binp->cntlst) - 1;
    if (rpp->minsamples > (indcnt + 1))
      return(0);
    break;
  case dm_chk_mincpgcrit:
    passcpgcrit = 0;
    cpgno = 1;   /* since stray counts in offset 0 */
    while (cpgno <= binp->maxccnt)
      {
      cpgcp = cntp->cpgcntarr+cpgno;
      if ((rpp->cntthreshold > 0) &&
           ((cpgcp->metcnt+cpgcp->unmetcnt) >= rpp->cntthreshold))
        passcpgcrit++;
      cpgno++;
      }
    if ((passcpgcrit < rpp->cpgcntscrit) && (passcpgcrit < binp->maxccnt) &&
          (rpp->cpgcntscrit > 0))  /* must allow maxccpg to be less than criterion */
      return(0);
    break;
  case dm_chk_folddiff:
    cp = binp->cntlst;
    minmethprop = 1.0;
    maxmethprop = 0.0;
    while (cp != NULL)
      {
      totcnt = cp->methcnt + cp->unmethcnt;
      if (totcnt != 0)
        {
        mprop =(float) cp->methcnt/totcnt;
        if (mprop > maxmethprop)
          maxmethprop = mprop;
        if (mprop < minmethprop)
          minmethprop = mprop;
        }
      cp = cp->nxtelt;
      }
    if (!((minmethprop == 0.0) ||
         ((rpp->folddiff > 1.0) && (maxmethprop/minmethprop >= rpp->folddiff)) ||
         (((rpp->folddiff != 0.0) && (rpp->folddiff < 1.0)) &&
            (maxmethprop/minmethprop <= rpp->folddiff))))
      return(0);
    break;
  default:
    return(0);
    break;
  }
/* made it right through */
return(1);
}

int dm_chkcnts4mask(DM_RUNPARS *rpp,
                    int chkmask,
                    DM_CPG_BIN *binp,
                    DM_METH_CNTS *cntp)
/* depending on chkmask, apply any indicated tests to cntp
returning 0 if failed */
{
DM_CNTCHK_TYPE cct;

cct = dm_chk_cntcrit;
while (cct <= dm_chk_folddiff)
  {
  if (chkmask & 1 << cct)
    if (!dm_chkcnts4critrn(rpp,cct,binp,cntp))
      return(0);
  cct++;
  }
/* made it right through */
return(1);
}

int dm_cntsarevalid(DM_RUNPARS *rpp,
                    DM_CPG_BIN *binp,
                    DM_METH_CNTS *cntp)
/* depending on rpp->binchkmask, apply any indicated tests to cntp
returning 0 if failed */
{
return(dm_chkcnts4mask(rpp,rpp->binchkmask,binp,cntp));
}

double dm_chisq4cntlst(DM_RUNPARS *rpp,
                       DM_CPG_BIN *binp,
                       DM_METH_CNTS *clstp,
                       int *df)
/* traverse clstp, performing Chi sq evaluation
of the meth and unmeth counts, based on the total
proportions.  return degrees of freedom to df*/
{
int methtot;
int unmethtot;
DM_METH_CNTS *cp;
float methprop;
float unmethprop;
float exp;
float thistot;
double chisq;
float ominuse;
int dfcnt;

(void) dm_summethnunmeth(clstp,&methtot,&unmethtot);
chisq = 0.0;
methprop = (float) methtot/(methtot+unmethtot);
unmethprop = 1.0 - methprop;
cp = clstp;
dfcnt = 0;
while (cp != NULL)
  {
  if (dm_chkcnts4mask(rpp,(rpp->binchkmask&(~1<<dm_chk_indcnt)),binp,cp))
    {  /* avoid pr-based calculation here */
    thistot = (float) (cp->methcnt + cp->unmethcnt);
    if (thistot > 0)
      {
      exp = methprop * thistot;
      ominuse = (float) cp->methcnt - exp;
      if (exp > 0.0)
        chisq += ominuse*ominuse/exp;
      exp = unmethprop * thistot;
      ominuse = (float) cp->unmethcnt - exp;
      if (exp > 0.0)
        chisq += ominuse*ominuse/exp;
      dfcnt++;
      }
    cp->cntused = 1;
    }
  cp = cp->nxtelt;
  }
if (df != NULL)
  *df = dfcnt - 1;
return(chisq);
}

double dm_chisq4cnt2x2(DM_METH_CNTS *c1p,
                       DM_METH_CNTS *c2p)
/* performing Chi sq evaluation for c1p & c2p
of the meth and unmeth counts, based on the total
proportions. */
{
int methtot;
int unmethtot;
float methprop;
float unmethprop;
float exp;
float thistot;
double chisq;
float ominuse;

methtot = c1p->methcnt + c2p->methcnt;
unmethtot = c1p->unmethcnt + c2p->unmethcnt;
chisq = 0.0;
methprop = (float) methtot/(methtot+unmethtot);
unmethprop = 1.0 - methprop;
thistot = (float) (c1p->methcnt + c1p->unmethcnt);
exp = methprop * thistot;
ominuse = (float) c1p->methcnt - exp;
if (exp > 0.0)
  chisq += ominuse*ominuse/exp;
exp = unmethprop * thistot;
ominuse = (float) c1p->unmethcnt - exp;
if (exp > 0.0)
  chisq += ominuse*ominuse/exp;
ominuse = (float) c2p->methcnt - exp;
if (exp > 0.0)
  chisq += ominuse*ominuse/exp;
exp = unmethprop * thistot;
ominuse = (float) c2p->unmethcnt - exp;
if (exp > 0.0)
  chisq += ominuse*ominuse/exp;
return(chisq);
}

int dm_chkpercpgcnts(DM_RUNPARS *rpp,
                     DM_CPG_BIN *binp,
                     DM_METH_CNTS *dmcp)
/* check that the counts in this count element
conform to min and max hits/cpg criteria */
{
float hitspercpg;

if ((binp != NULL) && (binp->cpgcnt > 0) && (dmcp != NULL))
  {
  hitspercpg = (float) (dmcp->methcnt + dmcp->unmethcnt)/binp->cpgcnt;
  return(((rpp->minhitspercpg == 0.0) || (hitspercpg >= rpp->minhitspercpg)) &&
         ((rpp->maxhitspercpg == 0.0) || (hitspercpg < rpp->maxhitspercpg)));
  }
else
  return(0);
}

double dm_stddev4nonngtv(DM_RUNPARS *rpp,
                         DM_CPG_BIN *binp,
                         DM_METH_CNTS *clstp,
                         double (* proptnfn)(DM_METH_CNTS *xp),
                         double *mean,
                         int *np)
/* scan clstp methylation proportions, returning
std deviation.  proptnfn() returns a double value
representing the proportion for the
meth vs unmeth counts, -1 if non-valid  */
{
double sumx;
double sumxsq;
int n;
DM_METH_CNTS *cp;
double proptn;

sumx = sumxsq = 0.0;
n = 0;
cp = clstp;
while (cp != NULL)
  {
  if ((dm_chkpercpgcnts(rpp,binp,cp)) && dm_cntsarevalid(rpp,binp,cp))
    {
    proptn = (* proptnfn)(cp);
    if (proptn >= 0.0)
      {
      sumx += proptn;
      sumxsq += proptn*proptn;
      cp->cntused = 1;
      n++;
      }
    }
  cp = cp->nxtelt;
  }
if (np != NULL)
  *np = n;
if ((mean != NULL) && (n > 0))
  *mean = sumx/(double) n;
if (n > 1)
  return(sqrt((sumxsq-(sumx*sumx/((double) n)))/((double) n-1)));
else
  return(-1);
}

char *dm_cpg_or_c(DM_RUNPARS *rpp)
  /* For tidiness of listings for non-CpG work */
{
if (rpp->allcs)
  return("C");
else
  return("CpG");
}

void dm_commonhdr(FILE *ofl,
                  DM_RUNPARS *rpars)
/* the sections of header for various outputs */
{
switch (rpars->omode)
  {
  case dm_out_anova:
  case dm_out_anova_gtr:
  case dm_out_anova_more:
    fprintf(ofl,"#Chr%sStart%sEnd%sLen%s%ss",rpars->listdelmtr,
              rpars->listdelmtr,rpars->listdelmtr,rpars->listdelmtr,
              dm_cpg_or_c(rpars));
    break;
  default:
    fprintf(ofl,"#Chr%sStart%sEnd%sLen%s%ss%s+&-Hits%s+-hits/%s",rpars->listdelmtr,
              rpars->listdelmtr,rpars->listdelmtr,rpars->listdelmtr,
              dm_cpg_or_c(rpars),rpars->listdelmtr,rpars->listdelmtr,
              dm_cpg_or_c(rpars));
    break;
  }
}

void dm_headdmethout(FILE *ofl,
                     DM_RUNPARS *rpars)
/* header line, for differential output */
{
dm_commonhdr(ofl,rpars);
switch (rpars->omode)
  {
  case dm_out_cpg_detail:
  case dm_out_cpg_nzdetl:
    fputc('\n',ofl);
    break;
  case dm_out_sdev:
    fprintf(ofl,"%sMean%sStddev%sNvalid\n",rpars->listdelmtr,rpars->listdelmtr,
              rpars->listdelmtr);
    break;
  case dm_out_anova:
    fprintf(ofl,"%sPr%sTest\n",rpars->listdelmtr,rpars->listdelmtr);
    break;
  case dm_out_anova_gtr:
  case dm_out_anova_more:
    fprintf(ofl,"%sPr%sTest%s>Meth%sSample_counts%sPropMeth\n",rpars->listdelmtr,
              rpars->listdelmtr,rpars->listdelmtr,rpars->listdelmtr,
              rpars->listdelmtr);
    break;
  case dm_out_chisq:
  case dm_out_chiforce:
    fprintf(ofl,"%sPr%sTest\n",rpars->listdelmtr,rpars->listdelmtr);
    break;
  default:
    fprintf(ofl,"%sPr%sTest",rpars->listdelmtr,rpars->listdelmtr);
    if (rpars->maxgrp > 0)
      fprintf(ofl,"%s>Meth%sPropMeth%sFold_diff",rpars->listdelmtr,rpars->listdelmtr,
                rpars->listdelmtr);
    fputc('\n',ofl);
    break;
  }
}

int dm_chkfemethcriteria(DM_RUNPARS *rpp,
                         DM_CPG_BIN *binp,
                         DM_METH_CNTS *c1p,
                         DM_METH_CNTS *c2p)
/* apply the checks for diff methylation for
Fisher's Exact test */
{
/* WGBS doesn-t set cpgcnts, so test this separately */
if (rpp->binwidth > 0)
  return(((c1p->methcnt + c1p->unmethcnt) > 0) &&
           ((c2p->methcnt + c2p->unmethcnt) > 0) &&
           dm_cntsarevalid(rpp,binp,c1p) && dm_cntsarevalid(rpp,binp,c2p));
else
  return(((binp->cpgcnt >= 5) && (c1p->unmethcnt >= 0) && (c2p->unmethcnt >= 0)
          && dm_cntsarevalid(rpp,binp,c1p) && dm_cntsarevalid(rpp,binp,c2p)));
}

double dm_checkprval(double pr)
  /* check for NaN - return 0 if
found */
{
if (isnan(pr))
  return(0.0);
else
  return(pr);
}

void dm_scandiffmethcnts(FILE *ofl,
                         DM_RUNPARS *rpars,
                         RBC_CHRNO chrno,
                         DM_CPG_BIN *binp,
                         DM_METH_CNTS *cnt1p,
                         DM_METH_CNTS *cnt2p,
                         int p1no,
                         int p2no)
/* scan counts cntNp (N=1 & 2) which belong to the
same bin (binp).  Perform tests of
Li, et al. (2010) PLOSBiology,11,e1000533 for
differential methylation - print results to ofl */
{
float c1prop;
float c2prop;
char *gtrgrp;

if (dm_chkfemethcriteria(rpars,binp,cnt1p,cnt2p))
  {
  dm_startoputline(ofl,rpars,binp,chrno);
  if (rpars->displycnts)
    fprintf(ofl,"%ss\t%d: %d+ %d-\t%d: %d+ %d-",p1no,dm_cpg_or_c(rpars),
              cnt1p->methcnt,cnt1p->unmethcnt,p2no,
              cnt2p->methcnt,cnt2p->unmethcnt);
  fprintf(ofl,"\t%g\tFE",dm_checkprval(fishers_exact_t1(cnt1p->methcnt,cnt1p->unmethcnt,
                                                         cnt2p->methcnt,cnt2p->unmethcnt)));
  if (cnt1p->smplgroup != cnt2p->smplgroup)  /* separate groups, indicate > meth group */
    {
    c1prop = dm_methpropnoftot(cnt1p);
    c2prop = dm_methpropnoftot(cnt2p);
    if ((c1prop >= 0.0) && (c2prop >= 0.0))
      if (c1prop == c2prop)
        gtrgrp = "=";
      else
        if (((c1prop > c2prop) && (cnt1p->smplgroup == 1)) ||
            ((c1prop < c2prop) && (cnt2p->smplgroup == 1)))
          gtrgrp = *(rpars->groupidlist);
        else
          gtrgrp = *(rpars->groupidlist+1);
    else
      gtrgrp = "-";
    fprintf(ofl,"%s%s",rpars->listdelmtr,gtrgrp);
    if ((c1prop <= 0.0) || (c2prop <= 0.0))
      {
      if ((c1prop <= 0.0) && (c2prop > 0.0))
        fprintf(ofl,"%s%s=-,%s=%.4f",rpars->listdelmtr,
              *(rpars->groupidlist),*(rpars->groupidlist+1),
              c2prop);
      else
        if ((c2prop <= 0.0) && (c1prop > 0.0))
          fprintf(ofl,"%s%s=%.4f,%s=-",rpars->listdelmtr,
                *(rpars->groupidlist),c1prop,*(rpars->groupidlist+1));
      fprintf(ofl,"%s-",rpars->listdelmtr);
      }
    else
      {
      fprintf(ofl,"%s%s=%.4f,%s=%.4f",rpars->listdelmtr,
                *(rpars->groupidlist),
                (cnt1p->smplgroup==1?c1prop:c2prop),
		*(rpars->groupidlist+1),
                (cnt2p->smplgroup!=1?c2prop:c1prop));
      if (c2prop > c1prop)
        fprintf(ofl,"%s%.2f",rpars->listdelmtr,c2prop/c1prop);
      else
        fprintf(ofl,"%s%.2f",rpars->listdelmtr,c1prop/c2prop);
      }
    }
  fputc('\n',ofl);
  }
}

double dm_rawmethunmethratio(DM_METH_CNTS *cp)
  /* return the proportion of meth to unmeth counts, -1 if
unmeth is zero */
{
if ((cp == NULL) || (cp->unmethcnt <= 0))
  return(-1.0);
else
  return((double) cp->methcnt/cp->unmethcnt);
}

double dm_methproportn(DM_METH_CNTS *cp)
  /* return proportionof meth to total counts, -1 if
total is zero */
{
int tcnts;

if ((cp == NULL) || ((tcnts = cp->unmethcnt + cp->methcnt) <= 0))
  return(-1.0);
else
  return(((double) cp->methcnt/tcnts));
}

void dm_prtusedcnts(FILE *ofl,
                    DM_METH_CNTS *clst,
                    char *delmtr)
/* traverse clst, printing used counts, prefix with
element No. */
{
DM_METH_CNTS *cp;
int eno;

cp = clst;
eno = 1;
while (cp != NULL)
  {
  if (cp->cntused)
    fprintf(ofl,"%s%d:%d+/%d-",delmtr,eno,cp->methcnt,cp->unmethcnt);
  cp = cp->nxtelt;
  eno++;
  }
}

double dm_minpr4chiorpair(DM_RUNPARS *rpp,
                          DM_CPG_BIN *binp,
                          int *dfcnt)
/* apply chisq/FE tests for list of binp, returning
the minimum probability found for pairwise, else the
chisq Pr.  Return 1.0 for failure.
if dfcnt is non-null, return either the df or the
count of pairwise FE tests to it */
{
DM_METH_CNTS *c1p;
DM_METH_CNTS *c2p;
double minpr;
double chisq;
int df;

if ((dm_minlstcontingency(binp->cntlst) >= 5.0) &&
      ((chisq = dm_chisq4cntlst(rpp,binp,binp->cntlst,&df)) > 0.0) &&
      (df > 0))
  {
  if (dfcnt != NULL)
    *dfcnt = df;
  return(1.0 - chi_sq_pr(df,chisq));
  }
else
  {
  minpr = 1.0;
  c1p = binp->cntlst;
  df = 0;
  while (c1p != NULL)
    {
    c2p = c1p->nxtelt;
    while (c2p != NULL)
      {
      if (dm_chkfemethcriteria(rpp,binp,c1p,c2p))
        {
        df++;
        if ((chisq = fishers_exact_t1(c1p->methcnt,c1p->unmethcnt,
                                        c2p->methcnt,c2p->unmethcnt)) < minpr)
          minpr = chisq;
        }
      c2p = c2p->nxtelt;
      }
    c1p = c1p->nxtelt;
    }
  if (dfcnt != NULL)
    *dfcnt = df;
  return(minpr);
  }
}

double dm_minpr4validfe(DM_RUNPARS *rpp,
                        DM_CPG_BIN *binp,
                        int *dfcnt)
/* apply pairwise FE tests for list of binp, returning
the minimum probability found.  Return 1.0 for failure.
if dfcnt is non-null, return the count of pairwise FE tests to it */
{
DM_METH_CNTS *c1p;
DM_METH_CNTS *c2p;
double minpr;
int df;
double fepr;

minpr = 1.0;
c1p = binp->cntlst;
df = 0;
while (c1p != NULL)
  {
  c2p = c1p->nxtelt;
  while (c2p != NULL)
    {
    if (dm_chkfemethcriteria(rpp,binp,c1p,c2p))
      {
      df++;
      if ((fepr = dm_checkprval(fishers_exact_t1(c1p->methcnt,c1p->unmethcnt,
                                                   c2p->methcnt,c2p->unmethcnt))) < minpr)
        minpr = fepr;
      }
    c2p = c2p->nxtelt;
    }
  c1p = c1p->nxtelt;
  }
if (dfcnt != NULL)
  *dfcnt = df;
return(minpr);
}

double dm_chi4validsamples(DM_RUNPARS *rpp,
                           DM_CPG_BIN *binp,
                           double *chisq,
                           int *df)
/* check all samples for binp for valid contingency, rejecting
those that don't, return chisq pr for remaining valid samples.
df returns the degrees of freedom: 0 if didn't work */
{
int methtot;
int unmethtot;
float methprop;
DM_METH_CNTS *mcp;
int vcnt;
int prvcnt;
int vhitssum;
float mexp;
float thistot;
float unmexp;
double chisum;
int dfcnt;
float ominuse;

mcp = binp->cntlst;
while (mcp != NULL)  /* count all, set cntused for all */
  {
  mcp->cntused = dm_chkcnts4mask(rpp,(rpp->binchkmask&(~1<<dm_chk_indcnt)),binp,mcp);
               /* avoid pr-based check calculation here */
  mcp = mcp->nxtelt;
  }
vcnt = dm_sumusedmethnunmeth(binp->cntlst,&methtot,&unmethtot);
prvcnt = 0;
chisum = 0.0;
dfcnt = 0;
while ((vcnt > 1) && (vcnt != prvcnt))
  {
  prvcnt = vcnt;
  if ((vhitssum = methtot + unmethtot) > 0)
    {
    methprop = (float) methtot/vhitssum;
    mcp = binp->cntlst;
    chisum = 0.0;
    dfcnt = 0;
    while (mcp != NULL)
      {
      if (mcp->cntused)
        {
        thistot = (float) (mcp->methcnt + mcp->unmethcnt);
        mexp = methprop * thistot;
        unmexp = thistot - mexp;
        if ((mexp < 5.0) || (unmexp < 5.0) || !dm_chkpercpgcnts(rpp,binp,mcp))
          mcp->cntused = 0;
        else
          {
          ominuse = (float) mcp->methcnt - mexp;
          chisum += ominuse*ominuse/mexp;
          ominuse = (float) mcp->unmethcnt - unmexp;
          chisum += ominuse*ominuse/unmexp;
          dfcnt++;
          }
        }
      mcp = mcp->nxtelt;
      }
    vcnt = dm_sumusedmethnunmeth(binp->cntlst,&methtot,&unmethtot);
    }
  else
    vcnt = 0;
  }
if (df != NULL)
  *df = dfcnt - 1;
if (chisq != NULL)
  *chisq = chisum;
if (dfcnt > 1)
  return(1.0 - chi_sq_pr(dfcnt-1,chisum));
else
  return(1.0);
}

double dm_chipr4validsamples(DM_RUNPARS *rpp,
                             DM_CPG_BIN *binp,
                             int *df)
/* check all samples for binp for valid contingency, rejecting
those that don't, return chisq pr for remaining valid samples.
df returns the degrees of freedom: 0 if didn't work */
{
return(dm_chi4validsamples(rpp,binp,NULL,df));
}

int dm_binmeetsthresholds(DM_RUNPARS *rpp,
                          DM_CPG_BIN *binp)
/* run thru the various threshold checks for bin binp.
(rpp->binprfun)() returns the minimum pr for binx, depending
on required method */
{
float hitspercpg;
float pr;
int df;
int tothits;

if ((binp->cpgcnt <= 0) && (rpp->binwidth > 0))
/* WGBS omits cpg counts, so return true */
  return(1);
else
  if ((binp->cpgcnt < rpp->mincpgs) || (binp->cpgcnt <= 0))
    return(0);
  else
    {
    tothits = dm_totbincnts(binp->cntlst,1);
    hitspercpg = (float) tothits/binp->cpgcnt;
    if (rpp->binprfun != NULL)
      {
      pr = dm_checkprval((*rpp->binprfun)(rpp,binp,&df));
      return(((rpp->cntthreshold <= 0) || (tothits >= rpp->cntthreshold)) &&
  /*           ((rpp->minhitspercpg == 0.0) || (hitspercpg >= rpp->minhitspercpg)) &&
             ((rpp->maxhitspercpg == 0.0) || (hitspercpg < rpp->maxhitspercpg)) && */
               ((rpp->prthreshold == 1.0) || (pr <= rpp->prthreshold)) &&
               ((df+1) >= rpp->minsamples));
      }
    else
      return((rpp->cntthreshold <= 0) || (tothits >= rpp->cntthreshold));
    }
}

double dm_totalss4bin(DM_RUNPARS *rpp,
                      DM_CPG_BIN *binp,
                      int *tdf)
/* return total ss and (if non-NULL) total
df in tdf.  return -1.0 if can't compute */
{
double sx;
double ssx;
int n;
DM_METH_CNTS *cntp;
double xmeth;

n = 0;
sx = ssx = 0.0;
cntp = binp->cntlst;
while (cntp != NULL)
  {
  if (dm_chkcnts4mask(rpp,rpp->binchkmask,binp,cntp) &&((cntp->methcnt > 0) || (cntp->unmethcnt > 0)))
    {
    xmeth = (double) cntp->methcnt/(cntp->methcnt+cntp->unmethcnt);
    n++;
    sx += xmeth;
    ssx += xmeth*xmeth;
    }
  cntp = cntp->nxtelt;
  }
if (tdf != NULL)
  *tdf = n - 1;
if (n > 0)
  return(ssx - sx*sx/n);
else
  return(-1.0);

}

void dm_clear_anova_dat(DM_ANOVA_DAT *adp)
  /* zero the values of adp */
{
int i;

adp->grandtotal = 0.0;
for (i=0; i< adp->maxgroup; i++)
  {
  *(adp->sx + i) = *(adp->ssx + i) = 0.0;
  *(adp->grpcnts + i) = 0;
  }
}

double dm_treatss4bin(DM_RUNPARS *rpp,
                      DM_CPG_BIN *binp,
                      DM_ANOVA_DAT *adp,
                      int *treatdf)
/* examine counts and, depending on group, accumulate
ss for treatment, if treatdf is non-NULL, then assign
treatment df to it */
{
DM_METH_CNTS *cntp;
double xmeth;
int gp;
int tdf;
double treatss;
int n;
double retval;

dm_clear_anova_dat(adp);
n = 0;
cntp = binp->cntlst;
while (cntp != NULL)
  {
  if (dm_chkcnts4mask(rpp,rpp->binchkmask,binp,cntp) &&
       ((cntp->methcnt > 0) || (cntp->unmethcnt > 0)))
    {
    xmeth = (double) cntp->methcnt/(cntp->methcnt+cntp->unmethcnt);
/* fprintf(stdout,"%d:%d:%.4f\t",n+1,cntp->smplgroup,xmeth); */
    *(adp->sx+cntp->smplgroup-1) += xmeth;
    *(adp->ssx+cntp->smplgroup-1) += xmeth*xmeth;
    *(adp->grpcnts+cntp->smplgroup-1) += 1;
    n++;
    adp->grandtotal += xmeth;
    }
  cntp = cntp->nxtelt;
  }
tdf = 0;
for (gp = 0; gp < adp->maxgroup; gp++)
  if (*(adp->grpcnts+gp) > 0)
    tdf++;
if (treatdf != NULL)
  *treatdf = tdf - 1;
if (tdf > 1)
  {
  treatss = 0.0;
  for (gp = 0; gp < adp->maxgroup; gp++)
    if (*(adp->grpcnts+gp) > 0)
      treatss += *(adp->sx+gp)**(adp->sx+gp)/(*(adp->grpcnts+gp));
  retval = treatss - adp->grandtotal*adp->grandtotal/n;
  }
else
  retval = -1.0;
return(retval);
}

int dm_groupcnt4bin(DM_RUNPARS *rpp,
                    DM_CPG_BIN *binp,
                    DM_ANOVA_DAT *adp)
  /* scan methcnts in binp and return number of different
groups which can be used in anova.  */
{
DM_METH_CNTS *cntp;
int garrp;
int validcnt;

dm_clear_anova_dat(adp);
cntp = binp->cntlst;
while (cntp != NULL)
  {
  if (dm_chkcnts4mask(rpp,rpp->binchkmask,binp,cntp) &&
        ((cntp->methcnt > 0) || (cntp->unmethcnt > 0)))
    *(adp->grpcnts+cntp->smplgroup-1) += 1;
  cntp = cntp->nxtelt;
  }
validcnt = 0;
for (garrp = 0; garrp < adp->maxgroup; garrp++)
  if (*(adp->grpcnts+garrp) > 1)
    validcnt++;
return(validcnt);
}

#define CF_ITERLIMIT 100
#ifdef DBL_EPSILON
#define CF_PRECLIMIT DBL_EPSILON
#else
#define CF_PRECLIMIT 3.0e-7
#endif

double dm_greaterdbl(double x,
                     double y)
/* return the greater of x,y */
{
if (x < y)
  return(y);
else
  return(x);
}

double dm_betainccontfract(double a,
                           double b,
                           double x)
/* use a modified Lentz method to derive
the incomplete beta function (regularised) by
continued fraction. Modified from Numerical Recipes
in C.  Return -1.0 for error states */
{
int i;
int ix2;
double tmp;
double c;
double d;
double del;
double h;
double aplusb;
double am1;
double qp1;
double plimit;

plimit = CF_PRECLIMIT;
aplusb = a + b;
qp1 = a + 1.0;
am1 = a - 1.0;
c = 1.0;
d = dm_greaterdbl((1.0 - aplusb*x/qp1),DBL_MIN);
h = d = 1.0/d;
for (i = 1; i <= CF_ITERLIMIT; i++)
  {
  ix2 = i * 2;
  tmp = i*(b-i)*x/((am1+ix2)*(a+ix2));
  d = dm_greaterdbl((1.0 + tmp*d),DBL_MIN);
  c = dm_greaterdbl((1.0 + tmp/c),DBL_MIN);
  d = 1.0/d;
  h *= d*c;
  tmp = -(a + i)*(aplusb + i)*x/((a + ix2)*(qp1 + ix2));
  d = dm_greaterdbl((1.0 + tmp * d),DBL_MIN);
  c = dm_greaterdbl((1.0 + tmp/c),DBL_MIN);
  d = 1.0/d;
  del = d*c;
  h *= del;
  if (fabs(del-1.0) < CF_PRECLIMIT)
    break;
  }
if (i > CF_ITERLIMIT)
  return(-1.0);
else
  return(h);
}

double dm_incomplbeta(double a,
                      double b,
                      double x)
/* Returns the incomplete beta function Ix(a,b).
return -1.0 for erroneous values */
{
double betincval;
double bx;

if ((x < 0.0) || (x > 1.0))
  return(-1.0);
else
  {
  if ((x == 0.0) || (x == 1.0))
    bx = 0.0;
  else
    bx = exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
    {
    betincval = dm_betainccontfract(a,b,x)/a;
    if (betincval >= 0.0)
      return(bx*betincval);
    else
      return(-1.0);
    }
  else
    {
    betincval = dm_betainccontfract(b,a,(1.0 - x))/b;
    if (betincval >= 0.0)
      return(1.0 - bx*betincval);
    else
      return(-1.0);
    }
  }
}

double dm_pr_f_dist(double x,
                    int df1,
                    int df2)
/* Return pr(F=x,df1,df2). upper tail.
returns -1.0 for erroneous data */
{
return(dm_incomplbeta(df2/2.0,df1/2.0,(df2/(df1*x+df2))));
}

float dm_grpmeth4binchk(DM_CPG_BIN *binp,
                        int *totlcnts,
                        int grp,
                        DM_RUNPARS *rpp,
                        int (*isvalid)(DM_RUNPARS *xrpp,
                                       DM_CPG_BIN *xbinp,
                                       DM_METH_CNTS *xcntp))
/* scan binp for entries matching grp and
return meth proportion for it.  If totlcnts is nonNULL
then return total to it.  Return -1.0 if
can't calculate proportion. Function (*isvalid) checks
if counts meet criteria */
{
int gmeth;
int gtot;
DM_METH_CNTS *blp;

gmeth = gtot = 0;
blp = binp->cntlst;
while (blp != NULL)
  {
  if ((blp->smplgroup == grp) && (*isvalid)(rpp,binp,blp))
    {
    gmeth += blp->methcnt;
    gtot += blp->methcnt + blp->unmethcnt;
    }
  blp = blp->nxtelt;
  }
if (totlcnts != NULL)
  *totlcnts = gtot;
if (gtot > 0)
  return((float)gmeth/(float)gtot);
else
  return(-1.0);
}

float dm_ratiogrpmeths4abin(DM_CPG_BIN *binp,
                            int *totlcnts)
/* work out the ratio of methylation proportion
for each of the two groups in binp - return
-1.0 if can't calculate.  If totlcnts nonNULL
then return the total counts to it */
{
int g1meth;
int g1tot;
int g2meth;
int g2tot;
DM_METH_CNTS *blp;

g1meth = g1tot = g2meth = g2tot = 0;
blp = binp->cntlst;
while (blp != NULL)
  {
  if (blp->smplgroup)
    {
    g2meth += blp->methcnt;
    g2tot += blp->methcnt + blp->unmethcnt;
    }
  else
    {
    g1meth += blp->methcnt;
    g1tot += blp->methcnt + blp->unmethcnt;
    }
  blp = blp->nxtelt;
  }
if (totlcnts != NULL)
  *totlcnts = g1tot + g2tot;
if ((g1tot > 0) && (g2meth > 0))
  return(((float)g1meth*(float)g2tot)/((float)g1tot*(float)g2meth));
else
  return(-1.0);
}

char *dm_meth2grpstr(DM_RUNPARS *rpp,
                     int mxgrp,
                     float meth)
/* return string for mratio and mxgrp */
{
if (meth == -1.0)
  return("??");
else
  if (*(rpp->groupidlist+mxgrp) != NULL)
    return(*(rpp->groupidlist+mxgrp));
  else
    return("<null>");
}

double dm_anovapr4bin(DM_RUNPARS *rparsp,
                      DM_CPG_BIN *binp,
                      int *xtotldf,
                      int *xtreatdf,
                      double *xf_val)
/* perform Anova on this bin. Return probability.  If
xf_val, xtotldf & xtreatdf are non-NULL, then return those
parameters to them.  Return -1.0 for invalid, when the other
passed parameters will be undefined. */
{
int totldf;
int treatdf;
double totlss;
double treatss;
double f_val;

if (dm_groupcnt4bin(rparsp,binp,rparsp->anova_dat) > 1)
  {
  totlss = dm_totalss4bin(rparsp,binp,&totldf);
  treatss = dm_treatss4bin(rparsp,binp,rparsp->anova_dat,&treatdf);
  if ((totlss > 0.0) && (treatss > 0.0))
    {
    f_val = (treatss/treatdf)*((totldf - treatdf)/(totlss - treatss));
    if (xtotldf != NULL)
      *xtotldf = totldf;
    if (xtreatdf != NULL)
      *xtreatdf = treatdf;
    if (xf_val != NULL)
      *xf_val = f_val;
    return(dm_checkprval(dm_pr_f_dist(f_val,treatdf,(totldf-treatdf))));
    }
  }
return(-1.0);
}

double dm_anovapr4validsampls(DM_RUNPARS *rpp,
                              DM_CPG_BIN *binp,
                              int *totldf)
/* return the Pr and optionally totldf for
anova on binp. the Total df are returned to
a non-NULL totldf */
{
return(dm_anovapr4bin(rpp,binp,totldf,NULL,NULL));
}

char *dm_int2letterstring(int cnt)
  /* return a character string to conform
to cnt, "a","b","c",.."z","aa","ab".. for
1,2,3... */
{
int chrcnt;
int cpy;
char *sp;

cpy = cnt;
chrcnt = 0;
while (cpy > 0)
  {
  chrcnt++;
  cpy /= 26;
  }
sp = &istring[chrcnt];
*sp = '\0';
sp--;
cpy = cnt;
while (cpy > 0)
  {
  *sp = 'a' + cpy%26 - 1;
  cpy /= 26;
  sp--;
  }
return(&istring[0]);
}

void dm_diffcntsscan4abin(DM_OUTMODE omod,
                          DM_RUNPARS *rparsp,
                          FILE *ofl,
                          int chrno,
                          DM_CPG_BIN *binp)
/* scan binp and do pairwise diff meth cnts */
{
DM_METH_CNTS *c1p;
DM_METH_CNTS *c2p;
int ocnt;
int c1no;
int c2no;
double chisq;
double chisqpr;
double *prvals;
int nprvals;
double minpr;
int cp;
int totldf;
int treatdf;
double f_val;
double fpr;
double sdev;
double mean;
int scnt;
int gtrmethgrp;
float maxmeth;
float xmeth;
int comma;
char *gmethstr;
int smpl;

switch (omod)
  {
  case dm_out_bestpairchi:
  case dm_out_allpairchi:
    if (dm_binmeetsthresholds(rparsp,binp))
      {
/* try to run with best of FE or chi, giving min pairwise PR or chiPR */
/* first try chi */
      chisq = -1.0;
      if (dm_minlstcontingency(binp->cntlst) >= 5.0)
        chisq = dm_chisq4cntlst(rparsp,binp,binp->cntlst,&c1no);
      if ((chisq > 0.0) && (c1no > 0))
        {
        chisqpr = 1.0 - dm_checkprval(chi_sq_pr(c1no,chisq));
        dm_startoputline(ofl,rparsp,binp,chrno);
        fprintf(ofl,"%s%4g%sChi_%.4f_%ddf",rparsp->listdelmtr,dm_checkprval(chisqpr),
                  rparsp->listdelmtr,chisq,c1no);
        if (rparsp->displycnts)
          dm_prtusedcnts(ofl,binp->cntlst,rparsp->listdelmtr);
        fputc('\n',ofl);
        }
      else
        chisq = -1.0;
      if (chisq < 0.0)   /* need to do the pairwise thing */
        {
        c1p = binp->cntlst;
        nprvals = 0;
        while (c1p != NULL)
          {
          c2p = c1p->nxtelt;
          while (c2p != NULL)
            {
            if (dm_chkfemethcriteria(rparsp,binp,c1p,c2p))
              nprvals++;
            c2p = c2p->nxtelt;
            }
          c1p = c1p->nxtelt;
          }
        if (nprvals > 0)
          {
          prvals = (double *)getmemory(sizeof(double)*nprvals,"Pr list");
          c1p = binp->cntlst;
          nprvals = 0;
          while (c1p != NULL)
            {
            c2p = c1p->nxtelt;
            while (c2p != NULL)
              {
              if (dm_chkfemethcriteria(rparsp,binp,c1p,c2p))
                {
                *(prvals+nprvals) = dm_checkprval(fishers_exact_t1(c1p->methcnt,c1p->unmethcnt,
                                                                     c2p->methcnt,c2p->unmethcnt));
                c1p->cntused = c2p->cntused = 1;
                nprvals++;
                }
              c2p = c2p->nxtelt;
              }
            c1p = c1p->nxtelt;
            }
          minpr = dm_ascendsort(prvals,nprvals);
          dm_startoputline(ofl,rparsp,binp,chrno);
          fprintf(ofl,"%s%4g%sFE_%d",rparsp->listdelmtr,minpr,rparsp->listdelmtr,
                    nprvals);
          if (rparsp->omode == dm_out_allpairchi)
            for (ocnt = 1; ocnt < nprvals; ocnt++)
              fprintf(ofl,"%s%g",rparsp->listdelmtr,*(prvals+ocnt));
          if (rparsp->displycnts)
            dm_prtusedcnts(ofl,binp->cntlst,rparsp->listdelmtr);
          fputc('\n',ofl);
          memfree(prvals);
          }
        }
      }
    break;
  case dm_out_chisq:
    if (dm_binmeetsthresholds(rparsp,binp))
      {
      if (dm_minlstcontingency(binp->cntlst) >= 5.0)
        {
        chisq = dm_chisq4cntlst(rparsp,binp,binp->cntlst,&c1no);
        if ((c1no > 0) && (chisq > 0.0))
          {
          chisqpr = 1.0 - dm_checkprval(chi_sq_pr(c1no,chisq));
          dm_startoputline(ofl,rparsp,binp,chrno);
          fprintf(ofl,"%s%4g%sChi_%.4f_%ddf\n",rparsp->listdelmtr,chisqpr,
                    rparsp->listdelmtr,chisq,c1no);
          }
        }
      else
        dm_diffcntsscan4abin(dm_out_fisherforce,rparsp,ofl,chrno,binp);
      }
    break;
  case dm_out_chiforce:
    if (dm_binmeetsthresholds(rparsp,binp))
      {
      chisqpr = dm_checkprval(dm_chi4validsamples(rparsp,binp,&chisq,&c1no));
      if ((c1no > 0) && (chisq > 0.0))
        {
        dm_startoputline(ofl,rparsp,binp,chrno);
        fprintf(ofl,"%s%4g%sChi_%.4f_%ddf\n",rparsp->listdelmtr,chisqpr,
                  rparsp->listdelmtr,chisq,c1no);
        }
      }
    break;
  case dm_out_cpg_nzdetl:
    if (dm_totbincnts(binp->cntlst,1))
      dm_diffcntsscan4abin(dm_out_cpg_detail,rparsp,ofl,chrno,binp);
    break;
  case dm_out_cpg_detail:
    if (dm_binmeetsthresholds(rparsp,binp))
      {
      dm_startoputline(ofl,rparsp,binp,chrno);
      fputc('\n',ofl);
      c1p = binp->cntlst;
      c1no = 1;
      while (c1p != NULL)
        {
        if (dm_cntsarevalid(rparsp,binp,c1p)) /* (c1p->methcnt > 0) || (c1p->unmethcnt > 0)) */
          {
          cp = 1;
          fprintf(ofl,"Sample %d",c1no);
          if ((c1p->cpgcntarr->metcnt > 0) || (c1p->cpgcntarr->unmetcnt > 0))
            fprintf(ofl,"%sStray:%d+%d-",rparsp->listdelmtr,
                      c1p->cpgcntarr->metcnt,c1p->cpgcntarr->unmetcnt);
          while (cp <= binp->maxccnt)
            {
            fprintf(ofl,"%s%d:%d+%d-",rparsp->listdelmtr,(c1p->cpgcntarr+cp)->cpgpos,
                      (c1p->cpgcntarr+cp)->metcnt,(c1p->cpgcntarr+cp)->unmetcnt);
            cp++;
            }
          fputc('\n',ofl);
          }
        c1no++;
        c1p = c1p->nxtelt;
        }
      }
    break;
  case dm_out_sdev:
    if (dm_binmeetsthresholds(rparsp,binp))
      {
      c1p = binp->cntlst;
      sdev = dm_stddev4nonngtv(rparsp,binp,binp->cntlst,dm_methproportn,
                                 &mean,&scnt);
      if (scnt > 1)
        {
        dm_startoputline(ofl,rparsp,binp,chrno);
        fprintf(ofl,"%s%.4f%s%.4f%s%d\n",rparsp->listdelmtr,mean,rparsp->listdelmtr,
                  sdev,rparsp->listdelmtr,scnt);
        }
      }
    break;
  case dm_out_anova:
  case dm_out_anova_gtr:
  case dm_out_anova_more:
    if (((fpr = dm_anovapr4bin(rparsp,binp,&totldf,&treatdf,&f_val)) >= 0.0) &&
         ((rparsp->prthreshold == 1.0) || ((float)fpr <= rparsp->prthreshold)))
      {
      dm_startoputline(ofl,rparsp,binp,chrno);
      fprintf(ofl,"%s%0g%sF(%d,%d)=%.2f",rparsp->listdelmtr,fpr,
                rparsp->listdelmtr,treatdf,(totldf-treatdf),f_val);
      if ((rparsp->omode == dm_out_anova_gtr) || (rparsp->omode == dm_out_anova_more))
        {
        maxmeth = -1.0;
        gtrmethgrp = -1;
        for (cp = 0; cp < rparsp->anova_dat->maxgroup; cp++)
          if ((xmeth = dm_grpmeth4binchk(binp,NULL,cp+1,rparsp,dm_cntsarevalid)) > maxmeth)
            {
            maxmeth = xmeth;
            gtrmethgrp = cp;
            }
        gmethstr = dm_meth2grpstr(rparsp,gtrmethgrp,maxmeth);
        fprintf(ofl,"%s%s%s",rparsp->listdelmtr,
                  gmethstr,
                  rparsp->listdelmtr);
        for (cp = 0; cp < rparsp->anova_dat->maxgroup; cp++)
          fprintf(ofl,"%s%s=%d",(cp>0?",":""),
                    ((*(rparsp->groupidlist+cp)!=NULL)?*(rparsp->groupidlist+cp):"<null>"),
                    *(rparsp->anova_dat->grpcnts+cp));
        fputs(rparsp->listdelmtr,ofl);
        comma = 0;
        for (cp=0; cp < rparsp->anova_dat->maxgroup; cp++)
          if ((xmeth = dm_grpmeth4binchk(binp,NULL,cp+1,rparsp,dm_cntsarevalid)) >= 0.0)
            {
            fprintf(ofl,"%s%s=%.4f",(comma?",":""),
                      ((*(rparsp->groupidlist+cp)!=NULL)?*(rparsp->groupidlist+cp):"<null>"),
                      xmeth);
            comma++;
            }
        if (rparsp->omode == dm_out_anova_more)
          {
          for (cp = 0; cp <= rparsp->anova_dat->maxgroup;cp++)
            {
            if ((cp < rparsp->anova_dat->maxgroup) &&
                  (dm_grpmeth4binchk(binp,NULL,cp+1,rparsp,dm_cntsarevalid) >= 0.0))
              fputc(';',ofl);
            c1p = binp->cntlst;
            smpl = 1;
            comma = 0;
            while (c1p != NULL)
              {
              if (c1p->smplgroup == cp + 1)
                {
                if (dm_chkcnts4mask(rparsp,rparsp->binchkmask,binp,c1p) &&
                      ((c1p->methcnt > 0) || (c1p->unmethcnt > 0)))
                  {
                  fprintf(ofl,"%s%s%s=%.4f",(comma?",":""),
                            ((*(rparsp->groupidlist+cp)!=NULL)?*(rparsp->groupidlist+cp):"<null>"),
                            dm_int2letterstring(smpl),
                            (float)c1p->methcnt/((float) c1p->methcnt+(float)c1p->unmethcnt));
                  comma++;
                  }
                smpl++;
                }
              c1p = c1p->nxtelt;
              }
            }
          }
        }
      fputc('\n',ofl);
      }
    break;
  case dm_out_pairwise:
  case dm_out_fisherforce:
  default:
    if (dm_binmeetsthresholds(rparsp,binp))
      {
      if ((dm_minlstcontingency(binp->cntlst) >= 5.0) && /* should use chi */
           (omod != dm_out_fisherforce))
        dm_diffcntsscan4abin(dm_out_chiforce,rparsp,ofl,chrno,binp);
      else
        {      
        c1p = binp->cntlst;
        c1no = 1;
        while (c1p != NULL)
          {
          c2p = c1p->nxtelt;
          c2no = c1no + 1;
          while (c2p != NULL)
            {
            dm_scandiffmethcnts(ofl,rparsp,chrno,binp,c1p,c2p,c1no,c2no);
            c2p = c2p->nxtelt;
            c2no++;
            }
          c1p = c1p->nxtelt;
          c1no++;
          }
        }
      }
    break;
  }
}

void dm_scanbins4achromosome(FILE *ofl,
                             DM_RUNPARS *rparsp,
                             RBC_CHRNO chrno)
/* don't check chrno: assume valid;  go through all bins
in chromosome chrno and check ... */
{
DM_CPG_BIN *bp;

bp = (rparsp->chrbininfo+chrno)->binlst;
while (bp != NULL)
  {
  dm_diffcntsscan4abin(rparsp->omode,rparsp,ofl,chrno,bp);
  bp = bp->nxtbin;
  }
}

void dm_scanallchrombins(FILE *ofl,
                         DM_RUNPARS *rpars)
/* go thru all chromosomes, scan bins */
{
RBC_CHRNO chrno;

for (chrno = 0; chrno < rpars->maxchrno; chrno++)
  dm_scanbins4achromosome(ofl,rpars,chrno);
}

int bc_amlgmatelocpgbins(DM_CHR_INFO *cbinlst,
                         RBC_CHRNO maxchr,
                         int minfrag,
                         int maxfrag,
                         int cpgmin)
/* scan cbinlst for bins which may fail a CpG count
but which could be amalgamated with adjacent bins to
pass the CpG count criterion implied by cpgmin.
Amalgamate valid bins.  Note that this may leave bins
that exceed length criteriam but if these are
made by joining valid adjacent bins, then the
size selection during restricted representation library
preparation should have left valid reads for those
regions.  Return number of bins pruned or amalgamated */
{
int abincnt;
RBC_CHRNO chrno;
DM_CPG_BIN *bp;
DM_CPG_BIN *nxt;
DM_METH_CNTS *cntp;

abincnt = 0;
for (chrno = Chr1; chrno <= maxchr; chrno++)
  {
  abincnt += bc_prunebins4len(&(cbinlst+chrno-1)->binlst,minfrag,maxfrag);
  bp = (cbinlst+chrno-1)->binlst;
  while (bp != NULL)
    {
    nxt = bp->nxtbin;
    if (bc_adjacentbins(bp,nxt) &&
          (bc_binundercpglmt(bp,cpgmin) || bc_binundercpglmt(nxt,cpgmin))
             && bc_cpglmttest(cpgmin,(bp->cpgcnt + nxt->cpgcnt)))
      {  /* checks that bins are adjacent, one or other fails cpg criterion & together they pass it */
/*          if (debuglevel > RBC_dbg_none)
            fprintf(stdout,"C%s amalgamating %d..%d (%d) & %d..%d (%d) with %d+%dCpGs\n",
                      rbc_chrno2str(chrno,1),bp->spos,(bp->spos+bp->binlen-1),bp->binlen,
                      nxt->spos,(nxt->spos+nxt->binlen-1),nxt->binlen,bp->cpgcnt,nxt->cpgcnt); */
      bp->binlen += nxt->binlen;
      bp->cpgcnt += nxt->cpgcnt;
      cntp = nxt->cntlst;
      while (cntp != NULL)
        {
        (void) dm_appndcntelt(&bp->cntlst,cntp->methcnt,cntp->unmethcnt,cntp->smplgroup);
        cntp = cntp->nxtelt;
        }
      rbc_delcntbin(nxt,&(cbinlst+chrno-1)->binlst);
      abincnt++;
      }
    bp = bp->nxtbin;
    }
/* now prune failed cpg bins */
  abincnt += bc_prunebins4cpgmin(&(cbinlst+chrno-1)->binlst,cpgmin);
  }
return(abincnt);
}

int dm_jointestsucceedalways(DM_RUNPARS *rpp,
                             DM_CPG_BIN *binp1,
                             DM_CPG_BIN *binp2)
{
return(1);
}

int dm_joinadjacntfrags(DM_RUNPARS *rpp,
                        DM_CHR_INFO *cbinlst,
                        int (*jointest)(DM_RUNPARS *rpp,
                                        DM_CPG_BIN *binp1,
                                        DM_CPG_BIN *binp2))
  /* scan cbinlst for adjacent RRGS fragment bins.
Join them if they satisfy jointest(), returning
count of number so treated. */
{
int abincnt;
RBC_CHRNO chrno;
DM_CPG_BIN *bp;
DM_CPG_BIN *nxt;
DM_METH_CNTS *cntp;

abincnt = 0;
for (chrno = Chr1; chrno <= rpp->maxchrno; chrno++)
  {
  bp = (cbinlst+chrno-1)->binlst;
  while (bp != NULL)
    {
    nxt = bp->nxtbin;
    if (bc_adjacentbins(bp,nxt) && (*jointest)(rpp,bp,nxt))
      {
/*          if (debuglevel > RBC_dbg_none)
            fprintf(stdout,"C%s amalgamating %d..%d (%d) & %d..%d (%d) with %d+%dCpGs\n",
                      rbc_chrno2str(chrno,1),bp->spos,(bp->spos+bp->binlen-1),bp->binlen,
                      nxt->spos,(nxt->spos+nxt->binlen-1),nxt->binlen,bp->cpgcnt,nxt->cpgcnt); */
      bp->binlen += nxt->binlen;
      bp->cpgcnt += nxt->cpgcnt;
      cntp = nxt->cntlst;
      while (cntp != NULL)
        {
        (void) dm_appndcntelt(&bp->cntlst,cntp->methcnt,cntp->unmethcnt,cntp->smplgroup);
        cntp = cntp->nxtelt;
        }
      rbc_delcntbin(nxt,&(cbinlst+chrno-1)->binlst);
      abincnt++;
      }
    else
      bp = bp->nxtbin;
    }
  }
return(abincnt);
}

void dm_listoutdata(FILE *ofl,
                    DM_RUNPARS *rpp)
/* list bins to ofl */
{
RBC_CHRNO chrno;
DM_CPG_BIN *bp;
DM_METH_CNTS *clp;
int bcnt;
int vbincnt;

dm_commonhdr(ofl,rpp);
fputc('\n',ofl);
for (chrno = 0; chrno < rpp->maxchrno; chrno++)
  {
  bp = (rpp->chrbininfo+chrno)->binlst;
  while (bp != NULL)
    {
    if (rpp->omode == dm_out_listnz)
      {
      clp = bp->cntlst;
      vbincnt = bcnt = 0;
      while (clp != NULL)
        {
        bcnt += clp->methcnt + clp->unmethcnt;
        if (dm_cntsarevalid(rpp,bp,clp))
          vbincnt++;
        clp = clp->nxtelt;
        }
      }
    if ((rpp->omode == dm_out_list) ||
           ((rpp->omode == dm_out_listnz) && (bcnt > 0) && (vbincnt > 0)))
      {
      fprintf(ofl,"%s%s%d%s%d%s%d%s%d",
               wlu_intval2str(rpp->chridlu,chrno),
               rpp->listdelmtr,bp->spos,rpp->listdelmtr,(bp->spos+bp->binlen-1),
               rpp->listdelmtr,bp->binlen,rpp->listdelmtr,bp->cpgcnt);
      clp = bp->cntlst;
      while (clp != NULL)
        {
        if (dm_cntsarevalid(rpp,bp,clp))
          fprintf(ofl,"%s%d+/%d-",rpp->listdelmtr,clp->methcnt,clp->unmethcnt);
        else
          fprintf(ofl,"%s-",rpp->listdelmtr);
        clp = clp->nxtelt;
        }
      fputc('\n',ofl);
      }
    bp = bp->nxtbin;
    }
  }
}

int dm_scanchrinfofl(FILE *cifile,
                     DM_CHR_INFO *cinfo,
                     WRD_LUSTRCT *chridlu)
/* scan lines in cifile as chromosomeID Filenam.
return number.  If cinfo is non-NULL, then put
values into fields. '#' causes remainder of line
to be treated as comment (unless quoted). */
{
int lno;
char *cid;
char *chrfile;
char nc;
char *bp;
int onquote;
char *target;
char prvchr;
int incmt;

cid = (char *) getmemory(MAXPATHLEN + 1,"chrid");
chrfile = (char *) getmemory(MAXPATHLEN + 1,"chrflnam");
target = bp = cid;
*chrfile = *bp = '\0';
lno = 0;
incmt = onquote = 0;
prvchr = '\0';
while ((nc = fgetc(cifile)) != EOF)
  switch (nc)
    {
    case '\n':
      if ((strlen(cid) > 0) && (strlen(chrfile) > 0))
        {
        if (cinfo != NULL)
          {
          (cinfo+lno)->chrid = bas_strdup(cid);
          (cinfo+lno)->chrflname = bas_strdup(chrfile);
          }
        if (chridlu != NULL)
          wlu_addwrd(chridlu,cid,(int) lno, NULL);
        *chrfile = *cid = '\0';
        lno++;
        }
/* reset to start */
      target = bp = cid;
      incmt = onquote = 0;
      prvchr = '\0';
      break;
    case '"':    /* start or end of quoted string */
      if (!incmt)
        onquote = !onquote;
      break;
    case ' ':     /* white space */
    case '\t':
      if (onquote)
        {
        if ((target == cid) && (!incmt))
          bas_appchr(cid,&bp,nc,MAXPATHLEN);
        else
          if ((target == chrfile) && (!incmt))
            bas_appchr(chrfile,&bp,nc,MAXPATHLEN);
        }
      else 
        {
        if (!incmt)
          {
          if (!isblank(prvchr))
            if (target == cid)
              target = bp = chrfile;
            else
              target = NULL;
          prvchr = nc;
          }
        }
      break;
    case '#':    /* start of comment field */
      if (!onquote)
        {
        incmt = 1;
        target = NULL;
        }
    default:
      if ((target == cid) && !incmt)
        bas_appchr(cid,&bp,nc,MAXPATHLEN);
      else
        if ((target == chrfile) && !incmt)
          bas_appchr(chrfile,&bp,nc,MAXPATHLEN);
      prvchr = nc;
      break;
    }
memfree(chrfile);
memfree(cid);
return(lno);
}      

void dm_bldchrfilenames(char *hdrstr,
                        DM_RUNPARS *rpp)
/* create set of file names based on hdrstr.
put into rpp->chrbininfo array (pre-existing) */
{
int cno;
char *sqfnam;
int nblen;

sqfnam = (char *) getmemory((nblen = (int) strlen(hdrstr) + 16),"Sq filename buff");
for (cno = 0; cno < rpp->maxchrno; cno++)
  {
  (rpp->chrbininfo+cno)->chrid = bas_strdup(any_chrno2str(rpp->maxchrno,rpp->havexy,(cno+1),1));
  wlu_addwrd(rpp->chridlu,(rpp->chrbininfo+cno)->chrid,(int) cno,NULL);
  snprintf(sqfnam,nblen-1,"%s%s.fa",hdrstr,(rpp->chrbininfo+cno)->chrid);
  (rpp->chrbininfo+cno)->chrflname = bas_strdup(sqfnam);
  }
memfree(sqfnam);
}

DM_SITES *dm_appndsite(DM_SITES **sitlst,
                       char *site,
                       int cutoffset)
/* dup site and append an element for it to the end of sitlst.
return address of new element in case useful */
{
DM_SITES *prev, *end_ptr;

if (sitlst != NULL)
  {
  prev = end_ptr = *sitlst;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtsite;
    }
  end_ptr = (DM_SITES *) getmemory(sizeof(DM_SITES),"Site_elt");
  end_ptr->nxtsite = NULL;
  end_ptr->sitestr = bas_strdup(site);
  end_ptr->offset = cutoffset;
  if (*sitlst == NULL)
    *sitlst = end_ptr;
  else
    prev->nxtsite = end_ptr;
  return(end_ptr);
  }
else
  return(NULL);
}

DM_SITES *dm_appndsitestr(DM_SITES **sitlst,
                          char *recogstr)
/* recogstr is a restriction recognition site
with (optional) cleavage position indicated with
'^'.  append the site/offset derived from this to
*sitlst */
{
char *scopy;
char *sp;
char *dp;
int offst;

sp = dp = scopy = bas_strdup(recogstr);
offst = 0;
while (*sp != '\0')
  {
  if (*sp == '^')
    offst = (int) (sp - scopy);
  else
    {
    *dp = *sp;
    dp++;
    }
  sp++;
  }
*dp = '\0';
return(dm_appndsite(sitlst,scopy,offst));
}

int dm_readrstsites(FILE *rfl,
                    DM_SITES **rstlst)
/* read successive lines from rfl, check
length > 0 and process as restriction sites.
return number of sites found */
{
char *lbuf;
char *lp;
int scnt;

lbuf = (char *) getmemory(FILENAME_MAX + 1,"site buff");
scnt = 0;
while ((lp = fgets(lbuf,FILENAME_MAX,rfl)) != NULL)
  if (strlen(lp) > 0)
    {
    while (*lp != '\0')
      {
      if (*lp == '\n')
        *lp = '\0';
      lp++;
      }
    (void) dm_appndsitestr(rstlst,lbuf);
    scnt++;
    }
memfree(lbuf);
return(scnt);
}

int dm_chk_range_sanity(DM_RUNPARS *rpp,
                        int chrno,
                        int rstart,
                        int rstop)
/* check rstart & rstop valid for chrno.  die if not, but
return 1 for OK */
{
if (rstart <= rstop)
  return(1);
else
  {
  err_msg_die("region start %d > stop %d for Chr %d\n",rstart,rstop,chrno);
  return(0);
  }
}

int dm_readbinfile(FILE *binfl,
                   DM_RUNPARS *rpp)
/* read lines from binfl expecting 'chr\tstart\tstop\n'
as a definition of bin positions for chr.  Return
number of such bins */
{
int binc;
char *line;
char *linecpy;
char *cpyfree;
char **tokns;
int tcnt;
int chrno;
int bstart;
int bstop;
char **lp;

binc = 0;
line = (char *) getmemory(FILENAME_MAX,"bin line");
tokns = (char **) getmemory(3*sizeof(char *),"bin tokenlist");
while (fgets(line,FILENAME_MAX-1,binfl) != NULL)
  {
  linecpy = cpyfree = bas_strdup(line);
  tcnt = 0;
  for (lp = tokns; ((*lp = strsep(&linecpy," \t\n")) != NULL);)
    {
    if (**lp != '\0')
      {
      tcnt++;
      if (++lp >= (tokns+3))
        break;
      }
    }
  if ((chrno = wlu_chkwrd(rpp->chridlu,*tokns)) < 0)
    err_msg_die("Invalid chromosome ID '%s'\n",*tokns);
  else
    {
    if ((rpp->uchrno == NULL) || (chrno == wlu_chkwrd(rpp->chridlu,rpp->uchrno)))
      {
      bstart = (int) strtol(*(tokns+1),NULL,10);
      bstop = (int) strtol(*(tokns+2),NULL,10);
      if (dm_chk_range_sanity(rpp,chrno,bstart,bstop))
        {
        (void) dm_appndbin(&(rpp->chrbininfo+chrno)->binlst,bstart,bstop-bstart+1,0,rpp);
        binc++;
        }
      }
    }
  memfree(cpyfree);
  }
return(binc);
}

int main(int argc,
         char *argv[])
{
int ap;
char op;
int ecnts;
FILE *srcfl;
int chrcnt;
RBC_CHRNO chrno;
char *fendp;
int cpgcnt;
FS_FSMSTRCT *cpgfsmp;    /* ptr to CpG searching fsm */
int minfrag;
int maxfrag;
int amlgmate;
FILE *dfile;
DM_FNAM_ELT *fnamlst;
DM_FNAM_ELT *fnamptr;
DM_RUNPARS rpars;
DM_SITES *sitp;
int groupno;
char *groupid;
int i;

bas_initmalfl("diffmeth_mllcdgb.txt");
debuglevel = RBC_dbg_none;
ecnts = 0;
srcfl = NULL;
cpgfsmp = fs_initnewfsm(4,1,FS_inv_ignor);
minfrag = maxfrag = cpgcnt = 0;
amlgmate = 0;
dfile = NULL;
fnamlst = NULL;
rpars.maxchrno = ChrY;
rpars.havexy = 1;
rpars.folddiff = 2.0;
rpars.cntthreshold = 1;
rpars.minhitspercpg = rpars.maxhitspercpg = 0.0;
rpars.omode = dm_out_pairwise;
rpars.displycnts = 0;
rpars.chrbininfo = NULL;
rpars.binclusterfreq = DM_LIST_FREQ_DEF;
rpars.listdelmtr = "\t";
rpars.chrprefix = "";
rpars.cntcpgsas1 = 1;
rpars.prthreshold = 1.0;
rpars.minsamples = 2;
rpars.mincpgs = 1;
rpars.minmethprop = 0.0;
rpars.joinfrags = 0;
rpars.needseq = 0;
rpars.sambuflen = DEF_SAMBUFLEN;
rpars.uchrno = NULL;
rpars.cachechrno = -1;
rpars.currbin = NULL;
rpars.samfst3p2prv = 0;
rpars.binwidth = 0;
rpars.binchkmask = 0;
rpars.cpgdetails = 0;
rpars.cpgcntscrit = 0;
rpars.binprfun = NULL;
rpars.anova_dat = NULL;
rpars.maxgrp = 0;
rpars.genomsqhdrstr = NULL;
rpars.chrinfoflnam = NULL;
rpars.chrinfofl = NULL;
rpars.srcmode = DM_src_cpglist;
rpars.compressed = 0;
#ifndef NO_ZLIB
rpars.bamrdp = NULL;
#endif
rpars.rstfile = NULL;
rpars.rstflnam = NULL;
rpars.sitelist = NULL;
rpars.allcs = 0;
rpars.groupidlist = NULL;
rpars.binstyle = DM_bin_restenz;
rpars.binfile = NULL;
for (ap = 1; ap < argc; ap++)
  {
/*  fprintf(stdout,"argv[%d]=\"%s\"\n",ap,argv[ap]); */
  if (*argv[ap] == '-')   /* an option */
    switch (op = *(argv[ap]+1))
      {
      case 's':
      case 'r':
      case 'S':
      case 'R':
        if (strlen(argv[ap]) > 2)      /* group number specified */
          {
          groupno = (int) strtol((argv[ap]+2),&fendp,10);
          if (fendp == (argv[ap] + 2))
            err_msg_die("Can't convert group number for %s\n",argv[ap]);
          }
        else
          switch (op)
            {
            case 's':
            case 'S':
              groupno = 2;
              break;
            case 'r':
            case 'R':
            default:
              groupno = 1;
              break;
            }
        rpars.maxgrp = imax(rpars.maxgrp,(groupno-1));
        groupid = argv[ap] + 1;
        switch (op)
          {
          case 's':    /* CpG list input */
          case 'r':    /* CpG list input */
            if (++ap > argc)
              err_msg_die("-%c needs file name\n",op);
            else
              (void) dm_appndfnam(&fnamlst,argv[ap],rpars.srcmode,groupno,groupid);
            break;
          case 'S':   /* sam file input */
          case 'R':   /* sam file input */
            if (++ap > argc)
              err_msg_die("-%c needs file name\n",op);
            else
              {
              if (rpars.compressed)
                (void) dm_appndfnam(&fnamlst,argv[ap],DM_src_bamfile,groupno,groupid);
              else
                (void) dm_appndfnam(&fnamlst,argv[ap],DM_src_samfile,groupno,groupid);
              rpars.needseq = 1;
              }
            break;
          }
        break;
      case 'g':
        if (++ap > argc)
          err_msg_die("-%c needs header string\n",op);
        else
          rpars.genomsqhdrstr = bas_strdup(argv[ap]);
        break;
      case 'G':   /* file name for chromosome spec file */
        if (++ap > argc)
          err_msg_die("-%c needs file name\n",op);
        else
          {
          rpars.chrinfoflnam = bas_strdup(argv[ap]);
          if ((rpars.chrinfofl = fopen(rpars.chrinfoflnam,"r")) == NULL)
            err_msg_die("Can't open chromosome information file '%s' for -%c\n",
                          rpars.chrinfoflnam,op);
          }
        break;
      case 'c':   /* a chromosome */
        if (++ap > argc)
          err_msg_die("-%c needs a chromosome identifier (%s..%s,%s,%s)\n",op,
                        any_chrno2str(rpars.maxchrno,rpars.havexy,1,1),
                        any_chrno2str(rpars.maxchrno,rpars.havexy,rpars.maxchrno-2,1),
                        any_chrno2str(rpars.maxchrno,rpars.havexy,rpars.maxchrno-1,1),
                        any_chrno2str(rpars.maxchrno,rpars.havexy,rpars.maxchrno,1));
        else
          rpars.uchrno = bas_strdup(argv[ap]);
        break;
      case 'p':         /* diff methylated regions */
      case 'P':         /* force Fisher's exact */
      case 'x':         /* chi-sq for valid fragments, else pairwise FE */
      case 'X':         /* force chi-sq for valid fragments */
      case 'q':         /* best pr of chi/fisher's exact */
      case 'Q':         /* all FE prs or chi */
      case 'e':         /* details for individual CpGs */
      case 'E':         /* details for individual nonzero CpGs */
      case 'l':         /* give listing */
      case 'L':         /* list non zero bins */
      case 'a':         /* anova */
      case 'A':         /* anova, show gtr meth group */
      case 'B':         /* anova, even more detail */
      case 'D':         /* sdev */
        switch (op)
          {
          case 'P':
            rpars.omode = dm_out_fisherforce;
            rpars.binprfun = dm_minpr4validfe;
            break;
          case 'x':
            rpars.omode = dm_out_chisq;
            rpars.binprfun = dm_chipr4validsamples;
            break;
          case 'X':
            rpars.omode = dm_out_chiforce;
            rpars.binprfun = dm_chipr4validsamples;
            break;
          case 'q':
            rpars.omode = dm_out_bestpairchi;
            rpars.binprfun = dm_minpr4chiorpair;
            break;
          case 'Q':
            rpars.omode = dm_out_allpairchi;
            rpars.binprfun = dm_minpr4chiorpair;
            break;
          case 'e':
            rpars.omode = dm_out_cpg_detail;
/*             rpars.binprfun = dm_minpr4chiorpair; */
            rpars.cpgdetails = 1;
            break;
          case 'E':
            rpars.omode = dm_out_cpg_nzdetl;
            rpars.cpgdetails = 1;
            break;
          case 'l':         /* give listing */
            rpars.omode = dm_out_list;
            break;
          case 'L':         /* list non zero bins */
            rpars.omode = dm_out_listnz;
            break;
          case 'D':         /* sdev proportions */
            rpars.omode = dm_out_sdev;
            break;
          case 'a':
          case 'A':
          case 'B':
            switch (op)
              {
              case 'A':
                rpars.omode = dm_out_anova_gtr;
                break;
              case 'B':
                rpars.omode = dm_out_anova_more;
                break;
              case 'a':
              default:
                rpars.omode = dm_out_anova;
                break;
              }
            rpars.binprfun = dm_anovapr4validsampls;
            break;
          default:
          case 'p':
            rpars.omode = dm_out_pairwise;
            rpars.binprfun = dm_minpr4validfe;
            break;
          }
        if (++ap < argc)
          {
          if (*argv[ap] != '-')
            bc_commalst2ints(argv[ap],&minfrag,&maxfrag);
          else
            ap--;
          }
        break;
      case 'm':       /* treat Complementary Cs as separate */
        rpars.cntcpgsas1 = 0;
        break;
      case 'f':       /* change default fold methylation test */
        if (++ap > argc)
          err_msg_die("-%c needs float value\n",op);
        else
          {
          rpars.folddiff = strtof(argv[ap],NULL);
          rpars.binchkmask = rpars.binchkmask | 1 << dm_chk_folddiff;
          }
        break;
      case 'C':       /* -C set minimum CpG no for bins */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          {
          rpars.mincpgs = (int) strtol(argv[ap],&fendp,10);
          rpars.binchkmask = rpars.binchkmask | 1 << dm_chk_cpgcnt;
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (rpars.mincpgs < 0)
              err_msg_die("Invalid %s min No.'%s'\n",dm_cpg_or_c(&rpars),argv[ap]);
          }
        break;
/*      case 'A':
        amlgmate = 1;
       break; */
      case 't':       /* change threshold for including bins in stats */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          {
          rpars.cntthreshold = (int) strtol(argv[ap],&fendp,10);
          if (!(rpars.binchkmask & 1 << dm_chk_mincpgcrit))  /* see -F */
            rpars.binchkmask = rpars.binchkmask | 1 << dm_chk_cntcrit;
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (rpars.cntthreshold < 0)
              err_msg_die("Invalid count threshold '%s'\n",argv[ap]);
          }
        break;
      case 'T':       /* change threshold for including hits/CpG in stats */
        if (++ap > argc)
          err_msg_die("-%c needs float value\n",op);
        else
          {
          rpars.minhitspercpg = (float) strtof(argv[ap],&fendp);
          rpars.binchkmask = rpars.binchkmask | 1 << dm_chk_cpghits;
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert float '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 'u':       /* change max threshold for including hits/CpG in stats */
        if (++ap > argc)
          err_msg_die("-%c needs float value\n",op);
        else
          {
          rpars.maxhitspercpg = (float) strtof(argv[ap],&fendp);
          rpars.binchkmask = rpars.binchkmask | 1 << dm_chk_maxcpghits;
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert float '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 'U':    /* Set pr threshold for output */
        if (++ap > argc)
          err_msg_die("-%c needs float value\n",op);
        else
          {
          rpars.prthreshold = (float) strtof(argv[ap],&fendp);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert float '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 'd':       /* enable/disable counts for each individual */
        rpars.displycnts = !rpars.displycnts;
        break;
      case 'k':       /* change max chromosome number */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          rpars.maxchrno = (int) strtol(argv[ap],NULL,10);
        break;
      case 'Y':       /* disable XY chromosome labelling */
        rpars.havexy = 0;
        break;
      case 'K':       /* change default bin cluster frequency */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          rpars.binclusterfreq = (int) strtol(argv[ap],NULL,10);
        break;
      case 'H':         /* chromosome prefix */
        if (++ap > argc)
          err_msg_die("-%c needs string value\n",op);
        else
          rpars.chrprefix = bas_strdup(argv[ap]);
        break;
      case 'I':       /* Minimum sample count for chi */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          {
          rpars.minsamples = (int) strtol(argv[ap],&fendp,10);
          rpars.binchkmask = rpars.binchkmask | 1 << dm_chk_indcnt;
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (rpars.cntthreshold < 0)
              err_msg_die("Invalid min sample count '%s'\n",argv[ap]);
          }
        break;
      case 'M':       /* minmethylation */
        if (++ap > argc)
          err_msg_die("-%c needs float value\n",op);
        else
          {
          rpars.minmethprop = (float) strtof(argv[ap],&fendp);
          rpars.binchkmask = rpars.binchkmask | 1 << dm_chk_minmeth;
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert float '%s' for -%c\n",argv[ap],op);
          else
            if ((rpars.minmethprop < 0.0) || (rpars.minmethprop > 1.0))
              err_msg_die("Invalid minimum methylation proportion: %f",rpars.minmethprop);
          }
        break;
      case 'F':        /* validity criterion is n cpghits exceed cntthreshold  */
        if (++ap > argc)
          err_msg_die("-%c needs integer value\n",op);
        else
          {
          rpars.cpgcntscrit = (int) strtol(argv[ap],&fendp,10);
          rpars.binchkmask = rpars.binchkmask | 1 << dm_chk_mincpgcrit;
/* Don't want cnt threshold applied in usual form */
          rpars.binchkmask = rpars.binchkmask & ~(1 << dm_chk_cntcrit);
          rpars.cpgdetails = 1;
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert integer '%s' for -%c\n",argv[ap],op);
          else
            if (rpars.cpgcntscrit < 0)
              err_msg_die("Invalid %s count criterion: %d",dm_cpg_or_c(&rpars),
                             rpars.cpgcntscrit);
          }
        break;
      case 'j':        /* join fragments */
        rpars.joinfrags = 1;
        break;
      case 'N':
        rpars.samfst3p2prv = 1;
        rpars.cpgdetails = 1;
        break;
      case 'W':       /* fixed width bins */
        if (++ap > argc)
          err_msg_die("-%c needs bin width\n",op);
        else
          {
          rpars.binwidth = (int) strtol(argv[ap],&fendp,10);
          rpars.binstyle = DM_bin_fixed;
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert int '%s' for -%c\n",argv[ap],op);
          if (rpars.binwidth == 0)
            err_msg_die("Invalid fixed bin width '%s'\n",argv[ap]);
          }
        break;
      case 'w':        /* ctlist input */
        rpars.srcmode = DM_src_ctlist;
        break;
      case 'z':        /* -R & -S compressed (bam) */
#ifdef NO_ZLIB
        err_msg_die("Compiled with NO_ZLIB set: can't read bam files\n");
#else
        rpars.compressed = 1;
#endif
        break;
      case 'Z':        /* -R & -S not compressed (sam) */
        rpars.compressed = 0;
        break;
      case 'n':        /* -n nuclease site file */
        if (++ap > argc)
          err_msg_die("-%c needs nuclease site file name\n",op);
        else
          {
          rpars.rstflnam = bas_strdup(argv[ap]);
          if ((rpars.rstfile = fopen(rpars.rstflnam,"r")) == NULL)
            err_msg_die("Can't open nuclease site file '%s' for -%c\n",
                          rpars.rstflnam,op);
          }
        break;
      case 'b':      /* sam file bufferlength */
        if (++ap > argc)
          err_msg_die("-%c needs bufferlength\n",op);
        else
          {
          rpars.sambuflen = (int) strtol(argv[ap],&fendp,10);
          if (fendp == argv[ap])    /* failed to read */
            err_msg_die("Can't convert int '%s' for -%c\n",argv[ap],op);
          }
        break;
      case 'J':    /* all Cs; non-CpG methylation */
        rpars.allcs = 1;
        break;
      case 'y':    /* user bin file */
        if (++ap >= argc)
          err_msg_die("-%c needs bin file name\n",op);
        else
          if ((rpars.binfile = fopen(argv[ap],"r")) == NULL)
            err_msg_die("Can't open bin data file '%s'\n",argv[ap]);
          else
            rpars.binstyle = DM_bin_userlist;
        break;
      case 'h':
        say_usage(stdout,argv[0]);
        exit(0);
        break;
      default:
        err_msg("Unknown Option: '%s\n",argv[ap]);
        say_usage(stderr,argv[0]);
        exit(1);
        break;
      }
  }
/* should be ready to go now */
/* build FSM */
if (rpars.binwidth <= 0)
  {
  if (rpars.rstfile != NULL)
    {
    (void) dm_readrstsites(rpars.rstfile,&rpars.sitelist);
    fclose(rpars.rstfile);
    }
  else
    (void) dm_appndsitestr(&rpars.sitelist,"C^CGG");
  sitp = rpars.sitelist;
  while (sitp != NULL)
    {
    fs_adddatprs(cpgfsmp,sitp->sitestr,(void *) sitp->offset);
    sitp = sitp->nxtsite;
    }
  }
if (rpars.allcs)
  fs_adddatprs(cpgfsmp,"C",NULL);
else
  fs_adddatprs(cpgfsmp,"CG",NULL);
(void) fs_bldfsm(cpgfsmp,WLU_CASEIND,0,0,tr_int2nares,fs_chkinpstr,fs_shed2lurec);
if (rpars.chrinfofl != NULL)
  rpars.maxchrno = dm_scanchrinfofl(rpars.chrinfofl,NULL,NULL);
/* create chromosome information array */
rpars.chrbininfo = (DM_CHR_INFO *) getmemory(sizeof(DM_CHR_INFO)*rpars.maxchrno,
                                         "Chr bin info");
/* create chromosome lookup structure and init it */
rpars.chridlu = (WRD_LUSTRCT *) getmemory(sizeof(WRD_LUSTRCT),"ChrID lu struct");
wlu_initlustrct(rpars.chridlu,WLU_CASEIND,-1);
/* init these */
for (chrno = 0; chrno < rpars.maxchrno; chrno++)
  {
  (rpars.chrbininfo+chrno)->chrlen = 0;
  (rpars.chrbininfo+chrno)->chromseq = NULL;
  (rpars.chrbininfo+chrno)->binlst = NULL;
  (rpars.chrbininfo+chrno)->beltlst = NULL;
  (rpars.chrbininfo+chrno)->chrflname = NULL;
  (rpars.chrbininfo+chrno)->chrid = NULL;
  }
if ((rpars.genomsqhdrstr != NULL) && (rpars.chrinfofl == NULL))
  dm_bldchrfilenames(rpars.genomsqhdrstr,&rpars);
else
  if (rpars.chrinfofl != NULL)
    {
    rewind(rpars.chrinfofl);
    switch (rpars.binstyle)
      {
      case DM_bin_userlist:
        break;        
      case DM_bin_fixed:
      case DM_bin_restenz:
      default:
        break;
      }
    (void) dm_scanchrinfofl(rpars.chrinfofl,rpars.chrbininfo,rpars.chridlu);
    fclose(rpars.chrinfofl);
    }
  else
    {
    err_msg("No sequence file information via -g or -G options\n");
    say_usage(stderr,argv[0]);
    exit(1);
    }
if (rpars.maxchrno > 0)
  {
  if (rpars.binstyle == DM_bin_userlist)
    (void) dm_readbinfile(rpars.binfile,&rpars);
  chrcnt = dm_redrepgenomsqs(&rpars,cpgfsmp);
  if (debuglevel > RBC_dbg_none)
    fprintf(stdout,"%d chromosomes read\n",chrcnt);
  if (amlgmate)
    (void) bc_amlgmatelocpgbins(rpars.chrbininfo,rpars.maxchrno,minfrag,maxfrag,
                                  rpars.mincpgs);
  else
    switch (rpars.binstyle)
      {
      case DM_bin_fixed:
      case DM_bin_userlist:
        break;  /* user-defined bins shouldn't fail */
      case DM_bin_restenz:
      default:
        (void) bc_losefailedbins(&rpars,rpars.chrbininfo,rpars.maxchrno,minfrag,maxfrag,
                                   rpars.mincpgs);
        break;
      }
  if (rpars.joinfrags)
    (void) dm_joinadjacntfrags(&rpars,rpars.chrbininfo,dm_jointestsucceedalways);
  fnamptr = fnamlst;
  for (chrno = 0; chrno < rpars.maxchrno; chrno++)
    if ((rpars.uchrno == NULL) || (wlu_chkwrd(rpars.chridlu,rpars.uchrno) == chrno))
      (void) dm_scanbinlsts2eltlsts(&(rpars.chrbininfo+chrno)->beltlst,
                                      (rpars.chrbininfo+chrno)->binlst,
                                      rpars.binclusterfreq);
/* create an array of groupids, and null it */
  rpars.groupidlist = (char **) getmemory(sizeof(char *)*(rpars.maxgrp+1),"GroupIDlist");
  for (i=0;i<=rpars.maxgrp;i++)
    *(rpars.groupidlist+i) = NULL;
  while (fnamptr != NULL)
    {
    switch (fnamptr->srcmode)
      {
#ifndef NO_ZLIB
      case DM_src_bamfile:
        if ((rpars.bamrdp = bf_opnchkbamflmsg(fnamptr->flnam,err_msg_die)) == NULL)
          err_msg_die("Can't open BAM file %s\n",fnamptr->flnam);
        else
          {
          ecnts = dm_readnxtsrcfl(rpars.chrbininfo,NULL,&rpars,fnamptr->srcmode,fnamptr->sgroup);
          bf_disposerundata(rpars.bamrdp);
          rpars.bamrdp = NULL;
          }
        break;
#endif
      case DM_src_samfile:
      case DM_src_cpglist:
      case DM_src_ctlist:
      default:
        if ((srcfl = fopen(fnamptr->flnam,"r")) == NULL)
          err_msg("Can't open %s file %s\n",
                    ((fnamptr->srcmode==DM_src_cpglist)?"CpG":"SAM"),fnamptr->flnam);
        else
          {
          ecnts = dm_readnxtsrcfl(rpars.chrbininfo,srcfl,&rpars,fnamptr->srcmode,fnamptr->sgroup);
          fclose(srcfl);
          }
        break;
      }
    if (*(rpars.groupidlist+fnamptr->sgroup-1) == NULL)
      *(rpars.groupidlist+fnamptr->sgroup-1) = fnamptr->groupid;
    fnamptr = fnamptr->nxtfelt;
    }
  dm_resetcntlststohead(rpars.chrbininfo,&rpars);
  switch (rpars.omode)
    {
    case dm_out_none:
      break;
    case dm_out_list:
    case dm_out_listnz:
      dm_listoutdata(stdout,&rpars);
      break;
    case dm_out_anova:
    case dm_out_anova_gtr:     /* need to allocate storage for structures */
    case dm_out_anova_more:
      rpars.anova_dat = (DM_ANOVA_DAT *) getmemory(sizeof(DM_ANOVA_DAT),"AnovaDat");
      rpars.anova_dat->maxgroup = rpars.maxgrp + 1;
      rpars.anova_dat->ssx = (double *) getmemory(sizeof(double)*rpars.anova_dat->maxgroup,"AnovaSSX");
      rpars.anova_dat->sx = (double *) getmemory(sizeof(double)*rpars.anova_dat->maxgroup,"AnovaSx");
      rpars.anova_dat->grpcnts = (int *) getmemory(sizeof(int)*rpars.anova_dat->maxgroup,"AnovaGrpCnts");
    default:
      dm_headdmethout(stdout,&rpars);
      dm_scanallchrombins(stdout,&rpars);
      break;
    }
  exit(0);
  }
else
  err_msg_die("Need genome sequence file header string\n");
exit(0);
}
