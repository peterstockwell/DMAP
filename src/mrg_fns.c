/* mrg_fns.c: some routines to support mkrrgenome and related
programs */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqmat_fns.h"
#include "mrg_fns.h"
#include "rmapbsbed2cpg.h"
#include "rbs_fns.h"

MRG_REGN_ELT *mrg_appndrgnelt(MRG_REGN_ELT **lstrt,
                             int rgstart,
                             int rgstop,
                             int cpgs)
/* create and append a new element to *lstrt,
 Return address of new element */
{
MRG_REGN_ELT *prev, *end_ptr;

if (lstrt != NULL)
  {
  prev = end_ptr = *lstrt;
  while (end_ptr != NULL)
    {
    prev = end_ptr;
    end_ptr = end_ptr->nxtregn;
    }
  end_ptr = (MRG_REGN_ELT *) getmemory(sizeof(MRG_REGN_ELT),"region elt");
  end_ptr->nxtregn = NULL;
  end_ptr->rstart = rgstart;
  end_ptr->rstop = rgstop;
  end_ptr->cpgcnt = cpgs;
  if (*lstrt == NULL)
    {
    *lstrt = end_ptr;
    end_ptr->prvregn = NULL;
    }
  else
    {
    prev->nxtregn = end_ptr;
    end_ptr->prvregn = prev;
    }
  return(end_ptr);
  }
else
  return(NULL);
}

void mrg_delregnelt(MRG_REGN_ELT *ep,
                    MRG_REGN_ELT **lstrt)
/* delete ep from list *lstrt */
{
MRG_REGN_ELT *pt;

if (ep != NULL)
  {
  if ((pt = ep->prvregn) == NULL)
    *lstrt = ep->nxtregn;
  else
    pt->nxtregn = ep->nxtregn;
  if ((pt = ep->nxtregn) != NULL)
    pt->prvregn = ep->prvregn;
  memfree(ep);
  }
}

void mrg_clrallregnelts(MRG_REGN_ELT **lstrt)
  /* iteratively delete all of lstrt */
{
while (*lstrt != NULL)
  mrg_delregnelt(*lstrt,lstrt);
}

int mrg_cntregnelts(MRG_REGN_ELT *clst)
  /* recursively count list elements */
{
if (clst == NULL)
  return(0);
else
  return(mrg_cntregnelts(clst->nxtregn) + 1);
}

int mrg_sumregnelts(MRG_REGN_ELT *rlst)
  /* recursively sum lengths of elements */
{
if (rlst == NULL)
  return(0);
else
  return(mrg_sumregnelts(rlst->nxtregn) + rlst->rstop - rlst->rstart + 1);
}

int mrg_sumcpgregnelts(MRG_REGN_ELT *rlst)
  /* recursively ... */
{
if (rlst == NULL)
  return(0);
else
  return(mrg_sumcpgregnelts(rlst->nxtregn) + rlst->cpgcnt);
}

MRG_REGN_ELT *mrg_regnelt4pos(MRG_REGN_ELT *rlst,
                            int pos)
/* return the region element corresponding to pos, if any */
{
MRG_REGN_ELT *rp;

rp = rlst;
while (rp != NULL)
  if (rbc_intinrng(rp->rstart,pos,rp->rstop))
    return(rp);
  else
    rp = rp->nxtregn;
/* fell off end, reurn NULL */
return(NULL);
}

