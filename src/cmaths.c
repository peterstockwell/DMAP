/* cmaths.c: C routines for mathematical operations */
/* contents: chi_sq  - chi-square density value for given chi and df
             chi_sq_pr - cumulative probability for given chi and df */

/*
 *  Revised by GHJ, 1-Sep-2005, to use isnan() to test for limits for
 *  sums and products. There are smarter ways of doing this, but I
 *  just want this to work...
 *  This doesn't warrant the effort, but it seems to me that it'd be
 *  good if this code made use of some of the maths libraries out
 *  there. If nothing else, their code is polished to make good use
 *  of the CPU. 
 *  The "correct" test using __STDC_VERSION__ doesn't work as gcc is
 *  cautious since they don't _quite_ implement all of C99, hence use
 *  of __HAS_ISNAN__, below.
 */

/* revised, PAS, Nov-2011 - use tgamma() rather than calculate
own. - shortens the thing appreciably */

#include <math.h>
#include <limits.h>
#include <float.h>

#ifdef __sun
#include <floatingpoint.h>      /* to allow tracing of FP exceptions */
#endif

/* Let us use isnan() to do math checking on this machine  GHJ 1-Sep-2005 */
#define __HAS_ISNAN__

int chksuminrng(double term1,
                double term2)
/* check that adding term1 to term2 will not exceed double floating range
of machine */
{
#ifdef __HAS_ISNAN__
return (! isnan(term1+term2));
#else
if ((term1 > 0.0) && (term2 > 0.0))
  return(((DBL_MAX - term2) >= term1));
else
  if ((term1 < 0.0) && (term2 < 0.0))
    return(chksuminrng(-term1,-term2));
  else
    return(1);
#endif
}

int chkprodinrng(double term1,
                 double term2)
/* check that term1*term2 will not exceed integral range of machine */
{
#ifdef __HAS_ISNAN__
  return (! isnan( term1 * term2 ) );
#else

/*
 Testing code added GHJ 1-Sep-2005
 C99-style check. See Harbison&Steele for details.
 Note this effectively checks *after* an overflow, so code needs to be reworked...?
 Probably should apply C99 maths constants that are meant to make this simple,
 rather than the "work it out" logic this and chksuminrng() use.
if (isnan(term1) || isnan(term2))
{
  printf("GHJ: Gotcha!\n");
  return(1);
}
*/

if ((term1 == 0.0) || (term2 == 0.0))
  return(1);
else
  if ((term1 > 0.0) && (term2 > 0.0))
    if (term2 <= 1.0)
      return(1);
    else
      return((term1 <= DBL_MAX/term2));
  else
    if (term1 < 0.0)
      return(chkprodinrng(-term1,term2));
    else
      return(chkprodinrng(term1,-term2));
#endif
}

int chisummatn(double x,
               int nu,
               double *sret)
/* to derive the summation term for chi-sq distribution.  Add successive
  terms iteratively until difference between successive sums is below
  precision of machine.  Check for numerical overflow and return
  *sret = DBL_MAX and 0 if it occurs */
{
double prvsum;
double newtrm;
double sum;
double newmlt;

prvsum = 0.0;
newtrm = x/(nu + 2);
sum = 1.0 + newtrm;
while (sum != prvsum)
  {
  prvsum = sum;
  nu += 2;
  if (chkprodinrng(newtrm,(newmlt = x /((double)(nu + 2)))))
    newtrm *= newmlt;
  else
    {
    if (sum >= 0.0)
      *sret = DBL_MAX;
    else
      *sret = -DBL_MAX;
    return(0);
    }
  if (chksuminrng(sum,newtrm))
    sum += newtrm;
  else
    {
    if (sum >= 0.0)
      *sret = DBL_MAX;
    else
      *sret = -DBL_MAX;
    return(0);
    }
  }
*sret = sum;
return(1);
}

double chi_sq(int nu,
              double chi)
/* return chi_sq(x) for nu df
  Ref: Based on Abramowitz and Stegun, "Handbook of Mathematical Functions",
  National Bureau of Standards, 1970 */
{
return(pow((double) chi,(((double) nu)/2.0 - 1.0)) * exp(-((double) chi/2.0))
       /(pow(2.0,(((double)nu)/2.0)) * tgamma((double) nu/2.0)));
}

double chi_sq_pr(int nu,
                 double chi)
/* returns integral probability for chi value on nu degrees of freedom
  Ref: Based on Abramowitz and Stegun, "Handbook of Mathematical Functions",
  National Bureau of Standards, 1970 */
{
double csp;
double csm;

if (chisummatn(chi,nu,&csm))
  if ((csp = (chi * chi_sq(nu,chi) * csm /(((double)nu)/2.0))) > 1.0) 
    return(1.0);
  else
    return(csp);
else
  return(1.0);
}

#ifdef __sun

void setaIEEEhndlr(char *hndlr)
/* set handler hndlr to abort */
{
int rv;

if ((rv = ieee_handler("set",hndlr,SIGFPE_ABORT)) != 0)
  {
  fprintf(stderr,"Invalid IEEE %s handler call %d\n",hndlr,rv);
  exit(1);
  }
}

void setIEEEhandlr()
  /* set IEEE handler to abort on FP exceptions on Sun/Solaris */
{
setaIEEEhndlr("invalid");
/* setaIEEEhndlr("inexact"); */
setaIEEEhndlr("overflow");
setaIEEEhndlr("division");
}

#endif
