/* cmaths.h: C routine headers for mathematical operations */
/* contents: chi_sq  - chi-square density value for given chi and df
             chi_sq_pr - cumulative probability for given chi and df
             setIEEEhandlr - set sudden death exception handler for sun */

double chi_sq(int nu,
              double chi);
/* return chi_sq(x) for nu df
  Ref: Based on Abramowitz and Stegun, "Handbook of Mathematical Functions",
  National Bureau of Standards, 1970 */

double chi_sq_pr(int nu,
                 double chi);
/* returns integral probability for chi value on nu degrees of freedom
  Ref: Based on Abramowitz and Stegun, "Handbook of Mathematical Functions",
  National Bureau of Standards, 1970  */

#ifdef __sun

void setIEEEhandlr();
  /* set IEEE handler to abort on FP exceptions on Sun/Solaris */

#endif
