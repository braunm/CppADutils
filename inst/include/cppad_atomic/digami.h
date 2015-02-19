#include <Rcpp.h>

int digami(double *d__, double *x, double *p, double *gplog,
	   double * gp1log, double *psip, double *psip1,
	   double *psidp, double *psidp1, int *ifault)
{
    /* Initialized data */

    static double e = 1e-6f;
    static double oflo = 1e30f;
    static double tmax = 100.f;
    static double vsmall = 1e-30f;
    static double zero = 0.f;
    static double one = 1.f;
    static double two = 2.f;

    /* System generated locals */
    double r__1;

    /* Builtin functions */
    //   double log(doubledouble), exp(doubledouble);

    /* Local variables */
    static double a, b, c__, f;
    static int i__;
    static double s;
    static int i2;
    static double s0, an, cp, dp[6], pn[6], pm1, cpc, dfp, dpp[6], cpp, dsp, 
	    dfpp, dspp, term, xlog, tmaxp;


/*     ALGORITHM AS 187  APPL. STATIST. (1982) VOL.31, NO.3 */

/*     Computes derivatives of the incomplete gamma integral for positive */
/*     parameters, X, P, using a series expansion if P > X or X <= 1, and */
/*     a continued fraction expansion otherwise. */

/*     Calculation of D(4) in line 60 corrected 5 October 1993. */

/*     N.B. The user must input values of the incomplete gamma, digamma */
/*          and trigamma functions.  These can be obtained using AS 239 */
/*          (or 32), AS 103 and AS 121 respectively. */

    /* Parameter adjustments */
    --d__;

    /* Function Body */

    *ifault = 0;

/*     Derivatives with respect to X */

    pm1 = *p - one;
    xlog = log(*x);
    d__[1] = exp(-(*gplog) + pm1 * xlog - *x);
    d__[2] = d__[1] * (pm1 / *x - one);
    d__[5] = d__[1] * (xlog - *psip);

/*     Derivatives with respect to P */

    if (*x > one && *x >= *p) {
	goto L30;
    }

/*     Series expansion */

    f = exp(*p * xlog - *gp1log - *x);
    dfp = f * (xlog - *psip1);
    dfpp = dfp * dfp / f - f * *psidp1;

    tmaxp = tmax + *p;
    c__ = one;
    s = one;
    cp = zero;
    cpp = zero;
    dsp = zero;
    dspp = zero;
    a = *p;
L1:
    a += one;
    cpc = cp / c__;
    cp = cpc - one / a;
/* Computing 2nd power */
    r__1 = a;
    cpp = cpp / c__ - cpc * cpc + one / (r__1 * r__1);
    c__ = c__ * *x / a;
    cp *= c__;
    cpp = cpp * c__ + cp * cp / c__;
    s += c__;
    dsp += cp;
    dspp += cpp;
    if (a > tmaxp) {
	goto L1001;
    }
    if (c__ > e * s) {
	goto L1;
    }
    d__[6] = s * f;
    d__[3] = s * dfp + f * dsp;
    d__[4] = s * dfpp + two * dfp * dsp + f * dspp;
    return 0;

/*     Continued fraction expansion */

L30:
    f = exp(*p * xlog - *gplog - *x);
    dfp = f * (xlog - *psip);
    dfpp = dfp * dfp / f - f * *psidp;

    a = pm1;
    b = *x + one - a;
    term = zero;
    pn[0] = one;
    pn[1] = *x;
    pn[2] = *x + one;
    pn[3] = *x * b;
    s0 = pn[2] / pn[3];
    for (i__ = 1; i__ <= 4; ++i__) {
	dp[i__ - 1] = zero;
	dpp[i__ - 1] = zero;
/* L31: */
    }
    dp[3] = -(*x);

L32:
    a -= one;
    b += two;
    term += one;
    an = a * term;
    pn[4] = b * pn[2] + an * pn[0];
    pn[5] = b * pn[3] + an * pn[1];
    dp[4] = b * dp[2] - pn[2] + an * dp[0] + pn[0] * term;
    dp[5] = b * dp[3] - pn[3] + an * dp[1] + pn[1] * term;
    dpp[4] = b * dpp[2] + an * dpp[0] + two * (term * dp[0] - dp[2]);
    dpp[5] = b * dpp[3] + an * dpp[1] + two * (term * dp[1] - dp[3]);

    if (abs(pn[5]) < vsmall) {
	goto L35;
    }
    s = pn[4] / pn[5];
    c__ = (r__1 = s - s0, abs(r__1));
    if (c__ * *p > e) {
	goto L34;
    }
    if (c__ <= e * s) {
	goto L42;
    }

L34:
    s0 = s;
L35:
    for (i__ = 1; i__ <= 4; ++i__) {
	i2 = i__ + 2;
	dp[i__ - 1] = dp[i2 - 1];
	dpp[i__ - 1] = dpp[i2 - 1];
	pn[i__ - 1] = pn[i2 - 1];
/* L36: */
    }

    if (term > tmax) {
	goto L1001;
    }
    if (abs(pn[4]) < oflo) {
	goto L32;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	dp[i__ - 1] /= oflo;
	dpp[i__ - 1] /= oflo;
	pn[i__ - 1] /= oflo;
/* L41: */
    }
    goto L32;

L42:
    d__[6] = one - f * s;
    dsp = (dp[4] - s * dp[5]) / pn[5];
    dspp = (dpp[4] - s * dpp[5] - two * dsp * dp[5]) / pn[5];
    d__[3] = -f * dsp - s * dfp;
    d__[4] = -f * dspp - two * dsp * dfp - s * dfpp;
    return 0;

/*     Set fault indicator */

L1001:
    *ifault = 1;
    return 0;
} /* digami_ */
