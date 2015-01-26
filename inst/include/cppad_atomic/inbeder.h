/* Routine to compute 1st and 2nd
   derivatives of the regularized incomplete
   beta function.

Converted (using f2c) from Fortran code from:

Boik, Robert J. and James F. Robinson-Cox,
"Derivatives of the Incomplete Beta Function"
Journal of Statistical Software 3(1), March 1998.

*/

/*
Declares 3 private member functions for incbeta:
-- inbeder
-- derconf_
-- subd_

Included in definition of incbeta_cl
*/



/* Table of constant values */

// static int c__1 = 1;
// static double c_b4 = 2.;
// static double c_b5 = 10.;
// static int c__5 = 5;
// static int c__6 = 6;


int inbeder(double *x, double *p, double *q, 
	    double *psi, double *der, int *nappx, double *errapx, 
	    int *ifault)
{
    /* Initialized data */
  using std::abs;
  using std::min;
  using std::max;

    static double zero = 0.;
    static double one = 1.;
    static double err = 1e-12;
    static int maxappx = 200;
    static int minappx = 3;
    static double prmin = 1e-24;

    /* System generated locals */
    double d__1, d__2, d__3;

    /* Builtin functions */
    double log(double), exp(double);

    /* Local variables */
    static double c__[6], d__, e;
    static int i__, n;
    static double w, c0, d1, x1, an[6], bn[6], pa, pb;
    static int ii;
    static double dr[6], pp, rn, qq, pr, an1[6], an2[6], pa1, pb1, bn1[6],
	     bn2[6], pab, dan[6], dbn[6];
    static int iii;
    static double omx, pab1, logx1, prmax, logomx, der_old__[6];
    // extern /* Subroutine */ int derconf_(double *, double *, 
    //	    double *, int *, double *, double *);


/*             x: Input argument -- value to which beta function is integrated */
/*             p,q: Input arguments -- beta shape parameters */
/*             psi: input arguments -- vector of length 7 */
/*                                    psi(1) = log[Beta(p,q)] */
/*                                    psi(2) = digamma(p) */
/*                                    psi(3) = trigamma(p) */
/*                                    psi(4) = digamma(q) */
/*                                    psi(5) = trigamma(q) */
/*                                    psi(6) = digamma(p+q) */
/*                                    psi(7) = trigamma(p+q) */
/*             der: output -- vector of length 6 */
/*                            der(1) = I (incomplete beta function) */
/*                            der(2) = Ip */
/*                            der(3) = Ipp */
/*                            der(4) = Iq */
/*                            der(5) = Iqq */
/*                            der(6) = Ipq */
/*             nappx: output -- highest order approximant evaluated */
/*             errapx: output -- approximate size of maximum absolute error */
/*                               of computed derivatives */
/*             ifault: output -- error indicator */
/*                               ifault = 0: no error */
/*                               ifault = 1: x outside of (0,1) */
/*                               ifault = 2: p less that 0 */
/*                               ifault = 3: q less than 0 */
/*                               ifault = 4: derivatives set to 0 because I=0 */
/*                                           or I=1 */
/*                               ifault = 5: evaluation stopped after maxappx */
/*                                           terms */

    /* Parameter adjustments */
    --der;
    --psi;

    /* Function Body */
    prmax = one - err;

/*          Initialize derivative vectors */
/*          and check for admissability of input arguments */

    for (i__ = 1; i__ <= 6; ++i__) {
	der_old__[i__ - 1] = zero;
	c__[i__ - 1] = zero;
	dr[i__ - 1] = one;
	der[i__] = zero;
	an2[i__ - 1] = zero;
	bn2[i__ - 1] = zero;
	an1[i__ - 1] = zero;
/* L10: */
	bn1[i__ - 1] = zero;
    }
    an1[0] = one;
    bn1[0] = one;
    an2[0] = one;
    pab = psi[6];
    pab1 = psi[7];
    *ifault = 0;
    *nappx = 0;


    if (*x <= zero || *x >= one) {
	*ifault = 1;
	return 0;
    }
    if (*p <= zero) {
	*ifault = 2;
	return 0;
    }
    if (*q <= zero) {
	*ifault = 3;
	return 0;
    }
 

/*          Use I(x,p,q) = 1- I(1-x,q,p) if x > p/(p+q) */

    if (*x <= *p / (*p + *q)) {
	ii = 1;
	x1 = *x;
	omx = one - *x;
	pp = *p;
	qq = *q;
	pa = psi[2];
	pb = psi[4];
	pa1 = psi[3];
	pb1 = psi[5];
    } else {
	ii = 2;
	x1 = one - *x;
	omx = *x;
	pp = *q;
	qq = *p;
	pa = psi[4];
	pb = psi[2];
	pa1 = psi[5];
	pb1 = psi[3];
    }
    w = x1 / omx;
    logx1 = log(x1);
    logomx = log(omx);

/*          Compute derivatives of K(x,p,q) = x^p(1-x)^(q-1)/[p beta(p,q)] */

    c__[0] = pp * logx1 + (qq - 1) * logomx - psi[1] - log(pp);
    c0 = exp(c__[0]);
    c__[1] = logx1 - one / pp - pa + pab;
/* Computing 2nd power */
    d__1 = c__[1];
/* Computing 2nd power */
    d__2 = pp;
    c__[2] = d__1 * d__1 + one / (d__2 * d__2) - pa1 + pab1;
    c__[3] = logomx - pb + pab;
/* Computing 2nd power */
    d__1 = c__[3];
    c__[4] = d__1 * d__1 - pb1 + pab1;
    c__[5] = c__[1] * c__[3] + pab1;

/*          Set counter and begin iteration */

    n = 0;
L20:
    ++n;

/*          Compute derivatives of an and bn with respect to p and/or q */

    derconf_(&pp, &qq, &w, &n, an, bn);

/*          Use forward recurrance relations to compute An, Bn, */
/*          and their derivatives */

    dan[0] = an[0] * an2[0] + bn[0] * an1[0];
    dbn[0] = an[0] * bn2[0] + bn[0] * bn1[0];
    dan[1] = an[1] * an2[0] + an[0] * an2[1] + bn[1] * an1[0] + bn[0] * an1[1]
	    ;
    dbn[1] = an[1] * bn2[0] + an[0] * bn2[1] + bn[1] * bn1[0] + bn[0] * bn1[1]
	    ;
    dan[2] = an[2] * an2[0] + an[1] * 2 * an2[1] + an[0] * an2[2] + bn[2] * 
	    an1[0] + bn[1] * 2 * an1[1] + bn[0] * an1[2];
    dbn[2] = an[2] * bn2[0] + an[1] * 2 * bn2[1] + an[0] * bn2[2] + bn[2] * 
	    bn1[0] + bn[1] * 2 * bn1[1] + bn[0] * bn1[2];
    dan[3] = an[3] * an2[0] + an[0] * an2[3] + bn[3] * an1[0] + bn[0] * an1[3]
	    ;
    dbn[3] = an[3] * bn2[0] + an[0] * bn2[3] + bn[3] * bn1[0] + bn[0] * bn1[3]
	    ;
    dan[4] = an[4] * an2[0] + an[3] * 2 * an2[3] + an[0] * an2[4] + bn[4] * 
	    an1[0] + bn[3] * 2 * an1[3] + bn[0] * an1[4];
    dbn[4] = an[4] * bn2[0] + an[3] * 2 * bn2[3] + an[0] * bn2[4] + bn[4] * 
	    bn1[0] + bn[3] * 2 * bn1[3] + bn[0] * bn1[4];
    dan[5] = an[5] * an2[0] + an[1] * an2[3] + an[3] * an2[1] + an[0] * an2[5]
	     + bn[5] * an1[0] + bn[1] * an1[3] + bn[3] * an1[1] + bn[0] * an1[
	    5];
    dbn[5] = an[5] * bn2[0] + an[1] * bn2[3] + an[3] * bn2[1] + an[0] * bn2[5]
	     + bn[5] * bn1[0] + bn[1] * bn1[3] + bn[3] * bn1[1] + bn[0] * bn1[
	    5];

/*          Scale derivatives to prevent overflow */

    rn = dan[0];
    iii = 1;
    if (abs(dbn[0]) > abs(dan[0])) {
	rn = dbn[0];
	iii = 2;
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	an1[i__ - 1] /= rn;
/* L33: */
	bn1[i__ - 1] /= rn;
    }
    for (i__ = 2; i__ <= 6; ++i__) {
	dan[i__ - 1] /= rn;
/* L34: */
	dbn[i__ - 1] /= rn;
    }
    if (iii == 1) {
	dbn[0] /= dan[0];
	dan[0] = one;
    } else {
	dan[0] /= dbn[0];
	dbn[0] = one;
    }

/*          Compute components of derivatives of the nth approximant */

    dr[0] = dan[0] / dbn[0];
    rn = dr[0];
    dr[1] = (dan[1] - rn * dbn[1]) / dbn[0];
/* Computing 2nd power */
    d__1 = dbn[1];
/* Computing 2nd power */
    d__2 = dbn[0];
    dr[2] = (dan[1] * -2 * dbn[1] + rn * 2 * (d__1 * d__1)) / (d__2 * d__2) + 
	    (dan[2] - rn * dbn[2]) / dbn[0];
    dr[3] = (dan[3] - rn * dbn[3]) / dbn[0];
/* Computing 2nd power */
    d__1 = dbn[3];
/* Computing 2nd power */
    d__2 = dbn[0];
    dr[4] = (dan[3] * -2 * dbn[3] + rn * 2 * (d__1 * d__1)) / (d__2 * d__2) + 
	    (dan[4] - rn * dbn[4]) / dbn[0];
/* Computing 2nd power */
    d__1 = dbn[0];
    dr[5] = (-dan[1] * dbn[3] - dan[3] * dbn[1] + rn * 2 * dbn[1] * dbn[3]) / 
	    (d__1 * d__1) + (dan[5] - rn * dbn[5]) / dbn[0];

/*          Save terms corresponding to approximants n-1 and n-2 */

    for (i__ = 1; i__ <= 6; ++i__) {
	an2[i__ - 1] = an1[i__ - 1];
	an1[i__ - 1] = dan[i__ - 1];
	bn2[i__ - 1] = bn1[i__ - 1];
/* L30: */
	bn1[i__ - 1] = dbn[i__ - 1];
    }

/*          Check if I < prmin or I > prmax */

    if (dr[0] > zero) {
	pr = exp(c__[0] + log(dr[0]));
    } else {
	pr = zero;
    }
    der[1] = pr;
    if (pr < prmin || pr > prmax) {
	*errapx = (d__1 = der_old__[0] - pr, abs(d__1));
	if (*errapx <= err) {
	    *ifault = 4;
	    for (i__ = 2; i__ <= 6; ++i__) {
/* L72: */
		der[i__] = zero;
	    }
	    goto L75;
	}
    }

/*          Compute nth approximants */

    der[2] = pr * c__[1] + c0 * dr[1];
    der[3] = pr * c__[2] + c0 * 2 * c__[1] * dr[1] + c0 * dr[2];
    der[4] = pr * c__[3] + c0 * dr[3];
    der[5] = pr * c__[4] + c0 * 2 * c__[3] * dr[3] + c0 * dr[4];
    der[6] = pr * c__[5] + c0 * c__[3] * dr[1] + c0 * c__[1] * dr[3] + c0 * 
	    dr[5];

/*          Check for convergence, check for maximum and minimum iterations. */

    d__ = zero;
    *errapx = zero;
    for (i__ = 1; i__ <= 6; ++i__) {
/* Computing MAX */
	d__2 = err, d__3 = (d__1 = der[i__], abs(d__1));
	d1 = max(d__2,d__3);
	e = (d__1 = der_old__[i__ - 1] - der[i__], abs(d__1));
	d1 = e / d1;
	if (d1 > d__) {
	    d__ = d1;
	}
	if (e > *errapx) {
	    *errapx = e;
	}
/* L92: */
	der_old__[i__ - 1] = der[i__];
    }
    if (n < minappx) {
	d__ = one;
    }
    if (n >= maxappx) {
	d__ = zero;
	*ifault = 5;
    }
    if (d__ >= err) {
	goto L20;
    }
L75:

/*          Adjust results if I(x,p,q) = 1- I(1-x,q,p) was used */

    if (ii == 2) {
	der[1] = one - der[1];
	c0 = der[2];
	der[2] = -der[4];
	der[4] = -c0;
	c0 = der[3];
	der[3] = -der[5];
	der[5] = -c0;
	der[6] = -der[6];
    }
    *nappx = n;
    return 0;
} /* inbeder_ */




/* Subroutine */ int derconf_(double *p, double *q, double *w, 
	int *n, double *an, double *bn)
{
    /* Initialized data */

    static double zero = 0.;
    static double one = 1.;
    static double two = 2.;

    static double f, t1, t2, t3, t4;
 

/*          Compute derivatives of an and bn with respect to p and/or q */

    /* Parameter adjustments */
    --bn;
    --an;

    /* Function Body */
    f = *w * *q / *p;
    if (*n == 1) {
	t1 = one - one / (*p + one);
	t2 = one - one / *q;
	t3 = one - two / (*p + two);
	t4 = one - two / *q;
	an[1] = t1 * t2 * f;
	an[2] = -an[1] / (*p + one);
	an[3] = -two * an[2] / (*p + one);
	an[4] = t1 * f / *q;
	an[5] = zero;
	an[6] = -an[4] / (*p + one);
	bn[1] = one - t3 * t4 * f;
	bn[2] = t3 * t4 * f / (*p + two);
	bn[3] = -two * bn[2] / (*p + two);
	bn[4] = -t3 * f / *q;
	bn[5] = zero;
	bn[6] = -bn[4] / (*p + two);
    } else {
	subd_(n, p, q, &f, &an[1], &bn[1]);
    }
    return 0;
} /* derconf_ */




/* Subroutine */ int subd_(int *in, double *p, double *q, 
	double *f, double *an, double *bn)
{
    /* Initialized data */

    static double c36 = 36.;
    static double c38 = 38.;
    static double c40 = 40.;
    static double c44 = 44.;
    static double c48 = 48.;
    static double c50 = 50.;
    static double c51 = 51.;
    static double c52 = 52.;
    static double c53 = 53.;
    static double c54 = 54.;
    static double c60 = 60.;
    static double c64 = 64.;
    static double c65 = 65.;
    static double c69 = 69.;
    static double c70 = 70.;
    static double c72 = 72.;
    static double c77 = 77.;
    static double c80 = 80.;
    static double c87 = 87.;
    static double c88 = 88.;
    static double c96 = 96.;
    static double c104 = 104.;
    static double c128 = 128.;
    static double c130 = 130.;
    static double c144 = 144.;
    static double c155 = 155.;
    static double c192 = 192.;
    static double c224 = 224.;
    static double c240 = 240.;
    static double c288 = 288.;
    static double c0 = 0.;
    static double c1 = 1.;
    static double c2 = 2.;
    static double c3 = 3.;
    static double c4 = 4.;
    static double c5 = 5.;
    static double c6 = 6.;
    static double c7 = 7.;
    static double c8 = 8.;
    static double c9 = 9.;
    static double c10 = 10.;
    static double c11 = 11.;
    static double c12 = 12.;
    static double c13 = 13.;
    static double c14 = 14.;
    static double c16 = 16.;
    static double c18 = 18.;
    static double c19 = 19.;
    static double c20 = 20.;
    static double c22 = 22.;
    static double c24 = 24.;
    static double c25 = 25.;
    static double c26 = 26.;
    static double c27 = 27.;
    static double c28 = 28.;
    static double c32 = 32.;
    static double c35 = 35.;

    /* System generated locals */
    double d__1;

    /* Local variables */
    static double n, t2, t3, t5, t7, t8, t9, t10, t11, t12, t13, t14, t15,
	     t17, t19, t20, t22, t23, t24, t27, t28, t30, t32, t33, t36, t37, 
	    t38, t39, t41, t43, t47, t49, t50, t51, t52, t54, t55, t57, t59, 
	    t62, t63, t65, t69, t70, t73, t74, t76, t79, t81, t82, t84, t89, 
	    t91, t92, t96, t97, t98, t100, t102, t104, t105, t107, t108, t110,
	     t113, t116, t117, t118, t119, t120, t122, t123, t126, t132, t133,
	     t135, t136, t137, t138, t139, t141, t142, t143, t144, t145, t149,
	     t151, t152, t155, t156, t161, t162, t163, t164, t165, t167, t172,
	     t175, t179, t181, t182, t184, t186, t192, t198, t199, t200, t201,
	     t207, t210, t211, t216, t217, t218, t220, t222, t229, t230, t231,
	     t232, t233, t234, t236, t241, t242, t243, t245, t251, t256, t257,
	     t258, t267, t268, t269, t270, t271, t272, t276, t277, t278, t279,
	     t281, t286, t287, t288, t292, t293, t294, t295, t296, t297, t299,
	     t304, t305, t308, t309, t311, t315, t317, t319, t320, t321, t322,
	     t323, t324, t328, t329, t330, t336, t337, t340, t341, t343, t350,
	     t357, t358, t363, t364, t367, t368, t378, t383, t384, t387, t388,
	     t390, t397, t410, t413, t414, t422, t425, t426, t436, t440, t449,
	     t456, t458, t468, t469, t477, t479, t485, t489, t492, t500, t507,
	     t512, t517, t521, t526, t527, t528, t531, t532, t533, t537, t538,
	     t540, t541, t544, t546, t548, t550, t551, t552, t553, t554, t555,
	     t557, t559, t562, t563, t564, t567, t568, t572, t574, t578, t580,
	     t582, t583, t598, t608, t610;


/*          Compute derivatives of an and bn with respect to p and/or q */
/*          when n > 1 */

    /* Parameter adjustments */
    --bn;
    --an;

    /* Function Body */
    n = (double) (*in);
/* Computing 2nd power */
    d__1 = *f;
    t2 = d__1 * d__1;
    t3 = c2 * n - c2;
    t5 = *p * *q;
    t7 = c1 / (t3 * *q + t5);
    t8 = t2 * t7;
/* Computing 2nd power */
    d__1 = n;
    t9 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = t9;
    t10 = d__1 * d__1;
    t11 = t2 * t10;
    t12 = c4 * n - c2;
/* Computing 2nd power */
    d__1 = *q;
    t13 = d__1 * d__1;
    t14 = t12 * t13;
    t15 = *p * t13;
    t17 = c1 / (t14 + c2 * t15);
    t19 = t9 * n;
    t20 = t19 * t2;
    t22 = c1 / (*p + c2 * n - c1);
    t23 = t20 * t22;
    t24 = c2 * n - c1;
    t27 = c1 / (t24 * *q + t5);
    t28 = t20 * t27;
    t30 = t10 * n * t2;
    t32 = n * t2;
    t33 = c2 * n - c3;
    t36 = c1 / (t33 * t13 + t15);
    t37 = t32 * t36;
    t38 = t9 * t2;
    t39 = c1 / t13;
    t41 = t32 * t39;
    t43 = (-c8 + c4 * n) * n;
    t47 = c1 / (c4 + t43 + (c4 * n - c4 + *p) * *p);
    t49 = t38 * t17;
    t50 = t38 * t47;
    t51 = t20 * t47;
    t52 = c1 / *q;
    t54 = t2 * t47;
    t55 = t32 * t47;
    t57 = c1 / (c2 * *p + c4 * n - c6);
    t59 = c4 * t8 - c3 * t11 * t17 - c4 * t23 - t28 - c4 * t30 * t27 + c9 * 
	    t37 - t38 * t39 + t41 + c4 * t11 * t47 - t49 + c24 * t50 - c16 * 
	    t51 - t2 * t52 + c4 * t54 - c16 * t55 - c53 * t38 * t57;
    t62 = c1 / (*p + c2 * n - c2);
    t63 = t32 * t62;
    t65 = c1 / (c2 * *p + c4 * n - c2);
    t69 = t2 / (*p + c2 * n - c3);
    t70 = t69 * t19;
    t73 = c1 / (t3 * t13 + t15);
    t74 = t11 * t73;
    t76 = t10 * t9 * t2;
    t79 = c1 / (t24 * t13 + t15);
    t81 = t2 * t62;
    t82 = c4 + t43;
    t84 = c4 * n - c4;
    t89 = c1 / (t82 * t13 + (t84 * t13 + t15) * *p);
    t91 = t20 * t36;
    t92 = t11 * t27;
    t96 = t20 * t89;
    t97 = t20 * t7;
    t98 = t12 * *q;
    t100 = c1 / (t98 + c2 * t5);
    t102 = c51 * t32 * t57 - c24 * t63 + c5 * t38 * t65 + c12 * t70 + c40 * 
	    t74 + c2 * t76 * t79 + c8 * t81 + c4 * t76 * t89 + c52 * t91 + c6 
	    * t92 - c2 * t69 * t10 - c8 * t20 * t62 + c2 * t11 * t22 - c16 * 
	    t96 - c64 * t97 + t32 * t100;
    t104 = t38 * t62;
    t105 = t30 * t36;
    t107 = c4 * n - c6;
    t108 = t107 * *q;
    t110 = c1 / (t108 + c2 * t5);
    t113 = t38 * t73;
    t116 = c1 / (t33 * *q + t5);
    t117 = t11 * t116;
    t118 = t20 * t116;
    t119 = t30 * t79;
    t120 = t32 * t73;
    t122 = t20 * t73;
    t123 = t20 * t79;
    t126 = c24 * t104 + c14 * t105 + t32 * t52 + c87 * t32 * t110 - c9 * t69 
	    - c12 * t30 * t73 + c24 * t113 - c26 * t117 + c65 * t118 - c2 * 
	    t119 - c4 * t120 + c4 * t30 * t116 - c48 * t122 + c2 * t123 - c2 *
	     t76 * t36 - c3 * t38 * t100;
    t132 = c1 / (t82 * *q + (t84 * *q + t5) * *p);
    t133 = t20 * t132;
    t135 = t38 * t89;
    t136 = t11 * t89;
    t137 = t30 * t89;
    t138 = t11 * t132;
    t139 = t107 * t13;
    t141 = c1 / (t139 + c2 * t15);
    t142 = t38 * t141;
    t143 = t32 * t132;
    t144 = t32 * t7;
    t145 = t38 * t7;
    t149 = t38 * t132;
    t151 = t2 * t116;
    t152 = -c48 * t133 - c8 * t30 * t132 + c4 * t135 + c24 * t136 - c16 * 
	    t137 + c32 * t138 - c69 * t142 - c8 * t143 - c32 * t144 + c72 * 
	    t145 - t32 * t65 + c20 * t11 * t7 - c77 * t11 * t141 + c32 * t149 
	    - c155 * t38 * t110 - c9 * t151;
    t155 = t84 * n;
    t156 = c1 + t155;
    t161 = c1 / (t156 * t13 + (t14 + t15) * *p);
    t162 = t30 * t161;
    t163 = -c8 + c8 * n;
    t164 = t163 * n;
    t165 = c2 + t164;
    t167 = -c4 + c8 * n;
    t172 = c1 / (t165 * t13 + (t167 * t13 + c2 * t15) * *p);
    t175 = (-c24 + c8 * n) * n;
    t179 = c1 / (c18 + t175 + (-c12 + c8 * n + c2 * *p) * *p);
    t181 = t20 * t161;
    t182 = t38 * t22;
    t184 = (c24 + t175) * n;
    t186 = (-c24 + c12 * n) * n;
    t192 = c1 / (-c8 + t184 + (c12 + t186 + (-c6 + c6 * n + *p) * *p) * *p);
    t198 = c1 / (t156 * *q + (t98 + t5) * *p);
    t199 = t11 * t198;
    t200 = t20 * t192;
    t201 = -c4 * t8 + c2 * t162 + c3 * t11 * t172 - c51 * t32 * t179 + c2 * 
	    t23 + c4 * t28 - c2 * t181 - c3 * t182 - c8 * t11 * t192 - c6 * 
	    t199 + c32 * t200 - c6 * t37;
    t207 = c1 / (t165 * *q + (t167 * *q + c2 * t5) * *p);
    t210 = (-c12 + c4 * n) * n;
    t211 = c9 + t210;
    t216 = c1 / (t211 * t13 + (t139 + t15) * *p);
    t217 = t32 * t216;
    t218 = -c8 + t184;
    t220 = c12 + t186;
    t222 = -c6 + c6 * n;
    t229 = c1 / (t218 * t13 + (t220 * t13 + (t222 * t13 + t15) * *p) * *p);
    t230 = t11 * t229;
    t231 = t20 * t216;
    t232 = t69 * n;
    t233 = t30 * t216;
    t234 = c18 + t175;
    t236 = -c12 + c8 * n;
    t241 = c1 / (t234 * t13 + (t236 * t13 + c2 * t15) * *p);
    t242 = t38 * t241;
    t243 = c3 * t38 * t207 - c36 * t50 + c12 * t51 - c12 * t54 - c9 * t217 + 
	    c36 * t55 + c12 * t63 - c48 * t230 - c52 * t231 - c13 * t232 - 
	    c14 * t233 + c69 * t242;
    t245 = t32 * t192;
    t251 = c1 / (t234 * *q + (t236 * *q + c2 * t5) * *p);
    t256 = c1 / (c1 + t155 + (c4 * n - c2 + *p) * *p);
    t257 = t20 * t256;
    t258 = c32 * t245 - c2 * t70 - c10 * t74 - c6 * t81 - c22 * t91 - c4 * 
	    t92 + c60 * t96 + c16 * t97 - c6 * t104 - c87 * t32 * t251 - c2 * 
	    t105 + c4 * t257;
    t267 = c1 / (t218 * *q + (t220 * *q + (t222 * *q + t5) * *p) * *p);
    t268 = t11 * t267;
    t269 = t11 * t79;
    t270 = t30 * t229;
    t271 = t32 * t267;
    t272 = c6 * t69 - c64 * t268 - c18 * t113 + c4 * t117 - c20 * t118 - t269 
	    + c32 * t270 + c2 * t119 + c4 * t120 + c24 * t122 - c2 * t123 + 
	    c16 * t271;
    t276 = t32 * t27;
    t277 = t69 * t9;
    t278 = t38 * t116;
    t279 = t38 * t192;
    t281 = c77 * t11 * t241 - t276 + c88 * t133 - c28 * t135 - c52 * t136 + 
	    c16 * t137 + c9 * t277 + c35 * t278 - c28 * t138 - c48 * t279 + 
	    c40 * t143 + c155 * t38 * t251;
    t286 = c1 / (t211 * *q + (t108 + t5) * *p);
    t287 = t20 * t286;
    t288 = t2 * t192;
    t292 = c1 / (c9 + t210 + (c4 * n - c6 + *p) * *p);
    t293 = t2 * t292;
    t294 = t2 * t286;
    t295 = t20 * t267;
    t296 = t2 * t132;
    t297 = t32 * t89;
    t299 = c24 * t144 - c36 * t145 - c96 * t149 - c65 * t287 + c6 * t151 - c8 
	    * t288 + c9 * t293 + c9 * t294 + c96 * t295 - c4 * t296 + c4 * 
	    t297 - c4 * t30 * t286;
    t304 = t11 * t286;
    t305 = t32 * t116;
    t308 = t38 * t267;
    t309 = t11 * t36;
    t311 = t38 * t79;
    t315 = c1 / (c2 + t164 + (-c4 + c8 * n + c2 * *p) * *p);
    t317 = c2 * t11 * t292 - t32 * t207 - c2 * t11 * t256 + c26 * t304 - c25 *
	     t305 + c4 * t30 * t198 + c16 * t30 * t267 - c64 * t308 + c11 * 
	    t309 - c8 * t76 * t229 + t311 - c5 * t38 * t315;
    t319 = t32 * t22;
    t320 = t20 * t198;
    t321 = t20 * t292;
    t322 = t38 * t229;
    t323 = t38 * t27;
    t324 = t20 * t229;
    t328 = t38 * t36;
    t329 = t38 * t172;
    t330 = t32 * t315 + t319 + t320 - c12 * t321 - c8 * t322 + t323 + c32 * 
	    t324 - c2 * t76 * t161 + c2 * t76 * t216 + c53 * t38 * t179 + c19 
	    * t328 + t329;
    t336 = (c6 + t236 * n) * n;
    t337 = -c1 + t336;
    t340 = (-c12 + c12 * n) * n;
    t341 = c3 + t340;
    t343 = -c3 + c6 * n;
    t350 = c1 / (t337 * t13 + (t341 * t13 + (t343 * t13 + t15) * *p) * *p);
    t357 = (-c64 + (c96 + (-c64 + c16 * n) * n) * n) * n;
    t358 = c16 + t357;
    t363 = (c96 + (-c96 + c32 * n) * n) * n;
    t364 = -c32 + t363;
    t367 = (-c48 + c24 * n) * n;
    t368 = c24 + t367;
    t378 = c1 / (t358 * *q + (t364 * *q + (t368 * *q + (t163 * *q + t5) * *p) 
	    * *p) * *p);
    t383 = (c54 + (-c36 + c8 * n) * n) * n;
    t384 = -c27 + t383;
    t387 = (-c36 + c12 * n) * n;
    t388 = c27 + t387;
    t390 = -c9 + c6 * n;
    t397 = c1 / (t384 * *q + (t388 * *q + (t390 * *q + t5) * *p) * *p);
    t410 = c1 / (t358 * t13 + (t364 * t13 + (t368 * t13 + (t163 * t13 + t15) *
	     *p) * *p) * *p);
    t413 = t32 * t286;
    t414 = -c3 * t11 * t350 + c2 * t8 - c4 * t162 - c2 * t28 + c4 * t181 + 
	    t182 + c8 * t199 - c288 * t20 * t378 - c32 * t200 + c2 * t37 + c8 
	    * t30 * t397 + c24 * t38 * t410 + c4 * t76 * t350 + c50 * t413 + 
	    c14 * t50;
    t422 = c1 / (c16 + t357 + (-c32 + t363 + (c24 + t367 + (-c8 + c8 * n + *p)
	     * *p) * *p) * *p);
    t425 = t32 * t229;
    t426 = -c96 * t20 * t422 + c14 * t54 + c12 * t217 - c28 * t55 - c96 * t30 
	    * t410 - c2 * t63 + c128 * t230 + c44 * t231 + c3 * t232 + c4 * 
	    t233 - c96 * t245 - c8 * t425 + c2 * t81 + c4 * t91 - c52 * t96;
    t436 = c1 / (t337 * *q + (t341 * *q + (t343 * *q + t5) * *p) * *p);
    t440 = c12 * t11 * t436 - c4 * t257 - c2 * t69 + c72 * t268 + c6 * t113 - 
	    c18 * t38 * t292 + c2 * t118 + t269 + c144 * t11 * t410 - c40 * 
	    t270 - c2 * t120 - c4 * t122 - c96 * t271 + t276 - c36 * t133;
    t449 = c1 / (t384 * t13 + (t388 * t13 + (t390 * t13 + t15) * *p) * *p);
    t456 = c1 / (-c1 + t336 + (c3 + t340 + (-c3 + c6 * n + *p) * *p) * *p);
    t458 = c38 * t135 + c22 * t136 - t277 - c69 * t38 * t449 - c7 * t278 + 
	    c96 * t279 - c52 * t143 - c8 * t144 + c6 * t145 + c80 * t149 + 
	    c40 * t287 - c2 * t151 + c32 * t288 - c12 * t293 + c4 * t11 * 
	    t456 - c12 * t294;
    t468 = -c224 * t295 + c8 * t296 - c8 * t297 + c24 * t76 * t410 - c8 * 
	    t304 - c52 * t11 * t397 + c7 * t305 + c87 * t32 * t397 + c240 * 
	    t308 - t309 - t311 + c104 * t20 * t449 - c8 * t20 * t456 + c5 * 
	    t38 * t456 + t32 * t436;
    t469 = t38 * t198;
    t477 = -c2 * t469 + c26 * t32 * t292 - t319 - c8 * t320 + c4 * t321 + c64 
	    * t322 + t323 - c144 * t324 - c5 * t328 - c18 * t2 * t397 + c24 * 
	    t2 * t422 + c130 * t20 * t397 - c155 * t38 * t397 - c96 * t32 * 
	    t422 - c96 * t20 * t410;
    t479 = t11 * t161;
    t485 = c1 / (-c27 + t383 + (c27 + t387 + (c6 * n - c9 + *p) * *p) * *p);
    t489 = t32 * t198;
    t492 = t2 * t267;
    t500 = c2 * t479 + c51 * t32 * t485 - c77 * t11 * t449 + c28 * t30 * t449 
	    + c2 * t489 + c18 * t32 * t449 - c2 * t38 * t161 + c8 * t492 - 
	    t38 * t350 + c4 * t20 * t350 - c8 * t30 * t436 - c4 * t30 * t350 
	    + c6 * t38 * t256 - c3 * t38 * t436 + c24 * t11 * t422;
    t507 = t38 * t286;
    t512 = t11 * t216;
    t517 = -c2 * t32 * t256 - c4 * t11 * t485 - c53 * t38 * t485 + c144 * t38 
	    * t422 + c24 * t20 * t485 - c18 * t2 * t485 - c70 * t507 - t32 * 
	    t456 - c48 * t32 * t378 - c48 * t30 * t378 + c192 * t38 * t378 - 
	    c22 * t512 + c192 * t11 * t378 - c38 * t38 * t216 - c2 * t20 * 
	    t436 - c4 * t76 * t449;
    t521 = c16 * t8 - c8 * t28 + t41 - c3 * t49 + c20 * t74 + c65 * t91 + c4 *
	     t92 - c48 * t96 - c16 * t97 + c4 * t105 + c72 * t113 - c4 * t117 
	    + c24 * t118 + c6 * t269 - c4 * t119 - c32 * t120 - c64 * t122 - 
	    t123 - t276 - c32 * t133;
    t526 = t2 * t73;
    t527 = t2 * t36;
    t528 = c48 * t149 - c18 * t151 + c8 * t296 - c8 * t297 + c51 * t305 - c26 
	    * t309 + c5 * t323 + t32 * t17 + c87 * t32 * t141 + c4 * t526 - 
	    c9 * t527;
    t531 = t2 * t89;
    t532 = t32 * t79;
    t533 = -c32 * t96 + c8 * t136 + c48 * t135 - c48 * t120 + c24 * t91 + c48 
	    * t113 - c16 * t122 - c53 * t328 - c18 * t527 - c8 * t123 + c16 * 
	    t526 + c8 * t531 + c5 * t311 - t532 + c51 * t37 - c4 * t309 + c4 *
	     t269 - c32 * t297;
    t537 = c9 * t2 * t216 - c87 * t32 * t241 - t32 * t172 - c12 * t8 + c4 * 
	    t162 + c4 * t28 + t181 - c4 * t199 - c25 * t37 - c51 * t413 - c64 
	    * t230 - c65 * t231 - c4 * t233 + c155 * t242 + c16 * t425;
    t538 = -c20 * t91 + c88 * t96 - c16 * t268 - c36 * t113 - c4 * t118 - c4 *
	     t269 + c16 * t270 + c24 * t120 + c16 * t122 + c4 * t123 + c64 * 
	    t271 + c2 * t276 + c24 * t133 - c96 * t135 - c28 * t136 + c18 * 
	    t278;
    t540 = c72 * t143 + c24 * t144 - c12 * t145 - c72 * t149 - c24 * t287 + 
	    c12 * t151 + c18 * t294 + c64 * t295 - c24 * t296 + c40 * t297 + 
	    c4 * t304 - c26 * t305 - c96 * t308 + c4 * t309 + t311;
    t541 = -c5 * t469 + c8 * t320 - c64 * t322 - c6 * t323 + c96 * t324 + c35 
	    * t328 + c3 * t329 - c6 * t479 + t489 - c16 * t492 + c53 * t507 + 
	    c26 * t512 - c4 * t526 + c6 * t527 - t532 - c4 * t531;
    t544 = t9 * *f;
    t546 = c1 / (*p + c2 * n);
    t548 = *q * n;
    t550 = c1 / (t5 + c2 * t548);
    t551 = t544 * t550;
    t552 = t544 * t7;
    t553 = n * *f;
    t554 = t553 * t7;
    t555 = t19 * *f;
    t557 = *f * t62;
    t559 = t557 * n;
    t562 = c1 - *f + c2 * t544 * t546 - c2 * t551 - c4 * t552 + c2 * t554 + 
	    c2 * t555 * t7 - c2 * t557 - c2 * t557 * t9 + c4 * t559 - c2 * 
	    t555 * t550 + c2 * t553 * t52;
    t563 = t553 * t550;
    t564 = t553 * t132;
    t567 = t544 * t132;
    t568 = *f * t47;
    t572 = c1 / (c4 * t9 + (c4 * n + *p) * *p);
    t574 = *q * t9;
    t578 = c1 / (c4 * t574 + (c4 * t548 + t5) * *p);
    t580 = t544 * t578;
    t582 = t553 * t47;
    t583 = -t563 - c2 * t564 + c2 * t544 * t47 - c2 * t555 * t132 + c4 * t567 
	    + c2 * t568 - c2 * t544 * t572 + c2 * t555 * t578 - t551 + c2 * 
	    t580 + t552 - t554 + t557 - t559 + t553 * t546 - c4 * t582;
    t598 = c1 / (c8 * *q * t19 + (c12 * t574 + (c6 * t548 + t5) * *p) * *p);
    t608 = c2 * t564 - c2 * t567 - c2 * t568 + c2 * t580 + c2 * t553 * t578 + 
	    c4 * t544 / (c8 * t19 + (c12 * t9 + (c6 * n + *p) * *p) * *p) - 
	    c4 * t555 * t598 - c2 * t553 * t572 + c8 * t553 * t192 - c4 * *f *
	     t192 - c4 * t544 * t192 + c4 * t553 * t267 - c8 * t544 * t267 + 
	    c4 * t555 * t267 - c4 * t544 * t598 + c2 * t582;
    t610 = *f * t7;
    an[1] = t59 + t102 + t126 + t152;
    an[2] = t201 + t243 + t258 + t272 + t281 + t299 + t317 + t330;
    an[3] = t414 + t426 + t440 + t458 + t468 + t477 + t500 + t517;
    an[4] = t521 + c32 * t135 + c32 * t136 - c8 * t137 - t2 * t39 - c53 * 
	    t278 + c8 * t138 - c155 * t142 - c32 * t143 - c48 * t144 + c48 * 
	    t145 + t528;
    an[5] = t533;
    an[6] = t537 + t538 + t540 + t541;
    bn[1] = t562;
    bn[2] = t583;
    bn[3] = t608;
    bn[4] = -(*f) * t52 - c2 * t552 + c4 * t554 - c2 * t610 + c2 * t551;
    bn[5] = c0;
    bn[6] = c2 * t567 - t554 - c4 * t564 + t563 - c2 * t580 + t610 + c2 * *f *
	     t132;
    return 0;
} /* subd_ */
