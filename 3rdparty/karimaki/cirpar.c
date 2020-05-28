/* cirpar.f -- translated by f2c (version 20191129).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* CMZ :  8.14/00 01/07/98  20.38.08  by  Stephen Burke */
/* CMZU:  8.13/00 24/04/98  16.03.21  by  C.Niebuhr */
/* CMZ :  8.04/00 27/06/96  20.27.36  by  Stephen Burke */
/* CMZ :  7.02/00 12/09/95  16.37.24  by  Gaby Raedel */
/* -- Author :  Volker Blobel */

/* The descriptions inserted below will appear asis with the */
/* initial asterisk removed in the WWW writeup. HTML commands */
/* (links etc.) may be included in the description if desired. */
/* The regions can be expanded to as many comment lines as */
/* required up to a maximum of 100 from *HTMLP to *HTMLE inclusive. */

/* HTMLP     : Describe the Purpose of the routine */



/* HTMLI     : Describe the Input variables to the routine */



/* HTMLO     : Describe the Output of the routine */



/* HTMLE     : Terminates the HTML documentation */

/* Subroutine */ int cirpar_(real *x, real *y, integer *n, real *tr, real *
	cdf)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    static doublereal equiv_9[10];

    /* Builtin functions */
    double atan2(doublereal, doublereal), sin(doublereal), cos(doublereal), 
	    sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal u, q1, q2, r2, cx, cy;
    static real xi, yi;
#define sw (equiv_9)
#define sx (equiv_9 + 1)
#define sy (equiv_9 + 2)
    static doublereal cr2, cr4;
#define sr2 (equiv_9 + 8)
#define sr4 (equiv_9 + 9)
    static doublereal dca, fka;
    static real dln;
    static doublereal phi, dlt;
    static real spa, spb, spc, dtr;
    static doublereal cxx, cxy, cyy, rnv;
    static real spt;
    static doublereal sqt;
#define sxx (equiv_9 + 3)
#define sxy (equiv_9 + 4)
#define syy (equiv_9 + 5)
    static doublereal cxr2, cyr2;
#define sxr2 (equiv_9 + 6)
#define syr2 (equiv_9 + 7)
#define zero (equiv_9)
    static doublereal cosphi, sinphi;


/*     Circle parameter from N points (X,Y) */

/*     Input: */
/*            X(.), Y(.) = arrays of N data points in a plane */
/*            W(.)       = arrays of N weights (=1/sigma**2) */

/*     Result: */
/*     TR(1) =RNV = inverse radius */
/*     TR(2) =DCA = distance of closest approach to 0.0, 0.0 */
/*     TR(3) =PHI = phi angle, all in H1 convention */
/*            CDF = chi square, divided by (N-3) */

/* -- Author :  Volker Blobel, Algorithm by Veikko Karimaeki */
/*     Veikko Karimaeki, Fast Code to fit circular Arcs, University of */
/*     Helsinki Report HU-SEFT-1991-10 */
/*     ... */
    /* Parameter adjustments */
    --tr;
    --y;
    --x;

    /* Function Body */
    *cdf = 0.f;
/*     ------------------------------------------------------------------ */
/*     fit to shifted reference point (to reduce round-off errors) */
/*     ------------------------------------------------------------------ */
    if (*n < 3) {
	goto L100;
    }
/*     zero all sums */
    for (i__ = 1; i__ <= 10; ++i__) {
/* L10: */
	zero[i__ - 1] = 0.;
    }
/*     shift to first point */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*sw += 1.f;
	*sx += (doublereal) x[i__] - (doublereal) x[1];
/* L20: */
	*sy += (doublereal) y[i__] - (doublereal) y[1];
    }
    cx = *sx / *sw;
    cy = *sy / *sw;
/*     form sums ... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*     ... using shifted coordinates */
	xi = (doublereal) x[i__] - (doublereal) x[1];
	yi = (doublereal) y[i__] - (doublereal) y[1];
	r2 = xi * xi + yi * yi;
	*sxx += (xi - cx) * (xi - cx);
	*sxy += (xi - cx) * (yi - cy);
	*syy += (yi - cy) * (yi - cy);
	*sxr2 += xi * r2;
	*syr2 += yi * r2;
	*sr2 += r2;
/* L30: */
	*sr4 += r2 * r2;
    }
/*     calculate averages and moments */
    cxx = *sxx / *sw;
    cxy = *sxy / *sw;
    cyy = *syy / *sw;
    *sxx += *sw * cx * cx;
    *sxy += *sw * cx * cy;
    *syy += *sw * cy * cy;
    cxr2 = (*sxr2 - cx * *sr2) / *sw;
    cyr2 = (*syr2 - cy * *sr2) / *sw;
    cr2 = *sr2 / *sw;
    cr4 = (*sr4 - *sr2 * cr2) / *sw;
/*     calculate circle parameter */
    q1 = cr4 * cxy - cxr2 * cyr2;
    q2 = cr4 * (cxx - cyy) - cxr2 * cxr2 + cyr2 * cyr2;
    phi = atan2(q1 * 2., q2) * .5;
    sinphi = sin(phi);
    cosphi = cos(phi);
/*     compare PHI with initial track direction */
/*     IF(COSPHI*CX+SINPHI*CY.LT.0.0) THEN */
    if (cosphi * ((doublereal) x[1] + cx) + sinphi * ((doublereal) y[1] + cy) 
	    < 0.f) {
/*         reverse direction */
	if (phi < 0.f) {
	    phi += 3.141592654f;
	} else {
	    phi += -3.141592654f;
	}
	cosphi = -cosphi;
	sinphi = -sinphi;
    }
    fka = (sinphi * cxr2 - cosphi * cyr2) / cr4;
    dlt = -fka * cr2 + sinphi * cx - cosphi * cy;
    sqt = sqrt(1. - dlt * 4. * fka);
    rnv = fka * 2. / sqt;
    dca = dlt * 2. / (sqt + 1.);
    u = rnv * dca + 1.;
/*     ------------------------------------------------------------------ */
/*     calculate estimate for standard deviation */
/*     ------------------------------------------------------------------ */
    if (*n > 3) {
/*     estimate for chi square/(n-3) (note: fit is to unweighted data!).. */
/* Computing 2nd power */
	d__1 = rnv * dca + 1.;
	*cdf = *sw * (d__1 * d__1) * (sinphi * sinphi * cxx - sinphi * 2. * 
		cosphi * cxy + cosphi * cosphi * cyy - fka * fka * cr4) / (
		doublereal) (*n - 3);
    }
/*     ------------------------------------------------------------------ */
/*     shift to reference point (0,0) */
/*     ------------------------------------------------------------------ */
    dtr = (doublereal) x[1] * sinphi - (doublereal) y[1] * cosphi + dca;
    dln = (doublereal) x[1] * cosphi + (doublereal) y[1] * sinphi;
    spa = dtr * 2.f + rnv * (dtr * dtr + dln * dln);
    spb = rnv * (doublereal) x[1] + u * sinphi;
    spc = -rnv * (doublereal) y[1] + u * cosphi;
    spt = sqrt(rnv * spa + 1.f);
/*     shifted parameters */
    tr[1] = rnv;
    tr[2] = spa / (spt + 1.f);
    tr[3] = atan2(spb, spc);
/*     ------------------------------------------------------------------ */
/*     change sign of RNV */
/*     ------------------------------------------------------------------ */
    tr[1] = -tr[1];
L100:
    return 0;
} /* cirpar_ */

#undef zero
#undef syr2
#undef sxr2
#undef syy
#undef sxy
#undef sxx
#undef sr4
#undef sr2
#undef sy
#undef sx
#undef sw


