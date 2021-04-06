/* slauc1.f -- translated by f2c (version 20041007).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "rrqr.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

logical slauc1(integer *k, real *x, real *smin, real *w, real *gamma, real *
	thresh)
{
    /* System generated locals */
    integer i__1;
    logical ret_val;

    /* Local variables */
    static real sine;
    static real cosine, sminpr;


/*     This code is part of a release of the package for computing */
/*     rank-revealing QR Factorizations written by: */
/*     ================================================================== */
/*     Christian H. Bischof        and   Gregorio Quintana-Orti */
/*     Math. and Comp. Sci. Div.         Departamento de Informatica */
/*     Argonne National Lab.             Universidad Jaime I */
/*     Argonne, IL 60439                 Campus P. Roja, 12071 Castellon */
/*     USA                               Spain */
/*     bischof@mcs.anl.gov               gquintan@inf.uji.es */
/*     ================================================================== */
/*     $Revision: 1.84 $ */
/*     $Date: 96/12/30 16:59:16 $ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLAUC1 applies incremental condition estimation to determine whether */
/*  the K-th column of A, stored in vector W, would be acceptable as a pivot */
/*  column with respect to the threshold THRESH. */


/*  Arguments */
/*  ========= */

/*  K       (input) INTEGER */
/*          Length of vector X. */

/*  X       (input/output) REAL array, dimension ( K ) */
/*          On entry, X(1:K-1) contains an approximate smallest left */
/*          singular vector of the upper triangle of A(1:k-1,1:k-1). */
/*          On exit, if SLAUC1 == .TRUE., X contains an approximate */
/*          smallest left singular vector of the upper triangle of */
/*          A(1:k,1:k); if SLAUC1 == .FALSE., X is unchanged. */

/*  SMIN    (input/output) REAL */
/*          On entry, an estimate for the smallest singular value of the */
/*          upper triangle of A(1:k-1,1:k-1). */
/*          On exit, if SLAUC1 == .TRUE., SMIN is an estimate of the */
/*          smallest singular value of the upper triangle of  A(1:k,1:k); */
/*          if SLAUC1 == .FALSE., SMIN is unchanged. */

/*  W       (input) REAL array, dimension ( K-1 ) */
/*          The K-th column of matrix A excluding the diagonal element. */

/*  GAMMA   (input) REAL */
/*          Diagonal entry in k-th column of A if column k were to */
/*          be accepted. */

/*  THRESH  (input) REAL */
/*          If the approximate smallest singular value for A(1:K,1:K) */
/*          is smaller than THRESH, the kth column is rejected. */

/*  (SLAUC1) (output) LOGICAL */
/*          If the k-th column of A is found acceptable, SLAUC1 */
/*          returns .TRUE., otherwise it returns .FALSE. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */

/*     .. */
/*     .. Executable Statements .. */


/*     Try to use diagonal element as condition estimator */

    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
    if (*thresh > dabs(*gamma)) {
	ret_val = FALSE_;
	return ret_val;
    }

/*     Use incremental condition estimation to determine an estimate */
/*     SMINPR and an approximate singular vector [SINE*X,COSINE]' */
/*     for A(K,K). */

    i__1 = *k - 1;
    slaic1(&c__2, &i__1, &x[1], smin, &w[1], gamma, &sminpr, &sine, &cosine);
    if (*thresh > sminpr) {
	ret_val = FALSE_;
    } else {
	i__1 = *k - 1;
	sscal(&i__1, &sine, &x[1], &c__1);
	x[*k] = cosine;
	*smin = sminpr;
	ret_val = TRUE_;
    }
    return ret_val;

/*     End of SLAUC1 */

} /* slauc1_ */

