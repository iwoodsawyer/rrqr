/* ztrrnk.f -- translated by f2c (version 20041007).
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

/* Subroutine */ int ztrrnk(integer *n, doublecomplex *r__, integer *ldr,
	doublereal *rcond, integer *rank, doublecomplex *work, integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublecomplex c1, c2, s1, s2;
    static doublereal smin, smax;
    static doublereal sminpr, smaxpr;


/*     $Revision: 1.42 $ */
/*     $Date: 96/12/30 16:59:45 $ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */

/*  Purpose */
/*  ======= */

/*  ZTRRNK computes an estimate for the numerical rank of a */
/*  triangular n-by-n matrix R. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          Number of rows and columns of the matrix R.  N >= 0. */

/*  R       (input) COMPLEX*16 array, dimension (LDR,N) */
/*          On entry, the n by n matrix R. */

/*  LDR     (input) INTEGER */
/*          The leading dimension of the array R. LDR >= max(1,N). */

/*  RCOND   (input) DOUBLE PRECISION */
/*          Threshold value for the numerical rank. */

/*  RANK    (output) INTEGER */
/*          Numerical rank for threshold RCOND. */

/*  WORK    (workspace) COMPLEX*16 array, dimension (2*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  Successful exit. */
/*          < 0:  If INFO = -i, the i-th argument had an illegal value. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*ldr < max(1,*n)) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("ZTRRNK", &i__1, 6);
	return 0;
    }

/*     Determine RANK using incremental condition estimation. */

    work[1].r = 1., work[1].i = 0.;
    i__1 = *n + 1;
    work[i__1].r = 1., work[i__1].i = 0.;
    smax = z_abs(&r__[r_dim1 + 1]);
    smin = smax;
    if (z_abs(&r__[r_dim1 + 1]) == 0.) {
	*rank = 0;
	goto L30;
    } else {
	*rank = 1;
    }

L10:
    if (*rank < *n) {
	i__ = *rank + 1;
	zlaic1(&c__2, rank, &work[1], &smin, &r__[i__ * r_dim1 + 1], &r__[
		i__ + i__ * r_dim1], &sminpr, &s1, &c1);
	zlaic1(&c__1, rank, &work[*n + 1], &smax, &r__[i__ * r_dim1 + 1], &
		r__[i__ + i__ * r_dim1], &smaxpr, &s2, &c2);

	if (smaxpr * *rcond <= sminpr) {
	    i__1 = *rank;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		i__3 = i__;
		z__1.r = s1.r * work[i__3].r - s1.i * work[i__3].i, z__1.i =
			s1.r * work[i__3].i + s1.i * work[i__3].r;
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
		i__2 = *n + i__;
		i__3 = *n + i__;
		z__1.r = s2.r * work[i__3].r - s2.i * work[i__3].i, z__1.i =
			s2.r * work[i__3].i + s2.i * work[i__3].r;
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
/* L20: */
	    }
	    i__1 = *rank + 1;
	    work[i__1].r = c1.r, work[i__1].i = c1.i;
	    i__1 = *n + *rank + 1;
	    work[i__1].r = c2.r, work[i__1].i = c2.i;
	    smin = sminpr;
	    smax = smaxpr;
	    ++(*rank);
	    goto L10;
	}
    }
L30:

    return 0;

/*     End of ZTRRNK */

} /* ztrrnk_ */

