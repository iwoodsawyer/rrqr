/* dgeqpw.f -- translated by f2c (version 20041007).
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

static integer c__1 = 1;

/* Subroutine */ int dgeqpw(integer *m, integer *lwsize, integer *nb,
	integer *offset, integer *lacptd, doublereal *a, integer *lda,
	integer *jpvt, doublereal *ircond, doublereal *x, doublereal *smin,
	doublereal *mxnm, doublereal *tau, doublereal *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__, k, i1;
    static doublereal akk, temp, smax;
    static doublereal temp2, gamma;
    static integer lastk;
    static integer pvtidx;


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
/*     $Date: 96/12/30 16:59:12 $ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */


/*  Purpose */
/*  ======= */

/*  DGEQPW applies one block step of the Householder QR factorization */
/*  algorithm with restricted pivoting. It is called by DGEQPB. */

/*  Let A be the partial QR factorization of an M by (OFFSET+LWSIZE) */
/*  matrix C, i.e. we have computed an orthogonal matrix Q1 and a */
/*  permutation matrix P1 such that */
/*            C * P1 = Q1 * A */
/*  and A(:,1:OFFSET) is upper triangular. Let us denote A(:,1:OFFSET) */
/*  by B. Then in addition let */
/*  X be an approximate smallest left singular vector of B in the sense */
/*  that */
/*       sigma_min(B) ~ twonorm(B'*X) = SMIN */
/*  and */
/*       sigma_max(B) ~ ((offset)**(1./3.))*MXNM = SMAX */
/*  with */
/*       cond_no(B) ~ SMAX/SMIN <= 1/IRCOND */

/*  Then DGEQP2 tries to identify NB columns in */
/*   A(:,OFFSET+1:OFFSET+LWSIZE) such that */
/*       cond_no([B,D]) < 1/IRCOND */
/*  where D are the KB columns of A(:,OFFSET+1:OFFSET+LWSIZE) that were */
/*  considered independent with respect to the threshold 1/IRCOND. */

/*  On exit, */
/*       C * P2 = Q2 * A */
/*  is again a partial QR factorization of C, but columns */
/*  OFFSET+1:OFFSET+LACPTD of A have been reduced via */
/*  a series of elementary reflectors to upper */
/*  trapezoidal form. Further */
/*       sigma_min(A(:,1:OFFSET+LACPTD)) */
/*                ~ twonorm(A(:,1:OFFSET+LACPTD)'*x) = SMIN */
/*  and */
/*       sigma_max(A(:,1:OFFSET+LACPTD)) ~ sqrt(OFFSET+LACPTD)*MXNM = SMAX */
/*  with */
/*       cond_no(A(:,1:OFFSET+LACPTD)) */
/*                   ~ SMAX/SMIN <= 1/IRCOND. */

/*  In the ideal case, LACPTD = NB, that is, */
/*  we found NB independent columns in the window consisting of */
/*  the first LWSIZE columns of A. */


/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A. */

/*  LWSIZE  (input) INTEGER */
/*          The size of the pivot window in A. */

/*  NB      (input) INTEGER */
/*          The number of independent columns one would like to identify. */
/*          This equals the desired blocksize in DGEQPB. */

/*  OFFSET  (input) INTEGER */
/*          The number of rows and columns of A that need not be updated. */

/*  LACPTD  (output) INTEGER */
/*          The number of columns in A(:,OFFSET+LWSIZE) that were */
/*          accepted as linearly independent. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,OFFSET+LWSIZE) */
/*          On entry, the upper triangle of A(:,1:OFFSET) contains the */
/*          partially completed triangular factor R; the elements below */
/*          the diagonal, with the array TAU, represent the orthogonal */
/*          matrix Q1 as a product of elementary reflectors. */
/*          On exit, the upper triangle of A(:,OFFSET+LACPTD) contains */
/*          the partially completed upper triangular factor R; the */
/*          elements below the diagonal, with the array TAU, represent */
/*          the orthogonal matrix Q2 as a product of elementary */
/*          reflectors. */
/*          A(OFFSET:M,LACPTD+1:LWSIZE) has been updated by the product */
/*          of these elementary reflectors. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. LDA >= M. */

/*  JPVT    (input/output) INTEGER array, dimension (OFFSET+LWSIZE) */
/*          On entry and exit, jpvt(i) = k if the i-th column */
/*          of A was the k-th column of C. */

/*  IRCOND  (input) DOUBLE PRECISION */
/*          1/IRCOND is the threshold for the condition number. */

/*  X       (input/output) DOUBLE PRECISION array, dimension (OFFSET+NB) */
/*          On entry, X(1:OFFSET) is an approximate left nullvector of */
/*          the upper triangle of A(1:OFFSET,1:OFFSET). */
/*          On exit, X(1:OFFSET+LACPTD) is an approximate left */
/*          nullvector of the matrix in the upper triangle of */
/*          A(1:OFFSET+LACPTD,1:OFFSET+LACPTD). */

/*  SMIN    (input/output) DOUBLE PRECISION */
/*          On entry, SMIN is an estimate for the smallest singular */
/*          value of the upper triangle of A(1:OFFSET,1:OFFSET). */
/*          On exit, SMIN is an estimate for the smallest singular */
/*          value of the matrix in the upper triangle of */
/*          A(1:OFFSET+LACPTD,1:OFFSET+LACPTD). */

/*  MXNM    (input) DOUBLE PRECISION */
/*          The norm of the largest column in matrix A. */

/*  TAU     (output) DOUBLE PRECISION array, dimension (OFFSET+LWSIZE) */
/*          On exit, TAU(1:OFFSET+LACPTD) contains details of */
/*          the orthogonal matrix Q2. */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*LWSIZE) */

/*  ================================================================ */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */


/*     Initialize partial column norms (stored in the first */
/*     LWSIZE entries of WORK) and exact column norms (stored */
/*     in the second LWSIZE entries of WORK) for the first */
/*     batch of columns. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --x;
    --tau;
    --work;

    /* Function Body */
    i__1 = *lwsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m - *offset;
	work[i__] = dnrm2(&i__2, &a[*offset + 1 + (*offset + i__) * a_dim1],
		&c__1);
	work[*lwsize + i__] = work[i__];
/* L10: */
    }

/*     ************* */
/*     * Main loop * */
/*     ************* */

/* Computing MIN */
    i__1 = *m, i__2 = *offset + *lwsize;
    lastk = min(i__1,i__2);
    *lacptd = 0;
L1000:
    if (*lacptd == *nb) {
	goto L2000;
    }

/*        Determine pivot candidate. */
/*        ========================= */

    i__1 = *lwsize - *lacptd;
    pvtidx = *offset + *lacptd + idamax(&i__1, &work[*lacptd + 1], &c__1);
    k = *offset + *lacptd + 1;

/*        Exchange current column and pivot column. */

    if (pvtidx != k) {
	dswap(m, &a[pvtidx * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
	i1 = jpvt[pvtidx];
	jpvt[pvtidx] = jpvt[k];
	jpvt[k] = i1;
	temp = work[pvtidx - *offset];
	work[pvtidx - *offset] = work[k - *offset];
	work[k - *offset] = temp;
	temp = work[pvtidx - *offset + *lwsize];
	work[pvtidx - *offset + *lwsize] = work[k + *lwsize - *offset];
	work[k + *lwsize - *offset] = temp;
    }

/*        Determine (offset+lacptd+1)st diagonal element GAMMA of */
/*        matrix A if elementary reflector were applied. */

    if (a[k + k * a_dim1] == 0.) {
	gamma = -work[k - *offset];
    } else {
	gamma = -d_sign(&work[k - *offset], &a[k + k * a_dim1]);
    }

/*        Update estimate for largest singular value. */

    smax = dlasmx(&k) * *mxnm;

/*        Is candidate pivot column acceptable ? */
/*        ===================================== */

    d__1 = smax * *ircond;
    if (dlauc1(&k, &x[1], smin, &a[k * a_dim1 + 1], &gamma, &d__1)) {

/*           Pivot candidate was accepted. */
/*           ============================ */

	++(*lacptd);

/*           Generate Householder vector. */

	if (k < *m) {
	    i__1 = *m - k + 1;
	    dlarfg(&i__1, &a[k + k * a_dim1], &a[k + 1 + k * a_dim1], &c__1,
		    &tau[k]);
	} else {
	    dlarfg(&c__1, &a[*m + k * a_dim1], &a[*m + k * a_dim1], &c__1, &
		    tau[k]);
	}

/*           Apply Householder reflection to A(k:m,k+1:lwsize). */

	if (*lacptd < *lwsize) {
	    akk = a[k + k * a_dim1];
	    a[k + k * a_dim1] = 1.;
	    i__1 = *m - k + 1;
	    i__2 = *lwsize - *lacptd;
	    dlarf("Left", &i__1, &i__2, &a[k + k * a_dim1], &c__1, &tau[k], &
		    a[k + (k + 1) * a_dim1], lda, &work[(*lwsize << 1) + 1]);
	    a[k + k * a_dim1] = akk;
	}

/*           Update partial column norms. */

	if (k < lastk) {
	    i__1 = *lwsize;
	    for (i__ = *lacptd + 1; i__ <= i__1; ++i__) {
		if (work[i__] != 0.) {
/* Computing 2nd power */
		    d__2 = (d__1 = a[k + (*offset + i__) * a_dim1], abs(d__1))
			     / work[i__];
		    temp = 1. - d__2 * d__2;
		    temp = max(temp,0.);
/* Computing 2nd power */
		    d__1 = work[i__] / work[i__ + *lwsize];
		    temp2 = temp * .05f * (d__1 * d__1) + 1.;
		    if (temp2 == 1.) {
			i__2 = *m - k;
			work[i__] = dnrm2(&i__2, &a[k + 1 + (*offset + i__) *
				 a_dim1], &c__1);
			work[i__ + *lwsize] = work[i__];
		    } else {
			work[i__] *= sqrt(temp);
		    }
		}
/* L20: */
	    }
	}
    } else {

/*           Reject all remaining columns in pivot window. */
/*           ============================================ */

	goto L2000;
    }

/*     End while. */

    goto L1000;
L2000:
    return 0;

/*     End of DGEQPW */

} /* dgeqpw_ */

