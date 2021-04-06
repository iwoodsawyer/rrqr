/* dgeqpc.f -- translated by f2c (version 20041007).
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
static integer c__2 = 2;

/* Subroutine */ int dgeqpc(integer *job, integer *m, integer *n, integer *k,
	 doublereal *a, integer *lda, doublereal *c__, integer *ldc, integer *
	dsrd, integer *offset, doublereal *ircond, integer *lacptd, integer *
	jpvt, doublereal *tau, doublereal *x, doublereal *svlues, doublereal *
	mxnm, doublereal *work, integer *lwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, mn;
    static doublereal aii;
    static integer pvt, info;
    static doublereal sine, temp, smin, smax;
    static doublereal temp2;
    static integer lasti, itemp;
    static doublereal cosine;
    static doublereal sminpr, smaxpr;


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
/*     $Date: 96/12/30 16:59:11 $ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose: */
/*  ======= */

/*  DGEQPC continues a partial QR factorization of A. If */
/*  A(1:OFFSET,1:OFFSET) has been reduced to upper triangular */
/*  form, then DGEQPC applies the traditional column pivoting */
/*  strategy to identify DSRD more independent columns of A with */
/*  the restriction that the condition number of the leading */
/*  triangle of A should not be larger than 1/IRCOND.  If */
/*  LACPTD ( <= DSRD) such columns are found, then the condition */
/*  number of */
/*     A(1:OFFSET+LACPTD,1:OFFSET+LACPTD) is less than 1/IRCOND. */
/*  If LACPTD < DSRD, then the QR factorization of A is completed, */
/*  otherwise only DSRD new steps were performed. */

/*  Arguments: */
/*  ========= */

/*  JOB     (input) INTEGER */
/*          The job to do: */
/*          = 1: The orthogonal transformations needed in the */
/*               triangularization are only applied to matrix A. */
/*               Thus, matrix C is not updated. */
/*          = 2: The same orthogonal transformations needed in the */
/*               triangularization of matrix A are applied to */
/*               matrix C from the left. */
/*               That is, if Q'*A*P=R, then C := Q'*C. */
/*               In this case, matrix C is m-by-k. */
/*          = 3: The transpose of the orthogonal transformations needed */
/*               in the triangularization of matrix A are applied */
/*               to matrix C from the right. */
/*               That is, if Q'*A*P=R, then C := C*Q. */
/*               In this case, matrix C is k-by-m. */
/*          In these three cases, the permutations are always stored */
/*          in vector JPVT. */

/*  M       (input) INTEGER */
/*          The number of rows of matrices A. M >= 0. */
/*          If JOB=2, M is the number of rows of matrix C. */
/*          If JOB=3, M is the number of columns of matrix C. */

/*  N       (input) INTEGER */
/*          The number of columns of matrix A.  N >= 0. */

/*  K       (input) INTEGER */
/*          It defines the dimension of matrix C. K >= 0. */
/*          If JOB=2, K is the number of columns of matrix C. */
/*          If JOB=3, K is the number of rows of matrix C. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the m by n matrix A. */
/*          On exit, the upper triangle of the array contains the */
/*          min(m,n) by n upper trapezoidal matrix R; the lower triangle */
/*          array is filled with zeros. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of array A. LDA >= max(1,M). */

/*  C       (input/output) DOUBLE PRECISION array, dimension */
/*                ( LDC, K ) if JOB=2. */
/*                ( LDC, M ) if JOB=3. */
/*          If argument JOB asks, all the orthogonal transformations */
/*          applied to matrix A are also applied to matrix C. */

/*  LDC     (input) INTEGER */
/*          The leading dimension of array C. */
/*          If JOB=2, then LDC >= MAX(1,M). */
/*          If JOB=3, then LDC >= MAX(1,K). */

/*  DSRD    (input) INTEGER */
/*          The number of independent columns one would like to */
/*          extract. */

/*  OFFSET  (input) INTEGER */
/*          A(1:OFFSET,1:OFFSET) has already been factored. */
/*          OFFSET >= 0. */

/*  IRCOND  (input) DOUBLE PRECISION */
/*          1/IRCOND is threshold for condition number. */

/*  LACPTD  (output) INTEGER */
/*          The number of additional columns that were identified */
/*          as independent. */

/*  JPVT    (input/output) INTEGER array, dimension (N) */
/*          If JPVT(I) = K, then the Ith column of the permuted */
/*          A was the Kth column of the original A. */

/*  TAU     (input/output) DOUBLE PRECISION array, dimension (MIN(M,N)) */
/*          Further details of the orthogonal matrix Q (see A). */

/*  X       (input/output) DOUBLE PRECISION array, dimension (MIN(M,N)) */
/*          On entry: X(1:OFFSET) contains an approximate smallest */
/*          left singular vector of A(1:OFFSET,1:OFFSET). */
/*          On exit: X(1:OFFSET+LACPTD) contains an approximate */
/*          smallest left singular vector of */
/*          A(1:OFFSET+LACPTD,1:OFFSET+LACPTD). */

/*  SVLUES  (input/output) DOUBLE PRECISION array, dimension(4) */
/*          The estimates of the singular values. */
/*          On entry: SVLUES(1) = sigma_max(A(1:M,1:N)) */
/*                    SVLUES(2) = sigma_min(A(1:OFFSET,1:OFFSET)) */
/*          On exit: SVLUES(1) = sigma_max(A(1:M,1:N)) */
/*                   SVLUES(2) = sigma_r(B) */
/*                   SVLUES(3) = sigma(min(r+1,min(m,n)))(B) */
/*                   SVLUES(4) = sigma_min(A) */
/*          where r = OFFSET+LACPTD and B = A(1:r,1:r) */

/*  MXNM    (input/output) DOUBLE PRECISION */
/*          On entry: norm of largest column in A(1:OFFSET,1:OFFSET) */
/*          On exit: norm of largest column in */
/*                   A(1:J,1:J) where J = OFFSET+LACPTD */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */

/*  LWORK   (input) INTEGER */
/*             MAX( 1, 3*N, N*NB )                          if JOB=1, or */
/*             MAX( 1, 2*N + MAX( N, K ), MAX( N, K)*NB )   otherwise. */
/*          where NB is the maximum of blocksize used within xGEQRF and */
/*          blocksize used within xORMQR. */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(1) H(2) . . . H(n) */

/*  Each H(i) has the form */

/*     H = I - tau * v * v' */

/*  where tau is a real scalar, and v is a real vector with */
/*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i). */

/*  The matrix P is represented in jpvt as follows: If */
/*     jpvt(j) = i */
/*  then the jth column of P is the ith canonical unit vector. */

/*  ===================================================================== */

/*     .. Parameters .. */

/*     Indices into the `svlues' array. */

/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --jpvt;
    --tau;
    --x;
    --svlues;
    --work;

    /* Function Body */
    mn = min(*m,*n);
    *lacptd = 0;
    if (*offset > 0) {
	smax = svlues[1];
	smin = svlues[2];
    }

/*     Initialize partial column norms. The first n entries of */
/*     work store the exact column norms. */

    i__1 = *n;
    for (i__ = *offset + 1; i__ <= i__1; ++i__) {
	i__2 = *m - *offset;
	work[i__] = dnrm2(&i__2, &a[*offset + 1 + i__ * a_dim1], &c__1);
	work[*n + i__] = work[i__];
/* L10: */
    }

/*     Compute factorization. */

/* Computing MIN */
    i__1 = mn, i__2 = *offset + *dsrd;
    lasti = min(i__1,i__2);
    i__1 = lasti;
    for (i__ = *offset + 1; i__ <= i__1; ++i__) {

/*        Determine ith pivot column and swap if necessary. */

	i__2 = *n - i__ + 1;
	pvt = i__ - 1 + idamax(&i__2, &work[i__], &c__1);
	if (pvt != i__) {
	    dswap(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		    c__1);
	    itemp = jpvt[pvt];
	    jpvt[pvt] = jpvt[i__];
	    jpvt[i__] = itemp;
	    work[pvt] = work[i__];
	    work[*n + pvt] = work[*n + i__];
	}

/*        Generate elementary reflector H(i). */

	if (i__ < *m) {
	    i__2 = *m - i__ + 1;
	    dlarfg(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * a_dim1],
		     &c__1, &tau[i__]);
	} else {
	    dlarfg(&c__1, &a[*m + *m * a_dim1], &a[*m + *m * a_dim1], &c__1,
		    &tau[*m]);
	}

/*        Apply elementary reflector H(I) to the corresponding block */
/*        of matrices A and C. */

	aii = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = 1.;
	if (i__ < *n) {

/*           Apply H(I) to A(I:M,I+1:N) from the left. */

	    i__2 = *m - i__ + 1;
	    i__3 = *n - i__;
	    dlarf("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &tau[
		    i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[(*n << 1)
		    + 1]);
	}
	if (*job == 2 && *k > 0) {

/*           Apply H(I) to C(I:M,1:K) from the left. */

	    i__2 = *m - i__ + 1;
	    dlarf("Left", &i__2, k, &a[i__ + i__ * a_dim1], &c__1, &tau[i__],
		     &c__[i__ + c_dim1], ldc, &work[(*n << 1) + 1]);
	} else if (*job == 3 && *k > 0) {

/*           Apply the transpose of H(I) to C(1:K,I:M) from the right. */

	    i__2 = *m - i__ + 1;
	    dlarf("Right", k, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[i__]
		    , &c__[i__ * c_dim1 + 1], ldc, &work[(*n << 1) + 1]);
	}
	a[i__ + i__ * a_dim1] = aii;

/*        Update partial column norms. */

	if (i__ < lasti) {
	    i__2 = *n;
	    for (j = i__ + 1; j <= i__2; ++j) {
		if (work[j] != 0.) {
/* Computing 2nd power */
		    d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)) / work[j];
		    temp = 1. - d__2 * d__2;
		    temp = max(temp,0.);
/* Computing 2nd power */
		    d__1 = work[j] / work[*n + j];
		    temp2 = temp * .05f * (d__1 * d__1) + 1.;
		    if (temp2 == 1.) {
			i__3 = *m - i__;
			work[j] = dnrm2(&i__3, &a[i__ + 1 + j * a_dim1], &
				c__1);
			work[*n + j] = work[j];
		    } else {
			work[j] *= sqrt(temp);
		    }
		}
/* L30: */
	    }
	}

/*        Check new column for independence. */

	if (i__ == 1) {
	    *mxnm = (d__1 = a[a_dim1 + 1], abs(d__1));
	    smin = *mxnm;
	    smax = *mxnm;
	    x[1] = 1.;
	    if (*mxnm > 0.) {
		*lacptd = 1;
	    } else {
		svlues[3] = smin;
		goto L50;
	    }
	} else {
	    smaxpr = dlasmx(&i__) * *mxnm;
	    d__1 = smaxpr * *ircond;
	    if (dlauc1(&i__, &x[1], &smin, &a[i__ * a_dim1 + 1], &a[i__ +
		    i__ * a_dim1], &d__1)) {

/*              Column accepted. */

		smax = smaxpr;
		++(*lacptd);
	    } else {

/*              Column rejected. */

		goto L50;
	    }
	}
/* L20: */
    }

L50:
    svlues[1] = smax;
    svlues[2] = smin;
    if (*lacptd == *dsrd) {

/*        DSRD independent columns have been found. */

	svlues[3] = smin;
	svlues[4] = smin;
    } else {

/*        All remaining columns rejected. */

	i__ = *offset + *lacptd + 1;
	if (i__ < mn) {

/*           Factor remaining columns. */

	    i__1 = *m - i__;
	    i__2 = *n - i__;
	    dgeqrf(&i__1, &i__2, &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &tau[
		    i__ + 1], &work[1], lwork, &info);

/*           Apply the transformations computed in DGEQRF to the */
/*           corresponding part of matrix C. */

	    if (*job == 2 && *k > 0) {

/*              Apply them to C(I+1:M,1:K) from the left. */

		i__1 = *m - i__;
		i__2 = mn - i__;
		dormqr("Left", "Transpose", &i__1, k, &i__2, &a[i__ + 1 + (
			i__ + 1) * a_dim1], lda, &tau[i__ + 1], &c__[i__ + 1
			+ c_dim1], ldc, &work[1], lwork, &info);
	    } else if (*job == 3 && *k > 0) {

/*              Apply the transpose of them to C(1:K,I+1:M) from the */
/*              right. */

		i__1 = *m - i__;
		i__2 = mn - i__;
		dormqr("Right", "No Transpose", k, &i__1, &i__2, &a[i__ + 1
			+ (i__ + 1) * a_dim1], lda, &tau[i__ + 1], &c__[(i__
			+ 1) * c_dim1 + 1], ldc, &work[1], lwork, &info);
	    }
	}

/*        Use incremental condition estimation to get an estimate */
/*        of the smallest singular value. */

/* Computing MAX */
	i__1 = 2, i__2 = *offset + *lacptd + 1;
	i__3 = mn;
	for (i__ = max(i__1,i__2); i__ <= i__3; ++i__) {
	    i__1 = i__ - 1;
	    dlaic1(&c__2, &i__1, &x[1], &smin, &a[i__ * a_dim1 + 1], &a[i__
		    + i__ * a_dim1], &sminpr, &sine, &cosine);
	    i__1 = i__ - 1;
	    dscal(&i__1, &sine, &x[1], &c__1);
	    x[i__] = cosine;
	    smin = sminpr;
	    if (i__ == *offset + *lacptd + 1) {
		svlues[3] = smin;
	    }
/* L60: */
	}
	svlues[4] = smin;
    }
    return 0;

/*     End of DGEQPC */

} /* dgeqpc_ */

