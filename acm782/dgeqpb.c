/* dgeqpb.f -- translated by f2c (version 20041007).
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
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__3 = 3;

/* Subroutine */ int dgeqpb(integer *job, integer *m, integer *n, integer *k,
	 doublereal *a, integer *lda, doublereal *c__, integer *ldc, integer *
	jpvt, doublereal *ircond, doublereal *orcond, integer *rank,
	doublereal *svlues, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, kb, nb, kk, mn;
    static doublereal smin, mxnm;
    static logical block;
    static doublereal rcond;
    static integer itemp;
    static integer wkmin, mvidx, strej, wsize;
    static integer accptd;
    static integer lacptd;
    static integer nllity, lwsize;


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
/*     $Date: 96/12/30 16:59:10 $ */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEQPB computes a QR factorization */
/*       A*P = Q*[ R11 R12 ] */
/*               [  0  R22 ] */
/*  of a real m by n matrix A. The permutation P is */
/*  chosen with the goal to reveal the rank of A by a */
/*  suitably dimensioned trailing submatrix R22 with norm(R22) */
/*  being small. */

/*  Arguments */
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

/*  JPVT    (output) INTEGER array, dimension (N) */
/*          JPVT(I) = K <==> Column K of A has been permuted into */
/*                           position I in AP. */
/*          JPVT(1:RANK) contains the indices of the columns considered */
/*          linearly independent. */
/*          JPVT(RANK+1:N) contains the indices of the columns considered */
/*          linearly dependent from the previous ones. */

/*  IRCOND  (input) DOUBLE PRECISION */
/*          1/IRCOND specifies an upper bound on the condition number */
/*          of R11. If IRCOND == 0, IRCOND = machine precision is chosen */
/*          as default. IRCOND must be >= 0. */

/*  ORCOND  (output) DOUBLE PRECISION */
/*          1/ORCOND is an estimate for the condition number of R11. */

/*  RANK    (output) INTEGER */
/*          RANK is an estimate for the numerical rank of A with respect */
/*          to the threshold 1/IRCOND in the sense that */
/*               RANK = arg_max(cond_no(R(1:r,1:r))<1/IRCOND) */
/*          This may be an underestimate of the rank if the leading */
/*          columns were not well-conditioned. */

/*  SVLUES  (output) DOUBLE PRECISION array, dimension (4) */
/*          On exit, SVLUES contains estimates of some of the singular */
/*          values of the triangular factor R. */
/*          SVLUES(1): largest singular value of R(1:RANK,1:RANK) */
/*          SVLUES(2): smallest singular value of R(1:RANK,1:RANK) */
/*          SVLUES(3): smallest singular value of R(1:RANK+1,1:RANK+1) */
/*          SVLUES(4): smallest singular value of R */
/*          If the triangular factorization is a rank-revealing one */
/*          (which will be the case if the leading columns were well- */
/*          conditioned), then SVLUES(1) will also be an estimate for */
/*          the largest singular value of A, SVLUES(2) and SVLUES(3) */
/*          will be estimates for the RANK-th and (RANK+1)-st singular */
/*          value of A, and SVLUES(4) will be an estimate for the */
/*          smallest singular value of A. */
/*          By examining these values, one can confirm that the rank is */
/*          well defined with respect to the threshold chosen. */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK) */
/*          On exit: WORK(1) is the size of the storage array needed */
/*                   for optimal performance */

/*  LWORK   (input) INTEGER */
/*          The dimension of array WORK. */
/*          If JOB=1: */
/*             The unblocked strategy requires that: */
/*                 LWORK >= 2*MN+3*N. */
/*             The block algorithm requires that: */
/*                 LWORK >= 2*MN+N*NB. */
/*          If JOB<>1: */
/*             The unblocked strategy requires that: */
/*                 LWORK >= 2*MN+2*N+MAX(K,N). */
/*             The block algorithm requires that: */
/*                 LWORK >= 2*MN+NB*NB+NB*MAX(K,N). */
/*          Where MN = min(M,N) and NB is the maximum of blocksize */
/*          used within xGEQRF and blocksize used within xORMQR. */
/*          In both cases, the minimum required workspace is the */
/*          one for the unblocked strategy. */

/*  INFO    (output) INTEGER */
/*          = 0: Successful exit */
/*          < 0: If INFO = -i, the i-th argument had an illegal value */

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
    --svlues;
    --work;

    /* Function Body */
    mn = min(*m,*n);

/*     Compute the minimum required workspace size. */

    if (*job == 1) {
	wkmin = (mn << 1) + *n * 3;
    } else {
	wkmin = (mn << 1) + (*n << 1) + max(*n,*k);
    }

/*     Test input arguments */
/*     ==================== */

    *info = 0;
    if (*job < 1 || *job > 3) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*k < 0) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -6;
    } else if (*job == 1 && *ldc < 1 || *job == 2 && *ldc < max(1,*m) || *job
	    == 3 && *ldc < max(1,*k)) {
	*info = -8;
    } else if (*ircond < 0.) {
	*info = -10;
    } else if (*lwork < max(1,wkmin)) {
	*info = -15;
    }
    if (*info == 0 || *info == -15) {

/*        Compute the optimal workspace size. */

	if (*job == 1) {
	    nb = ilaenv(&c__1, "DGEQRF", " ", m, n, &c__0, &c__0, (ftnlen)6,
		    (ftnlen)1);
/* Computing MAX */
	    i__1 = *n * 3, i__2 = *n * nb;
	    wsize = (mn << 1) + max(i__1,i__2);
	} else {
/* Computing MAX */
	    i__1 = ilaenv(&c__1, "DGEQRF", " ", m, n, &c__0, &c__0, (ftnlen)
		    6, (ftnlen)1), i__2 = ilaenv(&c__1, "DORMQR", " ", m, n,
		    &c__0, &c__0, (ftnlen)6, (ftnlen)1);
	    nb = max(i__1,i__2);
/* Computing MAX */
	    i__1 = (mn << 1) + (*n << 1) + max(*n,*k), i__2 = (mn << 1) + nb *
		     nb + nb * max(*n,*k);
	    wsize = max(i__1,i__2);
	}
	work[1] = (doublereal) wsize;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("DGEQPB", &i__1);
	return 0;
    }

/*     Initialization of vector JPVT. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jpvt[j] = j;
/* L70: */
    }

/*     Quick return if possible */

    if (mn == 0) {
	*rank = 0;
	*orcond = 0.;
	svlues[1] = 0.;
	svlues[2] = 0.;
	svlues[3] = 0.;
	svlues[4] = 0.;
	return 0;
    }

/*     Check whether Threshold for condition number was supplied. */
/*     If not, choose machine precision as default for RCOND. */

    if (*ircond > 0.) {
	rcond = *ircond;
    } else {
	rcond = dlamch("Epsilon");
    }

/*     Determine the allowed block size for the given workspace */
/*     and whether to use blocked code at all. */

    if (*lwork < wsize) {
	if (*job == 1) {
	    nb = (*lwork - (mn << 1)) / *n;
	} else {
/* Computing 2nd power */
	    i__1 = max(*k,*n);
	    itemp = (integer) sqrt((doublereal) (i__1 * i__1 + (*lwork << 2)
		    - (mn << 3)));
	    nb = (itemp - max(*k,*n)) / 2;
	}
    }

    block = nb > 1 && nb >= ilaenv(&c__2, "DGEQRF", " ", m, n, &c__0, &c__0,
	    (ftnlen)6, (ftnlen)1) && mn >= ilaenv(&c__3, "DGEQRF", " ", m, n,
	     &c__0, &c__0, (ftnlen)6, (ftnlen)1);

/*     The size of the pivot window is chosen to be NB + NLLITY */
/*     for the blocked algorithm. */

/* Computing MIN */
/* Computing MAX */
    i__3 = 10, i__4 = nb / 2 + *n * 5 / 100;
    i__1 = mn, i__2 = max(i__3,i__4);
    nllity = min(i__1,i__2);

/*     *************************************************** */
/*     * Move column with largest residual norm up front * */
/*     *************************************************** */

    i__1 = *lwork - (mn << 1);
    dgeqpc(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &c__1, &
	    c__0, &rcond, &lacptd, &jpvt[1], &work[1], &work[mn + 1], &svlues[
	    1], &mxnm, &work[(mn << 1) + 1], &i__1);
    if (lacptd == 1) {
	if (lacptd == mn) {
	    *rank = 1;
	    *orcond = svlues[2] / svlues[1];
	    goto L30;
	} else {
	    smin = svlues[2];
	}
    } else {
	*rank = 0;
	*orcond = 0.;
	svlues[1] = 0.;
	svlues[2] = 0.;
	svlues[3] = 0.;
	svlues[4] = 0.;
	goto L30;
    }

/*     **************************** */
/*     * Factor remaining columns * */
/*     **************************** */

    if (block) {

/*        *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */
/*        * Using blocked code with restricted pivoting strategy  * */
/*        *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* */

	strej = *n + 1;
	kk = 2;

L10:
	if (kk >= strej || kk > mn) {
	    goto L20;
	}

/*           invariant: A(:,KK) is the first column in currently */
/*           considered block column. */

/* Computing MIN */
/* Computing MIN */
	i__3 = mn + 1;
	i__1 = nb, i__2 = min(i__3,strej) - kk;
	kb = min(i__1,i__2);

/*           The goal now is to find "KB" independent columns */
/*           among the remaining STREJ-KK not yet rejected columns. */

/* Computing MIN */
	i__1 = strej - kk, i__2 = kb + nllity;
	lwsize = min(i__1,i__2);
	i__1 = kk - 1;
	dgeqpw(m, &lwsize, &kb, &i__1, &lacptd, &a[a_offset], lda, &jpvt[1],
		&rcond, &work[mn + 1], &smin, &mxnm, &work[1], &work[(mn << 1)
		 + 1]);
	if (lacptd > 0) {

/*              Accumulate Householder vectors in a block reflector. */

	    i__1 = *m - kk + 1;
	    dlarft("Forward", "Columnwise", &i__1, &lacptd, &a[kk + kk *
		    a_dim1], lda, &work[kk], &work[(mn << 1) + 1], &lacptd);

/*              Apply block reflector to A(KK:M,KK+LWSIZE:N). */

	    if (kk + lwsize <= *n) {
		i__1 = *m - kk + 1;
		i__2 = *n - kk - lwsize + 1;
		i__3 = *n - kk - lwsize + 1;
		dlarfb("Left", "Transpose", "Forward", "Columnwise", &i__1, &
			i__2, &lacptd, &a[kk + kk * a_dim1], lda, &work[(mn <<
			 1) + 1], &lacptd, &a[kk + (kk + lwsize) * a_dim1],
			lda, &work[(mn << 1) + lacptd * lacptd + 1], &i__3);
	    }

/*              Apply block reflector to the corresponding part */
/*              of matrix C. */

	    if (*job == 2 && *k > 0) {

/*                 Apply it to matrix C(KK:M,1:K) from the left. */

		i__1 = *m - kk + 1;
		dlarfb("Left", "Transpose", "Forward", "Columnwise", &i__1,
			k, &lacptd, &a[kk + kk * a_dim1], lda, &work[(mn << 1)
			 + 1], &lacptd, &c__[kk + c_dim1], ldc, &work[(mn <<
			1) + lacptd * lacptd + 1], k);
	    } else if (*job == 3 && *k > 0) {

/*                 Apply the transpose of it to matrix C(1:K,KK:M) */
/*                 from the right. */

		i__1 = *m - kk + 1;
		dlarfb("Right", "No Transpose", "Forward", "Columnwise", k, &
			i__1, &lacptd, &a[kk + kk * a_dim1], lda, &work[(mn <<
			 1) + 1], &lacptd, &c__[kk * c_dim1 + 1], ldc, &work[(
			mn << 1) + lacptd * lacptd + 1], k);
	    }
	}

/*           Move rejected columns to the end if there is space. */

	if (lacptd < kb) {
	    if (strej <= kk + lwsize) {
		strej = kk + lacptd;
	    } else {
		mvidx = strej;
/* Computing MIN */
		i__2 = kk + lwsize - 1, i__3 = strej - lwsize + lacptd - 1;
		i__1 = min(i__2,i__3);
		for (i__ = kk + lacptd; i__ <= i__1; ++i__) {
		    --mvidx;
		    dswap(m, &a[i__ * a_dim1 + 1], &c__1, &a[mvidx * a_dim1
			    + 1], &c__1);
		    itemp = jpvt[i__];
		    jpvt[i__] = jpvt[mvidx];
		    jpvt[mvidx] = itemp;
/* L40: */
		}
		strej = mvidx;
	    }
	}
	kk += lacptd;
	goto L10;
L20:
	accptd = kk - 1;
	svlues[1] = dlasmx(&accptd) * mxnm;
	svlues[2] = smin;
	if (accptd < mn) {

/*           Process rejected columns. */

	    i__1 = mn - kk + 1;
	    i__2 = kk - 1;
	    i__3 = *lwork - (mn << 1);
	    dgeqpc(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &
		    i__1, &i__2, &rcond, &lacptd, &jpvt[1], &work[1], &work[
		    mn + 1], &svlues[1], &mxnm, &work[(mn << 1) + 1], &i__3);
	    *rank = accptd + lacptd;
	} else {
	    *rank = accptd;
	    svlues[3] = smin;
	    svlues[4] = smin;
	}
    } else {

/*        *-*-*-*-*-*-*-*-*-*-*-*-* */
/*        * using unblocked code  * */
/*        *-*-*-*-*-*-*-*-*-*-*-*-* */

	accptd = 1;
	i__1 = mn - accptd;
	i__2 = *lwork - (mn << 1);
	dgeqpc(job, m, n, k, &a[a_offset], lda, &c__[c_offset], ldc, &i__1, &
		accptd, &rcond, &lacptd, &jpvt[1], &work[1], &work[mn + 1], &
		svlues[1], &mxnm, &work[(mn << 1) + 1], &i__2);
	*rank = accptd + lacptd;

    }
    *orcond = svlues[2] / svlues[1];

/*     Nullify the lower part of matrix A. */

L30:
    i__1 = mn;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = 0.;
/* L60: */
	}
/* L50: */
    }

    work[1] = (doublereal) wsize;
    return 0;

/*     End of DGEQPB */

} /* dgeqpb_ */

